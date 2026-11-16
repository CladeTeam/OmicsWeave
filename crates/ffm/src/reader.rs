//! Random access to an FFM file over any [`ReadAt`] source.

use std::ops::Deref;
use std::sync::Arc;

use omicsweave_core::ReadAt;
use omicsweave_core::cache::{BlockCache, BlockKey, CacheLease, SourceToken};

use crate::alphabet::Alphabet;
use crate::block::Block;
use crate::error::{Result, invalid};
use crate::layout::{Cursor, DirEntry, HEADER_LEN, Header, Span, parse_directory, unseal};
use crate::meta::{IdTable, MotifIndex};
use crate::options::{ReaderOptions, allocate};
use crate::schema::Schema;
use crate::stream::next_stream;
use crate::view::MotifView;

/// Block-cache kind of FFM blocks; `block` is the block number.
const CACHE_KIND: u32 = 7;

/// A decoded block, owned or leased from a shared cache.
#[derive(Clone)]
enum BlockHandle {
    Owned(Arc<Block>),
    Leased(CacheLease<Block>),
}

impl Deref for BlockHandle {
    type Target = Block;

    fn deref(&self) -> &Block {
        match self {
            BlockHandle::Owned(b) => b,
            BlockHandle::Leased(b) => b,
        }
    }
}

/// Opens an FFM file once and serves motifs, IDs and annotations by motif
/// number or ID. All lookups are O(1) (expected, for IDs); a block is read
/// whole and decoded once per shared-cache residency. Without a cache,
/// random accesses decode independently; [`Self::motifs`] holds the current
/// block for the duration of one sequential traversal.
pub struct FfmReader<R> {
    source: R,
    header: Header,
    alphabet: Alphabet,
    schema: Arc<Schema>,
    ids: IdTable,
    directory: Vec<DirEntry>,
    index: MotifIndex,
    cache: Option<Arc<BlockCache>>,
    namespace: SourceToken,
    options: ReaderOptions,
}

impl<R: ReadAt> FfmReader<R> {
    /// Reads the header and all metadata (two reads), checks them and
    /// builds the ID lookup table.
    pub fn open(source: R) -> Result<Self> {
        Self::open_with_options(source, ReaderOptions::default())
    }

    pub fn open_with_options(source: R, options: ReaderOptions) -> Result<Self> {
        let file_len = source.len();
        let mut head = [0u8; HEADER_LEN];
        if file_len < HEADER_LEN as u64 {
            return Err(invalid("FFM header", "not an FFM file (too short)"));
        }
        source.read_exact_at(&mut head, 0)?;
        let header = Header::parse(&head, file_len)?;
        let meta_span = header.metadata();
        let meta_len = usize::try_from(meta_span.len)
            .map_err(|_| invalid("FFM metadata", "length does not fit memory address space"))?;
        let mut meta = allocate(meta_len, options.max_metadata_bytes, "FFM metadata")?;
        source.read_exact_at(&mut meta, meta_span.offset)?;
        let section = |s: Span| {
            let a = (s.offset - meta_span.offset) as usize;
            &meta[a..a + s.len as usize]
        };

        let alphabet = Alphabet::parse(unseal(section(header.alphabet), "FFM alphabet")?)?;
        let schema = Schema::parse(unseal(section(header.schema), "FFM schema")?)?;

        let mut r = Cursor::new(section(header.ids), "FFM IDs");
        let ids_stream = next_stream(&mut r, "FFM IDs")?;
        let min_ids = (header.motif_count + 1) * 4 + header.motif_count;
        if (ids_stream.decoded_len as u64) < min_ids {
            return Err(invalid("FFM IDs", "stream too short for motif count"));
        }
        let ids = ids_stream.bytes("FFM IDs", options.max_stream_bytes)?;
        r.finish()?;
        let ids = IdTable::parse_stored(
            &ids,
            header.motif_count as usize,
            options.max_resident_bytes,
        )?;

        let directory = parse_directory(section(header.directory), &header, file_len)?;

        let mut r = Cursor::new(section(header.index), "FFM index");
        let index_stream = next_stream(&mut r, "FFM index")?;
        if index_stream.decoded_len as u64 != 4 * header.motif_count {
            return Err(invalid("FFM index", "length is not 4 x motif count"));
        }
        let columns = index_stream.numbers::<u32>("FFM index", options.max_stream_bytes)?;
        r.finish()?;
        if columns.len() as u64 != 4 * header.motif_count {
            return Err(invalid("FFM index", "length is not 4 x motif count"));
        }
        let index = MotifIndex::new(columns, &directory, !schema.is_empty())?;

        Ok(FfmReader {
            source,
            header,
            alphabet,
            schema: Arc::new(schema),
            ids,
            directory,
            index,
            cache: None,
            namespace: SourceToken::new(),
            options,
        })
    }

    /// Number of motifs.
    pub fn len(&self) -> usize {
        self.header.motif_count as usize
    }

    pub fn is_empty(&self) -> bool {
        self.header.motif_count == 0
    }

    pub fn block_count(&self) -> usize {
        self.directory.len()
    }

    /// Minor format version of the file.
    pub fn minor_version(&self) -> u16 {
        self.header.minor
    }

    pub fn alphabet(&self) -> &Alphabet {
        &self.alphabet
    }

    pub fn schema(&self) -> &Schema {
        &self.schema
    }

    /// The ID of motif `i`.
    pub fn id(&self, i: usize) -> Option<&str> {
        self.ids.get(i)
    }

    /// The motif number of an ID.
    pub fn find(&self, id: &str) -> Option<usize> {
        self.ids.find(id)
    }

    /// Shares decoded blocks through `cache`; `None` decodes independently.
    /// Cached sessions require source contents to remain unchanged, as with
    /// other readers using [`BlockCache`].
    pub fn set_cache(&mut self, cache: Option<Arc<BlockCache>>) {
        self.cache = cache;
    }

    pub fn cache(&self) -> Option<&Arc<BlockCache>> {
        self.cache.as_ref()
    }

    pub fn cache_namespace(&self) -> SourceToken {
        self.namespace
    }

    /// The model of motif `i`.
    pub fn motif(&self, i: usize) -> Result<MotifRef> {
        self.check_motif(i)?;
        let b = self.index.block(i);
        let block = self.block(b)?;
        Ok(MotifRef {
            block,
            record: i - self.index.block_first(b),
            alphabet_len: self.alphabet.len(),
        })
    }

    /// The annotation of motif `i`; the block's annotation stream is
    /// decompressed on first use.
    pub fn annotation(&self, i: usize) -> Result<AnnotationRef> {
        self.check_motif(i)?;
        let b = self.index.block(i);
        let block = self.block(b)?;
        block.annotation(&self.schema)?;
        Ok(AnnotationRef {
            block,
            record: i - self.index.block_first(b),
            schema: self.schema.clone(),
        })
    }

    /// Every motif in order; each block is decoded once.
    pub fn motifs(&self) -> impl Iterator<Item = Result<MotifRef>> + '_ {
        let mut current: Option<(usize, BlockHandle)> = None;
        (0..self.len()).map(move |i| {
            let b = self.index.block(i);
            if current.as_ref().is_none_or(|(at, _)| *at != b) {
                current = Some((b, self.block(b)?));
            }
            Ok(MotifRef {
                block: current.as_ref().expect("current traversal block").1.clone(),
                record: i - self.index.block_first(b),
                alphabet_len: self.alphabet.len(),
            })
        })
    }

    /// Consumes the reader, fully validates annotations and packs all models
    /// into resident arenas. The result retains neither source nor block cache.
    pub fn into_set(self) -> Result<crate::MotifSet> {
        let base = self.ids.heap_bytes() + self.alphabet.heap_bytes() + self.schema.heap_bytes();
        let limit = self
            .options
            .max_resident_bytes
            .checked_sub(base)
            .ok_or_else(|| invalid("FFM collection", "metadata exceeds resident limit"))?;
        let mut arena = crate::collection::MotifArena::default();
        for (b, entry) in self.directory.iter().enumerate() {
            let block = self.decode(b)?;
            block.annotation(&self.schema)?;
            for j in 0..entry.record_count as usize {
                let model = MotifView::from_block(&block, j, self.alphabet.len());
                arena.push_view(
                    model,
                    block.annotation_of(&self.schema, j)?,
                    self.schema.len(),
                    limit,
                )?;
            }
        }
        Ok(crate::MotifSet::from_parts(
            self.alphabet,
            Arc::unwrap_or_clone(self.schema),
            self.ids,
            arena,
        ))
    }

    fn check_motif(&self, i: usize) -> Result<()> {
        if i >= self.len() {
            return Err(invalid(
                "FFM reader",
                format!("motif {i} out of {} motifs", self.len()),
            ));
        }
        Ok(())
    }

    fn block(&self, b: usize) -> Result<BlockHandle> {
        if let Some(cache) = &self.cache {
            let key = BlockKey::new(self.namespace, CACHE_KIND, [0, 0], b as u64, 0);
            let lease = cache.get_or_load(key, || {
                let block = self.decode(b)?;
                let heap = block.heap_bytes();
                Ok((block, heap))
            })?;
            return Ok(BlockHandle::Leased(lease));
        }
        Ok(BlockHandle::Owned(Arc::new(self.decode(b)?)))
    }

    fn decode(&self, b: usize) -> Result<Block> {
        let e = self.directory[b];
        let mut bytes = allocate(e.length as usize, self.options.max_block_bytes, "FFM block")?;
        self.source.read_exact_at(&mut bytes, e.offset)?;
        let starts = self.index.block_offsets(b, e.record_count);
        Block::decode(
            &bytes,
            &starts,
            self.alphabet.len(),
            self.options.max_stream_bytes,
        )
    }
}

/// A motif's model, holding its decoded block.
#[derive(Clone)]
pub struct MotifRef {
    block: BlockHandle,
    record: usize,
    alphabet_len: usize,
}

impl MotifRef {
    pub fn view(&self) -> MotifView<'_> {
        MotifView::from_block(&self.block, self.record, self.alphabet_len)
    }
}

/// A motif's annotation values, one per schema key; empty means missing.
#[derive(Clone)]
pub struct AnnotationRef {
    block: BlockHandle,
    record: usize,
    schema: Arc<Schema>,
}

impl AnnotationRef {
    /// The record's values joined by tabs.
    pub fn as_str(&self) -> &str {
        self.block
            .annotation_of(&self.schema, self.record)
            .expect("annotation decoded when the reference was made")
    }

    pub fn values(&self) -> impl Iterator<Item = &str> {
        let keyed = !self.schema.is_empty();
        keyed
            .then(|| self.as_str().split('\t'))
            .into_iter()
            .flatten()
    }

    /// The value of key `i`.
    pub fn get(&self, i: usize) -> Option<&str> {
        self.values().nth(i)
    }

    /// The value of the key named `name`.
    pub fn get_by_name(&self, name: &str) -> Option<&str> {
        self.get(self.schema.index(name)?)
    }

    pub fn value(&self, i: usize) -> Option<crate::AnnotationValue<'_>> {
        self.schema.value(i, self.get(i)?)
    }
}
