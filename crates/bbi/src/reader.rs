//! [`BbiReader`]: the standard access to a bigWig or bigBed file, composing
//! the crate's steps over a caller-provided byte source.

use std::ops::Range;
use std::sync::Arc;

use omicsweave_core::cache::{BlockCache, BlockKey, CacheLease, SourceToken};

use omicsweave_core::orient::orient_values;
use omicsweave_core::{
    BinAccumulator, BinLayout, BinSource, BinSpec, IntervalEntry, IntervalRead, IntervalSize,
    IntervalsOut, Locus, ReadAt, Region, Span, TrackEntry, TrackRead, ZoomPiece, check_output,
};

use crate::block::{BedBlock, OwnedBedBlock, OwnedWigBlock, WigBlock, ZoomBlock};
use crate::chroms::{ChromId, Chroms};
use crate::error::{Error, Result};
use crate::header::{HEADER_LEN, Header, Kind, TotalSummary, ZoomLevel};
use crate::index::{BlockRef, Index, NodeInfo};
use crate::inflate::Inflater;
use crate::layout::{Layout, Level};

const READ_BATCH_LEN: u64 = 1 << 20;

#[derive(Debug, Clone)]
struct ReadRun {
    blocks: Range<usize>,
    file: Range<u64>,
    arena_at: usize,
}

fn grow(buf: &mut Vec<u8>, need: usize) -> Result<()> {
    if need > buf.len() {
        buf.try_reserve(need - buf.len())
            .map_err(|e| Error::invalid("BBI buffer", e.to_string()))?;
        buf.resize(need, 0);
    }
    Ok(())
}

/// A bigWig ([`TrackRead`]) or bigBed ([`IntervalRead`]) file. Opening reads
/// the header, chromosome tree and base index once and keeps them; each
/// query reads and decodes only the data blocks it overlaps, once.
///
/// Variable-size results stay in the reader until the next query: track
/// entries as decoded, bigBed features as positions into the decompressed
/// blocks, which are kept rather than copied. Writing them out is the only
/// copy.
pub struct BbiReader<S> {
    source: S,
    header: Header,
    zooms: Arc<Vec<ZoomLevel>>,
    layout: Arc<Layout>,
    chroms: Arc<Chroms>,
    index_base: u64,
    index: Arc<Vec<u8>>,
    index_nodes: Arc<Vec<NodeInfo>>,
    inflater: Option<Inflater>,
    blocks: Vec<BlockRef>,
    raw: Vec<u8>,
    raw_run: Option<ReadRun>,
    inflated: Vec<u8>,
    fill_scratch: Vec<usize>,
    /// Reused staging storage; track_ready marks a completed query.
    track_entries: Vec<TrackEntry>,
    track_ready: bool,
    /// Decompressed blocks of the last interval query, `arena[..used]`; the
    /// buffer only grows, so later queries reuse it without clearing.
    arena: Vec<u8>,
    used: usize,
    arena_run: Option<ReadRun>,
    intervals: Vec<StagedInterval>,
    interval_size: Option<IntervalSize>,
    cache: Option<Arc<BlockCache>>,
    namespace: SourceToken,
    cached_blocks: Vec<CacheLease<CachedBlock>>,
    cached_rows: Vec<CachedInterval>,
    cached_intervals: bool,
    /// Per-bin work state of bin queries.
    bins: BinAccumulator,
    /// Zoom-level R-trees, read on first use.
    zoom_indexes: Vec<Option<Arc<LevelIndex>>>,
}

/// A zoom level's R-tree, read and prepared once.
struct LevelIndex {
    base: u64,
    bytes: Vec<u8>,
    nodes: Vec<NodeInfo>,
}

enum CachedData {
    Wig(OwnedWigBlock),
    Bed(OwnedBedBlock),
    /// Decompressed zoom records.
    Zoom(Vec<u8>),
}

struct CachedBlock {
    key: BlockKey,
    data: CachedData,
}

/// A feature of the last interval query; its annotation is `arena[note]`.
#[derive(Debug, Clone)]
struct StagedInterval {
    start: u32,
    end: u32,
    note: Range<usize>,
}

struct CachedInterval {
    entry: StagedInterval,
    block: usize,
}

fn read<S: ReadAt>(source: &S, range: std::ops::Range<u64>) -> Result<Vec<u8>> {
    if range.start > range.end || range.end > source.len() {
        return Err(Error::OutOfRange {
            what: "BBI read",
            offset: range.start,
            lo: 0,
            hi: source.len(),
        });
    }
    let len = usize::try_from(range.end - range.start)
        .map_err(|_| Error::invalid("BBI read", "range exceeds usize"))?;
    let mut buf = Vec::new();
    grow(&mut buf, len)?;
    source.read_exact_at(&mut buf, range.start)?;
    Ok(buf)
}

impl<S: ReadAt> BbiReader<S> {
    pub fn open(source: S) -> Result<BbiReader<S>> {
        let header = Header::parse(&read(&source, 0..HEADER_LEN as u64)?)?;
        let zooms = ZoomLevel::parse_all(&header, &read(&source, header.zoom_headers())?)?;
        let layout = Layout::new(&header, &zooms, source.len())?;
        let tree = layout.chrom_tree();
        let chroms = Chroms::parse(&read(&source, tree.clone())?, tree.start)?;
        let index_range = layout.index(Level::Base).expect("base index");
        let index = read(&source, index_range.clone())?;
        let view = Index::new(&index, index_range.start)?;
        let mut index_nodes = Vec::new();
        view.prepare(&mut index_nodes)?;
        let levels = zooms.len();
        Ok(BbiReader {
            source,
            inflated: Vec::new(),
            fill_scratch: Vec::new(),
            header,
            zooms: Arc::new(zooms),
            layout: Arc::new(layout),
            chroms: Arc::new(chroms),
            index_base: index_range.start,
            index: Arc::new(index),
            index_nodes: Arc::new(index_nodes),
            inflater: None,
            blocks: Vec::new(),
            raw: Vec::new(),
            raw_run: None,
            track_entries: Vec::new(),
            track_ready: false,
            arena: Vec::new(),
            used: 0,
            arena_run: None,
            intervals: Vec::new(),
            interval_size: None,
            cache: None,
            namespace: SourceToken::new(),
            cached_blocks: Vec::new(),
            cached_rows: Vec::new(),
            cached_intervals: false,
            bins: BinAccumulator::new(),
            zoom_indexes: vec![None; levels],
        })
    }

    /// Attaches caller-owned decoded block storage. Source contents and
    /// indexes must remain unchanged for this reader and its clones.
    pub fn set_cache(&mut self, cache: Option<Arc<BlockCache>>) -> Result<()> {
        self.cache = cache;
        Ok(())
    }

    pub fn cache_namespace(&self) -> SourceToken {
        self.namespace
    }

    fn cache_for_query(&self) -> Option<Arc<BlockCache>> {
        let cache = self.cache.as_ref()?;
        let (bytes, entries) = cache.limits();
        // Preparing a block that cannot be retained adds a full scan and a
        // copy on every miss. Conservatively bypass that work when the
        // declared ceiling (or selected raw block size) cannot fit together
        // with record offsets. This is an admission choice, never an input
        // rejection; the ordinary path still decodes the actual bytes.
        let ceiling = if self.header.compressed() {
            self.header.block_buffer_len() as u64
        } else {
            self.blocks
                .iter()
                .map(|block| block.size)
                .max()
                .unwrap_or(0)
        };
        let estimate = ceiling
            .saturating_mul(2)
            .saturating_add(size_of::<CachedBlock>() as u64 + 256);
        (bytes as u64 >= estimate && entries != 0).then(|| Arc::clone(cache))
    }

    /// `level` is 0 for the full-resolution data and `1 + i` for zoom
    /// level `i`.
    fn block_key_at(&self, block: BlockRef, level: u64) -> BlockKey {
        BlockKey::new(
            self.namespace,
            4,
            [self.header.kind as u64, level],
            block.offset,
            block.size,
        )
    }

    fn block_key(&self, block: BlockRef) -> BlockKey {
        self.block_key_at(block, 0)
    }

    /// Lists the full-resolution blocks of a locus' spans, each once, in
    /// file order; the caller controls key retention and eviction
    /// independently of any output window.
    pub fn describe_cache(&mut self, locus: &Locus<'_>, out: &mut Vec<BlockKey>) -> Result<()> {
        out.clear();
        locus.check()?;
        let chrom = match self.chroms.id(locus.chrom) {
            Some(chrom) => chrom,
            None if self.header.kind == Kind::BigBed => return Ok(()),
            None => return Err(Error::UnknownChrom(locus.chrom.to_owned())),
        };
        if !locus.is_empty() {
            self.find_locus_blocks(None, chrom, locus.spans())?;
            out.extend(self.blocks.iter().map(|&block| self.block_key(block)));
        }
        Ok(())
    }

    /// Synchronously prepares dependencies in the external cache without
    /// allocating a final dense signal or interval result.
    pub fn prefetch_cache(&mut self, locus: &Locus<'_>) -> Result<()> {
        if self.cache.is_none() {
            return Err(Error::invalid("BBI cache", "no external cache is attached"));
        }
        let mut keys = Vec::new();
        self.describe_cache(locus, &mut keys)?;
        let Some(cache) = self.cache_for_query() else {
            return Ok(());
        };
        for i in 0..keys.len() {
            self.load_cached(i, &cache)?;
        }
        Ok(())
    }

    fn load_cached(&mut self, i: usize, cache: &BlockCache) -> Result<CacheLease<CachedBlock>> {
        let key = self.block_key(self.blocks[i]);
        let kind = self.header.kind;
        let block = cache.get_or_load(key, || {
            let bytes = self.load(i)?.to_vec();
            let (data, heap) = match kind {
                Kind::BigWig => {
                    let owned = OwnedWigBlock::new(bytes)?;
                    let heap = owned.heap_bytes();
                    (CachedData::Wig(owned), heap)
                }
                Kind::BigBed => {
                    let owned = OwnedBedBlock::new(bytes)?;
                    let heap = owned.heap_bytes();
                    (CachedData::Bed(owned), heap)
                }
            };
            Ok((CachedBlock { key, data }, heap))
        })?;
        if block.key != key {
            return Err(Error::invalid(
                "BBI cache",
                "cached block identity mismatch",
            ));
        }
        Ok(block)
    }

    /// Block `i` of `self.blocks`, a block of zoom level `level`, through
    /// the external cache.
    fn load_cached_zoom(
        &mut self,
        i: usize,
        level: usize,
        cache: &BlockCache,
    ) -> Result<CacheLease<CachedBlock>> {
        let key = self.block_key_at(self.blocks[i], 1 + level as u64);
        let block = cache.get_or_load(key, || {
            let bytes = self.load(i)?.to_vec();
            ZoomBlock::parse(&bytes)?;
            let heap = bytes.capacity();
            Ok((
                CachedBlock {
                    key,
                    data: CachedData::Zoom(bytes),
                },
                heap,
            ))
        })?;
        if block.key != key {
            return Err(Error::invalid(
                "BBI cache",
                "cached block identity mismatch",
            ));
        }
        Ok(block)
    }

    pub fn header(&self) -> &Header {
        &self.header
    }

    pub fn zooms(&self) -> &[ZoomLevel] {
        &self.zooms
    }

    pub fn layout(&self) -> &Layout {
        &self.layout
    }

    pub fn chroms(&self) -> &Chroms {
        &self.chroms
    }

    /// Statistics over every covered base, `None` when the file stores none.
    pub fn total_summary(&self) -> Result<Option<TotalSummary>> {
        match self.layout.total_summary() {
            Some(range) => Ok(Some(TotalSummary::parse(&read(&self.source, range)?)?)),
            None => Ok(None),
        }
    }

    /// The autoSql schema of a bigBed file, `None` when it stores none.
    pub fn auto_sql(&self) -> Result<Option<String>> {
        let Some(range) = self.layout.auto_sql() else {
            return Ok(None);
        };
        // The section can include arbitrary padding. Stop at the string's
        // terminator so only schema bytes need to be retained and decoded.
        let mut schema = Vec::new();
        let mut offset = range.start;
        while offset < range.end {
            let len = (range.end - offset).min(4096) as usize;
            let at = schema.len();
            schema
                .try_reserve(len)
                .map_err(|e| Error::invalid("autoSql", e.to_string()))?;
            schema.resize(at + len, 0);
            self.source.read_exact_at(&mut schema[at..], offset)?;
            if let Some(end) = memchr::memchr(0, &schema[at..]) {
                schema.truncate(at + end);
                return String::from_utf8(schema)
                    .map(Some)
                    .map_err(|e| Error::invalid("autoSql", e.utf8_error().to_string()));
            }
            offset += len as u64;
        }
        Err(Error::invalid("autoSql", "missing NUL terminator"))
    }

    fn expect(&self, kind: Kind) -> Result<()> {
        if self.header.kind != kind {
            return Err(Error::Unsupported(match kind {
                Kind::BigWig => "track access to a bigBed file",
                Kind::BigBed => "interval access to a bigWig file",
            }));
        }
        Ok(())
    }

    /// Collects the blocks overlapping `[start, end)` into `self.blocks`.
    fn find_blocks(&mut self, chrom: ChromId, start: u32, end: u32) -> Result<()> {
        self.blocks.clear();
        self.raw_run = None;
        self.arena_run = None;
        // SAFETY: open prepared these nodes for this exact index and base.
        // Both buffers and the base are private and never modified after open;
        // queries only borrow them, and their lifetimes cover the local view.
        let index = unsafe {
            Index::new(&self.index, self.index_base)?.with_prepared_nodes(&self.index_nodes)
        };
        for block in index.query(chrom, start, end) {
            self.blocks.push(block?);
        }
        Ok(())
    }

    /// Collects the blocks of the base data (`level: None`) or of a zoom
    /// level overlapping any of `spans` into `self.blocks`, each once, in
    /// file order.
    fn find_locus_blocks(
        &mut self,
        level: Option<usize>,
        chrom: ChromId,
        spans: &[Span],
    ) -> Result<()> {
        if level.is_none() && spans.len() == 1 {
            return self.find_blocks(chrom, spans[0].start, spans[0].end);
        }
        self.blocks.clear();
        self.raw_run = None;
        self.arena_run = None;
        let tree = match level {
            None => None,
            Some(level) => Some(self.zoom_index(level)?),
        };
        // SAFETY: open (or zoom_index) prepared these nodes for this exact
        // index and base; neither changes afterwards, and both outlive the
        // local view.
        let index = unsafe {
            match &tree {
                None => {
                    Index::new(&self.index, self.index_base)?.with_prepared_nodes(&self.index_nodes)
                }
                Some(tree) => Index::new(&tree.bytes, tree.base)?.with_prepared_nodes(&tree.nodes),
            }
        };
        for span in spans.iter().filter(|s| !s.is_empty()) {
            for block in index.query(chrom, span.start, span.end) {
                self.blocks.push(block?);
            }
        }
        if spans.len() > 1 {
            self.blocks.sort_unstable_by_key(|b| b.offset);
            self.blocks.dedup();
        }
        Ok(())
    }

    /// The R-tree of zoom level `level`, read on first use.
    fn zoom_index(&mut self, level: usize) -> Result<Arc<LevelIndex>> {
        if let Some(tree) = &self.zoom_indexes[level] {
            return Ok(Arc::clone(tree));
        }
        let range = self
            .layout
            .index(Level::Zoom(level))
            .ok_or_else(|| Error::invalid("BBI zoom", "zoom level without an index"))?;
        let bytes = read(&self.source, range.clone())?;
        let mut nodes = Vec::new();
        Index::new(&bytes, range.start)?.prepare(&mut nodes)?;
        let tree = Arc::new(LevelIndex {
            base: range.start,
            bytes,
            nodes,
        });
        self.zoom_indexes[level] = Some(Arc::clone(&tree));
        Ok(tree)
    }

    /// The zoom level libBigWig's `bwStats` summarizes bins of `resolution`
    /// bases from: the coarsest whose reduction is at most half a bin.
    pub fn zoom_level_for(&self, resolution: u32) -> Option<usize> {
        let half = resolution / 2;
        let mut best: Option<(u32, usize)> = None;
        for (i, z) in self.zooms.iter().enumerate() {
            if z.reduction <= half {
                let diff = half - z.reduction;
                if best.is_none_or(|(d, _)| diff < d) {
                    best = Some((diff, i));
                }
            }
        }
        best.map(|(_, i)| i)
    }

    fn plan_run(&self, i: usize) -> Result<ReadRun> {
        let first = self.blocks[i];
        let end_of = |block: BlockRef| -> Result<u64> {
            let end = block
                .offset
                .checked_add(block.size)
                .filter(|&end| end <= self.source.len());
            end.ok_or(Error::OutOfRange {
                what: "BBI data block",
                offset: block.offset,
                lo: 0,
                hi: self.source.len(),
            })
        };
        let mut end = end_of(first)?;
        let mut next = i + 1;
        while next < self.blocks.len() {
            let block = self.blocks[next];
            if block.offset != end {
                break;
            }
            let block_end = end_of(block)?;
            if block_end - first.offset > READ_BATCH_LEN {
                break;
            }
            end = block_end;
            next += 1;
        }
        usize::try_from(end - first.offset)
            .map_err(|_| Error::invalid("BBI read", "range exceeds usize"))?;
        Ok(ReadRun {
            blocks: i..next,
            file: first.offset..end,
            arena_at: 0,
        })
    }

    fn raw_range(&mut self, i: usize) -> Result<Range<usize>> {
        if !self
            .raw_run
            .as_ref()
            .is_some_and(|run| run.blocks.contains(&i))
        {
            self.raw_run = None;
            let run = self.plan_run(i)?;
            let len = (run.file.end - run.file.start) as usize;
            grow(&mut self.raw, len)?;
            self.source
                .read_exact_at(&mut self.raw[..len], run.file.start)?;
            self.raw_run = Some(run);
        }
        let run = self.raw_run.as_ref().expect("read run populated");
        let block = self.blocks[i];
        let lo = (block.offset - run.file.start) as usize;
        Ok(lo..lo + block.size as usize)
    }

    /// Reads block `i` of `self.blocks` and returns its decompressed bytes.
    fn load(&mut self, i: usize) -> Result<&[u8]> {
        let raw = self.raw_range(i)?;
        if self.header.compressed() {
            let n = self
                .inflater
                .get_or_insert_with(Inflater::new)
                .inflate_growing(
                    &self.raw[raw],
                    &mut self.inflated,
                    0,
                    self.header.block_buffer_len(),
                )?;
            Ok(&self.inflated[..n])
        } else {
            Ok(&self.raw[raw])
        }
    }

    fn track_chrom(&self, region: &Region<'_>) -> Result<ChromId> {
        region.check()?;
        self.chroms
            .id(region.chrom)
            .ok_or_else(|| Error::UnknownChrom(region.chrom.to_owned()))
    }

    /// Visits each track interval overlapping `region`, clipped to it.
    fn each_track_entry(
        &mut self,
        region: Region<'_>,
        mut f: impl FnMut(TrackEntry),
    ) -> Result<()> {
        self.expect(Kind::BigWig)?;
        let chrom = self.track_chrom(&region)?;
        if region.is_empty() {
            return Ok(());
        }
        self.find_blocks(chrom, region.start, region.end)?;
        if let Some(cache) = self.cache_for_query() {
            for i in 0..self.blocks.len() {
                let data = self.load_cached(i, &cache)?;
                let CachedData::Wig(block) = &data.data else {
                    return Err(Error::invalid("BBI cache", "expected a bigWig block"));
                };
                for (s, e, value) in block.view().overlapping(chrom, region.start, region.end) {
                    f(TrackEntry {
                        start: s.max(region.start),
                        end: e.min(region.end),
                        value,
                    });
                }
            }
            return Ok(());
        }
        for i in 0..self.blocks.len() {
            let block = WigBlock::parse(self.load(i)?)?;
            for (s, e, value) in block.overlapping(chrom, region.start, region.end) {
                f(TrackEntry {
                    start: s.max(region.start),
                    end: e.min(region.end),
                    value,
                });
            }
        }
        Ok(())
    }

    /// Reads block `i` of `self.blocks`, decompressed, onto the end of the
    /// arena and returns where it landed.
    fn append_block(&mut self, i: usize) -> Result<Range<usize>> {
        if !self.header.compressed() {
            if !self
                .arena_run
                .as_ref()
                .is_some_and(|run| run.blocks.contains(&i))
            {
                self.arena_run = None;
                let mut run = self.plan_run(i)?;
                let len = (run.file.end - run.file.start) as usize;
                let end = self
                    .used
                    .checked_add(len)
                    .ok_or_else(|| Error::invalid("BBI arena", "length overflows usize"))?;
                grow(&mut self.arena, end)?;
                self.source
                    .read_exact_at(&mut self.arena[self.used..end], run.file.start)?;
                run.arena_at = self.used;
                self.used = end;
                self.arena_run = Some(run);
            }
            let run = self.arena_run.as_ref().expect("arena run populated");
            let block = self.blocks[i];
            let at = run.arena_at + (block.offset - run.file.start) as usize;
            return Ok(at..at + block.size as usize);
        }
        let raw = self.raw_range(i)?;
        let at = self.used;
        let len = self
            .inflater
            .get_or_insert_with(Inflater::new)
            .inflate_growing(
                &self.raw[raw],
                &mut self.arena,
                at,
                self.header.block_buffer_len(),
            )?;
        self.used = at + len;
        Ok(at..at + len)
    }

    /// Decodes the features overlapping `region` into the staged state.
    fn stage_intervals(&mut self, region: Region<'_>) -> Result<IntervalSize> {
        self.expect(Kind::BigBed)?;
        region.check()?;
        let mut size = IntervalSize::default();
        let Some(chrom) = self.chroms.id(region.chrom) else {
            return Ok(size);
        };
        if region.is_empty() {
            return Ok(size);
        }
        self.find_blocks(chrom, region.start, region.end)?;
        if let Some(cache) = self.cache_for_query() {
            self.cached_intervals = true;
            for i in 0..self.blocks.len() {
                let data = self.load_cached(i, &cache)?;
                let CachedData::Bed(block) = &data.data else {
                    return Err(Error::invalid("BBI cache", "expected a bigBed block"));
                };
                let base = block.bytes().as_ptr() as usize;
                let block_index = self.cached_blocks.len();
                for entry in block.view().overlapping(chrom, region.start, region.end) {
                    let at = entry.rest.as_ptr() as usize - base;
                    self.cached_rows.push(CachedInterval {
                        entry: StagedInterval {
                            start: entry.start,
                            end: entry.end,
                            note: at..at + entry.rest.len(),
                        },
                        block: block_index,
                    });
                    size.entries += 1;
                    size.annotation_bytes = size
                        .annotation_bytes
                        .checked_add(entry.rest.len())
                        .ok_or_else(|| {
                            Error::invalid("BBI query", "annotation size overflows usize")
                        })?;
                }
                self.cached_blocks.push(data);
            }
        } else {
            for i in 0..self.blocks.len() {
                let range = self.append_block(i)?;
                let (arena, intervals) = (&self.arena, &mut self.intervals);
                let base = arena.as_ptr() as usize;
                for entry in BedBlock::parse(&arena[range])?.entries() {
                    let entry = entry?;
                    if entry.overlaps(chrom, region.start, region.end) {
                        let at = entry.rest.as_ptr() as usize - base;
                        intervals.push(StagedInterval {
                            start: entry.start,
                            end: entry.end,
                            note: at..at + entry.rest.len(),
                        });
                        size.entries += 1;
                        size.annotation_bytes = size
                            .annotation_bytes
                            .checked_add(entry.rest.len())
                            .ok_or_else(|| {
                                Error::invalid("BBI query", "annotation size overflows usize")
                            })?;
                    }
                }
            }
        }
        if size.entries > u32::MAX as usize {
            return Err(Error::invalid("BBI query", "more than 2^32-1 features"));
        }
        // The high-level interface promises start order even when blocks or
        // their records are unordered. Stable sorting retains stored order
        // for equal starts and keeps each annotation attached to its record.
        if !self.intervals.is_sorted_by_key(|entry| entry.start) {
            self.intervals.sort_by_key(|entry| entry.start);
        }
        if !self.cached_rows.is_sorted_by_key(|row| row.entry.start) {
            self.cached_rows.sort_by_key(|row| row.entry.start);
        }
        Ok(size)
    }
}

impl<S: ReadAt> BbiReader<S> {
    /// Feeds the full-resolution data of `locus` to `acc`, every block
    /// decoded once and its records visited in file order.
    fn accumulate_exact(
        &mut self,
        chrom: ChromId,
        locus: &Locus<'_>,
        acc: &mut BinAccumulator,
    ) -> Result<()> {
        self.find_locus_blocks(None, chrom, locus.spans())?;
        let pieces: Vec<(u32, Region<'_>)> = locus.pieces().collect();
        let feed = |block: &WigBlock<'_>, acc: &mut BinAccumulator| {
            for &(at, r) in &pieces {
                for (s, e, value) in block.overlapping(chrom, r.start, r.end) {
                    let (s, e) = (s.max(r.start), e.min(r.end));
                    if s < e {
                        acc.add(at + (s - r.start), at + (e - r.start), value);
                    }
                }
            }
        };
        if let Some(cache) = self.cache_for_query() {
            for i in 0..self.blocks.len() {
                let data = self.load_cached(i, &cache)?;
                let CachedData::Wig(block) = &data.data else {
                    return Err(Error::invalid("BBI cache", "expected a bigWig block"));
                };
                feed(&block.view(), acc);
            }
            return Ok(());
        }
        for i in 0..self.blocks.len() {
            feed(&WigBlock::parse(self.load(i)?)?, acc);
        }
        Ok(())
    }

    /// Feeds the records of zoom level `level` overlapping `locus` to
    /// `acc`, every block decoded once.
    fn accumulate_zoom(
        &mut self,
        level: usize,
        chrom: ChromId,
        locus: &Locus<'_>,
        acc: &mut BinAccumulator,
    ) -> Result<()> {
        self.find_locus_blocks(Some(level), chrom, locus.spans())?;
        let pieces: Vec<(u32, Region<'_>)> = locus.pieces().collect();
        let feed = |block: &ZoomBlock<'_>, acc: &mut BinAccumulator| {
            for &(at, r) in &pieces {
                for z in block.overlapping(chrom, r.start, r.end) {
                    let (s, e) = (z.start.max(r.start), z.end.min(r.end));
                    if s < e {
                        acc.add_zoom(ZoomPiece {
                            start: at + (s - r.start),
                            end: at + (e - r.start),
                            record_len: z.end - z.start,
                            bases: z.valid,
                            min: z.min,
                            max: z.max,
                            sum: z.sum,
                            squares: z.sum_squares,
                        });
                    }
                }
            }
        };
        if let Some(cache) = self.cache_for_query() {
            for i in 0..self.blocks.len() {
                let data = self.load_cached_zoom(i, level, &cache)?;
                let CachedData::Zoom(bytes) = &data.data else {
                    return Err(Error::invalid("BBI cache", "expected a zoom block"));
                };
                feed(&ZoomBlock::parse(bytes)?, acc);
            }
            return Ok(());
        }
        for i in 0..self.blocks.len() {
            feed(&ZoomBlock::parse(self.load(i)?)?, acc);
        }
        Ok(())
    }
}

impl<S: ReadAt> TrackRead for BbiReader<S> {
    fn chrom_len(&self, chrom: &str) -> Result<u32> {
        let id = self
            .chroms
            .id(chrom)
            .ok_or_else(|| Error::UnknownChrom(chrom.to_owned()))?;
        Ok(self.chroms.get(id).expect("id from the table").len)
    }

    fn values(&mut self, locus: &Locus<'_>, out: &mut [f32]) -> Result<()> {
        self.expect(Kind::BigWig)?;
        locus.check()?;
        let chrom = self.track_chrom(&locus.hull())?;
        check_output("values", out.len(), locus.len() as usize)?;
        out.fill(f32::NAN);
        if locus.is_empty() {
            return Ok(());
        }
        self.find_locus_blocks(None, chrom, locus.spans())?;
        let mut scratch = std::mem::take(&mut self.fill_scratch);
        let pieces: Vec<(u32, Region<'_>)> = locus.pieces().collect();
        let result = (|| {
            let mut fill = |block: &WigBlock<'_>, scratch: &mut Vec<usize>| -> Result<()> {
                for &(at, r) in &pieces {
                    let at = at as usize;
                    block.fill_reusing(
                        chrom,
                        r.start,
                        &mut out[at..at + r.len() as usize],
                        scratch,
                    )?;
                }
                Ok(())
            };
            if let Some(cache) = self.cache_for_query() {
                for i in 0..self.blocks.len() {
                    let data = self.load_cached(i, &cache)?;
                    let CachedData::Wig(block) = &data.data else {
                        return Err(Error::invalid("BBI cache", "expected a bigWig block"));
                    };
                    fill(&block.view(), &mut scratch)?;
                }
                return Ok(());
            }
            for i in 0..self.blocks.len() {
                fill(&WigBlock::parse(self.load(i)?)?, &mut scratch)?;
            }
            Ok(())
        })();
        self.fill_scratch = scratch;
        result?;
        orient_values(locus.strand, out);
        Ok(())
    }

    fn bins(&mut self, locus: &Locus<'_>, spec: &BinSpec<'_>, out: &mut [f32]) -> Result<()> {
        self.expect(Kind::BigWig)?;
        locus.check()?;
        let layout = BinLayout::new(locus, spec)?;
        let chrom = self.track_chrom(&locus.hull())?;
        check_output("bins", out.len(), spec.output_len(locus.len()))?;
        if out.is_empty() {
            return Ok(());
        }
        let level = match spec.source {
            BinSource::Exact => None,
            BinSource::Zoom => self.zoom_level_for(spec.resolution),
        };
        let mut acc = std::mem::take(&mut self.bins);
        acc.reset(layout, spec.stats, level.is_some());
        let result = match level {
            None => self.accumulate_exact(chrom, locus, &mut acc),
            Some(level) => self.accumulate_zoom(level, chrom, locus, &mut acc),
        };
        if result.is_ok() {
            acc.finish(spec.stats, out);
        }
        self.bins = acc;
        result
    }

    fn query_entries(&mut self, region: Region<'_>) -> Result<usize> {
        self.track_ready = false;
        let mut staged = std::mem::take(&mut self.track_entries);
        staged.clear();
        let result = self.each_track_entry(region, |entry| staged.push(entry));
        let n = staged.len();
        self.track_entries = staged;
        result?;
        if !self.track_entries.is_sorted_by_key(|entry| entry.start) {
            self.track_entries.sort_by_key(|entry| entry.start);
        }
        self.track_ready = true;
        Ok(n)
    }

    fn write_entries(&self, out: &mut [TrackEntry]) -> Result<()> {
        if !self.track_ready {
            return Err(Error::NotQueried);
        }
        let staged = &self.track_entries;
        check_output("track entries", out.len(), staged.len())?;
        out.copy_from_slice(staged);
        Ok(())
    }
}

impl<S: ReadAt> IntervalRead for BbiReader<S> {
    fn query(&mut self, region: Region<'_>) -> Result<IntervalSize> {
        self.interval_size = None;
        self.intervals.clear();
        self.cached_blocks.clear();
        self.cached_rows.clear();
        self.cached_intervals = false;
        self.used = 0;
        let size = self.stage_intervals(region)?;
        self.interval_size = Some(size);
        Ok(size)
    }

    fn write(&self, out: IntervalsOut<'_>) -> Result<()> {
        let size = self.interval_size.ok_or(Error::NotQueried)?;
        out.check(size)?;
        out.annotations.offsets[0] = 0;
        let mut byte = 0;
        if self.cached_intervals {
            for (row, cached) in self.cached_rows.iter().enumerate() {
                let staged = &cached.entry;
                let CachedData::Bed(block) = &self.cached_blocks[cached.block].data else {
                    return Err(Error::invalid("BBI cache", "expected a bigBed block"));
                };
                let note = &block.bytes()[staged.note.clone()];
                let end = byte + note.len();
                out.annotations.data[byte..end].copy_from_slice(note);
                out.entries[row] = IntervalEntry {
                    start: staged.start,
                    end: staged.end,
                    annotation: row as u32,
                };
                out.annotations.offsets[row + 1] = end as i64;
                byte = end;
            }
            debug_assert_eq!(byte, size.annotation_bytes);
            return Ok(());
        }
        for (row, staged) in self.intervals.iter().enumerate() {
            // SAFETY: A successful query records exactly size.entries ranges
            // inside the retained arena, and their lengths sum to annotation_bytes.
            // check verified all destination lengths before writing. Staged data
            // stays immutable until the next query; separate output borrows cannot
            // alias it. The query also bounds the annotation row index to u32.
            unsafe {
                let note = self.arena.get_unchecked(staged.note.clone());
                let end = byte + note.len();
                out.annotations
                    .data
                    .get_unchecked_mut(byte..end)
                    .copy_from_slice(note);
                *out.entries.get_unchecked_mut(row) = IntervalEntry {
                    start: staged.start,
                    end: staged.end,
                    annotation: row as u32,
                };
                *out.annotations.offsets.get_unchecked_mut(row + 1) = end as i64;
                byte = end;
            }
        }
        debug_assert_eq!(byte, size.annotation_bytes);
        Ok(())
    }
}

impl BbiReader<omicsweave_core::FileSource> {
    /// Shares immutable metadata and the external cache while keeping all
    /// staging and decoder scratch independent for the new execution context.
    pub fn try_clone(&self) -> Result<Self> {
        Ok(Self {
            source: self.source.try_clone()?,
            header: self.header,
            zooms: Arc::clone(&self.zooms),
            layout: Arc::clone(&self.layout),
            chroms: Arc::clone(&self.chroms),
            index_base: self.index_base,
            index: Arc::clone(&self.index),
            index_nodes: Arc::clone(&self.index_nodes),
            inflater: None,
            blocks: Vec::new(),
            raw: Vec::new(),
            raw_run: None,
            inflated: Vec::new(),
            fill_scratch: Vec::new(),
            track_entries: Vec::new(),
            track_ready: false,
            arena: Vec::new(),
            used: 0,
            arena_run: None,
            intervals: Vec::new(),
            interval_size: None,
            cache: self.cache.as_ref().map(Arc::clone),
            namespace: self.namespace,
            cached_blocks: Vec::new(),
            cached_rows: Vec::new(),
            cached_intervals: false,
            bins: BinAccumulator::new(),
            zoom_indexes: self.zoom_indexes.clone(),
        })
    }
}
