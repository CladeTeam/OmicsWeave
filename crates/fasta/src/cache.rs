//! Optional externally owned block data. Ordinary readers bypass this module.

use std::ops::Range;
use std::sync::Arc;

use omicsweave_core::cache::{BlockCache, BlockKey, SourceToken};
use omicsweave_core::{FileSource, Locus, ReadAt, Region, check_output};

use crate::bgzf::{Inflater, MAX_BLOCK_LEN};
use crate::error::{Error, Result};
use crate::fai::TextDecoder;
use crate::gzi::BlockRef;
use crate::reader::{BgzfRead, FastaReader, TwoBitReader, entry, grow, push_block, read_len};

const BGZF_KIND: u32 = 1;
const TWOBIT_KIND: u32 = 2;

struct BgzfText {
    key: BlockKey,
    bytes: Vec<u8>,
}

struct TwoBitTile {
    key: BlockKey,
    bytes: Vec<u8>,
}

fn bgzf_key(namespace: SourceToken, block: &BlockRef) -> BlockKey {
    BlockKey::new(namespace, BGZF_KIND, [0, 0], block.compressed.start, 0)
}

fn twobit_key(namespace: SourceToken, record: usize, tile: u32, bases: u32) -> BlockKey {
    BlockKey::new(
        namespace,
        TWOBIT_KIND,
        [record as u64, 0],
        tile as u64,
        bases as u64,
    )
}

impl<S: ReadAt> FastaReader<S> {
    /// Describes the external cache entries this locus uses, in file
    /// order and each once, clearing `out` first. Keys refer to GZI-indexed
    /// member offsets; the implicit first entry may include leading empty
    /// members, validated on load.
    pub fn describe_cache(&self, locus: &Locus<'_>, out: &mut Vec<BlockKey>) -> Result<()> {
        out.clear();
        if self.gzi.is_none() {
            return Err(Error::Unsupported(
                "plain FASTA does not use a decoded-block cache",
            ));
        }
        locus.check()?;
        locus.check_within(entry(&self.fai, locus.chrom)?.len)?;
        for (_, region) in locus.pieces() {
            self.describe_region(region, out)?;
        }
        out.dedup();
        Ok(())
    }

    fn describe_region(&self, region: Region<'_>, out: &mut Vec<BlockKey>) -> Result<()> {
        let gzi = self.gzi.as_ref().ok_or(Error::Unsupported(
            "plain FASTA does not use a decoded-block cache",
        ))?;
        let entry = entry(&self.fai, region.chrom)?;
        region.check_within(entry.len)?;
        if !region.is_empty() {
            let range = entry.text_range(region.start, region.end)?;
            for block in gzi.blocks(range, self.source.len()) {
                out.push(bgzf_key(self.namespace, &block));
            }
        }
        Ok(())
    }

    pub(crate) fn sequence_cached(
        &mut self,
        region: Region<'_>,
        out: &mut [u8],
        cache: &BlockCache,
    ) -> Result<()> {
        self.cached_blocks(region, Some(out), cache)
    }

    /// Loads the locus' decoded BGZF blocks into the attached external
    /// cache without allocating or producing a sequence output. Admission
    /// still obeys the external cache's budget and policy.
    pub fn prefetch_cache(&mut self, locus: &Locus<'_>) -> Result<()> {
        let cache = self
            .cache
            .as_ref()
            .map(Arc::clone)
            .ok_or_else(|| Error::invalid("FASTA cache", "no external cache is attached"))?;
        locus.check()?;
        for (_, region) in locus.pieces() {
            self.cached_blocks(region, None, &cache)?;
        }
        Ok(())
    }

    fn cached_blocks(
        &mut self,
        region: Region<'_>,
        mut out: Option<&mut [u8]>,
        cache: &BlockCache,
    ) -> Result<()> {
        let entry = entry(&self.fai, region.chrom)?;
        region.check_within(entry.len)?;
        if let Some(out) = &out {
            check_output("sequence", out.len(), region.len() as usize)?;
        }
        if region.is_empty() {
            return Ok(());
        }
        let range = entry.text_range(region.start, region.end)?;
        let gzi = self.gzi.as_deref().ok_or(Error::Unsupported(
            "plain FASTA does not use a decoded-block cache",
        ))?;
        let mut read = BgzfRead {
            source: &self.source,
            gzi,
            inflater: self.inflater.get_or_insert_with(Inflater::new),
            raw: &mut self.raw,
            text: &mut self.text,
        };
        let indices = gzi.indices(range.clone());
        let file_len = self.source.len();
        let mut raw_run: Option<Range<u64>> = None;
        let mut decoder = out
            .as_ref()
            .map(|_| TextDecoder::new(entry, region.start, region.end))
            .transpose()?;
        let mut text_cursor = read.block(indices.start, file_len).text_start;
        for i in indices.clone() {
            let block = read.block(i, file_len);
            if block.text_start != text_cursor {
                return Err(Error::invalid(
                    "gzi",
                    "text offsets disagree with decompressed blocks",
                ));
            }
            let key = bgzf_key(self.namespace, &block);
            let data = cache.get_or_load(key, || {
                let count = read_len(&block.compressed, file_len)?;
                if !raw_run.as_ref().is_some_and(|run| {
                    run.start <= block.compressed.start && block.compressed.end <= run.end
                }) {
                    // A miss reads one bounded adjacent run, so following
                    // misses reuse compressed bytes instead of issuing one
                    // positional read per block. Hits perform no reads.
                    let start = block.compressed.start;
                    let mut end = block.compressed.end;
                    for j in i + 1..indices.end {
                        let next = read.block(j, file_len);
                        if next.compressed.start != end
                            || next.compressed.end < next.compressed.start
                            || next.compressed.end.saturating_sub(start) > 1 << 20
                        {
                            break;
                        }
                        end = next.compressed.end;
                    }
                    let run = start..end;
                    let len = read_len(&run, file_len)?;
                    grow(read.raw, len)?;
                    read.source.read_exact_at(&mut read.raw[..len], start)?;
                    raw_run = Some(run);
                }
                let run = raw_run.as_ref().expect("miss populated compressed run");
                let at = (block.compressed.start - run.start) as usize;
                grow(read.text, MAX_BLOCK_LEN)?;
                let (used, mut n) = read
                    .inflater
                    .inflate_block(&read.raw[at..at + count], &mut read.text[..MAX_BLOCK_LEN])?;
                if n == 0 && block.text_len != Some(0) {
                    n = read.after_empty_members(block.compressed.start + used as u64, file_len)?;
                }
                if block.text_len.is_some_and(|expected| expected != n as u64) {
                    return Err(Error::invalid(
                        "gzi",
                        "indexed block lengths differ from BGZF data",
                    ));
                }
                let bytes = read.text[..n].to_vec();
                let charge = bytes.capacity();
                Ok((BgzfText { key, bytes }, charge))
            })?;
            if data.key != key || block.text_len.is_some_and(|n| n != data.bytes.len() as u64) {
                return Err(Error::invalid(
                    "FASTA cache",
                    "entry does not match indexed block",
                ));
            }
            if let (Some(decoder), Some(out)) = (decoder.as_mut(), out.as_deref_mut()) {
                push_block(decoder, &data.bytes, block.text_start, &range, out)?;
            }
            text_cursor = block
                .text_start
                .checked_add(data.bytes.len() as u64)
                .ok_or_else(|| Error::invalid("gzi", "text offset overflows u64"))?;
        }
        if text_cursor < range.end {
            return Err(Error::Truncated {
                what: "FASTA text",
                need: range.end,
                have: text_cursor,
            });
        }
        if let Some(decoder) = decoder {
            decoder.finish(region.len() as usize)?;
        }
        Ok(())
    }
}

impl<S: ReadAt> TwoBitReader<S> {
    /// Describes the fixed decoded tiles this locus uses, in order and each
    /// once, clearing `out`.
    pub fn describe_cache(&self, locus: &Locus<'_>, out: &mut Vec<BlockKey>) -> Result<()> {
        out.clear();
        locus.check()?;
        for span in locus.spans() {
            self.describe_region(Region::new(locus.chrom, span.start, span.end), out)?;
        }
        out.dedup();
        Ok(())
    }

    fn describe_region(&self, region: Region<'_>, out: &mut Vec<BlockKey>) -> Result<()> {
        let record = self
            .directory
            .index(region.chrom)
            .ok_or_else(|| Error::UnknownChrom(region.chrom.to_owned()))?;
        region.check_within(self.records[record].len())?;
        if !region.is_empty() {
            for tile in region.start / self.tile_bases..=(region.end - 1) / self.tile_bases {
                out.push(twobit_key(self.namespace, record, tile, self.tile_bases));
            }
        }
        Ok(())
    }

    pub(crate) fn sequence_cached(
        &mut self,
        region: Region<'_>,
        out: &mut [u8],
        cache: &BlockCache,
    ) -> Result<()> {
        self.cached_tiles(region, Some(out), cache)
    }

    /// Loads fixed decoded tiles into the attached external cache without
    /// constructing any query result. Tiles larger than the cache budget
    /// are skipped; later queries use the normal direct path.
    pub fn prefetch_cache(&mut self, locus: &Locus<'_>) -> Result<()> {
        let cache = self
            .cache
            .as_ref()
            .map(Arc::clone)
            .ok_or_else(|| Error::invalid("2bit cache", "no external cache is attached"))?;
        locus.check()?;
        for (_, region) in locus.pieces() {
            self.cached_tiles(region, None, &cache)?;
        }
        Ok(())
    }

    fn cached_tiles(
        &mut self,
        region: Region<'_>,
        mut out: Option<&mut [u8]>,
        cache: &BlockCache,
    ) -> Result<()> {
        let record_index = self
            .directory
            .index(region.chrom)
            .ok_or_else(|| Error::UnknownChrom(region.chrom.to_owned()))?;
        let record = &self.records[record_index];
        region.check_within(record.len())?;
        if let Some(out) = &out {
            check_output("sequence", out.len(), region.len() as usize)?;
        }
        if region.is_empty() {
            return Ok(());
        }
        let (bytes_budget, entry_budget) = cache.limits();
        let tile_len = self.tile_bases.min(record.len()) as usize;
        if entry_budget == 0 || bytes_budget <= tile_len.saturating_add(size_of::<TwoBitTile>()) {
            let Some(out) = out else {
                return Ok(());
            };
            // A cache too small for any full tile must not turn a tiny query
            // into a full-tile read/decode on every admission rejection.
            let range = record.packed_range(region.start, region.end)?;
            let n = read_len(&range, self.source.len())?;
            if n == 0 {
                out.fill(b'N');
                return Ok(());
            }
            grow(&mut self.packed, n)?;
            self.source
                .read_exact_at(&mut self.packed[..n], range.start)?;
            return record.decode(region.start, &self.packed[..n], out, self.mask);
        }
        for tile in region.start / self.tile_bases..=(region.end - 1) / self.tile_bases {
            let start = tile * self.tile_bases;
            let end = record.len().min(start.saturating_add(self.tile_bases));
            let key = twobit_key(self.namespace, record_index, tile, self.tile_bases);
            let data = cache.get_or_load(key, || {
                let range = record.packed_range(start, end)?;
                let n = read_len(&range, self.source.len())?;
                grow(&mut self.packed, n)?;
                if n != 0 {
                    self.source
                        .read_exact_at(&mut self.packed[..n], range.start)?;
                }
                let mut bytes = vec![0; (end - start) as usize];
                record.decode(start, &self.packed[..n], &mut bytes, false)?;
                let charge = bytes.capacity();
                Ok((TwoBitTile { key, bytes }, charge))
            })?;
            if data.key != key || data.bytes.len() != (end - start) as usize {
                return Err(Error::invalid(
                    "2bit cache",
                    "entry does not match decoded tile",
                ));
            }
            let lo = start.max(region.start);
            let hi = end.min(region.end);
            if let Some(out) = out.as_deref_mut() {
                out[(lo - region.start) as usize..(hi - region.start) as usize]
                    .copy_from_slice(&data.bytes[(lo - start) as usize..(hi - start) as usize]);
            }
        }
        if let Some(out) = out.filter(|_| self.mask) {
            record.mask_output(region.start, out, true);
        }
        Ok(())
    }
}

impl FastaReader<FileSource> {
    /// Shares the immutable source session, indexes and external cache;
    /// the new execution context has independent decompression scratch.
    pub fn try_clone(&self) -> Result<Self> {
        Ok(Self {
            source: self.source.try_clone()?,
            fai: Arc::clone(&self.fai),
            gzi: self.gzi.as_ref().map(Arc::clone),
            inflater: None,
            raw: Vec::new(),
            text: Vec::new(),
            hull: Vec::new(),
            cache: self.cache.as_ref().map(Arc::clone),
            namespace: self.namespace,
        })
    }
}

impl TwoBitReader<FileSource> {
    /// Shares the immutable source session, metadata and external cache;
    /// packed input scratch is independent in the new execution context.
    pub fn try_clone(&self) -> Result<Self> {
        Ok(Self {
            source: self.source.try_clone()?,
            directory: Arc::clone(&self.directory),
            records: Arc::clone(&self.records),
            mask: self.mask,
            packed: Vec::new(),
            cache: self.cache.as_ref().map(Arc::clone),
            namespace: self.namespace,
            tile_bases: self.tile_bases,
        })
    }
}
