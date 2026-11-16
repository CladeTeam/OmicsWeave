//! Standard [`SequenceRead`] access, composing the crate's steps over a
//! caller-provided byte source.

use std::iter::FusedIterator;
use std::ops::Range;
use std::sync::Arc;

use omicsweave_core::cache::{BlockCache, SourceToken};
use omicsweave_core::orient::orient_bases;
use omicsweave_core::{Locus, ReadAt, Region, SequenceRead, check_output};

use crate::bgzf::{Inflater, MAX_BLOCK_LEN};
use crate::error::{Error, Result};
use crate::fai::{Fai, FaiEntry, TextDecoder};
use crate::gzi::{BlockRef, Gzi};
use crate::twobit::{Directory, HEADER_LEN, Parsed, Record, TwoBitHeader};

const READ_BATCH_LEN: u64 = 1 << 20;

pub(crate) fn read_len(range: &Range<u64>, file_len: u64) -> Result<usize> {
    if range.start > range.end || range.end > file_len {
        return Err(Error::OutOfRange {
            what: "sequence bytes",
            offset: range.start,
            lo: 0,
            hi: file_len,
        });
    }
    usize::try_from(range.end - range.start)
        .map_err(|_| Error::invalid("sequence bytes", "range exceeds usize"))
}

pub(crate) fn grow(buf: &mut Vec<u8>, need: usize) -> Result<()> {
    if need > buf.len() {
        buf.try_reserve(need - buf.len())
            .map_err(|e| Error::invalid("sequence buffer", e.to_string()))?;
        buf.resize(need, 0);
    }
    Ok(())
}

pub(crate) fn entry<'f>(fai: &'f Fai, chrom: &str) -> Result<&'f FaiEntry> {
    fai.get(chrom)
        .ok_or_else(|| Error::UnknownChrom(chrom.to_owned()))
}

/// A FASTA file with its `.fai`: plain text, or bgzip when its `.gzi` is
/// given. Bases come back as stored.
pub struct FastaReader<S> {
    pub(crate) source: S,
    pub(crate) fai: Arc<Fai>,
    pub(crate) gzi: Option<Arc<Gzi>>,
    pub(crate) inflater: Option<Inflater>,
    pub(crate) raw: Vec<u8>,
    pub(crate) text: Vec<u8>,
    /// Bases of a run of spans that share BGZF blocks, decoded once.
    pub(crate) hull: Vec<u8>,
    pub(crate) cache: Option<Arc<BlockCache>>,
    pub(crate) namespace: SourceToken,
}

impl<S: ReadAt> FastaReader<S> {
    /// `gzi` is required for bgzip files and must be `None` for plain text.
    pub fn new(source: S, fai: Fai, gzi: Option<Gzi>) -> FastaReader<S> {
        FastaReader {
            source,
            fai: Arc::new(fai),
            gzi: gzi.map(Arc::new),
            inflater: None,
            raw: Vec::new(),
            text: Vec::new(),
            hull: Vec::new(),
            cache: None,
            namespace: SourceToken::new(),
        }
    }

    pub fn fai(&self) -> &Fai {
        &self.fai
    }

    pub fn gzi(&self) -> Option<&Gzi> {
        self.gzi.as_deref()
    }

    /// Attaches an externally managed decoded-block cache. The source and
    /// indexes must remain unchanged while this reader or its clones exist.
    /// Plain FASTA is read directly and does not accept a block cache.
    pub fn set_cache(&mut self, cache: Option<Arc<BlockCache>>) -> Result<()> {
        if cache.is_some() && self.gzi.is_none() {
            return Err(Error::Unsupported(
                "plain FASTA does not use a decoded-block cache",
            ));
        }
        self.cache = cache;
        Ok(())
    }

    /// The immutable-source session used to isolate external cache entries.
    pub fn cache_namespace(&self) -> SourceToken {
        self.namespace
    }
}

impl<S: ReadAt> SequenceRead for FastaReader<S> {
    fn chrom_len(&self, chrom: &str) -> Result<u32> {
        Ok(entry(&self.fai, chrom)?.len)
    }

    fn sequence(&mut self, locus: &Locus<'_>, out: &mut [u8]) -> Result<()> {
        locus.check()?;
        let entry_len = entry(&self.fai, locus.chrom)?.len;
        locus.check_within(entry_len)?;
        check_output("sequence", out.len(), locus.len() as usize)?;
        if locus.spans().len() == 1 {
            self.read_region(locus.hull(), out)?;
        } else {
            self.read_spans(locus, out)?;
        }
        orient_bases(locus.strand, out);
        Ok(())
    }
}

impl<S: ReadAt> FastaReader<S> {
    /// The forward concatenation of a multi-span locus. In bgzip files,
    /// consecutive spans that share BGZF blocks are decoded in bounded
    /// runs. A run boundary may revisit one block, keeping a sparse locus
    /// from materializing a chromosome-sized hull in working memory.
    fn read_spans(&mut self, locus: &Locus<'_>, out: &mut [u8]) -> Result<()> {
        if let Some(cache) = self.cache.as_ref().map(Arc::clone) {
            return self.sequence_cached_locus(locus, out, &cache);
        }
        let Some(gzi) = self.gzi.as_ref().map(Arc::clone) else {
            for (at, region) in locus.pieces() {
                let at = at as usize;
                self.read_region(region, &mut out[at..at + region.len() as usize])?;
            }
            return Ok(());
        };
        let fai = Arc::clone(&self.fai);
        let entry = entry(&fai, locus.chrom)?;
        let spans = locus.spans();
        let mut i = 0;
        let mut output_at = 0usize;
        let mut pending = None;
        while i < spans.len() {
            if spans[i].is_empty() {
                i += 1;
                continue;
            }
            let first = spans[i];
            let blocks = match pending.take() {
                Some(blocks) => blocks,
                None => gzi.indices(entry.text_range(first.start, first.end)?),
            };
            let mut block_end = blocks.end;
            let mut hull_end = first.end;
            let mut count = 1;
            let mut j = i + 1;
            while j < spans.len() {
                let span = spans[j];
                if span.is_empty() {
                    j += 1;
                    continue;
                }
                let blocks = gzi.indices(entry.text_range(span.start, span.end)?);
                if blocks.start >= block_end || (span.end - first.start) as u64 > READ_BATCH_LEN {
                    pending = Some(blocks);
                    break;
                }
                hull_end = span.end;
                block_end = blocks.end;
                count += 1;
                j += 1;
            }
            if count == 1 {
                let len = first.len() as usize;
                self.read_region(
                    Region::new(locus.chrom, first.start, first.end),
                    &mut out[output_at..output_at + len],
                )?;
                output_at += len;
            } else {
                let hull = Region::new(locus.chrom, first.start, hull_end);
                let mut bases = std::mem::take(&mut self.hull);
                let len = hull.len() as usize;
                let result = (|| {
                    if len > bases.len() {
                        bases
                            .try_reserve_exact(len - bases.len())
                            .map_err(|e| Error::invalid("sequence buffer", e.to_string()))?;
                        bases.resize(len, 0);
                    }
                    self.read_region(hull, &mut bases[..len])
                })();
                if result.is_ok() {
                    for span in &spans[i..j] {
                        if span.is_empty() {
                            continue;
                        }
                        let lo = (span.start - hull.start) as usize;
                        let len = span.len() as usize;
                        out[output_at..output_at + len].copy_from_slice(&bases[lo..lo + len]);
                        output_at += len;
                    }
                }
                self.hull = bases;
                result?;
            }
            i = j;
        }
        Ok(())
    }
}

impl<S: ReadAt> FastaReader<S> {
    /// The bases of one region, as stored.
    fn read_region(&mut self, region: Region<'_>, out: &mut [u8]) -> Result<()> {
        if let Some(cache) = self.cache.as_ref().map(Arc::clone) {
            return self.sequence_cached(region, out, &cache);
        }
        let FastaReader {
            source,
            fai,
            gzi,
            inflater,
            raw,
            text,
            ..
        } = self;
        let entry = entry(fai, region.chrom)?;
        region.check_within(entry.len)?;
        check_output("sequence", out.len(), region.len() as usize)?;
        if out.is_empty() {
            return Ok(());
        }
        let range = entry.text_range(region.start, region.end)?;
        match gzi {
            None => {
                let len = read_len(&range, source.len())?;
                if len == out.len() {
                    return source.read_exact_at(out, range.start);
                }
                grow(raw, len)?;
                source.read_exact_at(&mut raw[..len], range.start)?;
                entry.decode(region.start, &raw[..len], out)
            }
            Some(gzi) => BgzfRead {
                source,
                gzi,
                inflater: inflater.get_or_insert_with(Inflater::new),
                raw,
                text,
            }
            .decode(entry, region, range, out),
        }
    }
}

pub(crate) struct BgzfRead<'a, S> {
    pub(crate) source: &'a S,
    pub(crate) gzi: &'a Gzi,
    pub(crate) inflater: &'a mut Inflater,
    pub(crate) raw: &'a mut Vec<u8>,
    pub(crate) text: &'a mut Vec<u8>,
}

impl<S: ReadAt> BgzfRead<'_, S> {
    pub(crate) fn block(&self, index: usize, file_len: u64) -> BlockRef {
        let mut block = self.gzi.block(index, file_len);
        // An index may skip empty members. Only the member at this offset
        // is needed, and its BSIZE cannot exceed 64 KiB.
        block.compressed.end = block
            .compressed
            .end
            .min(block.compressed.start.saturating_add(MAX_BLOCK_LEN as u64));
        block
    }

    fn decode(
        &mut self,
        entry: &FaiEntry,
        region: Region<'_>,
        range: Range<u64>,
        out: &mut [u8],
    ) -> Result<()> {
        let file_len = self.source.len();
        let indices = self.gzi.indices(range.clone());
        let mut i = indices.start;
        let mut decoder = TextDecoder::new(entry, region.start, region.end)?;
        let mut text_cursor = self.block(i, file_len).text_start;
        let contiguous = entry.line_bases == entry.line_width;
        while i < indices.end {
            let first = self.block(i, file_len);
            let mut end = first.compressed.end;
            read_len(&first.compressed, file_len)?;
            let mut last = i + 1;
            while last < indices.end {
                let next = self.block(last, file_len);
                if next.compressed.end < next.compressed.start
                    || next.compressed.start != end
                    || next.compressed.end.saturating_sub(first.compressed.start) > READ_BATCH_LEN
                {
                    break;
                }
                end = next.compressed.end;
                last += 1;
            }
            let batch = first.compressed.start..end;
            let len = read_len(&batch, file_len)?;
            grow(self.raw, len)?;
            self.source
                .read_exact_at(&mut self.raw[..len], batch.start)?;
            for k in i..last {
                let b = self.block(k, file_len);
                if b.text_start != text_cursor {
                    return Err(Error::invalid(
                        "gzi",
                        "text offsets disagree with decompressed blocks",
                    ));
                }
                let lo = (b.compressed.start - batch.start) as usize;
                let hi = (b.compressed.end - batch.start) as usize;
                let encoded = self
                    .raw
                    .get(lo..hi)
                    .ok_or_else(|| Error::invalid("gzi", "block outside read batch"))?;
                let direct = contiguous
                    && b.text_start >= range.start
                    && b.text_len.is_some_and(|n| {
                        b.text_start
                            .checked_add(n)
                            .is_some_and(|end| end <= range.end)
                    });
                let (used, n) = if direct {
                    let at = (b.text_start - range.start) as usize;
                    let count = b.text_len.expect("direct block has known length") as usize;
                    if at != decoder.written() {
                        return Err(Error::invalid("gzi", "noncontiguous sequence bytes"));
                    }
                    let dst = out
                        .get_mut(at..at + count)
                        .ok_or_else(|| Error::invalid("gzi", "block outside output range"))?;
                    let result = self.inflater.inflate_block(encoded, dst)?;
                    decoder.advance_contiguous(result.1);
                    result
                } else {
                    grow(self.text, MAX_BLOCK_LEN)?;
                    let result = self
                        .inflater
                        .inflate_block(encoded, &mut self.text[..MAX_BLOCK_LEN])?;
                    push_block(
                        &mut decoder,
                        &self.text[..result.1],
                        b.text_start,
                        &range,
                        out,
                    )?;
                    result
                };
                let n = if n == 0 && b.text_len != Some(0) {
                    let n = self.after_empty_members(b.compressed.start + used as u64, file_len)?;
                    push_block(&mut decoder, &self.text[..n], b.text_start, &range, out)?;
                    n
                } else {
                    n
                };
                if b.text_len.is_some_and(|expected| expected != n as u64) {
                    return Err(Error::invalid(
                        "gzi",
                        "indexed block lengths differ from BGZF data",
                    ));
                }
                text_cursor = b
                    .text_start
                    .checked_add(n as u64)
                    .ok_or_else(|| Error::invalid("gzi", "text offset overflows u64"))?;
            }
            i = last;
        }
        decoder.finish(out.len())
    }
}

impl<S: ReadAt> BgzfRead<'_, S> {
    /// Inflates into `text` the first non-empty member at or after
    /// `offset` and returns its length. htslib indexes a file that begins
    /// with empty members at offset 0 all the same, so the first indexed
    /// block can start with them.
    pub(crate) fn after_empty_members(&mut self, mut offset: u64, file_len: u64) -> Result<usize> {
        let mut raw = Vec::new();
        grow(self.text, MAX_BLOCK_LEN)?;
        let mut filled = (file_len.saturating_sub(offset)).min(MAX_BLOCK_LEN as u64) as usize;
        if filled == 0 {
            return Err(Error::invalid("gzi", "indexed block holds no data"));
        }
        grow(&mut raw, filled)?;
        self.source.read_exact_at(&mut raw[..filled], offset)?;
        let mut window_start = offset;
        loop {
            let at = (offset - window_start) as usize;
            match self
                .inflater
                .inflate_block(&raw[at..filled], &mut self.text[..MAX_BLOCK_LEN])
            {
                Ok((used, n)) => {
                    if n > 0 {
                        return Ok(n);
                    }
                    offset += used as u64;
                }
                Err(Error::Truncated {
                    what: "bgzf header" | "bgzf block",
                    ..
                }) if at != 0 => {
                    // Consume complete members from one bounded window.
                    // A partial header or block is retried at its physical
                    // offset; inflate_block has not written output yet.
                    let len = (file_len.saturating_sub(offset)).min(MAX_BLOCK_LEN as u64) as usize;
                    if len == 0 {
                        return Err(Error::invalid("gzi", "indexed block holds no data"));
                    }
                    grow(&mut raw, len)?;
                    self.source.read_exact_at(&mut raw[..len], offset)?;
                    window_start = offset;
                    filled = len;
                }
                Err(error) => return Err(error),
            }
        }
    }
}

/// Hands the part of an inflated block (`text`, starting at text offset
/// `text_start`) that falls in `range` to the decoder.
pub(crate) fn push_block(
    decoder: &mut TextDecoder<'_>,
    text: &[u8],
    text_start: u64,
    range: &Range<u64>,
    out: &mut [u8],
) -> Result<()> {
    let block_end = text_start
        .checked_add(text.len() as u64)
        .ok_or_else(|| Error::invalid("gzi", "text offset overflows u64"))?;
    let lo = range.start.max(text_start);
    let hi = range.end.min(block_end);
    if lo < hi {
        decoder.push(
            &text[(lo - text_start) as usize..(hi - text_start) as usize],
            out,
        )?;
    }
    Ok(())
}

/// A 2bit file. Opening reads the header, directory and every record header
/// (N and mask blocks) once; each query reads only its packed bases. Bases
/// come back uppercase with `N` runs; with `mask`, soft-masked runs are
/// lowercase, as `twoBitToFa` writes them.
pub struct TwoBitReader<S> {
    pub(crate) source: S,
    pub(crate) directory: Arc<Directory>,
    pub(crate) records: Arc<[Record]>,
    pub(crate) mask: bool,
    pub(crate) packed: Vec<u8>,
    pub(crate) cache: Option<Arc<BlockCache>>,
    pub(crate) namespace: SourceToken,
    pub(crate) tile_bases: u32,
}

fn read<S: ReadAt>(source: &S, offset: u64, len: u64) -> Result<Vec<u8>> {
    let len = len.min(source.len().saturating_sub(offset));
    let len = usize::try_from(len)
        .map_err(|_| Error::invalid("2bit metadata", "read length exceeds usize"))?;
    let mut buf = Vec::new();
    grow(&mut buf, len)?;
    source.read_exact_at(&mut buf, offset)?;
    Ok(buf)
}

impl<S: ReadAt> TwoBitReader<S> {
    pub fn open(source: S, mask: bool) -> Result<TwoBitReader<S>> {
        let header = TwoBitHeader::parse(&read(&source, 0, HEADER_LEN as u64)?)?;
        let bound = header.directory_bound(source.len());
        let buf = read(&source, bound.start, bound.end - bound.start)?;
        let directory = match Directory::parse(&header, &buf)? {
            Parsed::Done(d) => d,
            Parsed::Need(n) => {
                return Err(Error::Truncated {
                    what: "2bit directory",
                    need: n as u64,
                    have: buf.len() as u64,
                });
            }
        };
        let mut prefix = Vec::new();
        let records = directory
            .entries()
            .iter()
            .map(|&(_, offset)| {
                // Every record has a 16-byte header even without N/mask
                // blocks. Fetch that minimum once rather than reading its
                // three count/reserved fields in separate I/O rounds.
                let mut want = 16;
                let mut filled = 0;
                loop {
                    let available = source.len().saturating_sub(offset);
                    let target = (want as u64).min(available) as usize;
                    grow(&mut prefix, target)?;
                    let mut at = offset
                        .checked_add(filled as u64)
                        .ok_or_else(|| Error::invalid("2bit record", "offset overflows u64"))?;
                    if at >= bound.start && at < bound.end {
                        let n = (target - filled).min((bound.end - at) as usize);
                        prefix[filled..filled + n].copy_from_slice(
                            &buf[(at - bound.start) as usize..(at - bound.start) as usize + n],
                        );
                        filled += n;
                        at += n as u64;
                    }
                    if filled < target {
                        source.read_exact_at(&mut prefix[filled..target], at)?;
                        filled = target;
                    }
                    match Record::parse(&prefix[..filled], offset)? {
                        Parsed::Done(record) => {
                            record.check_file_len(source.len())?;
                            return Ok(record);
                        }
                        Parsed::Need(n) if n > filled && filled == want => want = n,
                        Parsed::Need(n) => {
                            return Err(Error::Truncated {
                                what: "2bit record",
                                need: n as u64,
                                have: filled as u64,
                            });
                        }
                    }
                }
            })
            .collect::<Result<Vec<_>>>()?;
        Ok(TwoBitReader {
            directory: Arc::new(directory),
            source,
            records: records.into(),
            mask,
            packed: Vec::new(),
            cache: None,
            namespace: SourceToken::new(),
            tile_bases: crate::DEFAULT_CACHE_TILE_BASES,
        })
    }

    /// Sequence names in file order.
    pub fn names(
        &self,
    ) -> impl ExactSizeIterator<Item = &str> + DoubleEndedIterator + FusedIterator {
        self.directory
            .entries()
            .iter()
            .map(|(name, _)| name.as_str())
    }

    fn record(&self, chrom: &str) -> Result<&Record> {
        self.directory
            .index(chrom)
            .map(|i| &self.records[i])
            .ok_or_else(|| Error::UnknownChrom(chrom.to_owned()))
    }

    /// Attaches an externally managed fixed-tile cache. The input must
    /// remain unchanged while this reader or its clones exist.
    pub fn set_cache(&mut self, cache: Option<Arc<BlockCache>>) -> Result<()> {
        self.cache = cache;
        Ok(())
    }

    /// Chooses the tile size explicitly, in bases. Different tile layouts
    /// have distinct cache keys. Sizes must be powers of two, at least four.
    pub fn set_cache_with_tile(
        &mut self,
        cache: Option<Arc<BlockCache>>,
        tile_bases: u32,
    ) -> Result<()> {
        if tile_bases < 4 || !tile_bases.is_power_of_two() {
            return Err(Error::invalid(
                "2bit cache",
                "tile size must be a power of two >= 4",
            ));
        }
        self.tile_bases = tile_bases;
        self.cache = cache;
        Ok(())
    }

    pub fn cache_namespace(&self) -> SourceToken {
        self.namespace
    }
}

impl<S: ReadAt> SequenceRead for TwoBitReader<S> {
    fn chrom_len(&self, chrom: &str) -> Result<u32> {
        Ok(self.record(chrom)?.len())
    }

    /// Spans are unpacked one by one: 2bit is not compressed, so no
    /// decoding repeats.
    fn sequence(&mut self, locus: &Locus<'_>, out: &mut [u8]) -> Result<()> {
        locus.check()?;
        locus.check_within(self.record(locus.chrom)?.len())?;
        check_output("sequence", out.len(), locus.len() as usize)?;
        if locus.spans().len() == 1 {
            self.read_region(locus.hull(), out)?;
        } else if let Some(cache) = self.cache.as_ref().map(Arc::clone) {
            self.sequence_cached_locus(locus, out, &cache)?;
        } else {
            for (at, region) in locus.pieces() {
                let at = at as usize;
                self.read_region(region, &mut out[at..at + region.len() as usize])?;
            }
        }
        orient_bases(locus.strand, out);
        Ok(())
    }
}

impl<S: ReadAt> TwoBitReader<S> {
    /// The bases of one region, as stored.
    fn read_region(&mut self, region: Region<'_>, out: &mut [u8]) -> Result<()> {
        if let Some(cache) = self.cache.as_ref().map(Arc::clone) {
            return self.sequence_cached(region, out, &cache);
        }
        let TwoBitReader {
            source,
            directory,
            records,
            mask,
            packed,
            ..
        } = self;
        let record = directory
            .index(region.chrom)
            .map(|i| &records[i])
            .ok_or_else(|| Error::UnknownChrom(region.chrom.to_owned()))?;
        region.check_within(record.len())?;
        check_output("sequence", out.len(), region.len() as usize)?;
        if out.is_empty() {
            return Ok(());
        }
        let range = record.packed_range(region.start, region.end)?;
        let len = read_len(&range, source.len())?;
        if len == 0 {
            out.fill(b'N');
            return Ok(());
        }
        grow(packed, len)?;
        source.read_exact_at(&mut packed[..len], range.start)?;
        record.decode(region.start, &packed[..len], out, *mask)
    }
}
