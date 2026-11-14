//! [`BbiReader`]: the standard access to a bigWig or bigBed file, composing
//! the crate's steps over a caller-provided byte source.

use std::ops::Range;

use omicsweave_core::{
    IntervalEntry, IntervalRead, IntervalSize, IntervalsOut, ReadAt, Region, TrackEntry, TrackRead,
    check_output,
};

use crate::block::{BedBlock, WigBlock};
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
    zooms: Vec<ZoomLevel>,
    layout: Layout,
    chroms: Chroms,
    index_base: u64,
    index: Vec<u8>,
    index_nodes: Vec<NodeInfo>,
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
}

/// A feature of the last interval query; its annotation is `arena[note]`.
#[derive(Debug, Clone)]
struct StagedInterval {
    start: u32,
    end: u32,
    note: Range<usize>,
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
        Ok(BbiReader {
            source,
            inflated: Vec::new(),
            fill_scratch: Vec::new(),
            header,
            zooms,
            layout,
            chroms,
            index_base: index_range.start,
            index,
            index_nodes,
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
        })
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
        if size.entries > u32::MAX as usize {
            return Err(Error::invalid("BBI query", "more than 2^32-1 features"));
        }
        // The high-level interface promises start order even when blocks or
        // their records are unordered. Stable sorting retains stored order
        // for equal starts and keeps each annotation attached to its record.
        if !self.intervals.is_sorted_by_key(|entry| entry.start) {
            self.intervals.sort_by_key(|entry| entry.start);
        }
        Ok(size)
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

    fn values(&mut self, region: Region<'_>, out: &mut [f32]) -> Result<()> {
        self.expect(Kind::BigWig)?;
        let chrom = self.track_chrom(&region)?;
        check_output("values", out.len(), region.len() as usize)?;
        out.fill(f32::NAN);
        if region.is_empty() {
            return Ok(());
        }
        self.find_blocks(chrom, region.start, region.end)?;
        let mut scratch = std::mem::take(&mut self.fill_scratch);
        let result = (|| {
            for i in 0..self.blocks.len() {
                WigBlock::parse(self.load(i)?)?.fill_reusing(
                    chrom,
                    region.start,
                    out,
                    &mut scratch,
                )?;
            }
            Ok(())
        })();
        self.fill_scratch = scratch;
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
