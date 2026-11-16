//! [`BigWigWriter`] and [`BigBedWriter`]: write standard BBI files (version
//! 4) in one streaming pass, and return the [`BbiSidecar`] of what they wrote.
//!
//! See [`WriteOptions`] for the layout and the zoom levels written.
//!
//! A writer cuts items into blocks and zoom summaries on the caller's
//! thread, compresses each block using an [`Exec`], and writes the
//! compressed blocks back in the order they were cut. Offsets, the R-trees
//! and the sidecar come from that ordered commit, so the file and sidecar
//! are the same for every executor and thread count.

use std::cell::RefCell;
use std::cmp::Reverse;
use std::collections::{BinaryHeap, VecDeque};
use std::io::{Seek, SeekFrom, Write};

use omicsweave_core::exec::{Exec, Ordered};

use crate::chroms::{CHROM_TREE_MAGIC, Chrom, ChromId};
use crate::error::{Result, invalid};
use crate::header::{
    BIGBED_MAGIC, BIGWIG_MAGIC, HEADER_LEN, Kind, TOTAL_SUMMARY_LEN, ZOOM_HEADER_LEN,
};
use crate::index::INDEX_MAGIC;
use crate::sidecar::{BbiSidecar, BlockInfo, LevelBlocks};

/// File format version written.
pub const BBI_WRITE_VERSION: u16 = 4;
/// Items whose mean span sets the first derived zoom reduction.
pub const ZOOM_SAMPLE: usize = 10_000;

const EXTENSION_LEN: usize = 64;
const INDEX_HEADER_LEN: usize = 48;
const RTREE_LEAF_LEN: usize = 32;
const RTREE_BRANCH_LEN: usize = 24;
const ZOOM_RECORD_LEN: usize = 32;
const WIG_HEADER_LEN: usize = 24;
const BEDGRAPH: u8 = 1;

/// How a file is laid out and compressed.
///
/// Sections follow the UCSC tools: header, zoom headers, autoSql (bigBed),
/// total summary, extension header, chromosome B+ tree, data, R-tree index,
/// each zoom level's data and index, trailing magic. Tree nodes are padded to
/// the tree's block size, as `bedGraphToBigWig` and `bedToBigBed` pad them.
///
/// Chromosome ids are the positions of the chromosomes as given; the B+
/// tree is keyed by name. Items must come sorted by (chromosome id, start).
/// Base-level blocks never span chromosomes; zoom blocks may.
///
/// **Zoom levels.** Unless given explicitly, the first reduction is ten
/// times the mean item span of the first [`ZOOM_SAMPLE`] items (at least
/// 10), each further level four times the previous, at most
/// [`WriteOptions::zoom_levels`] levels and none past the longest
/// chromosome; a derived level is kept only while it has fewer summaries
/// than the level below. Summaries start at the first covered base and span
/// `reduction` bases (clipped to the chromosome), as `bbiAddRangeToSummary`
/// makes them. bigWig summaries are over item values; bigBed summaries, and
/// the bigBed total summary, are over coverage depth. Zoom data is kept in
/// memory, compressed, until [`BigWigWriter::finish`]: about the size of the
/// file's zoom sections.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct WriteOptions {
    /// Items per data block; `None` takes the UCSC default (1024 for bigWig,
    /// 512 for bigBed). Zoom blocks hold as many summaries.
    pub items_per_slot: Option<u32>,
    /// Children per tree node (R-trees and the chromosome tree), 2 to 65535.
    pub block_size: u32,
    /// zlib-compress blocks.
    pub compress: bool,
    /// libdeflate level, 0 to 12.
    pub level: i32,
    /// Most zoom levels derived.
    pub zoom_levels: usize,
    /// Explicit zoom reductions, written as given (each above 0), instead of
    /// derived ones.
    pub zoom_reductions: Option<Vec<u32>>,
}

impl Default for WriteOptions {
    /// The UCSC defaults: block size 256, compressed, up to 10 zoom levels;
    /// level 6 is zlib's default.
    fn default() -> WriteOptions {
        WriteOptions {
            items_per_slot: None,
            block_size: 256,
            compress: true,
            level: 6,
            zoom_levels: 10,
            zoom_reductions: None,
        }
    }
}

/// The record layout of a bigBed file.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct BedSchema {
    /// Fields per record, the three coordinates included.
    pub field_count: u16,
    /// Leading fields that follow the BED standard.
    pub defined_field_count: u16,
    /// The autoSql description, stored NUL-terminated; `None` stores none.
    pub auto_sql: Option<String>,
}

impl BedSchema {
    /// Plain BED with `fields` standard fields and no autoSql.
    pub fn bed(fields: u16) -> BedSchema {
        BedSchema {
            field_count: fields,
            defined_field_count: fields,
            auto_sql: None,
        }
    }
}

/// Little-endian appends.
#[derive(Default)]
struct Bytes(Vec<u8>);

impl Bytes {
    fn u8(&mut self, v: u8) {
        self.0.push(v);
    }
    fn u16(&mut self, v: u16) {
        self.0.extend_from_slice(&v.to_le_bytes());
    }
    fn u32(&mut self, v: u32) {
        self.0.extend_from_slice(&v.to_le_bytes());
    }
    fn u64(&mut self, v: u64) {
        self.0.extend_from_slice(&v.to_le_bytes());
    }
    fn f32(&mut self, v: f32) {
        self.0.extend_from_slice(&v.to_le_bytes());
    }
    fn f64(&mut self, v: f64) {
        self.0.extend_from_slice(&v.to_le_bytes());
    }
    fn zeros(&mut self, n: usize) {
        self.0.resize(self.0.len() + n, 0);
    }
}

type Point = (ChromId, u32);

/// A zoom summary or total summary being accumulated.
#[derive(Debug, Clone, Copy)]
struct Acc {
    chrom: ChromId,
    start: u32,
    end: u32,
    valid: u64,
    min: f64,
    max: f64,
    sum: f64,
    squares: f64,
}

impl Acc {
    fn new(chrom: ChromId, start: u32, end: u32) -> Acc {
        Acc {
            chrom,
            start,
            end,
            valid: 0,
            min: f64::INFINITY,
            max: f64::NEG_INFINITY,
            sum: 0.0,
            squares: 0.0,
        }
    }

    fn add(&mut self, bases: u32, value: f64) {
        let n = bases as f64;
        self.valid += bases as u64;
        self.min = self.min.min(value);
        self.max = self.max.max(value);
        self.sum += n * value;
        self.squares += n * value * value;
    }
}

/// A block being filled: its uncompressed bytes, span and item count.
#[derive(Default)]
struct Pending {
    raw: Vec<u8>,
    items: u32,
    start: Point,
    end: Point,
}

impl Pending {
    fn extend_end(&mut self, end: Point) {
        if self.items == 0 || end > self.end {
            self.end = end;
        }
    }
}

/// The bytes to store for `raw`: zlib at libdeflate `level`, or `raw`
/// itself when uncompressed.
fn pack(level: Option<i32>, raw: Vec<u8>) -> Result<Vec<u8>> {
    let Some(level) = level else {
        return Ok(raw);
    };
    let mut packed = Vec::new();
    let n = compress(level, &raw, &mut packed)?;
    packed.truncate(n);
    Ok(packed)
}

/// Keeps the full output bound initialized so repeated inline blocks can
/// overwrite the same buffer without allocating or zeroing its tail again.
fn compress(level: i32, raw: &[u8], packed: &mut Vec<u8>) -> Result<usize> {
    thread_local! {
        // A libdeflate compressor allocates its matchfinder tables (hundreds
        // of KiB) when created. Jobs run on whichever pool thread is free and
        // carry no state, so each thread keeps one compressor per level it
        // has used instead of building one per block.
        static COMPRESSORS: RefCell<Vec<(i32, libdeflater::Compressor)>> =
            const { RefCell::new(Vec::new()) };
    }
    COMPRESSORS.with_borrow_mut(|compressors| {
        let i = match compressors.iter().position(|(l, _)| *l == level) {
            Some(i) => i,
            None => {
                let lvl = libdeflater::CompressionLvl::new(level)
                    .map_err(|_| invalid("BBI writer", "compression level must be 0 to 12"))?;
                compressors.push((level, libdeflater::Compressor::new(lvl)));
                compressors.len() - 1
            }
        };
        let c = &mut compressors[i].1;
        let bound = c.zlib_compress_bound(raw.len());
        if packed.len() < bound {
            packed.resize(bound, 0);
        }
        c.zlib_compress(raw, packed)
            .map_err(|_| invalid("BBI writer", "compression failed"))
    })
}

/// One zoom level: the summary being accumulated and the block being
/// filled.
struct ZoomAcc {
    reduction: u32,
    current: Option<Acc>,
    block: Pending,
    records: u64,
}

impl ZoomAcc {
    fn new(reduction: u32) -> ZoomAcc {
        ZoomAcc {
            reduction,
            current: None,
            block: Pending::default(),
            records: 0,
        }
    }
}

/// A zoom level's compressed blocks, kept until the file's zoom section.
#[derive(Default)]
struct ZoomOut {
    data: Vec<u8>,
    /// `offset` is relative to `data` until written.
    blocks: Vec<BlockInfo>,
}

/// A block cut and submitted, waiting for its compressed bytes: its level
/// (`None` for base data) and its sidecar entry, without offset and size.
struct Sealed {
    level: Option<usize>,
    info: BlockInfo,
}

/// The ordered half of a writer: compression jobs in flight, and where
/// their results go. Blocks are committed in the order they were sealed, so
/// offsets do not depend on how the jobs ran.
struct Sink<W: Write> {
    out: W,
    position: u64,
    failed: bool,
    /// Inline blocks borrow their Pending bytes directly; only pool jobs
    /// need owned buffers and an ordered queue.
    exec: Option<Box<dyn Ordered<Vec<u8>>>>,
    inline_packed: Vec<u8>,
    level: Option<i32>,
    max_inflated: usize,
    sealed: VecDeque<Sealed>,
    blocks: Vec<BlockInfo>,
    zooms: Vec<ZoomOut>,
}

impl<W: Write> Sink<W> {
    fn write(&mut self, bytes: &[u8]) -> Result<u64> {
        if self.failed {
            return Err(invalid("BBI writer", "an earlier write or block failed"));
        }
        let at = self.position;
        if let Err(e) = self.out.write_all(bytes) {
            self.failed = true;
            return Err(e.into());
        }
        self.position += bytes.len() as u64;
        Ok(at)
    }

    fn check(&self) -> Result<()> {
        if self.failed {
            return Err(invalid("BBI writer", "an earlier write or block failed"));
        }
        Ok(())
    }

    /// Takes the next result in order. A failed job fails the writer, as
    /// the blocks after it can no longer be placed; so does a refused
    /// submission.
    fn next(&mut self, wait: bool) -> Result<Option<Vec<u8>>> {
        match &mut self.exec {
            Some(exec) => exec.next(wait).inspect_err(|_| self.failed = true),
            None => Ok(None),
        }
    }

    /// Seals `block` of `level`, leaving it empty. Inline work borrows the
    /// bytes; pool work owns them until its ordered result is committed.
    fn seal(&mut self, level: Option<usize>, block: &mut Pending) -> Result<()> {
        self.check()?;
        self.max_inflated = self.max_inflated.max(block.raw.len());
        let info = BlockInfo {
            start: block.start,
            end: block.end,
            offset: 0,
            size: 0,
            inflated: block.raw.len() as u32,
            items: block.items,
        };
        let sealed = Sealed { level, info };
        if self.exec.is_none() {
            let result = if let Some(level) = self.level {
                let mut packed = std::mem::take(&mut self.inline_packed);
                let result = compress(level, &block.raw, &mut packed)
                    .and_then(|n| self.commit_bytes(sealed, &packed[..n]));
                self.inline_packed = packed;
                result
            } else {
                self.commit_bytes(sealed, &block.raw)
            };
            if result.is_err() {
                self.failed = true;
            }
            result?;
            block.raw.clear();
            block.items = 0;
            return Ok(());
        }
        let capacity = block.raw.len();
        let raw = std::mem::replace(&mut block.raw, Vec::with_capacity(capacity));
        block.items = 0;
        while self.exec.as_ref().expect("pool queue").full() {
            let packed = self
                .next(true)?
                .ok_or_else(|| invalid("BBI writer", "compression job lost"))?;
            self.commit(packed)?;
        }
        let level_of = self.level;
        self.exec
            .as_mut()
            .expect("pool queue")
            .submit(Box::new(move || pack(level_of, raw)))
            .inspect_err(|_| self.failed = true)?;
        self.sealed.push_back(sealed);
        while let Some(packed) = self.next(false)? {
            self.commit(packed)?;
        }
        Ok(())
    }

    /// Places the oldest sealed block's stored bytes.
    fn commit(&mut self, packed: Vec<u8>) -> Result<()> {
        let sealed = self.sealed.pop_front().expect("results follow seals");
        self.commit_bytes(sealed, &packed)
    }

    fn commit_bytes(&mut self, sealed: Sealed, packed: &[u8]) -> Result<()> {
        let Sealed { level, mut info } = sealed;
        info.size = packed.len() as u32;
        match level {
            None => {
                info.offset = self.write(packed)?;
                self.blocks.push(info);
            }
            Some(l) => {
                let z = &mut self.zooms[l];
                info.offset = z.data.len() as u64;
                z.data.extend_from_slice(packed);
                z.blocks.push(info);
            }
        }
        Ok(())
    }

    /// Commits every block in flight.
    fn drain(&mut self) -> Result<()> {
        self.check()?;
        while let Some(packed) = self.next(true)? {
            self.commit(packed)?;
        }
        Ok(())
    }
}

/// What both writers share: output, layout, base blocks, zoom levels and
/// the total summary.
struct Core<W: Write + Seek> {
    sink: Sink<W>,
    kind: Kind,
    chroms: Vec<Chrom>,
    items_per_slot: u32,
    block_size: u32,
    max_zooms: usize,
    explicit_zooms: bool,
    zoom_slots: usize,
    schema: Option<BedSchema>,
    auto_sql_offset: u64,
    summary_offset: u64,
    extension_offset: u64,
    chrom_tree_offset: u64,
    data_offset: u64,
    block: Pending,
    /// Items added, for the zoom level rule.
    items: u64,
    total: Option<Acc>,
    /// Item spans sampled for the first reduction, and the zoom input
    /// buffered until the reductions are fixed.
    sample: (u64, u64),
    buffered: Vec<(ChromId, u32, u32, f64)>,
    zooms: Option<Vec<ZoomAcc>>,
}

impl<W: Write + Seek> Core<W> {
    fn create(
        out: W,
        kind: Kind,
        chroms: Vec<Chrom>,
        schema: Option<BedSchema>,
        options: WriteOptions,
        exec: Exec<'_>,
    ) -> Result<Core<W>> {
        if chroms.is_empty() || chroms.len() > u32::MAX as usize {
            return Err(invalid("BBI writer", "needs 1 to 2^32 - 1 chromosomes"));
        }
        let mut names = std::collections::HashSet::new();
        for c in &chroms {
            if c.name.is_empty() || c.name.contains('\0') || !names.insert(c.name.as_str()) {
                return Err(invalid(
                    "BBI writer",
                    format!("bad or duplicate chromosome name {:?}", c.name),
                ));
            }
        }
        let items_per_slot = options.items_per_slot.unwrap_or(match kind {
            Kind::BigWig => 1024,
            Kind::BigBed => 512,
        });
        if items_per_slot == 0 || items_per_slot > u16::MAX as u32 {
            return Err(invalid("BBI writer", "items per slot must be 1 to 65535"));
        }
        if !(2..=u16::MAX as u32).contains(&options.block_size) {
            return Err(invalid("BBI writer", "block size must be 2 to 65535"));
        }
        let level = if options.compress {
            libdeflater::CompressionLvl::new(options.level)
                .map_err(|_| invalid("BBI writer", "compression level must be 0 to 12"))?;
            Some(options.level)
        } else {
            None
        };
        let (zooms, zoom_slots) = match &options.zoom_reductions {
            Some(r) => {
                if r.contains(&0) || r.len() > u16::MAX as usize {
                    return Err(invalid("BBI writer", "zoom reductions must be positive"));
                }
                (
                    Some(r.iter().map(|&r| ZoomAcc::new(r)).collect::<Vec<_>>()),
                    r.len(),
                )
            }
            None if options.zoom_levels == 0 => (Some(Vec::new()), 0),
            None => (None, options.zoom_levels.min(u16::MAX as usize)),
        };
        if let Some(s) = &schema {
            if s.field_count < 3 || s.defined_field_count > s.field_count {
                return Err(invalid("bigBed schema", "bad field counts"));
            }
            if s.auto_sql.as_ref().is_some_and(|a| a.contains('\0')) {
                return Err(invalid("bigBed schema", "autoSql contains NUL"));
            }
        }
        let zoom_data = zooms.as_ref().map_or(0, |z| z.len());
        let mut core = Core {
            sink: Sink {
                out,
                position: 0,
                failed: false,
                exec: match exec {
                    Exec::Inline => None,
                    Exec::Pool { .. } => Some(exec.queue()),
                },
                inline_packed: Vec::new(),
                level,
                max_inflated: 0,
                sealed: VecDeque::new(),
                blocks: Vec::new(),
                zooms: (0..zoom_data).map(|_| ZoomOut::default()).collect(),
            },
            kind,
            chroms,
            items_per_slot,
            block_size: options.block_size,
            max_zooms: zoom_slots,
            explicit_zooms: options.zoom_reductions.is_some(),
            zoom_slots,
            schema,
            auto_sql_offset: 0,
            summary_offset: 0,
            extension_offset: 0,
            chrom_tree_offset: 0,
            data_offset: 0,
            block: Pending::default(),
            items: 0,
            total: None,
            sample: (0, 0),
            buffered: Vec::new(),
            zooms,
        };
        core.write_head()?;
        Ok(core)
    }

    fn write(&mut self, bytes: &[u8]) -> Result<u64> {
        self.sink.write(bytes)
    }

    /// Everything before the data blocks; offsets and sizes known only at
    /// the end are patched by `finish`.
    fn write_head(&mut self) -> Result<()> {
        let mut b = Bytes::default();
        b.zeros(HEADER_LEN + self.zoom_slots * ZOOM_HEADER_LEN);
        if let Some(sql) = self.schema.as_ref().and_then(|s| s.auto_sql.as_ref()) {
            self.auto_sql_offset = b.0.len() as u64;
            b.0.extend_from_slice(sql.as_bytes());
            b.u8(0);
        }
        self.summary_offset = b.0.len() as u64;
        b.zeros(TOTAL_SUMMARY_LEN);
        self.extension_offset = b.0.len() as u64;
        b.u16(EXTENSION_LEN as u16);
        b.zeros(EXTENSION_LEN - 2);
        self.chrom_tree_offset = b.0.len() as u64;
        chrom_tree(
            &mut b,
            &self.chroms,
            self.block_size,
            self.chrom_tree_offset,
        );
        self.data_offset = b.0.len() as u64;
        b.u64(0); // data count, patched
        self.write(&b.0)?;
        Ok(())
    }

    fn check_item(&self, chrom: ChromId, start: u32, end: u32, last: Option<Point>) -> Result<()> {
        let Some(c) = self.chroms.get(chrom as usize) else {
            return Err(invalid(
                "BBI item",
                format!("chromosome id {chrom} out of range"),
            ));
        };
        if start > end || end > c.len {
            return Err(invalid(
                "BBI item",
                format!("{}:{start}-{end} outside 0..{}", c.name, c.len),
            ));
        }
        if last.is_some_and(|l| (chrom, start) < l) {
            return Err(invalid(
                "BBI item",
                "items must be sorted by chromosome and start",
            ));
        }
        Ok(())
    }

    /// Starts a new base block when the current one is full or on another
    /// chromosome.
    fn room(&mut self, chrom: ChromId) -> Result<()> {
        if self.block.items > 0
            && (self.block.items == self.items_per_slot || self.block.start.0 != chrom)
        {
            self.flush_block()?;
        }
        Ok(())
    }

    fn flush_block(&mut self) -> Result<()> {
        if self.block.items == 0 {
            return Ok(());
        }
        let block = &mut self.block;
        if self.kind == Kind::BigWig {
            let head = &mut block.raw[..WIG_HEADER_LEN];
            for (dst, value) in head[..20].chunks_exact_mut(4).zip([
                block.start.0,
                block.start.1,
                block.end.1,
                0,
                0,
            ]) {
                dst.copy_from_slice(&value.to_le_bytes());
            }
            head[20] = BEDGRAPH;
            head[21] = 0;
            head[22..24].copy_from_slice(&(block.items as u16).to_le_bytes());
        }
        self.sink.seal(None, block)
    }

    /// Opens the base block for an item if it is empty.
    fn begin(&mut self, chrom: ChromId, start: u32) {
        if self.block.items == 0 {
            self.block.start = (chrom, start);
            self.block.raw.clear();
            if self.kind == Kind::BigWig {
                self.block.raw.resize(WIG_HEADER_LEN, 0);
            }
        }
    }

    /// Feeds a run of `value` over `[start, end)` to the total summary and
    /// the zoom levels (fixing the reductions once enough items are seen).
    fn summarize(&mut self, chrom: ChromId, start: u32, end: u32, value: f64) -> Result<()> {
        if start >= end {
            return Ok(());
        }
        let t = self.total.get_or_insert(Acc::new(0, 0, 0));
        t.add(end - start, value);
        if self.zooms.is_none() {
            self.buffered.push((chrom, start, end, value));
            if self.sample.0 >= ZOOM_SAMPLE as u64 {
                self.fix_zooms()?;
            }
            return Ok(());
        }
        self.zoom_add(chrom, start, end, value)
    }

    fn fix_zooms(&mut self) -> Result<()> {
        let (n, span) = self.sample;
        let mut reductions = Vec::new();
        if let Some(mean) = span.checked_div(n) {
            let longest = self.chroms.iter().map(|c| c.len).max().unwrap_or(0);
            let mut r = mean.saturating_mul(10).max(10);
            while reductions.len() < self.max_zooms && r <= longest as u64 {
                reductions.push(r as u32);
                r = r.saturating_mul(4);
            }
        }
        self.sink.zooms = (0..reductions.len()).map(|_| ZoomOut::default()).collect();
        self.zooms = Some(reductions.into_iter().map(ZoomAcc::new).collect());
        for (c, s, e, v) in std::mem::take(&mut self.buffered) {
            self.zoom_add(c, s, e, v)?;
        }
        Ok(())
    }

    fn zoom_add(&mut self, chrom: ChromId, start: u32, end: u32, value: f64) -> Result<()> {
        let len = self.chroms[chrom as usize].len;
        let zooms = self.zooms.as_mut().expect("reductions fixed");
        for (level, z) in zooms.iter_mut().enumerate() {
            let mut s = start;
            while s < end {
                let fresh = match &z.current {
                    Some(a) => a.chrom != chrom || s >= a.end,
                    None => true,
                };
                if fresh {
                    if let Some(a) = z.current.take() {
                        push_summary(z, level, &a, self.items_per_slot, &mut self.sink)?;
                    }
                    z.current = Some(Acc::new(chrom, s, s.saturating_add(z.reduction).min(len)));
                }
                let a = z.current.as_mut().expect("summary open");
                let e = end.min(a.end);
                a.add(e - s, value);
                s = e;
            }
        }
        Ok(())
    }

    fn finish(mut self) -> Result<(W, BbiSidecar)> {
        self.flush_block()?;
        if self.zooms.is_none() {
            self.fix_zooms()?;
        }
        let mut zooms = self.zooms.take().expect("reductions fixed");
        for (level, z) in zooms.iter_mut().enumerate() {
            if let Some(a) = z.current.take() {
                push_summary(z, level, &a, self.items_per_slot, &mut self.sink)?;
            }
            seal_zoom_block(z, level, &mut self.sink)?;
        }
        self.sink.drain()?;
        let mut zoom_data = std::mem::take(&mut self.sink.zooms);
        if !self.explicit_zooms {
            let mut below = self.items;
            let keep = zooms
                .iter()
                .take_while(|z| {
                    let fewer = z.records > 0 && z.records < below;
                    below = z.records;
                    fewer
                })
                .count();
            zooms.truncate(keep);
            zoom_data.truncate(keep);
        }

        let index_offset = self.sink.position;
        let mut b = Bytes::default();
        rtree(
            &mut b,
            &self.sink.blocks,
            self.block_size,
            self.items_per_slot,
            index_offset,
            index_offset,
        );
        self.write(&b.0)?;

        let mut levels = vec![LevelBlocks {
            reduction: 0,
            blocks: std::mem::take(&mut self.sink.blocks),
        }];
        let mut zoom_heads = Vec::with_capacity(zooms.len());
        for (z, data) in zooms.iter().zip(zoom_data) {
            let data_offset = self.write(&(data.blocks.len() as u32).to_le_bytes())?;
            let base = self.write(&data.data)?;
            let blocks: Vec<BlockInfo> = data
                .blocks
                .iter()
                .map(|b| BlockInfo {
                    offset: b.offset + base,
                    ..*b
                })
                .collect();
            let index_offset = self.sink.position;
            let mut b = Bytes::default();
            rtree(
                &mut b,
                &blocks,
                self.block_size,
                self.items_per_slot,
                index_offset,
                index_offset,
            );
            self.write(&b.0)?;
            zoom_heads.push((z.reduction, data_offset, index_offset));
            levels.push(LevelBlocks {
                reduction: z.reduction,
                blocks,
            });
        }
        let magic = match self.kind {
            Kind::BigWig => BIGWIG_MAGIC,
            Kind::BigBed => BIGBED_MAGIC,
        };
        self.write(&magic.to_le_bytes())?;
        let end = self.sink.position;

        let uncompress_buf_size = if self.sink.level.is_some() {
            u32::try_from(self.sink.max_inflated)
                .map_err(|_| invalid("BBI writer", "block exceeds 4 GiB"))?
        } else {
            0
        };
        let (fields, defined) = self
            .schema
            .as_ref()
            .map_or((0, 0), |s| (s.field_count, s.defined_field_count));
        let mut h = Bytes::default();
        h.u32(magic);
        h.u16(BBI_WRITE_VERSION);
        h.u16(zoom_heads.len() as u16);
        h.u64(self.chrom_tree_offset);
        h.u64(self.data_offset);
        h.u64(index_offset);
        h.u16(fields);
        h.u16(defined);
        h.u64(self.auto_sql_offset);
        h.u64(self.summary_offset);
        h.u32(uncompress_buf_size);
        h.u64(self.extension_offset);
        for &(reduction, data, index) in &zoom_heads {
            h.u32(reduction);
            h.u32(0);
            h.u64(data);
            h.u64(index);
        }
        let mut s = Bytes::default();
        match self.total {
            Some(t) => {
                s.u64(t.valid);
                s.f64(t.min);
                s.f64(t.max);
                s.f64(t.sum);
                s.f64(t.squares);
            }
            None => s.zeros(TOTAL_SUMMARY_LEN),
        }
        let count = match self.kind {
            Kind::BigWig => levels[0].blocks.len() as u64,
            Kind::BigBed => self.items,
        };
        let patch = |out: &mut W, at: u64, bytes: &[u8]| -> std::io::Result<()> {
            out.seek(SeekFrom::Start(at))?;
            out.write_all(bytes)
        };
        let out = &mut self.sink.out;
        let patched = (|| {
            patch(out, 0, &h.0)?;
            patch(out, self.summary_offset, &s.0)?;
            patch(out, self.data_offset, &count.to_le_bytes())?;
            out.seek(SeekFrom::Start(end))?;
            out.flush()
        })();
        patched?;
        Ok((
            self.sink.out,
            BbiSidecar {
                kind: self.kind,
                uncompress_buf_size,
                levels,
            },
        ))
    }
}

/// Appends a finished summary to its level's block, sealing the block when
/// it is full.
fn push_summary<W: Write>(
    z: &mut ZoomAcc,
    level: usize,
    a: &Acc,
    items_per_slot: u32,
    sink: &mut Sink<W>,
) -> Result<()> {
    if z.block.items == items_per_slot {
        seal_zoom_block(z, level, sink)?;
    }
    if z.block.items == 0 {
        z.block.start = (a.chrom, a.start);
        z.block.raw.clear();
    }
    z.block.extend_end((a.chrom, a.end));
    let valid = u32::try_from(a.valid).map_err(|_| invalid("zoom summary", "too many bases"))?;
    let mut b = Bytes(std::mem::take(&mut z.block.raw));
    b.u32(a.chrom);
    b.u32(a.start);
    b.u32(a.end);
    b.u32(valid);
    b.f32(a.min as f32);
    b.f32(a.max as f32);
    b.f32(a.sum as f32);
    b.f32(a.squares as f32);
    z.block.raw = b.0;
    z.block.items += 1;
    z.records += 1;
    Ok(())
}

fn seal_zoom_block<W: Write>(z: &mut ZoomAcc, level: usize, sink: &mut Sink<W>) -> Result<()> {
    if z.block.items == 0 {
        return Ok(());
    }
    debug_assert_eq!(z.block.raw.len(), z.block.items as usize * ZOOM_RECORD_LEN);
    sink.seal(Some(level), &mut z.block)
}

/// Nodes per tree level, leaves first, for `n` items and fanout `b`.
fn level_sizes(n: usize, b: usize) -> Vec<usize> {
    let mut sizes = vec![n.div_ceil(b).max(1)];
    while *sizes.last().unwrap() > 1 {
        sizes.push(sizes.last().unwrap().div_ceil(b));
    }
    sizes
}

/// Writes a tree root level first; `leaf(b, i)` and `branch(b, level,
/// node)` append one item each (a branch item without its child offset,
/// which is appended here). Nodes are padded to `fanout` items.
#[allow(clippy::too_many_arguments)]
fn tree(
    b: &mut Bytes,
    n: usize,
    fanout: usize,
    base: u64,
    leaf_len: usize,
    branch_len: usize,
    leaf: &dyn Fn(&mut Bytes, usize),
    branch: &dyn Fn(&mut Bytes, usize, usize),
) {
    let sizes = level_sizes(n, fanout);
    let node_len = |level: usize| 4 + fanout * if level == 0 { leaf_len } else { branch_len };
    // File offset of each level's first node, root level first.
    let mut starts = vec![0u64; sizes.len()];
    let mut at = base;
    for level in (0..sizes.len()).rev() {
        starts[level] = at;
        at += (sizes[level] * node_len(level)) as u64;
    }
    for level in (0..sizes.len()).rev() {
        let below = if level == 0 { n } else { sizes[level - 1] };
        for node in 0..sizes[level] {
            let first = node * fanout;
            let count = below.saturating_sub(first).min(fanout);
            b.u8((level == 0) as u8);
            b.u8(0);
            b.u16(count as u16);
            for i in first..first + count {
                if level == 0 {
                    leaf(b, i);
                } else {
                    branch(b, level - 1, i);
                    b.u64(starts[level - 1] + (i * node_len(level - 1)) as u64);
                }
            }
            let item = if level == 0 { leaf_len } else { branch_len };
            b.zeros((fanout - count) * item);
        }
    }
}

/// The chromosome B+ tree, keyed by name.
fn chrom_tree(b: &mut Bytes, chroms: &[Chrom], block_size: u32, base: u64) {
    let mut order: Vec<usize> = (0..chroms.len()).collect();
    order.sort_by(|&x, &y| chroms[x].name.as_bytes().cmp(chroms[y].name.as_bytes()));
    let key_size = chroms.iter().map(|c| c.name.len()).max().unwrap_or(1);
    let fanout = (block_size as usize).min(chroms.len()).max(1);
    b.u32(CHROM_TREE_MAGIC);
    b.u32(fanout as u32);
    b.u32(key_size as u32);
    b.u32(8);
    b.u64(chroms.len() as u64);
    b.u64(0);
    let key = |b: &mut Bytes, i: usize| {
        let name = chroms[order[i]].name.as_bytes();
        b.0.extend_from_slice(name);
        b.zeros(key_size - name.len());
    };
    // A node at `level` (0: leaves) covers fanout^(level+1) chromosomes.
    let span = |level: usize| fanout.pow(level as u32 + 1);
    tree(
        b,
        chroms.len(),
        fanout,
        base + 32,
        key_size + 8,
        key_size + 8,
        &|b, i| {
            key(b, i);
            b.u32(order[i] as u32);
            b.u32(chroms[order[i]].len);
        },
        &|b, level, node| key(b, node * span(level)),
    );
}

/// An R-tree over `blocks` (sorted by start), whose header is at `base`;
/// `data_end` is the end of the data it indexes.
fn rtree(
    b: &mut Bytes,
    blocks: &[BlockInfo],
    block_size: u32,
    items_per_slot: u32,
    base: u64,
    data_end: u64,
) {
    let fanout = (block_size as usize).min(blocks.len()).max(1);
    let sizes = level_sizes(blocks.len(), fanout);
    // Span of every node, leaves first: start of its first block, the
    // largest end among its blocks.
    let mut spans: Vec<Vec<(Point, Point)>> = Vec::with_capacity(sizes.len());
    let leaf_spans: Vec<(Point, Point)> = blocks
        .chunks(fanout)
        .map(|c| (c[0].start, c.iter().map(|x| x.end).max().unwrap()))
        .collect();
    spans.push(leaf_spans);
    for level in 1..sizes.len() {
        let below = &spans[level - 1];
        let up = below
            .chunks(fanout)
            .map(|c| (c[0].0, c.iter().map(|x| x.1).max().unwrap()))
            .collect();
        spans.push(up);
    }
    let (first, last) = match spans.last().and_then(|s| s.first()) {
        Some(&s) => s,
        None => ((0, 0), (0, 0)),
    };
    b.u32(INDEX_MAGIC);
    b.u32(fanout as u32);
    b.u64(blocks.len() as u64);
    b.u32(first.0);
    b.u32(first.1);
    b.u32(last.0);
    b.u32(last.1);
    b.u64(data_end);
    b.u32(items_per_slot);
    b.u32(0);
    tree(
        b,
        blocks.len(),
        fanout,
        base + INDEX_HEADER_LEN as u64,
        RTREE_LEAF_LEN,
        RTREE_BRANCH_LEN,
        &|b, i| {
            let x = &blocks[i];
            b.u32(x.start.0);
            b.u32(x.start.1);
            b.u32(x.end.0);
            b.u32(x.end.1);
            b.u64(x.offset);
            b.u64(x.size as u64);
        },
        &|b, level, node| {
            let (s, e) = spans[level][node];
            b.u32(s.0);
            b.u32(s.1);
            b.u32(e.0);
            b.u32(e.1);
        },
    );
}

/// Writes a bigWig file of bedGraph items; see [`WriteOptions`].
pub struct BigWigWriter<W: Write + Seek> {
    core: Core<W>,
    last: Option<Point>,
}

impl<W: Write + Seek> BigWigWriter<W> {
    /// Writes the head of a file to `out`, which must start empty at
    /// position 0 (wrap a file in a `BufWriter`). Blocks are compressed on
    /// the caller's thread.
    pub fn create(out: W, chroms: Vec<Chrom>, options: WriteOptions) -> Result<BigWigWriter<W>> {
        BigWigWriter::create_with(out, chroms, options, Exec::Inline)
    }

    /// [`BigWigWriter::create`], compressing blocks on `exec`; the file and
    /// sidecar are the same for every executor. With [`Exec::Pool`], at
    /// most `limit` blocks are in flight, and an output or compression
    /// error may surface on a later call than the one that cut the block.
    /// Such an error, or a job the executor refuses (its pool was dropped),
    /// fails every later call.
    pub fn create_with(
        out: W,
        chroms: Vec<Chrom>,
        options: WriteOptions,
        exec: Exec<'_>,
    ) -> Result<BigWigWriter<W>> {
        Ok(BigWigWriter {
            core: Core::create(out, Kind::BigWig, chroms, None, options, exec)?,
            last: None,
        })
    }

    /// Adds `value` over `[start, end)` of chromosome `chrom`. Items come
    /// sorted by chromosome id and start and do not overlap; empty and
    /// non-finite items are rejected.
    pub fn add(&mut self, chrom: ChromId, start: u32, end: u32, value: f32) -> Result<()> {
        let c = &mut self.core;
        c.sink.check()?;
        c.check_item(chrom, start, end, self.last)?;
        if start == end || !value.is_finite() {
            return Err(invalid(
                "bigWig item",
                "items must be non-empty with finite values",
            ));
        }
        self.last = Some((chrom, end));
        c.room(chrom)?;
        c.begin(chrom, start);
        let mut b = Bytes(std::mem::take(&mut c.block.raw));
        b.u32(start);
        b.u32(end);
        b.f32(value);
        c.block.raw = b.0;
        c.block.extend_end((chrom, end));
        c.block.items += 1;
        c.items += 1;
        if c.sample.0 < ZOOM_SAMPLE as u64 {
            c.sample.0 += 1;
            c.sample.1 += (end - start) as u64;
        }
        c.summarize(chrom, start, end, value as f64)
    }

    /// Writes the indexes, zoom levels and header; returns the output and
    /// the sidecar of the file written.
    pub fn finish(self) -> Result<(W, BbiSidecar)> {
        self.core.finish()
    }
}

/// Writes a bigBed file; see [`WriteOptions`].
pub struct BigBedWriter<W: Write + Seek> {
    core: Core<W>,
    last: Option<Point>,
    /// Coverage sweep of the current chromosome: ends of the features
    /// covering the sweep position, and that position.
    ends: BinaryHeap<Reverse<u32>>,
    at: (ChromId, u32),
}

impl<W: Write + Seek> BigBedWriter<W> {
    /// Writes the head of a file to `out`, which must start empty at
    /// position 0. Blocks are compressed on the caller's thread.
    pub fn create(
        out: W,
        chroms: Vec<Chrom>,
        schema: BedSchema,
        options: WriteOptions,
    ) -> Result<BigBedWriter<W>> {
        BigBedWriter::create_with(out, chroms, schema, options, Exec::Inline)
    }

    /// [`BigBedWriter::create`], compressing blocks on `exec`, as
    /// [`BigWigWriter::create_with`] does.
    pub fn create_with(
        out: W,
        chroms: Vec<Chrom>,
        schema: BedSchema,
        options: WriteOptions,
        exec: Exec<'_>,
    ) -> Result<BigBedWriter<W>> {
        Ok(BigBedWriter {
            core: Core::create(out, Kind::BigBed, chroms, Some(schema), options, exec)?,
            last: None,
            ends: BinaryHeap::new(),
            at: (0, 0),
        })
    }

    /// Adds a feature over `[start, end)` of chromosome `chrom`; `rest` is
    /// its fields after the third, tab-separated, without a NUL. Features
    /// come sorted by chromosome id and start.
    pub fn add(&mut self, chrom: ChromId, start: u32, end: u32, rest: &[u8]) -> Result<()> {
        self.core.sink.check()?;
        self.core.check_item(chrom, start, end, self.last)?;
        if rest.contains(&0) {
            return Err(invalid("bigBed item", "fields contain NUL"));
        }
        self.last = Some((chrom, start));
        self.sweep(chrom, start)?;
        if start < end {
            self.ends.push(Reverse(end));
        }
        let c = &mut self.core;
        c.room(chrom)?;
        c.begin(chrom, start);
        let mut b = Bytes(std::mem::take(&mut c.block.raw));
        b.u32(chrom);
        b.u32(start);
        b.u32(end);
        b.0.extend_from_slice(rest);
        b.u8(0);
        c.block.raw = b.0;
        c.block.extend_end((chrom, end));
        c.block.items += 1;
        c.items += 1;
        if c.sample.0 < ZOOM_SAMPLE as u64 {
            c.sample.0 += 1;
            c.sample.1 += (end - start) as u64;
        }
        Ok(())
    }

    /// Emits the coverage of the current chromosome up to `(chrom, pos)`.
    fn sweep(&mut self, chrom: ChromId, pos: u32) -> Result<()> {
        if chrom != self.at.0 {
            self.sweep_to(u32::MAX)?;
            self.at = (chrom, 0);
        }
        self.sweep_to(pos)
    }

    fn sweep_to(&mut self, pos: u32) -> Result<()> {
        let chrom = self.at.0;
        loop {
            let depth = self.ends.len();
            let next = match self.ends.peek() {
                Some(&Reverse(e)) if e <= pos => e,
                _ => break,
            };
            if self.at.1 < next {
                self.core.summarize(chrom, self.at.1, next, depth as f64)?;
                self.at.1 = next;
            }
            self.ends.pop();
        }
        if !self.ends.is_empty() && self.at.1 < pos {
            let depth = self.ends.len();
            self.core.summarize(chrom, self.at.1, pos, depth as f64)?;
        }
        if pos != u32::MAX {
            self.at.1 = self.at.1.max(pos);
        }
        Ok(())
    }

    /// Writes the indexes, zoom levels and header; returns the output and
    /// the sidecar of the file written.
    pub fn finish(mut self) -> Result<(W, BbiSidecar)> {
        self.sweep_to(u32::MAX)?;
        self.core.finish()
    }
}
