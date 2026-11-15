//! [`HicWriter`]: writes a `.hic` file (v8 or v9) from binned maps the
//! caller provides.
//!
//! Layout, as the Juicer tools write it: header, blocks, matrix records,
//! footer (master index, raw expected vectors), normalized expected vectors,
//! normalization vector index. No expected or normalization vectors are
//! stored: their lists are written empty, so readers find a valid file with
//! raw counts only.
//!
//! Maps are cut into blocks on the caller's thread; each block is encoded
//! and compressed as a job of an [`Exec`] and committed (written, indexed,
//! counted into its map's statistics) in the order it was added, so the
//! output does not depend on how the jobs ran.

use std::cell::RefCell;
use std::collections::{BTreeMap, HashMap, HashSet, VecDeque};
use std::io::{Seek, SeekFrom, Write};
use std::sync::Arc;

use omicsweave_core::exec::{Exec, Ordered};
use omicsweave_core::{ContactUnit, Error, Result};

use crate::block::{Deflater, Record, encode_block};
use crate::matrix::{BlockEntry, Grid};
use crate::meta::{Chrom, DEPTH_BASE_ATTRIBUTE, Depth, MAX_VERSION, MIN_VERSION};
use crate::sidecar::{BlockSize, HicSidecar, MatrixBlocks};

/// Small resolution lists fit in a few comparisons without a hash lookup.
const LINEAR_RESOLUTIONS: usize = 8;

/// What a new file declares before any map is written.
#[derive(Debug, Clone, PartialEq)]
pub struct HicHeader {
    /// 8 or 9.
    pub version: i32,
    pub genome: String,
    /// Key-value pairs; `v9-depth-base` selects the v9 block layout.
    pub attributes: Vec<(String, String)>,
    pub chroms: Vec<Chrom>,
    pub bp_resolutions: Vec<u32>,
    pub frag_resolutions: Vec<u32>,
    /// Restriction site positions (sorted) per chromosome; one list per
    /// chromosome whenever there are fragment resolutions. Fragment `i` of a
    /// chromosome runs from site `i - 1` to site `i`.
    pub sites: Vec<Vec<u32>>,
}

impl HicHeader {
    /// A v9 header with base-pair resolutions only.
    pub fn new(genome: &str, chroms: Vec<Chrom>, bp_resolutions: Vec<u32>) -> HicHeader {
        HicHeader {
            version: MAX_VERSION,
            genome: genome.to_owned(),
            attributes: vec![("software".into(), "omicsweave".into())],
            chroms,
            bp_resolutions,
            frag_resolutions: Vec::new(),
            sites: Vec::new(),
        }
    }

    pub fn depth(&self) -> Depth {
        Depth::from_attribute(
            self.attributes
                .iter()
                .find(|(k, _)| k == DEPTH_BASE_ATTRIBUTE)
                .map(|(_, v)| v.as_str()),
        )
    }

    pub fn resolutions(&self, unit: ContactUnit) -> &[u32] {
        match unit {
            ContactUnit::Bp => &self.bp_resolutions,
            ContactUnit::Frag => &self.frag_resolutions,
        }
    }

    /// A chromosome's span in `unit`: its length, or its fragment count.
    pub fn extent(&self, chrom: u32, unit: ContactUnit) -> u64 {
        match unit {
            ContactUnit::Bp => self.chroms[chrom as usize].len,
            ContactUnit::Frag => self.sites[chrom as usize].len() as u64,
        }
    }

    /// Bins of a chromosome at `bin_size` (as the Juicer tools count them).
    pub fn bins(&self, chrom: u32, unit: ContactUnit, bin_size: u32) -> u64 {
        self.extent(chrom, unit) / bin_size as u64 + 1
    }

    pub(crate) fn check(&self) -> Result<()> {
        let bad = |detail: String| Err(Error::invalid("hic header", detail));
        if !(MIN_VERSION..=MAX_VERSION).contains(&self.version) {
            return Err(Error::Unsupported(
                "writing .hic versions other than 8 and 9",
            ));
        }
        if self.chroms.len() > i32::MAX as usize {
            return bad("too many chromosomes".into());
        }
        let mut names = HashSet::new();
        for c in &self.chroms {
            if !names.insert(&c.name) {
                return bad(format!("duplicate chromosome name {:?}", c.name));
            }
            if c.name.contains('\0') {
                return bad(format!("chromosome name {:?} contains NUL", c.name));
            }
            if (self.version < 9 && c.len > i32::MAX as u64) || c.len > i64::MAX as u64 {
                return bad(format!("{} too long for v{}", c.name, self.version));
            }
        }
        for (k, v) in &self.attributes {
            if k.contains('\0') || v.contains('\0') {
                return bad("attribute contains NUL".into());
            }
        }
        if self.genome.contains('\0') {
            return bad("genome contains NUL".into());
        }
        for unit in [ContactUnit::Bp, ContactUnit::Frag] {
            if self.resolutions(unit).len() > i32::MAX as usize {
                return bad(format!("too many {} resolutions", unit.as_str()));
            }
            let mut seen = HashSet::new();
            for &r in self.resolutions(unit) {
                if r == 0 || r > i32::MAX as u32 {
                    return bad(format!("bad bin size {r}"));
                }
                if !seen.insert(r) {
                    return bad(format!("duplicate {} resolution {r}", unit.as_str()));
                }
            }
        }
        if !self.frag_resolutions.is_empty() {
            if self.sites.len() != self.chroms.len() {
                return bad("fragment resolutions need one site list per chromosome".into());
            }
            for s in &self.sites {
                if s.windows(2).any(|w| w[0] > w[1])
                    || s.last().is_some_and(|&p| p > i32::MAX as u32)
                {
                    return bad("site positions must be sorted and below 2^31".into());
                }
            }
        }
        for unit in [ContactUnit::Bp, ContactUnit::Frag] {
            for &res in self.resolutions(unit) {
                for c in 0..self.chroms.len() as u32 {
                    if self.bins(c, unit, res) > i32::MAX as u64 + 1 {
                        return bad(format!(
                            "chromosome {c} at {res} {} exceeds signed bin coordinates",
                            unit.as_str()
                        ));
                    }
                }
            }
        }
        Ok(())
    }
}

/// Runs of equal values, `(value, count)`, sorted by [`f32::total_cmp`].
type Runs = Vec<(f32, u32)>;

/// The runs of a block's values (fewer than 2^31 of them).
fn runs_of(records: &[Record]) -> Runs {
    let mut values: Vec<f32> = records.iter().map(|r| r.value).collect();
    values.sort_unstable_by(f32::total_cmp);
    values
        .chunk_by(|a, b| a.total_cmp(b).is_eq())
        .map(|c| (c[0], c.len() as u32))
        .collect()
}

/// The nonzero values of a map, as a multiset for its percentiles.
///
/// Raw contact counts repeat heavily (mostly small integers), so values are
/// kept as runs: each block's job hands over its own, and they are merged
/// into the map's whenever the unmerged ones outnumber them. Memory follows
/// the number of distinct values, not of cells, the serial work the number
/// of runs, and the percentiles stay exact.
#[derive(Default)]
struct Values {
    fresh: Runs,
    /// Equal values may span several runs only when one count would
    /// overflow.
    runs: Runs,
    count: u64,
}

impl Values {
    const MIN_FRESH: usize = 1 << 12;

    fn add(&mut self, runs: &[(f32, u32)]) {
        self.fresh.extend_from_slice(runs);
        self.count += runs.iter().map(|&(_, n)| n as u64).sum::<u64>();
        if self.fresh.len() >= self.runs.len().max(Self::MIN_FRESH) {
            self.fold();
        }
    }

    fn fold(&mut self) {
        if self.fresh.is_empty() {
            return;
        }
        self.fresh.sort_unstable_by(|a, b| a.0.total_cmp(&b.0));
        let mut fresh = self.fresh.iter().copied().peekable();
        let mut old = self.runs.iter().copied().peekable();
        let mut runs = Vec::with_capacity(self.runs.len() + 16);
        let mut put = |v: f32, n: u32| {
            let mut n = n as u64;
            if let Some(last) = runs
                .last_mut()
                .filter(|l: &&mut (f32, u32)| l.0.total_cmp(&v).is_eq())
            {
                let room = ((u32::MAX - last.1) as u64).min(n);
                last.1 += room as u32;
                n -= room;
            }
            if n > 0 {
                runs.push((v, n as u32));
            }
        };
        loop {
            let next = match (old.peek(), fresh.peek()) {
                (Some(a), Some(b)) if a.0.total_cmp(&b.0).is_le() => old.next(),
                (Some(_), Some(_)) => fresh.next(),
                (Some(_), None) => old.next(),
                (None, Some(_)) => fresh.next(),
                (None, None) => break,
            };
            let (v, n) = next.expect("peeked");
            put(v, n);
        }
        self.runs = runs;
        self.fresh.clear();
    }

    /// The `i`-th smallest value; call after [`Values::fold`].
    fn nth(&self, i: u64) -> f32 {
        let mut seen = 0u64;
        for &(v, n) in &self.runs {
            seen += n as u64;
            if i < seen {
                return v;
            }
        }
        unreachable!("rank {i} past {seen} values")
    }

    /// Apache Commons' default percentile.
    fn percentile(&mut self, p: f64) -> f32 {
        self.fold();
        let n = self.count;
        if n == 0 {
            return 0.0;
        }
        let pos = p * (n + 1) as f64 / 100.0;
        let value = if pos < 1.0 {
            self.nth(0) as f64
        } else if pos >= n as f64 {
            self.nth(n - 1) as f64
        } else {
            let at = pos.floor() as u64;
            let (lower, upper) = (self.nth(at - 1) as f64, self.nth(at) as f64);
            lower + (pos - pos.floor()) * (upper - lower)
        };
        value as f32
    }
}

/// What a matrix record stores about its map, accumulated block by block
/// as blocks are committed (in block order, so the float sum is the same
/// for every executor).
#[derive(Default)]
struct Stats {
    sum: f64,
    cells: u64,
    values: Values,
}

impl Stats {
    fn add(&mut self, intra: bool, records: &[Record], runs: &[(f32, u32)]) {
        for r in records {
            // A stored off-diagonal cell stands for both halves.
            let weight = if intra && r.x != r.y { 2.0 } else { 1.0 };
            self.sum += weight * r.value as f64;
        }
        self.cells += records.len() as u64;
        self.values.add(runs);
    }

    /// sum, occupied cells, 5th and 95th percentiles, as stored.
    fn summary(mut self) -> [f32; 4] {
        [
            self.sum as f32,
            self.cells as f32,
            self.values.percentile(5.0),
            self.values.percentile(95.0),
        ]
    }
}

/// A stored resolution of a chromosome pair's matrix.
struct Zoom {
    unit: ContactUnit,
    index: i32,
    intra: bool,
    bin_size: u32,
    grid: Grid,
    blocks: Vec<BlockEntry>,
    /// Per block: inflated bytes and records, for the sidecar.
    sizes: Vec<(u32, u32)>,
    /// Accumulating until the map has ended and its last block is
    /// committed.
    stats: Option<Stats>,
    summary: [f32; 4],
    in_flight: usize,
    ended: bool,
}

impl Zoom {
    fn settle(&mut self) {
        if self.ended
            && self.in_flight == 0
            && let Some(stats) = self.stats.take()
        {
            self.summary = stats.summary();
        }
    }
}

/// A block's cells, merged, nonzero and sorted by (y, x), the runs of
/// their values, and its stored (compressed) bytes; empty when no cell is
/// left.
struct Encoded {
    records: Vec<Record>,
    runs: Runs,
    packed: Vec<u8>,
    inflated: u32,
}

/// A thread's encoding buffer and zlib states, one per level used.
struct Scratch {
    raw: Vec<u8>,
    deflaters: Vec<(u32, Deflater)>,
}

thread_local! {
    // zlib state (about 256 KiB at the default level) and the encoding
    // buffer outlive a block: jobs run on whichever pool thread is free and
    // carry no state, so each thread keeps its own instead of allocating
    // them per block.
    static SCRATCH: RefCell<Scratch> = const {
        RefCell::new(Scratch {
            raw: Vec::new(),
            deflaters: Vec::new(),
        })
    };
}

/// The pure half of writing a block, run as a job: checks the records
/// against the block when `check` is given, merges cells (in input order),
/// encodes and compresses. Without `check`, `write_whole` supplies records
/// already sorted by (y, x) with unique cells.
fn encode(
    version: i32,
    level: u32,
    check: Option<(MatrixLayout, u32)>,
    mut records: Vec<Record>,
) -> Result<Encoded> {
    if let Some((layout, number)) = check {
        for r in &records {
            layout.check(r, number)?;
        }
    }
    records.retain(|r| r.value != 0.0);
    if check.is_some() {
        // Stable, so records of one cell add up in the order given.
        records.sort_by_key(|r| (r.y, r.x));
        records.dedup_by(|b, a| {
            let same = a.x == b.x && a.y == b.y;
            if same {
                a.value += b.value;
            }
            same
        });
    }
    if records.iter().any(|r| !r.value.is_finite()) {
        return Err(Error::invalid("hic map", "merged count is not finite"));
    }
    records.retain(|r| r.value != 0.0);
    if records.is_empty() {
        return Ok(Encoded {
            records,
            runs: Vec::new(),
            packed: Vec::new(),
            inflated: 0,
        });
    }
    let level = level.min(9);
    let (packed, inflated) = SCRATCH.with_borrow_mut(|scratch| -> Result<_> {
        let Scratch { raw, deflaters } = scratch;
        encode_block(version, &records, raw)?;
        let i = match deflaters.iter().position(|(l, _)| *l == level) {
            Some(i) => i,
            None => {
                deflaters.push((level, Deflater::new(level)));
                deflaters.len() - 1
            }
        };
        let mut packed = Vec::new();
        deflaters[i].1.deflate(raw, &mut packed)?;
        let inflated = u32::try_from(raw.len())
            .map_err(|_| Error::invalid("hic map", "block exceeds 4 GiB"))?;
        Ok((packed, inflated))
    })?;
    Ok(Encoded {
        runs: runs_of(&records),
        records,
        packed,
        inflated,
    })
}

/// Little-endian serialization into a byte buffer.
#[derive(Default)]
struct Bytes(Vec<u8>);

impl Bytes {
    fn i32(&mut self, v: i32) {
        self.0.extend_from_slice(&v.to_le_bytes());
    }
    fn i64(&mut self, v: i64) {
        self.0.extend_from_slice(&v.to_le_bytes());
    }
    fn f32(&mut self, v: f32) {
        self.0.extend_from_slice(&v.to_le_bytes());
    }
    fn str(&mut self, s: &str) {
        self.0.extend_from_slice(s.as_bytes());
        self.0.push(0);
    }
    /// `int` before v9, `long` from v9.
    fn count(&mut self, version: i32, v: u64) {
        if version > 8 {
            self.i64(v as i64)
        } else {
            self.i32(v as i32)
        }
    }
}

/// Bins per block side the default layout gives maps within a chromosome:
/// a 131-kb window at 1 kb then touches about 2 blocks of the v9 diagonal
/// layout and decodes about 2.5 times the contacts it uses (the Juicer
/// grid's ~1000-bin blocks decode about 15 times as many).
pub const DEFAULT_INTRA_BLOCK_BINS: u32 = 128;

/// How a map is cut into blocks.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum BlockBins {
    /// The grid the Juicer tools choose: about 1000 bins per block side.
    Juicer,
    /// This many bins per block side, at every resolution. Any reader of
    /// the format takes a map's block size from the file.
    Fixed(u32),
}

/// Block sizes of the maps a writer writes: within a chromosome, and
/// between chromosomes.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct BlockLayout {
    pub intra: BlockBins,
    pub inter: BlockBins,
}

impl Default for BlockLayout {
    /// [`DEFAULT_INTRA_BLOCK_BINS`]-bin blocks within a chromosome, where
    /// loaders read windows along the diagonal; the Juicer grid between
    /// chromosomes, whose sparse maps would only fragment into tiny blocks.
    fn default() -> BlockLayout {
        BlockLayout {
            intra: BlockBins::Fixed(DEFAULT_INTRA_BLOCK_BINS),
            inter: BlockBins::Juicer,
        }
    }
}

impl BlockLayout {
    /// The Juicer tools' grid for every map.
    pub fn juicer() -> BlockLayout {
        BlockLayout {
            intra: BlockBins::Juicer,
            inter: BlockBins::Juicer,
        }
    }

    fn of(&self, intra: bool) -> BlockBins {
        if intra { self.intra } else { self.inter }
    }

    fn check(&self) -> Result<()> {
        if [self.intra, self.inter].contains(&BlockBins::Fixed(0)) {
            return Err(Error::invalid("hic block grid", "zero bins per block"));
        }
        Ok(())
    }
}

/// Where the cells of one map go: how [`HicWriter::write_matrix`] orients
/// a record and which block holds it. Callers streaming a map with
/// [`HicWriter::add_block`] use it to partition records by block.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct MatrixLayout {
    /// As stored: `chr1 <= chr2`.
    chr1: u32,
    chr2: u32,
    /// The chromosomes were given as (`chr2`, `chr1`).
    flip: bool,
    unit: ContactUnit,
    bin_size: u32,
    index: i32,
    blocks: BlockBins,
    grid: Grid,
    bins: (u64, u64),
    version: i32,
    depth: Depth,
}

impl MatrixLayout {
    /// The stored chromosome pair, `chr1 <= chr2`.
    pub fn chroms(&self) -> (u32, u32) {
        (self.chr1, self.chr2)
    }

    /// Whether the chromosomes were given in the other order, so that
    /// [`MatrixLayout::place`] swaps x and y.
    pub fn flipped(&self) -> bool {
        self.flip
    }

    pub fn unit(&self) -> ContactUnit {
        self.unit
    }

    pub fn bin_size(&self) -> u32 {
        self.bin_size
    }

    pub fn grid(&self) -> Grid {
        self.grid
    }

    /// Bins of the stored `chr1` (x) and `chr2` (y).
    pub fn bins(&self) -> (u64, u64) {
        self.bins
    }

    fn intra(&self) -> bool {
        self.chr1 == self.chr2
    }

    /// The block number and stored (x, y) of the cell at bin `x` of the
    /// first and bin `y` of the second chromosome, as given to
    /// [`HicWriter::matrix_layout`]: swapped to the stored chromosome order,
    /// then, within a chromosome, to the upper half (`x <= y`). Bins outside
    /// the map are placed too; [`HicWriter::add_block`] rejects them.
    pub fn place(&self, x: u32, y: u32) -> (u32, u32, u32) {
        let (mut x, mut y) = if self.flip { (y, x) } else { (x, y) };
        if self.intra() && x > y {
            std::mem::swap(&mut x, &mut y);
        }
        (self.number(x, y), x, y)
    }

    /// The block of a stored cell, saturated to `u32`.
    fn number(&self, x: u32, y: u32) -> u32 {
        let n = self
            .grid
            .block_number(self.version, self.depth, self.intra(), x, y);
        u32::try_from(n).unwrap_or(u32::MAX)
    }

    /// Checks a stored record of block `number`.
    fn check(&self, r: &Record, number: u32) -> Result<()> {
        if !r.value.is_finite() {
            return Err(Error::invalid("hic map", "non-finite value"));
        }
        let (n1, n2) = self.bins;
        if r.x as u64 >= n1 || r.y as u64 >= n2 {
            return Err(Error::invalid(
                "hic map",
                format!("bin ({}, {}) outside the {n1} x {n2} map", r.x, r.y),
            ));
        }
        if self.intra() && r.x > r.y {
            return Err(Error::invalid(
                "hic map",
                format!("bin ({}, {}) below the diagonal", r.x, r.y),
            ));
        }
        if self.number(r.x, r.y) != number {
            return Err(Error::invalid(
                "hic map",
                format!("bin ({}, {}) is not in block {number}", r.x, r.y),
            ));
        }
        Ok(())
    }
}

/// The layouts of a file's maps: what [`HicWriter::matrix_layout`] needs
/// of the header and block sizes, cheap to clone and shareable across
/// threads, so cells can be placed in blocks away from the writer.
#[derive(Debug, Clone)]
pub struct MatrixLayouts {
    inner: Arc<Layouts>,
}

#[derive(Debug)]
struct Layouts {
    chroms: usize,
    /// Per unit (BP, FRAG): each chromosome's span, and the resolutions.
    extents: [Vec<u64>; 2],
    resolutions: [Vec<u32>; 2],
    /// Resolution indexes, for lists too long to search.
    ids: HashMap<(ContactUnit, u32), i32>,
    version: i32,
    depth: Depth,
    blocks: BlockLayout,
}

fn unit_slot(unit: ContactUnit) -> usize {
    match unit {
        ContactUnit::Bp => 0,
        ContactUnit::Frag => 1,
    }
}

impl MatrixLayouts {
    /// The layouts of a file with `header` whose maps are cut as `blocks`
    /// says; fails where [`HicWriter::create_with`] would.
    pub fn new(header: &HicHeader, blocks: BlockLayout) -> Result<MatrixLayouts> {
        header.check()?;
        blocks.check()?;
        Ok(MatrixLayouts::derive(header, blocks))
    }

    /// From a checked header and block layout.
    fn derive(header: &HicHeader, blocks: BlockLayout) -> MatrixLayouts {
        let units = [ContactUnit::Bp, ContactUnit::Frag];
        let extents = units.map(|unit| {
            if header.resolutions(unit).is_empty() {
                Vec::new()
            } else {
                (0..header.chroms.len() as u32)
                    .map(|c| header.extent(c, unit))
                    .collect()
            }
        });
        let ids = units
            .into_iter()
            .filter(|&unit| header.resolutions(unit).len() > LINEAR_RESOLUTIONS)
            .flat_map(|unit| {
                header
                    .resolutions(unit)
                    .iter()
                    .enumerate()
                    .map(move |(index, &res)| ((unit, res), index as i32))
            })
            .collect();
        MatrixLayouts {
            inner: Arc::new(Layouts {
                chroms: header.chroms.len(),
                extents,
                resolutions: units.map(|unit| header.resolutions(unit).to_vec()),
                ids,
                version: header.version,
                depth: header.depth(),
                blocks,
            }),
        }
    }

    /// How maps are cut into blocks.
    pub fn blocks(&self) -> BlockLayout {
        self.inner.blocks
    }

    pub(crate) fn resolution_index(&self, unit: ContactUnit, bin_size: u32) -> Result<i32> {
        let l = &self.inner;
        let resolutions = &l.resolutions[unit_slot(unit)];
        let index = if resolutions.len() <= LINEAR_RESOLUTIONS {
            resolutions
                .iter()
                .position(|&r| r == bin_size)
                .map(|i| i as i32)
        } else {
            l.ids.get(&(unit, bin_size)).copied()
        };
        index.ok_or_else(|| {
            Error::invalid(
                "hic map",
                format!("{bin_size} {} is not a declared resolution", unit.as_str()),
            )
        })
    }

    /// The layout of the map of `chr1` (x) and `chr2` (y) at a resolution,
    /// cut into blocks as [`MatrixLayouts::blocks`] says.
    pub fn layout(
        &self,
        chr1: u32,
        chr2: u32,
        unit: ContactUnit,
        bin_size: u32,
    ) -> Result<MatrixLayout> {
        let blocks = self.inner.blocks.of(chr1 == chr2);
        self.layout_with_blocks(chr1, chr2, unit, bin_size, blocks)
    }

    /// [`MatrixLayouts::layout`] with this map's own block size.
    pub fn layout_with_blocks(
        &self,
        chr1: u32,
        chr2: u32,
        unit: ContactUnit,
        bin_size: u32,
        blocks: BlockBins,
    ) -> Result<MatrixLayout> {
        let l = &self.inner;
        for c in [chr1, chr2] {
            if c as usize >= l.chroms {
                return Err(Error::invalid(
                    "hic map",
                    format!("chromosome index {c} out of range"),
                ));
            }
        }
        let index = self.resolution_index(unit, bin_size)?;
        let (c1, c2, flip) = if chr1 <= chr2 {
            (chr1, chr2, false)
        } else {
            (chr2, chr1, true)
        };
        let extents = &l.extents[unit_slot(unit)];
        let (e1, e2) = (extents[c1 as usize], extents[c2 as usize]);
        let grid = match blocks {
            BlockBins::Juicer => Grid::juicer(unit, bin_size, e1.max(e2), c1 == c2),
            BlockBins::Fixed(bins) => Grid::with_block_bins(bin_size, e1.max(e2), bins)?,
        };
        // As the Juicer tools count bins.
        let bins = |extent: u64| extent / bin_size as u64 + 1;
        Ok(MatrixLayout {
            chr1: c1,
            chr2: c2,
            flip,
            unit,
            bin_size,
            index,
            blocks,
            grid,
            bins: (bins(e1), bins(e2)),
            version: l.version,
            depth: l.depth,
        })
    }
}

/// The map being streamed.
struct Open {
    layout: MatrixLayout,
    /// Its position among its chromosome pair's resolutions.
    zoom: usize,
    last: Option<u32>,
}

/// A block submitted and not yet committed.
struct Flight {
    pair: (u32, u32),
    zoom: usize,
    number: u32,
}

/// Writes a `.hic` file to `W`, which must start empty at position 0 (wrap
/// a file in a `BufWriter`). Maps may be given in any order, whole
/// ([`HicWriter::write_matrix`]) or block by block
/// ([`HicWriter::begin_matrix`], [`HicWriter::add_block`],
/// [`HicWriter::end_matrix`]); [`HicWriter::finish`] writes the indexes.
///
/// Blocks are encoded and compressed as jobs of the writer's [`Exec`] and
/// written in the order they were added, so the file and sidecar are the
/// same for every executor. An output I/O failure, a block its job
/// rejects, or a job the executor refuses (its pool was dropped) makes
/// subsequent writes and `finish` fail: a partial write can leave the
/// output position unknown, and the blocks before a lost one are already
/// written. With [`Exec::Pool`] such an error may surface on a
/// later call than the one that added the block.
pub struct HicWriter<W: Write + Seek> {
    out: W,
    header: HicHeader,
    position: u64,
    failed: bool,
    nvi_field: u64,
    level: u32,
    keyed: Vec<(u64, Record)>,
    matrices: BTreeMap<(u32, u32), Vec<Zoom>>,
    written: HashSet<(u32, u32, ContactUnit, u32)>,
    layouts: MatrixLayouts,
    open: Option<Open>,
    exec: Box<dyn Ordered<Encoded>>,
    flights: VecDeque<Flight>,
}

impl<W: Write + Seek> HicWriter<W> {
    /// Checks and writes the header. `level` is the zlib level of blocks
    /// (the Juicer tools use 6); maps are cut into blocks as
    /// [`BlockLayout::default`] says.
    pub fn create(out: W, header: HicHeader, level: u32) -> Result<HicWriter<W>> {
        HicWriter::create_with(out, header, level, BlockLayout::default())
    }

    /// [`HicWriter::create`] with the block sizes of `layout`.
    pub fn create_with(
        out: W,
        header: HicHeader,
        level: u32,
        layout: BlockLayout,
    ) -> Result<HicWriter<W>> {
        HicWriter::create_exec(out, header, level, layout, Exec::Inline)
    }

    /// [`HicWriter::create_with`], encoding and compressing blocks on
    /// `exec`; with [`Exec::Pool`], at most `limit` blocks are in flight.
    pub fn create_exec(
        mut out: W,
        header: HicHeader,
        level: u32,
        layout: BlockLayout,
        exec: Exec<'_>,
    ) -> Result<HicWriter<W>> {
        header.check()?;
        layout.check()?;
        let v = header.version;
        let mut b = Bytes::default();
        b.0.extend_from_slice(b"HIC\0");
        b.i32(v);
        b.i64(0); // footer position
        b.str(&header.genome);
        let nvi_field = b.0.len() as u64;
        if v > 8 {
            b.i64(0);
            b.i64(0);
        }
        b.i32(header.attributes.len() as i32);
        for (k, val) in &header.attributes {
            b.str(k);
            b.str(val);
        }
        b.i32(header.chroms.len() as i32);
        for c in &header.chroms {
            b.str(&c.name);
            b.count(v, c.len);
        }
        for list in [&header.bp_resolutions, &header.frag_resolutions] {
            b.i32(list.len() as i32);
            for &r in list.iter() {
                b.i32(r as i32);
            }
        }
        if !header.frag_resolutions.is_empty() {
            for sites in &header.sites {
                b.i32(sites.len() as i32);
                for &s in sites {
                    b.i32(s as i32);
                }
            }
        }
        out.write_all(&b.0)?;
        let layouts = MatrixLayouts::derive(&header, layout);
        Ok(HicWriter {
            out,
            header,
            position: b.0.len() as u64,
            failed: false,
            nvi_field,
            level,
            keyed: Vec::new(),
            matrices: BTreeMap::new(),
            written: HashSet::new(),
            layouts,
            open: None,
            exec: exec.queue(),
            flights: VecDeque::new(),
        })
    }

    /// How new maps are cut into blocks.
    pub fn layout(&self) -> BlockLayout {
        self.layouts.blocks()
    }

    /// The layouts of this writer's maps, to place cells on other threads.
    pub fn layouts(&self) -> MatrixLayouts {
        self.layouts.clone()
    }

    pub fn header(&self) -> &HicHeader {
        &self.header
    }

    fn check_writable(&self) -> Result<()> {
        if self.failed {
            Err(Error::invalid(
                "hic writer",
                "an earlier write or block failed",
            ))
        } else {
            Ok(())
        }
    }

    fn write(&mut self, bytes: &[u8]) -> Result<u64> {
        self.check_writable()?;
        let at = self.position;
        let end = at
            .checked_add(bytes.len() as u64)
            .filter(|&end| end <= i64::MAX as u64)
            .ok_or_else(|| Error::invalid("hic writer", "file position exceeds i64"))?;
        if let Err(error) = self.out.write_all(bytes) {
            self.failed = true;
            return Err(error.into());
        }
        self.position = end;
        Ok(at)
    }

    /// The layout of the map of `chr1` (x) and `chr2` (y) at a resolution,
    /// cut into blocks as the writer's [`BlockLayout`] says
    /// ([`MatrixLayouts::layout`]).
    pub fn matrix_layout(
        &self,
        chr1: u32,
        chr2: u32,
        unit: ContactUnit,
        bin_size: u32,
    ) -> Result<MatrixLayout> {
        self.layouts.layout(chr1, chr2, unit, bin_size)
    }

    /// [`HicWriter::matrix_layout`] with this map's own block size.
    pub fn matrix_layout_with_blocks(
        &self,
        chr1: u32,
        chr2: u32,
        unit: ContactUnit,
        bin_size: u32,
        blocks: BlockBins,
    ) -> Result<MatrixLayout> {
        self.layouts
            .layout_with_blocks(chr1, chr2, unit, bin_size, blocks)
    }

    /// Fails unless a map of `layout` may begin now.
    fn check_new(&self, layout: &MatrixLayout) -> Result<()> {
        self.check_writable()?;
        if self.open.is_some() {
            return Err(Error::invalid("hic map", "another map has not ended"));
        }
        let (given1, given2) = if layout.flip {
            (layout.chr2, layout.chr1)
        } else {
            (layout.chr1, layout.chr2)
        };
        let own = self.matrix_layout_with_blocks(
            given1,
            given2,
            layout.unit,
            layout.bin_size,
            layout.blocks,
        )?;
        if own != *layout {
            return Err(Error::invalid(
                "hic map",
                "layout does not belong to this writer",
            ));
        }
        let (c1, c2) = (layout.chr1, layout.chr2);
        if self
            .written
            .contains(&(c1, c2, layout.unit, layout.bin_size))
        {
            return Err(Error::invalid(
                "hic map",
                format!(
                    "{c1}_{c2} at {} {} written twice",
                    layout.bin_size,
                    layout.unit.as_str()
                ),
            ));
        }
        Ok(())
    }

    /// Writes one chromosome pair's map at one resolution, cut into blocks
    /// as the writer's [`BlockLayout`] says. `records` hold x bins on `chr1`
    /// and y bins on `chr2`, in any order and orientation (a map within a
    /// chromosome may give either half, or both: records of the same stored
    /// cell add up, in the order given). Zero values are not stored.
    /// Invalid records fail the call before anything of the map is written.
    pub fn write_matrix(
        &mut self,
        chr1: u32,
        chr2: u32,
        unit: ContactUnit,
        bin_size: u32,
        records: &[Record],
    ) -> Result<()> {
        let blocks = self.layouts.blocks().of(chr1 == chr2);
        self.write_matrix_with_blocks(chr1, chr2, unit, bin_size, records, blocks)
    }

    /// [`HicWriter::write_matrix`] with this map's own block size.
    pub fn write_matrix_with_blocks(
        &mut self,
        chr1: u32,
        chr2: u32,
        unit: ContactUnit,
        bin_size: u32,
        records: &[Record],
        blocks: BlockBins,
    ) -> Result<()> {
        self.check_writable()?;
        let layout = self.matrix_layout_with_blocks(chr1, chr2, unit, bin_size, blocks)?;
        self.check_new(&layout)?;
        let mut keyed = std::mem::take(&mut self.keyed);
        let result = self.write_whole(&layout, records, &mut keyed);
        keyed.clear();
        self.keyed = keyed;
        result
    }

    /// Validates and merges a whole map, so a bad record fails before the
    /// map begins, then streams it block by block.
    fn write_whole(
        &mut self,
        layout: &MatrixLayout,
        records: &[Record],
        keyed: &mut Vec<(u64, Record)>,
    ) -> Result<()> {
        let (n1, n2) = layout.bins;
        keyed.clear();
        keyed.reserve(records.len());
        for r in records {
            if !r.value.is_finite() {
                return Err(Error::invalid("hic map", "non-finite value"));
            }
            let (number, x, y) = layout.place(r.x, r.y);
            if x as u64 >= n1 || y as u64 >= n2 {
                return Err(Error::invalid(
                    "hic map",
                    format!("bin ({x}, {y}) outside the {n1} x {n2} map"),
                ));
            }
            if r.value != 0.0 {
                let value = r.value;
                keyed.push((number as u64, Record { x, y, value }));
            }
        }
        keyed.sort_by_key(|&(n, r)| (n, r.y, r.x));
        // Same cells are adjacent now: add them up.
        keyed.dedup_by(|b, a| {
            let same = a.1.x == b.1.x && a.1.y == b.1.y;
            if same {
                a.1.value += b.1.value;
            }
            same
        });
        if keyed.iter().any(|(_, r)| !r.value.is_finite()) {
            return Err(Error::invalid("hic map", "merged count is not finite"));
        }
        if keyed.last().is_some_and(|&(n, _)| n > i32::MAX as u64) {
            return Err(Error::invalid("hic map", "block number exceeds i32"));
        }
        self.begin_matrix(layout)?;
        for chunk in keyed.chunk_by(|a, b| a.0 == b.0) {
            let records = chunk.iter().map(|&(_, r)| r).collect();
            self.submit_block(chunk[0].0 as u32, records, false)?;
        }
        self.end_matrix()
    }

    /// Begins streaming the map of `layout` (from this writer's
    /// [`HicWriter::matrix_layout`]): its blocks follow with
    /// [`HicWriter::add_block`] and it ends with [`HicWriter::end_matrix`].
    /// One map streams at a time; the result is the file
    /// [`HicWriter::write_matrix`] writes for the same records.
    pub fn begin_matrix(&mut self, layout: &MatrixLayout) -> Result<()> {
        self.check_new(layout)?;
        let key = (layout.chr1, layout.chr2, layout.unit, layout.bin_size);
        let zooms = self.matrices.entry((layout.chr1, layout.chr2)).or_default();
        zooms.push(Zoom {
            unit: layout.unit,
            index: layout.index,
            intra: layout.intra(),
            bin_size: layout.bin_size,
            grid: layout.grid,
            blocks: Vec::new(),
            sizes: Vec::new(),
            stats: Some(Stats::default()),
            summary: [0.0; 4],
            in_flight: 0,
            ended: false,
        });
        self.open = Some(Open {
            layout: *layout,
            zoom: zooms.len() - 1,
            last: None,
        });
        self.written.insert(key);
        Ok(())
    }

    /// Adds block `number` of the map being streamed: all its records, as
    /// [`MatrixLayout::place`] orients them, in any order. Block numbers
    /// increase strictly. Records of one cell add up in the order given,
    /// zero values are not stored, and an empty block is skipped. A wrong
    /// block number or a missing map fails the call and changes nothing;
    /// records outside the block, outside the map or non-finite (also once
    /// merged) fail the writer.
    pub fn add_block(&mut self, number: u32, records: Vec<Record>) -> Result<()> {
        self.submit_block(number, records, true)
    }

    fn submit_block(&mut self, number: u32, records: Vec<Record>, check: bool) -> Result<()> {
        self.check_writable()?;
        let Some(open) = &self.open else {
            return Err(Error::invalid("hic map", "no map has begun"));
        };
        if number > i32::MAX as u32 {
            return Err(Error::invalid("hic map", "block number exceeds i32"));
        }
        if open.last.is_some_and(|last| number <= last) {
            return Err(Error::invalid(
                "hic map",
                format!("block {number} after block {}", open.last.unwrap_or(0)),
            ));
        }
        if records.is_empty() {
            return Ok(());
        }
        let layout = open.layout;
        let flight = Flight {
            pair: (layout.chr1, layout.chr2),
            zoom: open.zoom,
            number,
        };
        if let Some(open) = &mut self.open {
            open.last = Some(number);
        }
        self.zoom_mut(flight.pair, flight.zoom).in_flight += 1;
        while self.exec.full() {
            let done = self
                .next(true)?
                .ok_or_else(|| Error::invalid("hic writer", "block job lost"))?;
            self.commit(done)?;
        }
        let (version, level) = (self.header.version, self.level);
        let check = check.then_some((layout, number));
        self.exec
            .submit(Box::new(move || encode(version, level, check, records)))
            .inspect_err(|_| self.failed = true)?;
        self.flights.push_back(flight);
        while let Some(done) = self.next(false)? {
            self.commit(done)?;
        }
        Ok(())
    }

    /// Ends the map being streamed. Its blocks may still be in flight; they
    /// are written in turn.
    pub fn end_matrix(&mut self) -> Result<()> {
        self.check_writable()?;
        let Some(open) = self.open.take() else {
            return Err(Error::invalid("hic map", "no map has begun"));
        };
        let zoom = self.zoom_mut((open.layout.chr1, open.layout.chr2), open.zoom);
        zoom.ended = true;
        zoom.settle();
        Ok(())
    }

    fn zoom_mut(&mut self, pair: (u32, u32), zoom: usize) -> &mut Zoom {
        &mut self.matrices.get_mut(&pair).expect("map begun")[zoom]
    }

    /// The next block result in order; a rejected block fails the writer.
    fn next(&mut self, wait: bool) -> Result<Option<Encoded>> {
        self.exec.next(wait).inspect_err(|_| self.failed = true)
    }

    /// Writes the oldest block in flight and adds it to its map.
    fn commit(&mut self, done: Encoded) -> Result<()> {
        let flight = self.flights.pop_front().expect("results follow blocks");
        let position = if done.records.is_empty() {
            None
        } else {
            let size = u32::try_from(done.packed.len())
                .ok()
                .filter(|&s| s <= i32::MAX as u32)
                .ok_or_else(|| Error::invalid("hic map", "block exceeds 2 GiB"))
                .inspect_err(|_| self.failed = true)?;
            Some((self.write(&done.packed)?, size))
        };
        let zoom = self.zoom_mut(flight.pair, flight.zoom);
        zoom.in_flight -= 1;
        let intra = zoom.intra;
        if let Some(stats) = &mut zoom.stats {
            stats.add(intra, &done.records, &done.runs);
        }
        if let Some((position, size)) = position {
            zoom.blocks.push(BlockEntry {
                number: flight.number,
                position,
                size,
            });
            zoom.sizes.push((done.inflated, done.records.len() as u32));
        }
        zoom.settle();
        Ok(())
    }

    /// Writes the matrix records, footer and indexes, and returns the output.
    pub fn finish(self) -> Result<W> {
        Ok(self.finish_with_sidecar()?.0)
    }

    /// [`HicWriter::finish`], also returning the sidecar of the file
    /// written: every block's inflated size and record count, as
    /// [`HicSidecar::build`] would measure them.
    pub fn finish_with_sidecar(mut self) -> Result<(W, HicSidecar)> {
        self.check_writable()?;
        if self.open.is_some() {
            return Err(Error::invalid("hic writer", "a map has not ended"));
        }
        while let Some(done) = self.next(true)? {
            self.commit(done)?;
        }
        let v = self.header.version;
        let mut records = Bytes::default();
        let mut master = Vec::new();
        let base = self.position;
        for (&(c1, c2), zooms) in &self.matrices {
            let start = records.0.len();
            records.i32(c1 as i32);
            records.i32(c2 as i32);
            records.i32(zooms.len() as i32);
            for z in zooms {
                debug_assert!(z.stats.is_none(), "every map settled");
                let [sum, cells, percent5, percent95] = z.summary;
                records.str(z.unit.as_str());
                records.i32(z.index);
                records.f32(sum);
                records.f32(cells);
                records.f32(percent5);
                records.f32(percent95);
                records.i32(z.bin_size as i32);
                records.i32(z.grid.block_bins as i32);
                records.i32(z.grid.columns as i32);
                records.i32(z.blocks.len() as i32);
                for blk in &z.blocks {
                    records.i32(blk.number as i32);
                    records.i64(blk.position as i64);
                    records.i32(blk.size as i32);
                }
            }
            let size = records.0.len() - start;
            if size > i32::MAX as usize {
                return Err(Error::invalid("hic matrix", "record exceeds 2 GiB"));
            }
            master.push((c1, c2, base + start as u64, size as i32));
        }
        self.write(&records.0)?;

        let mut v5 = Bytes::default();
        v5.i32(master.len() as i32);
        for (c1, c2, position, size) in master {
            v5.str(&format!("{c1}_{c2}"));
            v5.i64(position as i64);
            v5.i32(size);
        }
        v5.i32(0); // raw expected vectors
        let mut footer = Bytes::default();
        footer.count(v, v5.0.len() as u64);
        footer.0.extend_from_slice(&v5.0);
        footer.i32(0); // normalized expected vectors
        let footer_position = self.write(&footer.0)?;

        let mut nvi = Bytes::default();
        nvi.i32(0); // normalization vectors
        let nvi_position = self.write(&nvi.0)?;

        self.out.seek(SeekFrom::Start(8))?;
        self.out.write_all(&footer_position.to_le_bytes())?;
        if v > 8 {
            self.out.seek(SeekFrom::Start(self.nvi_field))?;
            self.out.write_all(&nvi_position.to_le_bytes())?;
            self.out.write_all(&(nvi.0.len() as u64).to_le_bytes())?;
        }
        self.out.seek(SeekFrom::Start(self.position))?;
        self.out.flush()?;
        let mut matrices: Vec<MatrixBlocks> = self
            .matrices
            .iter()
            .flat_map(|(&(chr1, chr2), zooms)| {
                zooms.iter().map(move |z| MatrixBlocks {
                    chr1,
                    chr2,
                    unit: z.unit,
                    bin_size: z.bin_size,
                    block_bin_count: z.grid.block_bins,
                    block_column_count: z.grid.columns,
                    blocks: z
                        .blocks
                        .iter()
                        .zip(&z.sizes)
                        .map(|(b, &(inflated, records))| BlockSize {
                            number: b.number,
                            position: b.position,
                            size: b.size,
                            inflated,
                            records,
                        })
                        .collect(),
                })
            })
            .collect();
        matrices.sort_by_key(|m| (m.chr1, m.chr2, m.unit, m.bin_size));
        let sidecar = HicSidecar {
            version: self.header.version,
            matrices,
        };
        Ok((self.out, sidecar))
    }
}
