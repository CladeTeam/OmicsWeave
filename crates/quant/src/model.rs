//! Products as data, independent of the files they are written to:
//! coverage [`Track`]s, [`CountMatrix`] rows over a [`FeatureTable`] and a
//! [`BarcodeTable`], ATAC [`FragmentTable`]s and per-barcode
//! [`BarcodeMetrics`]. Large products are streams read once, in order, by a
//! sink (their records may live in spill files); small ones are plain
//! tables.
//!
//! Methylation, contact and allele products add per-cytosine
//! [`Cytosine`] streams, sorted [`Contact`] streams and per-variant
//! [`VariantTable`]s; their region, gene and variant matrices are
//! [`CountMatrix`]es with named layers.
//!
//! Splicing and labeling products fit the same shapes: a junction × row
//! matrix is a [`CountMatrix`] whose [`FeatureTable`] rows are junctions;
//! spliced / unspliced / ambiguous and labeled counts are layers of a gene
//! matrix; the junctions themselves are a [`JunctionTable`]; per-read
//! conversions are a [`Histogram`] and mismatches by cycle a
//! [`MismatchTable`].

use omicsweave_convert::stage::Merge;
use omicsweave_core::{Error, Result};
use smallvec::SmallVec;

use crate::accum::counts::{Entry, Fold};
use crate::spool::{Fixed, Spooled};

/// The features a matrix counts over: one row per column of the matrix.
#[derive(Debug, Clone, Default, PartialEq)]
pub struct FeatureTable {
    /// The 10x feature type (`Gene Expression`, `Peaks`, ...).
    pub kind: String,
    pub id: Vec<String>,
    pub name: Vec<String>,
    pub chrom: Vec<String>,
    /// 0-based, half-open.
    pub start: Vec<u32>,
    pub end: Vec<u32>,
    /// `+`, `-` or `.`.
    pub strand: Vec<char>,
    /// Bases counted over (exonic length for genes).
    pub length: Vec<u64>,
}

impl FeatureTable {
    pub fn new(kind: &str) -> FeatureTable {
        FeatureTable {
            kind: kind.to_owned(),
            ..FeatureTable::default()
        }
    }

    #[allow(clippy::too_many_arguments)]
    pub fn push(
        &mut self,
        id: &str,
        name: &str,
        chrom: &str,
        start: u32,
        end: u32,
        strand: char,
        length: u64,
    ) {
        self.id.push(id.to_owned());
        self.name.push(name.to_owned());
        self.chrom.push(chrom.to_owned());
        self.start.push(start);
        self.end.push(end);
        self.strand.push(strand);
        self.length.push(length);
    }

    pub fn len(&self) -> usize {
        self.id.len()
    }

    pub fn is_empty(&self) -> bool {
        self.id.is_empty()
    }

    pub fn heap_bytes(&self) -> u64 {
        let s: usize = self
            .id
            .iter()
            .chain(&self.name)
            .chain(&self.chrom)
            .map(|x| x.capacity() + 24)
            .sum();
        (s + self.len() * 20) as u64
    }
}

/// The columns of a matrix as `features.tsv.gz` lists them. A
/// [`FeatureTable`] holds them as strings; long lists (every site of a
/// variant panel) produce their lines on demand.
pub trait FeatureList: Send + Sync {
    fn len(&self) -> usize;
    fn is_empty(&self) -> bool {
        self.len() == 0
    }
    /// Appends feature `i`'s id.
    fn id(&self, i: usize, out: &mut Vec<u8>);
    /// Appends feature `i`'s line: `id`, `name`, `type`, `chrom`, `start`,
    /// `end` (0-based half-open), tab-separated, without the newline.
    fn line(&self, i: usize, out: &mut Vec<u8>);
    /// The features as a table, when they are one (featureCounts' table
    /// needs their coordinates).
    fn table(&self) -> Option<&FeatureTable> {
        None
    }
}

impl FeatureList for FeatureTable {
    fn len(&self) -> usize {
        self.id.len()
    }

    fn id(&self, i: usize, out: &mut Vec<u8>) {
        out.extend_from_slice(self.id[i].as_bytes());
    }

    fn line(&self, i: usize, out: &mut Vec<u8>) {
        out.extend_from_slice(self.id[i].as_bytes());
        out.push(b'\t');
        out.extend_from_slice(self.name[i].as_bytes());
        out.push(b'\t');
        out.extend_from_slice(self.kind.as_bytes());
        out.push(b'\t');
        out.extend_from_slice(self.chrom[i].as_bytes());
        out.push(b'\t');
        crate::sink::push_u64(out, self.start[i] as u64);
        out.push(b'\t');
        crate::sink::push_u64(out, self.end[i] as u64);
    }

    fn table(&self) -> Option<&FeatureTable> {
        Some(self)
    }
}

/// The rows of a matrix: barcodes, read groups or input files, in row
/// order.
#[derive(Debug, Clone, Default, PartialEq, Eq)]
pub struct BarcodeTable {
    pub names: Vec<Vec<u8>>,
}

impl BarcodeTable {
    pub fn len(&self) -> usize {
        self.names.len()
    }

    pub fn is_empty(&self) -> bool {
        self.names.is_empty()
    }
}

/// A run of constant coverage on one reference: `value` over
/// `[start, end)`, before scaling.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Run {
    pub tid: u32,
    pub start: u32,
    pub end: u32,
    pub value: u32,
}

impl Fixed for Run {
    const SIZE: usize = 16;
    fn put(&self, o: &mut [u8]) {
        o[0..4].copy_from_slice(&self.tid.to_le_bytes());
        o[4..8].copy_from_slice(&self.start.to_le_bytes());
        o[8..12].copy_from_slice(&self.end.to_le_bytes());
        o[12..16].copy_from_slice(&self.value.to_le_bytes());
    }
    fn get(b: &[u8]) -> Run {
        let u = |i: usize| u32::from_le_bytes(b[i..i + 4].try_into().unwrap());
        Run {
            tid: u(0),
            start: u(4),
            end: u(8),
            value: u(12),
        }
    }
}

/// Coverage of one strand (or both) over the references of a BAM: maximal
/// runs of equal nonzero value, by reference and position.
pub struct Track {
    /// The references, by tid.
    pub chroms: Vec<(String, u32)>,
    /// `None` for both strands, `Some(false)` forward, `Some(true)` reverse.
    pub strand: Option<bool>,
    pub runs: Spooled<Run>,
    /// Multiplies every value when written.
    pub scale: f64,
}

/// A count matrix, read once row by row: rows × features, each nonzero cell
/// with one value per layer.
pub struct CountMatrix {
    pub rows: BarcodeTable,
    pub features: std::sync::Arc<dyn FeatureList>,
    /// Layer names; the first is the main matrix.
    pub layers: Vec<String>,
    /// Every layer's file is named after it (`methylated.cpk`,
    /// `total.cpk`); otherwise the first is `matrix.*`.
    pub named_layers: bool,
    pub(crate) cells: Cells,
}

/// The sorted cell entries of a matrix and how to fold them.
pub(crate) struct Cells {
    pub merge: Option<Merge<Entry>>,
    pub fold: Fold,
    pub head: Option<Entry>,
}

/// One row of a [`CountMatrix`]: `(column, values)` with ascending
/// columns, `values` holding one count per layer.
#[derive(Debug, Clone, Default, PartialEq, Eq)]
pub struct MatrixRow {
    pub row: u32,
    pub cols: Vec<u32>,
    /// `values[layer * cols.len() + i]`.
    pub values: Vec<u32>,
}

impl CountMatrix {
    /// The next row with any nonzero cell, in row order (rows without one
    /// are skipped).
    pub fn next_row(&mut self, out: &mut MatrixRow) -> Result<bool> {
        let layers = self.layers.len();
        self.cells.next_row(layers, out)
    }

    /// Every cell, for tests and small matrices: `(row, col, values)`.
    pub fn collect(mut self) -> Result<Vec<(u32, u32, Vec<u32>)>> {
        let mut out = Vec::new();
        let mut r = MatrixRow::default();
        let n = self.layers.len();
        while self.next_row(&mut r)? {
            for (i, &c) in r.cols.iter().enumerate() {
                let v = (0..n).map(|l| r.values[l * r.cols.len() + i]).collect();
                out.push((r.row, c, v));
            }
        }
        Ok(out)
    }
}

impl Cells {
    fn pull(&mut self) -> Result<Option<Entry>> {
        if let Some(e) = self.head.take() {
            return Ok(Some(e));
        }
        match &mut self.merge {
            Some(m) => m.next_record(),
            None => Ok(None),
        }
    }

    fn next_row(&mut self, layers: usize, out: &mut MatrixRow) -> Result<bool> {
        out.cols.clear();
        out.values.clear();
        if layers == 0 {
            return Err(Error::invalid(
                "matrix",
                "a count matrix needs at least one layer",
            ));
        }
        let Some(first) = self.pull()? else {
            return Ok(false);
        };
        out.row = first.row;
        // Per column of the row: one value per layer.
        let mut cell: SmallVec<[u32; 4]> = smallvec::smallvec![0; layers];
        let mut col = first.col;
        let mut group = (first.layer, first.sub);
        let mut acc = first.value;
        let mut flat = Vec::new();
        let fold = self.fold;
        let close_sub = |cell: &mut [u32], layer: u16, acc: u32| {
            if let Some((l, n)) = fold.close(layer, acc) {
                let l = l.min(layers - 1);
                cell[l] = cell[l].saturating_add(n);
            }
        };
        loop {
            let next = self.pull()?;
            match next {
                Some(e) if e.row == out.row && e.col == col && (e.layer, e.sub) == group => {
                    acc = fold.combine(e.layer, acc, e.value);
                }
                Some(e) if e.row == out.row && e.col == col => {
                    close_sub(&mut cell, group.0, acc);
                    group = (e.layer, e.sub);
                    acc = e.value;
                }
                other => {
                    close_sub(&mut cell, group.0, acc);
                    if cell.iter().any(|&v| v != 0) {
                        out.cols.push(col);
                        flat.extend_from_slice(&cell);
                    }
                    cell.fill(0);
                    match other {
                        Some(e) if e.row == out.row => {
                            col = e.col;
                            group = (e.layer, e.sub);
                            acc = e.value;
                        }
                        Some(e) => {
                            self.head = Some(e);
                            break;
                        }
                        None => break,
                    }
                }
            }
        }
        for l in 0..layers {
            out.values.extend(flat.chunks_exact(layers).map(|v| v[l]));
        }
        Ok(true)
    }
}

/// An ATAC fragment: `[start, end)` of reference `tid` (Tn5-corrected),
/// barcode row `row`, seen `count` times (duplicates collapsed).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Frag {
    pub tid: u32,
    pub start: u32,
    pub end: u32,
    pub row: u32,
    pub count: u32,
}

impl Fixed for Frag {
    const SIZE: usize = 20;
    fn put(&self, o: &mut [u8]) {
        for (i, v) in [self.tid, self.start, self.end, self.row, self.count]
            .into_iter()
            .enumerate()
        {
            o[4 * i..4 * i + 4].copy_from_slice(&v.to_le_bytes());
        }
    }
    fn get(b: &[u8]) -> Frag {
        let u = |i: usize| u32::from_le_bytes(b[4 * i..4 * i + 4].try_into().unwrap());
        Frag {
            tid: u(0),
            start: u(1),
            end: u(2),
            row: u(3),
            count: u(4),
        }
    }
}

/// The fragments of a single-cell ATAC library, sorted by (reference,
/// start, end, barcode row).
pub struct FragmentTable {
    pub chroms: Vec<(String, u32)>,
    pub barcodes: BarcodeTable,
    pub fragments: Spooled<Frag>,
}

/// Per-barcode statistics: one column of counts per name, one value per
/// barcode row.
#[derive(Debug, Clone, Default, PartialEq, Eq)]
pub struct BarcodeMetrics {
    pub barcodes: BarcodeTable,
    pub columns: Vec<(String, Vec<u64>)>,
}

/// Counts over a two-dimensional integer domain (for example reads by
/// number of T and number of T>C conversions).
#[derive(Debug, Clone, Default, PartialEq, Eq)]
pub struct Histogram {
    pub axes: [String; 2],
    /// The nonzero cells `(x, y, count)`, ascending.
    pub cells: Vec<(u32, u32, u64)>,
}

/// Aligned bases by sequencing cycle: for each mate (read 1 or unpaired,
/// read 2) and cycle (position from the read's 5' end as sequenced, hard
/// clips included), counts by transcript-sense reference base and read
/// base (`A`, `C`, `G`, `T`).
#[derive(Debug, Clone, Default, PartialEq, Eq)]
pub struct MismatchTable {
    /// `cycles[mate][cycle][reference][read]`.
    pub cycles: [Vec<[[u64; 4]; 4]>; 2],
}

impl MismatchTable {
    /// Counts over every mate and cycle: `[reference][read]`.
    pub fn totals(&self) -> [[u64; 4]; 4] {
        let mut t = [[0u64; 4]; 4];
        for c in self.cycles.iter().flatten() {
            for r in 0..4 {
                for q in 0..4 {
                    t[r][q] += c[r][q];
                }
            }
        }
        t
    }
}

/// A splice junction: an intron `[start, end)` of reference `tid`, with
/// STAR's `SJ.out.tab` columns.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Junction {
    pub tid: u32,
    pub start: u32,
    pub end: u32,
    /// 0 undefined, 1 `+`, 2 `-`.
    pub strand: u8,
    /// 0 non-canonical, 1 `GT/AG`, 2 `CT/AC`, 3 `GC/AG`, 4 `CT/GC`,
    /// 5 `AT/AC`, 6 `GT/AT`.
    pub motif: u8,
    pub annotated: bool,
    /// Reads (pairs) with `NH` 1 crossing it, and with `NH` > 1.
    pub unique: u32,
    pub multi: u32,
    /// The longest overhang of any read crossing it: the shorter of the
    /// two aligned segments next to the intron.
    pub max_overhang: u32,
}

/// Every junction seen, sorted by (reference, start, end).
#[derive(Debug, Clone, Default, PartialEq, Eq)]
pub struct JunctionTable {
    pub chroms: Vec<(String, u32)>,
    pub junctions: Vec<Junction>,
}

/// Calls at one cytosine: reference `tid`, position (the `C`, or the `G`
/// of a reverse-strand cytosine), context and strand (`code`: context
/// index 0 CpG, 1 CHG, 2 CHH, bit 2 set on the reverse strand), and the
/// methylated and unmethylated calls there.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Cytosine {
    pub tid: u32,
    pub pos: u32,
    pub code: u8,
    pub meth: u32,
    pub unmeth: u32,
}

impl Cytosine {
    pub fn reverse(&self) -> bool {
        self.code & 4 != 0
    }

    pub fn context(&self) -> usize {
        (self.code & 3) as usize
    }
}

impl Fixed for Cytosine {
    const SIZE: usize = 17;
    fn put(&self, o: &mut [u8]) {
        o[0..4].copy_from_slice(&self.tid.to_le_bytes());
        o[4..8].copy_from_slice(&self.pos.to_le_bytes());
        o[8] = self.code;
        o[9..13].copy_from_slice(&self.meth.to_le_bytes());
        o[13..17].copy_from_slice(&self.unmeth.to_le_bytes());
    }
    fn get(b: &[u8]) -> Cytosine {
        let u = |i: usize| u32::from_le_bytes(b[i..i + 4].try_into().unwrap());
        Cytosine {
            tid: u(0),
            pos: u(4),
            code: b[8],
            meth: u(9),
            unmeth: u(13),
        }
    }
}

/// Methylation calls by context, mate (read 1 or unpaired, read 2) and
/// position in the read as sequenced (0-based here, from the 5' end, hard
/// clips included): `[methylated, unmethylated]`.
#[derive(Debug, Clone, Default, PartialEq, Eq)]
pub struct Mbias {
    /// `counts[mate * 3 + context][position]`.
    pub counts: [Vec<[u64; 2]>; 6],
}

impl Mbias {
    #[inline]
    pub fn add(&mut self, mate: usize, context: usize, position: usize, meth: bool) {
        let v = &mut self.counts[mate * 3 + context];
        if v.len() <= position {
            v.resize(position + 1, [0; 2]);
        }
        v[position][!meth as usize] += 1;
    }

    pub fn merge(&mut self, o: &Mbias) {
        for (mine, theirs) in self.counts.iter_mut().zip(&o.counts) {
            if mine.len() < theirs.len() {
                mine.resize(theirs.len(), [0; 2]);
            }
            for (a, b) in mine.iter_mut().zip(theirs) {
                a[0] += b[0];
                a[1] += b[1];
            }
        }
    }

    /// `{context: {read1|read2: {methylated: [...], unmethylated: [...]}}}`
    /// by 1-based position, for contexts and mates with calls.
    pub fn json(&self) -> crate::sink::json::Json {
        use crate::sink::json::Json;
        let mut out = Vec::new();
        for (c, name) in crate::accum::methylation::CONTEXTS.iter().enumerate() {
            let mut mates = Vec::new();
            for (m, mate) in ["read1", "read2"].iter().enumerate() {
                let v = &self.counts[m * 3 + c];
                if v.is_empty() {
                    continue;
                }
                let col = |k: usize| Json::Arr(v.iter().map(|x| Json::from(x[k])).collect());
                mates.push((
                    *mate,
                    Json::obj([("methylated", col(0)), ("unmethylated", col(1))]),
                ));
            }
            if !mates.is_empty() {
                out.push((*name, Json::obj(mates)));
            }
        }
        Json::obj(out)
    }
}

/// A Hi-C contact: the two 5' ends of a read pair, ordered so that
/// `(c1, p1) <= (c2, p2)` by reference index then position (1-based, as
/// 4DN pairs write them), with strands (`true` reverse), MAPQs (`255`
/// unknown), the read name and the row of its barcode (`u32::MAX`
/// without one). `input` and `ordinal` order equal contacts by the record
/// they came from.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Contact {
    pub c1: u32,
    pub p1: u32,
    pub c2: u32,
    pub p2: u32,
    pub s1: bool,
    pub s2: bool,
    pub mapq1: u8,
    pub mapq2: u8,
    pub row: u32,
    pub input: u32,
    pub ordinal: u64,
    pub name: Vec<u8>,
}

impl Contact {
    /// The fields that make two contacts duplicates.
    pub fn ends(&self) -> (u32, u32, u32, u32, bool, bool, u32) {
        (
            self.c1, self.c2, self.p1, self.p2, self.s1, self.s2, self.row,
        )
    }
}

/// The contacts of a library, sorted as 4DN pairs are (`chr1`, `chr2`,
/// `pos1`, `pos2`, in reference order), with what labels them.
pub struct ContactTable {
    pub chroms: Vec<(String, u32)>,
    /// Barcode rows (empty without a barcode tag).
    pub barcodes: BarcodeTable,
    pub contacts: Merge<Contact>,
}

/// A biallelic SNV of a site panel.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Site {
    /// The panel's contig.
    pub contig: u32,
    /// 0-based.
    pub pos: u32,
    pub reference: u8,
    pub alt: u8,
    /// The allele (0 REF, 1 ALT) of haplotype 1 of the phasing sample at a
    /// phased heterozygous genotype; `None` otherwise.
    pub hap1: Option<u8>,
}

/// Sites of a panel with their counts: per site, the folded `ref`, `alt`
/// and `other` units over all rows.
pub struct VariantTable {
    pub contigs: std::sync::Arc<Vec<String>>,
    pub sites: std::sync::Arc<Vec<Site>>,
    /// `[ref, alt, other]` per site.
    pub counts: Vec<[u64; 3]>,
}
