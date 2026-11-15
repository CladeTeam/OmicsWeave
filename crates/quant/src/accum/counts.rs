//! Count matrices: rows (barcodes, read groups or input files) × features,
//! as one accumulator of composable dimensions.
//!
//! - **Features** (columns): any [`FeatureSpace`] — genes by exons or by
//!   span, promoter or gene-activity windows, bigBed peaks, genome bins.
//! - **Unit**: what is counted ([`CountUnit`]).
//! - **Row key**: a tag (`CB`, ...), the read group, or the input file
//!   ([`RowKey`]), optionally restricted to a whitelist.
//! - **Deduplication**: none, by UMI, or by fragment position
//!   ([`Dedup`]).
//! - **Assignment**: featureCounts' rules — a unit overlapping features
//!   of more than one column is ambiguous unless `multi_overlap`, strand
//!   checked per [`Strandedness`], multi-mappers (`NH` > 1) excluded
//!   unless `multi_mapping`.
//! - **Layers** of gene matrices, beside the counts: spliced / unspliced /
//!   ambiguous molecules ([`Counts::velocity`], see [`super::velocity`])
//!   and labeled units ([`Counts::labeled`], see [`super::labeling`]).
//!
//! Batches emit cell entries `(row, column, layer, sub, value)` — `sub` is
//! the UMI or fragment a value belongs to — collapsed per batch; the total
//! sorts them out of core ([`Sorter`]) and the matrix folds equal keys when
//! it is read ([`Fold`]), so memory is bounded by the sort buffer, not by
//! the number of nonzero cells.

use std::collections::HashMap;
use std::io::Read;
use std::path::PathBuf;
use std::sync::Arc;

use omicsweave_convert::stage::{Sorter, Spill, read_or_end};
use omicsweave_core::{Error, Result};
use omicsweave_gtf::Flank;

use super::coverage::Tally;
use super::labeling::{Conversion, Converted, convert};
use super::velocity::{self, GeneModels};
use super::{Accumulator, BatchInfo, Delivery, MemoryNeed, Needs, PlanInfo, RunCtx};
use crate::assign::{FeatureSpace, Query, SpaceId};
use crate::model::{BarcodeTable, Cells, CountMatrix, FeatureList, FeatureTable};
use crate::prelude::{
    Filter, Hits, MateKey, Mates, Pairs, ReadView, Strandedness, expects_mate, fnv, umi_key,
};
use crate::sink::json::Json;

/// What one count stands for.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum CountUnit {
    /// An alignment, assigned by its aligned blocks (featureCounts without
    /// `-p`; with [`Dedup::Umi`], a molecule as STARsolo / Cell Ranger
    /// count them).
    Read,
    /// A read pair, assigned by the aligned blocks of both mates together
    /// (featureCounts `-p --countReadPairs`); a mate whose pair cannot be
    /// completed counts alone, as do unpaired reads.
    Pair,
    /// A pair's fragment span from its leftmost mate (`[pos, pos +
    /// |TLEN|)`, Tn5-corrected with `tn5`), assigned by that span.
    Fragment { tn5: bool },
    /// Each Tn5 insertion site of a pair's corrected fragment (two per
    /// fragment), assigned on its own.
    CutSite,
}

/// What identifies a row.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum RowKey {
    /// The value of an optional field (`CB`, `CR`, ...).
    Tag([u8; 2]),
    /// The `RG` field, one row per `@RG` of the headers.
    ReadGroup,
    /// The input BAM, one row per input.
    File,
}

/// How repeated observations of one molecule are merged.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Dedup {
    None,
    /// One count per distinct UMI (exact `UB` value) per row and feature;
    /// reads without the UMI tag are not counted.
    Umi,
    /// One count per distinct fragment position (reference, start, end) per
    /// row and feature: duplicate fragments of a barcode count once.
    Position,
}

/// How cells with equal keys fold when the matrix is read.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Reduce {
    /// Values add up.
    Sum,
    /// Each distinct `sub` counts one.
    Distinct,
    /// Each distinct `sub` counts its largest value (a fragment with both
    /// cut sites in a feature counts two however often it was seen).
    MaxPerSub,
}

impl Reduce {
    #[inline]
    pub fn combine(self, a: u32, b: u32) -> u32 {
        match self {
            Reduce::Sum => a.saturating_add(b),
            Reduce::Distinct => 1,
            Reduce::MaxPerSub => a.max(b),
        }
    }

    #[inline]
    pub fn count(self, acc: u32) -> u32 {
        match self {
            Reduce::Distinct => 1,
            _ => acc,
        }
    }
}

/// Entries of this layer carry a molecule's velocity evidence bits.
pub(crate) const VELOCITY_LAYER: u16 = 1;
/// Entries of this layer mark labeled units.
pub(crate) const LABEL_LAYER: u16 = 2;

/// How the entries of a matrix fold into its layers: the counts by
/// [`Reduce`], two-bit evidence (velocity's exonic / intronic, a
/// molecule's haplotypes) OR'ed per molecule and classified into three
/// consecutive layers (first bit only, second only, both), labeled units by
/// [`Reduce`] into the labeled layer. With `direct`, an entry's layer is
/// its output layer and every layer folds by [`Reduce`].
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Fold {
    pub reduce: Reduce,
    /// The output layer of `spliced` (then `unspliced`, `ambiguous`): the
    /// classes of the evidence layer.
    pub spliced: Option<usize>,
    /// The output layer of `labeled`.
    pub labeled: Option<usize>,
    pub direct: bool,
}

impl Fold {
    pub fn plain(reduce: Reduce) -> Fold {
        Fold {
            reduce,
            spliced: None,
            labeled: None,
            direct: false,
        }
    }

    /// Entry layers are output layers, folded by `reduce`.
    pub fn direct(reduce: Reduce) -> Fold {
        Fold {
            direct: true,
            ..Fold::plain(reduce)
        }
    }

    #[inline]
    pub fn combine(&self, layer: u16, a: u32, b: u32) -> u32 {
        if layer == VELOCITY_LAYER && !self.direct {
            a | b
        } else {
            self.reduce.combine(a, b)
        }
    }

    /// The output layer and count of one folded `(layer, sub)` group.
    #[inline]
    pub fn close(&self, layer: u16, acc: u32) -> Option<(usize, u32)> {
        if self.direct {
            return Some((layer as usize, self.reduce.count(acc)));
        }
        match layer {
            VELOCITY_LAYER => Some((self.spliced? + velocity::class(acc), 1)),
            LABEL_LAYER => Some((self.labeled?, self.reduce.count(acc))),
            _ => Some((0, self.reduce.count(acc))),
        }
    }
}

/// A labeled layer: units whose conversions reach `min_conversions`.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Labeled {
    pub conversion: Conversion,
    pub min_conversions: u32,
}

impl Default for Labeled {
    /// One `T>C` or more, at base quality 20.
    fn default() -> Labeled {
        Labeled {
            conversion: Conversion::default(),
            min_conversions: 1,
        }
    }
}

/// Where the columns come from.
#[derive(Debug, Clone)]
pub enum FeatureSource {
    /// Genes by their merged exons (needs an annotation).
    GeneExons,
    /// Genes by their whole span (STARsolo `GeneFull`).
    GeneSpans,
    /// A window around each gene's TSS.
    Promoters(Flank),
    /// Each gene body extended upstream / downstream (gene activity).
    GeneActivity(Flank),
    /// Every feature of a bigBed file (peaks; BED converts to bigBed with
    /// `ow convert`).
    BigBed(PathBuf),
    /// Tiles of this many bases over every reference of the first input.
    Bins(u32),
    /// A space built by the caller.
    Space(Arc<FeatureSpace>),
}

/// A count matrix product; see the module documentation. The presets fill
/// in what a library type usually needs.
#[derive(Debug, Clone)]
pub struct Counts {
    pub name: String,
    pub features: FeatureSource,
    /// Only genes of these biotypes (`gene_type` / `gene_biotype`).
    pub biotypes: Option<Vec<String>>,
    pub unit: CountUnit,
    pub rows: RowKey,
    /// Only these row keys, one per line (plain or gzip), rows in this
    /// order.
    pub whitelist: Option<PathBuf>,
    pub umi_tag: [u8; 2],
    pub dedup: Dedup,
    pub strand: Strandedness,
    pub multi_overlap: bool,
    pub filter: Filter,
    /// Fragments longer than this are not counted (0: no limit).
    pub max_fragment: u32,
    /// Also write a featureCounts-style `counts.tsv` (features × rows).
    pub table: bool,
    /// Add the spliced, unspliced and ambiguous layers (gene features,
    /// reads or pairs; see [`super::velocity`]).
    pub velocity: bool,
    /// Add a labeled layer: units with at least this many conversions
    /// (see [`super::labeling`]). A unit's transcript strand is the read's
    /// under a stranded library, else that of the gene it is counted in.
    /// A pair's conversions are those of both mates (a position counted
    /// once); a molecule is labeled when any of its reads is.
    pub labeled: Option<Labeled>,
}

impl Counts {
    /// Single-cell RNA (STARsolo `Gene`, Cell Ranger): genes × UMIs × `CB`,
    /// the raw matrix over every barcode, `CB` / `UB` used as given,
    /// unique alignments only, forward-stranded (10x 3' and 5').
    pub fn rna_single_cell() -> Counts {
        Counts {
            name: "genes".into(),
            features: FeatureSource::GeneExons,
            biotypes: None,
            unit: CountUnit::Read,
            rows: RowKey::Tag(*b"CB"),
            whitelist: None,
            umi_tag: *b"UB",
            dedup: Dedup::Umi,
            strand: Strandedness::Forward,
            multi_overlap: false,
            filter: Filter {
                unique: true,
                ..Filter::new()
            },
            max_fragment: 0,
            table: false,
            velocity: false,
            labeled: None,
        }
    }

    /// Bulk RNA (featureCounts defaults): genes × reads (pairs with
    /// [`CountUnit::Pair`]) × input file, unstranded, unique alignments
    /// only, every mapped record counted (duplicates, secondary and
    /// supplementary alignments included; multi-mappers are told by `NH`).
    pub fn bulk_rna() -> Counts {
        Counts {
            name: "genes".into(),
            features: FeatureSource::GeneExons,
            biotypes: None,
            unit: CountUnit::Read,
            rows: RowKey::File,
            whitelist: None,
            umi_tag: *b"UB",
            dedup: Dedup::None,
            strand: Strandedness::Unstranded,
            multi_overlap: false,
            filter: Filter {
                require: 0,
                exclude: omicsweave_bam::flag::UNMAPPED,
                min_mapq: 0,
                unique: true,
            },
            max_fragment: 0,
            table: true,
            velocity: false,
            labeled: None,
        }
    }

    /// Single-cell ATAC: `features` × Tn5 cut sites × `CB`, duplicate
    /// fragments counted once.
    pub fn atac(name: &str, features: FeatureSource, filter: Filter) -> Counts {
        Counts {
            name: name.into(),
            features,
            biotypes: None,
            unit: CountUnit::CutSite,
            rows: RowKey::Tag(*b"CB"),
            whitelist: None,
            umi_tag: *b"UB",
            dedup: Dedup::Position,
            strand: Strandedness::Unstranded,
            multi_overlap: true,
            filter,
            max_fragment: 0,
            table: false,
            velocity: false,
            labeled: None,
        }
    }

    /// ChIP-seq / CUT&Tag bulk: `features` (bins or peaks) × fragments ×
    /// input file, duplicates excluded by flag, a fragment counted in every
    /// feature it overlaps.
    pub fn chip(features: FeatureSource) -> Counts {
        Counts {
            name: "chip".into(),
            features,
            biotypes: None,
            unit: CountUnit::Fragment { tn5: false },
            rows: RowKey::File,
            whitelist: None,
            umi_tag: *b"UB",
            dedup: Dedup::None,
            strand: Strandedness::Unstranded,
            multi_overlap: true,
            filter: Filter::new(),
            max_fragment: 0,
            table: true,
            velocity: false,
            labeled: None,
        }
    }
}

/// One cell observation. Sorted by `(row, col, layer, sub)`.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Entry {
    pub row: u32,
    pub col: u32,
    pub layer: u16,
    pub sub: u64,
    pub value: u32,
}

impl Entry {
    #[inline]
    fn key(&self) -> (u32, u32, u16, u64) {
        (self.row, self.col, self.layer, self.sub)
    }
}

impl Spill for Entry {
    type Key = (u32, u32, u16, u64);
    fn key(&self) -> Self::Key {
        Entry::key(self)
    }
    fn mem_size(&self) -> usize {
        0
    }
    fn encode(&self, out: &mut Vec<u8>) {
        out.extend_from_slice(&self.row.to_le_bytes());
        out.extend_from_slice(&self.col.to_le_bytes());
        out.extend_from_slice(&self.layer.to_le_bytes());
        out.extend_from_slice(&self.sub.to_le_bytes());
        out.extend_from_slice(&self.value.to_le_bytes());
    }
    fn decode(input: &mut impl Read) -> Result<Option<Entry>> {
        let mut b = [0u8; 22];
        if !read_or_end(input, &mut b)? {
            return Ok(None);
        }
        let u = |i: usize| u32::from_le_bytes(b[i..i + 4].try_into().unwrap());
        Ok(Some(Entry {
            row: u(0),
            col: u(4),
            layer: u16::from_le_bytes([b[8], b[9]]),
            sub: u64::from_le_bytes(b[10..18].try_into().unwrap()),
            value: u(18),
        }))
    }
}

/// Sorts `entries` and folds equal keys.
pub(crate) fn collapse(entries: &mut Vec<Entry>, fold: Fold) {
    entries.sort_unstable_by_key(Entry::key);
    let mut w = 0;
    for i in 0..entries.len() {
        if w > 0 && entries[w - 1].key() == entries[i].key() {
            entries[w - 1].value =
                fold.combine(entries[i].layer, entries[w - 1].value, entries[i].value);
        } else {
            entries[w] = entries[i];
            w += 1;
        }
    }
    entries.truncate(w);
}

/// How the rows of a stage are numbered.
#[derive(Debug)]
pub(crate) enum Rows {
    /// Tag values in order of first appearance.
    Tag([u8; 2]),
    /// Tag values of a whitelist, in its order.
    Listed([u8; 2], Arc<HashMap<Vec<u8>, u32>>, Arc<Vec<Vec<u8>>>),
    ReadGroup(Arc<HashMap<Vec<u8>, u32>>, Arc<Vec<Vec<u8>>>),
    File(Arc<Vec<Vec<u8>>>),
}

/// Row keys first seen in one batch, numbered locally.
#[derive(Default)]
pub(crate) struct LocalRows {
    ids: HashMap<Vec<u8>, u32>,
    names: Vec<Vec<u8>>,
}

impl LocalRows {
    #[inline]
    pub fn id(&mut self, key: &[u8]) -> u32 {
        if let Some(&i) = self.ids.get(key) {
            return i;
        }
        let i = self.names.len() as u32;
        self.ids.insert(key.to_vec(), i);
        self.names.push(key.to_vec());
        i
    }
}

/// Row keys of the whole run, numbered in order of first appearance.
#[derive(Default)]
pub(crate) struct GlobalRows {
    ids: HashMap<Vec<u8>, u32>,
    pub names: Vec<Vec<u8>>,
    bytes: u64,
}

impl GlobalRows {
    /// The global numbers of a batch's local rows.
    pub fn remap(&mut self, local: LocalRows) -> Vec<u32> {
        local
            .names
            .into_iter()
            .map(|name| {
                let n = self.names.len() as u32;
                *self.ids.entry(name).or_insert_with_key(|k| {
                    self.bytes += 2 * k.len() as u64 + 48;
                    self.names.push(k.clone());
                    n
                })
            })
            .collect()
    }

    pub fn bytes(&self) -> u64 {
        self.bytes
    }
}

impl Rows {
    /// Whether row numbers are local to a batch until merged.
    pub(crate) fn local(&self) -> bool {
        matches!(self, Rows::Tag(_))
    }

    #[inline]
    pub(crate) fn of(&self, r: &ReadView<'_>, local: &mut LocalRows) -> Option<u32> {
        match self {
            Rows::Tag(t) => r.tag_bytes(*t).map(|v| local.id(v)),
            Rows::Listed(t, ids, _) => ids.get(r.tag_bytes(*t)?).copied(),
            Rows::ReadGroup(ids, _) => ids.get(r.tag_bytes(*b"RG")?).copied(),
            Rows::File(_) => Some(r.input),
        }
    }

    pub fn fixed_names(&self) -> Option<&[Vec<u8>]> {
        match self {
            Rows::Tag(_) => None,
            Rows::Listed(_, _, n) | Rows::ReadGroup(_, n) | Rows::File(n) => Some(n),
        }
    }
}

/// Why units were or were not counted (featureCounts' summary).
#[derive(Debug, Clone, Copy, Default)]
pub struct Status {
    pub assigned: u64,
    pub no_features: u64,
    pub ambiguous: u64,
    pub no_row: u64,
    pub no_umi: u64,
    /// Mates whose pair is counted from the other mate, and pairs too long.
    pub not_leader: u64,
    pub too_long: u64,
    /// Velocity: reads over no gene span, over several, or held by no
    /// transcript of their gene.
    pub velocity_no_gene: u64,
    pub velocity_multi_gene: u64,
    pub velocity_no_model: u64,
    /// Labeling: units without a reference to compare with.
    pub no_reference: u64,
}

impl Status {
    fn add(&mut self, o: &Status) {
        self.assigned += o.assigned;
        self.no_features += o.no_features;
        self.ambiguous += o.ambiguous;
        self.no_row += o.no_row;
        self.no_umi += o.no_umi;
        self.not_leader += o.not_leader;
        self.too_long += o.too_long;
        self.velocity_no_gene += o.velocity_no_gene;
        self.velocity_multi_gene += o.velocity_multi_gene;
        self.velocity_no_model += o.velocity_no_model;
        self.no_reference += o.no_reference;
    }
}

/// A mate waiting for the other: its row, UMI, hits, and for a labeled
/// layer its conversions and transcript strand.
#[derive(Debug, Clone)]
pub(crate) struct Half {
    row: u32,
    sub: u64,
    hits: Hits,
    conv: Option<Converted>,
    rev: bool,
}

/// The velocity layers of a stage: the gene spans reads are assigned by,
/// and the gene models they are classified against.
pub struct Velocity {
    pub(crate) space: SpaceId,
    pub(crate) models: Arc<GeneModels>,
}

/// A [`Counts`] product resolved against a task: its space, rows and
/// fold.
pub struct CountStage {
    pub(crate) spec: Counts,
    pub(crate) space: SpaceId,
    pub(crate) features: Arc<FeatureTable>,
    pub(crate) rows: Rows,
    pub(crate) fold: Fold,
    pub(crate) layers: Vec<String>,
    pub(crate) velocity: Option<Velocity>,
    /// Fixed bytes: whitelist and read-group tables.
    pub(crate) table_bytes: u64,
}

pub struct CountBatch {
    entries: Vec<Entry>,
    local: LocalRows,
    pairs: Pairs<Half>,
    dangling: Vec<(MateKey, Half)>,
    status: Status,
    tally: Tally,
}

pub struct CountTotal<'p> {
    sorter: Sorter<'p, Entry>,
    rows: GlobalRows,
    mates: Mates<Half>,
    status: Status,
    tally: Tally,
    scratch: Vec<Entry>,
}

/// Visits the sorted, distinct assignment lists once, with each column's
/// number of contributing ends. Assignment produces this order directly.
fn visit_hits(mut a: &[u32], mut b: &[u32], mut visit: impl FnMut(u32, u32)) {
    while let (Some(&x), Some(&y)) = (a.first(), b.first()) {
        match x.cmp(&y) {
            std::cmp::Ordering::Less => {
                visit(x, 1);
                a = &a[1..];
            }
            std::cmp::Ordering::Greater => {
                visit(y, 1);
                b = &b[1..];
            }
            std::cmp::Ordering::Equal => {
                visit(x, 2);
                a = &a[1..];
                b = &b[1..];
            }
        }
    }
    for &c in a.iter().chain(b) {
        visit(c, 1);
    }
}

fn pair_hits(a: &[u32], b: &[u32], multi_overlap: bool) -> Hits {
    let mut hits = Hits::new();
    let mut best = 0;
    visit_hits(a, b, |col, ends| {
        if multi_overlap {
            hits.push(col);
        } else if ends >= best {
            if ends > best {
                hits.clear();
                best = ends;
            }
            hits.push(col);
        }
    });
    hits
}

impl CountStage {
    /// Counts `hits` for one unit: assigned, ambiguous or without
    /// features; an assigned unit is also marked labeled where its
    /// conversions `conv` (for a transcript on the reverse strand when
    /// `rev`) reach the threshold.
    #[allow(clippy::too_many_arguments)]
    #[inline]
    fn decide(
        &self,
        hits: &Hits,
        row: u32,
        sub: u64,
        conv: Option<&Converted>,
        rev: bool,
        out: &mut Vec<Entry>,
        st: &mut Status,
    ) {
        if hits.is_empty() {
            st.no_features += 1;
        } else if hits.len() > 1 && !self.spec.multi_overlap {
            st.ambiguous += 1;
        } else {
            st.assigned += 1;
            for &col in hits {
                out.push(Entry {
                    row,
                    col,
                    layer: 0,
                    sub,
                    value: 1,
                });
                if let (Some(l), Some(c)) = (&self.spec.labeled, conv) {
                    let reverse = match self.spec.strand {
                        Strandedness::Unstranded => match self.features.strand[col as usize] {
                            '+' => false,
                            '-' => true,
                            _ => rev,
                        },
                        _ => rev,
                    };
                    if c.count(reverse) >= l.min_conversions {
                        out.push(Entry {
                            row,
                            col,
                            layer: LABEL_LAYER,
                            sub,
                            value: 1,
                        });
                    }
                }
            }
        }
    }

    /// The velocity evidence of one read, under molecule `sub`.
    #[allow(clippy::too_many_arguments)]
    fn velocity(
        &self,
        v: &Velocity,
        r: &ReadView<'_>,
        row: u32,
        sub: u64,
        strand: Option<bool>,
        out: &mut Vec<Entry>,
        st: &mut Status,
    ) {
        let hits = r.assign(v.space, Query::Blocks, strand);
        match hits.len() {
            0 => st.velocity_no_gene += 1,
            1 => match v.models.classify(hits[0], r.blocks(false)) {
                0 => st.velocity_no_model += 1,
                bits => out.push(Entry {
                    row,
                    col: hits[0],
                    layer: VELOCITY_LAYER,
                    sub,
                    value: bits,
                }),
            },
            _ => st.velocity_multi_gene += 1,
        }
    }

    fn strand(&self, r: &ReadView<'_>) -> Option<bool> {
        match self.spec.strand {
            Strandedness::Unstranded => None,
            lib => Some(r.reverse_for(lib)),
        }
    }

    /// Counts a pair from both mates' hits. Without `multi_overlap`, as
    /// featureCounts votes: the features overlapped by the most ends (two
    /// before one) win, and the pair is ambiguous only when several tie;
    /// with it, every feature either end overlaps counts.
    fn pair(&self, a: Half, b: Half, out: &mut Vec<Entry>, st: &mut Status) {
        let hits = pair_hits(&a.hits, &b.hits, self.spec.multi_overlap);
        let conv = match (a.conv, b.conv) {
            (Some(mut x), Some(y)) => {
                x.join(&y);
                Some(x)
            }
            (x, None) => x,
            (None, y) => y,
        };
        self.decide(&hits, a.row, a.sub, conv.as_ref(), a.rev, out, st);
    }

    pub fn space(&self) -> SpaceId {
        self.space
    }
}

impl Accumulator for CountStage {
    type Batch = CountBatch;
    type Total<'p> = CountTotal<'p>;
    type Product = (CountMatrix, Status, Tally, u64);

    fn name(&self) -> &str {
        &self.spec.name
    }

    fn needs(&self) -> Needs {
        let mut n = Needs::CIGAR | Needs::ASSIGN;
        if !matches!(self.spec.rows, RowKey::File) || self.spec.filter.unique {
            n |= Needs::TAGS;
        }
        if self.spec.dedup == Dedup::Umi {
            n |= Needs::TAGS;
        }
        if self.spec.unit == CountUnit::Pair {
            n |= Needs::MATES;
        }
        if self.spec.labeled.is_some() {
            n |=
                Needs::TAGS | Needs::SEQUENCE | Needs::QUALITY | Needs::REFERENCE | Needs::VARIANTS;
        }
        n
    }

    fn estimate(&self, info: &PlanInfo) -> MemoryNeed {
        MemoryNeed {
            fixed: self.features.heap_bytes()
                + self.table_bytes
                + self.velocity.as_ref().map_or(0, |v| v.models.heap_bytes()),
            // About one entry of 32 bytes per record of a few hundred.
            per_batch: info.batch_bytes / 4,
            min_buffer: 1 << 20,
            weight: 4,
        }
    }

    fn batch(&self, _info: BatchInfo) -> CountBatch {
        CountBatch {
            entries: Vec::new(),
            local: LocalRows::default(),
            pairs: Pairs::default(),
            dangling: Vec::new(),
            status: Status::default(),
            tally: Tally::default(),
        }
    }

    fn observe(&self, b: &mut CountBatch, r: &ReadView<'_>) -> Result<()> {
        b.tally.seen += 1;
        if let Some(why) = self.spec.filter.reject(r) {
            b.tally.reject(why);
            return Ok(());
        }
        let fragment = match self.spec.unit {
            CountUnit::CutSite | CountUnit::Fragment { tn5: true } => {
                match r.atac_fragment(self.spec.max_fragment) {
                    Some(f) => Some(f),
                    None => {
                        if r.atac_fragment(0).is_some() {
                            b.status.too_long += 1;
                        } else {
                            b.status.not_leader += 1;
                        }
                        return Ok(());
                    }
                }
            }
            CountUnit::Fragment { tn5: false } => match r.fragment(self.spec.max_fragment) {
                Some(f) => Some(f),
                None if !r.is_paired() => None,
                None => {
                    if r.fragment(0).is_some() {
                        b.status.too_long += 1;
                    } else {
                        b.status.not_leader += 1;
                    }
                    return Ok(());
                }
            },
            _ => None,
        };
        let Some(row) = self.rows.of(r, &mut b.local) else {
            b.status.no_row += 1;
            return Ok(());
        };
        let sub = match self.spec.dedup {
            Dedup::None => 0,
            Dedup::Umi => match r.tag_bytes(self.spec.umi_tag) {
                Some(u) => umi_key(u),
                None => {
                    b.status.no_umi += 1;
                    return Ok(());
                }
            },
            Dedup::Position => match fragment {
                Some(f) => (f.start as u64) << 32 | f.end as u64,
                None => (r.pos() as u64) << 32 | r.end() as u64,
            },
        };
        let strand = self.strand(r);
        let conv = match &self.spec.labeled {
            Some(l) => {
                let c = convert(r, &l.conversion)?;
                if c.is_none() {
                    b.status.no_reference += 1;
                }
                c
            }
            None => None,
        };
        let rev = r.reverse_for(self.spec.strand);
        if let Some(v) = &self.velocity {
            // A molecule: the UMI; without one, the read, or the pair by
            // its name.
            let molecule = match (self.spec.dedup, self.spec.unit) {
                (Dedup::Umi, _) => sub,
                (_, CountUnit::Pair) => fnv(r.name()),
                _ => ((r.input as u64) << 48) ^ r.ordinal,
            };
            self.velocity(v, r, row, molecule, strand, &mut b.entries, &mut b.status);
        }
        match self.spec.unit {
            CountUnit::Read => {
                let hits = r.assign(self.space, Query::Blocks, strand);
                let st = &mut b.status;
                self.decide(&hits, row, sub, conv.as_ref(), rev, &mut b.entries, st);
            }
            CountUnit::Pair => {
                let hits = r.assign(self.space, Query::Blocks, strand);
                if expects_mate(r) {
                    let half = Half {
                        row,
                        sub,
                        hits,
                        conv,
                        rev,
                    };
                    if let Some((m, h)) = b.pairs.offer(r, half) {
                        self.pair(m, h, &mut b.entries, &mut b.status);
                    }
                } else {
                    let st = &mut b.status;
                    self.decide(&hits, row, sub, conv.as_ref(), rev, &mut b.entries, st);
                }
            }
            CountUnit::Fragment { tn5 } => {
                let q = match (fragment, tn5) {
                    (None, _) => Query::Span,
                    (Some(_), true) => Query::AtacFragment,
                    (Some(_), false) => Query::Fragment,
                };
                let hits = r.assign(self.space, q, strand);
                self.decide(&hits, row, sub, None, rev, &mut b.entries, &mut b.status);
            }
            CountUnit::CutSite => {
                let a = r.assign(self.space, Query::CutStart, None);
                let z = r.assign(self.space, Query::CutEnd, None);
                let assigned = [&a, &z]
                    .iter()
                    .any(|h| !h.is_empty() && (h.len() == 1 || self.spec.multi_overlap));
                if !assigned {
                    if a.is_empty() && z.is_empty() {
                        b.status.no_features += 1;
                    } else {
                        b.status.ambiguous += 1;
                    }
                    return Ok(());
                }
                b.status.assigned += 1;
                let eligible = |h: &Hits| h.len() == 1 || self.spec.multi_overlap;
                let a = if eligible(&a) { a.as_slice() } else { &[] };
                let z = if eligible(&z) { z.as_slice() } else { &[] };
                visit_hits(a, z, |col, n| {
                    b.entries.push(Entry {
                        row,
                        col,
                        layer: 0,
                        sub,
                        value: n,
                    });
                });
            }
        }
        Ok(())
    }

    fn seal(&self, b: &mut CountBatch) -> Result<()> {
        // Local row numbers map one to one onto global ones, so entries
        // collapsed here stay distinct once renumbered.
        collapse(&mut b.entries, self.fold);
        b.dangling = std::mem::take(&mut b.pairs).into_dangling();
        Ok(())
    }

    fn total<'p>(&'p self, rc: &RunCtx<'p>) -> Result<CountTotal<'p>> {
        Ok(CountTotal {
            sorter: Sorter::new(rc.pool, rc.buffer, &rc.tmp),
            rows: GlobalRows::default(),
            mates: Mates::new(std::mem::size_of::<Half>()),
            status: Status::default(),
            tally: Tally::default(),
            scratch: Vec::new(),
        })
    }

    fn merge<'p>(
        &'p self,
        t: &mut CountTotal<'p>,
        b: CountBatch,
        frontier: u64,
        rc: &RunCtx<'p>,
    ) -> Result<()> {
        t.status.add(&b.status);
        t.tally.add(&b.tally);
        let CountBatch {
            mut entries,
            local,
            dangling,
            ..
        } = b;
        let map = if self.rows.local() {
            let map = t.rows.remap(local);
            for e in &mut entries {
                e.row = map[e.row as usize];
            }
            Some(map)
        } else {
            None
        };
        let out = &mut t.scratch;
        out.clear();
        for (key, mut half) in dangling {
            if let Some(m) = &map {
                half.row = m[half.row as usize];
            }
            if let Some((a, z)) = t.mates.offer(key, half) {
                self.pair(a, z, out, &mut t.status);
            }
        }
        let st = &mut t.status;
        t.mates.expire(frontier, |h| {
            self.decide(&h.hits, h.row, h.sub, h.conv.as_ref(), h.rev, out, st);
        });
        if t.rows.bytes() + t.mates.bytes() > rc.buffer.max(64 << 20) {
            return Err(Error::invalid(
                "budget",
                format!(
                    "{}: {} row keys and {} mates waiting for their pair outgrow the memory \
                     budget; give a whitelist or raise the memory",
                    self.spec.name,
                    t.rows.names.len(),
                    t.mates.len()
                ),
            ));
        }
        entries.append(out);
        t.sorter.extend(entries)
    }

    fn end_input<'p>(
        &'p self,
        t: &mut CountTotal<'p>,
        _input: u32,
        _rc: &RunCtx<'p>,
    ) -> Result<()> {
        let mut out = Vec::new();
        let st = &mut t.status;
        t.mates
            .drain(|h| self.decide(&h.hits, h.row, h.sub, h.conv.as_ref(), h.rev, &mut out, st));
        t.sorter.extend(out)
    }

    fn finish<'p>(
        &'p self,
        t: CountTotal<'p>,
        _rc: &RunCtx<'p>,
    ) -> Result<(CountMatrix, Status, Tally, u64)> {
        let spilled = t.sorter.spilled();
        let rows = match self.rows.fixed_names() {
            Some(n) => n.to_vec(),
            None => t.rows.names,
        };
        let merge = t.sorter.finish()?;
        Ok((
            CountMatrix {
                rows: BarcodeTable { names: rows },
                features: Arc::clone(&self.features) as Arc<dyn FeatureList>,
                layers: self.layers.clone(),
                named_layers: false,
                cells: Cells {
                    merge: Some(merge),
                    fold: self.fold,
                    head: None,
                },
            },
            t.status,
            t.tally,
            spilled,
        ))
    }

    fn deliver(
        &self,
        (matrix, status, tally, spilled): (CountMatrix, Status, Tally, u64),
        d: &mut Delivery<'_>,
    ) -> Result<Json> {
        let n_rows = matrix.rows.len();
        let n_features = matrix.features.len();
        let written = crate::sink::matrix::write(matrix, &self.spec.name, self.spec.table, d)?;
        let unit = match self.spec.unit {
            CountUnit::Read if self.spec.dedup == Dedup::Umi => "molecule",
            CountUnit::Read => "read",
            CountUnit::Pair => "pair",
            CountUnit::Fragment { .. } => "fragment",
            CountUnit::CutSite => "cut-site",
        };
        let mut fields = vec![
            ("product", Json::from("counts")),
            ("unit", Json::from(unit)),
            ("rows", Json::from(n_rows as u64)),
            ("features", Json::from(n_features as u64)),
            ("nonzero", Json::from(written.nonzero)),
            ("assigned", Json::from(status.assigned)),
            ("unassigned_no_features", Json::from(status.no_features)),
            ("unassigned_ambiguity", Json::from(status.ambiguous)),
            ("unassigned_no_row_key", Json::from(status.no_row)),
            ("unassigned_no_umi", Json::from(status.no_umi)),
            (
                "mates_counted_with_their_pair",
                Json::from(status.not_leader),
            ),
            ("fragments_too_long", Json::from(status.too_long)),
            ("spilled_bytes", Json::from(spilled)),
        ];
        if self.layers.len() > 1 {
            let layers = self.layers.iter().map(|l| Json::from(l.clone()));
            fields.push(("layers", Json::Arr(layers.collect())));
        }
        if self.velocity.is_some() {
            fields.extend([
                (
                    "velocity_unassigned_no_gene",
                    Json::from(status.velocity_no_gene),
                ),
                (
                    "velocity_unassigned_multi_gene",
                    Json::from(status.velocity_multi_gene),
                ),
                (
                    "velocity_unassigned_no_model",
                    Json::from(status.velocity_no_model),
                ),
            ]);
        }
        if let Some(l) = &self.spec.labeled {
            fields.extend([
                ("conversion", Json::from(l.conversion.name())),
                ("min_conversions", Json::from(l.min_conversions as u64)),
                ("unlabeled_no_reference", Json::from(status.no_reference)),
            ]);
        }
        fields.extend(tally.json());
        fields.push((
            "files",
            Json::Arr(written.files.into_iter().map(Json::from).collect()),
        ));
        Ok(Json::obj(fields))
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn mate_votes_and_cut_counts_match_set_membership() {
        for a_mask in 0u32..128 {
            let a: Vec<_> = (0..7).filter(|i| a_mask & (1 << i) != 0).collect();
            for b_mask in 0u32..128 {
                let b: Vec<_> = (0..7).filter(|i| b_mask & (1 << i) != 0).collect();
                let all: Vec<_> = (0..7).filter(|c| a.contains(c) || b.contains(c)).collect();
                let common: Vec<_> = (0..7).filter(|c| a.contains(c) && b.contains(c)).collect();
                let best = if common.is_empty() { &all } else { &common };
                assert_eq!(pair_hits(&a, &b, false).as_slice(), best);
                assert_eq!(pair_hits(&a, &b, true).as_slice(), all);
                let mut counts = Vec::new();
                visit_hits(&a, &b, |c, n| counts.push((c, n)));
                let expected: Vec<_> = all
                    .iter()
                    .map(|&c| (c, a.contains(&c) as u32 + b.contains(&c) as u32))
                    .collect();
                assert_eq!(counts, expected);
            }
        }
    }

    #[test]
    fn collapse_folds_equal_keys() {
        let e = |row, col, sub, value| Entry {
            row,
            col,
            layer: 0,
            sub,
            value,
        };
        let mut v = vec![e(1, 2, 0, 1), e(0, 5, 0, 1), e(1, 2, 0, 3), e(1, 2, 9, 1)];
        collapse(&mut v, Fold::plain(Reduce::Sum));
        assert_eq!(v, vec![e(0, 5, 0, 1), e(1, 2, 0, 4), e(1, 2, 9, 1)]);
        let mut v = vec![e(1, 2, 7, 2), e(1, 2, 7, 1)];
        collapse(&mut v, Fold::plain(Reduce::MaxPerSub));
        assert_eq!(v, vec![e(1, 2, 7, 2)]);
        // Velocity evidence ORs per molecule and classifies when read.
        let fold = Fold {
            reduce: Reduce::Distinct,
            spliced: Some(1),
            labeled: Some(4),
            direct: false,
        };
        let vel = |sub, value| Entry {
            row: 0,
            col: 0,
            layer: VELOCITY_LAYER,
            sub,
            value,
        };
        let mut v = vec![vel(3, 1), vel(3, 2), vel(4, 1)];
        collapse(&mut v, fold);
        assert_eq!(v, vec![vel(3, 3), vel(4, 1)]);
        assert_eq!(fold.close(VELOCITY_LAYER, 3), Some((3, 1)));
        assert_eq!(fold.close(VELOCITY_LAYER, 1), Some((1, 1)));
        assert_eq!(fold.close(LABEL_LAYER, 7), Some((4, 1)));
    }
}
