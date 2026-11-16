//! Reducers: what a scan keeps from the windows it scores.
//!
//! A reducer declares what it needs ([`Need`]); the scan lowers the union
//! of all needs into the kernel's mode (every window, windows above a
//! cutoff, each sequence's best window) and hands each reducer batches of
//! scored windows. Reducers run thread-locally ([`Reducer::fork`]) and are
//! merged at the end; results must not depend on how the input was split.

use std::sync::Arc;
use std::sync::atomic::{AtomicU32, Ordering};

use omicsweave_core::{Result, Strand, Strings};

use crate::model::{StrandModel, Threshold, Variant};
use crate::program::Program;

/// What a reducer needs from the kernel.
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum Need {
    /// Every valid window.
    All,
    /// Windows passing a threshold.
    Above(Threshold),
    /// Each sequence's highest-scoring window per motif and strand.
    Best,
    /// Windows passing a threshold, handed to [`Reducer::sequence`] one
    /// sequence stretch at a time, in sequence and position order.
    Ordered(Threshold),
}

/// A [`Need`], optionally limited to one part of the program.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Demand {
    pub need: Need,
    /// Only motifs of this part; `None` for every motif.
    pub part: Option<usize>,
}

impl From<Need> for Demand {
    fn from(need: Need) -> Demand {
        Demand { need, part: None }
    }
}

/// A window passed to the ordered stage, with its reported statistics.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Match {
    pub score: f64,
    pub pvalue: f64,
    pub motif: u32,
    /// 0-based, half-open.
    pub start: u32,
    pub end: u32,
    pub strand: Strand,
    /// Background context symbols (see [`StrandModel::variant`]).
    pub context: u8,
    /// The integer score.
    pub int_score: i32,
}

/// The ordered stage's view of one stretch of a sequence: every window
/// starting in `windows` that passes an [`Need::Ordered`] threshold, sorted
/// by start, then motif, then strand (forward first). A long sequence comes
/// as several stretches, in order; `first` and `last` mark its ends.
pub struct SequenceChunk<'a> {
    pub sequence: u32,
    /// Window starts this stretch covers.
    pub windows: std::ops::Range<u32>,
    pub first: bool,
    pub last: bool,
    pub matches: &'a [Match],
}

/// One scored window. `context` is how many background context symbols it
/// had (see [`StrandModel::variant`]); `score` is its integer score.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct Window {
    pub sequence: u32,
    /// 0-based start in the sequence.
    pub start: u32,
    pub context: u8,
    pub score: i32,
}

/// A batch of windows of one motif on one strand, in no particular order;
/// a sequence's windows may be spread over several batches.
pub struct Windows<'a> {
    pub motif: u32,
    /// The program part of the motif.
    pub part: u32,
    pub strand: Strand,
    pub windows: &'a [Window],
    pub(crate) model: &'a StrandModel,
    pub(crate) width: u32,
}

impl Windows<'_> {
    /// Whether the batch's motif is in `part` (any part for `None`).
    pub fn in_part(&self, part: Option<usize>) -> bool {
        part.is_none_or(|p| p == self.part as usize)
    }

    pub fn width(&self) -> u32 {
        self.width
    }

    pub fn variant(&self, w: &Window) -> &Variant {
        self.model.variant(w.context as usize)
    }

    /// The real (unrounded) score.
    pub fn raw(&self, w: &Window) -> f64 {
        self.variant(w).int_model().raw(w.score)
    }

    /// The reported score.
    pub fn score(&self, w: &Window) -> f64 {
        self.variant(w).score(w.score)
    }

    /// The reported p-value.
    pub fn pvalue(&self, w: &Window) -> f64 {
        self.variant(w).pvalue(w.score)
    }

    pub fn meets(&self, w: &Window, threshold: Threshold) -> bool {
        self.variant(w).meets(w.score, threshold)
    }
}

/// The input of a scan.
pub struct ScanInfo<'a> {
    pub sequences: Strings<'a>,
}

/// A consumer of scored windows.
pub trait Reducer: Send + Sized {
    type Output;

    /// Appends what this reducer needs.
    fn needs(&self, out: &mut Vec<Demand>);

    /// Called once before the scan.
    fn prepare(&mut self, _program: &Program, _info: &ScanInfo<'_>) -> Result<()> {
        Ok(())
    }

    /// An empty reducer for another thread.
    fn fork(&self) -> Self;

    /// Unordered batches, on worker threads (thread-local reducers).
    fn windows(&mut self, _batch: &Windows<'_>) {}

    /// The ordered stage, on one thread with the reducer passed to the scan:
    /// every sequence stretch in input order, whether or not it has matches.
    fn sequence(&mut self, _program: &Program, _chunk: &SequenceChunk<'_>) {}

    fn merge(&mut self, other: Self);

    fn finish(self, program: &Program, info: &ScanInfo<'_>) -> Self::Output;
}

/// One reported window.
///
/// numpy: `np.dtype([("score", "<f8"), ("pvalue", "<f8"), ("sequence", "<u4"),
/// ("motif", "<u4"), ("start", "<u4"), ("end", "<u4"), ("strand", "u1")],
/// align=True)` (40 bytes).
#[repr(C)]
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct HitEntry {
    pub score: f64,
    pub pvalue: f64,
    pub sequence: u32,
    pub motif: u32,
    /// 0-based, half-open.
    pub start: u32,
    pub end: u32,
    pub strand: Strand,
}

const _: () = assert!(std::mem::size_of::<HitEntry>() == 40);

/// Every window passing a threshold, ordered by motif, sequence, start and
/// strand (forward first), as FIMO lists them.
#[derive(Clone, Debug)]
pub struct Hits {
    threshold: Threshold,
    part: Option<usize>,
    hits: Vec<HitEntry>,
}

impl Hits {
    pub fn new(threshold: Threshold) -> Hits {
        Hits {
            threshold,
            part: None,
            hits: Vec::new(),
        }
    }

    /// Only motifs of program part `part`.
    pub fn in_part(mut self, part: usize) -> Hits {
        self.part = Some(part);
        self
    }

    pub fn pvalue(p: f64) -> Hits {
        Hits::new(Threshold::PValue(p))
    }
}

impl Reducer for Hits {
    type Output = Vec<HitEntry>;

    fn needs(&self, out: &mut Vec<Demand>) {
        out.push(Demand {
            need: Need::Above(self.threshold),
            part: self.part,
        });
    }

    fn fork(&self) -> Hits {
        Hits {
            threshold: self.threshold,
            part: self.part,
            hits: Vec::new(),
        }
    }

    fn windows(&mut self, batch: &Windows<'_>) {
        if !batch.in_part(self.part) {
            return;
        }
        for w in batch.windows {
            if batch.meets(w, self.threshold) {
                self.hits.push(HitEntry {
                    score: batch.score(w),
                    pvalue: batch.pvalue(w),
                    sequence: w.sequence,
                    motif: batch.motif,
                    start: w.start,
                    end: w.start + batch.width,
                    strand: batch.strand,
                });
            }
        }
    }

    fn merge(&mut self, mut other: Hits) {
        self.hits.append(&mut other.hits);
    }

    fn finish(mut self, _: &Program, _: &ScanInfo<'_>) -> Vec<HitEntry> {
        self.hits
            .sort_by_key(|h| (h.motif, h.sequence, h.start, h.strand));
        self.hits
    }
}

/// A sequence's best window for one motif.
///
/// numpy: `np.dtype([("score", "<f8"), ("pvalue", "<f8"), ("sequence", "<u4"),
/// ("motif", "<u4"), ("start", "<u4"), ("strand", "u1")], align=True)`
/// (32 bytes).
#[repr(C)]
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct BestEntry {
    pub score: f64,
    pub pvalue: f64,
    pub sequence: u32,
    pub motif: u32,
    pub start: u32,
    pub strand: Strand,
}

const _: () = assert!(std::mem::size_of::<BestEntry>() == 32);

/// The highest-scoring window of every (sequence, motif) pair with at least
/// one valid window, over both strands; ties go to the earlier start, then
/// the forward strand. Ordered by sequence, then motif.
#[derive(Clone, Debug, Default)]
pub struct SequenceBest {
    part: Option<usize>,
    best: Vec<(f64, BestEntry)>,
}

impl SequenceBest {
    pub fn new() -> SequenceBest {
        SequenceBest::default()
    }

    /// Only motifs of program part `part`.
    pub fn in_part(mut self, part: usize) -> SequenceBest {
        self.part = Some(part);
        self
    }
}

fn better(a: &(f64, BestEntry), b: &(f64, BestEntry)) -> bool {
    a.0 > b.0 || (a.0 == b.0 && (a.1.start, a.1.strand) < (b.1.start, b.1.strand))
}

impl Reducer for SequenceBest {
    type Output = Vec<BestEntry>;

    fn needs(&self, out: &mut Vec<Demand>) {
        out.push(Demand {
            need: Need::Best,
            part: self.part,
        });
    }

    fn fork(&self) -> SequenceBest {
        SequenceBest {
            part: self.part,
            best: Vec::new(),
        }
    }

    fn windows(&mut self, batch: &Windows<'_>) {
        if !batch.in_part(self.part) {
            return;
        }
        // Batches are mostly sorted; runs of one sequence collapse here and
        // `finish` merges the rest.
        let mut current: Option<(f64, BestEntry)> = None;
        for w in batch.windows {
            let cand = (
                batch.raw(w),
                BestEntry {
                    score: batch.score(w),
                    pvalue: batch.pvalue(w),
                    sequence: w.sequence,
                    motif: batch.motif,
                    start: w.start,
                    strand: batch.strand,
                },
            );
            match &mut current {
                Some(c) if c.1.sequence == w.sequence => {
                    if better(&cand, c) {
                        *c = cand;
                    }
                }
                _ => {
                    if let Some(c) = current.take() {
                        self.best.push(c);
                    }
                    current = Some(cand);
                }
            }
        }
        if let Some(c) = current {
            self.best.push(c);
        }
    }

    fn merge(&mut self, mut other: SequenceBest) {
        self.best.append(&mut other.best);
    }

    fn finish(mut self, _: &Program, _: &ScanInfo<'_>) -> Vec<BestEntry> {
        self.best.sort_by_key(|(_, e)| (e.sequence, e.motif));
        let mut out: Vec<(f64, BestEntry)> = Vec::new();
        for c in self.best {
            match out.last_mut() {
                Some(last) if (last.1.sequence, last.1.motif) == (c.1.sequence, c.1.motif) => {
                    if better(&c, last) {
                        *last = c;
                    }
                }
                _ => out.push(c),
            }
        }
        out.into_iter().map(|(_, e)| e).collect()
    }
}

/// Windows passing a threshold, counted per (sequence, motif) over both
/// strands, as a sparse sequence × motif matrix.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct CountMatrix {
    /// `indptr[i] .. indptr[i + 1]` are sequence `i`'s entries.
    pub indptr: Vec<u64>,
    pub motifs: Vec<u32>,
    pub counts: Vec<u32>,
}

#[derive(Clone, Debug)]
pub struct SequenceCounts {
    threshold: Threshold,
    part: Option<usize>,
    counts: Vec<(u32, u32, u32)>,
}

impl SequenceCounts {
    pub fn new(threshold: Threshold) -> SequenceCounts {
        SequenceCounts {
            threshold,
            part: None,
            counts: Vec::new(),
        }
    }

    /// Only motifs of program part `part`.
    pub fn in_part(mut self, part: usize) -> SequenceCounts {
        self.part = Some(part);
        self
    }
}

impl Reducer for SequenceCounts {
    type Output = CountMatrix;

    fn needs(&self, out: &mut Vec<Demand>) {
        out.push(Demand {
            need: Need::Above(self.threshold),
            part: self.part,
        });
    }

    fn fork(&self) -> SequenceCounts {
        SequenceCounts {
            threshold: self.threshold,
            part: self.part,
            counts: Vec::new(),
        }
    }

    fn windows(&mut self, batch: &Windows<'_>) {
        if !batch.in_part(self.part) {
            return;
        }
        for w in batch.windows {
            if batch.meets(w, self.threshold) {
                match self.counts.last_mut() {
                    Some((s, m, c)) if *s == w.sequence && *m == batch.motif => *c += 1,
                    _ => self.counts.push((w.sequence, batch.motif, 1)),
                }
            }
        }
    }

    fn merge(&mut self, mut other: SequenceCounts) {
        self.counts.append(&mut other.counts);
    }

    fn finish(mut self, _: &Program, info: &ScanInfo<'_>) -> CountMatrix {
        self.counts.sort_unstable_by_key(|&(s, m, _)| (s, m));
        let mut indptr = vec![0u64; info.sequences.len() + 1];
        let mut motifs = Vec::new();
        let mut counts: Vec<u32> = Vec::new();
        let mut last: Option<(u32, u32)> = None;
        for (s, m, c) in self.counts {
            if last == Some((s, m)) {
                *counts.last_mut().unwrap() += c;
            } else {
                motifs.push(m);
                counts.push(c);
                indptr[s as usize + 1] += 1;
                last = Some((s, m));
            }
        }
        for i in 0..info.sequences.len() {
            indptr[i + 1] += indptr[i];
        }
        CountMatrix {
            indptr,
            motifs,
            counts,
        }
    }
}

/// Every window's reported score: sequence `i` occupies
/// `values[offsets[i] .. offsets[i + 1]]` as a `[motif][strand][position]`
/// array of `len_i` positions per (motif, strand), NaN where no window
/// scores. Strands are in the program's order (forward first).
#[derive(Clone, Debug, PartialEq)]
pub struct DenseScores {
    pub offsets: Vec<u64>,
    pub values: Vec<f32>,
    pub motifs: usize,
    pub strands: usize,
}

/// Collects [`DenseScores`].
#[derive(Clone, Debug, Default)]
pub struct Dense {
    cells: Option<Arc<Vec<AtomicU32>>>,
    offsets: Arc<Vec<u64>>,
    lengths: Arc<Vec<u32>>,
    motifs: usize,
    strands: usize,
}

impl Dense {
    pub fn new() -> Dense {
        Dense::default()
    }
}

impl Reducer for Dense {
    type Output = DenseScores;

    fn needs(&self, out: &mut Vec<Demand>) {
        out.push(Need::All.into());
    }

    fn prepare(&mut self, program: &Program, info: &ScanInfo<'_>) -> Result<()> {
        self.motifs = program.len();
        self.strands = if program.is_empty() {
            0
        } else {
            program.strand_list().len()
        };
        let per = (self.motifs * self.strands) as u64;
        let mut offsets = Vec::with_capacity(info.sequences.len() + 1);
        let mut lengths = Vec::with_capacity(info.sequences.len());
        let mut at = 0u64;
        offsets.push(0);
        for i in 0..info.sequences.len() {
            let len = info.sequences.get(i).len() as u64;
            lengths.push(len as u32);
            at += len * per;
            offsets.push(at);
        }
        let nan = f32::NAN.to_bits();
        self.cells = Some(Arc::new((0..at).map(|_| AtomicU32::new(nan)).collect()));
        self.offsets = Arc::new(offsets);
        self.lengths = Arc::new(lengths);
        Ok(())
    }

    fn fork(&self) -> Dense {
        self.clone()
    }

    fn windows(&mut self, batch: &Windows<'_>) {
        let Some(cells) = &self.cells else { return };
        let strand = match batch.strand {
            Strand::Reverse if self.strands == 2 => 1,
            _ => 0,
        };
        let row = batch.motif as u64 * self.strands as u64 + strand;
        for w in batch.windows {
            let len = self.lengths[w.sequence as usize] as u64;
            let at = self.offsets[w.sequence as usize] + row * len + w.start as u64;
            cells[at as usize].store((batch.score(w) as f32).to_bits(), Ordering::Relaxed);
        }
    }

    fn merge(&mut self, _other: Dense) {}

    fn finish(self, _: &Program, _: &ScanInfo<'_>) -> DenseScores {
        let values = match self.cells {
            Some(cells) => Arc::try_unwrap(cells)
                .map(|v| {
                    v.into_iter()
                        .map(|a| f32::from_bits(a.into_inner()))
                        .collect()
                })
                .unwrap_or_else(|shared| {
                    shared
                        .iter()
                        .map(|a| f32::from_bits(a.load(Ordering::Relaxed)))
                        .collect()
                }),
            None => Vec::new(),
        };
        DenseScores {
            offsets: Arc::unwrap_or_clone(self.offsets),
            values,
            motifs: self.motifs,
            strands: self.strands,
        }
    }
}

/// Reducers chosen at run time (all of one type, typically an enum): each
/// sees every batch and the output keeps their order.
impl<R: Reducer> Reducer for Vec<R> {
    type Output = Vec<R::Output>;

    fn needs(&self, out: &mut Vec<Demand>) {
        for r in self {
            r.needs(out);
        }
    }

    fn prepare(&mut self, program: &Program, info: &ScanInfo<'_>) -> Result<()> {
        for r in self.iter_mut() {
            r.prepare(program, info)?;
        }
        Ok(())
    }

    fn fork(&self) -> Self {
        self.iter().map(R::fork).collect()
    }

    fn windows(&mut self, batch: &Windows<'_>) {
        for r in self.iter_mut() {
            r.windows(batch);
        }
    }

    fn sequence(&mut self, program: &Program, chunk: &SequenceChunk<'_>) {
        for r in self.iter_mut() {
            r.sequence(program, chunk);
        }
    }

    fn merge(&mut self, other: Self) {
        for (r, o) in self.iter_mut().zip(other) {
            r.merge(o);
        }
    }

    fn finish(self, program: &Program, info: &ScanInfo<'_>) -> Self::Output {
        self.into_iter().map(|r| r.finish(program, info)).collect()
    }
}

macro_rules! tuple_reducer {
    ($($name:ident $idx:tt),+) => {
        impl<$($name: Reducer),+> Reducer for ($($name,)+) {
            type Output = ($($name::Output,)+);

            fn needs(&self, out: &mut Vec<Demand>) {
                $(self.$idx.needs(out);)+
            }

            fn prepare(&mut self, program: &Program, info: &ScanInfo<'_>) -> Result<()> {
                $(self.$idx.prepare(program, info)?;)+
                Ok(())
            }

            fn fork(&self) -> Self {
                ($(self.$idx.fork(),)+)
            }

            fn windows(&mut self, batch: &Windows<'_>) {
                $(self.$idx.windows(batch);)+
            }

            fn sequence(&mut self, program: &Program, chunk: &SequenceChunk<'_>) {
                $(self.$idx.sequence(program, chunk);)+
            }

            fn merge(&mut self, other: Self) {
                $(self.$idx.merge(other.$idx);)+
            }

            fn finish(self, program: &Program, info: &ScanInfo<'_>) -> Self::Output {
                ($(self.$idx.finish(program, info),)+)
            }
        }
    };
}

tuple_reducer!(A 0, B 1);
tuple_reducer!(A 0, B 1, C 2);
tuple_reducer!(A 0, B 1, C 2, D 3);
