//! MAST: sequences ranked by the combined significance of each motif's best
//! match, with a non-overlapping tiling of strong matches.
//!
//! The pieces are the framework's: [`scoring`] is a [`Scoring`] preset
//! (MAST's score units, scaling and p-value letters), and [`Mast`] is a
//! [`Reducer`] that takes each motif's best window per sequence
//! ([`Need::Best`]) and the ordered strong windows ([`Need::Ordered`]). How
//! window p-values become sequence p-values, how motifs combine and how
//! E-values are formed are the [`crate::stats`] interfaces in
//! [`MastOptions`]; MAST's own formulas are the defaults.

use std::collections::BTreeMap;
use std::io::{self, Write};
use std::sync::Arc;

use omicsweave_core::{Error, Result, Strand, Strings};
use omicsweave_ffm::{Alphabet, MotifView};

use crate::background::Background;
use crate::fimo::{format_e, format_g};
use crate::model::Threshold;
use crate::program::Program;
use crate::reduce::{Demand, Match, Need, Reducer, ScanInfo, SequenceChunk, Windows};
use crate::scoring::{
    Ambiguity, AmbiguityRule, BackgroundRule, LogScale, Quantizer, Rounding, Scoring, Strands,
    meme_round,
};
use crate::stats::{Adjust, BestOf, Combine, EValue, MemeBestOf, Product};

/// MAST's built-in nucleotide frequencies for p-values (`ACGT`).
const MAST_DNA: [f64; 4] = [0.281475655, 0.221785822, 0.228634607, 0.267048106];

/// MAST's built-in amino-acid frequencies for p-values, alphabetical.
const MAST_PROTEIN: [f64; 20] = [
    0.073091885,
    0.018145453,
    0.051687956,
    0.062278511,
    0.040243411,
    0.069259642,
    0.022405456,
    0.056227000,
    0.058435042,
    0.091621836,
    0.023044274,
    0.046032137,
    0.050623807,
    0.040715284,
    0.051846246,
    0.073729031,
    0.059352333,
    0.064298546,
    0.013328158,
    0.032649745,
];

const PROTEIN: &[u8] = b"ACDEFGHIKLMNPQRSTVWY";

/// How the two strands of a sequence are treated.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub enum StrandMode {
    /// Both strands together: a motif's best match is the better strand's
    /// and counts twice as many windows.
    Combine,
    /// Each strand scored and reported on its own.
    Separate,
    /// The forward strand only (or an alphabet without complement).
    Forward,
}

impl StrandMode {
    pub fn strands(self) -> Strands {
        match self {
            StrandMode::Combine | StrandMode::Separate => Strands::Both,
            StrandMode::Forward => Strands::Forward,
        }
    }
}

/// The letter probabilities MAST computes p-values with: `given` (a
/// background file's raw values) or MAST's built-in frequencies, each kept
/// to 8 significant digits and averaged with its complement when both
/// strands are scanned; not normalized.
pub fn null_letters(
    alphabet: &Alphabet,
    given: Option<&[f64]>,
    mode: StrandMode,
) -> Result<Vec<f64>> {
    let mut letters = match given {
        Some(g) if g.len() == alphabet.len() => g.to_vec(),
        Some(g) => {
            return Err(Error::invalid(
                "MAST letters",
                format!("{} letters for an alphabet of {}", g.len(), alphabet.len()),
            ));
        }
        None => match alphabet.symbols() {
            b"ACGT" | b"ACGU" => MAST_DNA.to_vec(),
            s if s == PROTEIN => MAST_PROTEIN.to_vec(),
            _ => {
                return Err(Error::invalid(
                    "MAST letters",
                    "built-in frequencies exist for nucleotides and amino acids; give letters",
                ));
            }
        },
    };
    for v in &mut letters {
        *v = meme_round(*v, 8);
    }
    if mode != StrandMode::Forward
        && let Some(comp) = alphabet.complement()
    {
        for i in 0..letters.len() {
            let c = comp[i] as usize;
            if comp[c] as usize == i && i < c {
                let avg = (letters[i] + letters[c]) / 2.0;
                letters[i] = avg;
                letters[c] = avg;
            }
        }
    }
    Ok(letters)
}

/// MAST's ambiguity codes: IUPAC nucleotides with `X` and `.` for `N`, or
/// `B Z J X` (with `*` and `.` for `X`) for amino acids.
pub fn ambiguity(alphabet: &Alphabet) -> Result<Ambiguity> {
    match alphabet.symbols() {
        b"ACGT" | b"ACGU" => Ambiguity::iupac(alphabet, AmbiguityRule::Mean)?
            .with_alias(b'X', b'N')?
            .with_alias(b'x', b'N')?
            .with_alias(b'.', b'N'),
        s if s == PROTEIN => {
            let sets: [(u8, Vec<u8>); 4] = [
                (b'B', vec![2, 11]),
                (b'Z', vec![3, 13]),
                (b'J', vec![7, 9]),
                (b'X', (0..20).collect()),
            ];
            let mut codes = Vec::new();
            for (b, set) in sets {
                for c in [b, b.to_ascii_lowercase()] {
                    if alphabet.lut()[c as usize] == omicsweave_ffm::INVALID {
                        codes.push((c, set.clone()));
                    }
                }
            }
            Ambiguity {
                rule: AmbiguityRule::Mean,
                codes,
                unknown: None,
            }
            .with_alias(b'*', b'X')?
            .with_alias(b'.', b'X')
        }
        _ => Ok(Ambiguity::skip()),
    }
}

/// MAST's scoring: pseudocount 0.01 over 100 sites by default (integer
/// part of `sites`), rows within 1e-5 of 1 kept, scores in hundredths of a
/// bit, `Quantizer::Mast { range: 100 }`, MEME's background preparation for
/// scores, [`null_letters`] for p-values with one table per motif, and
/// ambiguity codes scored by their weighted mean.
pub fn scoring(
    alphabet: &Alphabet,
    null_given: Option<&[f64]>,
    mode: StrandMode,
) -> Result<Scoring> {
    Ok(Scoring {
        pseudocount: 0.01,
        default_sites: 100.0,
        truncate_sites: true,
        normalize_tolerance: 0.00001,
        log: LogScale::Centibits,
        quantizer: Quantizer::Mast { range: 100 },
        rounding: Rounding::None,
        background: BackgroundRule::Meme,
        null_budget: 1 << 33,
        null_letters: Some(null_letters(alphabet, null_given, mode)?),
        shared_reverse_null: true,
        column_prior: None,
        ambiguity: ambiguity(alphabet)?,
    })
}

/// Compiles motifs with [`scoring`]. `background` is the score background
/// (MAST uses the motif file's unless given another); `null_given` the raw
/// letters of a background file for p-values.
pub fn compile(
    motifs: &[MotifView<'_>],
    alphabet: &Alphabet,
    background: &Background,
    null_given: Option<&[f64]>,
    mode: StrandMode,
    threads: Option<usize>,
) -> Result<Program> {
    let mode = if alphabet.complement().is_none() {
        StrandMode::Forward
    } else {
        mode
    };
    Program::compile(
        motifs,
        alphabet,
        background,
        &scoring(alphabet, null_given, mode)?,
        mode.strands(),
        threads,
    )
}

/// MAST's motif similarity: for each pair, the best mean Pearson
/// correlation of aligned columns of their integer scores (forward strand),
/// over every overlap of one motif's start along the other.
/// `out[i][j]` for the program's motifs `motifs[i]` and `motifs[j]`.
pub fn correlations(program: &Program, motifs: &[usize]) -> Vec<Vec<f64>> {
    let n = program.alphabet().len();
    let columns: Vec<Vec<&[i32]>> = motifs
        .iter()
        .map(|&m| {
            let int = program.motif(m).strands()[0].full().int_model();
            let mut cols: Vec<(i32, &[i32])> = int
                .terms
                .iter()
                .filter(|t| t.positions.len() == 1)
                .map(|t| (t.positions[0], &t.table[..n]))
                .collect();
            cols.sort_by_key(|c| c.0);
            cols.into_iter().map(|c| c.1).collect()
        })
        .collect();
    let means: Vec<Vec<f64>> = columns
        .iter()
        .map(|cols| {
            cols.iter()
                .map(|c| {
                    let mut sum = 0.0;
                    for &v in c.iter() {
                        sum += v as f64;
                    }
                    sum / n as f64
                })
                .collect()
        })
        .collect();
    let k = motifs.len();
    let mut out = vec![vec![0.0; k]; k];
    for i in 0..k {
        for j in 0..i {
            let mut best = f64::MIN;
            for (a, b) in [(i, j), (j, i)] {
                let (w1, w2) = (columns[a].len(), columns[b].len());
                for shift in 0..w2 {
                    let mut rsum = 0.0;
                    let mut l = 0;
                    let mut m = shift;
                    while l < w1 && m < w2 {
                        let (mut s1, mut s2, mut s3) = (0.0, 0.0, 0.0);
                        for (&cu, &cv) in columns[a][l].iter().zip(columns[b][m]) {
                            let u = cu as f64 - means[a][l];
                            let v = cv as f64 - means[b][m];
                            s1 += u * v;
                            s2 += u * u;
                            s3 += v * v;
                        }
                        let denom = (s2 * s3).sqrt();
                        rsum += if denom != 0.0 { s1 / denom } else { 1.0 };
                        l += 1;
                        m += 1;
                    }
                    best = best.max(rsum);
                }
            }
            let c = best / columns[i].len().min(columns[j].len()) as f64;
            out[i][j] = c;
            out[j][i] = c;
        }
    }
    out
}

/// MAST's redundancy rule: motif `i` is redundant when it correlates above
/// `max` with an earlier motif that is not itself redundant.
pub fn redundant(correlations: &[Vec<f64>], max: f64) -> Vec<bool> {
    let mut bad = vec![false; correlations.len()];
    for i in 1..correlations.len() {
        bad[i] = (0..i).any(|j| !bad[j] && correlations[i][j] > max);
    }
    bad
}

/// Statistics and thresholds of a MAST run.
#[derive(Clone)]
pub struct MastOptions {
    pub mode: StrandMode,
    /// Windows with position p-value below this enter the tiling (MAST's
    /// `-mt`; `-w` uses ten times it).
    pub hit_threshold: f64,
    /// Sequences (strands) with E-value at most this are ranked.
    pub evalue_threshold: f64,
    /// Sequence p-value of a motif's best window among its windows.
    pub best_of: Arc<dyn BestOf>,
    /// Combined p-value of a sequence over its motifs.
    pub combine: Arc<dyn Combine>,
    /// E-values from the combined p-values of every sequence (strand).
    pub adjust: Arc<dyn Adjust>,
}

impl MastOptions {
    /// MAST's defaults: `-mt 1e-4`, `-ev 10`, `p·n` / `1 − (1 − p)^n`,
    /// the product of p-values, E = p · number of sequences (strands).
    pub fn new(mode: StrandMode) -> MastOptions {
        MastOptions {
            mode,
            hit_threshold: 1e-4,
            evalue_threshold: 10.0,
            best_of: Arc::new(MemeBestOf),
            combine: Arc::new(Product),
            adjust: Arc::new(EValue::default()),
        }
    }
}

/// A match in a sequence's tiling.
#[derive(Clone, Debug, PartialEq)]
pub struct MastHit {
    pub motif: u32,
    /// 0-based, half-open.
    pub start: u32,
    pub end: u32,
    pub strand: Strand,
    /// The real score (hundredths of a bit under [`scoring`]).
    pub score: f64,
    pub pvalue: f64,
}

/// A motif's best window in a sequence.
#[derive(Clone, Debug, PartialEq)]
pub struct MotifBest {
    pub motif: u32,
    pub score: f64,
    /// Position p-value of the score.
    pub pvalue: f64,
    /// The p-value of the best of the sequence's windows.
    pub sequence_pvalue: f64,
}

/// One sequence (one strand of it with [`StrandMode::Separate`]).
#[derive(Clone, Debug, PartialEq)]
pub struct MastSequence {
    pub sequence: u32,
    /// `Some` with [`StrandMode::Separate`].
    pub strand: Option<Strand>,
    /// Combined p-value over the motifs that fit in the sequence.
    pub pvalue: f64,
    pub evalue: f64,
    pub best: Vec<MotifBest>,
    /// Non-overlapping strong matches, by start.
    pub hits: Vec<MastHit>,
}

/// Every sequence in input order, and the ranked ones.
#[derive(Clone, Debug, PartialEq)]
pub struct MastReport {
    pub sequences: Vec<MastSequence>,
    /// Indices into `sequences` within the E-value threshold, by p-value,
    /// then input order, then strand (forward first).
    pub ranked: Vec<usize>,
}

#[derive(Clone, Copy)]
struct BestScore {
    score: f64,
    pvalue: f64,
    start: u32,
}

impl BestScore {
    fn better_than(self, other: Self) -> bool {
        self.score > other.score || (self.score == other.score && self.start < other.start)
    }
}

/// The MAST reducer. Equal raw scores choose the earlier start, then the
/// forward strand, retaining that window's own p-value.
#[derive(Clone)]
pub struct Mast {
    options: MastOptions,
    part: Option<usize>,
    /// Best window statistics per (sequence, motif, strand) batch run.
    best: Vec<(u32, u32, Strand, BestScore)>,
    tilings: Vec<(u32, Option<Strand>, Vec<MastHit>)>,
    pending: Vec<Match>,
}

impl Mast {
    pub fn new(options: MastOptions) -> Mast {
        Mast {
            options,
            part: None,
            best: Vec::new(),
            tilings: Vec::new(),
            pending: Vec::new(),
        }
    }

    /// Only motifs of program part `part`.
    pub fn in_part(mut self, part: usize) -> Mast {
        self.part = Some(part);
        self
    }

    fn motifs(&self, program: &Program) -> std::ops::Range<usize> {
        match self.part {
            Some(p) => program.parts()[p].motifs(),
            None => 0..program.len(),
        }
    }
}

impl Reducer for Mast {
    type Output = MastReport;

    fn prepare(&mut self, program: &Program, _info: &ScanInfo<'_>) -> Result<()> {
        if self.part.is_some_and(|p| p >= program.parts().len()) {
            return Err(Error::invalid("MAST", "program part is out of bounds"));
        }
        Ok(())
    }

    fn needs(&self, out: &mut Vec<Demand>) {
        out.push(Demand {
            need: Need::Best,
            part: self.part,
        });
        out.push(Demand {
            need: Need::Ordered(Threshold::PValue(self.options.hit_threshold)),
            part: self.part,
        });
    }

    fn fork(&self) -> Mast {
        Mast {
            options: self.options.clone(),
            part: self.part,
            best: Vec::new(),
            tilings: Vec::new(),
            pending: Vec::new(),
        }
    }

    fn windows(&mut self, batch: &Windows<'_>) {
        if !batch.in_part(self.part) {
            return;
        }
        let mut run: Option<(u32, BestScore)> = None;
        for w in batch.windows {
            let candidate = BestScore {
                score: batch.raw(w),
                pvalue: batch.pvalue(w),
                start: w.start,
            };
            match &mut run {
                Some((s, best)) if *s == w.sequence => {
                    if candidate.better_than(*best) {
                        *best = candidate;
                    }
                }
                _ => {
                    if let Some((s, best)) = run.take() {
                        self.best.push((s, batch.motif, batch.strand, best));
                    }
                    run = Some((w.sequence, candidate));
                }
            }
        }
        if let Some((s, best)) = run {
            self.best.push((s, batch.motif, batch.strand, best));
        }
    }

    fn sequence(&mut self, program: &Program, chunk: &SequenceChunk<'_>) {
        let part = self.part;
        self.pending.extend(
            chunk
                .matches
                .iter()
                .filter(|m| part.is_none_or(|p| program.part_index(m.motif as usize) == p)),
        );
        if !chunk.last {
            return;
        }
        let len = chunk.windows.end;
        let motifs: Vec<usize> = self.motifs(program).collect();
        let threshold = self.options.hit_threshold;
        let strands: &[Option<Strand>] = match self.options.mode {
            StrandMode::Separate => &[Some(Strand::Forward), Some(Strand::Reverse)],
            _ => &[None],
        };
        for &strand in strands {
            let hits = tile(
                program,
                &self.pending,
                &motifs,
                len,
                threshold,
                self.options.mode,
                strand,
            );
            self.tilings.push((chunk.sequence, strand, hits));
        }
        self.pending.clear();
    }

    fn merge(&mut self, mut other: Mast) {
        self.best.append(&mut other.best);
        self.tilings.append(&mut other.tilings);
    }

    fn finish(mut self, program: &Program, info: &ScanInfo<'_>) -> MastReport {
        let mode = self.options.mode;
        let motifs: Vec<usize> = self.motifs(program).collect();
        self.best.sort_unstable_by_key(|&(s, m, st, _)| (s, m, st));
        let mut best: BTreeMap<(u32, u32, Strand), BestScore> = BTreeMap::new();
        for (s, m, st, v) in self.best {
            let e = best.entry((s, m, st)).or_insert(v);
            if v.better_than(*e) {
                *e = v;
            }
        }
        let mut tilings: BTreeMap<(u32, Option<Strand>), Vec<MastHit>> = BTreeMap::new();
        for (s, st, hits) in self.tilings {
            tilings.insert((s, st), hits);
        }
        let strands: &[Option<Strand>] = match mode {
            StrandMode::Separate => &[Some(Strand::Forward), Some(Strand::Reverse)],
            _ => &[None],
        };
        let mut sequences = Vec::new();
        for si in 0..info.sequences.len() {
            let len = info.sequences.get(si).len() as u64;
            for &strand in strands {
                let mut bests = Vec::new();
                let mut pvalues = Vec::new();
                for &mi in &motifs {
                    let model = program.motif(mi);
                    let w = model.width() as u64;
                    if w > len {
                        continue;
                    }
                    let get = |st: Strand| best.get(&(si as u32, mi as u32, st)).copied();
                    let score = match (mode, strand) {
                        (StrandMode::Separate, Some(st)) => get(st),
                        (StrandMode::Combine, _) => {
                            match (get(Strand::Forward), get(Strand::Reverse)) {
                                (Some(a), Some(b)) => Some(if b.better_than(a) { b } else { a }),
                                (a, b) => a.or(b),
                            }
                        }
                        _ => get(Strand::Forward),
                    };
                    let Some(score) = score else { continue };
                    let p = score.pvalue;
                    let mut n = (len - w + 1) as f64;
                    if mode == StrandMode::Combine {
                        n *= 2.0;
                    }
                    let sp = self.options.best_of.pvalue(p, n);
                    pvalues.push(sp);
                    bests.push(MotifBest {
                        motif: mi as u32,
                        score: score.score,
                        pvalue: p,
                        sequence_pvalue: sp,
                    });
                }
                sequences.push(MastSequence {
                    sequence: si as u32,
                    strand,
                    pvalue: self.options.combine.combine(&pvalues),
                    evalue: 0.0,
                    best: bests,
                    hits: tilings.remove(&(si as u32, strand)).unwrap_or_default(),
                });
            }
        }
        let pvalues: Vec<f64> = sequences.iter().map(|s| s.pvalue).collect();
        for (s, e) in sequences
            .iter_mut()
            .zip(self.options.adjust.adjust(&pvalues))
        {
            s.evalue = e;
        }
        // MAST first keeps a strand when p · (sequences read so far) is
        // below the threshold, then when its E-value is within it.
        let per = strands.len() as f64;
        let ev = self.options.evalue_threshold;
        let mut ranked: Vec<usize> = (0..sequences.len())
            .filter(|&i| {
                let s = &sequences[i];
                s.pvalue * (s.sequence as f64 + 1.0) * per < ev && s.evalue <= ev
            })
            .collect();
        let strand_rank = |s: Option<Strand>| match s {
            Some(Strand::Forward) => 0,
            None => 1,
            _ => 2,
        };
        ranked.sort_by(|&a, &b| {
            let (x, y) = (&sequences[a], &sequences[b]);
            x.pvalue
                .total_cmp(&y.pvalue)
                .then(x.sequence.cmp(&y.sequence))
                .then(strand_rank(x.strand).cmp(&strand_rank(y.strand)))
        });
        MastReport { sequences, ranked }
    }
}

/// MAST's tiling of one sequence (one strand with `Separate`): motifs in
/// order, windows by start; a window with p-value below `threshold` is kept
/// when its p-value is below the product of the p-values of the kept
/// matches it overlaps, which it then replaces.
fn tile(
    program: &Program,
    matches: &[Match],
    motifs: &[usize],
    len: u32,
    threshold: f64,
    mode: StrandMode,
    strand: Option<Strand>,
) -> Vec<MastHit> {
    let raw = |m: &Match| {
        program
            .motif(m.motif as usize)
            .strand(m.strand)
            .map_or(m.score, |s| {
                s.variant(m.context as usize).int_model().raw(m.int_score)
            })
    };
    // The best strand per (motif, start); matches arrive forward first.
    let mut best: BTreeMap<(u32, u32), &Match> = BTreeMap::new();
    for m in matches {
        let wanted = match (mode, strand) {
            (StrandMode::Separate, Some(st)) => m.strand == st,
            (StrandMode::Forward, _) => m.strand == Strand::Forward,
            _ => true,
        };
        if !wanted {
            continue;
        }
        best.entry((m.motif, m.start))
            .and_modify(|b| {
                if raw(m) > raw(b) {
                    *b = m;
                }
            })
            .or_insert(m);
    }
    struct Kept {
        width: u32,
        pvalue: f64,
        hit: MastHit,
    }
    let mut kept: BTreeMap<u32, Kept> = BTreeMap::new();
    let mut max_width = 0u32;
    for &mi in motifs {
        let width = program.motif(mi).width();
        if width > len {
            continue;
        }
        max_width = max_width.max(width);
        for (&(_, j), m) in best.range((mi as u32, 0)..(mi as u32 + 1, 0)) {
            let p = m.pvalue;
            if p >= threshold {
                continue;
            }
            let first = j.saturating_sub(max_width - 1);
            let last = (j + width).min(len);
            let mut prod = 1.0;
            let mut ok = true;
            for k in kept.range(j..last).map(|(_, k)| k) {
                prod *= k.pvalue;
                if p >= prod {
                    ok = false;
                    break;
                }
            }
            if ok {
                for (&k, kp) in kept.range(first..j) {
                    if kp.width > j - k {
                        prod *= kp.pvalue;
                        if p >= prod {
                            ok = false;
                            break;
                        }
                    }
                }
            }
            if !ok {
                continue;
            }
            let right: Vec<u32> = kept.range(j + 1..last).map(|(&k, _)| k).collect();
            let left: Vec<u32> = kept
                .range(first..j)
                .filter(|(k, kp)| kp.width > j - **k)
                .map(|(&k, _)| k)
                .collect();
            for k in right.into_iter().chain(left) {
                kept.remove(&k);
            }
            let raw = raw(m);
            kept.insert(
                j,
                Kept {
                    width,
                    pvalue: p,
                    hit: MastHit {
                        motif: mi as u32,
                        start: j,
                        end: j + width,
                        strand: m.strand,
                        score: raw,
                        pvalue: p,
                    },
                },
            );
        }
    }
    kept.into_values().map(|k| k.hit).collect()
}

/// Names and numbering the MAST text outputs print.
pub struct MastText<'a> {
    /// Per sequence.
    pub sequence_names: Strings<'a>,
    /// Per program motif: ID, alternate ID (empty for none) and the motif's
    /// number in its source file (1-based).
    pub motif_ids: &'a [&'a str],
    pub motif_alt_ids: &'a [&'a str],
    pub motif_numbers: &'a [u32],
    /// Whether strands are printed (`false` for an alphabet without
    /// complement).
    pub stranded: bool,
}

/// Writes every sequence's tiling as `mast -hit_list` does (after its two
/// `#` header lines): `name ±motif id alt start end score p-value`.
pub fn write_hit_list<W: Write>(
    mut out: W,
    report: &MastReport,
    text: &MastText<'_>,
) -> io::Result<()> {
    for s in &report.sequences {
        let name = text.sequence_names.get(s.sequence as usize);
        for h in &s.hits {
            let m = h.motif as usize;
            let strand = if !text.stranded {
                ""
            } else if h.strand == Strand::Reverse {
                "-"
            } else {
                "+"
            };
            let alt = match text.motif_alt_ids.get(m).copied().unwrap_or("") {
                "" => "-",
                a => a,
            };
            out.write_all(name)?;
            writeln!(
                out,
                " {strand}{} {} {alt} {} {} {:8.2} {:>8}",
                text.motif_numbers.get(m).copied().unwrap_or(m as u32 + 1),
                text.motif_ids[m],
                h.start + 1,
                h.end,
                h.score,
                format_e(h.pvalue, 2),
            )?;
        }
    }
    Ok(())
}

/// The combined p-value (`%.2e`) and E-value (`%.2g`) of a sequence as
/// `mast.xml` prints them.
pub fn format_scores(s: &MastSequence) -> (String, String) {
    (format_e(s.pvalue, 2), format_g(s.evalue, 2))
}
