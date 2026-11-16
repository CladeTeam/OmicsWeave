//! MCAST: clusters of motif matches found by dynamic programming over the
//! ordered hits of each sequence.
//!
//! # Model
//!
//! Every motif strand is a separate *state* (motif 1 forward, motif 1
//! reverse, motif 2 forward, …). A window of state `m` starting at `s` emits
//!
//! ```text
//! e_m(s) = log2(t) − log2(p)   when its p-value p ≤ t (a hit, ≥ 0)
//!        = log2(t)             otherwise (negative)
//! ```
//!
//! with `t` the hit threshold ([`ClusterOptions::pthresh`]) and MEME's
//! `log2(x) = ln(x) · 1.44269504`. A cluster is a chain of windows separated
//! by gaps; a gap of `L ≥ 1` letters costs `L · D / G` (`G` =
//! [`ClusterOptions::max_gap`]), and every cluster pays `D` once
//! ([`cluster_threshold`]), so a cluster's score is the sum of its windows'
//! emissions minus its gap costs and only clusters scoring above `D` are
//! reported. Two windows with no letter between them are joined only by
//! closing the first cluster, which costs `D`. Gaps longer than `G` cost more
//! than `D`, so they split clusters; a gap of exactly `G` ties and splits.
//!
//! The best set of non-overlapping clusters is found by a Viterbi pass over
//! the sequence with a background row (outside any cluster), a spacer row
//! (inside a gap) and one row per state. Ties prefer, in order: staying in
//! the background, extending the spacer, the earlier state. As in MCAST,
//! long sequences are cut into overlapping segments of
//! [`ClusterOptions::segment`] letters, each solved on its own; a cluster
//! that starts in the last [`ClusterOptions::overlap`] letters of a segment
//! is left to the next one.
//!
//! # Statistics
//!
//! Cluster p-values come from a null of cluster scores drawn by a
//! [`ClusterNull`] and fitted by a [`ClusterTail`]; E-values and q-values are
//! [`Adjust`] implementations over all clusters' p-values. MCAST's choices
//! are one implementation each: [`MemeSynthetic`] (random sequences of
//! random GC content), [`GcBinnedExponential`] (exponential tails in GC
//! bins), [`EValue`] and [`MemeQValue`] (Storey's π₀ with MEME's bootstrap).

use std::io::{self, Write};
use std::sync::Arc;

use omicsweave_core::{Error, Result, Strand, Strings};
use omicsweave_ffm::Alphabet;

use crate::background::Background;
use crate::fimo::format_g;
use crate::model::Threshold;
use crate::program::Program;
use crate::reduce::{Demand, Need, Reducer, ScanInfo, SequenceChunk};
use crate::scan::ScanOptions;
use crate::scoring::{Ambiguity, AmbiguityRule, LogScale, Rounding, Scoring, Strands};
use crate::stats::{Adjust, EValue, TailFit};
use crate::substitution::SubstitutionPrior;

/// MCAST's default hit threshold.
pub const DEFAULT_PTHRESH: f64 = 0.0005;
/// MCAST's default largest gap inside a cluster.
pub const DEFAULT_MAX_GAP: u32 = 50;
/// Letters the segments of one sequence share.
pub const DEFAULT_OVERLAP: usize = 1000;
/// MCAST's dynamic-programming budget: a segment holds `budget / rows`
/// letters, `rows` being the summed widths of all states plus 3.
pub const SEGMENT_BUDGET: usize = 10_000_000;
/// MEME's estimate of a hit's expected score, `1 / ln 2` rounded.
const EXPECTED_HIT_SCORE: f64 = 1.442;

/// MEME's truncated `log2(e)` and `ln(2)`.
#[allow(clippy::approx_constant)]
const MEME_LOG2_E: f64 = 1.44269504;
#[allow(clippy::approx_constant)]
const MEME_LN_2: f64 = 0.69314718;

/// MCAST's scoring policy under `background` (MCAST uses NRDB frequencies
/// unless given a background file; the choice is the caller's here).
///
/// Rows within 1e-5 of 1 are kept as given, no pseudocount; each column is
/// smoothed by [`SubstitutionPrior::mcast_dna`]; MEME bits; `Fimo { range:
/// 1000 }` scaling over the letters and the IUPAC codes, a code scoring the
/// sum of its members (`X`, `x` and `.` spell `N`); no rounding; MEME's
/// background preparation for both strands; p-value tables under the
/// prepared letters passed through MEME's `log2`/`exp2` round trip. Compile
/// with [`Strands::Both`].
pub fn scoring(background: &Background, alphabet: &Alphabet) -> Result<Scoring> {
    let mut s = Scoring::fimo();
    s.pseudocount = 0.0;
    s.rounding = Rounding::None;
    let prepared = s.prepare_background(background, alphabet, Strands::Both)?;
    let letters = prepared.letters();
    s.null_letters = Some(table_letters(letters));
    s.column_prior = Some(Arc::new(SubstitutionPrior::mcast_dna(letters)?));
    s.ambiguity = Ambiguity::iupac(alphabet, AmbiguityRule::Sum)?
        .with_alias(b'X', b'N')?
        .with_alias(b'x', b'N')?
        .with_alias(b'.', b'N')?
        .with_unknown(b'N')?;
    Ok(s)
}

/// Each letter through MEME's `exp2(log2(b))` round trip, the letters its
/// p-value tables are computed with.
pub fn table_letters(letters: &[f64]) -> Vec<f64> {
    letters
        .iter()
        .map(|&b| ((b.ln() * MEME_LOG2_E) * MEME_LN_2).exp())
        .collect()
}

/// MEME's probability of a best-of-`n`: `p·n` below 0.001 (strictly),
/// `1 − (1 − p)^n` otherwise.
pub fn meme_ev(p: f64, n: f64) -> f64 {
    if p * n < 0.001 {
        p * n
    } else {
        1.0 - (1.0 - p).powf(n)
    }
}

/// MCAST's cluster threshold `D` for `states` motif states:
///
/// ```text
/// q    = ev(t, G · states)            P(a hit within G letters)
/// h    = ev(t, states)                P(a hit at a position)
/// egap = Σ_{i=0..=G} i · h (1 − h)^i / q
/// D    = egcost · G · 1.442 / egap    (1e17 when egap = 0)
/// ```
///
/// so that a gap of `G` letters costs the expected score of `egcost` hits.
pub fn cluster_threshold(pthresh: f64, max_gap: u32, states: usize, egcost: f64) -> f64 {
    let n = states as f64;
    let q = meme_ev(pthresh, max_gap as f64 * n);
    let phit = meme_ev(pthresh, n);
    let mut egap = 0.0;
    let mut cum = phit;
    for i in 0..=max_gap {
        egap += i as f64 * cum;
        cum *= 1.0 - phit;
    }
    egap /= q;
    if egap == 0.0 {
        1e17
    } else {
        egcost * max_gap as f64 * EXPECTED_HIT_SCORE / egap
    }
}

/// How clusters are found.
#[derive(Clone, Debug, PartialEq)]
pub struct ClusterOptions {
    /// Windows with a p-value at most this are hits.
    pub pthresh: f64,
    /// The gap length whose cost equals the cluster threshold.
    pub max_gap: u32,
    /// The cost of a `max_gap` gap in expected hit scores.
    pub egcost: f64,
    /// The cluster threshold `D`; `None` uses [`cluster_threshold`].
    pub threshold: Option<f64>,
    /// Letters per dynamic-programming segment; `None` uses MCAST's
    /// [`SEGMENT_BUDGET`] rule.
    pub segment: Option<usize>,
    /// Letters consecutive segments share.
    pub overlap: usize,
    /// Only motifs of this program part.
    pub part: Option<usize>,
}

impl Default for ClusterOptions {
    fn default() -> ClusterOptions {
        ClusterOptions {
            pthresh: DEFAULT_PTHRESH,
            max_gap: DEFAULT_MAX_GAP,
            egcost: 1.0,
            threshold: None,
            segment: None,
            overlap: DEFAULT_OVERLAP,
            part: None,
        }
    }
}

impl ClusterOptions {
    fn check(&self) -> Result<()> {
        let bad = |d: String| Err(Error::invalid("MCAST options", d));
        if !(self.pthresh > 0.0 && self.pthresh <= 1.0) {
            return bad(format!("hit threshold {} is not in (0, 1]", self.pthresh));
        }
        if self.max_gap == 0 {
            return bad("max gap must be at least 1".into());
        }
        if !(self.egcost.is_finite() && self.egcost > 0.0) {
            return bad("egcost must be > 0".into());
        }
        if let Some(d) = self.threshold
            && !(d.is_finite() && d > 0.0)
        {
            return bad(format!("threshold {d} is not > 0"));
        }
        if let Some(s) = self.segment
            && s <= self.overlap
        {
            return bad(format!(
                "segment {s} is not longer than the overlap {}",
                self.overlap
            ));
        }
        Ok(())
    }

    /// The motif states of `program` this search uses: (motif, strand,
    /// width) in MCAST's order.
    pub fn states(&self, program: &Program) -> Vec<(u32, Strand, u32)> {
        let motifs = match self.part {
            Some(p) => program.parts().get(p).map_or(0..0, |p| p.motifs()),
            None => 0..program.len(),
        };
        let mut out = Vec::new();
        for m in motifs {
            let model = program.motif(m);
            for s in model.strands() {
                out.push((m as u32, s.strand(), model.width()));
            }
        }
        out
    }

    /// The cluster threshold `D` used with `program`.
    pub fn threshold_for(&self, program: &Program) -> f64 {
        self.threshold.unwrap_or_else(|| {
            cluster_threshold(
                self.pthresh,
                self.max_gap,
                self.states(program).len(),
                self.egcost,
            )
        })
    }

    /// Letters per segment with `program`.
    pub fn segment_for(&self, program: &Program) -> usize {
        self.segment
            .unwrap_or_else(|| {
                let rows: u64 = self.states(program).iter().map(|s| s.2 as u64).sum::<u64>() + 3;
                (SEGMENT_BUDGET as u64 / rows) as usize
            })
            .max(self.overlap + 1)
    }
}

/// One cluster.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Cluster {
    pub sequence: u32,
    /// 0-based, half-open.
    pub start: u32,
    pub end: u32,
    /// The cluster's score: hit scores minus gap costs.
    pub score: f64,
    /// `score − D` as the search computes it (the null's sample value).
    pub excess: f64,
    /// Windows in the cluster.
    pub hits: u32,
    /// GC fraction of the cluster and 500 letters on each side, within its
    /// segment (the segment's two padding columns count as non-GC).
    pub gc: f64,
    /// First letter and length of the segment the cluster was found in.
    pub segment_start: u32,
    pub segment_len: u32,
}

#[derive(Clone, Copy, Debug)]
struct Hit {
    start: u32,
    state: u32,
    score: f64,
}

#[derive(Clone, Copy, Debug)]
struct State {
    width: usize,
}

/// The cluster search as a [`Reducer`]: hits arrive in order through the
/// ordered stage; each sequence's segments are solved as soon as all their
/// hits are known. The output is in discovery order (sequence, then
/// position).
#[derive(Clone, Debug)]
pub struct Clusters {
    options: ClusterOptions,
    states: Vec<State>,
    /// State of (motif, strand index) as `base[motif] + strand index`.
    base: Vec<Option<u32>>,
    strand_order: Vec<Vec<Strand>>,
    lengths: Vec<u32>,
    d: f64,
    extend: f64,
    half_open: f64,
    miss: f64,
    log_t: f64,
    segment: usize,
    current: Option<Current>,
    out: Vec<Cluster>,
}

#[derive(Clone, Debug)]
struct Current {
    sequence: u32,
    len: u32,
    hits: Vec<Hit>,
    delivered: u32,
    complete: bool,
    /// First letter of the next segment.
    at: u32,
    /// Next padded column to look for a cluster in (MCAST's carried start).
    start_pos: i64,
    done: bool,
}

impl Clusters {
    pub fn new(options: ClusterOptions) -> Clusters {
        Clusters {
            options,
            states: Vec::new(),
            base: Vec::new(),
            strand_order: Vec::new(),
            lengths: Vec::new(),
            d: 0.0,
            extend: 0.0,
            half_open: 0.0,
            miss: 0.0,
            log_t: 0.0,
            segment: 0,
            current: None,
            out: Vec::new(),
        }
    }

    fn state_of(&self, motif: u32, strand: Strand) -> Option<u32> {
        let base = (*self.base.get(motif as usize)?)?;
        let k = self.strand_order[motif as usize]
            .iter()
            .position(|&s| s == strand)?;
        Some(base + k as u32)
    }

    fn advance(&mut self) {
        let Some(mut cur) = self.current.take() else {
            return;
        };
        while !cur.done {
            let n = (self.segment as u32).min(cur.len - cur.at);
            let complete = cur.at + n == cur.len;
            if !(cur.complete || cur.delivered >= cur.at + n) {
                break;
            }
            self.solve(&mut cur, n, complete);
            if complete {
                cur.done = true;
                break;
            }
            let remove = n - self.options.overlap as u32;
            cur.at += remove;
            cur.start_pos = (cur.start_pos - remove as i64).max(0);
            let at = cur.at;
            cur.hits.retain(|h| h.start >= at);
        }
        self.current = Some(cur);
    }

    /// Solves the segment of `n` letters at `cur.at` and records its
    /// clusters.
    fn solve(&mut self, cur: &mut Current, n: u32, complete: bool) {
        let a = cur.at;
        let n = n as usize;
        let lp = n + 2;
        let s_count = self.states.len();
        let neg = f64::NEG_INFINITY;
        // Hits per state as (padded start, score), ascending.
        let mut per: Vec<Vec<(usize, f64)>> = vec![Vec::new(); s_count];
        for h in &cur.hits {
            let w = self.states[h.state as usize].width;
            if h.start >= a && (h.start - a) as usize + w <= n {
                per[h.state as usize].push(((h.start - a) as usize + 1, h.score));
            }
        }
        let mut cursor = vec![0usize; s_count];
        let mut f0 = vec![0.0f64; lp + 1];
        let mut t0 = vec![0u32; lp + 1];
        let mut e0 = vec![0.0f64; lp + 1];
        let mut fsp = vec![neg; lp];
        let mut tsp = vec![0u32; lp];
        let mut prev = vec![neg; s_count];
        let mut next = vec![neg; s_count];
        let mut tm = vec![false; s_count * lp];
        for j in 1..=lp {
            let mut best = f0[j - 1];
            let mut arg = 0u32;
            let mut endv = 0.0;
            for (m, &v) in prev.iter().enumerate() {
                let c = v - self.d;
                if best < c {
                    best = c;
                    arg = m as u32 + 1;
                    endv = v;
                }
            }
            f0[j] = best;
            t0[j] = arg;
            e0[j] = endv;
            if j == lp {
                break;
            }
            let mut b = f0[j];
            let mut from = 0u32;
            let c = fsp[j - 1] - self.extend;
            if c > b {
                b = c;
                from = 1;
            }
            for (m, &v) in prev.iter().enumerate() {
                let c = v - self.half_open;
                if c > b {
                    b = c;
                    from = m as u32 + 2;
                }
            }
            fsp[j] = b;
            tsp[j] = from;
            for m in 0..s_count {
                let w = self.states[m].width;
                if j < w {
                    next[m] = neg;
                    continue;
                }
                let s = j + 1 - w;
                // A window reaching into the closing padding scores 0.
                let e = if s == lp - w {
                    0.0
                } else {
                    let list = &per[m];
                    let c = &mut cursor[m];
                    while *c < list.len() && list[*c].0 < s {
                        *c += 1;
                    }
                    if *c < list.len() && list[*c].0 == s {
                        list[*c].1
                    } else {
                        self.miss
                    }
                };
                let mut b = f0[s];
                let mut sp = false;
                let c = fsp[s - 1] - self.half_open;
                if c > b {
                    b = c;
                    sp = true;
                }
                next[m] = b + e;
                tm[m * lp + j] = sp;
            }
            std::mem::swap(&mut prev, &mut next);
        }

        // Trace back from the background row after the last column.
        #[derive(Clone, Copy)]
        enum Row {
            Background,
            Spacer,
            Motif(usize),
        }
        let entry = |c: usize| -> Row {
            match t0[c] {
                0 => Row::Background,
                k => Row::Motif(k as usize - 1),
            }
        };
        let mut trace = vec![0u32; lp + 1];
        let mut ends = vec![false; lp + 1];
        let mut i = lp as isize;
        let mut row = Row::Background;
        while i >= 0 {
            let iu = i as usize;
            match row {
                Row::Background => {
                    trace[iu] = 0;
                    let k = if iu >= 1 { t0[iu] } else { 0 };
                    i -= 1;
                    row = if k == 0 {
                        Row::Background
                    } else {
                        Row::Motif(k as usize - 1)
                    };
                }
                Row::Spacer => {
                    trace[iu] = 1;
                    let k = tsp[iu];
                    i -= 1;
                    row = match k {
                        0 => entry(iu),
                        1 => Row::Spacer,
                        k => Row::Motif(k as usize - 2),
                    };
                }
                Row::Motif(m) => {
                    let w = self.states[m].width;
                    for k in 0..w {
                        trace[iu - k] = m as u32 + 2;
                    }
                    ends[iu] = true;
                    let sp = tm[m * lp + iu];
                    i -= w as isize;
                    row = if sp {
                        Row::Spacer
                    } else {
                        entry((i + 1) as usize)
                    };
                }
            }
        }

        // Clusters are runs of non-background columns.
        let mut pos = cur.start_pos.max(0) as usize;
        loop {
            let mut found = None;
            let mut start = None;
            for (c, &t) in trace.iter().enumerate().take(lp).skip(pos) {
                match start {
                    None if t != 0 => start = Some(c),
                    Some(s) if t == 0 || (complete && c == lp - 1) => {
                        found = Some((s, c - 1));
                        break;
                    }
                    _ => {}
                }
            }
            let Some((st, en)) = found else { break };
            if !complete && st + self.options.overlap > lp {
                break;
            }
            if st as i64 == cur.start_pos {
                cur.start_pos = en as i64 + 1;
                pos = en + 1;
                continue;
            }
            cur.start_pos = en as i64 + 1;
            pos = en + 1;
            let raw = e0[en + 1] - f0[st];
            let excess = raw - self.d;
            self.out.push(Cluster {
                sequence: cur.sequence,
                start: a + st as u32 - 1,
                end: a + en as u32,
                score: excess + self.d,
                excess,
                hits: ends[st..=en].iter().filter(|&&e| e).count() as u32,
                gc: f64::NAN,
                segment_start: a,
                segment_len: n as u32,
            });
        }
    }

    fn close(&mut self) {
        if let Some(cur) = &mut self.current
            && !cur.done
        {
            cur.complete = true;
        }
        self.advance();
        self.current = None;
    }
}

/// GC fraction around `c` in `sequence`, as [`Cluster::gc`] defines it.
fn cluster_gc(c: &Cluster, sequence: &[u8]) -> f64 {
    let a = c.segment_start as usize;
    let n = c.segment_len as usize;
    let lp = n + 2;
    let st = c.start as usize - a + 1;
    let en = c.end as usize - a;
    let lo = st.saturating_sub(500);
    let hi = (en + 500).min(lp - 1);
    let mut gc = 0u64;
    for p in lo + 1..=hi {
        if (1..=n).contains(&p) {
            gc += matches!(sequence[a + p - 1], b'C' | b'G' | b'c' | b'g') as u64;
        }
    }
    gc as f64 / (hi - lo) as f64
}

impl Reducer for Clusters {
    type Output = Vec<Cluster>;

    fn needs(&self, out: &mut Vec<Demand>) {
        out.push(Demand {
            need: Need::Ordered(Threshold::PValue(self.options.pthresh)),
            part: self.options.part,
        });
    }

    fn prepare(&mut self, program: &Program, info: &ScanInfo<'_>) -> Result<()> {
        self.options.check()?;
        let states = self.options.states(program);
        if states.len() > u32::MAX as usize - 2 {
            return Err(Error::invalid("MCAST", "too many motif states"));
        }
        self.base = vec![None; program.len()];
        self.strand_order = vec![Vec::new(); program.len()];
        self.states.clear();
        for (k, &(m, strand, w)) in states.iter().enumerate() {
            if self.base[m as usize].is_none() {
                self.base[m as usize] = Some(k as u32);
            }
            self.strand_order[m as usize].push(strand);
            self.states.push(State { width: w as usize });
        }
        self.d = self.options.threshold_for(program);
        self.extend = self.d / self.options.max_gap as f64;
        self.half_open = self.extend / 2.0;
        self.log_t = LogScale::MemeBits.log(self.options.pthresh);
        self.miss = self.log_t;
        self.segment = self.options.segment_for(program);
        self.lengths = (0..info.sequences.len())
            .map(|i| info.sequences.get(i).len() as u32)
            .collect();
        Ok(())
    }

    fn fork(&self) -> Clusters {
        Clusters::new(self.options.clone())
    }

    fn sequence(&mut self, _program: &Program, chunk: &SequenceChunk<'_>) {
        if self.current.as_ref().map(|c| c.sequence) != Some(chunk.sequence) {
            self.close();
            self.current = Some(Current {
                sequence: chunk.sequence,
                len: self.lengths[chunk.sequence as usize],
                hits: Vec::new(),
                delivered: 0,
                complete: false,
                at: 0,
                start_pos: 0,
                done: false,
            });
        }
        let mut hits = Vec::with_capacity(chunk.matches.len());
        for m in chunk.matches {
            if let Some(state) = self.state_of(m.motif, m.strand) {
                hits.push(Hit {
                    start: m.start,
                    state,
                    score: self.log_t - LogScale::MemeBits.log(m.pvalue),
                });
            }
        }
        let cur = self.current.as_mut().expect("current sequence");
        cur.hits.extend(hits);
        cur.delivered = cur.delivered.max(chunk.windows.end);
        if chunk.last {
            cur.complete = true;
        }
        self.advance();
    }

    fn merge(&mut self, mut other: Clusters) {
        self.out.append(&mut other.out);
    }

    fn finish(mut self, _program: &Program, info: &ScanInfo<'_>) -> Vec<Cluster> {
        self.close();
        for c in &mut self.out {
            c.gc = cluster_gc(c, info.sequences.get(c.sequence as usize));
        }
        self.out
    }
}

/// A cluster search bound to a program, used to score data and null
/// sequences alike.
pub struct ClusterSearch<'a> {
    pub program: &'a Program,
    pub options: &'a ClusterOptions,
    pub scan: &'a ScanOptions,
}

impl ClusterSearch<'_> {
    /// The clusters of `sequences`, in discovery order.
    pub fn clusters(&self, sequences: Strings<'_>) -> Result<Vec<Cluster>> {
        self.program
            .scan(sequences, Clusters::new(self.options.clone()), self.scan)
    }

    /// The cluster threshold `D`.
    pub fn threshold(&self) -> f64 {
        self.options.threshold_for(self.program)
    }
}

/// One null cluster score.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct NullSample {
    /// Score above the cluster threshold.
    pub excess: f64,
    pub gc: f64,
    /// Letters the cluster spans.
    pub span: u32,
    /// Letters of the segment it was found in.
    pub segment: u32,
    /// Tie-break among equal GC contents (larger first).
    pub order: u64,
}

impl NullSample {
    fn of(c: &Cluster, order: u64) -> NullSample {
        NullSample {
            excess: c.excess,
            gc: c.gc,
            span: c.end - c.start,
            segment: c.segment_len,
            order,
        }
    }
}

/// A source of null cluster scores.
pub trait ClusterNull: Send + Sync {
    /// Null samples for `search`, given the clusters found in the data.
    fn sample(&self, search: &ClusterSearch<'_>, observed: &[Cluster]) -> Result<Vec<NullSample>>;
}

/// A fitted null: the p-value of a score excess at a GC content.
pub type ClusterPvalue = Box<dyn Fn(f64, f64) -> f64 + Send + Sync>;

/// A null cluster score distribution fitted to samples. `Ok(None)` when the
/// samples do not support a fit (no statistics are reported then).
pub trait ClusterTail: Send + Sync {
    fn fit(&self, samples: &[NullSample]) -> Result<Option<ClusterPvalue>>;
}

/// MT19937 with MEME's draws: `drand` from one output, `ldrand` from two.
struct Mt19937 {
    mt: [u32; 624],
    i: usize,
}

impl Mt19937 {
    /// The reference seeding (`init_genrand`).
    fn new(seed: u32) -> Mt19937 {
        let mut mt = [0u32; 624];
        mt[0] = seed;
        for i in 1..624 {
            mt[i] = 1812433253u32
                .wrapping_mul(mt[i - 1] ^ (mt[i - 1] >> 30))
                .wrapping_add(i as u32);
        }
        Mt19937 { mt, i: 624 }
    }

    /// The original 1998 seeding (`sgenrand`, multiplier 69069), which MEME
    /// uses for its unseeded global generator (seed 4357).
    fn knuth(seed: u32) -> Mt19937 {
        let mut mt = [0u32; 624];
        mt[0] = seed;
        for i in 1..624 {
            mt[i] = 69069u32.wrapping_mul(mt[i - 1]);
        }
        Mt19937 { mt, i: 624 }
    }

    fn next(&mut self) -> u32 {
        if self.i >= 624 {
            for k in 0..624 {
                let y = (self.mt[k] & 0x8000_0000) | (self.mt[(k + 1) % 624] & 0x7fff_ffff);
                let mag = if y & 1 == 1 { 0x9908_b0df } else { 0 };
                self.mt[k] = self.mt[(k + 397) % 624] ^ (y >> 1) ^ mag;
            }
            self.i = 0;
        }
        let mut y = self.mt[self.i];
        self.i += 1;
        y ^= y >> 11;
        y ^= (y << 7) & 0x9d2c_5680;
        y ^= (y << 15) & 0xefc6_0000;
        y ^ (y >> 18)
    }

    fn drand(&mut self) -> f64 {
        self.next() as f64 * (1.0 / 4294967296.0)
    }

    fn ldrand(&mut self) -> f64 {
        let a = self.next() as u64;
        let b = self.next() as u64;
        ((a << 32) | b) as f64 * (1.0 / 18446744073709551616.0)
    }
}

/// MCAST's null: rounds of one random sequence each, with A = T = (1 − g)/2
/// and C = G = g/2 for a GC content `g` drawn uniformly between the lowest
/// and highest GC content of the data's clusters, until at least
/// `min_rounds` rounds have produced `want` clusters. It gives up after
/// `give_up_after` rounds once fewer than `want · round / max_rounds`
/// clusters were found, and stops at `max_rounds`. At most `max_samples`
/// clusters are kept by reservoir sampling. One MT19937 stream seeded with
/// `seed` draws, in MCAST's order: the data's reservoir, then per round the
/// GC content, the sequence and the round's reservoir. Needs an `ACGT`
/// alphabet.
#[derive(Clone, Debug)]
pub struct MemeSynthetic {
    pub seed: u32,
    pub round_length: usize,
    pub min_rounds: usize,
    pub max_rounds: usize,
    pub give_up_after: usize,
    pub want: f64,
    pub max_samples: usize,
}

impl MemeSynthetic {
    /// MCAST's settings: 1e6 letters per round, 100 to 1000 rounds, 1e5
    /// clusters wanted and kept.
    pub fn mcast(seed: u32) -> MemeSynthetic {
        MemeSynthetic {
            seed,
            round_length: 1_000_000,
            min_rounds: 100,
            max_rounds: 1000,
            give_up_after: 10,
            want: 1e5,
            max_samples: 100_000,
        }
    }
}

/// Reservoir sampling as MCAST does it: the `seen`-th item (1-based)
/// replaces slot `⌊seen · u⌋` when that is below the capacity.
fn reservoir<T>(slots: &mut Vec<T>, cap: usize, seen: u64, item: T, rng: &mut Mt19937) -> bool {
    if slots.len() < cap {
        slots.push(item);
        return true;
    }
    let idx = (seen as f64 * rng.drand()) as i64;
    if idx >= 0 && (idx as usize) < cap {
        slots[idx as usize] = item;
        return true;
    }
    false
}

impl ClusterNull for MemeSynthetic {
    fn sample(&self, search: &ClusterSearch<'_>, observed: &[Cluster]) -> Result<Vec<NullSample>> {
        if search.program.alphabet().symbols() != b"ACGT" {
            return Err(Error::invalid("MCAST null", "needs the ACGT alphabet"));
        }
        let mut rng = Mt19937::new(self.seed);
        if observed.is_empty() {
            return Ok(Vec::new());
        }
        let mut gcs: Vec<f64> = Vec::new();
        for (k, c) in observed.iter().enumerate() {
            reservoir(&mut gcs, self.max_samples, k as u64 + 1, c.gc, &mut rng);
        }
        let (mut min_gc, mut max_gc) = (1.0f64, 0.0f64);
        for &g in &gcs {
            min_gc = min_gc.min(g);
            max_gc = max_gc.max(g);
        }
        let mut samples: Vec<NullSample> = Vec::new();
        let mut seen = 0u64;
        let mut serial = 0u64;
        let mut seq = vec![0u8; self.round_length];
        for round in 1..=self.max_rounds {
            let gc = min_gc + rng.drand() * (max_gc - min_gc);
            let mut p = [(1.0 - gc) / 2.0, gc / 2.0, gc / 2.0, (1.0 - gc) / 2.0];
            let total = p[0] + p[1] + p[2] + p[3];
            if total != 1.0 {
                p.iter_mut().for_each(|v| *v /= total);
            }
            let cum = [p[0], p[0] + p[1], p[0] + p[1] + p[2], 1.0];
            let _length = rng.ldrand();
            for b in &mut seq {
                let x = rng.ldrand();
                let (mut lo, mut hi) = (0usize, 3usize);
                while lo < hi {
                    let mid = lo + (hi - lo) / 2;
                    if x < cum[mid] {
                        hi = mid;
                    } else {
                        lo = mid + 1;
                    }
                }
                *b = b"ACGT"[lo];
            }
            let offsets = [0i64, seq.len() as i64];
            let found = search.clusters(Strings {
                offsets: &offsets,
                data: &seq,
            })?;
            for c in &found {
                seen += 1;
                let s = NullSample::of(c, serial);
                if reservoir(&mut samples, self.max_samples, seen, s, &mut rng) {
                    serial += 1;
                }
            }
            let got = seen as f64;
            if round >= self.min_rounds && got >= self.want {
                break;
            }
            if round >= self.give_up_after
                && got / self.want < round as f64 / self.max_rounds as f64
            {
                break;
            }
        }
        Ok(samples)
    }
}

/// A null from sequences drawn from a background model
/// ([`Background::sample`]): `sequences` sequences of `length` letters.
#[derive(Clone, Debug)]
pub struct BackgroundNull {
    pub background: Background,
    pub length: usize,
    pub sequences: usize,
    pub seed: u64,
}

impl ClusterNull for BackgroundNull {
    fn sample(&self, search: &ClusterSearch<'_>, _observed: &[Cluster]) -> Result<Vec<NullSample>> {
        let lengths = vec![self.length; self.sequences];
        let (offsets, data) =
            self.background
                .sample(search.program.alphabet(), &lengths, self.seed);
        let found = search.clusters(Strings {
            offsets: &offsets,
            data: &data,
        })?;
        Ok(found
            .iter()
            .enumerate()
            .map(|(i, c)| NullSample::of(c, i as u64))
            .collect())
    }
}

/// MCAST's tail: exponential tails `p = exp(−x / μ)` fitted in GC-content
/// bins, μ interpolated linearly between the bins' mean GC contents.
///
/// Samples whose segment leaves fewer than `2 · max_gap` letters around
/// them (counting the two padding columns) are dropped. The rest, sorted by
/// GC content, are cut into bins of at least `min_scores` samples, a bin
/// ending once its running mean GC moves by the ratio `max_ratio` (on the
/// AT side above GC 0.5) and enough samples remain; μ is the bin's mean
/// excess. No fit when fewer than `min_scores` samples remain, or fewer
/// than `min_scores` samples have p-values above `1 / (non-outliers)`.
#[derive(Clone, Debug)]
pub struct GcBinnedExponential {
    pub min_scores: usize,
    pub max_ratio: f64,
    pub max_gap: u32,
}

impl GcBinnedExponential {
    pub fn mcast(max_gap: u32) -> GcBinnedExponential {
        GcBinnedExponential {
            min_scores: 1000,
            max_ratio: 1.01,
            max_gap,
        }
    }
}

#[derive(Clone, Copy, Debug)]
struct Bin {
    min: f64,
    max: f64,
    mid: f64,
    mu: f64,
}

fn binned_pvalue(bins: &[Bin], s: f64, t: f64) -> f64 {
    let n = bins.len();
    let mut g = bins
        .iter()
        .position(|b| t >= b.min && t <= b.max)
        .unwrap_or(n - 1);
    let g1;
    if t > bins[g].mid {
        g1 = if g < n - 1 { g + 1 } else { g };
    } else if t < bins[g].mid {
        g1 = g;
        g = g.saturating_sub(1);
    } else {
        g1 = g;
    }
    let mu = if g == g1 {
        bins[g].mu
    } else {
        bins[g].mu + (t - bins[g].mid) * ((bins[g1].mu - bins[g].mu) / (bins[g1].mid - bins[g].mid))
    };
    (-s / mu).exp()
}

impl ClusterTail for GcBinnedExponential {
    fn fit(&self, samples: &[NullSample]) -> Result<Option<ClusterPvalue>> {
        let min = self.min_scores.max(1);
        let ranges = samples.len() / min + 1;
        let room = 2 * self.max_gap as i64;
        let mut valid: Vec<NullSample> = samples
            .iter()
            .filter(|s| s.segment as i64 + 2 - s.span as i64 >= room)
            .copied()
            .collect();
        if valid.is_empty() || valid.len() < min {
            return Ok(None);
        }
        valid.sort_by(|a, b| a.gc.total_cmp(&b.gc).then(b.order.cmp(&a.order)));
        let n = valid.len();
        let per = if ranges == 1 { n } else { min };
        let mut bins: Vec<Bin> = Vec::new();
        let mut low = 0usize;
        let mut gc_mid = valid[0].gc / 2.0;
        let mut gc_hi = 0.0f64;
        for _ in 0..ranges {
            let mut gc_sum = 0.0;
            let mut i = low;
            let mut count = 1usize;
            while i < n {
                let mid_new = gc_sum / count as f64;
                let ratio = if gc_mid <= 0.5 {
                    mid_new / gc_mid
                } else {
                    (1.0 - gc_mid) / (1.0 - mid_new)
                };
                gc_sum += valid[i].gc;
                if i == n - 1 || (count >= per && ratio >= self.max_ratio && i + per < n) {
                    let mut x_sum = 0.0;
                    for v in &valid[low..=i] {
                        x_sum += v.excess;
                    }
                    let mean = x_sum / count as f64;
                    let expected_max = -mean * (1.0 / count as f64).ln();
                    if mean > expected_max {
                        return Ok(None);
                    }
                    let mid = gc_sum / count as f64;
                    bins.push(Bin {
                        min: gc_hi,
                        max: valid[i].gc,
                        mid,
                        mu: mean,
                    });
                    gc_mid = mid;
                    gc_hi = valid[i].gc;
                    low = i + 1;
                    break;
                }
                i += 1;
                count += 1;
            }
            if i >= n - 1 {
                break;
            }
        }
        if let Some(last) = bins.last_mut() {
            last.max = 1.0;
        }
        // Enough samples must remain after removing the outliers.
        let mut pv: Vec<f64> = samples
            .iter()
            .map(|s| binned_pvalue(&bins, s.excess, s.gc))
            .collect();
        pv.sort_by(f64::total_cmp);
        let total = pv.len();
        let mut out = 0usize;
        while out < total {
            let thresh = 1.0 / (total - out) as f64;
            let i = pv.iter().take_while(|&&p| p < thresh).count();
            if i > out {
                out = i;
            } else {
                break;
            }
        }
        if total - out < min {
            return Ok(None);
        }
        Ok(Some(Box::new(move |s, gc| binned_pvalue(&bins, s, gc))))
    }
}

/// Any [`TailFit`] over the excesses, ignoring GC content.
pub struct Unconditional<T: TailFit>(pub T);

impl<T: TailFit> ClusterTail for Unconditional<T> {
    fn fit(&self, samples: &[NullSample]) -> Result<Option<ClusterPvalue>> {
        if samples.is_empty() {
            return Ok(None);
        }
        let x: Vec<f64> = samples.iter().map(|s| s.excess).collect();
        let f = self.0.fit(&x)?;
        Ok(Some(Box::new(move |s, _| f(s))))
    }
}

/// MEME's q-values: Benjamini–Hochberg with Storey's π₀.
///
/// With more than `min_values` p-values, π₀(λ) = #{p ≥ λ} / ((1 − λ) m) is
/// computed for λ = 0.5·k/`lambdas` (k = 1 ..= `lambdas`); the λ whose
/// bootstrap estimates (`bootstraps` resamples of `samples` p-values) have
/// the least squared error to the smallest π₀(λ) gives π₀, clamped to
/// [0, 1]. The bootstrap draws from MT19937 seeded as MEME's unseeded
/// global generator (`seed`, 4357). Otherwise π₀ = 1. Then
/// `q_(i) = min_{j ≥ i} min(1, π₀ p_(j) m / j)`.
#[derive(Clone, Copy, Debug)]
pub struct MemeQValue {
    pub bootstraps: usize,
    pub samples: usize,
    pub lambdas: usize,
    pub max_lambda: f64,
    pub min_values: usize,
    pub seed: u32,
}

impl Default for MemeQValue {
    fn default() -> MemeQValue {
        MemeQValue {
            bootstraps: 100,
            samples: 1000,
            lambdas: 100,
            max_lambda: 0.5,
            min_values: 100,
            seed: 4357,
        }
    }
}

fn pi_zero_at(lambda: f64, sorted: &[f64]) -> f64 {
    let start = if sorted.first().is_some_and(|&p| p >= lambda) {
        0
    } else {
        sorted.partition_point(|&p| p < lambda)
    };
    let pi = (sorted.len() - start) as f64 / ((1.0 - lambda) * sorted.len() as f64);
    if pi == 0.0 { f64::INFINITY } else { pi }
}

impl MemeQValue {
    /// Storey's π₀ of ascending p-values.
    pub fn pi_zero(&self, sorted: &[f64]) -> f64 {
        if sorted.is_empty() || self.lambdas == 0 {
            return 1.0;
        }
        let lambda = |k: usize| ((k + 1) as f64 / self.lambdas as f64) * self.max_lambda;
        let pis: Vec<f64> = (0..self.lambdas)
            .map(|k| pi_zero_at(lambda(k), sorted))
            .collect();
        let min = pis.iter().copied().fold(f64::INFINITY, f64::min);
        let mut mse = vec![0.0f64; self.lambdas];
        let mut rng = Mt19937::knuth(self.seed);
        let mut boot = vec![0.0f64; self.samples];
        for _ in 0..self.bootstraps {
            for b in &mut boot {
                let idx = (rng.ldrand() * sorted.len() as f64) as usize;
                *b = sorted[idx.min(sorted.len() - 1)];
            }
            boot.sort_by(f64::total_cmp);
            for (k, e) in mse.iter_mut().enumerate() {
                let err = pi_zero_at(lambda(k), &boot) - min;
                *e += err * err;
            }
        }
        let mut best = 0;
        for k in 1..self.lambdas {
            if mse[k] < mse[best] {
                best = k;
            }
        }
        pis.get(best).copied().unwrap_or(1.0).clamp(0.0, 1.0)
    }
}

impl Adjust for MemeQValue {
    fn adjust(&self, pvalues: &[f64]) -> Vec<f64> {
        let n = pvalues.len();
        let mut order: Vec<usize> = (0..n).collect();
        order.sort_by(|&a, &b| pvalues[a].total_cmp(&pvalues[b]));
        let sorted: Vec<f64> = order.iter().map(|&i| pvalues[i]).collect();
        let pi0 = if n > self.min_values && self.lambdas > 0 {
            self.pi_zero(&sorted)
        } else {
            1.0
        };
        let m = n as f32 as f64;
        let mut q: Vec<f64> = sorted
            .iter()
            .enumerate()
            .map(|(k, &p)| ((pi0 * p * m) / ((k + 1) as f32 as f64)).min(1.0))
            .collect();
        if n > 1 {
            let mut prev = q[n - 1];
            for k in (0..n - 1).rev() {
                if prev < q[k] {
                    q[k] = prev;
                }
                prev = q[k];
            }
        }
        let mut out = vec![0.0; n];
        for (k, &i) in order.iter().enumerate() {
            out[i] = q[k];
        }
        out
    }
}

/// How cluster statistics are computed.
pub struct Statistics {
    pub null: Box<dyn ClusterNull>,
    pub tail: Box<dyn ClusterTail>,
    pub evalue: Box<dyn Adjust>,
    pub qvalue: Box<dyn Adjust>,
}

impl Statistics {
    /// MCAST's statistics with its `--seed`.
    pub fn mcast(seed: u32, max_gap: u32) -> Statistics {
        Statistics {
            null: Box::new(MemeSynthetic::mcast(seed)),
            tail: Box::new(GcBinnedExponential::mcast(max_gap)),
            evalue: Box::new(EValue { tests: None }),
            qvalue: Box::new(MemeQValue::default()),
        }
    }
}

/// A cluster with its statistics (NaN when none were computed).
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ScoredCluster {
    pub cluster: Cluster,
    /// 1-based rank in discovery order (MCAST's `cluster-N`).
    pub id: u32,
    pub pvalue: f64,
    pub evalue: f64,
    pub qvalue: f64,
}

/// The clusters of a search, by increasing p-value (then decreasing score,
/// decreasing GC, discovery order); in discovery order without statistics.
#[derive(Clone, Debug)]
pub struct Report {
    pub clusters: Vec<ScoredCluster>,
    /// Whether p-, E- and q-values were computed.
    pub statistics: bool,
}

/// Finds the clusters of `sequences` and, with `statistics`, their p-, E-
/// and q-values. `program` should be compiled with [`scoring`] and both
/// strands to reproduce MCAST.
pub fn search(
    program: &Program,
    sequences: Strings<'_>,
    options: &ClusterOptions,
    statistics: Option<&Statistics>,
    scan: &ScanOptions,
) -> Result<Report> {
    let s = ClusterSearch {
        program,
        options,
        scan,
    };
    let found = s.clusters(sequences)?;
    s.report(&found, statistics)
}

impl ClusterSearch<'_> {
    /// `found` (clusters of this search, in discovery order) with their
    /// statistics, when `statistics` is given and its null supports a fit.
    pub fn report(&self, found: &[Cluster], statistics: Option<&Statistics>) -> Result<Report> {
        let mut clusters: Vec<ScoredCluster> = found
            .iter()
            .enumerate()
            .map(|(i, &c)| ScoredCluster {
                cluster: c,
                id: i as u32 + 1,
                pvalue: f64::NAN,
                evalue: f64::NAN,
                qvalue: f64::NAN,
            })
            .collect();
        let unscored = |clusters| {
            Ok(Report {
                clusters,
                statistics: false,
            })
        };
        let Some(stats) = statistics else {
            return unscored(clusters);
        };
        if clusters.is_empty() {
            return unscored(clusters);
        }
        let samples = stats.null.sample(self, found)?;
        let Some(tail) = stats.tail.fit(&samples)? else {
            return unscored(clusters);
        };
        let d = self.threshold();
        let p: Vec<f64> = clusters
            .iter()
            .map(|c| tail(c.cluster.score - d, c.cluster.gc))
            .collect();
        let e = stats.evalue.adjust(&p);
        let q = stats.qvalue.adjust(&p);
        for (i, c) in clusters.iter_mut().enumerate() {
            c.pvalue = p[i];
            c.evalue = e[i];
            c.qvalue = q[i];
        }
        clusters.sort_by(|a, b| {
            a.pvalue
                .total_cmp(&b.pvalue)
                .then(b.cluster.score.total_cmp(&a.cluster.score))
                .then(b.cluster.gc.total_cmp(&a.cluster.gc))
                .then(a.id.cmp(&b.id))
        });
        Ok(Report {
            clusters,
            statistics: true,
        })
    }
}

/// Which clusters [`write_text`] prints when statistics exist (every
/// cluster otherwise). MCAST's defaults: `E ≤ 10`; `--output-pthresh x`
/// means `p ≤ x` with no E limit; `--output-qthresh x` means `q ≤ x` and
/// `E ≤ 1`.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct OutputFilter {
    pub evalue: f64,
    pub pvalue: f64,
    pub qvalue: f64,
}

impl Default for OutputFilter {
    fn default() -> OutputFilter {
        OutputFilter {
            evalue: 10.0,
            pvalue: 1.0,
            qvalue: 1.0,
        }
    }
}

/// What [`write_text`] needs besides the report.
pub struct McastText<'a> {
    /// Sequence names, one row per scanned sequence.
    pub sequence_names: Strings<'a>,
    /// The scanned sequences.
    pub sequences: Strings<'a>,
    /// Added to every reported position; `None` for 0.
    pub origins: Option<&'a [u64]>,
    pub filter: OutputFilter,
}

fn printable(b: u8) -> u8 {
    match b.to_ascii_uppercase() {
        b'U' => b'T',
        c @ (b'A' | b'C' | b'G' | b'T' | b'R' | b'Y' | b'S' | b'W' | b'K' | b'M' | b'B' | b'D'
        | b'H' | b'V' | b'N') => c,
        _ => b'N',
    }
}

/// MCAST's `--text` table: a header and one line per reported cluster,
/// `cluster-N  name  start  stop  score(%g)  p(%.5g)  E(%.5g)  q(%.5g)
/// sequence` with 1-based inclusive positions and the matched letters in
/// upper case.
pub fn write_text<W: Write>(mut out: W, report: &Report, text: &McastText<'_>) -> io::Result<()> {
    out.write_all(
        b"pattern_name\tsequence_name\tstart\tstop\tscore\tp-value\tE-value\tq-value\tmatched_sequence\n",
    )?;
    let f = text.filter;
    let mut letters = Vec::new();
    for c in &report.clusters {
        if report.statistics
            && !(c.evalue <= f.evalue && c.pvalue <= f.pvalue && c.qvalue <= f.qvalue)
        {
            continue;
        }
        let k = &c.cluster;
        let s = k.sequence as usize;
        let origin = text.origins.and_then(|o| o.get(s).copied()).unwrap_or(0);
        write!(out, "cluster-{}\t", c.id)?;
        out.write_all(text.sequence_names.get(s))?;
        write!(
            out,
            "\t{}\t{}\t{}\t{}\t{}\t{}\t",
            origin + k.start as u64 + 1,
            origin + k.end as u64,
            format_g(k.score, 6),
            format_g(c.pvalue, 5),
            format_g(c.evalue, 5),
            format_g(c.qvalue, 5),
        )?;
        letters.clear();
        letters.extend(
            text.sequences.get(s)[k.start as usize..k.end as usize]
                .iter()
                .map(|&b| printable(b)),
        );
        out.write_all(&letters)?;
        out.write_all(b"\n")?;
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn mt19937_matches_the_reference() {
        let mut r = Mt19937::new(0);
        assert_eq!(
            [r.next(), r.next(), r.next()],
            [2357136044, 2546248239, 3071714933]
        );
        let mut r = Mt19937::new(0);
        assert_eq!(r.drand(), 0.548813502304256);
        assert_eq!(r.ldrand(), 0.5928446165554201);
        let mut r = Mt19937::knuth(4357);
        assert_eq!(r.ldrand(), 0.8173300600608487);
    }

    #[test]
    fn mcast_threshold_for_one_motif() {
        // 2 states, p 0.0005, max gap 60: MCAST's lex0 example.
        let d = cluster_threshold(0.0005, 60, 2, 1.0);
        assert!((d - 2.8677652116928782).abs() < 1e-12, "{d}");
        assert_eq!(meme_ev(0.0005, 2.0), 1.0 - (1.0 - 0.0005f64).powf(2.0));
    }

    #[test]
    fn qvalues_are_monotone_and_bounded() {
        let p: Vec<f64> = (0..300)
            .map(|i| ((i * 37) % 300) as f64 / 300.0 + 1e-4)
            .collect();
        let q = MemeQValue::default().adjust(&p);
        let mut pairs: Vec<(f64, f64)> = p.iter().copied().zip(q.iter().copied()).collect();
        pairs.sort_by(|a, b| a.0.total_cmp(&b.0));
        for w in pairs.windows(2) {
            assert!(w[0].1 <= w[1].1);
        }
        assert!(q.iter().all(|&v| (0.0..=1.0).contains(&v)));
        let small = MemeQValue::default().adjust(&[0.01, 0.02]);
        assert_eq!(small, vec![0.02, 0.02]);
    }

    #[test]
    fn pi_zero_handles_empty_samples_and_duplicate_cutoffs() {
        assert_eq!(MemeQValue::default().pi_zero(&[]), 1.0);
        assert_eq!(
            MemeQValue {
                lambdas: 0,
                ..MemeQValue::default()
            }
            .pi_zero(&[0.1, 0.5]),
            1.0
        );
        let sorted = [0.1, 0.2, 0.2, 0.5, 0.9];
        for lambda in [0.0, 0.2, 0.3, 0.9, 0.99] {
            let above = sorted.iter().filter(|&&p| p >= lambda).count();
            let expected = if above == 0 {
                f64::INFINITY
            } else {
                above as f64 / ((1.0 - lambda) * sorted.len() as f64)
            };
            assert_eq!(pi_zero_at(lambda, &sorted), expected);
        }
    }

    #[test]
    fn binned_tail_interpolates() {
        let samples: Vec<NullSample> = (0..3000)
            .map(|i| NullSample {
                excess: (i % 97) as f64 / 10.0,
                gc: 0.3 + (i as f64 / 3000.0) * 0.4,
                span: 30,
                segment: 100_000,
                order: i as u64,
            })
            .collect();
        let tail = GcBinnedExponential::mcast(50)
            .fit(&samples)
            .unwrap()
            .unwrap();
        let p = tail(5.0, 0.5);
        assert!(p > 0.0 && p < 1.0);
        assert!(tail(10.0, 0.5) < p);
        assert!(
            GcBinnedExponential::mcast(50)
                .fit(&samples[..10])
                .unwrap()
                .is_none()
        );
    }
}
