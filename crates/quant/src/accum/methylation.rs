//! DNA methylation (WGBS, EM-seq, single-cell bisulfite, long reads with
//! base modifications): per-read cytosine calls, merged into a
//! per-cytosine table, beta and coverage tracks, region × row matrices and
//! the M-bias of each read position.
//!
//! Definitions:
//! - A **call** is one cytosine of a read: its reference position (the `C`,
//!   or the `G` for a cytosine of the reverse strand), strand, context and
//!   state (methylated or not). Calls come from, per record
//!   ([`CallSource::Auto`]):
//!   - the Bismark `XM` tag (`Z`/`z` CpG, `X`/`x` CHG, `H`/`h` CHH,
//!     uppercase methylated; `U`/`u` and `.` are not calls), the strand from
//!     `XG` (`CT` forward, `GA` reverse);
//!   - else the `MM` / `ML` tags (SAM base modifications): each `C` of a
//!     `C+m` list is methylated when its probability (`(ML + 0.5) / 256`)
//!     reaches `min_probability`, unmethylated when it is at most
//!     `1 - min_probability`, and not called in between; unlisted `C`s are
//!     unmethylated, or not called under the `?` mode; the strand is the
//!     read's;
//!   - else the read against the reference (bwa-meth / MethylDackel): a
//!     read of the original top strand (`YD:Z:f`, `XG:Z:CT`, or without
//!     them read 1 forward / read 2 reverse) calls each reference `C`
//!     read as `C` (methylated) or `T` (unmethylated); one of the bottom
//!     strand each reference `G` read as `G` or `A`.
//! - The **context** of a cytosine is that of the reference (or, without
//!   one, `MD`-derived bases and then the read's own bases for `MM` calls):
//!   followed by `G`, CpG; else `H` then `G`, CHG; else `HH`, CHH (on the
//!   reverse strand read on the complement, leftwards). A cytosine whose
//!   next bases are unknown is not called.
//! - Only bases of quality `min_quality` or more (all bases of a record
//!   without qualities), and not at a known variant of the task, are
//!   called.
//! - **Overlapping mates** count once: for a pair whose mates lie within
//!   64 kb on one reference, the calls of read 2 inside read 1's aligned
//!   span are dropped (Bismark's `--no_overlap`).
//! - **M-bias**: calls by context, mate and position in the read as
//!   sequenced (1-based, hard clips included), before overlapping calls
//!   are dropped.
//!
//! Products: `<name>/cytosines.tsv.gz` (+ `.tbi`; chrom, start, end,
//! context, strand, methylated, unmethylated, 0-based half-open; with
//! `merge_cpg` the two cytosines of a CpG are one line `[C, C + 2)` with
//! strand `.`), `<name>/beta.bw` and `<name>/coverage.bw` over the same
//! records, `<name>/mbias.tsv` (also in `summary.json`), and region ×
//! row matrices `<name>/<matrix>/` with the layers `methylated` and
//! `total` (calls in the region's cytosines).

use std::collections::BTreeMap;
use std::path::PathBuf;
use std::sync::Arc;

use omicsweave_bam::flag;
use omicsweave_convert::stage::Sorter;
use omicsweave_core::{Error, Result};
use smallvec::SmallVec;

use super::counts::{
    Entry, FeatureSource, Fold, GlobalRows, LocalRows, Reduce, RowKey, Rows, collapse,
};
use super::coverage::Tally;
use super::{Accumulator, BatchInfo, Delivery, MemoryNeed, Needs, PlanInfo, RunCtx};
use crate::assign::SpaceId;
use crate::model::{BarcodeTable, Cells, CountMatrix, Cytosine, FeatureList, FeatureTable, Mbias};
use crate::prelude::{Filter, MateKey, Mates, Pairs, ReadView, Strandedness, complement};
use crate::sink::json::Json;
use crate::spool::{Spool, Spooled};

/// CpG context bit.
pub const CPG: u8 = 1;
/// CHG context bit.
pub const CHG: u8 = 2;
/// CHH context bit.
pub const CHH: u8 = 4;
/// Context names by index (bit `1 << index`).
pub const CONTEXTS: [&str; 3] = ["CpG", "CHG", "CHH"];

/// Mates further apart than this cannot overlap; they are not paired.
const MATE_WINDOW: u32 = 1 << 16;

/// Where a read's calls come from.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum CallSource {
    /// Per record: `XM` when present, else `MM` / `ML`, else the reference.
    #[default]
    Auto,
    /// The Bismark `XM` / `XG` tags.
    Bismark,
    /// The read against the reference (bwa-meth, gemBS, BISCUIT).
    Reference,
    /// The `MM` / `ML` base-modification tags (5mC, `C+m`).
    Modifications,
}

impl CallSource {
    pub fn from_name(name: &str) -> Result<CallSource> {
        Ok(match name {
            "auto" => CallSource::Auto,
            "bismark" => CallSource::Bismark,
            "reference" => CallSource::Reference,
            "modifications" => CallSource::Modifications,
            _ => {
                return Err(Error::invalid(
                    "methylation",
                    format!(
                        "unknown call source {name:?}: use auto, bismark, reference or \
                         modifications"
                    ),
                ));
            }
        })
    }

    pub fn name(self) -> &'static str {
        match self {
            CallSource::Auto => "auto",
            CallSource::Bismark => "bismark",
            CallSource::Reference => "reference",
            CallSource::Modifications => "modifications",
        }
    }
}

/// Context bits from names (`CpG`, `CHG`, `CHH`, any case).
pub fn contexts_from_names<S: AsRef<str>>(names: &[S]) -> Result<u8> {
    let mut bits = 0;
    for n in names {
        let n = n.as_ref();
        match CONTEXTS.iter().position(|c| c.eq_ignore_ascii_case(n)) {
            Some(i) => bits |= 1 << i,
            None => {
                return Err(Error::invalid(
                    "methylation",
                    format!("unknown context {n:?}: use CpG, CHG or CHH"),
                ));
            }
        }
    }
    if bits == 0 {
        return Err(Error::invalid("methylation", "no context to count"));
    }
    Ok(bits)
}

/// A region × row matrix of calls: layers `methylated` and `total`.
#[derive(Debug, Clone)]
pub struct MethylMatrix {
    /// The directory under the product's (`genes`, `bins`, ...).
    pub name: String,
    pub features: FeatureSource,
    pub biotypes: Option<Vec<String>>,
    pub rows: RowKey,
    /// Only these row keys, rows in this order.
    pub whitelist: Option<PathBuf>,
    /// Context bits counted.
    pub contexts: u8,
}

impl MethylMatrix {
    /// Calls of `contexts` in `features` per input file.
    pub fn new(name: &str, features: FeatureSource) -> MethylMatrix {
        MethylMatrix {
            name: name.into(),
            features,
            biotypes: None,
            rows: RowKey::File,
            whitelist: None,
            contexts: CPG,
        }
    }
}

/// The methylation product; see the module documentation.
#[derive(Debug, Clone)]
pub struct Methylation {
    /// The output directory (`methylation`).
    pub name: String,
    pub source: CallSource,
    /// Context bits of the table and tracks.
    pub contexts: u8,
    /// One line per CpG (both strands) in the table and tracks.
    pub merge_cpg: bool,
    pub filter: Filter,
    pub min_quality: u8,
    /// For `MM` / `ML` calls.
    pub min_probability: f64,
    /// Count the calls of overlapping mates once.
    pub overlap: bool,
    /// Write the per-cytosine table.
    pub table: bool,
    /// Write the beta and coverage bigWigs.
    pub bigwig: bool,
    pub matrices: Vec<MethylMatrix>,
}

impl Default for Methylation {
    /// MethylDackel's defaults: primary, not QC-failed or duplicate
    /// alignments of MAPQ 10 or more, bases of quality 5 or more; CpG
    /// only; calls by `XM`, `MM` or the reference; table and tracks.
    fn default() -> Methylation {
        Methylation {
            name: "methylation".into(),
            source: CallSource::Auto,
            contexts: CPG,
            merge_cpg: false,
            filter: Filter {
                require: 0,
                exclude: Filter::DEFAULT_EXCLUDE | flag::SUPPLEMENTARY,
                min_mapq: 10,
                unique: false,
            },
            min_quality: 5,
            min_probability: 0.5,
            overlap: true,
            table: true,
            bigwig: true,
            matrices: Vec::new(),
        }
    }
}

/// One call of a read. `code` is the context index with bit 2 set on the
/// reverse strand.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) struct Call {
    pub pos: u32,
    pub code: u8,
    pub meth: bool,
}

/// A call waiting with its read for the mate, with its matrix columns as
/// `(matrix, column)`.
#[derive(Debug, Clone)]
struct HCall {
    call: Call,
    cols: SmallVec<[(u8, u32); 2]>,
}

/// A mate waiting for the other.
#[derive(Debug, Clone)]
pub struct Half {
    tid: u32,
    span: (u32, u32),
    read1: bool,
    rows: SmallVec<[u32; 2]>,
    calls: Vec<HCall>,
}

impl Half {
    /// The key of its first call (`tid << 32 | pos`), if any.
    fn first(&self) -> Option<u64> {
        self.calls
            .iter()
            .map(|c| (self.tid as u64) << 32 | c.call.pos as u64)
            .min()
    }
}

/// A per-cytosine count of a batch or the total: `key` is `tid << 32 |
/// pos`.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
struct Cyt {
    key: u64,
    code: u8,
    meth: u32,
    unmeth: u32,
}

impl Cyt {
    fn order(&self) -> (u64, u8) {
        (self.key, self.code)
    }
}

fn collapse_cyt(v: &mut Vec<Cyt>) {
    v.sort_unstable_by_key(Cyt::order);
    let mut w = 0;
    for i in 0..v.len() {
        if w > 0 && v[w - 1].order() == v[i].order() {
            v[w - 1].meth += v[i].meth;
            v[w - 1].unmeth += v[i].unmeth;
        } else {
            v[w] = v[i];
            w += 1;
        }
    }
    v.truncate(w);
}

/// Why reads gave no calls, and what the calls were.
#[derive(Debug, Clone, Copy, Default)]
pub struct MethylStatus {
    /// Reads whose calls were taken.
    pub reads: u64,
    /// Reads without `XM` / `MM` and without a reference or `MD` tag.
    pub no_source: u64,
    /// Reads whose `XM`, `MM` / `ML` or `MN` tags do not fit the read.
    pub malformed: u64,
    /// Calls of read 2 inside read 1's span, dropped.
    pub overlap_dropped: u64,
    /// Calls counted, methylated and not, by context.
    pub meth: [u64; 3],
    pub unmeth: [u64; 3],
}

impl MethylStatus {
    fn add(&mut self, o: &MethylStatus) {
        self.reads += o.reads;
        self.no_source += o.no_source;
        self.malformed += o.malformed;
        self.overlap_dropped += o.overlap_dropped;
        for c in 0..3 {
            self.meth[c] += o.meth[c];
            self.unmeth[c] += o.unmeth[c];
        }
    }
}

/// A matrix of the stage, resolved: its space and rows.
pub(crate) struct MatrixStage {
    pub spec: MethylMatrix,
    pub space: SpaceId,
    pub features: Arc<FeatureTable>,
    pub rows: Rows,
}

/// A [`Methylation`] product resolved against a task.
pub struct MethylStage {
    pub(crate) spec: Methylation,
    pub(crate) matrices: Vec<MatrixStage>,
    /// Fixed bytes: whitelists.
    pub(crate) table_bytes: u64,
}

pub struct MethylBatch {
    cyt: Vec<Cyt>,
    entries: Vec<Vec<Entry>>,
    local: Vec<LocalRows>,
    pairs: Pairs<Half>,
    dangling: Vec<(MateKey, Half)>,
    mbias: Mbias,
    tally: Tally,
    status: MethylStatus,
    no_row: Vec<u64>,
    calls: Vec<(Call, u32)>,
}

pub struct MethylTotal<'p> {
    held: Vec<Cyt>,
    scratch: Vec<Cyt>,
    spool: Spool<Cytosine>,
    sorters: Vec<Sorter<'p, Entry>>,
    rows: Vec<GlobalRows>,
    mates: Mates<Half>,
    /// First-call keys of the halves waiting in `mates`, counted.
    waiting: BTreeMap<u64, u32>,
    mbias: Mbias,
    tally: Tally,
    status: MethylStatus,
    no_row: Vec<u64>,
    out: Vec<Vec<Entry>>,
}

/// The product of a stage: cytosines, matrices and statistics.
pub struct MethylProduct {
    cytosines: Spooled<Cytosine>,
    chroms: Vec<(String, u32)>,
    matrices: Vec<CountMatrix>,
    mbias: Mbias,
    tally: Tally,
    status: MethylStatus,
    no_row: Vec<u64>,
    spilled: u64,
}

/// The context index of a cytosine followed by `n1` then `n2` on its own
/// strand.
#[inline]
fn context(n1: u8, n2: u8) -> Option<u8> {
    let h = |b: u8| matches!(b, b'A' | b'C' | b'T');
    match (n1, n2) {
        (b'G', _) => Some(0),
        (a, b'G') if h(a) => Some(1),
        (a, b) if h(a) && h(b) => Some(2),
        _ => None,
    }
}

/// Whether the read stands for the original top strand of a bisulfite
/// library: `YD` (bwa-meth), `XG` (Bismark), else read 1 forward or read
/// 2 reverse.
fn top_strand(r: &ReadView<'_>) -> bool {
    match r.tag_bytes(*b"YD") {
        Some(b"f") => return true,
        Some(b"r") => return false,
        _ => {}
    }
    match r.tag_bytes(*b"XG") {
        Some(b"CT") => return true,
        Some(b"GA") => return false,
        _ => {}
    }
    !r.reverse_for(Strandedness::Forward)
}

/// What call extraction made of a read.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum Outcome {
    Called,
    NoSource,
    Malformed,
}

impl MethylStage {
    fn quality_ok(&self, qual: &[u8], q: usize) -> bool {
        match qual.first() {
            Some(0xff) | None => true,
            _ => qual.get(q).is_some_and(|&v| v >= self.spec.min_quality),
        }
    }

    /// The read's calls into `out` as `(call, query position)`.
    fn calls(&self, r: &ReadView<'_>, out: &mut Vec<(Call, u32)>) -> Result<Outcome> {
        out.clear();
        let source = match self.spec.source {
            CallSource::Auto => {
                if r.tag(*b"XM").is_some() {
                    CallSource::Bismark
                } else if r.tag(*b"MM").is_some() || r.tag(*b"Mm").is_some() {
                    CallSource::Modifications
                } else {
                    CallSource::Reference
                }
            }
            s => s,
        };
        match source {
            CallSource::Bismark => self.bismark(r, out),
            CallSource::Modifications => self.modifications(r, out),
            _ => self.reference(r, out),
        }
    }

    fn bismark(&self, r: &ReadView<'_>, out: &mut Vec<(Call, u32)>) -> Result<Outcome> {
        let Some(xm) = r.tag_bytes(*b"XM") else {
            return Ok(Outcome::NoSource);
        };
        if xm.len() != r.rec.l_seq() {
            return Ok(Outcome::Malformed);
        }
        let reverse = match r.tag_bytes(*b"XG") {
            Some(b"CT") => false,
            Some(b"GA") => true,
            _ => !top_strand(r),
        };
        let qual = r.rec.qual();
        let mask = r.masked()?;
        for &(rs, qs, len) in r.segments() {
            for i in 0..len {
                let q = (qs + i) as usize;
                let c = xm[q];
                let ctx = match c.to_ascii_lowercase() {
                    b'z' => 0,
                    b'x' => 1,
                    b'h' => 2,
                    _ => continue,
                };
                let pos = rs + i;
                if !self.quality_ok(qual, q) || mask.binary_search(&pos).is_ok() {
                    continue;
                }
                out.push((
                    Call {
                        pos,
                        code: ctx | (reverse as u8) << 2,
                        meth: c.is_ascii_uppercase(),
                    },
                    q as u32,
                ));
            }
        }
        Ok(Outcome::Called)
    }

    fn reference(&self, r: &ReadView<'_>, out: &mut Vec<(Call, u32)>) -> Result<Outcome> {
        let Some(bases) = r.aligned()? else {
            return Ok(Outcome::NoSource);
        };
        let genome = r.reference_bases()?.is_some();
        let top = top_strand(r);
        let (target, converted) = if top { (b'C', b'T') } else { (b'G', b'A') };
        // Reference bases around a position: the genome's, else those the
        // MD tag gives under neighbouring aligned bases.
        let near = |i: usize, delta: i64| -> Result<u8> {
            let p = bases[i].ref_pos as i64 + delta;
            if p < 0 {
                return Ok(b'N');
            }
            if genome {
                return Ok(r.reference_base(p as u32)?.unwrap_or(b'N'));
            }
            let j = i as i64 + delta;
            Ok(match bases.get(j.max(0) as usize) {
                Some(b) if j >= 0 && b.ref_pos as i64 == p => b.reference,
                _ => b'N',
            })
        };
        for (i, b) in bases.iter().enumerate() {
            if b.reference != target || b.masked {
                continue;
            }
            if b.qual != 255 && b.qual < self.spec.min_quality {
                continue;
            }
            let meth = if b.read == target {
                true
            } else if b.read == converted {
                false
            } else {
                continue;
            };
            let ctx = if top {
                context(near(i, 1)?, near(i, 2)?)
            } else {
                context(complement(near(i, -1)?), complement(near(i, -2)?))
            };
            let Some(ctx) = ctx else { continue };
            out.push((
                Call {
                    pos: b.ref_pos,
                    code: ctx | (!top as u8) << 2,
                    meth,
                },
                b.qpos,
            ));
        }
        Ok(Outcome::Called)
    }

    fn modifications(&self, r: &ReadView<'_>, out: &mut Vec<(Call, u32)>) -> Result<Outcome> {
        let mm = r.tag_bytes(*b"MM").or_else(|| r.tag_bytes(*b"Mm"));
        let ml = match r.tag(*b"ML").or_else(|| r.tag(*b"Ml")) {
            Some(crate::prelude::Aux::Array(b'C', v)) => Some(v),
            Some(_) => return Ok(Outcome::Malformed),
            None => None,
        };
        let (Some(mm), Some(ml)) = (mm, ml) else {
            return Ok(Outcome::NoSource);
        };
        let l_seq = r.rec.l_seq();
        if let Some(n) = r.tag(*b"MN").and_then(|a| a.int())
            && n as usize != l_seq
        {
            return Ok(Outcome::Malformed);
        }
        if r.hard_clips() != (0, 0) && r.tag(*b"MN").is_none() {
            return Ok(Outcome::Malformed);
        }
        let reverse = r.is_reverse();
        // Bases of the read as sequenced.
        let original = |k: usize| -> u8 {
            if reverse {
                complement(r.base(l_seq - 1 - k))
            } else {
                r.base(k)
            }
        };
        // Probabilities of 5mC per C of the read as sequenced: `Some(v)`
        // listed with ML value v, `None` unlisted.
        let mut ml_at = 0usize;
        let mut found: Option<(Vec<Option<u8>>, bool)> = None;
        for group in mm.split(|&c| c == b';') {
            if group.is_empty() {
                continue;
            }
            let mut fields = group.split(|&c| c == b',');
            let head = fields.next().unwrap_or(&[]);
            if head.len() < 3 {
                return Ok(Outcome::Malformed);
            }
            let (base, strand) = (head[0], head[1]);
            let mut codes = &head[2..];
            let mut explicit = false;
            if let Some(&last) = codes.last()
                && (last == b'?' || last == b'.')
            {
                explicit = last == b'?';
                codes = &codes[..codes.len() - 1];
            }
            let n_codes = if codes.iter().all(u8::is_ascii_digit) {
                1
            } else {
                codes.len()
            };
            let m_index = codes.iter().position(|&c| c == b'm');
            let skips: Vec<usize> = match fields
                .map(|f| std::str::from_utf8(f).ok()?.parse::<usize>().ok())
                .collect::<Option<Vec<usize>>>()
            {
                Some(v) => v,
                None => return Ok(Outcome::Malformed),
            };
            let take = skips.len() * n_codes;
            if ml_at + take > ml.len() {
                return Ok(Outcome::Malformed);
            }
            if base == b'C'
                && strand == b'+'
                && found.is_none()
                && let Some(mi) = m_index
            {
                let mut probs = Vec::new();
                let mut c_seen = 0usize;
                let mut listed = skips.iter().enumerate();
                let mut next = listed.next();
                let mut skip_left = next.map(|(_, &s)| s);
                for k in 0..l_seq {
                    if original(k) != b'C' {
                        continue;
                    }
                    c_seen += 1;
                    match (next, skip_left) {
                        (Some((li, _)), Some(0)) => {
                            probs.push(Some(ml[ml_at + li * n_codes + mi]));
                            next = listed.next();
                            skip_left = next.map(|(_, &s)| s);
                        }
                        (Some(_), Some(s)) => {
                            probs.push(None);
                            skip_left = Some(s - 1);
                        }
                        _ => probs.push(None),
                    }
                }
                if next.is_some() || c_seen != probs.len() {
                    return Ok(Outcome::Malformed);
                }
                found = Some((probs, explicit));
            }
            ml_at += take;
        }
        let Some((probs, explicit)) = found else {
            return Ok(Outcome::Called);
        };
        let p = self.spec.min_probability;
        let qual = r.rec.qual();
        let mask = r.masked()?;
        let genome = r.reference_bases()?.is_some();
        let mut c_index = 0usize;
        for k in 0..l_seq {
            if original(k) != b'C' {
                continue;
            }
            let prob = probs[c_index];
            c_index += 1;
            let meth = match prob {
                Some(v) => {
                    let x = (v as f64 + 0.5) / 256.0;
                    if x >= p {
                        true
                    } else if x <= 1.0 - p {
                        false
                    } else {
                        continue;
                    }
                }
                None if explicit => continue,
                None => false,
            };
            let q = if reverse { l_seq - 1 - k } else { k } as u32;
            let Some(pos) = r.reference_at(q) else {
                continue;
            };
            if !self.quality_ok(qual, q as usize) || mask.binary_search(&pos).is_ok() {
                continue;
            }
            let ctx = if genome {
                let at = |d: i64| -> Result<u8> {
                    let x = pos as i64 + d;
                    Ok(if x < 0 {
                        b'N'
                    } else {
                        r.reference_base(x as u32)?.unwrap_or(b'N')
                    })
                };
                if reverse {
                    if at(0)? != b'G' {
                        continue;
                    }
                    context(complement(at(-1)?), complement(at(-2)?))
                } else {
                    if at(0)? != b'C' {
                        continue;
                    }
                    context(at(1)?, at(2)?)
                }
            } else {
                let read = |d: i64| -> u8 {
                    let x = q as i64 + d;
                    if x < 0 || x as usize >= l_seq {
                        b'N'
                    } else {
                        r.base(x as usize)
                    }
                };
                if reverse {
                    context(complement(read(-1)), complement(read(-2)))
                } else {
                    context(read(1), read(2))
                }
            };
            let Some(ctx) = ctx else { continue };
            out.push((
                Call {
                    pos,
                    code: ctx | (reverse as u8) << 2,
                    meth,
                },
                q,
            ));
        }
        Ok(Outcome::Called)
    }

    /// Whether `r` waits for its mate so that overlapping calls count
    /// once (both mates decide alike).
    fn pairs(&self, r: &ReadView<'_>) -> bool {
        let f = r.flag();
        self.spec.overlap
            && f & flag::PAIRED != 0
            && f & (flag::UNMAPPED | flag::MATE_UNMAPPED) == 0
            && r.rec.mate_tid == r.rec.tid
            && r.rec.pos.abs_diff(r.rec.mate_pos) < MATE_WINDOW
    }

    /// The calls of a unit (a read, or a pair resolved) into the batch or
    /// total outputs.
    fn emit(
        &self,
        tid: u32,
        rows: &[u32],
        calls: impl Iterator<Item = HCall>,
        cyt: &mut Vec<Cyt>,
        entries: &mut [Vec<Entry>],
        st: &mut MethylStatus,
    ) {
        for h in calls {
            let c = h.call;
            let ctx = (c.code & 3) as usize;
            if c.meth {
                st.meth[ctx] += 1;
            } else {
                st.unmeth[ctx] += 1;
            }
            if self.spec.contexts & (1 << ctx) != 0 {
                cyt.push(Cyt {
                    key: (tid as u64) << 32 | c.pos as u64,
                    code: c.code,
                    meth: c.meth as u32,
                    unmeth: !c.meth as u32,
                });
            }
            for &(m, col) in &h.cols {
                let row = rows[m as usize];
                let e = &mut entries[m as usize];
                if c.meth {
                    e.push(Entry {
                        row,
                        col,
                        layer: 0,
                        sub: 0,
                        value: 1,
                    });
                }
                e.push(Entry {
                    row,
                    col,
                    layer: 1,
                    sub: 0,
                    value: 1,
                });
            }
        }
    }

    /// Both mates' calls, read 2's inside read 1's span dropped.
    fn resolve(
        &self,
        a: Half,
        b: Half,
        cyt: &mut Vec<Cyt>,
        entries: &mut [Vec<Entry>],
        st: &mut MethylStatus,
    ) {
        let (one, two) = if b.read1 && !a.read1 { (b, a) } else { (a, b) };
        let (s, e) = one.span;
        let before = two.calls.len();
        let kept: Vec<HCall> = two
            .calls
            .into_iter()
            .filter(|h| h.call.pos < s || h.call.pos >= e)
            .collect();
        st.overlap_dropped += (before - kept.len()) as u64;
        self.emit(one.tid, &one.rows, one.calls.into_iter(), cyt, entries, st);
        self.emit(two.tid, &two.rows, kept.into_iter(), cyt, entries, st);
    }

    fn emit_half(
        &self,
        h: Half,
        cyt: &mut Vec<Cyt>,
        entries: &mut [Vec<Entry>],
        st: &mut MethylStatus,
    ) {
        self.emit(h.tid, &h.rows, h.calls.into_iter(), cyt, entries, st);
    }
}

/// Folds the halves' first-call keys in and out of the waiting multiset.
fn wait(waiting: &mut BTreeMap<u64, u32>, key: Option<u64>) {
    if let Some(k) = key {
        *waiting.entry(k).or_insert(0) += 1;
    }
}

fn unwait(waiting: &mut BTreeMap<u64, u32>, key: Option<u64>) {
    if let Some(k) = key
        && let Some(n) = waiting.get_mut(&k)
    {
        *n -= 1;
        if *n == 0 {
            waiting.remove(&k);
        }
    }
}

impl MethylTotal<'_> {
    /// Merges `next` (sorted, collapsed) into the held cytosines and spools
    /// those before `limit`.
    fn advance(&mut self, next: &[Cyt], limit: u64) -> Result<()> {
        let held = std::mem::take(&mut self.held);
        let mut all = std::mem::take(&mut self.scratch);
        all.clear();
        let (mut i, mut j) = (0, 0);
        while i < held.len() || j < next.len() {
            let take_held =
                j == next.len() || (i < held.len() && held[i].order() <= next[j].order());
            let c = if take_held {
                i += 1;
                held[i - 1]
            } else {
                j += 1;
                next[j - 1]
            };
            match all.last_mut() {
                Some(last) if last.order() == c.order() => {
                    last.meth += c.meth;
                    last.unmeth += c.unmeth;
                }
                _ => all.push(c),
            }
        }
        let cut = all.partition_point(|c| c.key < limit);
        for c in &all[..cut] {
            self.spool.push(Cytosine {
                tid: (c.key >> 32) as u32,
                pos: c.key as u32,
                code: c.code,
                meth: c.meth,
                unmeth: c.unmeth,
            })?;
        }
        let mut rest = held;
        rest.clear();
        rest.extend_from_slice(&all[cut..]);
        self.held = rest;
        self.scratch = all;
        Ok(())
    }
}

impl Accumulator for MethylStage {
    type Batch = MethylBatch;
    type Total<'p> = MethylTotal<'p>;
    type Product = MethylProduct;

    fn name(&self) -> &str {
        &self.spec.name
    }

    fn needs(&self) -> Needs {
        let mut n = Needs::CIGAR
            | Needs::TAGS
            | Needs::SEQUENCE
            | Needs::QUALITY
            | Needs::REFERENCE
            | Needs::VARIANTS;
        if self.spec.overlap {
            n |= Needs::MATES;
        }
        if !self.matrices.is_empty() {
            n |= Needs::ASSIGN;
        }
        n
    }

    fn estimate(&self, info: &PlanInfo) -> MemoryNeed {
        let m = self.matrices.len() as u64;
        MemoryNeed {
            fixed: self.table_bytes
                + self
                    .matrices
                    .iter()
                    .map(|x| x.features.heap_bytes())
                    .sum::<u64>(),
            // Calls of a few bytes for most bases of a record.
            per_batch: info.batch_bytes / 2 + m * info.batch_bytes / 4,
            min_buffer: (2 << 20) + m * (1 << 20),
            weight: 1 + 4 * m as u32,
        }
    }

    fn batch(&self, _info: BatchInfo) -> MethylBatch {
        let m = self.matrices.len();
        MethylBatch {
            cyt: Vec::new(),
            entries: (0..m).map(|_| Vec::new()).collect(),
            local: (0..m).map(|_| LocalRows::default()).collect(),
            pairs: Pairs::default(),
            dangling: Vec::new(),
            mbias: Mbias::default(),
            tally: Tally::default(),
            status: MethylStatus::default(),
            no_row: vec![0; m],
            calls: Vec::new(),
        }
    }

    fn observe(&self, b: &mut MethylBatch, r: &ReadView<'_>) -> Result<()> {
        b.tally.seen += 1;
        if let Some(why) = self.spec.filter.reject(r) {
            b.tally.reject(why);
            return Ok(());
        }
        let mut calls = std::mem::take(&mut b.calls);
        let outcome = self.calls(r, &mut calls)?;
        match outcome {
            Outcome::NoSource => b.status.no_source += 1,
            Outcome::Malformed => b.status.malformed += 1,
            Outcome::Called => b.status.reads += 1,
        }
        // M-bias: every call by position in the read as sequenced.
        let mate = r.is_read2() as usize;
        let (lead, tail) = r.hard_clips();
        let l_seq = r.rec.l_seq() as u32;
        for &(c, q) in &calls {
            let cycle = if r.is_reverse() {
                tail + (l_seq - 1 - q)
            } else {
                lead + q
            };
            b.mbias
                .add(mate, (c.code & 3) as usize, cycle as usize, c.meth);
        }
        let mut rows: SmallVec<[u32; 2]> = SmallVec::new();
        for (m, st) in self.matrices.iter().enumerate() {
            match st.rows.of(r, &mut b.local[m]) {
                Some(row) => rows.push(row),
                None => {
                    rows.push(u32::MAX);
                    b.no_row[m] += 1;
                }
            }
        }
        let hcalls: Vec<HCall> = calls
            .iter()
            .map(|&(c, _)| {
                let mut cols = SmallVec::new();
                let bit = 1 << (c.code & 3);
                for (m, st) in self.matrices.iter().enumerate() {
                    if rows[m] != u32::MAX && st.spec.contexts & bit != 0 {
                        for col in r.assign_at(st.space, c.pos, c.pos + 1) {
                            cols.push((m as u8, col));
                        }
                    }
                }
                HCall { call: c, cols }
            })
            .collect();
        b.calls = calls;
        let half = Half {
            tid: r.tid(),
            span: (r.pos(), r.end()),
            read1: r.flag() & flag::READ1 != 0,
            rows,
            calls: hcalls,
        };
        if self.pairs(r) {
            if let Some((a, z)) = b.pairs.offer(r, half) {
                self.resolve(a, z, &mut b.cyt, &mut b.entries, &mut b.status);
            }
        } else {
            self.emit_half(half, &mut b.cyt, &mut b.entries, &mut b.status);
        }
        Ok(())
    }

    fn seal(&self, b: &mut MethylBatch) -> Result<()> {
        collapse_cyt(&mut b.cyt);
        for e in &mut b.entries {
            collapse(e, Fold::direct(Reduce::Sum));
        }
        b.dangling = std::mem::take(&mut b.pairs).into_dangling();
        b.calls = Vec::new();
        Ok(())
    }

    fn total<'p>(&'p self, rc: &RunCtx<'p>) -> Result<MethylTotal<'p>> {
        let m = self.matrices.len() as u64;
        // The spool takes a share like one matrix's.
        let share = rc.buffer / (1 + 4 * m);
        Ok(MethylTotal {
            held: Vec::new(),
            scratch: Vec::new(),
            spool: Spool::new(share, &rc.tmp),
            sorters: (0..m)
                .map(|_| Sorter::new(rc.pool, share * 4, &rc.tmp))
                .collect(),
            rows: (0..m).map(|_| GlobalRows::default()).collect(),
            mates: Mates::new(std::mem::size_of::<Half>() + 16 * std::mem::size_of::<HCall>()),
            waiting: BTreeMap::new(),
            mbias: Mbias::default(),
            tally: Tally::default(),
            status: MethylStatus::default(),
            no_row: vec![0; m as usize],
            out: (0..m).map(|_| Vec::new()).collect(),
        })
    }

    fn merge<'p>(
        &'p self,
        t: &mut MethylTotal<'p>,
        b: MethylBatch,
        frontier: u64,
        rc: &RunCtx<'p>,
    ) -> Result<()> {
        let MethylBatch {
            mut cyt,
            mut entries,
            local,
            dangling,
            mbias,
            tally,
            status,
            no_row,
            ..
        } = b;
        t.tally.add(&tally);
        t.status.add(&status);
        t.mbias.merge(&mbias);
        for (a, n) in t.no_row.iter_mut().zip(&no_row) {
            *a += n;
        }
        let maps: Vec<Option<Vec<u32>>> = local
            .into_iter()
            .enumerate()
            .map(|(m, l)| self.matrices[m].rows.local().then(|| t.rows[m].remap(l)))
            .collect();
        for (m, map) in maps.iter().enumerate() {
            if let Some(map) = map {
                for e in &mut entries[m] {
                    e.row = map[e.row as usize];
                }
            }
        }
        let remap = |h: &mut Half| {
            for (m, map) in maps.iter().enumerate() {
                if let Some(map) = map
                    && h.rows[m] != u32::MAX
                {
                    h.rows[m] = map[h.rows[m] as usize];
                }
            }
        };
        let mut extra = Vec::new();
        let mut out = std::mem::take(&mut t.out);
        for (key, mut half) in dangling {
            remap(&mut half);
            let first = half.first();
            match t.mates.offer(key, half) {
                Some((a, z)) => {
                    unwait(&mut t.waiting, a.first());
                    self.resolve(a, z, &mut extra, &mut out, &mut t.status);
                }
                None => wait(&mut t.waiting, first),
            }
        }
        {
            let (waiting, status) = (&mut t.waiting, &mut t.status);
            t.mates.expire(frontier, |h| {
                unwait(waiting, h.first());
                self.emit_half(h, &mut extra, &mut out, status);
            });
        }
        if t.mates.bytes() + t.rows.iter().map(GlobalRows::bytes).sum::<u64>()
            > rc.buffer.max(64 << 20)
        {
            return Err(Error::invalid(
                "budget",
                format!(
                    "{}: {} mates waiting for their pair and the row keys outgrow the memory \
                     budget; give a whitelist or raise the memory",
                    self.spec.name,
                    t.mates.len()
                ),
            ));
        }
        if !extra.is_empty() {
            cyt.extend(extra);
            collapse_cyt(&mut cyt);
        }
        // Later records start at or after the frontier's position, and so
        // do their calls; waiting halves hold calls from their first one.
        let mut limit = frontier.saturating_sub(1);
        if let Some((&k, _)) = t.waiting.iter().next() {
            limit = limit.min(k);
        }
        t.advance(&cyt, limit)?;
        for (m, mut e) in entries.into_iter().enumerate() {
            e.append(&mut out[m]);
            t.sorters[m].extend(e)?;
        }
        t.out = out;
        Ok(())
    }

    fn end_input<'p>(
        &'p self,
        t: &mut MethylTotal<'p>,
        _input: u32,
        _rc: &RunCtx<'p>,
    ) -> Result<()> {
        let mut extra = Vec::new();
        let mut out = std::mem::take(&mut t.out);
        let status = &mut t.status;
        t.mates
            .drain(|h| self.emit_half(h, &mut extra, &mut out, status));
        t.waiting.clear();
        collapse_cyt(&mut extra);
        t.advance(&extra, u64::MAX)?;
        for (m, o) in out.iter_mut().enumerate() {
            t.sorters[m].extend(std::mem::take(o))?;
        }
        t.out = out;
        Ok(())
    }

    fn finish<'p>(&'p self, t: MethylTotal<'p>, rc: &RunCtx<'p>) -> Result<MethylProduct> {
        let mut spilled = t.spool.spilled();
        let mut matrices = Vec::new();
        for ((st, sorter), rows) in self.matrices.iter().zip(t.sorters).zip(t.rows) {
            spilled += sorter.spilled();
            let names = match st.rows.fixed_names() {
                Some(n) => n.to_vec(),
                None => rows.names,
            };
            matrices.push(CountMatrix {
                rows: BarcodeTable { names },
                features: Arc::clone(&st.features) as Arc<dyn FeatureList>,
                layers: vec!["methylated".into(), "total".into()],
                named_layers: true,
                cells: Cells {
                    merge: Some(sorter.finish()?),
                    fold: Fold::direct(Reduce::Sum),
                    head: None,
                },
            });
        }
        Ok(MethylProduct {
            cytosines: t.spool.finish()?,
            chroms: rc.inputs[0].refs.clone(),
            matrices,
            mbias: t.mbias,
            tally: t.tally,
            status: t.status,
            no_row: t.no_row,
            spilled,
        })
    }

    fn deliver(&self, p: MethylProduct, d: &mut Delivery<'_>) -> Result<Json> {
        let name = &self.spec.name;
        let mut files = Vec::new();
        let tracks =
            (self.spec.bigwig).then(|| (format!("{name}/beta.bw"), format!("{name}/coverage.bw")));
        let table = self.spec.table.then(|| format!("{name}/cytosines.tsv.gz"));
        let written = crate::sink::methylation::write(
            p.cytosines,
            &p.chroms,
            self.spec.merge_cpg,
            table.as_deref(),
            tracks.as_ref().map(|(a, b)| (a.as_str(), b.as_str())),
            d,
        )?;
        if let Some(t) = &table {
            files.push(Json::from(t.clone()));
            files.push(Json::from(format!("{t}.tbi")));
        }
        if let Some((a, b)) = tracks {
            files.push(Json::from(a));
            files.push(Json::from(b));
        }
        let mbias_file = format!("{name}/mbias.tsv");
        let path = d.path(&mbias_file)?;
        crate::sink::table::write_mbias(&p.mbias, &path, d)?;
        files.push(Json::from(mbias_file));
        let mut matrices = Vec::new();
        for ((st, m), no_row) in self.matrices.iter().zip(p.matrices).zip(&p.no_row) {
            let (n_rows, n_features) = (m.rows.len(), m.features.len());
            let w = crate::sink::matrix::write(m, &format!("{name}/{}", st.spec.name), false, d)?;
            matrices.push(Json::obj([
                ("name", Json::from(st.spec.name.clone())),
                ("rows", Json::from(n_rows)),
                ("features", Json::from(n_features)),
                ("nonzero", Json::from(w.nonzero)),
                ("unassigned_no_row_key", Json::from(*no_row)),
                (
                    "files",
                    Json::Arr(w.files.into_iter().map(Json::from).collect()),
                ),
            ]));
        }
        let s = &p.status;
        let by_context =
            |v: &[u64; 3]| Json::obj(CONTEXTS.iter().zip(v).map(|(c, &n)| (*c, Json::from(n))));
        let contexts: Vec<Json> = (0..3)
            .filter(|i| self.spec.contexts & (1 << i) != 0)
            .map(|i| Json::from(CONTEXTS[i]))
            .collect();
        let mut fields = vec![
            ("product", Json::from("methylation")),
            ("source", Json::from(self.spec.source.name())),
            ("contexts", Json::Arr(contexts)),
            ("reads_called", Json::from(s.reads)),
            ("unassigned_no_source", Json::from(s.no_source)),
            ("unassigned_malformed_tags", Json::from(s.malformed)),
            ("overlap_calls_dropped", Json::from(s.overlap_dropped)),
            ("methylated_calls", by_context(&s.meth)),
            ("unmethylated_calls", by_context(&s.unmeth)),
            ("cytosines", Json::from(written.cytosines)),
            ("table_lines", Json::from(written.lines)),
            ("spilled_bytes", Json::from(p.spilled)),
            ("mbias", p.mbias.json()),
        ];
        if !matrices.is_empty() {
            fields.push(("matrices", Json::Arr(matrices)));
        }
        fields.extend(p.tally.json());
        fields.push(("files", Json::Arr(files)));
        Ok(Json::obj(fields))
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn contexts_follow_the_next_bases() {
        assert_eq!(context(b'G', b'A'), Some(0));
        assert_eq!(context(b'A', b'G'), Some(1));
        assert_eq!(context(b'T', b'C'), Some(2));
        assert_eq!(context(b'N', b'G'), None);
        assert_eq!(context(b'C', b'N'), None);
        assert_eq!(contexts_from_names(&["cpg", "CHH"]).unwrap(), CPG | CHH);
        assert!(contexts_from_names(&["CWG"]).is_err());
    }
}
