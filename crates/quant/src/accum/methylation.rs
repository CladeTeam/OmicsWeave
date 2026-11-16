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
//! Products: `<name>/cytosines.tsv.gz` (+ `.tbi` / `.csi`; chrom, start, end,
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
use omicsweave_bam::mods::{ModDecoder, ModSeq, ModTags};
use omicsweave_convert::stage::{Sorter, Spill, read_or_end};
use omicsweave_core::{Error, ModificationEntry, Result};
use omicsweave_gtf::Selection;
use smallvec::SmallVec;
use std::io::Read;

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
    /// Only the genes of this subset of the task's annotation.
    pub subset: Option<Arc<Selection>>,
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
            subset: None,
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

impl Spill for Cyt {
    type Key = (u64, u8);
    fn key(&self) -> Self::Key {
        self.order()
    }
    fn mem_size(&self) -> usize {
        0
    }
    fn encode(&self, out: &mut Vec<u8>) {
        out.extend_from_slice(&self.key.to_le_bytes());
        out.push(self.code);
        out.extend_from_slice(&self.meth.to_le_bytes());
        out.extend_from_slice(&self.unmeth.to_le_bytes());
    }
    fn decode(input: &mut impl Read) -> Result<Option<Self>> {
        let mut bytes = [0u8; 17];
        if !read_or_end(input, &mut bytes)? {
            return Ok(None);
        }
        Ok(Some(Self {
            key: u64::from_le_bytes(bytes[..8].try_into().unwrap()),
            code: bytes[8],
            meth: u32::from_le_bytes(bytes[9..13].try_into().unwrap()),
            unmeth: u32::from_le_bytes(bytes[13..].try_into().unwrap()),
        }))
    }
    const ABSORBS: bool = true;
    fn absorb(&mut self, next: Self) {
        self.meth += next.meth;
        self.unmeth += next.unmeth;
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
    mods: ModScratch,
}

/// Buffers of `MM` / `ML` decoding, reused from read to read.
#[derive(Default)]
struct ModScratch {
    decoder: ModDecoder,
    rows: Vec<ModificationEntry>,
}

pub struct MethylTotal<'p> {
    held: Vec<Cyt>,
    scratch: Vec<Cyt>,
    spool: Spool<Cytosine>,
    sorted: Option<Sorter<'p, Cyt>>,
    sorted_spilled: u64,
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
    fn calls(
        &self,
        r: &ReadView<'_>,
        out: &mut Vec<(Call, u32)>,
        mods: &mut ModScratch,
    ) -> Result<Outcome> {
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
            CallSource::Modifications => self.modifications(r, out, mods),
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
        if r.segments().iter().any(|&(rs, qs, len)| {
            rs.checked_add(len).is_none()
                || qs
                    .checked_add(len)
                    .is_none_or(|end| end as usize > xm.len() || end as usize > qual.len())
        }) {
            return Ok(Outcome::Malformed);
        }
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
        let genome = r.has_reference_context()?;
        let top = top_strand(r);
        let (target, converted) = if top { (b'C', b'T') } else { (b'G', b'A') };
        // Reference bases around a position: the genome's, else those the
        // MD tag gives under neighbouring aligned bases.
        let near = |i: usize, delta: i64| -> Result<u8> {
            let p = bases[i].ref_pos as i64 + delta;
            let Ok(p) = u32::try_from(p) else {
                return Ok(b'N');
            };
            if genome {
                return Ok(r.reference_context_base(p)?.unwrap_or(b'N'));
            }
            let j = i as i64 + delta;
            Ok(match bases.get(j.max(0) as usize) {
                Some(b) if j >= 0 && b.ref_pos == p => b.reference,
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

    fn modifications(
        &self,
        r: &ReadView<'_>,
        out: &mut Vec<(Call, u32)>,
        m: &mut ModScratch,
    ) -> Result<Outcome> {
        let (a, b) = r.rec.aux();
        let tags = match ModTags::find(&[a, b]) {
            Ok(Some(t)) if t.ml().is_some() => t,
            Ok(_) => return Ok(Outcome::NoSource),
            Err(_) => return Ok(Outcome::Malformed),
        };
        let l_seq = r.rec.l_seq();
        if let Some(n) = tags.mn()
            && n != l_seq as i64
        {
            return Ok(Outcome::Malformed);
        }
        if r.hard_clips() != (0, 0) && tags.mn().is_none() {
            return Ok(Outcome::Malformed);
        }
        let reverse = r.is_reverse();
        m.rows.clear();
        let seq = ModSeq::Packed {
            bytes: r.rec.packed_seq(),
            len: l_seq,
        };
        if m.decoder
            .decode(&tags, seq, reverse, &[b'm' as i32], &mut m.rows)
            .is_err()
        {
            return Ok(Outcome::Malformed);
        }
        let mut types = m
            .decoder
            .types()
            .iter()
            .filter(|t| t.canonical == b'C' && !t.minus && t.code == b'm' as i32);
        let Some(first) = types.next() else {
            return Ok(Outcome::Called);
        };
        let explicit = first.explicit;
        // Repeated types have no unique probability or skipped-base mode.
        if types.next().is_some() {
            return Ok(Outcome::Malformed);
        }
        let mut listed = m
            .rows
            .iter()
            .filter(|e| e.canonical == b'C' && e.flags & ModificationEntry::MINUS == 0);
        let mut next = if reverse {
            listed.next_back()
        } else {
            listed.next()
        };
        let packed = r.rec.packed_seq();
        let base = |q: usize| -> u8 {
            // SAFETY: callers pass q < l_seq; Record::parse validates that
            // packed contains ceil(l_seq / 2) bytes. A nibble is below 16.
            unsafe {
                let byte = *packed.get_unchecked(q / 2);
                *b"=ACMGRSVTWYHKDBN".get_unchecked(((byte >> (4 * (1 - q % 2))) & 15) as usize)
            }
        };
        // Bases of the read as sequenced.
        let original = |k: usize| -> u8 {
            if reverse {
                complement(base(l_seq - 1 - k))
            } else {
                base(k)
            }
        };
        let p = self.spec.min_probability;
        let qual = r.rec.qual();
        let mask = r.masked()?;
        let genome = r.has_reference_context()?;
        let mut segments = r.segments().iter();
        let mut segment = if reverse {
            segments.next_back()
        } else {
            segments.next()
        };
        // Walks the listed cytosines in the order the read was sequenced.
        for k in 0..l_seq {
            if original(k) != b'C' {
                continue;
            }
            let q = if reverse { l_seq - 1 - k } else { k } as u32;
            let prob = next.filter(|e| e.qpos as u32 == q).map(|e| {
                next = if reverse {
                    listed.next_back()
                } else {
                    listed.next()
                };
                e.probability
            });
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
            // Query positions are monotone, so each CIGAR segment is
            // visited once instead of binary-searching for each cytosine.
            while segment.is_some_and(|&(_, qs, len)| {
                if reverse {
                    q < qs
                } else {
                    q >= qs && q - qs >= len
                }
            }) {
                segment = if reverse {
                    segments.next_back()
                } else {
                    segments.next()
                };
            }
            let Some(&(rp, qs, len)) = segment else {
                continue;
            };
            let Some(delta) = q.checked_sub(qs).filter(|&d| d < len) else {
                continue;
            };
            let Some(pos) = rp.checked_add(delta) else {
                continue;
            };
            if !self.quality_ok(qual, q as usize) || mask.binary_search(&pos).is_ok() {
                continue;
            }
            let ctx = if genome {
                let at = |d: i64| -> Result<u8> {
                    let x = pos as i64 + d;
                    Ok(match u32::try_from(x) {
                        Ok(p) => r.reference_context_base(p)?.unwrap_or(b'N'),
                        Err(_) => b'N',
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
                        base(x as usize)
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

impl<'p> MethylTotal<'p> {
    /// Merges `next` (sorted, collapsed) into the held cytosines and spools
    /// those before `limit`.
    fn advance(&mut self, next: &[Cyt], limit: u64, rc: &RunCtx<'p>) -> Result<()> {
        let share = rc.buffer / (1 + 4 * self.sorters.len() as u64);
        let pending = self.held.len() - self.held.partition_point(|c| c.key < limit) + next.len()
            - next.partition_point(|c| c.key < limit);
        // Long alignments and distant mates can hold future calls over many
        // batches. Switch before repeatedly merging that growing suffix.
        if self.sorted.is_none()
            && pending.saturating_mul(size_of::<Cyt>()) as u64 > (share / 4).max(64 << 10)
        {
            self.spool.stream_to_disk()?;
            let mut sorter = Sorter::new(rc.pool, (share / 2).max(1), &rc.tmp);
            sorter.extend(std::mem::take(&mut self.held))?;
            self.scratch = Vec::new();
            self.sorted = Some(sorter);
        }
        if let Some(sorter) = &mut self.sorted {
            return sorter.extend(next.to_vec());
        }
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

    fn finish_sorted(&mut self) -> Result<()> {
        if let Some(sorter) = self.sorted.take() {
            self.sorted_spilled += sorter.spilled();
            let mut records = sorter.finish()?;
            let mut pending: Option<Cyt> = None;
            while let Some(c) = records.next_record()? {
                if let Some(last) = &mut pending
                    && last.order() == c.order()
                {
                    last.absorb(c);
                    continue;
                }
                if let Some(last) = pending.replace(c) {
                    self.spool.push(Cytosine {
                        tid: (last.key >> 32) as u32,
                        pos: last.key as u32,
                        code: last.code,
                        meth: last.meth,
                        unmeth: last.unmeth,
                    })?;
                }
            }
            if let Some(last) = pending {
                self.spool.push(Cytosine {
                    tid: (last.key >> 32) as u32,
                    pos: last.key as u32,
                    code: last.code,
                    meth: last.meth,
                    unmeth: last.unmeth,
                })?;
            }
        }
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
            mods: ModScratch::default(),
        }
    }

    fn observe(&self, b: &mut MethylBatch, r: &ReadView<'_>) -> Result<()> {
        b.tally.seen += 1;
        if let Some(why) = self.spec.filter.reject(r) {
            b.tally.reject(why);
            return Ok(());
        }
        let mut calls = std::mem::take(&mut b.calls);
        let outcome = self.calls(r, &mut calls, &mut b.mods)?;
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
            sorted: None,
            sorted_spilled: 0,
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
        t.advance(&cyt, limit, rc)?;
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
        rc: &RunCtx<'p>,
    ) -> Result<()> {
        let mut extra = Vec::new();
        let mut out = std::mem::take(&mut t.out);
        let status = &mut t.status;
        t.mates
            .drain(|h| self.emit_half(h, &mut extra, &mut out, status));
        t.waiting.clear();
        collapse_cyt(&mut extra);
        t.advance(&extra, u64::MAX, rc)?;
        t.finish_sorted()?;
        for (m, o) in out.iter_mut().enumerate() {
            t.sorters[m].extend(std::mem::take(o))?;
        }
        t.out = out;
        Ok(())
    }

    fn finish<'p>(&'p self, mut t: MethylTotal<'p>, rc: &RunCtx<'p>) -> Result<MethylProduct> {
        t.finish_sorted()?;
        let mut spilled = t.spool.spilled() + t.sorted_spilled;
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
        }
        if let Some(index) = &written.index {
            files.push(Json::from(
                index
                    .strip_prefix(d.dir)
                    .expect("index in output directory")
                    .to_string_lossy()
                    .into_owned(),
            ));
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
    fn repeated_modification_types_are_ambiguous() {
        use crate::assign::{BatchSpaces, Spaces};
        use crate::reference::BatchCtx;
        let spaces = Spaces::default();
        let ctx = BatchCtx {
            spaces: BatchSpaces::new(&spaces, 0),
            reference: None,
            variants: None,
        };
        let stage = MethylStage {
            spec: Methylation::default(),
            matrices: vec![],
            table_bytes: 0,
        };
        let mut raw = vec![0u8; 32];
        raw[8] = 2;
        raw[9] = 60;
        raw[12..14].copy_from_slice(&1u16.to_le_bytes());
        raw[16..20].copy_from_slice(&6i32.to_le_bytes());
        raw[20..24].copy_from_slice(&(-1i32).to_le_bytes());
        raw[24..28].copy_from_slice(&(-1i32).to_le_bytes());
        raw.extend_from_slice(b"r\0");
        raw.extend_from_slice(&(6u32 << 4).to_le_bytes());
        raw.extend_from_slice(&[0x24, 0x22, 0x42]);
        raw.extend_from_slice(&[30; 6]);
        raw.extend_from_slice(b"MMZC+m?,0;C+m.,1;\0MLBC");
        raw.extend_from_slice(&2u32.to_le_bytes());
        raw.extend_from_slice(&[255, 0]);
        let rec = omicsweave_bam::Record::parse(&raw, 1).unwrap();
        let view = ReadView::new(0, 0, rec, &ctx);
        let mut calls = Vec::new();
        let mut mods = ModScratch::default();
        assert_eq!(
            stage.modifications(&view, &mut calls, &mut mods).unwrap(),
            Outcome::Malformed
        );
        assert!(calls.is_empty());
    }

    #[test]
    fn repeated_spliced_calls_reuse_actual_fasta_windows() {
        use crate::assign::{BatchSpaces, Spaces};
        use crate::reference::{BatchCtx, RefWindow, Reference};
        use crate::source::Input;
        use std::cell::RefCell;
        let dir = tempfile::tempdir().unwrap();
        let path = dir.path().join("reference.fa");
        let len = 1_500_000u32;
        let mut fasta = b">chr1\n".to_vec();
        fasta.extend((0..len).map(|p| b"ACGT"[(p % 4) as usize]));
        fasta.push(b'\n');
        std::fs::write(&path, fasta).unwrap();
        std::fs::write(
            omicsweave_fasta::fai_path(&path),
            format!("chr1\t{len}\t6\t{len}\t{}\n", len + 1),
        )
        .unwrap();
        let inputs = [Input {
            path: dir.path().join("in.bam"),
            refs: vec![("chr1".into(), len)],
            read_groups: vec![],
            len: 0,
        }];
        let reference = Reference::open(&path, &inputs, None).unwrap();
        let spaces = Spaces::default();
        let ctx = BatchCtx {
            spaces: BatchSpaces::new(&spaces, 0),
            reference: Some(RefCell::new(RefWindow::new(&reference, 0))),
            variants: None,
        };
        let stage = MethylStage {
            spec: Methylation::default(),
            matrices: vec![],
            table_bytes: 0,
        };
        let mut raw = vec![0u8; 32];
        raw[8] = 2;
        raw[9] = 60;
        raw[12..14].copy_from_slice(&3u16.to_le_bytes());
        raw[16..20].copy_from_slice(&6i32.to_le_bytes());
        raw[20..24].copy_from_slice(&(-1i32).to_le_bytes());
        raw[24..28].copy_from_slice(&(-1i32).to_le_bytes());
        raw.extend_from_slice(b"r\0");
        for (n, op) in [(3u32, 0u32), (1_000_000, 3), (3, 0)] {
            raw.extend_from_slice(&((n << 4) | op).to_le_bytes());
        }
        let sequence = raw.len();
        raw.extend_from_slice(&[0; 3]);
        raw.extend_from_slice(&[30; 6]);
        let mut calls = Vec::new();
        for offset in 0..1000u32 {
            let pos = 100 + offset;
            raw[4..8].copy_from_slice(&(pos as i32).to_le_bytes());
            for i in 0..3 {
                raw[sequence + i] = ([1u8, 2, 4, 8][(pos as usize + 2 * i) % 4] << 4)
                    | [1u8, 2, 4, 8][(pos as usize + 2 * i + 1) % 4];
            }
            let rec = omicsweave_bam::Record::parse(&raw, 1).unwrap();
            let view = ReadView::new(0, offset as u64, rec, &ctx);
            calls.clear();
            assert!(matches!(
                stage.reference(&view, &mut calls).unwrap(),
                Outcome::Called
            ));
            let expected: Vec<_> = (0..6u32)
                .filter_map(|q| {
                    let p = pos + q + if q >= 3 { 1_000_000 } else { 0 };
                    (p % 4 == 1).then_some((p, 0, true, q))
                })
                .collect();
            assert_eq!(
                calls
                    .iter()
                    .map(|(c, q)| (c.pos, c.code, c.meth, *q))
                    .collect::<Vec<_>>(),
                expected
            );
        }
        let fetched = ctx.reference.as_ref().unwrap().borrow().fetched_bytes();
        assert_eq!(fetched, 1_048_576);
    }

    #[test]
    fn distant_calls_spill_without_remerging_the_whole_held_suffix() {
        use omicsweave_core::exec::Pool;
        use std::collections::BTreeMap;
        let dir = tempfile::tempdir().unwrap();
        let pool = Pool::new(2);
        let stage = MethylStage {
            spec: Methylation::default(),
            matrices: vec![],
            table_bytes: 0,
        };
        let small = RunCtx {
            pool: &pool,
            buffer: 1 << 20,
            tmp: dir.path().to_owned(),
            inputs: &[],
        };
        let large = RunCtx {
            pool: &pool,
            buffer: 128 << 20,
            tmp: dir.path().to_owned(),
            inputs: &[],
        };
        let mut a = stage.total(&small).unwrap();
        let mut b = stage.total(&large).unwrap();
        let mut want: BTreeMap<(u64, u8), (u32, u32)> = BTreeMap::new();
        for batch in 0..16u32 {
            let mut values: Vec<Cyt> = (0..2000u32)
                .map(|i| {
                    let pos = 100_000 + (i * 391 + batch * 127) % 20_000;
                    Cyt {
                        key: pos as u64,
                        code: if pos % 2 == 0 { 0 } else { 4 },
                        meth: i % 3,
                        unmeth: i % 2,
                    }
                })
                .collect();
            values.push(Cyt {
                key: 10,
                code: 0,
                meth: 1,
                unmeth: 0,
            });
            values.push(Cyt {
                key: 11,
                code: 4,
                meth: 0,
                unmeth: 1,
            });
            // A completed prefix is emitted only in the first batch.
            if batch != 0 {
                values.retain(|c| c.key >= 1000);
            }
            for c in &values {
                let e = want.entry(c.order()).or_default();
                e.0 += c.meth;
                e.1 += c.unmeth;
            }
            collapse_cyt(&mut values);
            a.advance(&values, 1000, &small).unwrap();
            b.advance(&values, 1000, &large).unwrap();
        }
        assert!(a.sorted.is_some());
        assert_eq!(a.held.capacity(), 0);
        assert_eq!(a.scratch.capacity(), 0);
        assert!(b.sorted.is_none());
        stage.end_input(&mut a, 0, &small).unwrap();
        stage.end_input(&mut b, 0, &large).unwrap();
        assert!(a.sorted.is_none());
        assert!(a.sorted_spilled > 0);
        let expected: Vec<_> = want
            .into_iter()
            .map(|((key, code), (meth, unmeth))| (key, code, meth, unmeth))
            .collect();
        let mut outputs = Vec::new();
        for spool in [a.spool, b.spool] {
            let mut got = Vec::new();
            spool
                .finish()
                .unwrap()
                .drain(|c| {
                    got.push((
                        ((c.tid as u64) << 32) | c.pos as u64,
                        c.code,
                        c.meth,
                        c.unmeth,
                    ));
                    Ok(())
                })
                .unwrap();
            assert_eq!(got, expected);
            outputs.push(got);
        }
        use omicsweave_convert::stage::Outputs;
        use omicsweave_core::{FileSource, Locus, Region, TrackRead};
        for merge in [false, true] {
            let mut tracks = Vec::new();
            let mut tables = Vec::new();
            for (mode, values) in outputs.iter().enumerate() {
                let path = dir.path().join(format!("{merge}-{mode}"));
                std::fs::create_dir_all(&path).unwrap();
                let mut files = Outputs::new();
                let mut delivery = Delivery {
                    outputs: &mut files,
                    dir: &path,
                    pool: &pool,
                    matrix: crate::task::MatrixFormat::Mtx,
                    tmp: dir.path(),
                    buffer: 4 << 20,
                };
                let cells = values
                    .iter()
                    .map(|&(key, code, meth, unmeth)| Cytosine {
                        tid: (key >> 32) as u32,
                        pos: key as u32,
                        code,
                        meth,
                        unmeth,
                    })
                    .collect();
                crate::sink::methylation::write(
                    Spooled::from_vec(cells),
                    &[("chr1".into(), 130_000)],
                    merge,
                    Some("calls.tsv.gz"),
                    Some(("beta.bw", "coverage.bw")),
                    &mut delivery,
                )
                .unwrap();
                files.commit().unwrap();
                tables.push(std::fs::read(path.join("calls.tsv.gz")).unwrap());
                let mut layers = Vec::new();
                for name in ["beta.bw", "coverage.bw"] {
                    let mut reader =
                        omicsweave_bbi::BbiReader::open(FileSource::open(path.join(name)).unwrap())
                            .unwrap();
                    let mut values = vec![0.0; 130_000];
                    reader
                        .values(&Locus::from(Region::new("chr1", 0, 130_000)), &mut values)
                        .unwrap();
                    layers.push(values.into_iter().map(f32::to_bits).collect::<Vec<_>>());
                }
                tracks.push(layers);
            }
            assert_eq!(tables[0], tables[1]);
            assert_eq!(tracks[0], tracks[1]);
        }
    }

    #[test]
    fn reference_calls_fetch_only_aligned_context_across_long_introns() {
        use crate::assign::{Query, SpaceId};
        use crate::prelude::{Context, Hits};
        use std::cell::RefCell;
        struct Genome(RefCell<Vec<(u32, u32)>>);
        impl Context for Genome {
            fn assign(&self, _: SpaceId, _: &ReadView<'_>, _: Query, _: Option<bool>) -> Hits {
                Hits::new()
            }
            fn assign_at(&self, _: SpaceId, _: &ReadView<'_>, _: u32, _: u32) -> Hits {
                Hits::new()
            }
            fn reference(
                &self,
                _: u32,
                start: u32,
                end: u32,
                with: &mut dyn FnMut(&[u8]),
            ) -> Result<bool> {
                self.0.borrow_mut().push((start, end));
                let bytes: Vec<_> = (start..end).map(|p| b"ACGT"[(p % 4) as usize]).collect();
                with(&bytes);
                Ok(true)
            }
            fn variants(&self, _: u32, _: u32, _: u32, _: &mut Vec<u32>) -> Result<()> {
                Ok(())
            }
        }
        let mut raw = vec![0u8; 32];
        raw[4..8].copy_from_slice(&100i32.to_le_bytes());
        raw[8] = 2;
        raw[9] = 60;
        raw[12..14].copy_from_slice(&3u16.to_le_bytes());
        raw[16..20].copy_from_slice(&6i32.to_le_bytes());
        raw[20..24].copy_from_slice(&(-1i32).to_le_bytes());
        raw[24..28].copy_from_slice(&(-1i32).to_le_bytes());
        raw.extend_from_slice(b"r\0");
        for (len, op) in [(3u32, 0u32), (1_000_000, 3), (3, 0)] {
            raw.extend_from_slice(&((len << 4) | op).to_le_bytes());
        }
        raw.extend_from_slice(&[0x12, 0x48, 0x12]);
        raw.extend_from_slice(&[30; 6]);
        let rec = omicsweave_bam::Record::parse(&raw, 1).unwrap();
        let genome = Genome(RefCell::new(Vec::new()));
        let view = ReadView::new(0, 0, rec, &genome);
        let stage = MethylStage {
            spec: Methylation::default(),
            matrices: vec![],
            table_bytes: 0,
        };
        let mut calls = Vec::new();
        assert!(matches!(
            stage.reference(&view, &mut calls).unwrap(),
            Outcome::Called
        ));
        assert_eq!(
            calls
                .iter()
                .map(|(c, q)| (c.pos, c.code, c.meth, *q))
                .collect::<Vec<_>>(),
            [(101, 0, true, 1), (1_000_105, 0, true, 5)]
        );
        let ranges = genome.0.borrow();
        assert_eq!(ranges.iter().map(|(s, e)| e - s).sum::<u32>(), 20);
        assert!(ranges.iter().all(|(s, e)| e - s <= 7));
    }

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
