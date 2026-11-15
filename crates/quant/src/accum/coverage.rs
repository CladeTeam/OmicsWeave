//! Coverage tracks: per-base depth of reads, read ends or fragments.
//!
//! A batch turns each counted read into intervals and records their
//! boundaries as sorted events. The total keeps the events at or past the
//! merge frontier (later reads may still start there) and sweeps the rest
//! into runs of constant depth, so memory follows the depth at the
//! frontier, not the genome; finished runs go to a [`Spool`].

use omicsweave_core::{Error, Result};

use super::{Accumulator, BatchInfo, Delivery, MemoryNeed, Needs, PlanInfo, RunCtx};
use crate::model::{Run, Track};
use crate::prelude::{Filter, ReadView, Rejected, Strandedness};
use crate::sink::json::Json;
use crate::spool::Spool;

/// What each counted read adds to the track.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum CoverageUnit {
    /// Its aligned bases (`M` / `=` / `X`), as `samtools depth` counts
    /// them; with `deletions`, also the bases of `D` operations
    /// (`samtools depth -J`). Introns (`N`) never count.
    Bases { deletions: bool },
    /// Its whole span, introns and deletions included.
    Span,
    /// The base at its 5' end.
    FivePrime,
    /// The base at its 3' end.
    ThreePrime,
    /// Its pair's fragment, once, from the leftmost mate (`[pos,
    /// pos + |TLEN|)`); fragments longer than `max_len` (0: no limit) are
    /// dropped. Unpaired reads count their span; mates whose pair is
    /// counted elsewhere add nothing.
    Fragment { max_len: u32 },
}

/// One track for both strands, or one per strand.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum StrandSplit {
    Merged,
    /// Forward and reverse tracks by the strand of the transcript (or
    /// fragment) the read stands for, under this library type: read 2
    /// counts as read 1's strand, a reverse-stranded library flips it.
    Split(Strandedness),
}

/// How values are scaled when written.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum Scale {
    /// Raw counts.
    None,
    /// Counts per million counted reads (fragments for
    /// [`CoverageUnit::Fragment`]).
    Cpm,
    Factor(f64),
}

/// Coverage as bigWig: `<name>.bw`, or `<name>.forward.bw` and
/// `<name>.reverse.bw`.
#[derive(Debug, Clone, PartialEq)]
pub struct Coverage {
    pub name: String,
    pub unit: CoverageUnit,
    pub strands: StrandSplit,
    pub filter: Filter,
    pub scale: Scale,
    /// Also write each bigWig's `.sidecar`.
    pub sidecar: bool,
}

impl Default for Coverage {
    /// Aligned bases of reads passing [`Filter::new`], both strands, raw.
    fn default() -> Coverage {
        Coverage {
            name: "coverage".into(),
            unit: CoverageUnit::Bases { deletions: false },
            strands: StrandSplit::Merged,
            filter: Filter::new(),
            scale: Scale::None,
            sidecar: false,
        }
    }
}

/// Counts of records by what the filter made of them.
#[derive(Debug, Clone, Copy, Default)]
pub struct Tally {
    pub seen: u64,
    pub flag: u64,
    pub mapq: u64,
    pub multi: u64,
}

impl Tally {
    pub fn reject(&mut self, why: Rejected) {
        match why {
            Rejected::Flag => self.flag += 1,
            Rejected::Mapq => self.mapq += 1,
            Rejected::MultiMapping => self.multi += 1,
        }
    }

    pub fn add(&mut self, o: &Tally) {
        self.seen += o.seen;
        self.flag += o.flag;
        self.mapq += o.mapq;
        self.multi += o.multi;
    }

    pub fn json(&self) -> Vec<(&'static str, Json)> {
        vec![
            ("records", Json::from(self.seen)),
            ("filtered_flag", Json::from(self.flag)),
            ("filtered_mapq", Json::from(self.mapq)),
            ("filtered_multimapping", Json::from(self.multi)),
        ]
    }
}

pub struct CoverageBatch {
    /// Per strand: `(tid << 32 | pos) << 1 | is_end`, sorted at seal.
    events: [Vec<u64>; 2],
    counted: u64,
    tally: Tally,
}

pub struct CoverageTotal {
    pending: [Vec<u64>; 2],
    depth: [i64; 2],
    last: [Option<u64>; 2],
    open: [Option<Run>; 2],
    spools: [Spool<Run>; 2],
    counted: u64,
    tally: Tally,
    scratch: Vec<u64>,
}

#[inline]
fn event(tid: u32, pos: u32, end: bool) -> u64 {
    ((tid as u64) << 32 | pos as u64) << 1 | end as u64
}

impl Coverage {
    fn split(&self) -> bool {
        matches!(self.strands, StrandSplit::Split(_))
    }

    fn strand_of(&self, r: &ReadView<'_>) -> usize {
        match self.strands {
            StrandSplit::Merged => 0,
            StrandSplit::Split(lib) => r.reverse_for(lib) as usize,
        }
    }
}

fn add(events: &mut Vec<u64>, tid: u32, s: u32, e: u32) {
    if e > s {
        events.push(event(tid, s, false));
        events.push(event(tid, e, true));
    }
}

/// Closes `run` into the spool, or extends the open run with it.
fn emit(open: &mut Option<Run>, spool: &mut Spool<Run>, run: Run) -> Result<()> {
    if let Some(o) = open {
        if o.tid == run.tid && o.end == run.start && o.value == run.value {
            o.end = run.end;
            return Ok(());
        }
        spool.push(*o)?;
    }
    *open = Some(run);
    Ok(())
}

impl CoverageTotal {
    /// Sweeps the events of strand `s` before `limit` (a `tid << 32 | pos`
    /// key) out of `pending` and `next` into runs.
    fn sweep(&mut self, s: usize, next: &[u64], limit: u64) -> Result<()> {
        let pending = std::mem::take(&mut self.pending[s]);
        let merged = &mut self.scratch;
        merged.clear();
        merged.reserve(pending.len() + next.len());
        let (mut i, mut j) = (0, 0);
        while i < pending.len() && j < next.len() {
            if pending[i] <= next[j] {
                merged.push(pending[i]);
                i += 1;
            } else {
                merged.push(next[j]);
                j += 1;
            }
        }
        merged.extend_from_slice(&pending[i..]);
        merged.extend_from_slice(&next[j..]);
        let cut = merged.partition_point(|&e| e >> 1 < limit);
        let (mut depth, mut last, mut open) = (self.depth[s], self.last[s], self.open[s]);
        for &e in &merged[..cut] {
            let p = e >> 1;
            if let Some(l) = last
                && l != p
                && depth > 0
            {
                let run = Run {
                    tid: (l >> 32) as u32,
                    start: l as u32,
                    end: if p >> 32 == l >> 32 {
                        p as u32
                    } else {
                        return Err(Error::invalid(
                            "coverage",
                            "an interval runs past the end of its reference",
                        ));
                    },
                    value: depth.min(u32::MAX as i64) as u32,
                };
                emit(&mut open, &mut self.spools[s], run)?;
            }
            depth += if e & 1 == 1 { -1 } else { 1 };
            last = Some(p);
        }
        self.depth[s] = depth;
        self.last[s] = last;
        self.open[s] = open;
        let mut rest = pending;
        rest.clear();
        rest.extend_from_slice(&merged[cut..]);
        self.pending[s] = rest;
        Ok(())
    }
}

impl Accumulator for Coverage {
    type Batch = CoverageBatch;
    type Total<'p> = CoverageTotal;
    type Product = (Vec<Track>, CoverageTotalStats);

    fn name(&self) -> &str {
        &self.name
    }

    fn needs(&self) -> Needs {
        let mut n = Needs::CIGAR;
        if self.filter.unique {
            n |= Needs::TAGS;
        }
        n
    }

    fn estimate(&self, info: &PlanInfo) -> MemoryNeed {
        MemoryNeed {
            fixed: 0,
            // Two events of 8 bytes per block, a few blocks per record of
            // a few hundred bytes.
            per_batch: info.batch_bytes / 4,
            min_buffer: 2 << 20,
            weight: 1,
        }
    }

    fn batch(&self, _info: BatchInfo) -> CoverageBatch {
        CoverageBatch {
            events: [Vec::new(), Vec::new()],
            counted: 0,
            tally: Tally::default(),
        }
    }

    fn observe(&self, b: &mut CoverageBatch, r: &ReadView<'_>) -> Result<()> {
        b.tally.seen += 1;
        if let Some(why) = self.filter.reject(r) {
            b.tally.reject(why);
            return Ok(());
        }
        let s = self.strand_of(r);
        let ev = &mut b.events[s];
        let tid = r.tid();
        match self.unit {
            CoverageUnit::Bases { deletions } => {
                for &(bs, be) in r.blocks(deletions) {
                    add(ev, tid, bs, be);
                }
            }
            CoverageUnit::Span => add(ev, tid, r.pos(), r.end()),
            CoverageUnit::FivePrime => {
                let p = r.five_prime();
                add(ev, tid, p, p + 1)
            }
            CoverageUnit::ThreePrime => {
                let p = r.three_prime();
                add(ev, tid, p, p + 1)
            }
            CoverageUnit::Fragment { max_len } => {
                if !r.is_paired() {
                    add(ev, tid, r.pos(), r.end());
                } else if let Some(f) = r.fragment(max_len) {
                    add(ev, tid, f.start, f.end);
                } else {
                    return Ok(());
                }
            }
        }
        b.counted += 1;
        Ok(())
    }

    fn seal(&self, b: &mut CoverageBatch) -> Result<()> {
        for e in &mut b.events {
            e.sort_unstable();
        }
        Ok(())
    }

    fn total<'p>(&'p self, rc: &RunCtx<'p>) -> Result<CoverageTotal> {
        let share = if self.split() {
            rc.buffer / 2
        } else {
            rc.buffer
        };
        Ok(CoverageTotal {
            pending: [Vec::new(), Vec::new()],
            depth: [0; 2],
            last: [None; 2],
            open: [None; 2],
            spools: [Spool::new(share, &rc.tmp), Spool::new(share, &rc.tmp)],
            counted: 0,
            tally: Tally::default(),
            scratch: Vec::new(),
        })
    }

    fn merge<'p>(
        &'p self,
        t: &mut CoverageTotal,
        b: CoverageBatch,
        frontier: u64,
        _rc: &RunCtx<'p>,
    ) -> Result<()> {
        t.counted += b.counted;
        t.tally.add(&b.tally);
        // The frontier counts positions from 1; events from 0.
        let limit = frontier.saturating_sub(1);
        for s in 0..2 {
            if !b.events[s].is_empty() || !t.pending[s].is_empty() {
                t.sweep(s, &b.events[s], limit)?;
            }
        }
        Ok(())
    }

    fn end_input<'p>(&'p self, t: &mut CoverageTotal, _input: u32, _rc: &RunCtx<'p>) -> Result<()> {
        for s in 0..2 {
            t.sweep(s, &[], u64::MAX)?;
            t.last[s] = None;
        }
        Ok(())
    }

    fn finish<'p>(
        &'p self,
        mut t: CoverageTotal,
        rc: &RunCtx<'p>,
    ) -> Result<(Vec<Track>, CoverageTotalStats)> {
        let chroms: Vec<(String, u32)> = rc.inputs[0]
            .refs
            .iter()
            .map(|r| (r.0.clone(), r.1))
            .collect();
        let scale = match self.scale {
            Scale::None => 1.0,
            Scale::Cpm => 1e6 / t.counted.max(1) as f64,
            Scale::Factor(f) => f,
        };
        let mut tracks = Vec::new();
        let spools = std::mem::replace(
            &mut t.spools,
            [Spool::new(0, &rc.tmp), Spool::new(0, &rc.tmp)],
        );
        for (s, mut spool) in spools.into_iter().enumerate() {
            if let Some(o) = t.open[s].take() {
                spool.push(o)?;
            }
            if s == 1 && !self.split() {
                continue;
            }
            tracks.push(Track {
                chroms: chroms.clone(),
                strand: self.split().then_some(s == 1),
                runs: spool.finish()?,
                scale,
            });
        }
        Ok((
            tracks,
            CoverageTotalStats {
                counted: t.counted,
                tally: t.tally,
                scale,
            },
        ))
    }

    fn deliver(
        &self,
        (tracks, stats): (Vec<Track>, CoverageTotalStats),
        d: &mut Delivery<'_>,
    ) -> Result<Json> {
        let mut files = Vec::new();
        let mut spilled = 0;
        for track in tracks {
            let file = match track.strand {
                None => format!("{}.bw", self.name),
                Some(false) => format!("{}.forward.bw", self.name),
                Some(true) => format!("{}.reverse.bw", self.name),
            };
            spilled += track.runs.spilled();
            let path = d.path(&file)?;
            let runs = crate::sink::bigwig::write(track, &path, self.sidecar, d)?;
            files.push(Json::obj([
                ("file", Json::from(file)),
                ("runs", Json::from(runs)),
            ]));
        }
        let unit = match self.unit {
            CoverageUnit::Bases { deletions: false } => "bases",
            CoverageUnit::Bases { deletions: true } => "bases+deletions",
            CoverageUnit::Span => "span",
            CoverageUnit::FivePrime => "five-prime",
            CoverageUnit::ThreePrime => "three-prime",
            CoverageUnit::Fragment { .. } => "fragment",
        };
        let mut fields = vec![
            ("product", Json::from("coverage")),
            ("unit", Json::from(unit)),
            ("counted", Json::from(stats.counted)),
            ("scale", Json::from(stats.scale)),
            ("spilled_bytes", Json::from(spilled)),
        ];
        fields.extend(stats.tally.json());
        fields.push(("files", Json::Arr(files)));
        Ok(Json::obj(fields))
    }
}

/// What a coverage pass counted.
#[derive(Debug, Clone, Copy)]
pub struct CoverageTotalStats {
    pub counted: u64,
    pub tally: Tally,
    pub scale: f64,
}
