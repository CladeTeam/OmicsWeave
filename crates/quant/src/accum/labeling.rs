//! Metabolic labeling (SLAM-seq, scSLAM-seq, TimeLapse-seq, 4sU): base
//! conversions of each read against the reference.
//!
//! Definitions:
//! - A read's **transcript strand** follows the library: read 2 counts as
//!   read 1's strand, a reverse-stranded library flips it (an unstranded
//!   library is read as forward-stranded here; gene matrices use the
//!   assigned gene's strand instead). A conversion `T>C` is a `T` of the
//!   transcript read as `C`: on the genome, `T>C` for reads of a forward
//!   transcript and `A>G` for reads of a reverse one.
//! - Only aligned bases (`M` / `=` / `X`) count, and among them only those
//!   with a known reference base and read base (`A`, `C`, `G`, `T`), a
//!   base quality of at least `min_quality` (bases of a record without
//!   qualities always pass), and no known variant at the position (the
//!   task's BCF). The reference is the task's genome, or the read's `MD`
//!   tag without one; a read with neither is not counted.
//! - **Coverage** of a read: its counted bases whose transcript base is the
//!   conversion's source (`T`); **conversions**: those read as its target
//!   (`C`).
//!
//! Products of [`Labeling`]: the histogram of reads by (coverage,
//! conversions) and the mismatch table by sequencing cycle (every
//! transcript-sense reference base against every read base), with the
//! overall rate of each of the twelve substitutions. Gene matrices get a
//! `labeled` layer from [`crate::Counts::labeled`].

use std::collections::HashMap;

use omicsweave_core::{Error, Result};
use smallvec::SmallVec;

use super::coverage::Tally;
use super::{Accumulator, BatchInfo, Delivery, MemoryNeed, Needs, PlanInfo, RunCtx};
use crate::model::{Histogram, MismatchTable};
use crate::prelude::{Filter, ReadView, Strandedness, complement};
use crate::sink::json::Json;

/// A base conversion in transcript sense, and which bases are trusted.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Conversion {
    /// The transcript's base (`T`).
    pub from: u8,
    /// What it is read as (`C`).
    pub to: u8,
    /// Bases of lower quality are not counted.
    pub min_quality: u8,
}

impl Default for Conversion {
    /// `T>C` (4sU), bases of quality 20 or more.
    fn default() -> Conversion {
        Conversion {
            from: b'T',
            to: b'C',
            min_quality: 20,
        }
    }
}

impl Conversion {
    /// `from>to` (`T>C`, `G>A`, ...).
    pub fn parse(text: &str, min_quality: u8) -> Result<Conversion> {
        let b = text.as_bytes();
        let ok = |c: u8| matches!(c, b'A' | b'C' | b'G' | b'T');
        if b.len() != 3 || b[1] != b'>' || !ok(b[0]) || !ok(b[2]) || b[0] == b[2] {
            return Err(Error::invalid(
                "conversion",
                format!("{text:?} is not a base conversion such as T>C"),
            ));
        }
        Ok(Conversion {
            from: b[0],
            to: b[2],
            min_quality,
        })
    }

    pub fn name(&self) -> String {
        format!("{}>{}", self.from as char, self.to as char)
    }
}

/// The conversions of one read (or pair), for a transcript on either
/// strand: index 0 the forward strand (`T>C` on the genome), 1 the reverse
/// (`A>G`).
#[derive(Debug, Clone, Default, PartialEq, Eq)]
pub struct Converted {
    pub coverage: [u32; 2],
    /// Reference positions of the conversions, ascending.
    pub positions: [SmallVec<[u32; 4]>; 2],
}

impl Converted {
    /// Conversions for a transcript on the reverse strand or not.
    #[inline]
    pub fn count(&self, reverse: bool) -> u32 {
        self.positions[reverse as usize].len() as u32
    }

    /// The other mate's conversions joined in, a position converted in both
    /// counted once.
    pub fn join(&mut self, o: &Converted) {
        for s in 0..2 {
            self.coverage[s] += o.coverage[s];
            let mut all = std::mem::take(&mut self.positions[s]);
            all.extend_from_slice(&o.positions[s]);
            all.sort_unstable();
            all.dedup();
            self.positions[s] = all;
        }
    }
}

/// Whether an aligned base counts under `c`.
#[inline]
fn trusted(b: &crate::prelude::AlignedBase, c: &Conversion) -> bool {
    !b.masked
        && (b.qual == 255 || b.qual >= c.min_quality)
        && matches!(b.reference, b'A' | b'C' | b'G' | b'T')
        && matches!(b.read, b'A' | b'C' | b'G' | b'T')
}

/// `r`'s conversions under `c` for both transcript strands; `None` when the
/// read has no reference to compare with.
pub fn convert(r: &ReadView<'_>, c: &Conversion) -> Result<Option<Converted>> {
    let Some(bases) = r.aligned()? else {
        return Ok(None);
    };
    let pairs = [(c.from, c.to), (complement(c.from), complement(c.to))];
    let mut out = Converted::default();
    for b in bases {
        if !trusted(b, c) {
            continue;
        }
        for (s, &(from, to)) in pairs.iter().enumerate() {
            if b.reference == from {
                out.coverage[s] += 1;
                if b.read == to {
                    out.positions[s].push(b.ref_pos);
                }
            }
        }
    }
    Ok(Some(out))
}

/// Per-read conversion statistics of a library; see the module
/// documentation.
#[derive(Debug, Clone, PartialEq)]
pub struct Labeling {
    /// The output directory (`labeling`).
    pub name: String,
    pub filter: Filter,
    /// The library type that gives each read's transcript strand.
    pub strand: Strandedness,
    pub conversion: Conversion,
}

impl Default for Labeling {
    /// Primary, unique alignments that are not duplicates or QC failures;
    /// forward-stranded (QuantSeq FWD, 10x); `T>C`, quality 20.
    fn default() -> Labeling {
        use omicsweave_bam::flag;
        Labeling {
            name: "labeling".into(),
            filter: Filter {
                require: 0,
                exclude: Filter::DEFAULT_EXCLUDE | flag::SUPPLEMENTARY,
                min_mapq: 0,
                unique: true,
            },
            strand: Strandedness::Forward,
            conversion: Conversion::default(),
        }
    }
}

const BASES: [u8; 4] = *b"ACGT";

#[inline]
fn base_index(b: u8) -> Option<usize> {
    BASES.iter().position(|&x| x == b)
}

/// Mismatch counts by (mate, cycle, reference base, read base), grown as
/// longer reads come.
#[derive(Default)]
struct Cycles {
    /// `counts[mate][cycle * 16 + ref * 4 + read]`.
    counts: [Vec<u64>; 2],
}

impl Cycles {
    #[inline]
    fn add(&mut self, mate: usize, cycle: usize, r: usize, q: usize) {
        let v = &mut self.counts[mate];
        let at = cycle * 16 + r * 4 + q;
        if v.len() <= at {
            v.resize((cycle + 1) * 16, 0);
        }
        v[at] += 1;
    }

    fn merge(&mut self, o: Cycles) {
        for (m, theirs) in o.counts.into_iter().enumerate() {
            let mine = &mut self.counts[m];
            if mine.len() < theirs.len() {
                mine.resize(theirs.len(), 0);
            }
            for (a, b) in mine.iter_mut().zip(theirs) {
                *a += b;
            }
        }
    }
}

pub struct LabelingBatch {
    hist: HashMap<(u32, u32), u64>,
    cycles: Cycles,
    tally: Tally,
    no_reference: u64,
}

pub struct LabelingTotal {
    hist: HashMap<(u32, u32), u64>,
    cycles: Cycles,
    tally: Tally,
    no_reference: u64,
}

impl Accumulator for Labeling {
    type Batch = LabelingBatch;
    type Total<'p> = LabelingTotal;
    type Product = (Histogram, MismatchTable, Tally, u64);

    fn name(&self) -> &str {
        &self.name
    }

    fn needs(&self) -> Needs {
        Needs::CIGAR
            | Needs::TAGS
            | Needs::SEQUENCE
            | Needs::QUALITY
            | Needs::REFERENCE
            | Needs::VARIANTS
    }

    fn estimate(&self, _info: &PlanInfo) -> MemoryNeed {
        MemoryNeed {
            fixed: 1 << 20,
            per_batch: 256 << 10,
            min_buffer: 0,
            weight: 0,
        }
    }

    fn batch(&self, _info: BatchInfo) -> LabelingBatch {
        LabelingBatch {
            hist: HashMap::new(),
            cycles: Cycles::default(),
            tally: Tally::default(),
            no_reference: 0,
        }
    }

    fn observe(&self, b: &mut LabelingBatch, r: &ReadView<'_>) -> Result<()> {
        b.tally.seen += 1;
        if let Some(why) = self.filter.reject(r) {
            b.tally.reject(why);
            return Ok(());
        }
        let Some(bases) = r.aligned()? else {
            b.no_reference += 1;
            return Ok(());
        };
        let rev = r.reverse_for(self.strand);
        let c = &self.conversion;
        let (from, to) = if rev {
            (complement(c.from), complement(c.to))
        } else {
            (c.from, c.to)
        };
        let mate = r.is_read2() as usize;
        let (lead, tail) = r.hard_clips();
        let l_seq = r.rec.l_seq() as u32;
        let stored_reverse = r.is_reverse();
        let (mut cov, mut conv) = (0u32, 0u32);
        for x in bases {
            if !trusted(x, c) {
                continue;
            }
            if x.reference == from {
                cov += 1;
                conv += (x.read == to) as u32;
            }
            let (rb, qb) = if rev {
                (complement(x.reference), complement(x.read))
            } else {
                (x.reference, x.read)
            };
            let cycle = if stored_reverse {
                tail + (l_seq - 1 - x.qpos)
            } else {
                lead + x.qpos
            };
            let ri = base_index(rb).expect("trusted bases are ACGT");
            let qi = base_index(qb).expect("trusted bases are ACGT");
            b.cycles.add(mate, cycle as usize, ri, qi);
        }
        *b.hist.entry((cov, conv)).or_insert(0) += 1;
        Ok(())
    }

    fn total<'p>(&'p self, _rc: &RunCtx<'p>) -> Result<LabelingTotal> {
        Ok(LabelingTotal {
            hist: HashMap::new(),
            cycles: Cycles::default(),
            tally: Tally::default(),
            no_reference: 0,
        })
    }

    fn merge<'p>(
        &'p self,
        t: &mut LabelingTotal,
        b: LabelingBatch,
        _frontier: u64,
        _rc: &RunCtx<'p>,
    ) -> Result<()> {
        for (k, v) in b.hist {
            *t.hist.entry(k).or_insert(0) += v;
        }
        t.cycles.merge(b.cycles);
        t.tally.add(&b.tally);
        t.no_reference += b.no_reference;
        Ok(())
    }

    fn finish<'p>(
        &'p self,
        t: LabelingTotal,
        _rc: &RunCtx<'p>,
    ) -> Result<(Histogram, MismatchTable, Tally, u64)> {
        let mut cells: Vec<(u32, u32, u64)> =
            t.hist.into_iter().map(|((a, b), v)| (a, b, v)).collect();
        cells.sort_unstable();
        let c = &self.conversion;
        let hist = Histogram {
            axes: [
                format!("{}_coverage", c.from as char),
                format!("{}_conversions", c.name()),
            ],
            cells,
        };
        let cycles = t.cycles.counts.map(|v| {
            v.chunks_exact(16)
                .map(|c| {
                    let mut a = [[0u64; 4]; 4];
                    for (i, x) in c.iter().enumerate() {
                        a[i / 4][i % 4] = *x;
                    }
                    a
                })
                .collect()
        });
        Ok((hist, MismatchTable { cycles }, t.tally, t.no_reference))
    }

    fn deliver(
        &self,
        (hist, table, tally, no_reference): (Histogram, MismatchTable, Tally, u64),
        d: &mut Delivery<'_>,
    ) -> Result<Json> {
        let conv_file = format!("{}/conversions.tsv.gz", self.name);
        let path = d.path(&conv_file)?;
        crate::sink::table::write_histogram(&hist, &path, d)?;
        let mm_file = format!("{}/mismatches.tsv", self.name);
        let path = d.path(&mm_file)?;
        crate::sink::table::write_mismatches(&table, &path, d)?;
        let reads: u64 = hist.cells.iter().map(|c| c.2).sum();
        let cov: u64 = hist.cells.iter().map(|c| c.0 as u64 * c.2).sum();
        let conv: u64 = hist.cells.iter().map(|c| c.1 as u64 * c.2).sum();
        let min = |k: u32| -> u64 { hist.cells.iter().filter(|c| c.1 >= k).map(|c| c.2).sum() };
        let totals = table.totals();
        let mut rates = Vec::new();
        for (ri, &rb) in BASES.iter().enumerate() {
            let covered: u64 = totals[ri].iter().sum();
            for (qi, &qb) in BASES.iter().enumerate() {
                if qi == ri {
                    continue;
                }
                let rate = if covered == 0 {
                    0.0
                } else {
                    totals[ri][qi] as f64 / covered as f64
                };
                rates.push((format!("{}>{}", rb as char, qb as char), Json::Num(rate)));
            }
        }
        let mut fields = vec![
            ("product", Json::from("labeling")),
            ("conversion", Json::from(self.conversion.name())),
            ("reads", Json::from(reads)),
            ("coverage", Json::from(cov)),
            ("conversions", Json::from(conv)),
            (
                "conversion_rate",
                Json::Num(if cov == 0 {
                    0.0
                } else {
                    conv as f64 / cov as f64
                }),
            ),
            ("reads_with_1_or_more", Json::from(min(1))),
            ("reads_with_2_or_more", Json::from(min(2))),
            ("unassigned_no_reference", Json::from(no_reference)),
            ("mismatch_rates", Json::obj(rates)),
        ];
        fields.extend(tally.json());
        fields.push((
            "files",
            Json::Arr(vec![Json::from(conv_file), Json::from(mm_file)]),
        ));
        Ok(Json::obj(fields))
    }
}
