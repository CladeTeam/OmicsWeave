//! Per-cytosine methylation → `cytosines.tsv.gz` (BGZF, tabix-indexed as
//! BED: chrom, start, end, context, strand, methylated, unmethylated;
//! 0-based half-open) and the `beta.bw` / `coverage.bw` tracks of the same
//! lines, written in one pass over the cytosines.

use std::fs::File;
use std::io::BufWriter;
use std::path::PathBuf;

use omicsweave_bbi::{BigWigWriter, Chrom, WriteOptions};
use omicsweave_core::exec::Exec;
use omicsweave_core::{Error, Result};

use super::fragments::IndexedBed;
use super::push_u64;
use crate::accum::Delivery;
use crate::accum::methylation::CONTEXTS;
use crate::model::Cytosine;
use crate::spool::Spooled;

/// What the sink wrote.
pub struct Written {
    /// Cytosines read.
    pub cytosines: u64,
    /// Lines of the table (CpGs merged).
    pub lines: u64,
    /// The table's index, when a table was requested.
    pub index: Option<PathBuf>,
}

/// One line: `[start, end)` with its strand (`+`, `-` or `.`).
struct Line {
    tid: u32,
    start: u32,
    end: u32,
    context: usize,
    strand: u8,
    meth: u64,
    unmeth: u64,
}

/// Merges the two cytosines of each CpG into one line when `merge`.
struct Lines {
    merge: bool,
    /// A forward CpG waiting for its reverse partner.
    pending: Option<Cytosine>,
}

fn single(c: &Cytosine) -> Result<Line> {
    let end = c
        .pos
        .checked_add(1)
        .ok_or_else(|| Error::invalid("methylation output", "cytosine end exceeds u32"))?;
    Ok(Line {
        tid: c.tid,
        start: c.pos,
        end,
        context: c.context(),
        strand: if c.reverse() { b'-' } else { b'+' },
        meth: u64::from(c.meth),
        unmeth: u64::from(c.unmeth),
    })
}

fn cpg(start: u32, tid: u32, meth: u64, unmeth: u64) -> Result<Line> {
    let end = start
        .checked_add(2)
        .ok_or_else(|| Error::invalid("methylation output", "CpG end exceeds u32"))?;
    Ok(Line {
        tid,
        start,
        end,
        context: 0,
        strand: b'.',
        meth,
        unmeth,
    })
}

impl Lines {
    fn push(&mut self, c: Cytosine, out: &mut impl FnMut(Line) -> Result<()>) -> Result<()> {
        if !self.merge {
            return out(single(&c)?);
        }
        let is_cpg = c.context() == 0;
        if let Some(p) = self.pending.take() {
            if is_cpg && c.reverse() && c.tid == p.tid && p.pos.checked_add(1) == Some(c.pos) {
                return out(cpg(
                    p.pos,
                    p.tid,
                    u64::from(p.meth) + u64::from(c.meth),
                    u64::from(p.unmeth) + u64::from(c.unmeth),
                )?);
            }
            out(cpg(p.pos, p.tid, u64::from(p.meth), u64::from(p.unmeth))?)?;
        }
        match (is_cpg, c.reverse()) {
            (true, false) => {
                self.pending = Some(c);
                Ok(())
            }
            (true, true) => out(cpg(
                c.pos.checked_sub(1).ok_or_else(|| {
                    Error::invalid(
                        "methylation output",
                        "reverse CpG has no preceding position",
                    )
                })?,
                c.tid,
                u64::from(c.meth),
                u64::from(c.unmeth),
            )?),
            _ => out(single(&c)?),
        }
    }

    fn finish(&mut self, out: &mut impl FnMut(Line) -> Result<()>) -> Result<()> {
        match self.pending.take() {
            Some(p) => out(cpg(p.pos, p.tid, u64::from(p.meth), u64::from(p.unmeth))?),
            None => Ok(()),
        }
    }
}

type Track = BigWigWriter<BufWriter<File>>;

/// Writes the table (`table`, relative to the output directory) and the
/// beta and coverage tracks (`tracks`) of `cytosines`.
pub fn write(
    cytosines: Spooled<Cytosine>,
    chroms: &[(String, u32)],
    merge: bool,
    table: Option<&str>,
    tracks: Option<(&str, &str)>,
    d: &mut Delivery<'_>,
) -> Result<Written> {
    let n = cytosines.len();
    let table_path = match table {
        Some(rel) => Some(d.path(rel)?),
        None => None,
    };
    let mut bed = match &table_path {
        Some(p) => Some(IndexedBed::create(d, p, "methylation", chroms)?),
        None => None,
    };
    let mut bigwigs: Option<[(Track, std::path::PathBuf); 2]> = match tracks {
        Some((beta, cov)) => {
            let mut open = |rel: &str| -> Result<(Track, std::path::PathBuf)> {
                let path = d.path(rel)?;
                let file = d.outputs.create(&path, "bigwig")?;
                let w = BigWigWriter::create_with(
                    BufWriter::with_capacity(1 << 20, file),
                    chroms
                        .iter()
                        .map(|(name, len)| Chrom {
                            name: name.clone(),
                            len: *len,
                        })
                        .collect(),
                    WriteOptions::default(),
                    Exec::Pool {
                        pool: d.pool,
                        limit: 4 * d.pool.threads(),
                    },
                )?;
                Ok((w, path))
            };
            Some([open(beta)?, open(cov)?])
        }
        None => None,
    };
    let mut lines = 0u64;
    // Track intervals must not overlap: a line starting inside the last
    // one (conflicting contexts at one position) is left out of them.
    let mut last: Option<(u32, u32)> = None;
    let mut text = Vec::with_capacity(128);
    let mut out = |l: Line| -> Result<()> {
        lines += 1;
        if let Some(b) = &mut bed {
            text.clear();
            text.extend_from_slice(chroms[l.tid as usize].0.as_bytes());
            text.push(b'\t');
            push_u64(&mut text, l.start as u64);
            text.push(b'\t');
            push_u64(&mut text, l.end as u64);
            text.push(b'\t');
            text.extend_from_slice(CONTEXTS[l.context].as_bytes());
            text.push(b'\t');
            text.push(l.strand);
            text.push(b'\t');
            push_u64(&mut text, l.meth);
            text.push(b'\t');
            push_u64(&mut text, l.unmeth);
            text.push(b'\n');
            b.push(l.tid, l.start, l.end, &text)?;
        }
        if let Some([(beta, _), (cov, _)]) = &mut bigwigs {
            let free = match last {
                Some((t, e)) => t != l.tid || l.start >= e,
                None => true,
            };
            let total = l.meth + l.unmeth;
            if free && total > 0 {
                beta.add(l.tid, l.start, l.end, (l.meth as f64 / total as f64) as f32)?;
                cov.add(l.tid, l.start, l.end, total as f32)?;
                last = Some((l.tid, l.end));
            }
        }
        Ok(())
    };
    let mut merger = Lines {
        merge,
        pending: None,
    };
    cytosines.drain(|c| merger.push(c, &mut out))?;
    merger.finish(&mut out)?;
    let index = match (bed, &table_path) {
        (Some(b), Some(p)) => Some(b.finish(p, d)?),
        _ => None,
    };
    if let Some(ws) = bigwigs {
        for (w, _) in ws {
            let (file, _) = w.finish()?;
            file.into_inner().map_err(|e| e.into_error())?;
        }
    }
    Ok(Written {
        cytosines: n,
        lines,
        index,
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn merged_cpg_counts_and_last_representable_end_do_not_wrap() {
        let forward = Cytosine {
            tid: 0,
            pos: u32::MAX - 2,
            code: 0,
            meth: u32::MAX,
            unmeth: u32::MAX,
        };
        let reverse = Cytosine {
            pos: u32::MAX - 1,
            code: 4,
            ..forward
        };
        let mut merger = Lines {
            merge: true,
            pending: None,
        };
        let mut rows = Vec::new();
        let mut output = |line| {
            rows.push(line);
            Ok(())
        };
        merger.push(forward, &mut output).unwrap();
        merger.push(reverse, &mut output).unwrap();
        merger.finish(&mut output).unwrap();
        assert_eq!(rows.len(), 1);
        assert_eq!((rows[0].start, rows[0].end), (u32::MAX - 2, u32::MAX));
        assert_eq!(
            (rows[0].meth, rows[0].unmeth),
            (u64::from(u32::MAX) * 2, u64::from(u32::MAX) * 2)
        );
    }

    #[test]
    fn unrepresentable_cytosine_and_cpg_ends_are_errors() {
        let c = Cytosine {
            tid: 0,
            pos: u32::MAX,
            code: 2,
            meth: 1,
            unmeth: 0,
        };
        assert!(single(&c).is_err());
        assert!(cpg(u32::MAX - 1, 0, 1, 0).is_err());
        let last = single(&Cytosine {
            pos: u32::MAX - 1,
            ..c
        })
        .unwrap();
        assert_eq!(last.end, u32::MAX);
        let mut merger = Lines {
            merge: true,
            pending: None,
        };
        assert!(
            merger
                .push(
                    Cytosine {
                        pos: 0,
                        code: 4,
                        ..c
                    },
                    &mut |_| Ok(())
                )
                .is_err()
        );
    }
}
