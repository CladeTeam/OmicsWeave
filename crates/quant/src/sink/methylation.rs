//! Per-cytosine methylation → `cytosines.tsv.gz` (BGZF, tabix-indexed as
//! BED: chrom, start, end, context, strand, methylated, unmethylated;
//! 0-based half-open) and the `beta.bw` / `coverage.bw` tracks of the same
//! lines, written in one pass over the cytosines.

use std::fs::File;
use std::io::BufWriter;

use omicsweave_bbi::{BigWigWriter, Chrom, WriteOptions};
use omicsweave_core::Result;
use omicsweave_core::exec::Exec;

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
}

/// One line: `[start, end)` with its strand (`+`, `-` or `.`).
struct Line {
    tid: u32,
    start: u32,
    end: u32,
    context: usize,
    strand: u8,
    meth: u32,
    unmeth: u32,
}

/// Merges the two cytosines of each CpG into one line when `merge`.
struct Lines {
    merge: bool,
    /// A forward CpG waiting for its reverse partner.
    pending: Option<Cytosine>,
}

fn single(c: &Cytosine) -> Line {
    Line {
        tid: c.tid,
        start: c.pos,
        end: c.pos + 1,
        context: c.context(),
        strand: if c.reverse() { b'-' } else { b'+' },
        meth: c.meth,
        unmeth: c.unmeth,
    }
}

fn cpg(start: u32, tid: u32, meth: u32, unmeth: u32) -> Line {
    Line {
        tid,
        start,
        end: start + 2,
        context: 0,
        strand: b'.',
        meth,
        unmeth,
    }
}

impl Lines {
    fn push(&mut self, c: Cytosine, out: &mut impl FnMut(Line) -> Result<()>) -> Result<()> {
        if !self.merge {
            return out(single(&c));
        }
        let is_cpg = c.context() == 0;
        if let Some(p) = self.pending.take() {
            if is_cpg && c.reverse() && c.tid == p.tid && c.pos == p.pos + 1 {
                return out(cpg(p.pos, p.tid, p.meth + c.meth, p.unmeth + c.unmeth));
            }
            out(cpg(p.pos, p.tid, p.meth, p.unmeth))?;
        }
        match (is_cpg, c.reverse()) {
            (true, false) => {
                self.pending = Some(c);
                Ok(())
            }
            (true, true) => out(cpg(c.pos.saturating_sub(1), c.tid, c.meth, c.unmeth)),
            _ => out(single(&c)),
        }
    }

    fn finish(&mut self, out: &mut impl FnMut(Line) -> Result<()>) -> Result<()> {
        match self.pending.take() {
            Some(p) => out(cpg(p.pos, p.tid, p.meth, p.unmeth)),
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
            push_u64(&mut text, l.meth as u64);
            text.push(b'\t');
            push_u64(&mut text, l.unmeth as u64);
            text.push(b'\n');
            b.push(l.tid, l.start, l.end, &text)?;
        }
        if let Some([(beta, _), (cov, _)]) = &mut bigwigs {
            let free = match last {
                Some((t, e)) => t != l.tid || l.start >= e,
                None => true,
            };
            let total = l.meth as u64 + l.unmeth as u64;
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
    if let (Some(b), Some(p)) = (bed, &table_path) {
        b.finish(p, d)?;
    }
    if let Some(ws) = bigwigs {
        for (w, _) in ws {
            let (file, _) = w.finish()?;
            file.into_inner().map_err(|e| e.into_error())?;
        }
    }
    Ok(Written {
        cytosines: n,
        lines,
    })
}
