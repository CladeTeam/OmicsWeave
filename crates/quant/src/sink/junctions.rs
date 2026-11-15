//! [`JunctionTable`] → `junctions.tab`, the nine columns of STAR's
//! `SJ.out.tab`: chromosome, first and last intron base (1-based), strand
//! (0 undefined, 1 `+`, 2 `-`), motif, annotated (0 / 1), unique reads,
//! multi-mapping reads, largest overhang.

use std::io::{BufWriter, Write};
use std::path::Path;

use omicsweave_core::Result;

use super::push_u64;
use crate::accum::Delivery;
use crate::model::JunctionTable;

pub fn write(table: &JunctionTable, path: &Path, d: &mut Delivery<'_>) -> Result<()> {
    let mut out = BufWriter::with_capacity(1 << 20, d.outputs.create(path, "junctions")?);
    let mut line = Vec::with_capacity(128);
    for j in &table.junctions {
        line.clear();
        line.extend_from_slice(table.chroms[j.tid as usize].0.as_bytes());
        for v in [
            j.start as u64 + 1,
            j.end as u64,
            j.strand as u64,
            j.motif as u64,
            j.annotated as u64,
            j.unique as u64,
            j.multi as u64,
            j.max_overhang as u64,
        ] {
            line.push(b'\t');
            push_u64(&mut line, v);
        }
        line.push(b'\n');
        out.write_all(&line)?;
    }
    out.flush()?;
    Ok(())
}
