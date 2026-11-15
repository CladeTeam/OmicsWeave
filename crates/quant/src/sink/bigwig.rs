//! [`Track`] → bigWig, with zoom levels and an optional sidecar.

use std::io::BufWriter;
use std::path::Path;

use omicsweave_bbi::{BigWigWriter, Chrom, WriteOptions};
use omicsweave_core::Result;
use omicsweave_core::exec::Exec;

use crate::accum::Delivery;
use crate::model::Track;

/// Writes `track` (every value times its scale) to `path`; returns the
/// number of runs written.
pub fn write(track: Track, path: &Path, sidecar: bool, d: &mut Delivery<'_>) -> Result<u64> {
    let file = d.outputs.create(path, "bigwig")?;
    let chroms = track
        .chroms
        .iter()
        .map(|(name, len)| Chrom {
            name: name.clone(),
            len: *len,
        })
        .collect();
    let mut w = BigWigWriter::create_with(
        BufWriter::with_capacity(1 << 20, file),
        chroms,
        WriteOptions::default(),
        Exec::Pool {
            pool: d.pool,
            limit: 4 * d.pool.threads(),
        },
    )?;
    let scale = track.scale;
    let mut n = 0u64;
    track.runs.drain(|r| {
        n += 1;
        let v = if scale == 1.0 {
            r.value as f32
        } else {
            (r.value as f64 * scale) as f32
        };
        w.add(r.tid, r.start, r.end, v)
    })?;
    let (out, car) = w.finish()?;
    out.into_inner().map_err(|e| e.into_error())?;
    if sidecar {
        let source = d.outputs.source(path)?;
        let file = d
            .outputs
            .create(&omicsweave_core::sidecar::sidecar_path(path), "sidecar")?;
        car.write_to(&source, file)?;
    }
    Ok(n)
}
