//! Sinks: products written as files.
//!
//! | Product | Files |
//! |---|---|
//! | [`Track`](crate::model::Track) | bigWig ([`bigwig`]), with an optional `.sidecar` |
//! | [`CountMatrix`](crate::model::CountMatrix) | a 10x directory: `barcodes.tsv.gz`, `features.tsv.gz`, `matrix.cpk` (CellPack), `matrix.mtx.gz` or `matrix.csv.gz`; optionally a featureCounts-style `counts.tsv` ([`matrix`]) |
//! | [`FragmentTable`](crate::model::FragmentTable) | 10x `fragments.tsv.gz` with its tabix `.tbi` ([`fragments`]) |
//! | [`BarcodeMetrics`](crate::model::BarcodeMetrics) | CSV ([`table`]) |
//! | [`JunctionTable`](crate::model::JunctionTable) | `junctions.tab`, STAR's `SJ.out.tab` columns ([`junctions`]) |
//! | [`Histogram`](crate::model::Histogram) | `conversions.tsv.gz` ([`table`]) |
//! | [`MismatchTable`](crate::model::MismatchTable) | `mismatches.tsv` ([`table`]) |
//! | [`Cytosine`](crate::model::Cytosine)s | `cytosines.tsv.gz` with its `.tbi`, `beta.bw`, `coverage.bw` ([`methylation`]) |
//! | [`Mbias`](crate::model::Mbias) | `mbias.tsv` ([`table`]) |
//! | [`ContactTable`](crate::model::ContactTable) | 4DN `contacts.pairs.gz` ([`pairs`]) and Juicer `contacts.hic` ([`hic`]) |
//! | [`VariantTable`](crate::model::VariantTable) | `variants.tsv.gz` ([`table`]) |
//! | the run | `summary.json` ([`json`]) |
//!
//! Compressed text is BGZF (gzip-compatible, compressed on the worker
//! pool). Every file is written under a temporary name and renamed into
//! place when the whole run succeeded.

pub mod bigwig;
pub mod fragments;
pub mod hic;
pub mod json;
pub mod junctions;
pub mod matrix;
pub mod methylation;
pub mod pairs;
pub mod table;

use std::fs::File;
use std::io::BufWriter;
use std::path::Path;

use omicsweave_bgzf::{Finished, Logical, Writer};
use omicsweave_core::Result;

use crate::accum::Delivery;

/// BGZF compression level of text outputs.
pub const TEXT_LEVEL: i32 = 6;

/// BGZF-compressed text written to an output.
pub struct GzText {
    w: Writer<BufWriter<File>>,
}

impl GzText {
    /// Creates `path` (kind `kind`) among the delivery's outputs.
    pub fn create(d: &mut Delivery<'_>, path: &Path, kind: &'static str) -> Result<GzText> {
        Self::create_mode(d, path, kind, true)
    }

    /// Writes text without retaining a directory of all physical BGZF blocks.
    /// Use [`GzText::close`] when no tabix/GZI index is needed.
    pub fn create_unindexed(
        d: &mut Delivery<'_>,
        path: &Path,
        kind: &'static str,
    ) -> Result<GzText> {
        Self::create_mode(d, path, kind, false)
    }

    fn create_mode(
        d: &mut Delivery<'_>,
        path: &Path,
        kind: &'static str,
        indexed: bool,
    ) -> Result<GzText> {
        let file = d.outputs.create(path, kind)?;
        let limit = 4 * d.pool.threads();
        let out = BufWriter::with_capacity(1 << 20, file);
        let w = if indexed {
            Writer::new(out, TEXT_LEVEL, Box::new(d.pool.ordered(limit)))?
        } else {
            Writer::new_streaming(out, TEXT_LEVEL, Box::new(d.pool.ordered(limit)))?
        };
        Ok(GzText { w })
    }

    #[inline]
    pub fn write(&mut self, text: &[u8]) -> Result<()> {
        self.w.write(text)
    }

    /// Position after what has been written.
    pub fn tell(&self) -> Logical {
        self.w.tell()
    }

    pub fn flush_block(&mut self) -> Result<()> {
        self.w.flush_block()
    }

    /// Returns the block directory created by [`GzText::create`].
    /// Unindexed text must be completed with [`GzText::close`] instead.
    pub fn finish(self) -> Result<Finished<BufWriter<File>>> {
        let done = self.w.finish()?;
        Ok(done)
    }

    /// Finishes and flushes the file.
    pub fn close(self) -> Result<()> {
        let out = self.w.finish_streaming()?;
        out.into_inner().map_err(|e| e.into_error())?;
        Ok(())
    }
}

/// `text` as a BGZF file at `path`.
pub fn write_gz(d: &mut Delivery<'_>, path: &Path, kind: &'static str, text: &[u8]) -> Result<()> {
    let mut g = GzText::create_unindexed(d, path, kind)?;
    g.write(text)?;
    g.close()
}

/// Appends the decimal digits of `v`.
#[inline]
pub(crate) fn push_u64(out: &mut Vec<u8>, mut v: u64) {
    let mut buf = [0u8; 20];
    let mut i = buf.len();
    loop {
        i -= 1;
        buf[i] = b'0' + (v % 10) as u8;
        v /= 10;
        if v == 0 {
            break;
        }
    }
    out.extend_from_slice(&buf[i..]);
}

#[cfg(test)]
mod tests {
    use super::*;
    use omicsweave_convert::stage::Outputs;
    use omicsweave_core::exec::Pool;
    use std::io::Read;

    #[test]
    fn unindexed_text_preserves_members_plaintext_and_logical_positions() {
        let dir = tempfile::tempdir().unwrap();
        let pool = Pool::new(2);
        for len in [0usize, 400_003] {
            let text: Vec<u8> = (0..len).map(|i| ((i * 37 + i / 31) % 251) as u8).collect();
            let mut results = Vec::new();
            let mut positions = Vec::new();
            for indexed in [true, false] {
                let path = dir.path().join(format!("text-{len}-{indexed}.gz"));
                let mut outputs = Outputs::new();
                let mut delivery = Delivery {
                    outputs: &mut outputs,
                    dir: dir.path(),
                    pool: &pool,
                    matrix: crate::task::MatrixFormat::Mtx,
                    tmp: dir.path(),
                    buffer: 4 << 20,
                };
                let mut g = if indexed {
                    GzText::create(&mut delivery, &path, "text")
                } else {
                    GzText::create_unindexed(&mut delivery, &path, "text")
                }
                .unwrap();
                let mut tells = vec![(g.tell(), 0u64)];
                let mut written = 0;
                for (i, chunk) in text.chunks(5000).enumerate() {
                    g.write(chunk).unwrap();
                    written += chunk.len();
                    if i % 3 == 0 {
                        g.flush_block().unwrap();
                    }
                    if i % 17 == 0 {
                        g.w.sync().unwrap();
                    }
                    tells.push((g.tell(), written as u64));
                }
                if indexed {
                    let done = g.finish().unwrap();
                    assert_eq!(done.blocks.text_len, len as u64);
                    for &(logical, text_at) in &tells {
                        let physical = done.blocks.resolve(logical);
                        let start = done
                            .blocks
                            .starts
                            .iter()
                            .find(|&&(c, _)| c == physical >> 16)
                            .map_or(done.blocks.text_len, |&(_, t)| t);
                        assert_eq!(start + (physical & 0xffff), text_at);
                    }
                    done.out.into_inner().unwrap();
                } else {
                    g.close().unwrap();
                }
                outputs.commit().unwrap();
                let data = std::fs::read(path).unwrap();
                assert!(data.ends_with(&omicsweave_bgzf::EOF_BLOCK));
                let mut plain = Vec::new();
                flate2::read::MultiGzDecoder::new(&data[..])
                    .read_to_end(&mut plain)
                    .unwrap();
                assert_eq!(plain, text);
                results.push(data);
                positions.push(tells);
            }
            assert_eq!(results[0], results[1]);
            assert_eq!(positions[0], positions[1]);
        }
    }
}
