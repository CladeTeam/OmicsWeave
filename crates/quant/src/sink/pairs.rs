//! Contacts → 4DN pairs (`.pairs.gz`, BGZF): the header (`## pairs format
//! v1.0`, `#sorted: chr1-chr2-pos1-pos2`, `#shape: upper triangle`,
//! `#genome_assembly`, a `#chromsize` line per reference in reference
//! order, `#columns`), then a line per contact: `readID chr1 pos1 chr2
//! pos2 strand1 strand2 pair_type mapq1 mapq2` (and `barcode`), positions
//! 1-based, pair type `UU`.

use std::fs::File;
use std::io::BufWriter;
use std::path::Path;

use omicsweave_bgzf::Finished;
use omicsweave_core::Result;

use super::{GzText, push_u64};
use crate::accum::Delivery;
use crate::model::Contact;

pub struct PairsWriter<'c> {
    g: GzText,
    chroms: &'c [(String, u32)],
    line: Vec<u8>,
}

impl<'c> PairsWriter<'c> {
    /// Creates `path` and writes the header.
    pub fn create(
        d: &mut Delivery<'_>,
        path: &Path,
        chroms: &'c [(String, u32)],
        genome: &str,
        barcode: bool,
    ) -> Result<PairsWriter<'c>> {
        Self::create_mode(d, path, chroms, genome, barcode, true)
    }

    /// Creates pairs text without retaining a directory for an unused index.
    pub fn create_unindexed(
        d: &mut Delivery<'_>,
        path: &Path,
        chroms: &'c [(String, u32)],
        genome: &str,
        barcode: bool,
    ) -> Result<PairsWriter<'c>> {
        Self::create_mode(d, path, chroms, genome, barcode, false)
    }

    fn create_mode(
        d: &mut Delivery<'_>,
        path: &Path,
        chroms: &'c [(String, u32)],
        genome: &str,
        barcode: bool,
        indexed: bool,
    ) -> Result<PairsWriter<'c>> {
        let mut g = if indexed {
            GzText::create(d, path, "pairs")?
        } else {
            GzText::create_unindexed(d, path, "pairs")?
        };
        let mut head = String::from(
            "## pairs format v1.0\n#sorted: chr1-chr2-pos1-pos2\n#shape: upper triangle\n",
        );
        head.push_str(&format!("#genome_assembly: {genome}\n"));
        for (n, l) in chroms {
            head.push_str(&format!("#chromsize: {n} {l}\n"));
        }
        head.push_str("#columns: readID chr1 pos1 chr2 pos2 strand1 strand2 pair_type mapq1 mapq2");
        if barcode {
            head.push_str(" barcode");
        }
        head.push('\n');
        g.write(head.as_bytes())?;
        Ok(PairsWriter {
            g,
            chroms,
            line: Vec::with_capacity(256),
        })
    }

    pub fn push(&mut self, c: &Contact, barcode: Option<&[u8]>) -> Result<()> {
        let l = &mut self.line;
        l.clear();
        if c.name.is_empty() {
            l.push(b'.');
        } else {
            l.extend_from_slice(&c.name);
        }
        let strand = |s: bool| if s { b'-' } else { b'+' };
        for (chrom, pos) in [(c.c1, c.p1), (c.c2, c.p2)] {
            l.push(b'\t');
            l.extend_from_slice(self.chroms[chrom as usize].0.as_bytes());
            l.push(b'\t');
            push_u64(l, pos as u64);
        }
        l.extend_from_slice(&[b'\t', strand(c.s1), b'\t', strand(c.s2)]);
        l.extend_from_slice(b"\tUU\t");
        push_u64(l, c.mapq1 as u64);
        l.push(b'\t');
        push_u64(l, c.mapq2 as u64);
        if let Some(b) = barcode {
            l.push(b'\t');
            l.extend_from_slice(b);
        }
        l.push(b'\n');
        self.g.write(l)
    }

    pub fn finish(self) -> Result<Finished<BufWriter<File>>> {
        self.g.finish()
    }

    /// Finishes and flushes pairs text when no block directory is consumed.
    pub fn close(self) -> Result<()> {
        self.g.close()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use omicsweave_convert::stage::Outputs;
    use omicsweave_core::exec::Pool;
    use std::io::Read;

    #[test]
    fn unindexed_pairs_preserve_the_complete_bgzf_stream() {
        let dir = tempfile::tempdir().unwrap();
        let pool = Pool::new(2);
        let chroms = vec![
            ("chr1".to_string(), 1_000_000),
            ("chr2".to_string(), 1_000_000),
        ];
        let mut results = Vec::new();
        for indexed in [true, false] {
            let path = dir.path().join(format!("pairs-{indexed}.gz"));
            let mut outputs = Outputs::new();
            let mut delivery = Delivery {
                outputs: &mut outputs,
                dir: dir.path(),
                pool: &pool,
                matrix: crate::task::MatrixFormat::Mtx,
                tmp: dir.path(),
                buffer: 4 << 20,
            };
            let mut w = if indexed {
                PairsWriter::create(&mut delivery, &path, &chroms, "test", true)
            } else {
                PairsWriter::create_unindexed(&mut delivery, &path, &chroms, "test", true)
            }
            .unwrap();
            for i in 0..2000u32 {
                let c = Contact {
                    c1: 0,
                    p1: i + 1,
                    c2: 1,
                    p2: i + 2,
                    s1: i % 2 != 0,
                    s2: i % 3 != 0,
                    mapq1: 60,
                    mapq2: 50,
                    row: 0,
                    input: 0,
                    ordinal: i as u64,
                    name: if i % 7 == 0 {
                        vec![]
                    } else {
                        format!("read{i}").into_bytes()
                    },
                };
                w.push(&c, Some(b"cellA")).unwrap();
            }
            if indexed {
                w.finish().unwrap().out.into_inner().unwrap();
            } else {
                w.close().unwrap();
            }
            outputs.commit().unwrap();
            let bytes = std::fs::read(path).unwrap();
            let mut plain = String::new();
            flate2::read::MultiGzDecoder::new(&bytes[..])
                .read_to_string(&mut plain)
                .unwrap();
            assert_eq!(plain.lines().filter(|l| !l.starts_with('#')).count(), 2000);
            assert!(plain.contains(".\tchr1\t1\tchr2\t2\t+\t+\tUU\t60\t50\tcellA\n"));
            results.push(bytes);
        }
        assert_eq!(results[0], results[1]);
    }
}
