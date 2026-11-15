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
        let mut g = GzText::create(d, path, "pairs")?;
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
}
