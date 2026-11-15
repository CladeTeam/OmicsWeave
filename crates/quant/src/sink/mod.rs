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
        let file = d.outputs.create(path, kind)?;
        let limit = 4 * d.pool.threads();
        Ok(GzText {
            w: Writer::new(
                BufWriter::with_capacity(1 << 20, file),
                TEXT_LEVEL,
                Box::new(d.pool.ordered(limit)),
            )?,
        })
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

    pub fn finish(self) -> Result<Finished<BufWriter<File>>> {
        let done = self.w.finish()?;
        Ok(done)
    }

    /// Finishes and flushes the file.
    pub fn close(self) -> Result<()> {
        let done = self.finish()?;
        done.out.into_inner().map_err(|e| e.into_error())?;
        Ok(())
    }
}

/// `text` as a BGZF file at `path`.
pub fn write_gz(d: &mut Delivery<'_>, path: &Path, kind: &'static str, text: &[u8]) -> Result<()> {
    let mut g = GzText::create(d, path, kind)?;
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
