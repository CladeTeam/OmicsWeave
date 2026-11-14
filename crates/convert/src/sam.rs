//! SAM → BAM (+ BAI / CSI).

use std::path::PathBuf;
use std::sync::Arc;

use omicsweave_bgzf::{IndexFormat, try_csi_levels};
use omicsweave_core::{Error, Result};
use omicsweave_sam::{BamEncoder, BamScratch, SamHeader, bam_header_bytes};

use crate::budget::{Budget, MemoryPlan};
use crate::conversion::{Context, Report};
use crate::format::Format;
use crate::pipeline::input::open_text;
use crate::pipeline::records::{self, Entry, Framing, LineEncoder, Spec};

/// Which index to write beside a BAM.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum BamIndex {
    None,
    /// `.bai` when every reference fits (below 2^29 bases), else `.csi`.
    Auto,
    Bai,
    /// `.csi` with bins of `2^min_shift` bases.
    Csi {
        min_shift: u32,
    },
}

/// SAM (plain, gzip or BGZF) to BAM, optionally indexed.
///
/// The uncompressed BAM stream is the one htslib writes for the same SAM
/// (`samtools view -b`), and records are cut into BGZF blocks where htslib
/// cuts them, so an index built here has the bins and chunks
/// `samtools index` builds for the file. Records keep their input order;
/// an index requires coordinate-sorted input. References missing from the
/// header are errors (htslib makes such reads unmapped).
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct SamToBam {
    pub src: PathBuf,
    pub dst: PathBuf,
    /// libdeflate level, 0 to 12.
    pub level: i32,
    pub index: BamIndex,
}

struct Bam(BamEncoder);

impl LineEncoder for Bam {
    type Scratch = BamScratch;

    fn scratch(&self) -> BamScratch {
        BamScratch::new()
    }

    fn encode(
        &self,
        scratch: &mut BamScratch,
        line: &[u8],
        line_no: u64,
        out: &mut Vec<u8>,
    ) -> Result<Option<Entry>> {
        let s = self.0.encode_line(scratch, line, line_no, out)?;
        Ok(Some(Entry {
            tid: s.tid,
            beg: s.beg,
            end: s.end,
            mapped: s.mapped,
        }))
    }
}

fn suffixed(path: &std::path::Path, suffix: &str) -> PathBuf {
    let mut s = path.as_os_str().to_owned();
    s.push(suffix);
    PathBuf::from(s)
}

impl SamToBam {
    pub fn new(src: impl Into<PathBuf>, dst: impl Into<PathBuf>) -> SamToBam {
        SamToBam {
            src: src.into(),
            dst: dst.into(),
            level: omicsweave_bgzf::DEFAULT_LEVEL,
            index: BamIndex::None,
        }
    }

    pub(crate) fn plan(&self, budget: &Budget) -> Result<MemoryPlan> {
        omicsweave_bgzf::check_level(self.level)?;
        Ok(records::plan(budget, 0))
    }

    fn index_format(&self, header: &SamHeader) -> Result<Option<(IndexFormat, &'static str)>> {
        let max_len = header.references().iter().map(|r| r.len).max().unwrap_or(0);
        let fits_bai = max_len <= 1 << 29;
        Ok(match self.index {
            BamIndex::None => None,
            BamIndex::Bai if !fits_bai => {
                return Err(Error::invalid(
                    "bam index",
                    format!("a reference of {max_len} bases does not fit a BAI index; use CSI"),
                ));
            }
            BamIndex::Bai => Some((IndexFormat::Bai, ".bai")),
            BamIndex::Auto if fits_bai => Some((IndexFormat::Bai, ".bai")),
            BamIndex::Auto => {
                let (min_shift, levels) = try_csi_levels(max_len, 14)?;
                Some((IndexFormat::Csi { min_shift, levels }, ".csi"))
            }
            BamIndex::Csi { min_shift } => {
                let (min_shift, levels) = try_csi_levels(max_len, min_shift)?;
                Some((IndexFormat::Csi { min_shift, levels }, ".csi"))
            }
        })
    }

    pub(crate) fn run(&self, ctx: &Context<'_>) -> Result<Report> {
        let s = records::batch_sizing(ctx.budget);
        let mut input = open_text(&self.src, &ctx.pool, s.limit)?;
        let header = SamHeader::read(&mut input.text)?;
        let index = self
            .index_format(&header)?
            .map(|(f, suffix)| (f, header.references().len(), suffixed(&self.dst, suffix)));
        let spec = Spec {
            dst: &self.dst,
            format: Format::Bam,
            level: self.level,
            header: bam_header_bytes(&header),
            framing: Framing::Record,
            encoder: Arc::new(Bam(BamEncoder::new(&header))),
            index,
            text: input.text,
            first_line: header.lines() + 1,
            read: input.read,
            file_len: input.file_len,
        };
        records::run(ctx, spec)
    }
}
