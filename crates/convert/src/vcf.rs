//! VCF → BCF (+ CSI).

use std::collections::HashSet;
use std::io::BufRead;
use std::path::PathBuf;
use std::sync::Arc;
use std::sync::atomic::Ordering;

use omicsweave_bgzf::{IndexFormat, try_csi_levels};
use omicsweave_core::{Error, Result};
use omicsweave_vcf::{BcfEncoder, BcfScratch, Header, bcf_header_bytes, read_vcf_header};

use crate::budget::text_batch;
use crate::budget::{Budget, MemoryPlan};
use crate::conversion::{Context, Report};
use crate::format::Format;
use crate::pipeline::input::{LineBatches, lines, open_text};
use crate::pipeline::records::{self, Entry, Framing, LineEncoder, Spec};

/// VCF (plain, gzip or BGZF) to BCF, optionally with its CSI index.
///
/// The uncompressed BCF stream is the one htslib writes for the same VCF
/// (`bcftools view -Ob`), cut into BGZF blocks as htslib cuts it, so the
/// index has the bins and chunks `bcftools index` builds for the file.
/// Contigs, FILTERs, INFO and FORMAT keys must be declared in the header:
/// BCF stores the header first, so undeclared names cannot be added later.
/// With [`VcfToBcf::add_missing_contigs`], a first pass over the CHROM
/// column declares the contigs the header lacks (many VCFs, ClinVar's
/// among them, have no `##contig` lines).
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct VcfToBcf {
    pub src: PathBuf,
    pub dst: PathBuf,
    /// libdeflate level, 0 to 12.
    pub level: i32,
    /// Also write `<dst>.csi` (records must then be sorted).
    pub index: bool,
    /// Declare contigs missing from the header, in order of first use, as
    /// `##contig=<ID=...>` lines before `#CHROM` (costs a pass over the
    /// input). Off, an undeclared contig is an error.
    pub add_missing_contigs: bool,
}

struct Bcf(BcfEncoder);

impl LineEncoder for Bcf {
    type Scratch = BcfScratch;

    fn scratch(&self) -> BcfScratch {
        BcfScratch::default()
    }

    fn encode(
        &self,
        scratch: &mut BcfScratch,
        line: &[u8],
        line_no: u64,
        out: &mut Vec<u8>,
    ) -> Result<Option<Entry>> {
        let s = self.0.encode_line(scratch, line, line_no, out)?;
        Ok(Some(Entry {
            tid: s.rid as i32,
            beg: s.beg,
            end: s.end,
            mapped: true,
        }))
    }
}

impl VcfToBcf {
    pub fn new(src: impl Into<PathBuf>, dst: impl Into<PathBuf>) -> VcfToBcf {
        VcfToBcf {
            src: src.into(),
            dst: dst.into(),
            level: omicsweave_bgzf::DEFAULT_LEVEL,
            index: false,
            add_missing_contigs: false,
        }
    }

    pub(crate) fn plan(&self, budget: &Budget) -> Result<MemoryPlan> {
        omicsweave_bgzf::check_level(self.level)?;
        Ok(records::plan(budget, 0))
    }

    pub(crate) fn run(&self, ctx: &Context<'_>) -> Result<Report> {
        let s = records::batch_sizing(ctx.budget);
        let mut input = open_text(&self.src, &ctx.pool, s.limit)?;
        let (mut header, consumed) = read_vcf_header(&mut input.text)?;
        if self.add_missing_contigs {
            header = self.with_missing_contigs(ctx, header)?;
        }
        let index = if self.index {
            // As bcftools: bins of 2^14 bases, levels for the longest contig
            // (2^31 - 1 when no contig has a length).
            let max_len = header
                .contigs
                .iter()
                .filter_map(|c| c.len)
                .max()
                .unwrap_or((1 << 31) - 1);
            let (min_shift, levels) = try_csi_levels(max_len, 14)?;
            let mut path = self.dst.as_os_str().to_owned();
            path.push(".csi");
            Some((
                IndexFormat::Csi { min_shift, levels },
                header.contigs.len(),
                PathBuf::from(path),
            ))
        } else {
            None
        };
        let spec = Spec {
            dst: &self.dst,
            format: Format::Bcf,
            level: self.level,
            header: bcf_header_bytes(&header)?,
            framing: Framing::Stream,
            encoder: Arc::new(Bcf(BcfEncoder::new(&header))),
            index,
            text: input.text,
            first_line: consumed + 1,
            read: input.read,
            file_len: input.file_len,
        };
        records::run(ctx, spec)
    }

    /// `header` with `##contig` lines for the contigs the data uses and it
    /// does not declare, inserted before `#CHROM` in order of first use.
    fn with_missing_contigs(&self, ctx: &Context<'_>, header: Header) -> Result<Header> {
        let declared: HashSet<&str> = header.contigs.iter().map(|c| c.name.as_str()).collect();
        let s = records::batch_sizing(ctx.budget);
        let mut input = open_text(&self.src, &ctx.pool, s.limit)?;
        let mut head = Vec::new();
        let chrom_line = loop {
            let mut line = Vec::new();
            if input.text.read_until(b'\n', &mut line)? == 0 {
                return Err(Error::invalid("vcf", "no #CHROM header line"));
            }
            if line.starts_with(b"#CHROM") {
                break line;
            }
            head.extend_from_slice(&line);
        };
        let (read, file_len) = (Arc::clone(&input.read), input.file_len);
        let mut batches = LineBatches::new(input.text, text_batch(ctx.budget), 1);
        let mut seen: HashSet<Vec<u8>> = HashSet::new();
        let mut missing: Vec<Vec<u8>> = Vec::new();
        while let Some(batch) = batches.next_batch()? {
            for (_, line) in lines(&batch) {
                let name = line.split(|&c| c == b'\t').next().unwrap_or(b"");
                let name = name.strip_suffix(b"\r").unwrap_or(name);
                if name.is_empty() || seen.contains(name) {
                    continue;
                }
                seen.insert(name.to_vec());
                let known = std::str::from_utf8(name).is_ok_and(|n| declared.contains(n));
                if !known {
                    missing.push(name.to_vec());
                }
            }
            ctx.progress("scan", read.load(Ordering::Relaxed), file_len);
        }
        if missing.is_empty() {
            return Ok(header);
        }
        for name in &missing {
            head.extend_from_slice(b"##contig=<ID=");
            head.extend_from_slice(name);
            head.extend_from_slice(b">\n");
        }
        head.extend_from_slice(&chrom_line);
        Ok(read_vcf_header(&mut &head[..])?.0)
    }
}
