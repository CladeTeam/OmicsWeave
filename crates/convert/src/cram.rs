//! CRAM → BAM (+ BAI / CSI).

use std::fs::File;
use std::io::{BufWriter, Read};
use std::path::{Path, PathBuf};
use std::sync::Arc;

use omicsweave_bgzf::{BinningIndex, IndexFormat, Writer, try_csi_levels};
use omicsweave_core::exec::Ordered;
use omicsweave_core::{Error, FileSource, Locus, Result, SequenceRead};
use omicsweave_cram::{BamRecords, Container, CramHeader, CramOptions, Scanner, SharedReference};
use omicsweave_fasta::{Fai, FastaReader, Gzi, TWOBIT_MAGIC, TwoBitReader, fai_path, gzi_path};
use omicsweave_sam::bam_header_bytes;

use crate::budget::{Budget, MemoryPlan, WRITE_BUF, sizing};
use crate::conversion::{Context, Report};
use crate::format::Format;
use crate::pipeline::input::with_path;
use crate::pipeline::output::Outputs;
use crate::pipeline::records::block_sizing;
use crate::sam::BamIndex;

/// CRAM (3.0 or 3.1) to BAM, optionally indexed.
///
/// The uncompressed BAM stream is the one htslib writes when it reads the
/// CRAM and writes BAM (`samtools view -b` without `@PG`): its header
/// rebuilt from the CRAM's, every record as htslib decodes it (MD and NM
/// generated where htslib generates them, reads without a stored name
/// named after the file), cut into BGZF blocks where htslib cuts them, so
/// an index built here has the bins and chunks `samtools index` builds.
/// Slices decode in parallel; slices compressed against an external
/// reference need `reference` (FASTA with `.fai`, bgzip FASTA with `.fai`
/// and `.gzi`, or 2bit).
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct CramToBam {
    pub src: PathBuf,
    pub dst: PathBuf,
    /// libdeflate level, 0 to 12.
    pub level: i32,
    pub index: BamIndex,
    pub reference: Option<PathBuf>,
}

/// Bytes a slice holds in flight: its share of the container and the BAM
/// records it decodes to (htslib writes 10 000 reads a slice).
const SLICE_UNIT: u64 = 24 << 20;

fn slice_sizing(budget: &Budget) -> crate::budget::Sizing {
    sizing(budget.threads, SLICE_UNIT, budget.memory / 2)
}

/// A reference sequence file, read by name.
pub(crate) enum ReferenceFile {
    Fasta(FastaReader<FileSource>),
    TwoBit(TwoBitReader<FileSource>),
}

impl SequenceRead for ReferenceFile {
    fn chrom_len(&self, chrom: &str) -> Result<u32> {
        match self {
            ReferenceFile::Fasta(r) => r.chrom_len(chrom),
            ReferenceFile::TwoBit(r) => r.chrom_len(chrom),
        }
    }

    fn sequence(&mut self, locus: &Locus<'_>, out: &mut [u8]) -> Result<()> {
        match self {
            ReferenceFile::Fasta(r) => r.sequence(locus, out),
            ReferenceFile::TwoBit(r) => r.sequence(locus, out),
        }
    }
}

fn read_index(path: &Path) -> Result<Vec<u8>> {
    std::fs::read(path).map_err(|e| Error::Io {
        kind: e.kind(),
        detail: format!(
            "{}: {e}; build the reference's indexes first",
            path.display()
        ),
    })
}

/// Opens the reference at `path` by its content: 2bit, bgzip FASTA (with
/// `.fai` and `.gzi` beside it) or FASTA (with `.fai`).
pub(crate) fn open_reference(path: &Path) -> Result<ReferenceFile> {
    let mut head = [0u8; 4];
    let n = File::open(path)
        .and_then(|mut f| f.read(&mut head))
        .map_err(|e| with_path(e, path))?;
    let source = FileSource::open(path).map_err(|e| match e {
        Error::Io { kind, detail } => Error::Io {
            kind,
            detail: format!("{}: {detail}", path.display()),
        },
        e => e,
    })?;
    if n == 4 {
        let magic = u32::from_le_bytes(head);
        if magic == TWOBIT_MAGIC || magic.swap_bytes() == TWOBIT_MAGIC {
            return Ok(ReferenceFile::TwoBit(TwoBitReader::open(source, false)?));
        }
    }
    let fai = Fai::parse(&read_index(&fai_path(path))?)?;
    let gzi = if head.starts_with(&[0x1f, 0x8b]) {
        Some(Gzi::parse(&read_index(&gzi_path(path))?)?)
    } else {
        None
    };
    Ok(ReferenceFile::Fasta(FastaReader::new(source, fai, gzi)))
}

fn suffixed(path: &Path, suffix: &str) -> PathBuf {
    let mut s = path.as_os_str().to_owned();
    s.push(suffix);
    PathBuf::from(s)
}

impl CramToBam {
    pub fn new(src: impl Into<PathBuf>, dst: impl Into<PathBuf>) -> CramToBam {
        CramToBam {
            src: src.into(),
            dst: dst.into(),
            level: omicsweave_bgzf::DEFAULT_LEVEL,
            index: BamIndex::None,
            reference: None,
        }
    }

    pub(crate) fn plan(&self, budget: &Budget) -> Result<MemoryPlan> {
        omicsweave_bgzf::check_level(self.level)?;
        Ok(MemoryPlan {
            fixed: 0,
            pipeline: slice_sizing(budget).per_queue
                + 2 * block_sizing(budget).per_queue
                + WRITE_BUF as u64,
            buffer: 0,
        })
    }

    fn index_format(&self, max_len: u64) -> Result<Option<(IndexFormat, &'static str)>> {
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
        let source = FileSource::open(&self.src).map_err(|e| match e {
            Error::Io { kind, detail } => Error::Io {
                kind,
                detail: format!("{}: {detail}", self.src.display()),
            },
            e => e,
        })?;
        let file_len = omicsweave_core::ReadAt::len(&source);
        let header = Arc::new(CramHeader::read(&source)?);
        let reference = match &self.reference {
            Some(p) => Some(Arc::new(SharedReference::new(open_reference(p)?))),
            None => None,
        };
        let options = Arc::new(CramOptions {
            // htslib names unnamed reads after the file.
            name_prefix: self
                .src
                .file_name()
                .map(|n| n.to_string_lossy().as_bytes().to_vec())
                .unwrap_or_default(),
            ..CramOptions::default()
        });
        let sam = header.sam();
        let max_len = sam.references().iter().map(|r| r.len).max().unwrap_or(0);
        let index = self
            .index_format(max_len)?
            .map(|(f, suffix)| (f, sam.references().len(), suffixed(&self.dst, suffix)));

        let mut outputs = Outputs::new();
        let file = outputs.create(&self.dst, Format::Bam)?;
        let blocks = block_sizing(ctx.budget);
        let out = BufWriter::with_capacity(WRITE_BUF, file);
        let exec: Box<dyn Ordered<Vec<u8>>> = Box::new(ctx.pool.ordered(blocks.limit));
        let mut writer = if index.is_some() {
            Writer::new(out, self.level, exec)?
        } else {
            Writer::new_streaming(out, self.level, exec)?
        };
        writer.write(&bam_header_bytes(sam))?;
        writer.flush_block()?;
        let mut bins = index
            .as_ref()
            .map(|(format, n_refs, _)| BinningIndex::new(*format, *n_refs, writer.tell()));

        let mut scan = Scanner::new(source, &header, options.verify_checksums);
        let mut exec = ctx
            .pool
            .ordered::<BamRecords>(slice_sizing(ctx.budget).limit);
        let mut records = 0u64;
        let mut take = |r: BamRecords,
                        writer: &mut Writer<_>,
                        bins: &mut Option<BinningIndex>|
         -> Result<()> {
            let mut at = 0usize;
            for &(len, e) in &r.records {
                let rec = &r.bytes[at..at + len as usize];
                at += len as usize;
                writer.flush_try(rec.len())?;
                if let Some(ix) = bins.as_mut() {
                    ix.amend_last(writer.tell());
                }
                writer.write(rec)?;
                if let Some(ix) = bins.as_mut() {
                    ix.push(e.tid, e.beg, e.end, writer.tell(), e.mapped)?;
                }
            }
            records += r.records.len() as u64;
            Ok(())
        };
        while let Some(container) = scan.next_container()? {
            let container: Arc<Container> = Arc::new(container);
            for i in 0..container.slices() {
                while exec.full() {
                    take(exec.next(true)?.expect("pending"), &mut writer, &mut bins)?;
                }
                let (c, h, r, o) = (
                    Arc::clone(&container),
                    Arc::clone(&header),
                    reference.clone(),
                    Arc::clone(&options),
                );
                exec.submit(Box::new(move || {
                    let mut out = BamRecords::default();
                    c.decode_bam(i, &h, r.as_deref(), &o, &mut out)?;
                    Ok(out)
                }))?;
                while let Some(r) = exec.next(false)? {
                    take(r, &mut writer, &mut bins)?;
                }
            }
            ctx.progress("decode", scan.position(), file_len);
        }
        while let Some(r) = exec.next(true)? {
            take(r, &mut writer, &mut bins)?;
        }
        writer.flush_block()?;
        let end = writer.tell();
        if let (Some(ix), Some((kind, _, path))) = (bins, &index) {
            let done = writer.finish()?;
            done.out.into_inner().map_err(|e| e.into_error())?;
            let bytes = ix.finish(end, &done.blocks)?;
            let format = match kind {
                IndexFormat::Bai => Format::Bai,
                IndexFormat::Csi { .. } => Format::Csi,
            };
            outputs.write(path, format, &bytes)?;
        } else {
            writer
                .finish_streaming()?
                .into_inner()
                .map_err(|e| e.into_error())?;
        }
        Ok(Report {
            files: outputs.commit()?,
            records,
            bytes_read: scan.position(),
            spilled: 0,
            elapsed: Default::default(),
            warnings: Vec::new(),
        })
    }
}
