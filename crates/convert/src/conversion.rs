//! A conversion from one file to another: choosing it, checking it against
//! a budget, and running it.

use std::path::{Path, PathBuf};
use std::time::{Duration, Instant};

use omicsweave_core::exec::Pool;
use omicsweave_core::{Error, Result};

use crate::bed::{BedToBigBed, ChromSizes};
use crate::budget::{self, Budget, MemoryPlan};
use crate::cool::CoolToHic;
use crate::fasta::{FastaToBgzipFasta, FastaToTwoBit};
use crate::format::Format;
use crate::sam::{BamIndex, SamToBam};
use crate::vcf::VcfToBcf;

/// A conversion of one source file into one destination file (and its
/// companions), with its options.
#[derive(Debug, Clone, PartialEq)]
pub enum Conversion {
    FastaToBgzipFasta(FastaToBgzipFasta),
    FastaToTwoBit(FastaToTwoBit),
    VcfToBcf(VcfToBcf),
    SamToBam(SamToBam),
    BedToBigBed(BedToBigBed),
    CoolToHic(CoolToHic),
}

/// Which index to write beside a sorted BCF or BAM.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Index {
    /// BCF: `.csi`. BAM: `.bai` when every reference fits, else `.csi`.
    Auto,
    Bai,
    Csi,
}

impl Index {
    /// The index called `name`: `auto`, `bai` or `csi`.
    pub fn from_name(name: &str) -> Result<Index> {
        match name {
            "auto" => Ok(Index::Auto),
            "bai" => Ok(Index::Bai),
            "csi" => Ok(Index::Csi),
            _ => Err(Error::invalid(
                "index",
                format!("unknown index {name:?}: use \"auto\", \"bai\" or \"csi\""),
            )),
        }
    }
}

/// The options of a conversion. Each applies to the conversions that have
/// it; [`Conversion::configure`] rejects one given to any other.
#[derive(Debug, Clone, PartialEq)]
pub struct Options {
    /// Compression level (bgzip, BCF, BAM: 0-12; bigBed, hic: zlib 1-9).
    pub level: Option<i32>,
    /// BCF, BAM: also write an index (records must be sorted).
    pub index: Option<Index>,
    /// BED to bigBed: chromosome lengths (`chrom.sizes`, `.fai` or
    /// `.2bit`) instead of each chromosome's largest interval end.
    pub chrom_sizes: Option<PathBuf>,
    /// `.mcool` to hic: the resolutions to convert instead of all.
    pub resolutions: Option<Vec<u32>>,
    /// bigBed, hic: also write the `.sidecar`.
    pub sidecar: bool,
    /// VCF to BCF: declare contigs the header lacks (one more pass).
    pub add_contigs: bool,
}

impl Default for Options {
    fn default() -> Options {
        Options {
            level: None,
            index: None,
            chrom_sizes: None,
            resolutions: None,
            sidecar: true,
            add_contigs: false,
        }
    }
}

fn not_for(option: &str, conversion: &Conversion) -> Error {
    Error::invalid(
        "convert",
        format!("option {option} does not apply to {}", conversion.name()),
    )
}

impl Conversion {
    /// The conversion of `src`, holding `from`, into `dst`, to hold `to`,
    /// with default options. A source is read in any of its storages, so
    /// `from` may be any of them.
    pub fn between(
        from: Format,
        to: Format,
        src: impl Into<PathBuf>,
        dst: impl Into<PathBuf>,
    ) -> Result<Conversion> {
        use Format::*;
        let (src, dst) = (src.into(), dst.into());
        Ok(match (from, to) {
            (Fasta | GzipFasta | BgzipFasta, BgzipFasta) => {
                Conversion::FastaToBgzipFasta(FastaToBgzipFasta::new(src, dst))
            }
            (Fasta | GzipFasta | BgzipFasta, TwoBit) => {
                Conversion::FastaToTwoBit(FastaToTwoBit::new(src, dst))
            }
            (Vcf | GzipVcf, Bcf) => Conversion::VcfToBcf(VcfToBcf::new(src, dst)),
            (Sam | GzipSam, Bam) => Conversion::SamToBam(SamToBam::new(src, dst)),
            (Bed | GzipBed, BigBed) => {
                Conversion::BedToBigBed(BedToBigBed::new(src, dst, ChromSizes::FromData))
            }
            (Cool | Mcool, Hic) => Conversion::CoolToHic(CoolToHic::new(src, dst)),
            _ => {
                return Err(Error::invalid(
                    "convert",
                    format!("no conversion from {from} to {to}"),
                ));
            }
        })
    }

    /// The conversion of `src` into `dst`, both formats told by their
    /// names ([`Format::from_extension`]).
    pub fn infer(src: impl Into<PathBuf>, dst: impl Into<PathBuf>) -> Result<Conversion> {
        let (src, dst) = (src.into(), dst.into());
        let (from, to) = (Format::from_extension(&src)?, Format::from_extension(&dst)?);
        Conversion::between(from, to, src, dst)
    }

    /// Applies `options`; an option this conversion does not have is an error.
    pub fn configure(&mut self, o: &Options) -> Result<()> {
        let takes = |option: bool, applies: bool, name: &str, c: &Conversion| {
            if option && !applies {
                Err(not_for(name, c))
            } else {
                Ok(())
            }
        };
        takes(
            o.level.is_some(),
            !matches!(self, Conversion::FastaToTwoBit(_)),
            "level",
            self,
        )?;
        takes(
            o.index.is_some(),
            matches!(self, Conversion::VcfToBcf(_) | Conversion::SamToBam(_)),
            "index",
            self,
        )?;
        takes(
            o.chrom_sizes.is_some(),
            matches!(self, Conversion::BedToBigBed(_)),
            "chrom_sizes",
            self,
        )?;
        takes(
            o.resolutions.is_some(),
            matches!(self, Conversion::CoolToHic(_)),
            "resolutions",
            self,
        )?;
        takes(
            o.add_contigs,
            matches!(self, Conversion::VcfToBcf(_)),
            "add_contigs",
            self,
        )?;
        match self {
            Conversion::FastaToBgzipFasta(c) => c.level = o.level.unwrap_or(c.level),
            Conversion::FastaToTwoBit(_) => {}
            Conversion::VcfToBcf(c) => {
                c.level = o.level.unwrap_or(c.level);
                c.add_missing_contigs = o.add_contigs;
                c.index = match o.index {
                    None => false,
                    Some(Index::Auto | Index::Csi) => true,
                    Some(Index::Bai) => {
                        return Err(Error::invalid(
                            "convert",
                            "a BCF is indexed by .csi, not .bai",
                        ));
                    }
                };
            }
            Conversion::SamToBam(c) => {
                c.level = o.level.unwrap_or(c.level);
                c.index = match o.index {
                    None => BamIndex::None,
                    Some(Index::Auto) => BamIndex::Auto,
                    Some(Index::Bai) => BamIndex::Bai,
                    Some(Index::Csi) => BamIndex::Csi { min_shift: 14 },
                };
            }
            Conversion::BedToBigBed(c) => {
                if let Some(p) = &o.chrom_sizes {
                    c.chrom_sizes = ChromSizes::File(p.clone());
                }
                c.options.level = o.level.unwrap_or(c.options.level);
                c.sidecar = o.sidecar;
            }
            Conversion::CoolToHic(c) => {
                if let Some(l) = o.level {
                    c.level = u32::try_from(l)
                        .map_err(|_| Error::invalid("convert", "hic level must be 1 to 9"))?;
                }
                c.resolutions = o.resolutions.clone();
                c.sidecar = o.sidecar;
            }
        }
        Ok(())
    }

    /// Every (source, destination) pair of formats there is a conversion for.
    pub fn routes() -> impl Iterator<Item = (Format, Format)> {
        Format::ALL.into_iter().flat_map(|from| {
            Format::ALL
                .into_iter()
                .filter(move |&to| Conversion::between(from, to, "", "").is_ok())
                .map(move |to| (from, to))
        })
    }

    /// What the conversion does, for messages: `FASTA to 2bit`, ...
    pub fn name(&self) -> &'static str {
        match self {
            Conversion::FastaToBgzipFasta(_) => "FASTA to bgzip FASTA",
            Conversion::FastaToTwoBit(_) => "FASTA to 2bit",
            Conversion::VcfToBcf(_) => "VCF to BCF",
            Conversion::SamToBam(_) => "SAM to BAM",
            Conversion::BedToBigBed(_) => "BED to bigBed",
            Conversion::CoolToHic(_) => "cool to hic",
        }
    }

    /// The format the conversion writes.
    pub fn target(&self) -> Format {
        match self {
            Conversion::FastaToBgzipFasta(_) => Format::BgzipFasta,
            Conversion::FastaToTwoBit(_) => Format::TwoBit,
            Conversion::VcfToBcf(_) => Format::Bcf,
            Conversion::SamToBam(_) => Format::Bam,
            Conversion::BedToBigBed(_) => Format::BigBed,
            Conversion::CoolToHic(_) => Format::Hic,
        }
    }

    pub fn src(&self) -> &Path {
        match self {
            Conversion::FastaToBgzipFasta(c) => &c.src,
            Conversion::FastaToTwoBit(c) => &c.src,
            Conversion::VcfToBcf(c) => &c.src,
            Conversion::SamToBam(c) => &c.src,
            Conversion::BedToBigBed(c) => &c.src,
            Conversion::CoolToHic(c) => &c.src,
        }
    }

    pub fn dst(&self) -> &Path {
        match self {
            Conversion::FastaToBgzipFasta(c) => &c.dst,
            Conversion::FastaToTwoBit(c) => &c.dst,
            Conversion::VcfToBcf(c) => &c.dst,
            Conversion::SamToBam(c) => &c.dst,
            Conversion::BedToBigBed(c) => &c.dst,
            Conversion::CoolToHic(c) => &c.dst,
        }
    }
}

/// A checked conversion with its memory division.
#[derive(Debug, Clone, PartialEq)]
pub struct Plan {
    pub conversion: Conversion,
    pub budget: Budget,
    pub memory: MemoryPlan,
}

/// Where a running conversion is.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Progress {
    /// What the conversion is doing (`"read"`, `"sort"`, `"write"`, ...).
    pub phase: &'static str,
    /// Input file bytes consumed so far.
    pub bytes_read: u64,
    /// Total input file bytes.
    pub bytes_total: u64,
}

/// What a finished conversion did.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Report {
    /// Files written with their formats, the destination first, then its
    /// companions.
    pub files: Vec<(Format, PathBuf)>,
    /// Records (sequences, variants, alignments, intervals, pixels)
    /// converted.
    pub records: u64,
    pub bytes_read: u64,
    /// Bytes written to spill files.
    pub spilled: u64,
    pub elapsed: Duration,
}

/// Progress callback.
pub type OnProgress<'a> = Option<&'a (dyn Fn(Progress) + Sync)>;

/// What a running conversion works with: its budget and memory division,
/// the worker pool, and where progress goes.
pub(crate) struct Context<'a> {
    pub budget: &'a Budget,
    pub memory: MemoryPlan,
    pub pool: Pool,
    progress: OnProgress<'a>,
    last: std::cell::Cell<Option<(Instant, &'static str)>>,
}

impl Context<'_> {
    /// Reports progress: on every change of phase, otherwise at most a few
    /// times a second.
    pub fn progress(&self, phase: &'static str, bytes_read: u64, bytes_total: u64) {
        let Some(f) = self.progress else { return };
        let now = Instant::now();
        if self
            .last
            .get()
            .is_some_and(|(t, p)| p == phase && now.duration_since(t) < Duration::from_millis(200))
        {
            return;
        }
        self.last.set(Some((now, phase)));
        f(Progress {
            phase,
            bytes_read,
            bytes_total,
        });
    }

    /// Directory for spill files.
    pub fn tmp_parent(&self, dst: &Path) -> PathBuf {
        match &self.budget.tmp_dir {
            Some(d) => d.clone(),
            None => dst
                .parent()
                .filter(|p| !p.as_os_str().is_empty())
                .map_or_else(|| PathBuf::from("."), Path::to_path_buf),
        }
    }
}

/// Checks `conversion` against `budget` and divides the memory.
pub fn plan(conversion: Conversion, budget: Budget) -> Result<Plan> {
    if budget.threads == 0 {
        return Err(Error::invalid("budget", "threads must be at least one"));
    }
    if !conversion.src().exists() {
        return Err(Error::Io {
            kind: std::io::ErrorKind::NotFound,
            detail: format!("{} does not exist", conversion.src().display()),
        });
    }
    if conversion.src() == conversion.dst() {
        return Err(Error::invalid(
            "convert",
            "source and destination are the same file",
        ));
    }
    let memory = match &conversion {
        Conversion::FastaToBgzipFasta(c) => c.plan(&budget)?,
        Conversion::FastaToTwoBit(c) => c.plan(&budget)?,
        Conversion::VcfToBcf(c) => c.plan(&budget)?,
        Conversion::SamToBam(c) => c.plan(&budget)?,
        Conversion::BedToBigBed(c) => c.plan(&budget)?,
        Conversion::CoolToHic(c) => c.plan(&budget)?,
    };
    let required = memory.fixed.checked_add(memory.pipeline).ok_or_else(|| {
        Error::invalid("budget", "fixed and pipeline memory exceed the size limit")
    })?;
    budget::require(&budget, required, "this conversion")?;
    Ok(Plan {
        conversion,
        budget,
        memory,
    })
}

/// Runs a plan.
pub fn run(plan: &Plan, progress: OnProgress<'_>) -> Result<Report> {
    let start = Instant::now();
    let ctx = Context {
        budget: &plan.budget,
        memory: plan.memory,
        pool: Pool::new(plan.budget.threads),
        progress,
        last: std::cell::Cell::new(None),
    };
    let mut report = match &plan.conversion {
        Conversion::FastaToBgzipFasta(c) => c.run(&ctx)?,
        Conversion::FastaToTwoBit(c) => c.run(&ctx)?,
        Conversion::VcfToBcf(c) => c.run(&ctx)?,
        Conversion::SamToBam(c) => c.run(&ctx)?,
        Conversion::BedToBigBed(c) => c.run(&ctx)?,
        Conversion::CoolToHic(c) => c.run(&ctx)?,
    };
    report.elapsed = start.elapsed();
    Ok(report)
}

/// [`plan`] then [`run`].
pub fn convert(conversion: Conversion, budget: Budget, progress: OnProgress<'_>) -> Result<Report> {
    run(&plan(conversion, budget)?, progress)
}
