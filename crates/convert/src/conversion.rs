//! A conversion from one file to another: choosing it, checking it against
//! a budget, and running it.

use std::path::{Path, PathBuf};
use std::time::{Duration, Instant};

use omicsweave_bbi::WriteOptions;
use omicsweave_core::exec::Pool;
use omicsweave_core::{Error, Result};

use crate::bamsort::BamSort;
use crate::bed::{BedToBigBed, ChromSizes};
use crate::bedmethyl::{BedMethylToBigBed, BedMethylToBigWig, Strand, Value};
use crate::budget::{self, Budget, MemoryPlan};
use crate::cool::CoolToHic;
use crate::cram::CramToBam;
use crate::fasta::{FastaToBgzipFasta, FastaToTwoBit};
use crate::format::Format;
use crate::pairs::{PairsDialect, PairsToHic};
use crate::sam::{BamIndex, SamToBam};
use crate::vcf::VcfToBcf;
use crate::wig::{BedGraphToBigWig, WigToBigWig};

/// A conversion of one source file into one destination file (and its
/// companions), with its options.
#[derive(Debug, Clone, PartialEq)]
pub enum Conversion {
    FastaToBgzipFasta(FastaToBgzipFasta),
    FastaToTwoBit(FastaToTwoBit),
    VcfToBcf(VcfToBcf),
    SamToBam(SamToBam),
    BamSort(BamSort),
    CramToBam(CramToBam),
    BedToBigBed(BedToBigBed),
    BedGraphToBigWig(BedGraphToBigWig),
    WigToBigWig(WigToBigWig),
    CoolToHic(CoolToHic),
    PairsToHic(PairsToHic),
    BedMethylToBigBed(BedMethylToBigBed),
    BedMethylToBigWig(BedMethylToBigWig),
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
    /// Compression level (bgzip, BCF, BAM: 0-12; bigBed, bigWig, hic:
    /// zlib 1-9).
    pub level: Option<i32>,
    /// BCF, BAM: also write an index (records must be sorted).
    pub index: Option<Index>,
    /// bigBed, bigWig: chromosome lengths (`chrom.sizes`, `.fai` or
    /// `.2bit`) instead of each chromosome's largest interval end. pairs to
    /// hic: the chromosomes, in order, instead of the `#chromsize` lines.
    pub chrom_sizes: Option<PathBuf>,
    /// `.mcool` to hic: the resolutions to convert instead of all. pairs to
    /// hic: the resolutions to bin at instead of the Juicer `pre` defaults.
    pub resolutions: Option<Vec<u32>>,
    /// bedMethyl to bigWig: the column stored ([`Value::from_name`]).
    pub value: Option<String>,
    /// bedMethyl to bigWig: the items of this modification code (column 4).
    pub mod_code: Option<String>,
    /// bedMethyl to bigWig: the items of this strand (`+`, `-`, `.`).
    pub strand: Option<String>,
    /// bigBed, bigWig, hic: also write the `.sidecar`.
    pub sidecar: bool,
    /// VCF to BCF: declare contigs the header lacks (one more pass).
    pub add_contigs: bool,
    /// BAM to BAM: sort the records by coordinate, the one conversion of a
    /// BAM into a BAM (which therefore requires it).
    pub sort: bool,
    /// CRAM to BAM: the reference the CRAM was compressed against (FASTA
    /// with `.fai`, bgzip FASTA with `.fai` and `.gzi`, or 2bit).
    pub reference: Option<PathBuf>,
    /// BAM `.csi` (`index` [`Index::Csi`]): bins of `2^min_shift` bases
    /// instead of `2^14`.
    pub min_shift: Option<u32>,
    /// CRAM to BAM: check block CRC32s and each slice's reference MD5, as
    /// htslib does unless told to ignore them.
    pub verify_checksums: bool,
    /// bigBed, bigWig: how blocks, the index and zoom levels are written;
    /// `level` still sets the compression level when given.
    pub bbi: Option<WriteOptions>,
}

impl Default for Options {
    fn default() -> Options {
        Options {
            level: None,
            index: None,
            chrom_sizes: None,
            resolutions: None,
            value: None,
            mod_code: None,
            strand: None,
            sidecar: false,
            add_contigs: false,
            sort: false,
            reference: None,
            min_shift: None,
            verify_checksums: true,
            bbi: None,
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
            (Bam, Bam) => Conversion::BamSort(BamSort::new(src, dst)),
            (Cram, Bam) => Conversion::CramToBam(CramToBam::new(src, dst)),
            // A bedGraph is BED with four fields.
            (Bed | GzipBed | BedGraph | GzipBedGraph, BigBed) => {
                Conversion::BedToBigBed(BedToBigBed::new(src, dst, ChromSizes::FromData))
            }
            (BedGraph | GzipBedGraph, BigWig) => {
                Conversion::BedGraphToBigWig(BedGraphToBigWig::new(src, dst, ChromSizes::FromData))
            }
            (Wig | GzipWig, BigWig) => {
                Conversion::WigToBigWig(WigToBigWig::new(src, dst, ChromSizes::FromData))
            }
            (Cool | Mcool, Hic) => Conversion::CoolToHic(CoolToHic::new(src, dst)),
            (Pairs | GzipPairs, Hic) => {
                Conversion::PairsToHic(PairsToHic::new(src, dst, PairsDialect::Pairs))
            }
            (JuicerShort | GzipJuicerShort, Hic) => {
                Conversion::PairsToHic(PairsToHic::new(src, dst, PairsDialect::Short))
            }
            (BedMethyl | GzipBedMethyl, BigBed) => Conversion::BedMethylToBigBed(
                BedMethylToBigBed::new(src, dst, ChromSizes::FromData),
            ),
            (BedMethyl | GzipBedMethyl, BigWig) => Conversion::BedMethylToBigWig(
                BedMethylToBigWig::new(src, dst, ChromSizes::FromData),
            ),
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
            matches!(
                self,
                Conversion::VcfToBcf(_)
                    | Conversion::SamToBam(_)
                    | Conversion::BamSort(_)
                    | Conversion::CramToBam(_)
            ),
            "index",
            self,
        )?;
        takes(
            o.chrom_sizes.is_some(),
            matches!(
                self,
                Conversion::BedToBigBed(_)
                    | Conversion::BedGraphToBigWig(_)
                    | Conversion::WigToBigWig(_)
                    | Conversion::PairsToHic(_)
                    | Conversion::BedMethylToBigBed(_)
                    | Conversion::BedMethylToBigWig(_)
            ),
            "chrom_sizes",
            self,
        )?;
        takes(
            o.resolutions.is_some(),
            matches!(self, Conversion::CoolToHic(_) | Conversion::PairsToHic(_)),
            "resolutions",
            self,
        )?;
        let methyl = matches!(self, Conversion::BedMethylToBigWig(_));
        takes(o.value.is_some(), methyl, "value", self)?;
        takes(o.mod_code.is_some(), methyl, "mod", self)?;
        takes(o.strand.is_some(), methyl, "strand", self)?;
        takes(
            o.add_contigs,
            matches!(self, Conversion::VcfToBcf(_)),
            "add_contigs",
            self,
        )?;
        takes(o.sort, matches!(self, Conversion::BamSort(_)), "sort", self)?;
        takes(
            o.reference.is_some(),
            matches!(self, Conversion::CramToBam(_)),
            "reference",
            self,
        )?;
        takes(
            !o.verify_checksums,
            matches!(self, Conversion::CramToBam(_)),
            "verify_checksums",
            self,
        )?;
        let bam = matches!(
            self,
            Conversion::SamToBam(_) | Conversion::BamSort(_) | Conversion::CramToBam(_)
        );
        takes(
            o.min_shift.is_some(),
            bam && o.index == Some(Index::Csi),
            "min_shift",
            self,
        )?;
        takes(
            o.bbi.is_some(),
            matches!(
                self,
                Conversion::BedToBigBed(_)
                    | Conversion::BedGraphToBigWig(_)
                    | Conversion::WigToBigWig(_)
                    | Conversion::BedMethylToBigBed(_)
                    | Conversion::BedMethylToBigWig(_)
            ),
            "bbi",
            self,
        )?;
        takes(
            o.sidecar,
            matches!(
                self,
                Conversion::BedToBigBed(_)
                    | Conversion::BedGraphToBigWig(_)
                    | Conversion::WigToBigWig(_)
                    | Conversion::BedMethylToBigBed(_)
                    | Conversion::BedMethylToBigWig(_)
                    | Conversion::CoolToHic(_)
                    | Conversion::PairsToHic(_)
            ),
            "sidecar",
            self,
        )?;
        let bam_index = |index: Option<Index>| match index {
            None => BamIndex::None,
            Some(Index::Auto) => BamIndex::Auto,
            Some(Index::Bai) => BamIndex::Bai,
            Some(Index::Csi) => BamIndex::Csi {
                min_shift: o.min_shift.unwrap_or(14),
            },
        };
        let bbi = |options: &mut WriteOptions| {
            if let Some(b) = &o.bbi {
                *options = b.clone();
            }
            options.level = o.level.unwrap_or(options.level);
        };
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
                c.index = bam_index(o.index);
            }
            Conversion::CramToBam(c) => {
                c.level = o.level.unwrap_or(c.level);
                c.reference = o.reference.clone();
                c.verify_checksums = o.verify_checksums;
                c.index = bam_index(o.index);
            }
            Conversion::BamSort(c) => {
                if !o.sort {
                    return Err(Error::invalid(
                        "convert",
                        if o.index.is_some() {
                            "an index needs a coordinate-sorted BAM: pass sort to sort this one \
                             (a BAM already sorted is indexed in place, not converted)"
                        } else {
                            "a BAM converts to a BAM only by sorting it: pass sort"
                        },
                    ));
                }
                c.level = o.level.unwrap_or(c.level);
                c.index = bam_index(o.index);
            }
            Conversion::BedToBigBed(c) => {
                if let Some(p) = &o.chrom_sizes {
                    c.chrom_sizes = ChromSizes::File(p.clone());
                }
                bbi(&mut c.options);
                c.sidecar = o.sidecar;
            }
            Conversion::BedGraphToBigWig(c) => {
                if let Some(p) = &o.chrom_sizes {
                    c.chrom_sizes = ChromSizes::File(p.clone());
                }
                bbi(&mut c.options);
                c.sidecar = o.sidecar;
            }
            Conversion::WigToBigWig(c) => {
                if let Some(p) = &o.chrom_sizes {
                    c.chrom_sizes = ChromSizes::File(p.clone());
                }
                bbi(&mut c.options);
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
            Conversion::PairsToHic(c) => {
                if let Some(l) = o.level {
                    c.level = u32::try_from(l)
                        .map_err(|_| Error::invalid("convert", "hic level must be 1 to 9"))?;
                }
                c.chrom_sizes = o.chrom_sizes.clone();
                c.resolutions = o.resolutions.clone();
                c.sidecar = o.sidecar;
            }
            Conversion::BedMethylToBigBed(c) => {
                if let Some(p) = &o.chrom_sizes {
                    c.chrom_sizes = ChromSizes::File(p.clone());
                }
                bbi(&mut c.options);
                c.sidecar = o.sidecar;
            }
            Conversion::BedMethylToBigWig(c) => {
                if let Some(p) = &o.chrom_sizes {
                    c.chrom_sizes = ChromSizes::File(p.clone());
                }
                if let Some(v) = &o.value {
                    c.value = Value::from_name(v)?;
                }
                c.mod_code = o.mod_code.clone();
                c.strand = o.strand.as_deref().map(Strand::from_name).transpose()?;
                bbi(&mut c.options);
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
            Conversion::BamSort(_) => "BAM to sorted BAM",
            Conversion::CramToBam(_) => "CRAM to BAM",
            Conversion::BedToBigBed(_) => "BED to bigBed",
            Conversion::BedGraphToBigWig(_) => "bedGraph to bigWig",
            Conversion::WigToBigWig(_) => "wig to bigWig",
            Conversion::CoolToHic(_) => "cool to hic",
            Conversion::PairsToHic(c) => match c.dialect {
                PairsDialect::Pairs => "pairs to hic",
                PairsDialect::Short => "Juicer short to hic",
            },
            Conversion::BedMethylToBigBed(_) => "bedMethyl to bigBed",
            Conversion::BedMethylToBigWig(_) => "bedMethyl to bigWig",
        }
    }

    /// The format the conversion writes.
    pub fn target(&self) -> Format {
        match self {
            Conversion::FastaToBgzipFasta(_) => Format::BgzipFasta,
            Conversion::FastaToTwoBit(_) => Format::TwoBit,
            Conversion::VcfToBcf(_) => Format::Bcf,
            Conversion::SamToBam(_) | Conversion::BamSort(_) | Conversion::CramToBam(_) => {
                Format::Bam
            }
            Conversion::BedToBigBed(_) | Conversion::BedMethylToBigBed(_) => Format::BigBed,
            Conversion::BedGraphToBigWig(_)
            | Conversion::WigToBigWig(_)
            | Conversion::BedMethylToBigWig(_) => Format::BigWig,
            Conversion::CoolToHic(_) | Conversion::PairsToHic(_) => Format::Hic,
        }
    }

    pub fn src(&self) -> &Path {
        match self {
            Conversion::FastaToBgzipFasta(c) => &c.src,
            Conversion::FastaToTwoBit(c) => &c.src,
            Conversion::VcfToBcf(c) => &c.src,
            Conversion::SamToBam(c) => &c.src,
            Conversion::BamSort(c) => &c.src,
            Conversion::CramToBam(c) => &c.src,
            Conversion::BedToBigBed(c) => &c.src,
            Conversion::BedGraphToBigWig(c) => &c.src,
            Conversion::WigToBigWig(c) => &c.src,
            Conversion::CoolToHic(c) => &c.src,
            Conversion::PairsToHic(c) => &c.src,
            Conversion::BedMethylToBigBed(c) => &c.src,
            Conversion::BedMethylToBigWig(c) => &c.src,
        }
    }

    pub fn dst(&self) -> &Path {
        match self {
            Conversion::FastaToBgzipFasta(c) => &c.dst,
            Conversion::FastaToTwoBit(c) => &c.dst,
            Conversion::VcfToBcf(c) => &c.dst,
            Conversion::SamToBam(c) => &c.dst,
            Conversion::BamSort(c) => &c.dst,
            Conversion::CramToBam(c) => &c.dst,
            Conversion::BedToBigBed(c) => &c.dst,
            Conversion::BedGraphToBigWig(c) => &c.dst,
            Conversion::WigToBigWig(c) => &c.dst,
            Conversion::CoolToHic(c) => &c.dst,
            Conversion::PairsToHic(c) => &c.dst,
            Conversion::BedMethylToBigBed(c) => &c.dst,
            Conversion::BedMethylToBigWig(c) => &c.dst,
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
    /// What the conversion passed over without failing (rows of
    /// chromosomes it does not know, say), one message each.
    pub warnings: Vec<String>,
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
    warnings: std::cell::RefCell<Vec<String>>,
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

    /// Notes something passed over, for [`Report::warnings`].
    pub fn warn(&self, message: String) {
        self.warnings.borrow_mut().push(message);
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
        Conversion::BamSort(c) => c.plan(&budget)?,
        Conversion::CramToBam(c) => c.plan(&budget)?,
        Conversion::BedToBigBed(c) => c.plan(&budget)?,
        Conversion::BedGraphToBigWig(c) => c.plan(&budget)?,
        Conversion::WigToBigWig(c) => c.plan(&budget)?,
        Conversion::CoolToHic(c) => c.plan(&budget)?,
        Conversion::PairsToHic(c) => c.plan(&budget)?,
        Conversion::BedMethylToBigBed(c) => c.plan(&budget)?,
        Conversion::BedMethylToBigWig(c) => c.plan(&budget)?,
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
        warnings: Default::default(),
    };
    let mut report = match &plan.conversion {
        Conversion::FastaToBgzipFasta(c) => c.run(&ctx)?,
        Conversion::FastaToTwoBit(c) => c.run(&ctx)?,
        Conversion::VcfToBcf(c) => c.run(&ctx)?,
        Conversion::SamToBam(c) => c.run(&ctx)?,
        Conversion::BamSort(c) => c.run(&ctx)?,
        Conversion::CramToBam(c) => c.run(&ctx)?,
        Conversion::BedToBigBed(c) => c.run(&ctx)?,
        Conversion::BedGraphToBigWig(c) => c.run(&ctx)?,
        Conversion::WigToBigWig(c) => c.run(&ctx)?,
        Conversion::CoolToHic(c) => c.run(&ctx)?,
        Conversion::PairsToHic(c) => c.run(&ctx)?,
        Conversion::BedMethylToBigBed(c) => c.run(&ctx)?,
        Conversion::BedMethylToBigWig(c) => c.run(&ctx)?,
    };
    report.elapsed = start.elapsed();
    report.warnings = ctx.warnings.take();
    Ok(report)
}

/// [`plan`] then [`run`].
pub fn convert(conversion: Conversion, budget: Budget, progress: OnProgress<'_>) -> Result<Report> {
    run(&plan(conversion, budget)?, progress)
}
