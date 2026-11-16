//! Streaming, parallel format conversion under a budget of threads and
//! memory.
//!
//! | [`Conversion`] | Source ([`Format`]) | Destination |
//! |---|---|---|
//! | [`FastaToBgzipFasta`] | FASTA (plain, gzip or BGZF) | bgzip FASTA + `.fai` + `.gzi` |
//! | [`FastaToTwoBit`] | FASTA (plain, gzip or BGZF) | `.2bit` |
//! | [`VcfToBcf`] | VCF (plain, gzip or BGZF) | `.bcf` (+ `.csi`) |
//! | [`SamToBam`] | SAM (plain, gzip or BGZF) | `.bam` (+ `.bai` / `.csi`) |
//! | [`BamSort`] | BAM | coordinate-sorted `.bam` (+ `.bai` / `.csi`) |
//! | [`CramToBam`] | CRAM 3.0 / 3.1 (+ a reference) | `.bam` (+ `.bai` / `.csi`) |
//! | [`BedToBigBed`] | BED (plain, gzip or BGZF) | `.bigbed` + `.sidecar` |
//! | [`BedGraphToBigWig`] | bedGraph (plain, gzip or BGZF) | `.bigwig` + `.sidecar` |
//! | [`WigToBigWig`] | wig (plain, gzip or BGZF) | `.bigwig` + `.sidecar` |
//! | [`CoolToHic`] | `.cool` / `.mcool` | `.hic` + `.sidecar` |
//! | [`PairsToHic`] | 4DN pairs, Juicer short (plain, gzip or BGZF) | `.hic` + `.sidecar` |
//! | [`BedMethylToBigBed`] | bedMethyl (plain, gzip or BGZF) | `.bigbed` + `.sidecar` |
//! | [`BedMethylToBigWig`] | bedMethyl (plain, gzip or BGZF) | `.bigwig` + `.sidecar` |
//! | [`MafToBigMaf`] | MAF (plain, gzip or BGZF) | bigMaf `.bigbed` + `.sidecar` |
//! | [`MotifToFfm`] | motif collections (MEME, JASPAR, TRANSFAC, HOMER, CIS-BP, HOCOMOCO, Cluster-Buster, BaMM, ProBound, GLAM2, FFM; one or more files) | `.ffm` |
//!
//! [`Conversion::between`] chooses the conversion for a pair of formats
//! ([`Format::detect`] tells a source's format by content,
//! [`Format::from_extension`] a destination's by name), [`Conversion::infer`]
//! for a pair of file names; [`Conversion::configure`] applies [`Options`].
//!
//! Every conversion is built from the same three kinds of stage
//! (`pipeline`): a serial reader that cuts the source into units (line
//! batches, BGZF blocks, pixel spans), pure tasks that decode or encode
//! units on a shared [`omicsweave_core::exec::Pool`], and a serial writer
//! that takes results back **in order**. Output therefore does not depend
//! on the number of threads, nor on whether a sort spilled to disk.
//!
//! [`plan`] checks a conversion against its [`Budget`] before anything is
//! read and reports how the memory is divided ([`MemoryPlan`]); [`run`]
//! executes the plan. Files are published at the end through unique
//! temporary files and atomic renames per file; several files are not a
//! filesystem transaction.

mod bamsort;
mod bed;
mod bedmethyl;
mod budget;
mod conversion;
mod cool;
mod cram;
mod fasta;
mod format;
mod maf;
mod motif;
mod pairs;
mod pipeline;
mod sam;
mod vcf;
mod wig;

pub use bamsort::BamSort;
pub use bed::{BedToBigBed, ChromSizes};
pub use bedmethyl::{BedMethylToBigBed, BedMethylToBigWig, Strand, Value};
pub use budget::{Budget, MemoryPlan, RSS_SLACK, human, parse_size};
pub use conversion::{
    Conversion, Index, OnProgress, Options, Plan, Progress, Report, convert, plan, run,
};
pub use cool::CoolToHic;
pub use cram::CramToBam;
pub use fasta::{FastaToBgzipFasta, FastaToTwoBit};
pub use format::Format;
pub use maf::{BIGMAF_AUTOSQL, MafToBigMaf};
pub use motif::{MotifToFfm, detect_all as detect_sources};
pub use pairs::{JUICER_RESOLUTIONS, PairsDialect, PairsToHic};
pub use sam::{BamIndex, SamToBam};
pub use vcf::VcfToBcf;
pub use wig::{BedGraphToBigWig, WigToBigWig};

/// Pipeline stages shared with other budgeted streaming tools (the
/// quantification layer): queue sizing, the external sorter, and
/// destinations published atomically.
pub mod stage {
    pub use crate::budget::{Sizing, require, sizing, text_batch};
    pub use crate::pipeline::output::{Outputs, TempDir};
    pub use crate::pipeline::sort::{Merge, Sorter, Spill, read_or_end};
}
