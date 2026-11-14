//! Streaming, parallel format conversion under a budget of threads and
//! memory.
//!
//! | [`Conversion`] | Source ([`Format`]) | Destination |
//! |---|---|---|
//! | [`FastaToBgzipFasta`] | FASTA (plain, gzip or BGZF) | bgzip FASTA + `.fai` + `.gzi` |
//! | [`FastaToTwoBit`] | FASTA (plain, gzip or BGZF) | `.2bit` |
//! | [`VcfToBcf`] | VCF (plain, gzip or BGZF) | `.bcf` (+ `.csi`) |
//! | [`SamToBam`] | SAM (plain, gzip or BGZF) | `.bam` (+ `.bai` / `.csi`) |
//! | [`BedToBigBed`] | BED (plain, gzip or BGZF) | `.bigbed` + `.sidecar` |
//! | [`CoolToHic`] | `.cool` / `.mcool` | `.hic` + `.sidecar` |
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

mod bed;
mod budget;
mod conversion;
mod cool;
mod fasta;
mod format;
mod pipeline;
mod sam;
mod vcf;

pub use bed::{BedToBigBed, ChromSizes};
pub use budget::{Budget, MemoryPlan, RSS_SLACK, human, parse_size};
pub use conversion::{
    Conversion, Index, OnProgress, Options, Plan, Progress, Report, convert, plan, run,
};
pub use cool::CoolToHic;
pub use fasta::{FastaToBgzipFasta, FastaToTwoBit};
pub use format::Format;
pub use sam::{BamIndex, SamToBam};
pub use vcf::VcfToBcf;
