//! Personal consensus sequences: `bcftools consensus` (1.23) in Rust,
//! producing the same sequences in memory.
//!
//! A [`Consensus`] reads a BCF by region through its CSI index
//! ([`omicsweave_vcf::BcfReader`]). Loading a region reads its reference,
//! through any [`omicsweave_core::SequenceRead`], and decodes the records
//! it takes, each once; it then builds sequences of it, one per query: a
//! sample's haplotype, IUPAC codes of samples, or every ALT. Loading once
//! and querying many samples or haplotypes is the batch path. Every
//! `bcftools consensus` option applies: haplotypes (`-H`), IUPAC codes
//! (`-I`), missing and absent characters (`-M`, `-a`), masks (`-m`,
//! `--mask-with`), marks (`--mark-*`), chains (`-c`), `--regions-overlap`
//! and `-i`/`-e` filters ([`RecordFilter`], applied as a region's records
//! are decoded). Sequences are written 5' to 3' of a strand: the reverse
//! complement for `-`.

mod chain;
mod consensus;
mod edit;
mod iupac;
mod mask;
mod records;
mod select;

pub use consensus::{Alleles, Consensus, ConsensusOptions, Mark};
pub use mask::{Mask, MaskWith};
pub use records::{Overlap, RecordFilter};
pub use select::is_haplotype;
