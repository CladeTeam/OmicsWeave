//! Personal consensus sequences: `bcftools consensus` (1.23) in Rust,
//! producing the same sequences in memory.
//!
//! A [`VcfStore`] reads a VCF or BCF file once into a resident form shared
//! by every builder. A [`Consensus`] loads a region — its reference, read
//! through any [`omicsweave_core::SequenceRead`], and the records it takes —
//! and then builds sequences of it, one per query: a sample's haplotype,
//! IUPAC codes of samples, or every ALT. Loading once and querying many
//! samples or haplotypes is the batch path. Every `bcftools consensus`
//! option applies: haplotypes (`-H`), IUPAC codes (`-I`), missing and
//! absent characters (`-M`, `-a`), masks (`-m`, `--mask-with`), marks
//! (`--mark-*`), chains (`-c`), `--regions-overlap` and `-i`/`-e` filters
//! ([`RecordFilter`], applied as the store is read).

mod chain;
mod consensus;
mod edit;
mod iupac;
mod mask;
mod select;
mod store;

pub use consensus::{Alleles, Consensus, ConsensusOptions, Mark};
pub use mask::{Mask, MaskWith};
pub use store::{Overlap, RecordFilter, VcfStore};
