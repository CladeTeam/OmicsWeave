//! Convert-specific CRAM input and independent htslib output.

mod fixtures;
mod oracle;

pub use fixtures::{SamSpec, random_sam, reference, write_fasta};
pub use oracle::{bam_of_cram, crai_of, cram_of};
