#![allow(dead_code, unused_imports)]

mod bcftools;
mod fixtures;
mod oracle;
mod random;

pub use bcftools::{Outcome, bcftools_consensus};
pub use fixtures::{Genome, bcf_of, consensus, open_bcf, to_indexed_bcf, to_indexed_raw_bcf};
pub use oracle::fasta_sequence;
pub use random::Rng;
