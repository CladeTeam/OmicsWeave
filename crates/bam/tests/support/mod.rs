//! Fixtures, input inspection and independent htslib oracles.
mod fixture;
mod htslib;
mod io;
pub mod modifications;
mod query;
mod sam;

pub use fixture::fixture;
pub use htslib::{Hts, Rec, bam_of, index, read_all};
pub use io::{blocks, parse_index, with_empty_blocks};
pub use query::{open, query_all, regions};
pub use sam::{Rng, SamSpec, random_sam};
