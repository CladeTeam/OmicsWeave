//! MM/ML/MN fixtures and adapters for the independent htslib oracle.
mod fixture;
mod oracle;
pub use fixture::{random_sam, spec_lines, write_bam};
pub use oracle::{HtsMod, as_hts, htslib_mods, our_mods, types_as_hts};
