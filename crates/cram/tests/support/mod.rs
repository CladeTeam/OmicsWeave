//! Shared fixtures and independent expectations for CRAM integration tests.

mod fixture;
mod htslib;
mod query;
mod sam;

pub use fixture::{Fixture, fasta_reader, fixture, reference, standard_refs, write_fasta};
pub use htslib::{Hts, Rec, crai_of, cram_of, read_all};
pub use query::{check_queries, only, open, query};
pub use sam::{SamSpec, random_sam};

use std::ffi::CString;
use std::path::Path;

/// A small deterministic generator (splitmix64).
pub struct Rng(u64);

impl Rng {
    pub fn new(seed: u64) -> Rng {
        Rng(seed)
    }

    pub fn next(&mut self) -> u64 {
        self.0 = self.0.wrapping_add(0x9E37_79B9_7F4A_7C15);
        let mut z = self.0;
        z = (z ^ (z >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
        z = (z ^ (z >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
        z ^ (z >> 31)
    }

    /// Uniform in `lo..hi`.
    pub fn range(&mut self, lo: u64, hi: u64) -> u64 {
        lo + self.next() % (hi - lo).max(1)
    }

    pub fn chance(&mut self, percent: u64) -> bool {
        self.next() % 100 < percent
    }
}

fn path_cstring(path: &Path) -> CString {
    CString::new(path.to_str().unwrap()).unwrap()
}

/// The text of a gzip file (one member or several).
pub fn gunzip(raw: &[u8]) -> Vec<u8> {
    use std::io::Read;
    let mut out = Vec::new();
    flate2::read::MultiGzDecoder::new(raw)
        .read_to_end(&mut out)
        .unwrap();
    out
}
