//! Shared helpers of the conversion tests.
#![allow(dead_code)]

use std::io::Write;
use std::path::Path;

use omicsweave_convert::Budget;

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

    pub fn pick<'a, T>(&mut self, items: &'a [T]) -> &'a T {
        &items[self.range(0, items.len() as u64) as usize]
    }
}

/// The three storages an input can have.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Storage {
    Plain,
    Gzip,
    Bgzf,
}

pub const STORAGES: [Storage; 3] = [Storage::Plain, Storage::Gzip, Storage::Bgzf];

/// Writes `text` to `path` in `storage`.
pub fn write_input(path: &Path, text: &[u8], storage: Storage) {
    let bytes = match storage {
        Storage::Plain => text.to_vec(),
        Storage::Gzip => {
            let mut e = flate2::write::GzEncoder::new(Vec::new(), flate2::Compression::new(3));
            e.write_all(text).unwrap();
            e.finish().unwrap()
        }
        Storage::Bgzf => {
            let mut w = omicsweave_bgzf::Writer::new(
                Vec::new(),
                1,
                Box::new(omicsweave_core::exec::Inline::new()),
            )
            .unwrap();
            w.write(text).unwrap();
            w.finish().unwrap().out
        }
    };
    std::fs::write(path, bytes).unwrap();
}

/// Budgets that exercise the inline and the threaded paths, roomy and
/// tight.
pub fn budgets() -> Vec<Budget> {
    vec![
        Budget::new(1, 1 << 30),
        Budget::new(6, 1 << 30),
        Budget::new(3, 6 << 20),
    ]
}

pub mod hts;

/// The text of a BGZF file.
pub fn inflate(path: &Path) -> Vec<u8> {
    use std::io::Read;
    let raw = std::fs::read(path).unwrap();
    let mut out = Vec::new();
    omicsweave_bgzf::Reader::new(&raw[..], Box::new(omicsweave_core::exec::Inline::new()))
        .read_to_end(&mut out)
        .unwrap();
    out
}
