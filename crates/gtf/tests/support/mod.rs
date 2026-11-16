//! Shared input fixtures, generation and independent base-set expectations.

pub mod base_sets;
pub mod fixtures;
pub mod generator;

use omicsweave_gtf::{Annotation, GtfLayout, ReadOptions};

pub fn read(text: &str) -> Annotation {
    Annotation::parse(text.as_bytes(), &ReadOptions::default()).unwrap()
}

pub fn gtf_text(ann: &Annotation, layout: GtfLayout) -> String {
    let mut out = Vec::new();
    ann.write_gtf(&mut out, layout).unwrap();
    String::from_utf8(out).unwrap()
}

pub const LEN: u32 = 3000;

pub fn random(seed: u64) -> Annotation {
    read(&generator::random_gtf(
        &mut generator::Rng(seed),
        3,
        40,
        LEN,
    ))
}
