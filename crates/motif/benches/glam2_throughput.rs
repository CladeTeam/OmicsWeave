//! Run on a worker with cargo bench -p omicsweave-motif --bench glam2_throughput.
#[allow(dead_code)]
#[path = "../tests/support/glam2.rs"]
mod glam2_support;
#[allow(dead_code)]
#[path = "../tests/support/rng.rs"]
mod rng;
use glam2_support::*;
use omicsweave_motif::Strands;
use omicsweave_motif::profile::{GLAM2_HITS, ProfileScoring, glam2_dna_alphabet};

fn main() {
    // Cargo passes --bench only for benchmarking; all-target test runs stay lightweight.
    if !std::env::args().any(|arg| arg == "--bench") {
        return;
    }
    let alphabet = glam2_dna_alphabet();
    let mut rng = Rng(4242);
    let motif_text = loop {
        let text = random_motif(&mut rng, b"acgt", 16, 20);
        if text.lines().nth(2).unwrap().matches('*').count() >= 12 {
            break text;
        }
    };
    let motif = parse_glam2_motif(&motif_text, &alphabet);
    let fasta_text = random_fasta(&mut rng, b"acgt", 400, 10_000);
    let fasta = Fasta::parse(&fasta_text);
    let scoring = ProfileScoring::glam2(&alphabet).unwrap();
    let p = program(
        std::slice::from_ref(&motif),
        &alphabet,
        &scoring,
        Strands::Both,
    )
    .unwrap();
    let bases = fasta.data.len();
    for threads in [Some(1), None] {
        let t = std::time::Instant::now();
        let text = hits_text(&p, &fasta, GLAM2_HITS, threads);
        let secs = t.elapsed().as_secs_f64();
        eprintln!(
            "threads {threads:?}: {bases} bases, width {}, both strands: {secs:.3} s ({:.1} Mb/s), {} lines",
            p.model(0).width(),
            bases as f64 / secs / 1e6,
            text.lines().count()
        );
    }
}
