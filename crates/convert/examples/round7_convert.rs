//! Shared-source baseline/candidate fixture; run only on an authorized worker.
#[allow(dead_code)]
#[path = "../src/pipeline/input.rs"]
mod input;

use std::hint::black_box;
use std::io::BufReader;
use std::time::{Duration, Instant};

fn bench(name: &str, mut operation: impl FnMut()) {
    operation();
    for sample in 0..5 {
        let start = Instant::now();
        let mut iterations = 0u64;
        while start.elapsed() < Duration::from_millis(100) {
            operation();
            iterations += 1;
        }
        println!(
            "{name} sample={sample} ns/iter={:.0} iterations={iterations}",
            start.elapsed().as_nanos() as f64 / iterations as f64
        );
    }
}

fn verify(name: &str, input: &[u8], capacity: usize, target: usize) {
    let reader = BufReader::with_capacity(capacity, input);
    let mut batches = input::LineBatches::new(reader, target, 1);
    let mut text = Vec::new();
    while let Some(batch) = batches.next_batch().unwrap() {
        text.extend_from_slice(&batch.text);
    }
    assert_eq!(text, input);
    println!(
        "fixture {name} bytes={} checksum={:016x}",
        text.len(),
        omicsweave_core::sidecar::checksum(&text)
    );
}

fn main() {
    let long = vec![b'x'; 1 << 20];
    verify("line_batches_1mib_single_line", &long, 4096, 1024);
    bench("line_batches_1mib_single_line", || {
        let reader = BufReader::with_capacity(4096, black_box(long.as_slice()));
        let mut batches = input::LineBatches::new(reader, 1024, 1);
        let mut bytes = 0;
        while let Some(batch) = batches.next_batch().unwrap() {
            bytes += batch.text.len();
            black_box(batch);
        }
        assert_eq!(black_box(bytes), long.len());
    });
    let ordinary = b"chr1\t123\t456\tfeature\n".repeat(50_000);
    verify("line_batches_ordinary_lines", &ordinary, 64 << 10, 64 << 10);
    bench("line_batches_ordinary_lines", || {
        let reader = BufReader::with_capacity(64 << 10, black_box(ordinary.as_slice()));
        let mut batches = input::LineBatches::new(reader, 64 << 10, 1);
        let mut bytes = 0;
        while let Some(batch) = batches.next_batch().unwrap() {
            bytes += batch.text.len();
            black_box(batch);
        }
        assert_eq!(black_box(bytes), ordinary.len());
    });
}
