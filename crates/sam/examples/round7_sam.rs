use std::hint::black_box;
use std::time::{Duration, Instant};

use omicsweave_sam::{BamEncoder, BamScratch, SamHeader};

fn measure(name: &str, encoder: &BamEncoder, line: &[u8]) {
    let mut scratch = BamScratch::new();
    let mut out = Vec::new();
    encoder
        .encode_line(&mut scratch, line, 1, &mut out)
        .unwrap();
    let expected = out.clone();
    println!(
        "fixture {name} bytes={} checksum={:016x}",
        expected.len(),
        omicsweave_core::sidecar::checksum(&expected),
    );
    for _ in 0..5 {
        let start = Instant::now();
        let mut iterations = 0u64;
        while start.elapsed() < Duration::from_millis(100) {
            out.clear();
            black_box(
                encoder
                    .encode_line(&mut scratch, black_box(line), 1, &mut out)
                    .unwrap(),
            );
            black_box(&out);
            iterations += 1;
        }
        let nanos = start.elapsed().as_nanos() as f64 / iterations as f64;
        assert_eq!(out, expected);
        println!("{name} {nanos:.3} ns/iter");
    }
}

fn main() {
    let mut header = &b"@SQ\tSN:chr1\tLN:2147483647\n"[..];
    let encoder = BamEncoder::new(&SamHeader::read(&mut header).unwrap());
    for n in [1usize, 128, 65_536] {
        let line = format!(
            "q\t0\tchr1\t1\t60\t{}\t*\t0\t0\t{}\t*\tNM:i:0\tMD:Z:{n}",
            "1M".repeat(n),
            "A".repeat(n),
        );
        measure(&format!("sam_cigar_{n}"), &encoder, line.as_bytes());
    }
    for n in [0usize, 16, 4096] {
        let line = format!(
            "q\t0\tchr1\t1\t60\t100S\t*\t0\t0\t{}\t*{}\tCG:B:I,1600\tXZ:Z:after",
            "A".repeat(100),
            "\tAA:Z:abcdefghijklmnopqrstuvwx".repeat(n),
        );
        measure(&format!("sam_cg_aux_{n}"), &encoder, line.as_bytes());
        let ordinary = line.replace("100S", "100M");
        measure(&format!("sam_plain_aux_{n}"), &encoder, ordinary.as_bytes());
    }
}
