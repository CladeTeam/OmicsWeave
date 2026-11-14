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
        "fixture {name} input_bytes={} input_checksum={:016x} output_bytes={} output_checksum={:016x}",
        line.len(),
        omicsweave_core::sidecar::checksum(line),
        out.len(),
        omicsweave_core::sidecar::checksum(&out),
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
    measure(
        "sam_plain",
        &encoder,
        b"q\t0\tchr1\t1\t60\t8M\t*\t0\t0\tACGTACGT\tIIIIIIII\tNM:i:0",
    );
    for (name, fake, real, before, after) in [
        ("sam_cg_small", 1, 1, 0, 0),
        ("sam_cg_inline", 1, 65_535, 0, 0),
        ("sam_cg_long", 1, 65_536, 0, 0),
        ("sam_cg_grow_aux", 1, 128, 131_072, 131_072),
        ("sam_cg_shrink_aux", 16, 65_536, 131_072, 131_072),
    ] {
        let mut line = b"q\t0\tchr1\t1\t60\t8S".to_vec();
        for _ in 1..fake {
            line.extend_from_slice(b"0N");
        }
        line.extend_from_slice(b"\t*\t0\t0\tACGTACGT\tIIIIIIII\tAA:Z:");
        line.extend(std::iter::repeat_n(b'a', before));
        line.extend_from_slice(b"\tCG:B:I,128");
        for _ in 1..real {
            line.extend_from_slice(b",2");
        }
        line.extend_from_slice(b"\tZZ:Z:");
        line.extend(std::iter::repeat_n(b'z', after));
        measure(name, &encoder, &line);
    }
}
