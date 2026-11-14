//! Compression boundary checks and a paired public-API workload.

use omicsweave_bgzf::{
    BLOCK_TEXT_LEN, EOF_BLOCK, MAX_BLOCK_LEN, Writer, check_level, deflate_block,
};
use omicsweave_core::exec::Inline;

fn text(len: usize, random: bool) -> Vec<u8> {
    let mut state = 0x1234_5678u32;
    (0..len)
        .map(|i| {
            if random {
                state ^= state << 13;
                state ^= state >> 17;
                state ^= state << 5;
                state as u8
            } else {
                b"ACGT\n"[i % 5]
            }
        })
        .collect()
}

fn checksum(bytes: &[u8]) -> u64 {
    bytes.iter().fold(0xcbf2_9ce4_8422_2325, |h, &b| {
        (h ^ u64::from(b)).wrapping_mul(0x100_0000_01b3)
    })
}

// Use libdeflate's gzip decoder directly: this independently checks the
// BGZF framing, CRC and ISIZE instead of the crate's own inflate path.
fn decode_gzip_block(raw: &[u8], expected: &[u8]) {
    assert!((28..=MAX_BLOCK_LEN).contains(&raw.len()));
    assert_eq!(&raw[..16], &EOF_BLOCK[..16]);
    assert_eq!(
        u16::from_le_bytes([raw[16], raw[17]]) as usize + 1,
        raw.len()
    );
    let mut decoded = vec![0; expected.len()];
    let n = libdeflater::Decompressor::new()
        .gzip_decompress(raw, &mut decoded)
        .unwrap();
    assert_eq!(n, expected.len());
    assert_eq!(decoded, expected);
}

#[test]
fn compression_levels_sizes_and_append_preserve_gzip_members() {
    for level in 0..=12 {
        let mut compressor = libdeflater::Compressor::new(check_level(level).unwrap());
        for len in [0, 1, 64, 4096, BLOCK_TEXT_LEN] {
            for random in [false, true] {
                let expected = text(len, random);
                let mut raw = b"preserved-prefix".to_vec();
                let base = raw.len();
                deflate_block(&mut compressor, &expected, &mut raw).unwrap();
                assert_eq!(&raw[..base], b"preserved-prefix");
                decode_gzip_block(&raw[base..], &expected);
            }
        }
    }
    let mut compressor = libdeflater::Compressor::new(check_level(6).unwrap());
    let mut raw = b"preserved-prefix".to_vec();
    assert!(deflate_block(&mut compressor, &text(BLOCK_TEXT_LEN + 1, false), &mut raw).is_err());
    assert_eq!(raw, b"preserved-prefix");
}

fn write_records(input: &[u8], blocks: usize) -> Vec<u8> {
    let mut writer = Writer::new(Vec::new(), 6, Box::new(Inline::new())).unwrap();
    for _ in 0..blocks {
        writer.write(input).unwrap();
        writer.flush_block().unwrap();
    }
    let done = writer.finish().unwrap();
    assert_eq!(done.blocks.starts.len(), blocks);
    assert_eq!(done.blocks.text_len, (input.len() * blocks) as u64);
    done.out
}

#[test]
fn short_flushed_records_and_full_records_are_independent_gzip_members() {
    for len in [1, 64, 4096, BLOCK_TEXT_LEN] {
        let input = text(len, false);
        let raw = write_records(&input, 3);
        let mut at = 0;
        for _ in 0..3 {
            let len = u16::from_le_bytes([raw[at + 16], raw[at + 17]]) as usize + 1;
            decode_gzip_block(&raw[at..at + len], &input);
            at += len;
        }
        assert_eq!(&raw[at..], &EOF_BLOCK);
    }
}

fn samples(name: &str, mut run: impl FnMut() -> Vec<u8>) {
    let expected = run();
    let fingerprint = checksum(&expected);
    for sample in 0..5 {
        let start = std::time::Instant::now();
        let mut iterations = 0u64;
        while iterations < 100 || start.elapsed().as_millis() < 100 {
            std::hint::black_box(run());
            iterations += 1;
        }
        let ns = start.elapsed().as_nanos();
        println!(
            "name={name} sample={sample} iterations={iterations} total_ns={ns} ns_per_iter={:.3} bytes={} checksum={fingerprint:016x}",
            ns as f64 / iterations as f64,
            expected.len()
        );
    }
    assert_eq!(run(), expected);
}

#[test]
#[ignore = "paired timing only on the designated CPU worker"]
fn audit_compression_workload() {
    for len in [0, 1, 64, 4096, BLOCK_TEXT_LEN] {
        for random in [false, true] {
            let input = text(len, random);
            let mut compressor = libdeflater::Compressor::new(check_level(6).unwrap());
            samples(&format!("deflate_{len}_{random}"), || {
                let mut raw = Vec::new();
                deflate_block(&mut compressor, std::hint::black_box(&input), &mut raw).unwrap();
                raw
            });
        }
    }
    for len in [1, 64, 4096, BLOCK_TEXT_LEN] {
        let input = text(len, false);
        samples(&format!("writer_16x{len}"), || write_records(&input, 16));
    }
}
