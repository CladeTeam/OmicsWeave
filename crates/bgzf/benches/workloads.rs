//! Bounded CPU workloads; run with `cargo bench -p omicsweave-bgzf --bench workloads`.
#[path = "../tests/support/wire.rs"]
mod wire;
use omicsweave_bgzf::{
    BLOCK_TEXT_LEN, BinningIndex, Blocks, IndexFormat, Logical, Writer, check_level, deflate_block,
    inflate_block_owned,
};
use omicsweave_core::exec::Inline;
use wire::{crc32, full_zero_block, stored};

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

fn compression() {
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

fn inflation() {
    for len in [0, 1, 4096, 65_536] {
        let text = vec![0; len];
        let raw = if len == 65_536 {
            full_zero_block()
        } else {
            stored(&text)
        };
        assert_eq!(inflate_block_owned(&raw).unwrap(), text);
        for sample in 0..5 {
            let mut rounds = 0u64;
            let begin = std::time::Instant::now();
            while rounds < 100 || begin.elapsed().as_millis() < 100 {
                let decoded = inflate_block_owned(std::hint::black_box(&raw)).unwrap();
                std::hint::black_box(decoded);
                rounds += 1;
            }
            let ns = begin.elapsed().as_nanos();
            println!(
                "name=owned_inflate_{len} sample={sample} iterations={rounds} total_ns={ns} ns_per_iter={:.3}",
                ns as f64 / rounds as f64
            );
        }
    }
}

fn indexing() {
    let blocks = Blocks {
        starts: vec![(0, 0)],
        end: 100,
        text_len: 8192,
    };
    for (name, width) in [("narrow", 1), ("overlap", 1i64 << 29)] {
        let mut push = std::time::Duration::ZERO;
        let mut rounds = 0u128;
        let mut bytes = Vec::new();
        while rounds < 5 || push < std::time::Duration::from_millis(100) {
            let mut index = BinningIndex::new(IndexFormat::Bai, 1, Logical::new(0, 0));
            let start = std::time::Instant::now();
            for i in 1..=8192 {
                index.push(0, 0, width, Logical::new(0, i), true).unwrap();
            }
            push += start.elapsed();
            bytes = index.finish(Logical::new(0, 8192), &blocks).unwrap();
            rounds += 1;
        }
        println!(
            "{name}: records_per_round=8192 windows={} rounds={} total_push_ns={} push_ns_per_round={} bytes={} crc32={:08x}",
            ((width - 1) >> 14) + 1,
            rounds,
            push.as_nanos(),
            push.as_nanos() / rounds,
            bytes.len(),
            crc32(&bytes),
        );
        std::hint::black_box(bytes);
    }
}

fn main() {
    compression();
    inflation();
    indexing();
}
