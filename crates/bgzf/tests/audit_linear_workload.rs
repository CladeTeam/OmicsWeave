use omicsweave_bgzf::{BinningIndex, Blocks, IndexFormat, Logical};

fn crc32(text: &[u8]) -> u32 {
    let mut crc = !0u32;
    for &b in text {
        crc ^= b as u32;
        for _ in 0..8 {
            crc = (crc >> 1) ^ (0xedb8_8320 & 0u32.wrapping_sub(crc & 1));
        }
    }
    !crc
}

#[test]
#[ignore = "manual timing workload: run only on the designated CPU pod"]
fn audit_linear_overlap_workload() {
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
