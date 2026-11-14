//! Portable paired benchmark: copy this exact file into both source trees.
use std::hint::black_box;
use std::time::{Duration, Instant};

use omicsweave_fasta::{Parsed, Record};

fn fixture(kind: &str) -> Vec<u8> {
    let len = 65536u32;
    let blocks: Vec<(u32, u32)> = match kind {
        "ordinary" => (0..8192).map(|i| (i * 4, 1)).collect(),
        "adjacent" => (0..len).map(|i| (i, 1)).collect(),
        "empty" => std::iter::once((0, 1))
            .chain(std::iter::repeat_n((1, 0), 65530))
            .chain(std::iter::once((1, 1)))
            .collect(),
        _ => unreachable!(),
    };
    let mut bytes = len.to_le_bytes().to_vec();
    for list in [&blocks[..], &blocks[..]] {
        bytes.extend_from_slice(&(list.len() as u32).to_le_bytes());
        for &(start, _) in list {
            bytes.extend_from_slice(&start.to_le_bytes());
        }
        for &(_, len) in list {
            bytes.extend_from_slice(&len.to_le_bytes());
        }
    }
    bytes.extend_from_slice(&0u32.to_le_bytes());
    bytes.resize(bytes.len() + len as usize / 4, 0x1B);
    bytes
}

fn parse(bytes: &[u8]) -> Record {
    let Parsed::Done(record) = Record::parse(bytes, 0).unwrap() else {
        panic!("complete fixture");
    };
    record
}

fn measure(name: &str, mut f: impl FnMut()) {
    for _ in 0..20 {
        f();
    }
    let mut samples = [0f64; 5];
    for sample in &mut samples {
        let start = Instant::now();
        let mut n = 0;
        while start.elapsed() < Duration::from_millis(100) {
            for _ in 0..16 {
                f();
            }
            n += 16;
        }
        *sample = start.elapsed().as_nanos() as f64 / n as f64;
    }
    samples.sort_by(f64::total_cmp);
    println!("{name}\t{:.2}", samples[2]);
}

fn main() {
    for kind in ["ordinary", "adjacent", "empty"] {
        let bytes = fixture(kind);
        let record = parse(&bytes);
        println!(
            "fixture\t{kind}\tbytes={}\tblocks={}",
            bytes.len(),
            record.n_blocks().len()
        );
        measure(&format!("runs.{kind}.parse"), || {
            black_box(parse(black_box(&bytes)));
        });
        for (start, end) in [(0, 2), (0, 3), (123, 187), (0, 65536)] {
            let range = record.packed_range(start, end).unwrap();
            let packed = &bytes[range.start as usize..range.end as usize];
            let mut out = vec![0; (end - start) as usize];
            record.decode(start, packed, &mut out, true).unwrap();
            let wanted: Vec<u8> = (start..end)
                .map(|p| match kind {
                    "adjacent" => b'N',
                    "empty" if p < 2 => b'N',
                    "ordinary" if p < 32768 && p % 4 == 0 => b'N',
                    _ => b"TCAG"[(p % 4) as usize],
                })
                .collect();
            assert_eq!(out, wanted);
            println!(
                "digest\t{kind}.{start}-{end}\t{:016x}",
                out.iter()
                    .fold(0u64, |h, &b| h.wrapping_mul(131).wrapping_add(b as u64))
            );
            measure(&format!("runs.{kind}.plan.{start}-{end}"), || {
                black_box(
                    record
                        .packed_range(black_box(start), black_box(end))
                        .unwrap(),
                );
            });
            measure(&format!("runs.{kind}.decode.{start}-{end}"), || {
                record
                    .decode(
                        black_box(start),
                        black_box(packed),
                        black_box(&mut out),
                        true,
                    )
                    .unwrap();
                black_box(&out);
            });
            assert_eq!(out, wanted);
        }
    }
}
