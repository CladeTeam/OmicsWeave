//! cargo bench -p omicsweave-bbi --bench paths
//! Cached in-memory paths; preparation and one-shot costs are reported separately.

use omicsweave_bbi::{BedBlock, INDEX_MAGIC, Index, NodeInfo, WigBlock, ZoomBlock};
use std::hint::black_box;
use std::time::{Duration, Instant};

fn measure(name: &str, mut run: impl FnMut()) {
    if std::env::var("BBI_BENCH_FILTER").is_ok_and(|filter| !name.contains(&filter)) {
        return;
    }
    for _ in 0..100 {
        run();
    }
    let mut samples = [0.0f64; 5];
    for sample in &mut samples {
        let start = Instant::now();
        let mut count = 0;
        while start.elapsed() < Duration::from_millis(100) {
            for _ in 0..32 {
                run();
            }
            count += 32;
        }
        *sample = start.elapsed().as_nanos() as f64 / count as f64;
    }
    samples.sort_by(f64::total_cmp);
    println!("{name}\t{:.2}", samples[2]);
}

fn wig(kind: u8, n: u16, step: u32, span: u32) -> Vec<u8> {
    let mut buf = Vec::new();
    for v in [0, 100, 100 + (n as u32 - 1) * step + span, step, span] {
        buf.extend(v.to_le_bytes());
    }
    buf.extend([kind, 0]);
    buf.extend(n.to_le_bytes());
    for i in 0..n as u32 {
        if kind != 3 {
            buf.extend((100 + i * step).to_le_bytes());
        }
        if kind == 1 {
            buf.extend((100 + i * step + span).to_le_bytes());
        }
        buf.extend((i as f32).to_le_bytes());
    }
    buf
}

fn leaf_index(n: u16) -> Vec<u8> {
    let mut buf = Vec::new();
    buf.extend(INDEX_MAGIC.to_le_bytes());
    buf.extend((n as u32).to_le_bytes());
    buf.extend((n as u64).to_le_bytes());
    buf.extend([0u8; 32]);
    buf.extend([1, 0]);
    buf.extend(n.to_le_bytes());
    for i in 0..n as u32 {
        for v in [0, i * 100, 0, i * 100 + 100] {
            buf.extend(v.to_le_bytes());
        }
        buf.extend((i as u64 * 32).to_le_bytes());
        buf.extend(32u64.to_le_bytes());
    }
    buf
}

fn branch_index(leaves: u16, per_leaf: u16) -> Vec<u8> {
    let mut buf = leaf_index(0)[..48].to_vec();
    buf[8..16].copy_from_slice(&(leaves as u64 * per_leaf as u64).to_le_bytes());
    buf[28..32].copy_from_slice(&(leaves as u32 * per_leaf as u32 * 100).to_le_bytes());
    buf.extend([0, 0]);
    buf.extend(leaves.to_le_bytes());
    let first = 48 + 4 + leaves as u64 * 24;
    for leaf in 0..leaves as u32 {
        let lo = leaf * per_leaf as u32 * 100;
        let hi = lo + per_leaf as u32 * 100;
        for value in [0, lo, 0, hi] {
            buf.extend(value.to_le_bytes());
        }
        buf.extend((first + leaf as u64 * (4 + per_leaf as u64 * 32)).to_le_bytes());
    }
    for leaf in 0..leaves as u32 {
        buf.extend([1, 0]);
        buf.extend(per_leaf.to_le_bytes());
        for i in leaf * per_leaf as u32..(leaf + 1) * per_leaf as u32 {
            for value in [0, i * 100, 0, i * 100 + 100] {
                buf.extend(value.to_le_bytes());
            }
            buf.extend((i as u64 * 32).to_le_bytes());
            buf.extend(32u64.to_le_bytes());
        }
    }
    buf
}

fn main() {
    println!("case\tns/op (median of 5 x 100ms, warmed caller buffers)");
    let tree = branch_index(64, 64);
    let raw = Index::new(&tree, 0).unwrap();
    let mut nodes = Vec::new();
    let prepared = raw.prepare(&mut nodes).unwrap();
    for (mode, index) in [("raw", raw), ("prepared", prepared)] {
        for (range, start, end) in [
            ("narrow", 200_000, 200_010),
            ("full", 0, 409_600),
            ("miss", 500_000, 500_010),
        ] {
            let expected: Vec<_> = (0..4096u64)
                .filter(|&i| i * 100 < end as u64 && i * 100 + 100 > start as u64)
                .map(|i| i * 32)
                .collect();
            let got: Vec<_> = index
                .query(0, start, end)
                .map(|r| r.unwrap().offset)
                .collect();
            assert_eq!(got, expected);
            measure(&format!("index.tree64.{mode}.{range}"), || {
                black_box(
                    black_box(index)
                        .query(0, black_box(start), black_box(end))
                        .map(|r| r.unwrap().offset)
                        .sum::<u64>(),
                );
            });
        }
    }
    let mut reused = leaf_index(512)[..48].to_vec();
    reused.extend([0, 0, 128, 0]);
    let child = 48 + 4 + 128 * 24;
    for _ in 0..128 {
        for n in [0u32, 0, 0, 100_000] {
            reused.extend(n.to_le_bytes());
        }
        reused.extend((child as u64).to_le_bytes());
    }
    reused.extend_from_slice(&leaf_index(512)[48..]);
    let reused = Index::new(&reused, 0).unwrap();
    let mut nodes = Vec::new();
    measure("index.reused_leaf.prepare", || {
        assert!(black_box(reused).prepare(black_box(&mut nodes)).is_err());
    });
    for n in [1u16, 64, 4096] {
        for (name, kind) in [("bedgraph", 1), ("varstep", 2), ("fixedstep", 3)] {
            let buf = wig(kind, n, 4, 2);
            let raw = WigBlock::parse(&buf).unwrap();
            let ready = raw.prepare();
            let start = 100 + n as u32 / 2 * 4;
            let prefix = format!("{name}.{n}");
            measure(&format!("{prefix}.prepare"), || {
                black_box(black_box(raw).prepare());
            });
            for (mode, block) in [("raw", raw), ("prepared", ready)] {
                let mut out = [f32::NAN; 64];
                measure(&format!("{prefix}.{mode}.fill.narrow"), || {
                    black_box(black_box(block).fill(0, black_box(start), black_box(&mut out)));
                });
                measure(&format!("{prefix}.{mode}.overlap.narrow"), || {
                    black_box(
                        black_box(block)
                            .overlapping(0, black_box(start), black_box(start + 64))
                            .fold(0.0, |s, (_, _, v)| s + v),
                    );
                });
            }
            let mut out = [f32::NAN; 64];
            measure(&format!("{prefix}.parse.fill.narrow"), || {
                black_box(WigBlock::parse(black_box(&buf)).unwrap().fill(
                    0,
                    black_box(start),
                    black_box(&mut out),
                ));
            });
        }
        let buf = leaf_index(n);
        let raw = Index::new(&buf, 0).unwrap();
        let mut nodes: Vec<NodeInfo> = Vec::new();
        measure(&format!("index.{n}.prepare"), || {
            black_box(raw.prepare(black_box(&mut nodes)).unwrap());
        });
        let ready = raw.prepare(&mut nodes).unwrap();
        for (mode, index) in [("raw", raw), ("prepared", ready)] {
            for (range, start, end) in [
                ("narrow", n as u32 / 2 * 100, n as u32 / 2 * 100 + 10),
                ("full", 0, n as u32 * 100),
                ("miss", n as u32 * 100 + 100, n as u32 * 100 + 110),
            ] {
                let expected: Vec<_> = (0..n as u64)
                    .filter(|&i| i * 100 < end as u64 && i * 100 + 100 > start as u64)
                    .map(|i| i * 32)
                    .collect();
                let got: Vec<_> = index
                    .query(0, start, end)
                    .map(|r| r.unwrap().offset)
                    .collect();
                assert_eq!(got, expected);
                measure(&format!("index.{n}.{mode}.{range}"), || {
                    black_box(
                        black_box(index)
                            .query(0, black_box(start), black_box(end))
                            .map(|r| r.unwrap().offset)
                            .sum::<u64>(),
                    );
                });
            }
        }
    }
    for (name, step, span) in [("same-start", 0, 8192), ("wide-span", 1, 8192)] {
        let buf = wig(3, 4096, step, span);
        let block = WigBlock::parse(&buf).unwrap();
        let mut out = [0.0; 64];
        measure(&format!("fixedstep.overlap.{name}"), || {
            black_box(black_box(block).fill(0, black_box(2148), black_box(&mut out)));
        });
    }
    for kind in [1, 2] {
        let buf = wig(kind, 4096, 1, 8192);
        let raw = WigBlock::parse(&buf).unwrap();
        for (mode, block) in [("raw", raw), ("prepared", raw.prepare())] {
            let mut out = [0.0; 64];
            measure(&format!("wig.overlap.{kind}.{mode}"), || {
                black_box(black_box(block).fill(0, black_box(2148), black_box(&mut out)));
            });
        }
        let mut arbitrary = buf.clone();
        if kind == 1 {
            for i in 0..4096 {
                let at = 24 + i * 12;
                let end = 100 + i as u32 + if i % 2 == 0 { 8192 } else { 32 };
                arbitrary[at + 4..at + 8].copy_from_slice(&end.to_le_bytes());
            }
            let block = WigBlock::parse(&arbitrary).unwrap().prepare();
            let mut out = [0.0; 64];
            let mut scratch = vec![0; 3 * (out.len() + 1)];
            measure("wig.overlap.arbitrary.fill", || {
                black_box(black_box(block).fill(0, black_box(2148), black_box(&mut out)));
            });
            measure("wig.overlap.arbitrary.scratch", || {
                black_box(
                    black_box(block)
                        .fill_with_scratch(
                            0,
                            black_box(2148),
                            black_box(&mut out),
                            black_box(&mut scratch),
                        )
                        .unwrap(),
                );
            });
        }
    }
    let mut zoom = Vec::new();
    for i in 0..4096u32 {
        for v in [0, i * 4, i * 4 + 4, 4, 1, 2, 3, 4] {
            zoom.extend(v.to_le_bytes());
        }
    }
    let raw = ZoomBlock::parse(&zoom).unwrap();
    measure("zoom.prepare", || {
        black_box(black_box(raw).prepare());
    });
    for (mode, block) in [("raw", raw), ("prepared", raw.prepare())] {
        measure(&format!("zoom.{mode}.narrow"), || {
            black_box(
                black_box(block)
                    .overlapping(0, black_box(8000), black_box(8064))
                    .map(|r| r.sum)
                    .sum::<f32>(),
            );
        });
    }
    for len in [8, 256] {
        let mut buf = Vec::new();
        for i in 0..4096u32 {
            for v in [0, i * 4, i * 4 + 2] {
                buf.extend(v.to_le_bytes());
            }
            buf.extend(std::iter::repeat_n(b'a', len));
            buf.push(0);
        }
        let raw = BedBlock::parse(&buf).unwrap();
        let mut offsets: Vec<u32> = Vec::new();
        measure(&format!("bed.{len}.prepare"), || {
            black_box(raw.prepare(black_box(&mut offsets)).unwrap());
        });
        let ready = raw.prepare(&mut offsets).unwrap();
        measure(&format!("bed.{len}.raw.narrow"), || {
            black_box(
                raw.entries()
                    .map(Result::unwrap)
                    .filter(|r| r.overlaps(0, black_box(8000), black_box(8064)))
                    .map(|r| r.rest.len())
                    .sum::<usize>(),
            );
        });
        measure(&format!("bed.{len}.prepared.narrow"), || {
            black_box(
                black_box(ready)
                    .overlapping(0, black_box(8000), black_box(8064))
                    .map(|r| r.rest.len())
                    .sum::<usize>(),
            );
        });
        measure(&format!("bed.{len}.raw.nth"), || {
            black_box(raw.entries().nth(black_box(2048)).unwrap().unwrap());
        });
        measure(&format!("bed.{len}.prepared.nth"), || {
            black_box(black_box(ready).entries().nth(black_box(2048)).unwrap());
        });
    }
}
