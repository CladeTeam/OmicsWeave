//! Writing on a worker pool and streaming maps block by block: for any
//! executor, and whether a map is given whole or block by block (cells of a
//! block in any order), the file and sidecar are byte for byte the same.

use crate::support as common;

use std::collections::{BTreeMap, HashSet};
use std::io::Cursor;
use std::sync::LazyLock;

use common::Rng;
use omicsweave_core::ContactUnit;
use omicsweave_core::exec::{Exec, Pool};
use omicsweave_hic::{
    BlockBins, BlockLayout, BuildOptions, Chrom, DEPTH_BASE_ATTRIBUTE, HicBuilder, HicHeader,
    HicReader, HicSidecar, HicWriter, MatrixLayouts, Record,
};
use proptest::prelude::*;

static POOL4: LazyLock<Pool> = LazyLock::new(|| Pool::new(4));
static POOL8: LazyLock<Pool> = LazyLock::new(|| Pool::new(8));

/// Inline, then pools with a single block in flight, a few, and many.
fn executors() -> [Exec<'static>; 5] {
    [
        Exec::Inline,
        Exec::Pool {
            pool: &POOL4,
            limit: 1,
        },
        Exec::Pool {
            pool: &POOL4,
            limit: 3,
        },
        Exec::Pool {
            pool: &POOL8,
            limit: 2,
        },
        Exec::Pool {
            pool: &POOL8,
            limit: 10_000,
        },
    ]
}

fn random_header(rng: &mut Rng) -> HicHeader {
    let chroms = (0..rng.range(1, 3))
        .map(|i| Chrom {
            name: format!("c{i}"),
            len: rng.range(1, 6_000),
        })
        .collect();
    let mut resolutions: Vec<u32> = [1, 3, 10, 50, 200, 1_000]
        .into_iter()
        .filter(|_| rng.chance(50))
        .collect();
    if resolutions.is_empty() {
        resolutions.push(25);
    }
    let mut header = HicHeader::new("g", chroms, resolutions);
    header.version = if rng.chance(30) { 8 } else { 9 };
    match rng.range(0, 3) {
        0 => {}
        1 => header
            .attributes
            .push((DEPTH_BASE_ATTRIBUTE.into(), "3".into())),
        2 => header
            .attributes
            .push((DEPTH_BASE_ATTRIBUTE.into(), "-2".into())),
        _ => header
            .attributes
            .push((DEPTH_BASE_ATTRIBUTE.into(), "0".into())),
    }
    if rng.chance(40) {
        header.frag_resolutions = vec![1, 4];
        header.sites = header
            .chroms
            .iter()
            .map(|c| {
                (1..c.len as u32)
                    .step_by(rng.range(5, 90) as usize)
                    .collect()
            })
            .collect();
    }
    header
}

/// A map as a caller gives it: chromosomes in either order, records in
/// either orientation, cells repeated, zeros included.
struct Map {
    chr1: u32,
    chr2: u32,
    unit: ContactUnit,
    bin_size: u32,
    blocks: Option<BlockBins>,
    records: Vec<Record>,
}

fn random_value(rng: &mut Rng) -> f32 {
    match rng.range(0, 9) {
        0 => 0.0,
        1 => rng.range(0, 1_000_000) as f32 / 64.0,
        2 => 40_000.0,
        _ => rng.range(1, 6) as f32,
    }
}

fn random_maps(rng: &mut Rng, header: &HicHeader) -> Vec<Map> {
    let n = header.chroms.len() as u32;
    let mut seen = HashSet::new();
    let mut maps = Vec::new();
    for unit in [ContactUnit::Bp, ContactUnit::Frag] {
        for &bin_size in header.resolutions(unit) {
            for _ in 0..rng.range(1, 4) {
                let (chr1, chr2) = (
                    rng.range(0, n as u64 - 1) as u32,
                    rng.range(0, n as u64 - 1) as u32,
                );
                if !seen.insert((chr1.min(chr2), chr1.max(chr2), unit, bin_size)) {
                    continue;
                }
                let (n1, n2) = (
                    header.bins(chr1, unit, bin_size) as u32,
                    header.bins(chr2, unit, bin_size) as u32,
                );
                let blocks = match rng.range(0, 3) {
                    0 => Some(BlockBins::Juicer),
                    1 => Some(BlockBins::Fixed(rng.range(1, 40) as u32)),
                    _ => None,
                };
                let mut records: Vec<Record> = (0..rng.range(0, 3_000))
                    .map(|_| Record {
                        x: rng.range(0, n1 as u64 - 1) as u32,
                        y: rng.range(0, n2 as u64 - 1) as u32,
                        value: random_value(rng),
                    })
                    .collect();
                for i in 0..records.len() / 8 {
                    let mut again = records[i * 7 % records.len()];
                    again.value = random_value(rng);
                    if chr1 == chr2 && rng.chance(50) {
                        std::mem::swap(&mut again.x, &mut again.y);
                    }
                    records.push(again);
                }
                maps.push(Map {
                    chr1,
                    chr2,
                    unit,
                    bin_size,
                    blocks,
                    records,
                });
            }
        }
    }
    maps
}

type Written = (Vec<u8>, HicSidecar);

/// Writes every map whole.
fn whole(header: &HicHeader, maps: &[Map], exec: Exec<'_>) -> Written {
    let out = Cursor::new(Vec::new());
    let mut w =
        HicWriter::create_exec(out, header.clone(), 1, BlockLayout::default(), exec).unwrap();
    for m in maps {
        match m.blocks {
            Some(b) => {
                w.write_matrix_with_blocks(m.chr1, m.chr2, m.unit, m.bin_size, &m.records, b)
            }
            None => w.write_matrix(m.chr1, m.chr2, m.unit, m.bin_size, &m.records),
        }
        .unwrap();
    }
    let (out, sidecar) = w.finish_with_sidecar().unwrap();
    (out.into_inner(), sidecar)
}

/// Streams every map block by block, shuffling each block's cells (records
/// of one cell keep their order, which decides how their floats add up).
fn streamed(header: &HicHeader, maps: &[Map], exec: Exec<'_>, seed: u64) -> Written {
    let out = Cursor::new(Vec::new());
    let mut w =
        HicWriter::create_exec(out, header.clone(), 1, BlockLayout::default(), exec).unwrap();
    let mut rng = Rng(seed);
    for m in maps {
        let layout = match m.blocks {
            Some(b) => w.matrix_layout_with_blocks(m.chr1, m.chr2, m.unit, m.bin_size, b),
            None => w.matrix_layout(m.chr1, m.chr2, m.unit, m.bin_size),
        }
        .unwrap();
        let mut blocks: BTreeMap<u32, Vec<Record>> = BTreeMap::new();
        for r in &m.records {
            let (number, x, y) = layout.place(r.x, r.y);
            blocks.entry(number).or_default().push(Record {
                x,
                y,
                value: r.value,
            });
        }
        w.begin_matrix(&layout).unwrap();
        for (number, mut records) in blocks {
            let salt = rng.next();
            records.sort_by_key(|r| {
                (((r.x as u64) << 32 | r.y as u64) ^ salt).wrapping_mul(0x9E37_79B9_7F4A_7C15)
            });
            w.add_block(number, records).unwrap();
        }
        w.end_matrix().unwrap();
    }
    let (out, sidecar) = w.finish_with_sidecar().unwrap();
    (out.into_inner(), sidecar)
}

fn same_everywhere(seed: u64) {
    let mut rng = Rng(seed);
    let header = random_header(&mut rng);
    let maps = random_maps(&mut rng, &header);
    let want = whole(&header, &maps, Exec::Inline);
    assert_eq!(want.1, HicSidecar::build(want.0.as_slice()).unwrap());
    for exec in executors() {
        let got = whole(&header, &maps, exec);
        assert!(got.0 == want.0, "whole maps differ");
        assert_eq!(got.1, want.1);
        let got = streamed(&header, &maps, exec, seed);
        assert!(got.0 == want.0, "streamed maps differ");
        assert_eq!(got.1, want.1);
    }
}

proptest! {
    #![proptest_config(ProptestConfig { cases: 48, failure_persistence: Some(Box::new(proptest::test_runner::FileFailurePersistence::Direct(concat!(env!("CARGO_MANIFEST_DIR"), "/tests/cases/execution.proptest-regressions")))), ..ProptestConfig::default() })]

    #[test]
    fn maps_are_the_same_on_every_executor_and_when_streamed(seed in any::<u64>()) {
        same_everywhere(seed);
    }

    #[test]
    fn builder_is_the_same_on_every_executor(seed in any::<u64>()) {
        let mut rng = Rng(seed);
        let header = random_header(&mut rng);
        let n = header.chroms.len() as u64;
        let contacts: Vec<_> = (0..rng.range(0, 4_000))
            .map(|_| {
                let (a, b) = (rng.range(0, n - 1) as u32, rng.range(0, n - 1) as u32);
                let pa = rng.range(0, header.chroms[a as usize].len + 5) as u32;
                let pb = rng.range(0, header.chroms[b as usize].len + 5) as u32;
                (a, pa, b, pb, random_value(&mut rng))
            })
            .collect();
        let build = |exec: Exec<'_>| {
            let mut builder = HicBuilder::new(header.clone()).unwrap();
            for &(a, pa, b, pb, v) in &contacts {
                builder.add(a, pa, b, pb, v).unwrap();
            }
            let options = BuildOptions { level: 1, blocks: BlockLayout::default() };
            let (out, sidecar) = builder.finish_exec(Cursor::new(Vec::new()), &options, exec).unwrap();
            (out.into_inner(), sidecar)
        };
        let want = build(Exec::Inline);
        for exec in executors() {
            let got = build(exec);
            prop_assert!(got.0 == want.0, "file differs");
            prop_assert_eq!(&got.1, &want.1);
        }
    }
}

/// The (sum, occupied cells, 5th, 95th percentile) stored for each
/// resolution of a chromosome pair's matrix record.
fn stored_summaries(bytes: &[u8], c1: u32, c2: u32) -> Vec<[f32; 4]> {
    let reader = HicReader::open(bytes).unwrap();
    let mut at = reader.footer().matrices[&(c1, c2)] as usize;
    let i32_at = |at: &mut usize| {
        let v = i32::from_le_bytes(bytes[*at..*at + 4].try_into().unwrap());
        *at += 4;
        v
    };
    at += 8;
    let n = i32_at(&mut at);
    let mut out = Vec::new();
    for _ in 0..n {
        at += bytes[at..].iter().position(|&b| b == 0).unwrap() + 1;
        at += 4;
        let f: Vec<f32> = (0..4)
            .map(|i| f32::from_le_bytes(bytes[at + 4 * i..at + 4 * i + 4].try_into().unwrap()))
            .collect();
        out.push([f[0], f[1], f[2], f[3]]);
        at += 16 + 12;
        let blocks = i32_at(&mut at);
        at += blocks as usize * 16;
    }
    out
}

/// Apache Commons' default percentile over all values at once.
fn percentile(sorted: &[f32], p: f64) -> f32 {
    let n = sorted.len();
    let pos = p * (n + 1) as f64 / 100.0;
    if pos < 1.0 {
        return sorted[0];
    }
    if pos >= n as f64 {
        return sorted[n - 1];
    }
    let at = pos.floor() as usize;
    let (lower, upper) = (sorted[at - 1] as f64, sorted[at] as f64);
    (lower + (pos - pos.floor()) * (upper - lower)) as f32
}

/// Enough cells that the writer folds values into runs several times, with
/// heavily repeated small counts and a tail of distinct ones; the stored
/// statistics match those computed over every value at once.
#[test]
fn map_statistics_match_all_values_at_once() {
    let header = HicHeader::new(
        "g",
        vec![Chrom {
            name: "s".into(),
            len: 700,
        }],
        vec![1],
    );
    let mut rng = Rng(7);
    let mut cells = BTreeMap::new();
    let mut records = Vec::new();
    for y in 0..=700u32 {
        for x in 0..=y {
            let value = match rng.range(0, 99) {
                0..=59 => 1.0,
                60..=89 => rng.range(2, 9) as f32,
                90..=97 => rng.range(10, 30_000) as f32,
                _ => rng.range(1, 1 << 30) as f32 / 1024.0,
            };
            cells.insert((x, y), value);
            records.push(Record { x, y, value });
        }
    }
    assert!(cells.len() > 3 << 16);
    let mut values: Vec<f32> = cells.values().copied().collect();
    values.sort_by(f32::total_cmp);
    let sum: f64 = cells
        .iter()
        .map(|(&(x, y), &v)| if x == y { v as f64 } else { 2.0 * v as f64 })
        .sum();
    for exec in [executors()[0], executors()[4]] {
        let out = Cursor::new(Vec::new());
        let mut w =
            HicWriter::create_exec(out, header.clone(), 1, BlockLayout::default(), exec).unwrap();
        w.write_matrix(0, 0, ContactUnit::Bp, 1, &records).unwrap();
        let bytes = w.finish().unwrap().into_inner();
        let [s, n, p5, p95] = stored_summaries(&bytes, 0, 0)[0];
        assert_eq!(n, values.len() as f32);
        assert_eq!(p5.to_bits(), percentile(&values, 5.0).to_bits());
        assert_eq!(p95.to_bits(), percentile(&values, 95.0).to_bits());
        assert!(((s as f64 - sum) / sum).abs() < 1e-6);
    }
}

fn shareable<T: Clone + Send + Sync + 'static>(_: &T) {}

/// Layouts made away from the writer, on other threads, are the writer's:
/// the same layout or the same error for every map, and a map streamed
/// with them writes the file `write_matrix` does.
#[test]
fn shared_layouts_are_the_writers() {
    for seed in 0..40u64 {
        let mut rng = Rng(seed);
        let header = random_header(&mut rng);
        let blocks = [BlockLayout::default(), BlockLayout::juicer()][seed as usize % 2];
        let out = Cursor::new(Vec::new());
        let w = HicWriter::create_with(out, header.clone(), 1, blocks).unwrap();
        let layouts = MatrixLayouts::new(&header, blocks).unwrap();
        shareable(&layouts);
        let n = header.chroms.len() as u32;
        let mut keys = Vec::new();
        for unit in [ContactUnit::Bp, ContactUnit::Frag] {
            for res in [1, 3, 4, 25, 1_000, 7] {
                for c1 in 0..=n {
                    for c2 in 0..=n {
                        for b in [None, Some(BlockBins::Juicer), Some(BlockBins::Fixed(3))] {
                            keys.push((c1, c2, unit, res, b));
                        }
                    }
                }
            }
        }
        let theirs: Vec<_> = std::thread::scope(|scope| {
            let shared = w.layouts();
            let from_new = layouts.clone();
            let keys = &keys;
            let one = scope.spawn(move || {
                keys.iter()
                    .map(|&(c1, c2, unit, res, b)| match b {
                        Some(b) => shared.layout_with_blocks(c1, c2, unit, res, b),
                        None => shared.layout(c1, c2, unit, res),
                    })
                    .collect::<Vec<_>>()
            });
            let two = scope.spawn(move || {
                keys.iter()
                    .map(|&(c1, c2, unit, res, b)| match b {
                        Some(b) => from_new.layout_with_blocks(c1, c2, unit, res, b),
                        None => from_new.layout(c1, c2, unit, res),
                    })
                    .collect::<Vec<_>>()
            });
            one.join()
                .unwrap()
                .into_iter()
                .zip(two.join().unwrap())
                .collect()
        });
        for (&(c1, c2, unit, res, b), (one, two)) in keys.iter().zip(theirs) {
            let own = match b {
                Some(b) => w.matrix_layout_with_blocks(c1, c2, unit, res, b),
                None => w.matrix_layout(c1, c2, unit, res),
            };
            match (own, one, two) {
                (Ok(own), Ok(one), Ok(two)) => {
                    assert_eq!(own, one);
                    assert_eq!(own, two);
                }
                (Err(_), Err(_), Err(_)) => {}
                other => panic!("{c1} {c2} {unit:?} {res} {b:?}: {other:?}"),
            }
        }

        // Cells placed on worker threads, streamed by the writer.
        let maps = random_maps(&mut rng, &header);
        let want = {
            let out = Cursor::new(Vec::new());
            let mut w = HicWriter::create_with(out, header.clone(), 1, blocks).unwrap();
            for m in &maps {
                w.write_matrix(m.chr1, m.chr2, m.unit, m.bin_size, &m.records)
                    .unwrap();
            }
            w.finish_with_sidecar().unwrap()
        };
        let placed: Vec<_> = std::thread::scope(|scope| {
            let handles: Vec<_> = maps
                .iter()
                .map(|m| {
                    let layouts = layouts.clone();
                    scope.spawn(move || {
                        let layout = layouts.layout(m.chr1, m.chr2, m.unit, m.bin_size).unwrap();
                        let mut blocks: BTreeMap<u32, Vec<Record>> = BTreeMap::new();
                        for r in &m.records {
                            let (number, x, y) = layout.place(r.x, r.y);
                            let value = r.value;
                            blocks
                                .entry(number)
                                .or_default()
                                .push(Record { x, y, value });
                        }
                        (layout, blocks)
                    })
                })
                .collect();
            handles.into_iter().map(|h| h.join().unwrap()).collect()
        });
        let out = Cursor::new(Vec::new());
        let mut w = HicWriter::create_with(out, header.clone(), 1, blocks).unwrap();
        for (layout, blocks) in placed {
            w.begin_matrix(&layout).unwrap();
            for (number, records) in blocks {
                w.add_block(number, records).unwrap();
            }
            w.end_matrix().unwrap();
        }
        let got = w.finish_with_sidecar().unwrap();
        assert!(got.0.into_inner() == want.0.into_inner());
        assert_eq!(got.1, want.1);
    }
}
