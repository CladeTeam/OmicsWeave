//! Differential tests against libBigWig: random bigWig files written by
//! libBigWig, read by both libBigWig and the crate.

mod common;

use std::path::Path;

use common::libbigwig::{Reader, Writer};
use common::{Opened, Rng, same};
use omicsweave_bbi::{Kind, Level};
use proptest::prelude::*;

/// Writes a random bigWig mixing all three section types and returns its
/// chromosomes.
fn write_random(path: &Path, rng: &mut Rng) -> Vec<(String, u32)> {
    let n = rng.range(1, 4) as usize;
    let chroms: Vec<(String, u32)> = (0..n)
        .map(|i| {
            let name = if rng.chance(30) {
                format!("scaffold_{i}_{}", rng.next() % 100_000)
            } else {
                format!("chr{}", i + 1)
            };
            (name, rng.range(1_000, 300_000) as u32)
        })
        .collect();
    let refs: Vec<(&str, u32)> = chroms.iter().map(|(n, l)| (n.as_str(), *l)).collect();
    let max_zooms = rng.range(0, 10) as i32;
    let mut w = Writer::create(path, max_zooms, &refs).expect("create");
    let mut written = false;

    for (name, len) in &chroms {
        let len = *len;
        let mut pos = rng.range(0, 2_000) as u32;
        for _ in 0..rng.range(1, 60) {
            if pos >= len {
                break;
            }
            let count = if rng.chance(5) {
                rng.range(1, 3_000)
            } else {
                rng.range(1, 300)
            };
            match rng.range(0, 2) {
                0 => {
                    let mut items = Vec::new();
                    for _ in 0..count {
                        if pos >= len {
                            break;
                        }
                        let end = (pos + rng.range(1, 50) as u32).min(len);
                        items.push((pos, end, rng.value()));
                        pos = end + rng.range(0, 30) as u32;
                    }
                    assert!(w.add_bedgraph(name, &items), "bedGraph");
                    written = true;
                }
                1 => {
                    let span = rng.range(1, 20) as u32;
                    let mut items = Vec::new();
                    for _ in 0..count {
                        if pos + span > len {
                            break;
                        }
                        items.push((pos, rng.value()));
                        pos += span + rng.range(0, 30) as u32;
                    }
                    if !items.is_empty() {
                        assert!(w.add_varstep(name, span, &items), "varStep");
                        written = true;
                    }
                }
                _ => {
                    let span = rng.range(1, 20) as u32;
                    let step = span + rng.range(0, 10) as u32;
                    let fit = if pos + span > len {
                        0
                    } else {
                        (len - pos - span) / step + 1
                    };
                    let k = (count as u32).min(fit);
                    if k > 0 {
                        let values: Vec<f32> = (0..k).map(|_| rng.value()).collect();
                        assert!(w.add_fixedstep(name, pos, span, step, &values), "fixedStep");
                        written = true;
                        pos += k * step + rng.range(0, 30) as u32;
                    }
                }
            }
        }
    }
    if !written {
        // Interval and zoom comparisons require a file with a data index.
        assert!(w.add_bedgraph(&chroms[0].0, &[(0, 1, 0.0)]));
    }
    drop(w);
    chroms
}

fn check_metadata(ours: &Opened, theirs: &Reader) {
    let h = theirs.header();
    assert_eq!(h.kind, 0);
    assert_eq!(ours.header.kind, Kind::BigWig);
    assert_eq!(ours.header.version, h.version);
    assert_eq!(ours.header.zoom_levels, h.zoom_levels);
    assert_eq!(ours.header.chrom_tree_offset, h.chrom_tree);
    assert_eq!(ours.header.full_data_offset, h.full_data);
    assert_eq!(ours.header.full_index_offset, h.full_index);
    assert_eq!(ours.header.auto_sql_offset, h.auto_sql);
    assert_eq!(ours.header.total_summary_offset, h.total_summary);
    assert_eq!(ours.header.uncompress_buf_size, h.buf_size);
    assert_eq!(ours.header.extension_offset, h.extension);

    let s = ours.summary.expect("total summary");
    assert_eq!(s.bases_covered, h.bases_covered);
    assert_eq!(
        (s.min, s.max, s.sum, s.sum_squares),
        (h.min, h.max, h.sum, h.sum_squares)
    );

    let zooms: Vec<(u32, u64, u64)> = ours
        .zooms
        .iter()
        .map(|z| (z.reduction, z.data_offset, z.index_offset))
        .collect();
    assert_eq!(zooms, theirs.zooms());

    let chroms: Vec<(String, u32)> = ours
        .chroms
        .iter()
        .map(|(_, c)| (c.name.clone(), c.len))
        .collect();
    assert_eq!(chroms, theirs.chroms());
}

fn check_query(ours: &mut Opened, theirs: &mut Reader, chrom: &str, start: u32, end: u32) {
    let a = ours.values(chrom, start, end);
    let b = theirs.values(chrom, start, end).expect("libBigWig values");
    assert_eq!(a.len(), b.len());
    if let Some(i) = (0..a.len()).find(|&i| !same(a[i], b[i])) {
        panic!(
            "{chrom}:{start}-{end}: base {} ours {} libBigWig {}",
            start as usize + i,
            a[i],
            b[i]
        );
    }
    let a = ours.intervals(chrom, start, end);
    let b = theirs
        .intervals(chrom, start, end)
        .expect("libBigWig intervals");
    assert_eq!(a, b, "{chrom}:{start}-{end} intervals");
}

/// Zoom records read by the crate must describe the base data libBigWig
/// reads: each record's covered bases, min and max equal those of the base
/// intervals inside its window, and the counts add up to the file total.
///
/// Record sums are not compared: libBigWig's writer resets the running sum
/// when a zoom buffer fills (`addIntervalValue`) before storing it, so the
/// sums it writes are not exact.
fn check_zooms(ours: &mut Opened, theirs: &mut Reader) {
    let h = theirs.header();
    let chroms: Vec<(u32, String, u32)> = ours
        .chroms
        .iter()
        .map(|(id, c)| (id, c.name.clone(), c.len))
        .collect();
    let base: Vec<Vec<(u32, u32, f32)>> = chroms
        .iter()
        .map(|(_, name, len)| theirs.intervals(name, 0, *len).expect("intervals"))
        .collect();
    for zoom in 0..ours.zooms.len() {
        let mut valid = 0u64;
        for ((id, _, len), items) in chroms.iter().zip(&base) {
            let records = ours.summaries(zoom, *id, 0, *len);
            for pair in records.windows(2) {
                assert!(
                    pair[0].end <= pair[1].start,
                    "zoom {zoom}: overlapping records {pair:?}"
                );
            }
            for r in &records {
                let first = items.partition_point(|&(_, e, _)| e <= r.start);
                let (mut covered, mut min, mut max) = (0u32, f32::INFINITY, f32::NEG_INFINITY);
                for &(s, e, v) in items[first..].iter().take_while(|&&(s, _, _)| s < r.end) {
                    covered += e.min(r.end) - s.max(r.start);
                    min = min.min(v);
                    max = max.max(v);
                }
                assert_eq!(r.valid, covered, "zoom {zoom}: covered bases of {r:?}");
                if covered > 0 {
                    assert_eq!((r.min, r.max), (min, max), "zoom {zoom}: range of {r:?}");
                }
                valid += r.valid as u64;
            }
        }
        assert_eq!(valid, h.bases_covered, "zoom {zoom} covered bases");
    }
}

proptest! {
    #![proptest_config(ProptestConfig { cases: 24, ..ProptestConfig::default() })]

    #[test]
    fn matches_libbigwig(seed in any::<u64>()) {
        let dir = tempfile::tempdir().unwrap();
        let path = dir.path().join("t.bw");
        let mut rng = Rng::new(seed);
        let chroms = write_random(&path, &mut rng);

        let mut ours = Opened::open(&path);
        let mut theirs = Reader::open(&path).expect("libBigWig open");
        check_metadata(&ours, &theirs);

        for (name, len) in &chroms {
            check_query(&mut ours, &mut theirs, name, 0, *len);
        }
        for _ in 0..40 {
            let (name, len) = &chroms[rng.range(0, chroms.len() as u64 - 1) as usize];
            let start = rng.range(0, *len as u64 + 50) as u32;
            let end = start + rng.range(1, 20_000) as u32;
            check_query(&mut ours, &mut theirs, name, start, end);
        }
        check_zooms(&mut ours, &mut theirs);
    }
}

#[test]
fn layout_bounds_every_index() {
    let dir = tempfile::tempdir().unwrap();
    let path = dir.path().join("t.bw");
    write_random(&path, &mut Rng::new(7));
    let ours = Opened::open(&path);
    let starts = ours.layout.section_starts();
    assert_eq!(*starts.last().unwrap(), ours.layout.file_len());
    let base = ours.layout.index(Level::Base).unwrap();
    assert_eq!(base.start, ours.header.full_index_offset);
    assert!(base.end <= ours.layout.file_len());
    for (i, z) in ours.zooms.iter().enumerate() {
        assert_eq!(
            ours.layout.index(Level::Zoom(i)).unwrap().start,
            z.index_offset
        );
    }
    assert!(ours.layout.index(Level::Zoom(ours.zooms.len())).is_none());
}

#[test]
fn standard_reader_reuses_prepared_index_across_windows() {
    use omicsweave_core::{Region, TrackEntry, TrackRead};
    for seed in [7, 23] {
        let dir = tempfile::tempdir().unwrap();
        let path = dir.path().join("reader.bw");
        let chroms = write_random(&path, &mut Rng::new(seed));
        let mut reader = omicsweave_bbi::BbiReader::open(std::fs::read(&path).unwrap()).unwrap();
        let mut reference = Reader::open(&path).unwrap();
        for _ in 0..2 {
            for (name, len) in &chroms {
                for start in [0, len / 2, len.saturating_sub(64), len + 3] {
                    let region = Region::new(name, start, start + 64);
                    let mut out = [f32::NAN; 64];
                    reader.values(region, &mut out).unwrap();
                    let expected = reference.values(name, start, start + 64).unwrap();
                    assert!(out.iter().zip(expected).all(|(&a, b)| same(a, b)));
                    let expected = reference.intervals(name, start, start + 64).unwrap();
                    let mut entries =
                        vec![TrackEntry::default(); reader.query_entries(region).unwrap()];
                    reader.write_entries(&mut entries).unwrap();
                    assert_eq!(
                        entries
                            .iter()
                            .map(|r| (r.start, r.end, r.value.to_bits()))
                            .collect::<Vec<_>>(),
                        expected
                            .into_iter()
                            .map(|(s, e, v)| (s.max(start), e.min(start + 64), v.to_bits()))
                            .collect::<Vec<_>>()
                    );
                }
            }
        }
    }
}
