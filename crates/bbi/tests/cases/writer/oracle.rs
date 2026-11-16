//! The writers against libBigWig and the crate's readers: random files
//! written by [`BigWigWriter`] and [`BigBedWriter`] read back by both, zoom
//! levels checked against per-base truth, and the writers' sidecars checked
//! against [`BbiSidecar::build`] on the written file.

use std::fs::File;
use std::io::BufWriter;

use crate::support::libbigwig::Reader;
use crate::support::{Opened, Rng, same};
use omicsweave_bbi::{
    BbiReader, BedSchema, BigBedWriter, BigWigWriter, Chrom, Level, Summary, WriteOptions,
};
use omicsweave_core::{
    FileSource, IntervalRead, IntervalsOut, Locus, Region, StringColumn, TrackRead,
};
use proptest::prelude::*;

fn random_chroms(rng: &mut Rng) -> Vec<Chrom> {
    (0..rng.range(1, 4))
        .map(|i| Chrom {
            name: if rng.chance(30) {
                format!("scaffold_{i}")
            } else {
                format!("chr{}", rng.range(1, 30) * 10 + i)
            },
            len: rng.range(500, 60_000) as u32,
        })
        .collect()
}

fn random_options(rng: &mut Rng) -> WriteOptions {
    WriteOptions {
        items_per_slot: rng.chance(50).then(|| rng.range(1, 64) as u32),
        block_size: rng.range(2, 9) as u32,
        compress: rng.chance(80),
        level: rng.range(1, 12) as i32,
        zoom_levels: rng.range(0, 10) as usize,
        zoom_reductions: rng.chance(20).then(|| {
            (0..rng.range(1, 4))
                .map(|_| rng.range(1, 5_000) as u32)
                .collect()
        }),
    }
}

/// Checks zoom summaries against per-base truth: `depth[c][i]` is `None`
/// where base `i` of chromosome `c` is uncovered, else its value. A level
/// with reduction `r` must cover every covered base exactly once, each
/// summary starting at the first covered base after the previous one and
/// spanning `r` bases clipped to the chromosome.
fn check_zoom(records: &[Summary], truth: &[Option<f64>], chrom: u32, reduction: u32) {
    let len = truth.len() as u32;
    let mut at = 0u32;
    for r in records {
        assert_eq!(r.chrom, chrom);
        let next = (at..len)
            .find(|&i| truth[i as usize].is_some())
            .expect("summary past coverage");
        assert_eq!(r.start, next, "summary start");
        assert_eq!(
            r.end,
            r.start.saturating_add(reduction).min(len),
            "summary end"
        );
        let window = &truth[r.start as usize..r.end as usize];
        let values: Vec<f64> = window.iter().flatten().copied().collect();
        assert_eq!(r.valid as usize, values.len(), "covered bases of {r:?}");
        let min = values.iter().copied().fold(f64::INFINITY, f64::min);
        let max = values.iter().copied().fold(f64::NEG_INFINITY, f64::max);
        assert_eq!((r.min, r.max), (min as f32, max as f32), "range of {r:?}");
        let sum: f64 = values.iter().sum();
        assert!(
            (r.sum as f64 - sum).abs() <= 1e-3 * sum.abs().max(1.0),
            "sum of {r:?}: {sum}"
        );
        at = r.end;
    }
    assert!(
        (at..len).all(|i| truth[i as usize].is_none()),
        "coverage after the last summary"
    );
}

use crate::support::sidecar_check::check_sidecar;

fn bigwig_case(seed: u64) {
    let mut rng = Rng::new(seed);
    let dir = tempfile::tempdir().unwrap();
    let path = dir.path().join("t.bw");
    let chroms = random_chroms(&mut rng);
    let options = random_options(&mut rng);
    let mut w = BigWigWriter::create(
        BufWriter::new(File::create(&path).unwrap()),
        chroms.clone(),
        options.clone(),
    )
    .unwrap();
    let mut items: Vec<Vec<(u32, u32, f32)>> = vec![Vec::new(); chroms.len()];
    for (id, c) in chroms.iter().enumerate() {
        if rng.chance(15) {
            continue;
        }
        let mut pos = rng.range(0, 300) as u32;
        let dense = rng.chance(50);
        while pos < c.len && items[id].len() < 3_000 {
            let end = (pos + rng.range(1, if dense { 5 } else { 80 }) as u32).min(c.len);
            let v = rng.value();
            w.add(id as u32, pos, end, v).unwrap();
            items[id].push((pos, end, v));
            pos = end + rng.range(0, if dense { 2 } else { 60 }) as u32;
        }
    }
    let (out, sidecar) = w.finish().unwrap();
    drop(out.into_inner().unwrap());

    let mut theirs = Reader::open(&path).expect("libBigWig opens the file");
    let h = theirs.header();
    assert_eq!((h.kind, h.version), (0, 4));
    let lib_chroms = theirs.chroms();
    assert_eq!(lib_chroms.len(), chroms.len());
    for (name, len) in &lib_chroms {
        assert!(chroms.iter().any(|c| &c.name == name && c.len == *len));
    }
    let (mut bases, mut min, mut max, mut sum, mut squares) =
        (0u64, f64::INFINITY, f64::NEG_INFINITY, 0.0, 0.0);
    for list in &items {
        for &(s, e, v) in list {
            let n = (e - s) as f64;
            let v = v as f64;
            bases += (e - s) as u64;
            min = min.min(v);
            max = max.max(v);
            sum += n * v;
            squares += n * v * v;
        }
    }
    assert_eq!(h.bases_covered, bases);
    if bases > 0 {
        assert_eq!(
            (h.min, h.max, h.sum, h.sum_squares),
            (min, max, sum, squares)
        );
    }

    let mut reader = BbiReader::open(FileSource::open(&path).unwrap()).unwrap();
    let mut ours = Opened::open(&path);
    assert_eq!(ours.zooms.len(), h.zoom_levels as usize);
    if let Some(r) = &options.zoom_reductions {
        assert_eq!(
            ours.zooms.iter().map(|z| z.reduction).collect::<Vec<_>>(),
            *r
        );
    }
    for (id, c) in chroms.iter().enumerate() {
        let mut truth = vec![None; c.len as usize];
        for &(s, e, v) in &items[id] {
            truth[s as usize..e as usize].fill(Some(v as f64));
        }
        let lib = theirs.values(&c.name, 0, c.len).expect("libBigWig values");
        let mut mine = vec![0.0; c.len as usize];
        reader
            .values(&Locus::from(Region::new(&c.name, 0, c.len)), &mut mine)
            .unwrap();
        for (i, t) in truth.iter().enumerate() {
            let want = t.map_or(f32::NAN, |v| v as f32);
            assert!(same(lib[i], want) && same(mine[i], want), "{}:{i}", c.name);
        }
        assert_eq!(
            theirs.intervals(&c.name, 0, c.len).unwrap_or_default(),
            items[id]
        );
        for z in 0..ours.zooms.len() {
            let records = ours.summaries(z, id as u32, 0, c.len);
            check_zoom(&records, &truth, id as u32, ours.zooms[z].reduction);
        }
    }
    let base = sidecar.level(Level::Base).unwrap();
    assert_eq!(
        base.blocks.iter().map(|b| b.items as usize).sum::<usize>(),
        items.iter().map(Vec::len).sum()
    );
    assert!(base.blocks.windows(2).all(|w| w[0].start.0 <= w[1].start.0));
    check_sidecar(&path, &sidecar);
}

fn bigbed_case(seed: u64) {
    let mut rng = Rng::new(seed);
    let dir = tempfile::tempdir().unwrap();
    let path = dir.path().join("t.bb");
    let chroms = random_chroms(&mut rng);
    let options = random_options(&mut rng);
    let schema = BedSchema {
        field_count: 5,
        defined_field_count: 5,
        auto_sql: rng.chance(50).then(|| "table t\n\"test\"\n(\n)".to_owned()),
    };
    let mut w = BigBedWriter::create(
        BufWriter::new(File::create(&path).unwrap()),
        chroms.clone(),
        schema.clone(),
        options,
    )
    .unwrap();
    let mut items: Vec<Vec<(u32, u32, Vec<u8>)>> = vec![Vec::new(); chroms.len()];
    for (id, c) in chroms.iter().enumerate() {
        if rng.chance(15) {
            continue;
        }
        let mut pos = rng.range(0, 300) as u32;
        while pos < c.len && items[id].len() < 1_500 {
            let max = if rng.chance(10) { 3_000 } else { 120 };
            let end = (pos + rng.range(1, max) as u32).min(c.len);
            let rest = format!("f{}\t{}", rng.next() % 1000, rng.range(0, 1000)).into_bytes();
            w.add(id as u32, pos, end, &rest).unwrap();
            items[id].push((pos, end, rest));
            pos += rng.range(0, 90) as u32;
        }
    }
    let (out, sidecar) = w.finish().unwrap();
    drop(out.into_inner().unwrap());

    let mut theirs = Reader::open(&path).expect("libBigWig opens the file");
    let h = theirs.header();
    assert_eq!(
        (h.kind, h.version, h.field_count, h.defined_field_count),
        (1, 4, 5, 5)
    );
    assert_eq!(theirs.auto_sql(), schema.auto_sql);
    let mut reader = BbiReader::open(FileSource::open(&path).unwrap()).unwrap();
    let mut ours = Opened::open(&path);
    let mut covered = 0u64;
    for (id, c) in chroms.iter().enumerate() {
        let lib = theirs.entries(&c.name, 0, c.len).unwrap_or_default();
        let got: Vec<_> = lib
            .iter()
            .map(|e| (e.start, e.end, e.rest.clone()))
            .collect();
        assert_eq!(got, items[id]);

        let size = reader.query(Region::new(&c.name, 0, c.len), false).unwrap();
        let mut entries = vec![Default::default(); size.entries];
        let mut offsets = vec![0i64; size.entries + 1];
        let mut data = vec![0u8; size.annotation_bytes];
        reader
            .write(IntervalsOut {
                entries: &mut entries,
                annotations: StringColumn {
                    offsets: &mut offsets,
                    data: &mut data,
                },
            })
            .unwrap();
        let mine: Vec<(u32, u32)> = entries.iter().map(|e| (e.start, e.end)).collect();
        let lib_spans: Vec<(u32, u32)> = lib.iter().map(|e| (e.start, e.end)).collect();
        assert_eq!(mine, lib_spans);

        let mut depth = vec![0u32; c.len as usize];
        for (s, e, _) in &items[id] {
            depth[*s as usize..*e as usize]
                .iter_mut()
                .for_each(|d| *d += 1);
        }
        let truth: Vec<Option<f64>> = depth.iter().map(|&d| (d > 0).then_some(d as f64)).collect();
        covered += truth.iter().flatten().count() as u64;
        for z in 0..ours.zooms.len() {
            let records = ours.summaries(z, id as u32, 0, c.len);
            check_zoom(&records, &truth, id as u32, ours.zooms[z].reduction);
        }
    }
    assert_eq!(h.bases_covered, covered);
    let records: usize = items.iter().map(Vec::len).sum();
    assert_eq!(
        sidecar
            .level(Level::Base)
            .unwrap()
            .blocks
            .iter()
            .map(|b| b.items as usize)
            .sum::<usize>(),
        records
    );
    check_sidecar(&path, &sidecar);
}

proptest! {
    #![proptest_config(ProptestConfig { cases: 48, ..ProptestConfig::default() })]

    #[test]
    fn bigwig_writer_matches_libbigwig(seed in any::<u64>()) {
        bigwig_case(seed);
    }

    #[test]
    fn bigbed_writer_matches_libbigwig(seed in any::<u64>()) {
        bigbed_case(seed);
    }
}
