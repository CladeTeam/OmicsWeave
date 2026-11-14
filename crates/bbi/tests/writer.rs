//! The writers against libBigWig and the crate's readers: random files
//! written by [`BigWigWriter`] and [`BigBedWriter`] read back by both, zoom
//! levels checked against per-base truth, and the writers' sidecars checked
//! against [`BbiSidecar::build`] on the written file.

mod common;

use std::fs::File;
use std::io::BufWriter;
use std::path::Path;

use common::libbigwig::{Reader, Writer};
use common::{Opened, Rng, same};
use omicsweave_bbi::{
    BbiReader, BbiSidecar, BedSchema, BigBedWriter, BigWigWriter, Chrom, Kind, Level, Summary,
    WriteOptions,
};
use omicsweave_core::sidecar::{SourceIdentity, sidecar_path};
use omicsweave_core::{
    Error, FileSource, IntervalRead, IntervalsOut, Region, StringColumn, TrackRead,
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

fn check_sidecar(path: &Path, written: &BbiSidecar) {
    let source = FileSource::open(path).unwrap();
    let built = BbiSidecar::build(&source).unwrap();
    assert_eq!(
        &built, written,
        "writer sidecar differs from one built from the file"
    );
    let identity = SourceIdentity::of(&source).unwrap();
    let bytes = written.to_bytes(identity);
    assert_eq!(&BbiSidecar::parse(&bytes, identity).unwrap(), written);
    let other = SourceIdentity {
        fingerprint: identity.fingerprint ^ 1,
        ..identity
    };
    assert_eq!(BbiSidecar::parse(&bytes, other), Err(Error::StaleSidecar));
    let mut damaged = bytes.clone();
    *damaged.last_mut().unwrap() ^= 1;
    assert!(matches!(
        BbiSidecar::parse(&damaged, identity),
        Err(Error::Invalid { .. })
    ));
}

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
            .values(Region::new(&c.name, 0, c.len), &mut mine)
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

        let size = reader.query(Region::new(&c.name, 0, c.len)).unwrap();
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

#[test]
fn sidecar_of_a_libbigwig_file_and_its_path() {
    let dir = tempfile::tempdir().unwrap();
    let path = dir.path().join("lib.bw");
    {
        let mut w = Writer::create(&path, 6, &[("chr1", 100_000), ("chr2", 50_000)]).unwrap();
        let items: Vec<(u32, u32, f32)> =
            (0..5_000).map(|i| (i * 20, i * 20 + 7, i as f32)).collect();
        assert!(w.add_bedgraph("chr1", &items));
        assert!(w.add_fixedstep("chr2", 100, 5, 10, &vec![1.5; 2_000]));
    }
    let built = BbiSidecar::build_file(&path).unwrap();
    assert_eq!(BbiSidecar::load(&path).unwrap(), built);
    assert!(sidecar_path(&path).ends_with("lib.bw.sidecar"));

    let mut ours = Opened::open(&path);
    assert_eq!(built.kind, Kind::BigWig);
    assert_eq!(built.levels.len(), 1 + ours.zooms.len());
    let total: u32 = built.levels[0].blocks.iter().map(|b| b.items).sum();
    assert_eq!(total, 5_000 + 2_000);
    for (i, level) in built.levels.iter().enumerate() {
        for b in &level.blocks {
            assert!(b.inflated <= built.uncompress_buf_size);
        }
        if i > 0 {
            let records: u32 = level.blocks.iter().map(|b| b.items).sum();
            let found: usize = (0..2)
                .map(|c| ours.summaries(i - 1, c, 0, u32::MAX).len())
                .sum();
            assert_eq!(records as usize, found);
        }
    }

    // Any change to the file makes the saved sidecar stale.
    let mut bytes = std::fs::read(&path).unwrap();
    bytes.push(0);
    std::fs::write(&path, &bytes).unwrap();
    assert_eq!(BbiSidecar::load(&path), Err(Error::StaleSidecar));
}

#[test]
fn writers_reject_unsorted_and_invalid_items() {
    let chroms = vec![
        Chrom {
            name: "a".into(),
            len: 1000,
        },
        Chrom {
            name: "b".into(),
            len: 1000,
        },
    ];
    let mut w = BigWigWriter::create(
        std::io::Cursor::new(Vec::new()),
        chroms.clone(),
        WriteOptions::default(),
    )
    .unwrap();
    w.add(1, 10, 20, 1.0).unwrap();
    assert!(w.add(1, 15, 30, 1.0).is_err(), "overlap");
    assert!(w.add(0, 0, 5, 1.0).is_err(), "earlier chromosome");
    assert!(w.add(1, 30, 30, 1.0).is_err(), "empty item");
    assert!(w.add(1, 30, 1001, 1.0).is_err(), "past the chromosome");
    assert!(w.add(1, 30, 40, f32::NAN).is_err(), "NaN");
    assert!(w.add(2, 0, 1, 1.0).is_err(), "unknown chromosome");

    let mut b = BigBedWriter::create(
        std::io::Cursor::new(Vec::new()),
        chroms.clone(),
        BedSchema::bed(3),
        WriteOptions::default(),
    )
    .unwrap();
    b.add(0, 10, 20, b"").unwrap();
    b.add(0, 10, 15, b"").unwrap();
    assert!(b.add(0, 5, 15, b"").is_err(), "unsorted");
    assert!(b.add(0, 30, 40, b"x\0y").is_err(), "NUL");

    let dup = vec![chroms[0].clone(), chroms[0].clone()];
    assert!(
        BigWigWriter::create(
            std::io::Cursor::new(Vec::new()),
            dup,
            WriteOptions::default()
        )
        .is_err()
    );
}

#[test]
fn empty_files_are_readable() {
    let dir = tempfile::tempdir().unwrap();
    for kind in [Kind::BigWig, Kind::BigBed] {
        let path = dir.path().join(format!("{kind:?}"));
        let out = BufWriter::new(File::create(&path).unwrap());
        let chroms = vec![Chrom {
            name: "chr1".into(),
            len: 10,
        }];
        let (out, sidecar) = match kind {
            Kind::BigWig => BigWigWriter::create(out, chroms, WriteOptions::default())
                .unwrap()
                .finish()
                .unwrap(),
            Kind::BigBed => {
                BigBedWriter::create(out, chroms, BedSchema::bed(3), WriteOptions::default())
                    .unwrap()
                    .finish()
                    .unwrap()
            }
        };
        drop(out.into_inner().unwrap());
        assert_eq!(sidecar.levels.len(), 1);
        assert!(sidecar.levels[0].blocks.is_empty());
        assert!(Reader::open(&path).is_some());
        let mut reader = BbiReader::open(FileSource::open(&path).unwrap()).unwrap();
        if kind == Kind::BigWig {
            let mut v = [0.0; 10];
            reader.values(Region::new("chr1", 0, 10), &mut v).unwrap();
            assert!(v.iter().all(|x| x.is_nan()));
        }
        check_sidecar(&path, &sidecar);
    }
}
