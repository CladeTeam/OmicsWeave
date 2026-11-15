//! The standard access of `BbiReader` against libBigWig, and bigBed against
//! a brute-force overlap scan of the records it was written from.

mod common;

use common::bigbed::{Record, Spec, write};
use common::libbigwig::{Reader, Writer};
use common::{Rng, same};
use omicsweave_bbi::BbiReader;
use omicsweave_core::{
    Error, IntervalEntry, IntervalRead, IntervalsOut, Locus, Region, StringColumn, Strings,
    TrackEntry, TrackRead,
};
use proptest::prelude::*;

#[test]
fn dense_values_preserve_last_record_wins_across_overlap_shapes() {
    for kind in ["disjoint", "monotone", "nested"] {
        let dir = tempfile::tempdir().unwrap();
        let path = dir.path().join("overlap.bw");
        let records: Vec<_> = (0..4096u32)
            .map(|i| {
                let (s, e) = match kind {
                    "disjoint" => (i * 4, i * 4 + 2),
                    "monotone" => (i, i + 8192),
                    _ => (i, if i % 2 == 0 { 16384 - i } else { i + 8192 }),
                };
                (s, e, if i % 127 == 0 { f32::NAN } else { i as f32 })
            })
            .collect();
        let mut writer = Writer::create(&path, 0, &[("s", 32768)]).unwrap();
        assert!(writer.add_bedgraph("s", &records));
        drop(writer);
        let mut reader = BbiReader::open(std::fs::read(path).unwrap()).unwrap();
        for (start, end) in [
            (5000, 9000),
            (0, 16384),
            (4095, 4096),
            (100, 100),
            (20000, 20100),
            (512, 768),
        ] {
            let region = Region::new("s", start, end);
            let mut out = vec![0.0; region.len() as usize];
            reader.values(&Locus::from(region), &mut out).unwrap();
            let mut expected = vec![f32::NAN; out.len()];
            for &(s, e, v) in &records {
                if s.max(start) < e.min(end) {
                    expected[(s.max(start) - start) as usize..(e.min(end) - start) as usize]
                        .fill(v);
                }
            }
            assert!(
                out.iter()
                    .zip(expected)
                    .all(|(a, b)| a.to_bits() == b.to_bits()),
                "{kind} {start}..{end}"
            );
        }
    }
}

fn random_bigwig(path: &std::path::Path, rng: &mut Rng) -> Vec<(String, u32)> {
    let chroms: Vec<(String, u32)> = (0..rng.range(1, 3))
        .map(|i| (format!("chr{i}"), rng.range(1_000, 100_000) as u32))
        .collect();
    let refs: Vec<(&str, u32)> = chroms.iter().map(|(n, l)| (n.as_str(), *l)).collect();
    let mut w = Writer::create(path, rng.range(0, 4) as i32, &refs).unwrap();
    for (name, len) in &chroms {
        let mut pos = 0u32;
        for _ in 0..rng.range(1, 40) {
            let mut items = Vec::new();
            for _ in 0..rng.range(1, 200) {
                pos += rng.range(0, 40) as u32;
                let end = (pos + rng.range(1, 60) as u32).min(*len);
                if pos >= end {
                    break;
                }
                items.push((pos, end, rng.value()));
                pos = end;
            }
            if !items.is_empty() {
                assert!(w.add_bedgraph(name, &items));
            }
        }
    }
    chroms
}

fn regions(rng: &mut Rng, chroms: &[(String, u32)]) -> Vec<(String, u32, u32)> {
    let mut out: Vec<(String, u32, u32)> = chroms.iter().map(|(n, l)| (n.clone(), 0, *l)).collect();
    for _ in 0..30 {
        let (name, len) = &chroms[rng.range(0, chroms.len() as u64 - 1) as usize];
        let start = rng.range(0, *len as u64 + 100) as u32;
        out.push((name.clone(), start, start + rng.range(0, 5_000) as u32));
    }
    out
}

proptest! {
    #![proptest_config(ProptestConfig { cases: 24, ..ProptestConfig::default() })]

    #[test]
    fn track_access_matches_libbigwig(seed in any::<u64>()) {
        let mut rng = Rng::new(seed);
        let dir = tempfile::tempdir().unwrap();
        let path = dir.path().join("t.bw");
        let chroms = random_bigwig(&path, &mut rng);
        let bytes = std::fs::read(&path).unwrap();
        let mut ours = BbiReader::open(bytes).unwrap();
        let mut theirs = Reader::open(&path).unwrap();

        for (chrom, start, end) in regions(&mut rng, &chroms) {
            let region = Region::new(&chrom, start, end);
            let mut values = vec![0f32; region.len() as usize];
            ours.values(&Locus::from(region), &mut values).unwrap();
            let n = ours.query_entries(region).unwrap();
            let mut entries = vec![TrackEntry::default(); n];
            ours.write_entries(&mut entries).unwrap();
            if region.is_empty() {
                prop_assert!(values.is_empty() && entries.is_empty());
                continue;
            }

            let expected_values = theirs.values(&chrom, start, end).unwrap();
            prop_assert!(values.iter().zip(&expected_values).all(|(a, b)| same(*a, *b)), "{}:{}-{}", chrom, start, end);
            let expected: Vec<TrackEntry> = theirs
                .intervals(&chrom, start, end)
                .unwrap()
                .into_iter()
                .map(|(s, e, value)| TrackEntry { start: s.max(start), end: e.min(end), value })
                .collect();
            prop_assert_eq!(&entries, &expected, "{}:{}-{}", chrom, start, end);
            prop_assert!(entries.iter().all(|e| e.start >= start && e.end <= end && e.start < e.end));
        }

        prop_assert_eq!(ours.chrom_len("chrUn"), Err(Error::UnknownChrom("chrUn".into())));
        let mut none: [f32; 0] = [];
        prop_assert!(matches!(ours.values(&Locus::from(Region::new("chrUn", 0, 0)), &mut none), Err(Error::UnknownChrom(_))));
        prop_assert!(matches!(ours.query(Region::new("chr0", 0, 10)), Err(Error::Unsupported(_))));
    }

    #[test]
    fn interval_access_matches_libbigwig_and_a_scan(seed in any::<u64>()) {
        let mut rng = Rng::new(seed);
        let spec = random_spec(&mut rng);
        prop_assume!(!spec.records.is_empty());
        let bigbed = write(&spec);
        let dir = tempfile::tempdir().unwrap();
        let path = dir.path().join("t.bb");
        std::fs::write(&path, &bigbed).unwrap();

        let mut ours = BbiReader::open(&bigbed[..]).unwrap();
        let mut theirs = Reader::open(&path).unwrap();

        let mut queries = regions(&mut rng, &spec.chroms);
        queries.push(("chrUn".into(), 0, 1_000));
        for (chrom, start, end) in queries {
            let region = Region::new(&chrom, start, end);
            let a = answer(&mut ours, region);
            let got: Vec<(u32, u32, Vec<u8>)> =
                a.0.iter().zip(&a.1).map(|(e, note)| (e.start, e.end, note.clone())).collect();
            prop_assert_eq!(&got, &scan(&spec, region), "{}:{}-{}", chrom, start, end);
            if region.is_empty() || chrom == "chrUn" {
                prop_assert!(a.0.is_empty());
                continue;
            }
            let expected: Vec<(u32, u32, Vec<u8>)> =
                theirs.entries(&chrom, start, end).unwrap().into_iter().map(|e| (e.start, e.end, e.rest)).collect();
            prop_assert_eq!(got, expected, "{}:{}-{}", chrom, start, end);
        }
        prop_assert!(matches!(ours.query_entries(Region::new(&spec.chroms[0].0, 0, 1)), Err(Error::Unsupported(_))));
    }
}

fn random_spec(rng: &mut Rng) -> Spec {
    let chroms: Vec<(String, u32)> = (0..rng.range(1, 5))
        .map(|i| (format!("chr{i}"), rng.range(10_000, 200_000) as u32))
        .collect();
    let extra = rng.range(0, 3) as usize;
    let mut records = Vec::new();
    for (id, (_, len)) in chroms.iter().enumerate() {
        let mut pos = 0u32;
        for _ in 0..rng.range(0, 300) {
            pos += rng.range(0, 300) as u32;
            if pos >= *len {
                break;
            }
            let end = (pos
                + if rng.chance(3) {
                    0
                } else {
                    rng.range(1, 2_000) as u32
                })
            .min(*len);
            let rest = (0..extra)
                .map(|f| format!("f{f}_{}", rng.next() % 1000))
                .collect::<Vec<_>>()
                .join("\t")
                .into_bytes();
            records.push(Record {
                chrom: id as u32,
                start: pos,
                end,
                rest,
            });
        }
    }
    Spec {
        chroms,
        records,
        items_per_block: rng.range(1, 32) as usize,
        split_chroms: rng.chance(50),
        chrom_fanout: rng.range(2, 6) as usize,
        index_fanout: rng.range(2, 8) as usize,
        auto_sql: None,
        compress: rng.chance(80),
        field_count: 3 + extra as u16,
        defined_field_count: 3,
    }
}

/// The records of `region` in file order: those starting before its end
/// and ending after its start (an empty region has none).
fn scan(spec: &Spec, region: Region<'_>) -> Vec<(u32, u32, Vec<u8>)> {
    let Some(id) = spec
        .chroms
        .iter()
        .position(|(name, _)| name == region.chrom)
    else {
        return Vec::new();
    };
    spec.records
        .iter()
        .filter(|r| {
            r.chrom == id as u32
                && !region.is_empty()
                && r.start < region.end
                && r.end > region.start
        })
        .map(|r| (r.start, r.end, r.rest.clone()))
        .collect()
}

/// A query's entries and their annotations.
fn answer(
    reader: &mut impl IntervalRead,
    region: Region<'_>,
) -> (Vec<IntervalEntry>, Vec<Vec<u8>>) {
    let size = reader.query(region).unwrap();
    let mut entries = vec![IntervalEntry::default(); size.entries];
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
    let strings = Strings {
        offsets: &offsets,
        data: &data,
    };
    let notes = entries
        .iter()
        .map(|e| strings.get(e.annotation as usize).to_vec())
        .collect();
    (entries, notes)
}
