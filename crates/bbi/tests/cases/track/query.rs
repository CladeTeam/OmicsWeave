use crate::support::libbigwig::{Reader, Writer};
use crate::support::{Rng, same};
use omicsweave_bbi::BbiReader;
use omicsweave_core::{Error, IntervalRead, Locus, Region, TrackEntry, TrackRead};
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

use crate::support::random::regions;

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
            let n = ours.query_entries(region, true).unwrap();
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
        prop_assert!(matches!(ours.query(Region::new("chr0", 0, 10), false), Err(Error::Unsupported(_))));
    }
}
