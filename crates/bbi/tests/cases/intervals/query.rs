use crate::support::Rng;
use crate::support::bigbed::{Record, Spec, write};
use crate::support::libbigwig::Reader;
use crate::support::random::regions;
use omicsweave_bbi::BbiReader;
use omicsweave_core::{
    Error, IntervalEntry, IntervalRead, IntervalsOut, Region, StringColumn, Strings, TrackRead,
};
use proptest::prelude::*;
proptest! {
    #![proptest_config(ProptestConfig { cases: 24, ..ProptestConfig::default() })]

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
        prop_assert!(matches!(ours.query_entries(Region::new(&spec.chroms[0].0, 0, 1), true), Err(Error::Unsupported(_))));
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
    let size = reader.query(region, false).unwrap();
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
