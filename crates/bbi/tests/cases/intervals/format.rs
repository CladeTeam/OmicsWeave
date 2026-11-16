//! Differential tests against libBigWig's bigBed reader on random files from
//! the test writer, including multi-level trees and blocks that span
//! chromosomes.

use crate::support::bigbed::{Record, Spec, write};
use crate::support::libbigwig::Reader;
use crate::support::{Opened, Rng};
use omicsweave_bbi::Kind;
use proptest::prelude::*;

fn random_spec(rng: &mut Rng) -> Spec {
    let n = rng.range(1, 12) as usize;
    let chroms: Vec<(String, u32)> = (0..n)
        .map(|i| {
            (
                format!("chr{}_{}", i, rng.next() % 1000),
                rng.range(10_000, 1_000_000) as u32,
            )
        })
        .collect();
    let extra = rng.range(0, 4) as usize;
    let mut records = Vec::new();
    for (id, (_, len)) in chroms.iter().enumerate() {
        let mut pos = rng.range(0, 1_000) as u32;
        for _ in 0..rng.range(0, 400) {
            if pos >= *len {
                break;
            }
            let width = if rng.chance(3) {
                0
            } else {
                rng.range(1, 500) as u32
            };
            let end = (pos + width).min(*len);
            let rest = if extra == 0 || rng.chance(10) {
                Vec::new()
            } else {
                (0..extra)
                    .map(|f| format!("f{f}_{}", rng.next() % 10_000))
                    .collect::<Vec<_>>()
                    .join("\t")
                    .into_bytes()
            };
            records.push(Record {
                chrom: id as u32,
                start: pos,
                end,
                rest,
            });
            // Overlapping records are legal in bigBed; only starts are sorted.
            pos += rng.range(0, 200) as u32;
        }
    }
    Spec {
        chroms,
        records,
        items_per_block: rng.range(1, 64) as usize,
        split_chroms: rng.chance(50),
        chrom_fanout: rng.range(2, 8) as usize,
        index_fanout: rng.range(2, 16) as usize,
        auto_sql: rng.chance(50).then(|| {
            "table t\n\"test\"\n(\nstring chrom;\nuint chromStart;\nuint chromEnd;\n)".to_owned()
        }),
        compress: rng.chance(80),
        field_count: 3 + extra as u16,
        defined_field_count: 3,
    }
}

fn check_query(ours: &mut Opened, theirs: &mut Reader, chrom: &str, start: u32, end: u32) {
    let a = ours.entries(chrom, start, end);
    let b: Vec<(u32, u32, Vec<u8>)> = theirs
        .entries(chrom, start, end)
        .expect("libBigWig entries")
        .into_iter()
        .map(|e| (e.start, e.end, e.rest))
        .collect();
    assert_eq!(a, b, "{chrom}:{start}-{end}");
}

proptest! {
    #![proptest_config(ProptestConfig { cases: 48, ..ProptestConfig::default() })]

    #[test]
    fn matches_libbigwig(seed in any::<u64>()) {
        let mut rng = Rng::new(seed);
        let spec = random_spec(&mut rng);
        prop_assume!(!spec.records.is_empty());
        let dir = tempfile::tempdir().unwrap();
        let path = dir.path().join("t.bb");
        std::fs::write(&path, write(&spec)).unwrap();

        let mut ours = Opened::open(&path);
        let mut theirs = Reader::open(&path).expect("libBigWig open");

        let h = theirs.header();
        prop_assert_eq!(h.kind, 1);
        prop_assert_eq!(ours.header.kind, Kind::BigBed);
        prop_assert_eq!(ours.header.field_count, h.field_count);
        prop_assert_eq!(ours.header.defined_field_count, h.defined_field_count);
        prop_assert_eq!(ours.header.uncompress_buf_size, h.buf_size);
        let s = ours.summary.expect("total summary");
        prop_assert_eq!(s.bases_covered, h.bases_covered);
        prop_assert_eq!(&ours.auto_sql, &spec.auto_sql);
        prop_assert_eq!(&ours.auto_sql, &theirs.auto_sql());
        let chroms: Vec<(String, u32)> = ours.chroms.iter().map(|(_, c)| (c.name.clone(), c.len)).collect();
        prop_assert_eq!(&chroms, &theirs.chroms());
        prop_assert_eq!(&chroms, &spec.chroms);

        for (name, len) in &spec.chroms {
            check_query(&mut ours, &mut theirs, name, 0, *len);
        }
        for _ in 0..60 {
            let (name, len) = &spec.chroms[rng.range(0, spec.chroms.len() as u64 - 1) as usize];
            let start = rng.range(0, *len as u64) as u32;
            let end = start + rng.range(1, 50_000) as u32;
            check_query(&mut ours, &mut theirs, name, start, end);
        }

        // Every record comes back from a whole-chromosome query, so the
        // index loses nothing.
        let total: usize = spec.chroms.iter().map(|(name, len)| ours.entries(name, 0, *len).len()).sum();
        let zero_width_at_zero = spec.records.iter().filter(|r| r.start == 0 && r.end == 0).count();
        prop_assert_eq!(total + zero_width_at_zero, spec.records.len());
    }
}
