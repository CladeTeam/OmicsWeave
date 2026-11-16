use crate::support::*;
use omicsweave_core::FileSource;
use omicsweave_core::cache::BlockCache;
use omicsweave_core::{
    AlignmentEntry, AlignmentFields, AlignmentRead, AlignmentsOut, Error, Region, StringColumn,
};
use omicsweave_cram::index_entries;
use std::sync::Arc;

#[test]
fn multi_reference_slices() {
    let refs: Vec<(String, u64)> = (0..40)
        .map(|i| (format!("ctg{i}"), 6_000 + 97 * i))
        .collect();
    let f = fixture(
        3,
        refs,
        SamSpec {
            per_ref: 12,
            gap: 100,
            long_cigars: false,
            unplaced: 10,
        },
    );
    for (name, opts) in [
        (
            "multi",
            &["version=3.1", "multi_seq_per_slice=1", "seqs_per_slice=150"][..],
        ),
        (
            "multi30",
            &[
                "version=3.0",
                "multi_seq_per_slice=1",
                "seqs_per_slice=1000",
                "embed_ref=1",
            ],
        ),
    ] {
        let (cram, crai) = f.cram(name, opts, true);
        assert!(check_queries(&f, &cram, &crai, true, 9) > 300, "{name}");
        let hts = Hts::open(&cram, &crai, Some(&f.fasta));
        let mut cached = open(&cram, &crai, Some(&f.fasta));
        cached.set_cache(Some(Arc::new(BlockCache::new(64 << 20, 4096))));
        for _ in 0..2 {
            for chrom in ["ctg0", "ctg4", "ctg17", "ctg39"] {
                let tid = hts.tid(chrom);
                let region = Region::new(chrom, 0, 2000);
                let expected = hts.query(tid, 0, 2000);
                assert_eq!(
                    query(&mut cached, tid, region, AlignmentFields::ALL),
                    expected
                );
            }
        }
        let want = String::from_utf8(gunzip(&std::fs::read(&crai).unwrap())).unwrap();
        let got: String = index_entries(FileSource::open(&cram).unwrap())
            .unwrap()
            .iter()
            .map(|e| e.line())
            .collect();
        assert_eq!(got, want, "{name}");
        assert!(
            want.lines().any(|l| l.starts_with("-1\t")),
            "{name}: unplaced entries"
        );
    }
}

fn entries_out(entries: &mut [AlignmentEntry]) -> AlignmentsOut<'_> {
    AlignmentsOut {
        entries,
        names: None,
        seqs: None,
        quals: None,
        cigars: None,
        aux: None,
        modifications: None,
    }
}

#[test]
fn writes_require_exact_buffers_and_keep_staged_records_after_failure() {
    let f = fixture(
        81,
        vec![("chr1".into(), 20_000)],
        SamSpec {
            per_ref: 30,
            gap: 20,
            long_cigars: false,
            unplaced: 0,
        },
    );
    let (cram, crai) = f.cram("output", &["version=3.1"], true);
    let mut reader = open(&cram, &crai, Some(&f.fasta));
    assert_eq!(reader.write(entries_out(&mut [])), Err(Error::NotQueried));
    let region = Region::new("chr1", 0, 20_000);
    let size = reader.query(region, AlignmentFields::NAME).unwrap();
    assert!(size.records > 0 && size.name_bytes > 0);
    let sentinel = AlignmentEntry {
        start: -123,
        ..AlignmentEntry::default()
    };
    for len in [size.records - 1, size.records + 1] {
        let mut entries = vec![sentinel; len];
        assert!(
            matches!(
                reader.write(entries_out(&mut entries)),
                Err(Error::OutputSize { .. })
            ),
            "entries={len}"
        );
        assert!(entries.iter().all(|&entry| entry == sentinel));
    }
    let mut entries = vec![sentinel; size.records];
    assert!(matches!(
        reader.write(entries_out(&mut entries)),
        Err(Error::Invalid { .. })
    ));
    assert!(entries.iter().all(|&entry| entry == sentinel));
    for (offset_len, data_len) in [
        (size.records, size.name_bytes),
        (size.records + 1, size.name_bytes - 1),
        (size.records + 1, size.name_bytes + 1),
    ] {
        let mut offsets = vec![-1; offset_len];
        let mut data = vec![0xfe; data_len];
        let mut out = entries_out(&mut entries);
        out.names = Some(StringColumn {
            offsets: &mut offsets,
            data: &mut data,
        });
        assert!(
            matches!(reader.write(out), Err(Error::OutputSize { .. })),
            "offsets={offset_len}, data={data_len}"
        );
        assert!(entries.iter().all(|&entry| entry == sentinel));
        assert!(offsets.iter().all(|&offset| offset == -1));
        assert!(data.iter().all(|&byte| byte == 0xfe));
    }
    let oracle = Hts::open(&cram, &crai, Some(&f.fasta)).query(0, 0, 20_000);
    let expected_names: Vec<_> = oracle
        .iter()
        .flat_map(|record| record.name.iter().copied())
        .collect();
    let mut expected_offsets = vec![0];
    for record in &oracle {
        expected_offsets.push(expected_offsets.last().unwrap() + record.name.len() as i64);
    }
    for _ in 0..2 {
        let mut offsets = vec![-1; size.records + 1];
        let mut data = vec![0xfe; size.name_bytes];
        let mut out = entries_out(&mut entries);
        out.names = Some(StringColumn {
            offsets: &mut offsets,
            data: &mut data,
        });
        reader.write(out).unwrap();
        assert_eq!(data, expected_names);
        assert_eq!(offsets, expected_offsets);
        assert_eq!(
            entries.iter().map(|entry| entry.start).collect::<Vec<_>>(),
            oracle.iter().map(|record| record.pos).collect::<Vec<_>>()
        );
    }
}

#[test]
fn failed_queries_invalidate_staged_output_and_reader_recovers() {
    let f = fixture(
        82,
        vec![("chr1".into(), 20_000)],
        SamSpec {
            per_ref: 30,
            gap: 20,
            long_cigars: false,
            unplaced: 0,
        },
    );
    let (cram, crai) = f.cram("state", &["version=3.0"], true);
    let mut reader = open(&cram, &crai, Some(&f.fasta));
    let region = Region::new("chr1", 0, 20_000);
    let expected = Hts::open(&cram, &crai, Some(&f.fasta)).query(0, 0, 20_000);
    for invalid in [Region::new("missing", 0, 10), Region::new("chr1", 10, 9)] {
        assert_eq!(
            query(&mut reader, 0, region, AlignmentFields::ALL),
            expected
        );
        let error = reader.query(invalid, AlignmentFields::NONE).unwrap_err();
        assert!(
            matches!(error, Error::UnknownChrom(_) | Error::InvalidRegion { .. }),
            "{invalid:?}: {error}"
        );
        assert_eq!(reader.write(entries_out(&mut [])), Err(Error::NotQueried));
        assert_eq!(
            query(&mut reader, 0, region, AlignmentFields::ALL),
            expected
        );
    }
    assert_eq!(
        reader
            .query(Region::new("chr1", 100, 100), AlignmentFields::NONE)
            .unwrap()
            .records,
        0
    );
    reader.write(entries_out(&mut [])).unwrap();
}
