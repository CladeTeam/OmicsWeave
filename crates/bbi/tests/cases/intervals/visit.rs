use std::{cell::Cell, rc::Rc, sync::Arc};

use crate::support::byte_fixture::{Failable, bed, container};
use omicsweave_bbi::{BIGBED_MAGIC, BbiReader};
use omicsweave_core::cache::BlockCache;
use omicsweave_core::{
    Error, IntervalEntry, IntervalRead, IntervalsOut, Region, StringColumn, Strings,
};

type Row = (u32, u32, Vec<u8>);

fn visited(reader: &BbiReader<impl omicsweave_core::ReadAt>) -> Vec<Row> {
    let mut rows = Vec::new();
    reader
        .visit_intervals(|start, end, note| {
            rows.push((start, end, note.to_vec()));
            Ok(())
        })
        .unwrap();
    rows
}

#[test]
fn visits_match_staged_output_with_clipping_and_cache_leases() {
    let blocks = [
        bed(&[(12, 20, b"late"), (0, 10, b"first"), (5, 14, b"")]),
        bed(&[(5, 15, b"\xff\x80\t"), (50, 60, b"miss")]),
    ];
    for compressed in [false, true] {
        for cached in [false, true] {
            let fail = Rc::new(Cell::new(false));
            let mut reader = BbiReader::open(Failable {
                bytes: container(BIGBED_MAGIC, &blocks, compressed),
                fail: fail.clone(),
            })
            .unwrap();
            let cache = Arc::new(BlockCache::new(1 << 20, 128));
            if cached {
                reader.set_cache(Some(cache.clone())).unwrap();
            }
            assert_eq!(
                reader.visit_intervals(|_, _, _| Ok(())),
                Err(Error::NotQueried)
            );
            for clip in [false, true] {
                fail.set(false);
                let size = reader.query(Region::new("s", 5, 15), clip).unwrap();
                let expected: Vec<Row> = [
                    (0u32, 10u32, b"first".as_slice()),
                    (5, 14, b"".as_slice()),
                    (5, 15, b"\xff\x80\t".as_slice()),
                    (12, 20, b"late".as_slice()),
                ]
                .into_iter()
                .map(|(start, end, note)| {
                    (
                        if clip { start.max(5) } else { start },
                        if clip { end.min(15) } else { end },
                        note.to_vec(),
                    )
                })
                .collect();
                let mut entries = vec![IntervalEntry::default(); size.entries];
                let mut offsets = vec![0; size.entries + 1];
                let mut data = vec![0; size.annotation_bytes];
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
                let copied: Vec<_> = entries
                    .iter()
                    .map(|row| {
                        (
                            row.start,
                            row.end,
                            strings.get(row.annotation as usize).to_vec(),
                        )
                    })
                    .collect();
                assert_eq!(
                    copied, expected,
                    "compressed={compressed}, cached={cached}, clip={clip}"
                );
                // Visiting and writing consume retained query state without new source reads.
                fail.set(true);
                cache.clear();
                assert_eq!(visited(&reader), expected);
                assert_eq!(visited(&reader), expected);
                assert!(
                    reader
                        .write(IntervalsOut {
                            entries: &mut [],
                            annotations: StringColumn {
                                offsets: &mut [0],
                                data: &mut []
                            },
                        })
                        .is_err()
                );
                assert_eq!(visited(&reader), expected);
            }
        }
    }
}

#[test]
fn callback_errors_stop_visiting_without_consuming_staged_rows() {
    let bytes = container(
        BIGBED_MAGIC,
        &[bed(&[(1, 3, b"a"), (4, 6, b"b"), (7, 9, b"c")])],
        false,
    );
    let mut reader = BbiReader::open(bytes).unwrap();
    reader.query(Region::new("s", 0, 10), false).unwrap();
    let expected = visited(&reader);
    let mut calls = 0;
    let error = Error::invalid("visitor", "stop after second row");
    assert_eq!(
        reader.visit_intervals(|_, _, _| {
            calls += 1;
            if calls == 2 {
                Err(error.clone())
            } else {
                Ok(())
            }
        }),
        Err(error)
    );
    assert_eq!(calls, 2);
    assert_eq!(visited(&reader), expected);
}

#[test]
fn empty_queries_visit_nothing_and_failed_queries_revoke_staging() {
    let fail = Rc::new(Cell::new(false));
    let mut reader = BbiReader::open(Failable {
        bytes: container(BIGBED_MAGIC, &[bed(&[(1, 3, b"a")])], true),
        fail: fail.clone(),
    })
    .unwrap();
    for region in [
        Region::new("s", 5, 5),
        Region::new("s", 10, 20),
        Region::new("missing", 0, 10),
    ] {
        assert_eq!(reader.query(region, false).unwrap().entries, 0);
        assert!(visited(&reader).is_empty());
    }
    reader.query(Region::new("s", 0, 10), false).unwrap();
    assert_eq!(visited(&reader), [(1, 3, b"a".to_vec())]);
    fail.set(true);
    assert!(reader.query(Region::new("s", 0, 10), false).is_err());
    let mut calls = 0;
    assert_eq!(
        reader.visit_intervals(|_, _, _| {
            calls += 1;
            Ok(())
        }),
        Err(Error::NotQueried)
    );
    assert_eq!(calls, 0);
    fail.set(false);
    reader.query(Region::new("s", 0, 10), false).unwrap();
    assert_eq!(visited(&reader), [(1, 3, b"a".to_vec())]);
}
