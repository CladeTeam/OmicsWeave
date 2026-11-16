use crate::support::byte_fixture::*;
use omicsweave_bbi::{BIGBED_MAGIC, BIGWIG_MAGIC, BbiReader};
use omicsweave_core::{
    IntervalEntry, IntervalRead, IntervalsOut, Region, StringColumn, TrackEntry, TrackRead,
};
use std::{cell::Cell, rc::Rc};
#[test]
fn high_level_bigbed_orders_all_hits_stably_across_blocks() {
    let records: Vec<(u32, u32, &[u8])> = vec![
        (90, 100, b"miss"),
        (20, 30, b"a"),
        (2, 8, b"b"),
        (20, 25, b"c"),
        (1, 3, b"d"),
    ];
    let blocks = [bed(&records[..3]), bed(&records[3..])];
    let mut expected: Vec<_> = records
        .iter()
        .copied()
        .filter(|r| r.0 < 40 && r.1 > 0)
        .collect();
    expected.sort_by_key(|r| r.0);
    for compressed in [false, true] {
        let mut reader = BbiReader::open(container(BIGBED_MAGIC, &blocks, compressed)).unwrap();
        let size = reader.query(Region::new("s", 0, 40), false).unwrap();
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
        assert_eq!(entries.len(), expected.len());
        for (row, &(start, end, note)) in expected.iter().enumerate() {
            assert_eq!(
                entries[row],
                IntervalEntry {
                    start,
                    end,
                    annotation: row as u32
                }
            );
            assert_eq!(
                &data[offsets[row] as usize..offsets[row + 1] as usize],
                note
            );
        }
    }
}

#[test]
fn malformed_bigbed_tail_invalidates_staged_output() {
    let mut bytes = bed(&[(90, 100, b"miss")]);
    bytes.push(0);
    let mut reader = BbiReader::open(container(BIGBED_MAGIC, &[bytes], false)).unwrap();
    assert!(reader.query(Region::new("s", 0, 1), false).is_err());
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
}

#[test]
fn failed_reads_revoke_staging_and_allow_recovery() {
    let fail = Rc::new(Cell::new(false));
    let mut track = BbiReader::open(Failable {
        bytes: container(BIGWIG_MAGIC, &[wig(&[(1, 3, 0x4040_0000)])], true),
        fail: fail.clone(),
    })
    .unwrap();
    let region = Region::new("s", 0, 40);
    assert_eq!(track.query_entries(region, true).unwrap(), 1);
    fail.set(true);
    assert!(track.query_entries(region, true).is_err());
    assert!(track.write_entries(&mut [TrackEntry::default()]).is_err());
    fail.set(false);
    assert_eq!(track.query_entries(region, true).unwrap(), 1);
    let mut entries = [TrackEntry::default()];
    track.write_entries(&mut entries).unwrap();
    assert_eq!(
        entries[0],
        TrackEntry {
            start: 1,
            end: 3,
            value: 3.0
        }
    );

    let mut intervals = BbiReader::open(Failable {
        bytes: container(BIGBED_MAGIC, &[bed(&[(1, 3, b"a")])], true),
        fail: fail.clone(),
    })
    .unwrap();
    assert_eq!(intervals.query(region, false).unwrap().entries, 1);
    fail.set(true);
    assert!(intervals.query(region, false).is_err());
    assert!(
        intervals
            .write(IntervalsOut {
                entries: &mut [IntervalEntry::default()],
                annotations: StringColumn {
                    offsets: &mut [0; 2],
                    data: &mut [0]
                },
            })
            .is_err()
    );
    fail.set(false);
    assert_eq!(intervals.query(region, false).unwrap().entries, 1);
    let mut entries = [IntervalEntry::default()];
    let mut offsets = [0; 2];
    let mut data = [0];
    intervals
        .write(IntervalsOut {
            entries: &mut entries,
            annotations: StringColumn {
                offsets: &mut offsets,
                data: &mut data,
            },
        })
        .unwrap();
    assert_eq!(
        entries[0],
        IntervalEntry {
            start: 1,
            end: 3,
            annotation: 0
        }
    );
    assert_eq!(offsets, [0, 1]);
    assert_eq!(&data, b"a");
}

#[test]
fn clip_cuts_entries_to_the_region_only_when_asked() {
    let source = |bytes| Failable {
        bytes,
        fail: Rc::new(Cell::new(false)),
    };
    let region = Region::new("s", 5, 15);
    let mut track = BbiReader::open(source(container(
        BIGWIG_MAGIC,
        &[wig(&[(0, 10, 0x3f80_0000), (12, 30, 0x4000_0000)])],
        false,
    )))
    .unwrap();
    for (clip, bounds) in [(true, [(5, 10), (12, 15)]), (false, [(0, 10), (12, 30)])] {
        let mut entries = vec![TrackEntry::default(); track.query_entries(region, clip).unwrap()];
        track.write_entries(&mut entries).unwrap();
        let got: Vec<_> = entries.iter().map(|e| (e.start, e.end)).collect();
        assert_eq!(got, bounds);
    }

    let mut features = BbiReader::open(source(container(
        BIGBED_MAGIC,
        &[bed(&[(0, 10, b"a"), (12, 30, b"b")])],
        false,
    )))
    .unwrap();
    for (clip, bounds) in [(false, [(0, 10), (12, 30)]), (true, [(5, 10), (12, 15)])] {
        let size = features.query(region, clip).unwrap();
        let mut entries = vec![IntervalEntry::default(); size.entries];
        let mut offsets = vec![0; size.entries + 1];
        let mut data = vec![0; size.annotation_bytes];
        features
            .write(IntervalsOut {
                entries: &mut entries,
                annotations: StringColumn {
                    offsets: &mut offsets,
                    data: &mut data,
                },
            })
            .unwrap();
        let got: Vec<_> = entries.iter().map(|e| (e.start, e.end)).collect();
        assert_eq!(got, bounds);
        assert_eq!(data, b"ab");
    }
}

mod recovery {
    use crate::support::compressed_fixture::{
        MutableSource, bed, bed_answer, compressed_container,
    };
    use omicsweave_bbi::{BIGBED_MAGIC, BbiReader};
    use omicsweave_core::{IntervalRead, IntervalsOut, Region, StringColumn};
    use std::{cell::RefCell, rc::Rc};
    #[test]
    fn bigbed_large_declared_ceiling_preserves_prefix_and_recovers_after_corruption() {
        let long = vec![b'x'; 180_000];
        let blocks = [bed(1, &long), bed(101, b"end")];
        let (bytes, ranges) = compressed_container(BIGBED_MAGIC, &blocks, u32::MAX, false);
        let source = MutableSource(Rc::new(RefCell::new(bytes)));
        let mut reader = BbiReader::open(source.clone()).unwrap();
        let full = Region::new("s", 0, 200);
        let short = Region::new("s", 100, 200);
        assert_eq!(
            bed_answer(&mut reader, full),
            [long.clone(), b"end".to_vec()]
        );
        assert_eq!(bed_answer(&mut reader, short), [b"end".to_vec()]);
        let original = source.0.borrow()[ranges[1].clone()].to_vec();
        source.0.borrow_mut()[ranges[1].clone()].fill(0);
        assert!(reader.query(full, false).is_err());
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
        source.0.borrow_mut()[ranges[1].clone()].copy_from_slice(&original);
        assert_eq!(bed_answer(&mut reader, full), [long, b"end".to_vec()]);
    }
}
