//! Independent byte-level oracle for staged BED output and annotation ordering.

use std::io::Write;

use omicsweave_bed::BedTable;
use omicsweave_core::{
    Error, IntervalEntry, IntervalRead, IntervalSize, IntervalsOut, Region, StringColumn,
};

#[derive(Clone)]
struct Record {
    chrom: &'static str,
    start: u32,
    end: u32,
    note: Vec<u8>,
}

fn render(records: &[Record]) -> Vec<u8> {
    let mut text = Vec::new();
    for record in records {
        write!(text, "{}\t{}\t{}\t", record.chrom, record.start, record.end).unwrap();
        text.extend_from_slice(&record.note);
        text.push(b'\n');
    }
    text
}

// File position is an explicit secondary key: no index/query code is reused.
fn expected(records: &[Record], region: Region<'_>) -> (Vec<IntervalEntry>, Vec<i64>, Vec<u8>) {
    let mut selected: Vec<_> = records
        .iter()
        .enumerate()
        .filter(|(_, r)| {
            region.start < region.end
                && r.chrom == region.chrom
                && r.start < region.end
                && r.end > region.start
        })
        .collect();
    selected.sort_unstable_by_key(|(file_row, r)| (r.start, *file_row));
    let mut entries = Vec::new();
    let mut offsets = vec![0];
    let mut data = Vec::new();
    for (_, r) in selected {
        entries.push(IntervalEntry {
            start: r.start,
            end: r.end,
            annotation: entries.len() as u32,
        });
        data.extend_from_slice(&r.note);
        offsets.push(data.len() as i64);
    }
    (entries, offsets, data)
}

fn assert_output(table: &BedTable, records: &[Record], region: Region<'_>, size: IntervalSize) {
    let (want_entries, want_offsets, want_data) = expected(records, region);
    assert_eq!(size.entries, want_entries.len());
    assert_eq!(size.annotation_bytes, want_data.len());
    let guard = IntervalEntry {
        start: 0xdead_beef,
        end: 0x1234_5678,
        annotation: u32::MAX,
    };
    let mut entries = vec![guard; size.entries + 2];
    let mut offsets = vec![-77; size.entries + 3];
    let mut data = vec![0xcc; size.annotation_bytes + 2];
    table
        .write(IntervalsOut {
            entries: &mut entries[1..size.entries + 1],
            annotations: StringColumn {
                offsets: &mut offsets[1..size.entries + 2],
                data: &mut data[1..size.annotation_bytes + 1],
            },
        })
        .unwrap();
    assert_eq!(entries.first(), Some(&guard));
    assert_eq!(entries.last(), Some(&guard));
    assert_eq!(offsets.first(), Some(&-77));
    assert_eq!(offsets.last(), Some(&-77));
    assert_eq!(data.first(), Some(&0xcc));
    assert_eq!(data.last(), Some(&0xcc));
    assert_eq!(&entries[1..size.entries + 1], want_entries);
    assert_eq!(&offsets[1..size.entries + 2], want_offsets);
    assert_eq!(&data[1..size.annotation_bytes + 1], want_data);
}

fn check_region(table: &mut BedTable, records: &[Record], region: Region<'_>) {
    let size = table.query(region).unwrap();
    assert_output(table, records, region, size);
    let snapshot = table.clone();
    // Invalidating the original cannot invalidate its independently cloned stage.
    assert!(table.query(Region::new("c", 2, 1)).is_err());
    let mut offsets = [719];
    assert_eq!(
        table.write(IntervalsOut {
            entries: &mut [],
            annotations: StringColumn {
                offsets: &mut offsets,
                data: &mut [],
            },
        }),
        Err(Error::NotQueried)
    );
    assert_eq!(offsets, [719]);
    assert_output(&snapshot, records, region, size);
    assert_output(&snapshot, records, region, size);
}

#[test]
fn byte_annotations_survive_every_small_file_order() {
    let records = [
        Record {
            chrom: "c",
            start: 0,
            end: 9,
            note: vec![0xff, 0],
        },
        Record {
            chrom: "c",
            start: 2,
            end: 4,
            note: vec![],
        },
        Record {
            chrom: "c",
            start: 2,
            end: 2,
            note: vec![0xfe, b'\t', 0x80],
        },
        Record {
            chrom: "d",
            start: 1,
            end: 8,
            note: vec![b'd', 0],
        },
        Record {
            chrom: "c",
            start: 5,
            end: 7,
            note: vec![0x80],
        },
    ];
    fn permutations(prefix: &mut Vec<Record>, rest: Vec<Record>) {
        if rest.is_empty() {
            let mut table = BedTable::load(&render(prefix)).unwrap();
            for start in 0..=9 {
                for end in start..=9 {
                    check_region(&mut table, prefix, Region::new("c", start, end));
                }
            }
            check_region(&mut table, prefix, Region::new("d", 0, 9));
            check_region(&mut table, prefix, Region::new("absent", 0, u32::MAX));
            return;
        }
        for i in 0..rest.len() {
            let mut remaining = rest.clone();
            prefix.push(remaining.remove(i));
            permutations(prefix, remaining);
            prefix.pop();
        }
    }
    permutations(&mut Vec::new(), records.to_vec());
}

#[test]
fn monotone_ends_do_not_certify_reordered_or_interleaved_annotations() {
    for count in [4095, 4096, 4097] {
        for layout in 0..3 {
            let mut records = Vec::new();
            for i in 0..count as u32 {
                // Layout 0 is contiguous/sorted, 1 interleaves chromosomes,
                // and 2 starts reversed while all ends remain monotone.
                let start = if layout == 2 { count as u32 - i } else { i / 2 };
                let note = if i % 5 == 0 {
                    vec![]
                } else {
                    vec![0xff, (i % 251) as u8, b'\t', 0x80]
                };
                // Keep raw note bytes outside BED line delimiters.
                let note: Vec<_> = note
                    .into_iter()
                    .map(|b| if b == b'\n' || b == b'\r' { 0 } else { b })
                    .collect();
                records.push(Record {
                    chrom: "c",
                    start,
                    end: u32::MAX,
                    note,
                });
                if layout == 1 {
                    records.push(Record {
                        chrom: "d",
                        start: i,
                        end: i + 1,
                        note: vec![b'x'; (i % 17) as usize],
                    });
                }
            }
            let mut table = BedTable::load(&render(&records)).unwrap();
            for (start, end) in [
                (0, u32::MAX),
                (0, 1),
                (17, 301),
                (u32::MAX - 1, u32::MAX),
                (7, 7),
            ] {
                check_region(&mut table, &records, Region::new("c", start, end));
            }
        }
    }
}

#[test]
fn gzip_empty_members_and_binary_annotation_splits_match_raw_oracle() {
    fn member(bytes: &[u8]) -> Vec<u8> {
        let mut gzip = flate2::write::GzEncoder::new(Vec::new(), flate2::Compression::fast());
        gzip.write_all(bytes).unwrap();
        gzip.finish().unwrap()
    }
    let records = vec![
        Record {
            chrom: "c",
            start: 1,
            end: 3,
            note: vec![0xff, 0, b'\t', 0xfe],
        },
        Record {
            chrom: "c",
            start: 2,
            end: 9,
            note: vec![],
        },
    ];
    let mut text = render(&records);
    text.pop();
    for split in 0..=text.len() {
        let raw: Vec<_> = [
            member(&[]),
            member(&text[..split]),
            member(&[]),
            member(&text[split..]),
            member(&[]),
        ]
        .into_iter()
        .flatten()
        .collect();
        let mut table = BedTable::load(&raw).unwrap();
        check_region(&mut table, &records, Region::new("c", 0, 10));
    }
}

#[test]
fn rejected_output_shapes_leave_all_buffers_untouched_and_stage_reusable() {
    let records: Vec<_> = (0..4)
        .map(|i| Record {
            chrom: "c",
            start: i,
            end: i + 4,
            note: vec![b'x'; i as usize],
        })
        .collect();
    let region = Region::new("c", 0, 10);
    let mut table = BedTable::load(&render(&records)).unwrap();
    let size = table.query(region).unwrap();
    for (n, o, b) in [
        (3, 5, 6),
        (5, 5, 6),
        (4, 4, 6),
        (4, 6, 6),
        (4, 5, 5),
        (4, 5, 7),
    ] {
        let sentinel = IntervalEntry {
            start: 11,
            end: 22,
            annotation: 33,
        };
        let mut entries = vec![sentinel; n];
        let mut offsets = vec![-3; o];
        let mut data = vec![0xff; b];
        assert!(matches!(
            table.write(IntervalsOut {
                entries: &mut entries,
                annotations: StringColumn {
                    offsets: &mut offsets,
                    data: &mut data
                },
            }),
            Err(Error::OutputSize { .. })
        ));
        assert_eq!(entries, vec![sentinel; n]);
        assert_eq!(offsets, vec![-3; o]);
        assert_eq!(data, vec![0xff; b]);
        assert_output(&table, &records, region, size);
    }
}

#[test]
fn consecutive_empty_annotations_need_offsets_but_no_bytes() {
    let records: Vec<_> = (0..8)
        .map(|i| Record {
            chrom: "c",
            start: i,
            end: i + 1,
            note: vec![],
        })
        .collect();
    let mut table = BedTable::load(&render(&records)).unwrap();
    for region in [
        Region::new("c", 0, 8),
        Region::new("c", 2, 6),
        Region::new("c", 3, 4),
        Region::new("c", 8, 9),
    ] {
        check_region(&mut table, &records, region);
    }
}
