//! Complete and partial coverage against an independent stable file-row oracle.
use std::io::Write;

use omicsweave_bed::BedTable;
use omicsweave_core::{Error, IntervalEntry, IntervalRead, IntervalsOut, Region, StringColumn};

#[derive(Clone)]
struct Row {
    chrom: &'static str,
    start: u32,
    end: u32,
    note: Vec<u8>,
}

fn render(rows: &[Row]) -> Vec<u8> {
    let mut bytes = Vec::new();
    for row in rows {
        write!(bytes, "{}\t{}\t{}\t", row.chrom, row.start, row.end).unwrap();
        bytes.extend_from_slice(&row.note);
        bytes.push(b'\n');
    }
    bytes
}

fn check(table: &mut BedTable, rows: &[Row], region: Region<'_>) {
    let mut selected: Vec<_> = rows
        .iter()
        .enumerate()
        .filter(|(_, row)| {
            region.start < region.end
                && row.chrom == region.chrom
                && row.start < region.end
                && row.end > region.start
        })
        .collect();
    selected.sort_unstable_by_key(|(file_row, row)| (row.start, *file_row));
    let mut expected_entries = Vec::new();
    let mut expected_offsets = vec![0];
    let mut expected_data = Vec::new();
    for (_, row) in selected {
        expected_entries.push(IntervalEntry {
            start: row.start,
            end: row.end,
            annotation: expected_entries.len() as u32,
        });
        expected_data.extend_from_slice(&row.note);
        expected_offsets.push(expected_data.len() as i64);
    }
    let size = table.query(region).unwrap();
    assert_eq!(size.entries, expected_entries.len());
    assert_eq!(size.annotation_bytes, expected_data.len());
    let guard = IntervalEntry {
        start: 765,
        end: 876,
        annotation: 987,
    };
    let mut entries = vec![guard; size.entries + 2];
    let mut offsets = vec![-991; size.entries + 3];
    let mut data = vec![0xab; size.annotation_bytes + 2];
    for _ in 0..2 {
        table
            .write(IntervalsOut {
                entries: &mut entries[1..size.entries + 1],
                annotations: StringColumn {
                    offsets: &mut offsets[1..size.entries + 2],
                    data: &mut data[1..size.annotation_bytes + 1],
                },
            })
            .unwrap();
        assert_eq!(&entries[1..size.entries + 1], expected_entries);
        assert_eq!(&offsets[1..size.entries + 2], expected_offsets);
        assert_eq!(&data[1..size.annotation_bytes + 1], expected_data);
        assert_eq!((entries[0], entries[size.entries + 1]), (guard, guard));
        assert_eq!((offsets[0], offsets[size.entries + 2]), (-991, -991));
        assert_eq!((data[0], data[size.annotation_bytes + 1]), (0xab, 0xab));
    }
}

#[test]
fn complete_and_partial_coverage_preserve_sorted_rows_at_all_boundaries() {
    for layout in ["contiguous", "interleaved", "reordered"] {
        for zero_end in [false, true] {
            let mut rows = vec![Row {
                chrom: "c",
                start: 0,
                end: u32::MAX,
                note: vec![0xff, 0, b'\t'],
            }];
            for i in 0..33 {
                rows.push(Row {
                    chrom: "c",
                    start: i * 3,
                    end: if zero_end && i == 0 { 0 } else { i * 3 + 5 },
                    note: if i % 3 == 0 {
                        vec![]
                    } else {
                        vec![b'x', i as u8 + 128]
                    },
                });
                if layout == "interleaved" {
                    rows.push(Row {
                        chrom: "d",
                        start: i,
                        end: i + 1,
                        note: vec![b'd'],
                    });
                }
            }
            if layout == "reordered" {
                rows.reverse();
            }
            let mut table = BedTable::load(&render(&rows)).unwrap();
            let boundaries = [0, 1, 4, 5, 6, 32, 95, 96, 97, 101, u32::MAX - 1, u32::MAX];
            for &start in &boundaries {
                for &end in &boundaries {
                    if start <= end {
                        check(&mut table, &rows, Region::new("c", start, end));
                    }
                }
            }
            check(&mut table, &rows, Region::new("d", 0, u32::MAX));
            check(&mut table, &rows, Region::new("absent", 0, u32::MAX));
        }
    }
}

#[test]
fn large_reordered_complete_spans_keep_annotation_pairing_and_recover() {
    let mut rows = vec![Row {
        chrom: "c",
        start: 0,
        end: u32::MAX,
        note: b"enclosing".to_vec(),
    }];
    for i in 0..4097u32 {
        let start = i.wrapping_mul(0x9e37_79b9) % 16_384;
        rows.push(Row {
            chrom: "c",
            start,
            end: start + 7,
            note: if i % 4 == 0 {
                vec![]
            } else {
                format!("row-{i}").into_bytes()
            },
        });
    }
    let mut table = BedTable::load(&render(&rows)).unwrap();
    for region in [
        Region::new("c", 0, u32::MAX),
        Region::new("c", 7, 16_384),
        Region::new("c", 0, u32::MAX),
        Region::new("c", 9_999, 10_000),
    ] {
        check(&mut table, &rows, region);
        let mut clone = table.clone();
        assert!(table.query(Region::new("c", 2, 1)).is_err());
        let mut offsets = [991];
        assert_eq!(
            table.write(IntervalsOut {
                entries: &mut [],
                annotations: StringColumn {
                    offsets: &mut offsets,
                    data: &mut []
                },
            }),
            Err(Error::NotQueried)
        );
        assert_eq!(offsets, [991]);
        check(&mut clone, &rows, region);
        check(&mut table, &rows, region);
    }
}
