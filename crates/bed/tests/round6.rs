//! Independent file-row oracle across annotation block boundaries.
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

fn fixture(n: u32, layout: &str, zero_length: bool) -> (Vec<Row>, Vec<u8>) {
    let mut rows = Vec::new();
    for i in 0..n {
        let start = i / 2 * 4;
        let mut note = Vec::new();
        if i % 5 != 0 {
            note.resize(i as usize % 23, 0xff);
            note.extend_from_slice(&[0, b'\t', 0x80]);
        }
        rows.push(Row {
            chrom: "c",
            start,
            end: start + if zero_length { 0 } else { 3 },
            note,
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
    let mut text = Vec::new();
    for row in &rows {
        write!(text, "{}\t{}\t{}\t", row.chrom, row.start, row.end).unwrap();
        text.extend_from_slice(&row.note);
        text.push(b'\n');
    }
    (rows, text)
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
    let sentinel = IntervalEntry {
        start: 123,
        end: 456,
        annotation: 789,
    };
    let mut entries = vec![sentinel; size.entries + 2];
    let mut offsets = vec![-97; size.entries + 3];
    let mut data = vec![0xa5; size.annotation_bytes + 2];
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
        assert_eq!(
            (entries[0], entries[size.entries + 1]),
            (sentinel, sentinel)
        );
        assert_eq!((offsets[0], offsets[size.entries + 2]), (-97, -97));
        assert_eq!((data[0], data[size.annotation_bytes + 1]), (0xa5, 0xa5));
    }
    let mut bad_offsets = [93];
    assert!(table.query(Region::new("c", 2, 1)).is_err());
    assert_eq!(
        table.write(IntervalsOut {
            entries: &mut [],
            annotations: StringColumn {
                offsets: &mut bad_offsets,
                data: &mut [],
            },
        }),
        Err(Error::NotQueried)
    );
    assert_eq!(bad_offsets, [93]);
}

#[test]
fn monotone_queries_match_file_rows_across_every_block_edge() {
    for n in [63, 64, 65, 127, 128, 129, 1025, 4097] {
        for layout in ["contiguous", "interleaved", "reordered"] {
            for zero_length in [false, true] {
                let (rows, text) = fixture(n, layout, zero_length);
                let mut table = BedTable::load(&text).unwrap();
                assert_eq!(table.chroms()[0], "c");
                let mut boundaries = vec![0, 1, 2, 3, 4, n * 2, n * 2 + 1];
                for boundary in [64u32, 128, 192, 1024, 4096] {
                    if boundary < n {
                        boundaries.extend([boundary * 2 - 1, boundary * 2, boundary * 2 + 1]);
                    }
                }
                for &start in &boundaries {
                    for &end in &boundaries {
                        if start <= end {
                            check(&mut table, &rows, Region::new("c", start, end));
                        }
                    }
                }
                check(&mut table, &rows, Region::new("d", 0, u32::MAX));
                check(&mut table, &rows, Region::new("missing", 0, u32::MAX));
                check(&mut table.clone(), &rows, Region::new("c", 1, n * 2));
                let mut encoder =
                    flate2::write::GzEncoder::new(Vec::new(), flate2::Compression::fast());
                encoder.write_all(&text).unwrap();
                let mut gzip = BedTable::load(&encoder.finish().unwrap()).unwrap();
                check(&mut gzip, &rows, Region::new("c", 1, n * 2));
            }
        }
    }
}
