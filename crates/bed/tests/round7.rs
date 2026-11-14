//! File-row oracle for single and multiple chromosome grouping at load time.
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

fn verify(table: &mut BedTable, rows: &[Row], region: Region<'_>) {
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
        start: 43,
        end: 57,
        annotation: 91,
    };
    let mut entries = vec![guard; size.entries + 2];
    let mut offsets = vec![-991; size.entries + 3];
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
        assert_eq!((entries[0], entries[size.entries + 1]), (guard, guard));
        assert_eq!((offsets[0], offsets[size.entries + 2]), (-991, -991));
        assert_eq!((data[0], data[size.annotation_bytes + 1]), (0xa5, 0xa5));
    }
}

#[test]
fn grouping_preserves_file_row_annotations_and_stable_starts() {
    for n in [0u32, 1, 2, 4095, 4096, 4097] {
        for order in ["sorted", "reverse", "permuted"] {
            for multiple in [false, true] {
                let mut rows: Vec<_> = (0..n)
                    .map(|i| {
                        let start = if order == "permuted" {
                            i.wrapping_mul(0x9e37_79b9) % 32_768
                        } else {
                            i / 3 * 4
                        };
                        Row {
                            chrom: "c",
                            start,
                            end: if i % 17 == 0 { u32::MAX } else { start + i % 7 },
                            note: if i % 4 == 0 {
                                vec![]
                            } else {
                                vec![0xff, 0, b'\t', 0x80 + (i % 31) as u8]
                            },
                        }
                    })
                    .collect();
                if order == "reverse" {
                    rows.reverse();
                }
                if multiple {
                    rows = rows
                        .into_iter()
                        .flat_map(|row| {
                            let other = Row {
                                chrom: "d",
                                start: row.start,
                                end: row.end,
                                note: vec![b'd'],
                            };
                            [row, other]
                        })
                        .collect();
                }
                let mut text = Vec::new();
                for row in &rows {
                    write!(text, "{}\t{}\t{}\t", row.chrom, row.start, row.end).unwrap();
                    text.extend_from_slice(&row.note);
                    text.push(b'\n');
                }
                let mut encoder =
                    flate2::write::GzEncoder::new(Vec::new(), flate2::Compression::fast());
                encoder.write_all(&text).unwrap();
                for bytes in [&text, &encoder.finish().unwrap()] {
                    let mut table = BedTable::load(bytes).unwrap();
                    assert_eq!(table.len(), rows.len());
                    assert_eq!(table.is_empty(), rows.is_empty());
                    for (i, row) in rows.iter().enumerate() {
                        assert_eq!(table.annotations().get(i), row.note);
                    }
                    for region in [
                        Region::new("c", 0, u32::MAX),
                        Region::new("c", 0, 1),
                        Region::new("c", 4096, 8192),
                        Region::new("c", u32::MAX - 1, u32::MAX),
                        Region::new("c", 10, 10),
                        Region::new("d", 0, u32::MAX),
                        Region::new("missing", 0, u32::MAX),
                    ] {
                        verify(&mut table, &rows, region);
                        verify(&mut table.clone(), &rows, region);
                    }
                    assert!(table.query(Region::new("c", 2, 1)).is_err());
                    let mut offsets = [123];
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
                    assert_eq!(offsets, [123]);
                    verify(&mut table, &rows, Region::new("c", 0, u32::MAX));
                }
            }
        }
    }
}
