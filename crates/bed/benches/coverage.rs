//! Complete nonmonotone spans, with fixture and output verification outside timing.
use omicsweave_bed::BedTable;
use omicsweave_core::{IntervalEntry, IntervalRead, IntervalsOut, Region, StringColumn};
use std::fmt::Write;
use std::hint::black_box;
use std::time::{Duration, Instant};

fn measure(name: &str, mut call: impl FnMut()) {
    if std::env::var("BED_BENCH_FILTER").is_ok_and(|s| !name.contains(&s)) {
        return;
    }
    call();
    let mut samples = [0.0; 5];
    for sample in &mut samples {
        let start = Instant::now();
        let mut count = 0;
        while start.elapsed() < Duration::from_millis(100) {
            for _ in 0..8 {
                call();
            }
            count += 8;
        }
        *sample = start.elapsed().as_nanos() as f64 / count as f64;
    }
    samples.sort_by(f64::total_cmp);
    println!("{name}\t{:.2}", samples[2]);
}

fn fixture(n: u32, layout: &str) -> (Vec<u8>, Vec<IntervalEntry>, Vec<i64>, Vec<u8>) {
    let mut rows: Vec<_> = (0..n)
        .map(|i| (i * 8, i * 8 + 3, format!("{i:08x}")))
        .collect();
    if layout == "reordered" {
        rows.reverse();
    }
    rows.insert(0, (0, u32::MAX, "enclosing".into()));
    let mut text = String::new();
    for (i, (start, end, note)) in rows.iter().enumerate() {
        writeln!(text, "c\t{start}\t{end}\t{note}").unwrap();
        if layout == "interleaved" {
            writeln!(text, "d\t{i}\t{}\tother", i + 1).unwrap();
        }
    }
    // Stable sorting preserves ties in their explicit input row order.
    rows.sort_by_key(|r| r.0);
    let mut expected = Vec::new();
    let mut offsets = vec![0];
    let mut data = Vec::new();
    for (start, end, note) in rows {
        expected.push(IntervalEntry {
            start,
            end,
            annotation: expected.len() as u32,
        });
        data.extend_from_slice(note.as_bytes());
        offsets.push(data.len() as i64);
    }
    (text.into_bytes(), expected, offsets, data)
}

fn main() {
    for n in [1024u32, 131_072] {
        for layout in ["contiguous", "interleaved", "reordered"] {
            let name = format!("bed.coverage.n{n}.{layout}");
            let (text, expected, expected_offsets, expected_data) = fixture(n, layout);
            let mut table = BedTable::load(&text).unwrap();
            let full = Region::new("c", 0, n * 8);
            let partial = Region::new("c", n * 6, n * 6 + 1);
            let size = table.query(full).unwrap();
            let mut entries = vec![IntervalEntry::default(); size.entries];
            let mut offsets = vec![0; size.entries + 1];
            let mut data = vec![0; size.annotation_bytes];
            table
                .write(IntervalsOut {
                    entries: &mut entries,
                    annotations: StringColumn {
                        offsets: &mut offsets,
                        data: &mut data,
                    },
                })
                .unwrap();
            assert_eq!(entries, expected);
            assert_eq!(offsets, expected_offsets);
            assert_eq!(data, expected_data);
            println!(
                "{name}.verified\t{}:{}:{}",
                entries.len(),
                offsets.last().unwrap(),
                data.iter().map(|&b| b as u64).sum::<u64>()
            );
            measure(&format!("{name}.load"), || {
                black_box(BedTable::load(black_box(&text)).unwrap());
            });
            measure(&format!("{name}.query.full"), || {
                black_box(table.query(black_box(full)).unwrap());
            });
            measure(&format!("{name}.write.full"), || {
                black_box(table.query(black_box(full)).unwrap());
                table
                    .write(IntervalsOut {
                        entries: black_box(&mut entries),
                        annotations: StringColumn {
                            offsets: black_box(&mut offsets),
                            data: black_box(&mut data),
                        },
                    })
                    .unwrap();
                black_box((&entries, &offsets, &data));
            });
            measure(&format!("{name}.query.partial"), || {
                black_box(table.query(black_box(partial)).unwrap());
            });
        }
    }
    block_cases();
    grouping_cases();
}

fn grouping_cases() {
    for n in [1024u32, 131_072] {
        for chroms in [1u32, 4] {
            for order in ["sorted", "reverse", "permuted"] {
                let name = format!("bed.grouping.n{n}.c{chroms}.{order}.load");
                if std::env::var("BED_BENCH_FILTER").is_ok_and(|s| !name.contains(&s)) {
                    continue;
                }
                let mut rows: Vec<_> = (0..n)
                    .map(|i| {
                        let start = if order == "permuted" {
                            i.wrapping_mul(0x9e37_79b9) % (n * 8)
                        } else {
                            i / 3 * 8
                        };
                        (
                            i % chroms,
                            start,
                            start + i % 13,
                            "x".repeat(i as usize % 23),
                        )
                    })
                    .collect();
                if order == "reverse" {
                    rows.reverse();
                }
                let mut text = String::new();
                for (chrom, start, end, note) in &rows {
                    writeln!(text, "c{chrom}\t{start}\t{end}\t{note}").unwrap();
                }
                let table = BedTable::load(text.as_bytes()).unwrap();
                assert_eq!(table.len(), rows.len());
                let mut checksum = 0u64;
                for chrom in 0..chroms {
                    let mut expected: Vec<_> = rows
                        .iter()
                        .enumerate()
                        .filter(|(_, row)| row.0 == chrom)
                        .collect();
                    expected.sort_unstable_by_key(|(file_row, row)| (row.1, *file_row));
                    let entries = table.entries(&format!("c{chrom}"));
                    assert_eq!(entries.len(), expected.len());
                    for (entry, (file_row, row)) in entries.iter().zip(expected) {
                        assert_eq!(
                            *entry,
                            IntervalEntry {
                                start: row.1,
                                end: row.2,
                                annotation: file_row as u32,
                            }
                        );
                        assert_eq!(table.annotations().get(file_row), row.3.as_bytes());
                        checksum = checksum.wrapping_add(
                            u64::from(entry.start)
                                + u64::from(entry.end)
                                + u64::from(entry.annotation)
                                + row.3.bytes().map(u64::from).sum::<u64>(),
                        );
                    }
                }
                println!("{name}.checksum\t{checksum}");
                measure(&name, || {
                    black_box(BedTable::load(black_box(text.as_bytes())).unwrap());
                });
            }
        }
    }
}

fn block_cases() {
    for n in [1024u32, 131_072] {
        for layout in ["contiguous", "interleaved", "reordered"] {
            let name = format!("bed.blocks.n{n}.{layout}");
            let mut rows: Vec<_> = (0..n)
                .map(|i| (i / 2 * 4, i / 2 * 4 + 3, vec![b'x'; i as usize % 23]))
                .collect();
            if layout == "reordered" {
                rows.reverse();
            }
            let mut text = Vec::new();
            for (i, (start, end, note)) in rows.iter().enumerate() {
                std::io::Write::write_fmt(&mut text, format_args!("c\t{start}\t{end}\t")).unwrap();
                text.extend_from_slice(note);
                text.push(b'\n');
                if layout == "interleaved" {
                    std::io::Write::write_fmt(
                        &mut text,
                        format_args!("d\t{i}\t{}\tother\n", i + 1),
                    )
                    .unwrap();
                }
            }
            rows.sort_by_key(|row| row.0);
            let mut table = BedTable::load(&text).unwrap();
            measure(&format!("{name}.load"), || {
                black_box(BedTable::load(black_box(&text)).unwrap());
            });
            for (width, region) in [
                ("narrow", Region::new("c", n, n + 1)),
                ("fringe", Region::new("c", 3, n * 2 - 3)),
                ("wide", Region::new("c", 0, n * 2)),
            ] {
                let mut expected_entries = Vec::new();
                let mut expected_offsets = vec![0];
                let mut expected_data = Vec::new();
                for (start, end, note) in &rows {
                    if *start < region.end && *end > region.start {
                        expected_entries.push(IntervalEntry {
                            start: *start,
                            end: *end,
                            annotation: expected_entries.len() as u32,
                        });
                        expected_data.extend_from_slice(note);
                        expected_offsets.push(expected_data.len() as i64);
                    }
                }
                let size = table.query(region).unwrap();
                assert_eq!(size.entries, expected_entries.len());
                assert_eq!(size.annotation_bytes, expected_data.len());
                let mut entries = vec![IntervalEntry::default(); size.entries];
                let mut offsets = vec![0; size.entries + 1];
                let mut data = vec![0; size.annotation_bytes];
                table
                    .write(IntervalsOut {
                        entries: &mut entries,
                        annotations: StringColumn {
                            offsets: &mut offsets,
                            data: &mut data,
                        },
                    })
                    .unwrap();
                assert_eq!(entries, expected_entries);
                assert_eq!(offsets, expected_offsets);
                assert_eq!(data, expected_data);
                println!(
                    "{name}.{width}.verified\t{}:{}:{}",
                    entries.len(),
                    data.len(),
                    data.iter().map(|&b| b as u64).sum::<u64>()
                );
                measure(&format!("{name}.query.{width}"), || {
                    black_box(table.query(black_box(region)).unwrap());
                });
                measure(&format!("{name}.write.{width}"), || {
                    black_box(table.query(black_box(region)).unwrap());
                    table
                        .write(IntervalsOut {
                            entries: black_box(&mut entries),
                            annotations: StringColumn {
                                offsets: black_box(&mut offsets),
                                data: black_box(&mut data),
                            },
                        })
                        .unwrap();
                    black_box((&entries, &offsets, &data));
                });
            }
        }
    }
}
