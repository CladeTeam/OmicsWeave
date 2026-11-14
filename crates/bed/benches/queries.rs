//! Synthetic in-memory queries and whole-file loading, with exact output checks.
use omicsweave_bed::BedTable;
use omicsweave_core::{IntervalEntry, IntervalRead, IntervalsOut, Region, StringColumn};
use std::{
    fmt::Write,
    hint::black_box,
    time::{Duration, Instant},
};

fn selected(name: &str) -> bool {
    !std::env::var("BED_BENCH_FILTER").is_ok_and(|filter| !name.contains(&filter))
}

fn measure(name: &str, batch: usize, mut run: impl FnMut()) {
    if !selected(name) {
        return;
    }
    run();
    let mut samples = [0.0f64; 5];
    for sample in &mut samples {
        let start = Instant::now();
        let mut count = 0;
        while start.elapsed() < Duration::from_millis(100) {
            for _ in 0..batch {
                run();
            }
            count += batch;
        }
        *sample = start.elapsed().as_nanos() as f64 / count as f64;
    }
    samples.sort_by(f64::total_cmp);
    println!("{name}\t{:.2}", samples[2]);
}

fn queries(name: &str, text: &[u8], region: Region<'_>, expected: &[IntervalEntry]) {
    queries_expected(
        name,
        text,
        region,
        expected,
        &(0..=expected.len() as i64).collect::<Vec<_>>(),
        &vec![b'x'; expected.len()],
    );
}

fn queries_expected(
    name: &str,
    text: &[u8],
    region: Region<'_>,
    expected: &[IntervalEntry],
    expected_offsets: &[i64],
    expected_data: &[u8],
) {
    let mut table = BedTable::load(text).unwrap();
    let size = table.query(region).unwrap();
    assert_eq!(size.entries, expected.len(), "{name}");
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
    assert_eq!(entries, expected, "{name}");
    assert_eq!(offsets, expected_offsets, "{name}");
    assert_eq!(data, expected_data, "{name}");
    measure(name, 32, || {
        black_box(table.query(black_box(region)).unwrap());
    });
    measure(&format!("{name}.write"), 8, || {
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
        black_box(&entries);
        black_box(&offsets);
        black_box(&data);
    });
}

fn main() {
    println!("case\tns/op (median of 5 x 100ms, warm buffers; loads include drop)");
    for sorted in [false, true] {
        let mut starts: Vec<_> = (0..131_072u32)
            .map(|i| {
                let mut x = i.wrapping_mul(0x9E37_79B9);
                x ^= x >> 16;
                x % 65_536 * 4
            })
            .collect();
        if sorted {
            starts.sort_unstable();
        }
        let mut text = String::new();
        for start in starts {
            writeln!(text, "s\t{start}\t{}\tx", start + 2).unwrap();
        }
        let name = format!("bed.load.sorted{sorted}");
        measure(&name, 1, || {
            black_box(BedTable::load(black_box(text.as_bytes())).unwrap());
        });
        let mut encoder = flate2::write::GzEncoder::new(Vec::new(), flate2::Compression::new(6));
        std::io::Write::write_all(&mut encoder, text.as_bytes()).unwrap();
        let gzip = encoder.finish().unwrap();
        assert_eq!(BedTable::load(&gzip).unwrap().len(), 131_072);
        measure(&format!("{name}.gzip"), 1, || {
            black_box(BedTable::load(black_box(&gzip)).unwrap());
        });
    }
    for n in [1_024u32, 131_072] {
        for nested in [false, true] {
            let mut text = String::new();
            let mut records = Vec::new();
            if nested {
                writeln!(text, "chr1\t0\t{}\tx", n * 10 + 100).unwrap();
                records.push((0, n * 10 + 100));
            }
            for i in 0..n {
                writeln!(text, "chr1\t{}\t{}\tx", i * 10, i * 10 + 5).unwrap();
                records.push((i * 10, i * 10 + 5));
            }
            for wide in [false, true] {
                let name = format!("bed.n{n}.nested{nested}.wide{wide}");
                let region = if wide {
                    Region::new("chr1", 0, n * 10)
                } else {
                    Region::new("chr1", n * 9, n * 9 + 10)
                };
                queries(&name, text.as_bytes(), region, &scan(&records, region));
            }
        }
    }
    for layout in ["interleaved", "duplicates", "sparse", "irregular"] {
        let n = 131_072u32;
        let mut text = String::new();
        let mut records = Vec::new();
        for i in 0..n {
            let chrom = if layout == "interleaved" { i % 4 } else { 0 };
            let start = if layout == "duplicates" {
                i / 16 * 10
            } else {
                i * 10
            };
            let end = start
                + if layout == "irregular" {
                    (i % 19) * 11
                } else {
                    5
                };
            writeln!(text, "chr{chrom}\t{start}\t{end}\tx").unwrap();
            if chrom == 0 {
                records.push((start, end));
            }
        }
        for wide in [false, true] {
            let start = if wide {
                0
            } else if layout == "duplicates" {
                n / 16 * 9
            } else {
                n * 9
            };
            let start = start + if !wide && layout == "sparse" { 6 } else { 0 };
            let end = if wide {
                n * 10
            } else {
                start + if layout == "sparse" { 1 } else { 10 }
            };
            let region = Region::new("chr0", start, end);
            queries(
                &format!("bed.{layout}.wide{wide}"),
                text.as_bytes(),
                region,
                &scan(&records, region),
            );
        }
    }
    for annotation_len in [0usize, 32, 256] {
        let annotation = "x".repeat(annotation_len);
        let n = 131_072u32;
        for nested in [false, true] {
            let mut text = String::new();
            let mut records = Vec::new();
            if nested {
                writeln!(text, "chr1\t0\t{}\t{annotation}", n * 10 + 100).unwrap();
                records.push((0, n * 10 + 100));
            }
            for i in 0..n {
                writeln!(text, "chr1\t{}\t{}\t{annotation}", i * 10, i * 10 + 5).unwrap();
                records.push((i * 10, i * 10 + 5));
            }
            for wide in [false, true] {
                let name = format!("bed.annotation{annotation_len}.nested{nested}.wide{wide}");
                if !selected(&name) && !selected(&format!("{name}.write")) {
                    continue;
                }
                let region = if wide {
                    Region::new("chr1", 0, n * 10)
                } else {
                    Region::new("chr1", n * 9, n * 9 + 10)
                };
                let expected = scan(&records, region);
                let offsets: Vec<_> = (0..=expected.len())
                    .map(|row| (row * annotation_len) as i64)
                    .collect();
                let data = vec![b'x'; expected.len() * annotation_len];
                queries_expected(&name, text.as_bytes(), region, &expected, &offsets, &data);
            }
        }
    }
}

fn scan(records: &[(u32, u32)], region: Region<'_>) -> Vec<IntervalEntry> {
    records
        .iter()
        .filter(|&&(start, end)| region.overlaps(start, end))
        .enumerate()
        .map(|(row, &(start, end))| IntervalEntry {
            start,
            end,
            annotation: row as u32,
        })
        .collect()
}
