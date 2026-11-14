//! Standard reader queries and staged writes, with source I/O counted separately.
#[path = "../tests/common/mod.rs"]
mod common;
use common::bigbed::{Record, Spec, write};
use omicsweave_bbi::BbiReader;
use omicsweave_core::{
    IntervalEntry, IntervalRead, IntervalsOut, ReadAt, Region, Result, StringColumn, TrackRead,
};
use std::{
    cell::Cell,
    hint::black_box,
    rc::Rc,
    time::{Duration, Instant},
};

#[derive(Clone)]
struct Counted<'a> {
    bytes: &'a [u8],
    calls: Rc<Cell<usize>>,
    read: Rc<Cell<usize>>,
}
impl ReadAt for Counted<'_> {
    fn len(&self) -> u64 {
        self.bytes.len() as u64
    }
    fn read_exact_at(&self, out: &mut [u8], at: u64) -> Result<()> {
        self.calls.set(self.calls.get() + 1);
        self.read.set(self.read.get() + out.len());
        self.bytes.read_exact_at(out, at)
    }
}
fn measure(name: &str, mut f: impl FnMut()) {
    for _ in 0..20 {
        f();
    }
    let mut samples = [0.0f64; 5];
    for sample in &mut samples {
        let t = Instant::now();
        let mut n = 0;
        while t.elapsed() < Duration::from_millis(100) {
            for _ in 0..16 {
                f();
            }
            n += 16;
        }
        *sample = t.elapsed().as_nanos() as f64 / n as f64;
    }
    samples.sort_by(f64::total_cmp);
    println!("{name}\t{:.2}", samples[2]);
}
fn main() {
    let large_spec = Spec {
        chroms: vec![("s".into(), 100_000)],
        records: (0..64)
            .map(|i| Record {
                chrom: 0,
                start: i * 16,
                end: i * 16 + 8,
                rest: vec![b'a'; 10_000],
            })
            .collect(),
        items_per_block: 8,
        split_chroms: true,
        chrom_fanout: 2,
        index_fanout: 16,
        auto_sql: None,
        compress: true,
        field_count: 4,
        defined_field_count: 3,
    };
    let large_bytes = write(&large_spec);
    let large_region = Region::new("s", 0, 100_000);
    let mut large_reader = BbiReader::open(&large_bytes[..]).unwrap();
    let size = large_reader.query(large_region).unwrap();
    assert_eq!(size.entries, 64);
    assert_eq!(size.annotation_bytes, 640_000);
    let mut rows = vec![IntervalEntry::default(); size.entries];
    let mut offsets = vec![0; size.entries + 1];
    let mut notes = vec![0; size.annotation_bytes];
    large_reader
        .write(IntervalsOut {
            entries: &mut rows,
            annotations: StringColumn {
                offsets: &mut offsets,
                data: &mut notes,
            },
        })
        .unwrap();
    for (i, row) in rows.iter().enumerate() {
        assert_eq!(
            (row.start, row.end, row.annotation),
            (i as u32 * 16, i as u32 * 16 + 8, i as u32)
        );
        assert_eq!(
            (offsets[i], offsets[i + 1]),
            (i as i64 * 10_000, (i as i64 + 1) * 10_000)
        );
    }
    assert!(notes.iter().all(|&b| b == b'a'));
    measure("bbi.interval.largeblocks.open_query.cold_buffers", || {
        let mut reader = BbiReader::open(black_box(&large_bytes[..])).unwrap();
        black_box(reader.query(black_box(large_region)).unwrap());
    });
    measure("bbi.interval.largeblocks.query.warm_buffers", || {
        black_box(large_reader.query(black_box(large_region)).unwrap());
    });
    for kind in ["disjoint", "monotone", "nested"] {
        let dir = tempfile::tempdir().unwrap();
        let path = dir.path().join("track.bw");
        let mut writer = common::libbigwig::Writer::create(&path, 0, &[("s", 100_000)]).unwrap();
        let records: Vec<_> = (0..4096)
            .map(|i| {
                let (s, e) = match kind {
                    "disjoint" => (i * 16, i * 16 + 8),
                    "monotone" => (i, 10_000 + i),
                    _ => (i, if i % 2 == 0 { 20_000 + i } else { 10_000 + i }),
                };
                (s, e, i as f32)
            })
            .collect();
        assert!(writer.add_bedgraph("s", &records));
        drop(writer);
        let bytes = std::fs::read(path).unwrap();
        let mut reader = BbiReader::open(bytes).unwrap();
        for (name, start, end) in [
            ("long", 5000, 9000),
            ("narrow", 5000, 5064),
            ("miss", 80000, 80064),
        ] {
            let region = Region::new("s", start, end);
            let mut values = vec![0.0; region.len() as usize];
            reader.values(region, &mut values).unwrap();
            let mut expected = vec![f32::NAN; values.len()];
            for &(s, e, value) in &records {
                let (lo, hi) = (s.max(start), e.min(end));
                if lo < hi {
                    expected[(lo - start) as usize..(hi - start) as usize].fill(value);
                }
            }
            assert!(
                values
                    .iter()
                    .zip(expected)
                    .all(|(a, b)| a.to_bits() == b.to_bits())
            );
            measure(&format!("bbi.values.{kind}.{name}"), || {
                reader
                    .values(black_box(region), black_box(&mut values))
                    .unwrap();
                black_box(&values);
            });
        }
    }
    for compress in [false, true] {
        let spec = Spec {
            chroms: vec![("s".into(), 100_000)],
            records: (0..4096)
                .map(|i| Record {
                    chrom: 0,
                    start: i * 16,
                    end: i * 16 + 8,
                    rest: vec![b'a'; 128],
                })
                .collect(),
            items_per_block: 16,
            split_chroms: true,
            chrom_fanout: 2,
            index_fanout: 16,
            auto_sql: None,
            compress,
            field_count: 4,
            defined_field_count: 3,
        };
        let bytes = write(&spec);
        let source = Counted {
            bytes: &bytes,
            calls: Rc::new(Cell::new(0)),
            read: Rc::new(Cell::new(0)),
        };
        let mut reader = BbiReader::open(source.clone()).unwrap();
        let region = Region::new("s", 0, 100_000);
        source.calls.set(0);
        source.read.set(0);
        let size = reader.query(region).unwrap();
        println!(
            "io\tbbi.interval.compressed{compress}\t{}\t{}",
            source.calls.get(),
            source.read.get()
        );
        let mut rows = vec![IntervalEntry::default(); size.entries];
        let mut offsets = vec![0; size.entries + 1];
        let mut data = vec![0; size.annotation_bytes];
        measure(&format!("bbi.interval.compressed{compress}.query"), || {
            black_box(reader.query(black_box(region)).unwrap());
        });
        measure(&format!("bbi.interval.compressed{compress}.write"), || {
            reader
                .write(IntervalsOut {
                    entries: black_box(&mut rows),
                    annotations: StringColumn {
                        offsets: black_box(&mut offsets),
                        data: black_box(&mut data),
                    },
                })
                .unwrap();
        });
        assert_eq!(rows.len(), 4096);
        assert_eq!(offsets[4096], size.annotation_bytes as i64);
    }
}
