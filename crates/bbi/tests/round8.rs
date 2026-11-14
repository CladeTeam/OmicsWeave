//! autoSql strings are bounded by their terminator, independently of padding.
use std::cell::Cell;
use std::hint::black_box;
use std::io::Cursor;
use std::time::{Duration, Instant};

use omicsweave_bbi::{BbiReader, BedSchema, BigBedWriter, Chrom, WriteOptions};
use omicsweave_core::sidecar::checksum;
use omicsweave_core::{ReadAt, Result};

fn fixture(schema: &[u8], padding: usize) -> Vec<u8> {
    let writer = BigBedWriter::create(
        Cursor::new(Vec::new()),
        vec![Chrom {
            name: "c".into(),
            len: 100,
        }],
        BedSchema::bed(3),
        WriteOptions {
            compress: false,
            zoom_reductions: Some(vec![]),
            ..WriteOptions::default()
        },
    )
    .unwrap();
    let mut bytes = writer.finish().unwrap().0.into_inner();
    let offset = bytes.len() as u64;
    bytes[36..44].copy_from_slice(&offset.to_le_bytes());
    bytes.extend_from_slice(schema);
    bytes.resize(bytes.len() + padding, b'p');
    bytes
}

struct CountReads<'a> {
    bytes: &'a [u8],
    read: Cell<usize>,
    largest: Cell<usize>,
}

impl ReadAt for CountReads<'_> {
    fn len(&self) -> u64 {
        self.bytes.len() as u64
    }

    fn read_exact_at(&self, out: &mut [u8], at: u64) -> Result<()> {
        self.read.set(self.read.get() + out.len());
        self.largest.set(self.largest.get().max(out.len()));
        self.bytes.read_exact_at(out, at)
    }
}

#[test]
fn auto_sql_reads_only_the_terminated_prefix_of_a_padded_section() {
    let bytes = fixture(b"table t\n\"test\"\n(\n)\0", 2 << 20);
    let source = CountReads {
        bytes: &bytes,
        read: Cell::new(0),
        largest: Cell::new(0),
    };
    let reader = BbiReader::open(&source).unwrap();
    source.read.set(0);
    source.largest.set(0);
    assert_eq!(
        reader.auto_sql().unwrap().as_deref(),
        Some("table t\n\"test\"\n(\n)")
    );
    assert!(source.read.get() <= 4096);
    assert!(source.largest.get() <= 4096);
}

#[test]
fn auto_sql_decodes_utf8_across_reads_and_terminators_at_read_edges() {
    for prefix in [0, 4095, 4096, 8191, 8192] {
        let mut schema = vec![b'x'; prefix];
        schema.extend_from_slice("染色体".as_bytes());
        let expected = String::from_utf8(schema.clone()).unwrap();
        schema.push(0);
        let bytes = fixture(&schema, 29);
        let reader = BbiReader::open(bytes).unwrap();
        assert_eq!(reader.auto_sql().unwrap(), Some(expected));
    }
    for len in [0, 4095, 4096, 8191] {
        let mut schema = vec![b'x'; len];
        schema.push(0);
        assert_eq!(
            BbiReader::open(fixture(&schema, 0))
                .unwrap()
                .auto_sql()
                .unwrap(),
            Some("x".repeat(len))
        );
    }
}

#[test]
fn auto_sql_rejects_unterminated_or_invalid_utf8_strings() {
    for schema in [vec![b'x'; 17], vec![b'x'; 8192], vec![0xff, 0]] {
        assert!(
            BbiReader::open(fixture(&schema, 0))
                .unwrap()
                .auto_sql()
                .is_err()
        );
    }
    assert_eq!(
        BbiReader::open(fixture(&[0, 0xff], 0))
            .unwrap()
            .auto_sql()
            .unwrap(),
        Some(String::new())
    );
}

#[test]
#[ignore = "root runs this identical source on baseline and candidate serially"]
fn paired_auto_sql_benchmark() {
    for padding in [0, 4096, 2 << 20] {
        let bytes = fixture(b"table t\n\"test\"\n(\n)\0", padding);
        let digest = checksum(&bytes);
        let source = CountReads {
            bytes: &bytes,
            read: Cell::new(0),
            largest: Cell::new(0),
        };
        let reader = BbiReader::open(&source).unwrap();
        source.read.set(0);
        assert_eq!(
            reader.auto_sql().unwrap().as_deref(),
            Some("table t\n\"test\"\n(\n)")
        );
        let read = source.read.get();
        for _ in 0..8 {
            black_box(reader.auto_sql().unwrap());
        }
        let mut samples = [0.0f64; 5];
        for sample in &mut samples {
            let start = Instant::now();
            let mut n = 0;
            while start.elapsed() < Duration::from_millis(100) {
                black_box(reader.auto_sql().unwrap());
                n += 1;
            }
            *sample = start.elapsed().as_nanos() as f64 / n as f64;
        }
        samples.sort_by(f64::total_cmp);
        println!(
            "reader.auto_sql.padding{padding}\tns={:.2}\tread={read}\tdigest={digest:016x}",
            samples[2]
        );
    }
}
