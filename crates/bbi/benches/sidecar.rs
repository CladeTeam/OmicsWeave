//! Sidecar construction over warm memory; files and expected metadata are
//! prepared before timing. Declarations remain ceilings, not output sizes.
use std::hint::black_box;
use std::io::Cursor;
use std::time::{Duration, Instant};

use omicsweave_bbi::{BbiSidecar, BedSchema, BigBedWriter, BigWigWriter, Chrom, WriteOptions};
use omicsweave_core::sidecar::checksum;

fn measure(name: &str, bytes: &[u8], expected: &BbiSidecar) {
    if std::env::var("BBI_SIDECAR_BENCH_FILTER").is_ok_and(|f| !name.contains(&f)) {
        return;
    }
    assert_eq!(&BbiSidecar::build(bytes).unwrap(), expected);
    println!("{name}.digest\t{:016x}", checksum(&expected.body()));
    let mut samples = [0.0; 5];
    for sample in &mut samples {
        let start = Instant::now();
        let mut count = 0;
        while start.elapsed() < Duration::from_millis(100) {
            black_box(BbiSidecar::build(black_box(bytes)).unwrap());
            count += 1;
        }
        *sample = start.elapsed().as_nanos() as f64 / count as f64;
    }
    samples.sort_by(f64::total_cmp);
    println!("{name}\t{:.2}", samples[2]);
    assert_eq!(&BbiSidecar::build(bytes).unwrap(), expected);
}

fn wig(n: u32, items: u32, compress: bool) -> (Vec<u8>, BbiSidecar) {
    let mut writer = BigWigWriter::create(
        Cursor::new(Vec::new()),
        vec![Chrom {
            name: "c".into(),
            len: n * 8 + 8,
        }],
        WriteOptions {
            items_per_slot: Some(items),
            block_size: 8,
            compress,
            zoom_reductions: Some(vec![128, 1024]),
            ..WriteOptions::default()
        },
    )
    .unwrap();
    for i in 0..n {
        writer.add(0, i * 8, i * 8 + 4, (i % 97) as f32).unwrap();
    }
    let (out, sidecar) = writer.finish().unwrap();
    (out.into_inner(), sidecar)
}

fn measure_writer(n: u32, items: u32, compress: bool) {
    let name = format!("writer.wig.{n}.items{items}.compressed{compress}");
    if std::env::var("BBI_SIDECAR_BENCH_FILTER").is_ok_and(|f| !name.contains(&f)) {
        return;
    }
    let expected = wig(n, items, compress);
    println!("{name}.digest\t{:016x}", checksum(&expected.0));
    let mut samples = [0.0; 5];
    for sample in &mut samples {
        let start = Instant::now();
        let mut count = 0;
        while start.elapsed() < Duration::from_millis(100) {
            black_box(wig(black_box(n), items, compress));
            count += 1;
        }
        *sample = start.elapsed().as_nanos() as f64 / count as f64;
    }
    samples.sort_by(f64::total_cmp);
    println!("{name}\t{:.2}", samples[2]);
    assert_eq!(wig(n, items, compress), expected);
}

fn main() {
    for compress in [false, true] {
        for items in [16, 1024] {
            measure_writer(4096, items, compress);
        }
    }
    let (bytes, expected) = wig(4096, 64, true);
    measure("sidecar.wig.small", &bytes, &expected);
    for ceiling in [16u32 << 20, 64 << 20] {
        let mut bytes = bytes.clone();
        bytes[52..56].copy_from_slice(&ceiling.to_le_bytes());
        let mut expected = expected.clone();
        expected.uncompress_buf_size = ceiling;
        measure(&format!("sidecar.wig.ceiling{ceiling}"), &bytes, &expected);
    }
    let (bytes, expected) = wig(12_000, 12_000, true);
    measure("sidecar.wig.large", &bytes, &expected);

    let mut writer = BigBedWriter::create(
        Cursor::new(Vec::new()),
        vec![Chrom {
            name: "c".into(),
            len: 100_000,
        }],
        BedSchema::bed(4),
        WriteOptions {
            items_per_slot: Some(8),
            block_size: 8,
            zoom_reductions: Some(vec![128]),
            ..WriteOptions::default()
        },
    )
    .unwrap();
    let note = vec![b'a'; 10_000];
    for i in 0..64 {
        writer.add(0, i * 16, i * 16 + 8, &note).unwrap();
    }
    let (out, sidecar) = writer.finish().unwrap();
    measure("sidecar.bed.large", &out.into_inner(), &sidecar);
}
