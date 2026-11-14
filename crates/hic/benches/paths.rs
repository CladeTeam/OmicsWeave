//! Warm caller buffers; synthetic in-memory fixtures, excluding setup.
use omicsweave_core::{ContactQuery, ContactRead, ContactUnit, Region, sidecar::checksum};
use omicsweave_hic::{
    BlockBins, BlockEntry, BlockLayout, Chrom, Depth, HicHeader, HicReader, HicSidecar, HicWriter,
    Record, ZoomData, decode_block, encode_block,
};
use std::{
    hint::black_box,
    io::Cursor,
    time::{Duration, Instant},
};

fn measure(name: &str, mut f: impl FnMut()) {
    if std::env::var("HIC_BENCH_FILTER").is_ok_and(|s| !name.contains(&s)) {
        return;
    }
    for _ in 0..5 {
        f();
    }
    let mut samples = [0.0f64; 5];
    for sample in &mut samples {
        let start = Instant::now();
        let mut n = 0;
        while start.elapsed() < Duration::from_millis(100) {
            f();
            n += 1;
        }
        *sample = start.elapsed().as_nanos() as f64 / n as f64;
    }
    samples.sort_by(f64::total_cmp);
    println!("{name}\t{:.2}", samples[2]);
}

fn whole_writer() {
    let header = HicHeader::new(
        "whole",
        vec![Chrom {
            name: "a".into(),
            len: 512,
        }],
        vec![1],
    );
    for side in [32, 256] {
        let records: Vec<_> = (0..side)
            .flat_map(|y| {
                (0..=y).map(move |x| Record {
                    x,
                    y,
                    value: (x % 7 + y % 11 + 1) as f32,
                })
            })
            .collect();
        let name = format!("hic.writer.whole.side{side}");
        let write = || {
            let mut writer = HicWriter::create_with(
                Cursor::new(Vec::new()),
                header.clone(),
                0,
                BlockLayout::juicer(),
            )
            .unwrap();
            writer
                .write_matrix(0, 0, ContactUnit::Bp, 1, black_box(&records))
                .unwrap();
            writer.finish().unwrap().into_inner()
        };
        println!("{name}.digest\t{:016x}", checksum(&write()));
        measure(&name, || {
            black_box(write());
        });
    }
}

fn main() {
    println!("case\tns/op (median of 5 x 100ms)");
    whole_writer();
    let zoom = ZoomData {
        chr1: 0,
        chr2: 1,
        unit: ContactUnit::Bp,
        bin_size: 1000,
        sum_counts: 1.0,
        block_bin_count: 1,
        block_column_count: 1024,
        blocks: (0..4096)
            .map(|i| BlockEntry {
                number: i,
                position: 0,
                size: 1,
            })
            .collect(),
    };
    let mut selected = Vec::new();
    for width in [8, 32, 256] {
        measure(&format!("hic.select.width{width}"), || {
            zoom.blocks_for_into(
                8,
                Depth::Log(2),
                black_box(0..width),
                black_box(0..width),
                black_box(&mut selected),
            );
        });
    }
    let wide_zoom = ZoomData {
        block_column_count: 16,
        blocks: (0..200_000)
            .filter(|number| number % 3 != 0)
            .map(|number| BlockEntry {
                number,
                position: 0,
                size: 1,
            })
            .collect(),
        ..zoom
    };
    for end in [200, u32::MAX] {
        measure(&format!("hic.select.widecols.end{end}"), || {
            wide_zoom.blocks_for_into(
                8,
                Depth::Log(2),
                black_box(3..40),
                black_box(100..end),
                black_box(&mut selected),
            );
            black_box(&selected);
        });
    }
    for dense in [false, true] {
        let records: Vec<_> = (0..64u32)
            .flat_map(|y| {
                (0..64u32)
                    .filter(move |x| dense || (x + y) % 7 == 0)
                    .map(move |x| Record {
                        x,
                        y,
                        value: (1 + x + y) as f32 / 3.0,
                    })
            })
            .collect();
        let mut encoded = Vec::new();
        encode_block(9, &records, &mut encoded).unwrap();
        measure(&format!("hic.decode.dense{dense}"), || {
            let mut sum = 0.0;
            decode_block(9, black_box(&encoded), |x, y, v| {
                sum += x as f32 + y as f32 + v
            })
            .unwrap();
            black_box(sum);
        });
        println!(
            "hic.encode.dense{dense}.digest\t{:016x}",
            checksum(&encoded)
        );
        measure(&format!("hic.encode.dense{dense}"), || {
            encode_block(9, black_box(&records), black_box(&mut encoded)).unwrap();
            black_box(&encoded);
        });
    }
    let header = HicHeader::new(
        "test",
        vec![Chrom {
            name: "s".into(),
            len: 512_000,
        }],
        vec![1000],
    );
    let mut writer =
        HicWriter::create_with(Cursor::new(Vec::new()), header, 1, BlockLayout::juicer()).unwrap();
    let records: Vec<_> = (0..512u32)
        .flat_map(|x| {
            (x..(x + 8).min(512)).map(move |y| Record {
                x,
                y,
                value: (x + y + 1) as f32,
            })
        })
        .collect();
    writer
        .write_matrix(0, 0, ContactUnit::Bp, 1000, &records)
        .unwrap();
    let bytes = writer.finish().unwrap().into_inner();
    let mut reader = HicReader::open(&bytes[..]).unwrap();
    let query = ContactQuery::new(
        Region::new("s", 0, 512_000),
        Region::new("s", 0, 512_000),
        1000,
    );
    let mut out = vec![0.0; 512 * 512];
    measure("hic.reader.matrix512", || {
        reader
            .matrix(black_box(query), black_box(&mut out))
            .unwrap();
        black_box(&out);
    });
    let header = HicHeader::new(
        "sidecar",
        vec![Chrom {
            name: "s".into(),
            len: 65_536,
        }],
        vec![1, 2, 4],
    );
    let mut writer = HicWriter::create_with(
        Cursor::new(Vec::new()),
        header,
        1,
        BlockLayout {
            intra: BlockBins::Fixed(16),
            inter: BlockBins::Fixed(16),
        },
    )
    .unwrap();
    let records: Vec<_> = (0..16_384)
        .map(|x| Record {
            x,
            y: x,
            value: (x % 31) as f32 + 0.5,
        })
        .collect();
    for res in [1, 2, 4] {
        writer
            .write_matrix(0, 0, ContactUnit::Bp, res, &records)
            .unwrap();
    }
    let (out, expected) = writer.finish_with_sidecar().unwrap();
    let bytes = out.into_inner();
    assert_eq!(HicSidecar::build(bytes.as_slice()).unwrap(), expected);
    println!(
        "hic.sidecar.adjacent.digest\t{:016x}",
        checksum(&expected.body())
    );
    measure("hic.sidecar.adjacent", || {
        black_box(HicSidecar::build(black_box(bytes.as_slice())).unwrap());
    });
}
