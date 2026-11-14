//! Builder finish, excluding contact collection; output is a warm Vec.
use omicsweave_core::sidecar::checksum;
use omicsweave_hic::{BuildOptions, Chrom, HicBuilder, HicHeader};
use std::hint::black_box;
use std::io::Cursor;
use std::time::{Duration, Instant};

fn fixture(chroms: u32, resolutions: u32, value: f32) -> HicBuilder {
    let header = HicHeader::new(
        "bench",
        (0..chroms)
            .map(|i| Chrom {
                name: format!("c{i}"),
                len: 100_000,
            })
            .collect(),
        (1..=resolutions).rev().collect(),
    );
    let mut builder = HicBuilder::new(header).unwrap();
    for a in 0..chroms {
        for b in a..chroms {
            builder.add(a, 100, b, 200, value).unwrap();
        }
    }
    builder
}

fn measure(chroms: u32, resolutions: u32, value: f32) {
    let name = format!(
        "hic.finish.chr{chroms}.res{resolutions}.zero{}",
        value == 0.0
    );
    if std::env::var("HIC_BENCH_FILTER").is_ok_and(|s| !name.contains(&s)) {
        return;
    }
    let options = BuildOptions {
        level: 1,
        blocks: omicsweave_hic::BlockLayout::juicer(),
    };
    let expected = fixture(chroms, resolutions, value)
        .finish(Cursor::new(Vec::new()), &options)
        .unwrap()
        .into_inner();
    println!("{name}.digest\t{:016x}", checksum(&expected));
    let mut samples = [0.0; 5];
    for sample in &mut samples {
        let mut elapsed = Duration::ZERO;
        let mut count = 0;
        while elapsed < Duration::from_millis(100) {
            let builder = fixture(chroms, resolutions, value);
            let start = Instant::now();
            let (out, sidecar) = builder
                .finish_with_sidecar(Cursor::new(Vec::new()), black_box(&options))
                .unwrap();
            elapsed += start.elapsed();
            black_box((out, sidecar));
            count += 1;
        }
        *sample = elapsed.as_nanos() as f64 / count as f64;
    }
    samples.sort_by(f64::total_cmp);
    println!("{name}\t{:.2}", samples[2]);
}

fn collect(resolutions: u32) {
    let name = format!("hic.collect.chr4.res{resolutions}.contacts2048");
    if std::env::var("HIC_BENCH_FILTER").is_ok_and(|s| !name.contains(&s)) {
        return;
    }
    let header = HicHeader::new(
        "bench",
        (0..4)
            .map(|i| Chrom {
                name: format!("c{i}"),
                len: 100_000,
            })
            .collect(),
        (1..=resolutions).rev().collect(),
    );
    let contacts: Vec<_> = (0..2048u32)
        .map(|i| {
            (
                i % 4,
                i * 17 % 100_001,
                i / 4 % 4,
                i * 31 % 100_001,
                (i % 7) as f32 + 0.5,
            )
        })
        .collect();
    let build = || {
        let mut builder = HicBuilder::new(header.clone()).unwrap();
        for &(a, p, b, q, value) in &contacts {
            builder.add(a, p, b, q, value).unwrap();
        }
        builder
    };
    let bytes = build()
        .finish(
            Cursor::new(Vec::new()),
            &BuildOptions {
                level: 1,
                blocks: omicsweave_hic::BlockLayout::juicer(),
            },
        )
        .unwrap()
        .into_inner();
    println!("{name}.digest\t{:016x}", checksum(&bytes));
    let mut samples = [0.0; 5];
    for sample in &mut samples {
        let start = Instant::now();
        let mut count = 0;
        while start.elapsed() < Duration::from_millis(100) {
            black_box(build());
            count += 1;
        }
        *sample = start.elapsed().as_nanos() as f64 / count as f64;
    }
    samples.sort_by(f64::total_cmp);
    println!("{name}\t{:.2}", samples[2]);
}

fn main() {
    println!("case\tns/op (median of 5 x 100ms, contact collection excluded)");
    for resolutions in [1, 8, 64, 256] {
        collect(resolutions);
    }
    for resolutions in [1, 2, 8, 64] {
        measure(4, resolutions, 1.0);
    }
    for resolutions in [64, 256] {
        measure(12, resolutions, 0.0);
    }
}
