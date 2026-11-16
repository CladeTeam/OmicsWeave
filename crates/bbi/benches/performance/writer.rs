use crate::allocator::counted;
use crate::writer_fixture::{Discard, bed, chroms, options, wig};
use omicsweave_bbi::BigWigWriter;
use omicsweave_core::sidecar::checksum;
use std::hint::black_box;
use std::time::{Duration, Instant};
fn measure<T>(name: &str, mut f: impl FnMut() -> T) {
    for _ in 0..8 {
        black_box(f());
    }
    let mut samples = [0.0; 5];
    for sample in &mut samples {
        let start = Instant::now();
        let mut n = 0;
        while start.elapsed() < Duration::from_millis(100) {
            black_box(f());
            n += 1;
        }
        *sample = start.elapsed().as_nanos() as f64 / n as f64;
    }
    samples.sort_by(f64::total_cmp);
    let (_, largest, allocated) = counted(|| black_box(f()));
    println!(
        "{name}\tns={:.2}\tlargest={largest}\tallocated={allocated}",
        samples[2]
    );
}

pub(crate) fn writer_benchmark() {
    for zooms in [0, 2] {
        measure(&format!("writer.create.zooms{zooms}"), || {
            BigWigWriter::create(Discard::default(), chroms(), options(false, zooms)).unwrap()
        });
        for n in [1, 9_999, 10_000, 20_000] {
            let fixture = wig(n, false, zooms);
            println!(
                "writer.wig.n{n}.zooms{zooms}.digest\t{:016x}",
                checksum(&fixture.0)
            );
            measure(&format!("writer.wig.n{n}.zooms{zooms}"), || {
                wig(black_box(n), false, zooms)
            });
        }
    }
    for n in [1, 9_999, 10_000, 20_000] {
        let fixture = bed(n, false);
        println!(
            "writer.bed.n{n}.zooms0.digest\t{:016x}",
            checksum(&fixture.0)
        );
        measure(&format!("writer.bed.n{n}.zooms0"), || {
            bed(black_box(n), false)
        });
    }
}

pub(crate) fn warm_flush_allocations() {
    for compress in [false, true] {
        let mut write_options = options(false, 0);
        write_options.compress = compress;
        write_options.items_per_slot = Some(16);
        let mut writer = BigWigWriter::create(Discard::default(), chroms(), write_options).unwrap();
        let warm = 2_050 * 16;
        for i in 0..warm {
            writer.add(0, i * 4, i * 4 + 2, (i % 31) as f32).unwrap();
        }
        // Retained block metadata has ample spare capacity after warm-up;
        // the counted phase isolates repeated block cutting and compression.
        let (_, largest, allocated) = counted(|| {
            for i in warm..warm + 512 * 16 {
                writer.add(0, i * 4, i * 4 + 2, (i % 31) as f32).unwrap();
            }
        });
        println!(
            "writer.warm512flushes.compressed{compress}\tlargest={largest}\tallocated={allocated}"
        );
        black_box(writer.finish().unwrap());
    }
}
