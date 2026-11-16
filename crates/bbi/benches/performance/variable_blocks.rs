use omicsweave_bbi::{BbiSidecar, BedSchema, BigBedWriter, Chrom, WriteOptions};
use omicsweave_core::exec::{Exec, Pool};
use omicsweave_core::sidecar::checksum;
use std::hint::black_box;
use std::io::Cursor;
use std::time::{Duration, Instant};
fn chroms() -> Vec<Chrom> {
    vec![Chrom {
        name: "c".into(),
        len: 100_000,
    }]
}

fn options(compress: bool, zooms: bool, slots: u32) -> WriteOptions {
    WriteOptions {
        compress,
        items_per_slot: Some(slots),
        zoom_reductions: Some(if zooms { vec![8, 128] } else { vec![] }),
        ..WriteOptions::default()
    }
}

fn bed(exec: Exec<'_>, compress: bool, zooms: bool, slots: u32) -> (Vec<u8>, BbiSidecar) {
    let mut writer = BigBedWriter::create_with(
        Cursor::new(Vec::new()),
        chroms(),
        BedSchema::bed(4),
        options(compress, zooms, slots),
        exec,
    )
    .unwrap();
    let notes: Vec<_> = [1, 8192, 2, 257]
        .into_iter()
        .map(|len| vec![b'x'; len])
        .collect();
    for i in 0..257 {
        writer
            .add(0, i * 4, i * 4 + 17, &notes[i as usize % notes.len()])
            .unwrap();
    }
    let (out, sidecar) = writer.finish().unwrap();
    (out.into_inner(), sidecar)
}

pub(crate) fn variable_block_benchmark() {
    let pool = Pool::new(2);
    for (mode, exec) in [
        ("inline", Exec::Inline),
        (
            "pool",
            Exec::Pool {
                pool: &pool,
                limit: 4,
            },
        ),
    ] {
        for compress in [false, true] {
            for zooms in [false, true] {
                for slots in [1, 64] {
                    let name =
                        format!("writer.bed.{mode}.compressed{compress}.zooms{zooms}.slots{slots}");
                    let expected = bed(exec, compress, zooms, slots);
                    let digest = checksum(&expected.0);
                    for _ in 0..4 {
                        black_box(bed(exec, compress, zooms, slots));
                    }
                    let mut samples = [0.0f64; 5];
                    for sample in &mut samples {
                        let start = Instant::now();
                        let mut n = 0;
                        while start.elapsed() < Duration::from_millis(100) {
                            black_box(bed(exec, compress, zooms, slots));
                            n += 1;
                        }
                        *sample = start.elapsed().as_nanos() as f64 / n as f64;
                    }
                    samples.sort_by(f64::total_cmp);
                    println!("{name}\tns={:.2}\tdigest={digest:016x}", samples[2]);
                }
            }
        }
    }
}
