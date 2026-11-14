//! Same-source baseline/candidate benchmark; fixtures and assertions are untimed.

#[path = "../tests/common/mod.rs"]
mod common;

use std::hint::black_box;
use std::time::{Duration, Instant};

use common::{MapSpec, Storage, write_cool};
use hdf5_metno_sys as sys;
use omicsweave_cool::{Cool, PixelBatch, PixelScratch};

fn measure(mut f: impl FnMut(), name: &str) {
    for trial in 0..5 {
        let start = Instant::now();
        let mut iterations = 0;
        while start.elapsed() < Duration::from_millis(100) {
            f();
            iterations += 1;
        }
        println!(
            "{name} trial={trial} ns/iter={:.1}",
            start.elapsed().as_nanos() as f64 / f64::from(iterations)
        );
    }
}

fn digest(batch: &PixelBatch) -> u64 {
    batch
        .bin1
        .iter()
        .zip(&batch.bin2)
        .zip(&batch.count)
        .fold(0u64, |sum, ((a, b), c)| {
            sum.wrapping_add(a.wrapping_mul(31) ^ b.wrapping_mul(17) ^ c.to_bits())
        })
}

fn main() {
    let dir = tempfile::tempdir().unwrap();
    let mut spec = MapSpec::new(&[("chr1", 51_200)], 100);
    spec.pixels = (0..32_768)
        .map(|i| (i / 256, i % 256 + 128, (i % 997 + 1) as f64))
        .collect();

    spec.storage = [Storage::plain(8); 3];
    let path = dir.path().join("index.cool");
    write_cool(&path, &spec);
    let verify = Cool::open(&path, "/").unwrap();
    let mut reference = PixelBatch::new();
    verify
        .read_pixels_library(0..verify.n_pixels(), &mut reference)
        .unwrap();
    let expected = digest(&reference);
    println!("cool_fixture digest={expected} pixels={}", reference.len());
    measure(
        || {
            let cool = Cool::open(black_box(&path), "/").unwrap();
            assert_eq!(cool.n_pixels(), 32_768);
            black_box(cool);
        },
        "cool_open_4096_chunks_per_column",
    );

    spec.storage = [Storage::gzip_shuffle(256); 3];
    let compressed_path = dir.path().join("compressed.cool");
    write_cool(&compressed_path, &spec);
    let compressed = Cool::open(&compressed_path, "/").unwrap();
    let mut scratch = PixelScratch::new();
    let mut batch = PixelBatch::new();
    compressed
        .read_pixels(0..32_768, &mut scratch, &mut batch)
        .unwrap();
    assert_eq!(batch, reference);
    assert_eq!(digest(&batch), expected);
    measure(
        || {
            compressed
                .read_pixels(black_box(0..32_768), &mut scratch, &mut batch)
                .unwrap();
            black_box(&batch);
        },
        "cool_gzip_shuffle_scan_32768",
    );

    // Each column has one full chunk whose two filters were skipped.
    spec.storage = [Storage::gzip_shuffle(32_768); 3];
    let skipped_path = dir.path().join("skipped.cool");
    write_cool(&skipped_path, &spec);
    {
        let file = hdf5_metno::File::open_rw(&skipped_path).unwrap();
        let columns = [
            (
                "pixels/bin1_id",
                spec.pixels
                    .iter()
                    .flat_map(|p| (p.0 as i64).to_ne_bytes())
                    .collect::<Vec<_>>(),
            ),
            (
                "pixels/bin2_id",
                spec.pixels
                    .iter()
                    .flat_map(|p| (p.1 as i64).to_ne_bytes())
                    .collect::<Vec<_>>(),
            ),
            (
                "pixels/count",
                spec.pixels
                    .iter()
                    .flat_map(|p| (p.2 as i32).to_ne_bytes())
                    .collect::<Vec<_>>(),
            ),
        ];
        for (name, raw) in columns {
            let ds = file.dataset(name).unwrap();
            let offset = [0u64];
            let _lock = sys::LOCK.lock();
            // SAFETY: each live 1-D dataset has one full 32768-row chunk.
            // raw contains exactly that many values in its native stored
            // numeric type. The single offset is aligned; mask 3 skips
            // the shuffle and deflate stages.
            let status = unsafe {
                sys::h5d::H5Dwrite_chunk(
                    ds.id(),
                    sys::h5p::H5P_DEFAULT,
                    3,
                    offset.as_ptr(),
                    raw.len(),
                    raw.as_ptr().cast(),
                )
            };
            assert!(status >= 0);
        }
    }
    let skipped = Cool::open(&skipped_path, "/").unwrap();
    skipped
        .read_pixels_library(137..138, &mut reference)
        .unwrap();
    skipped
        .read_pixels(137..138, &mut scratch, &mut batch)
        .unwrap();
    assert_eq!(batch, reference);
    println!("cool_skipped_fixture digest={}", digest(&batch));
    measure(
        || {
            skipped
                .read_pixels(black_box(137..138), &mut scratch, &mut batch)
                .unwrap();
            black_box(&batch);
        },
        "cool_skipped_filters_one_row_of_32768",
    );
}
