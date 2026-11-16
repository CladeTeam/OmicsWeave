//! Smoke check on a real `.mcool`, named by `OW_COOL_SMOKE` (read only):
//!
//! ```sh
//! OW_COOL_SMOKE=/path/to/file.mcool cargo test --release -p omicsweave-cool \
//!     --test real -- --ignored --nocapture
//! ```
//!
//! Every resolution's spans are spot-checked against libhdf5's reader; the
//! finest resolution (or `OW_COOL_SMOKE_RES`) is then streamed whole on
//! `OW_COOL_SMOKE_THREADS` threads (default 16) and its totals checked
//! against the file's `nnz` / `sum` attributes.

use std::path::PathBuf;
use std::sync::Mutex;
use std::sync::atomic::{AtomicUsize, Ordering};

use omicsweave_cool::{Cool, PixelBatch, PixelScratch, mcool_resolutions};

fn env_u64(name: &str) -> Option<u64> {
    std::env::var(name).ok().map(|v| v.parse().unwrap())
}

#[test]
#[ignore = "needs OW_COOL_SMOKE=<path to a real .mcool>"]
fn real_mcool() {
    let path = PathBuf::from(std::env::var("OW_COOL_SMOKE").expect("set OW_COOL_SMOKE"));
    let resolutions = mcool_resolutions(&path).unwrap();
    assert!(
        !resolutions.is_empty(),
        "{} has no resolutions",
        path.display()
    );
    println!("{}: resolutions {resolutions:?}", path.display());

    let mut scratch = PixelScratch::new();
    let (mut fast, mut library) = (PixelBatch::new(), PixelBatch::new());
    for &res in &resolutions {
        let cool = Cool::open(&path, &format!("/resolutions/{res}")).unwrap();
        let spans = cool.pixel_spans(1 << 20);
        println!(
            "res {res:>8}: {} chroms, {} bins, {} pixels, paths {:?}, {} spans",
            cool.chroms().len(),
            cool.n_bins(),
            cool.n_pixels(),
            cool.read_paths(),
            spans.len(),
        );
        let picks = [
            0,
            spans.len() / 3,
            spans.len() / 2,
            spans.len().saturating_sub(1),
        ];
        for &i in &picks {
            let Some(span) = spans.get(i) else { continue };
            cool.read_pixels(span.clone(), &mut scratch, &mut fast)
                .unwrap();
            cool.read_pixels_library(span.clone(), &mut library)
                .unwrap();
            assert_eq!(fast, library, "res {res} span {span:?}");
        }
    }

    let res = env_u64("OW_COOL_SMOKE_RES")
        .map(|r| u32::try_from(r).expect("OW_COOL_SMOKE_RES must fit u32"))
        .unwrap_or(resolutions[0]);
    let threads = usize::try_from(env_u64("OW_COOL_SMOKE_THREADS").unwrap_or(16))
        .expect("OW_COOL_SMOKE_THREADS must fit usize");
    assert!(threads > 0, "OW_COOL_SMOKE_THREADS must be positive");
    let group = format!("/resolutions/{res}");
    let cool = Cool::open(&path, &group).unwrap();
    let spans = cool.pixel_spans(1 << 20);

    let next = AtomicUsize::new(0);
    let totals = Mutex::new((0u64, 0f64));
    std::thread::scope(|scope| {
        for _ in 0..threads {
            scope.spawn(|| {
                let mut scratch = PixelScratch::new();
                let mut batch = PixelBatch::new();
                let (mut rows, mut sum) = (0u64, 0f64);
                loop {
                    let i = next.fetch_add(1, Ordering::Relaxed);
                    let Some(span) = spans.get(i) else { break };
                    cool.read_pixels(span.clone(), &mut scratch, &mut batch)
                        .unwrap();
                    assert!(
                        batch.bin1.windows(2).all(|w| w[0] <= w[1]),
                        "bin1 not sorted in {span:?}"
                    );
                    rows += batch.len() as u64;
                    sum += batch.count.iter().sum::<f64>();
                }
                let mut t = totals.lock().unwrap();
                t.0 += rows;
                t.1 += sum;
            });
        }
    });
    let (rows, sum) = totals.into_inner().unwrap();
    println!("res {res}: streamed {rows} pixels on {threads} threads, count sum {sum}");
    assert_eq!(rows, cool.n_pixels());

    let file = hdf5_metno::File::open(&path).unwrap();
    let grp = file.group(&group).unwrap();
    if let Ok(nnz) = grp.attr("nnz").and_then(|a| a.read_scalar::<i64>()) {
        assert_eq!(nnz as u64, rows);
    }
    if let Ok(total) = grp.attr("sum").and_then(|a| a.read_scalar::<f64>()) {
        assert!(
            (total - sum).abs() <= 1e-9 * total.abs().max(1.0),
            "{total} vs {sum}"
        );
    }
}
