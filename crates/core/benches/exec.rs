//! Ordered queue bookkeeping with inline nested jobs and external workers.
use omicsweave_core::exec::{Ordered, Pool};
use std::hint::black_box;
use std::sync::Arc;
use std::time::{Duration, Instant};

fn batch(pool: &Pool, depth: usize) -> usize {
    let mut queue = pool.ordered::<usize>(depth);
    let mut sum = 0;
    for i in 0..depth {
        queue.submit(Box::new(move || Ok(black_box(i)))).unwrap();
    }
    while let Some(value) = queue.next(true).unwrap() {
        sum += value;
    }
    assert_eq!(sum, depth * (depth - 1) / 2);
    black_box(sum)
}

fn measure(pool: &Pool, depth: usize, name: &str) {
    batch(pool, depth);
    let mut samples = [0.0; 5];
    for sample in &mut samples {
        let start = Instant::now();
        let mut count = 0;
        while start.elapsed() < Duration::from_millis(100) {
            black_box(batch(pool, depth));
            count += depth;
        }
        *sample = start.elapsed().as_nanos() as f64 / count as f64;
    }
    samples.sort_by(f64::total_cmp);
    println!("{name}.{depth}.ns_per_job\t{:.2}", samples[2]);
}

fn main() {
    let pool = Arc::new(Pool::new(4));
    for depth in [1, 8, 128, 4096] {
        let nested = Arc::clone(&pool);
        pool.run(move || {
            // Same-pool submissions run inline, isolating queue bookkeeping
            // from worker dispatch and OS scheduling.
            measure(&nested, depth, "exec.nested");
            Ok(())
        })
        .unwrap();
        measure(&pool, depth, "exec.external");
    }
}
