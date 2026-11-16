/// Sorted, unique pixels from `(bin1, bin2, count)` triples; with
/// `upper`, each pair is ordered so `bin1 <= bin2`.
pub fn normalize_pixels(raw: &[(u64, u64, f64)], upper: bool) -> Vec<(u64, u64, f64)> {
    let mut v: Vec<(u64, u64, f64)> = raw
        .iter()
        .map(|&(a, b, c)| if upper && a > b { (b, a, c) } else { (a, b, c) })
        .collect();
    v.sort_by_key(|x| (x.0, x.1));
    v.dedup_by(|x, y| x.0 == y.0 && x.1 == y.1);
    v
}

/// Dense-ish pseudo-random pixels over `n_bins` bins.
pub fn sample_pixels(n_bins: u64, n: usize, seed: u64, integer: bool) -> Vec<(u64, u64, f64)> {
    let mut state = seed
        .wrapping_mul(6364136223846793005)
        .wrapping_add(1442695040888963407);
    let mut next = || {
        state = state
            .wrapping_mul(6364136223846793005)
            .wrapping_add(1442695040888963407);
        state >> 33
    };
    let raw: Vec<(u64, u64, f64)> = (0..n)
        .map(|_| {
            let a = next() % n_bins;
            let b = next() % n_bins;
            let c = if integer {
                (next() % 1000 + 1) as f64
            } else {
                (next() % 100_000) as f64 / 7.0
            };
            (a, b, c)
        })
        .collect();
    normalize_pixels(&raw, true)
}
