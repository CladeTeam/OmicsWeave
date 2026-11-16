#![allow(dead_code)]
/// splitmix64, for deterministic fixtures from a proptest-chosen seed.
pub struct Rng(u64);

impl Rng {
    pub fn new(seed: u64) -> Rng {
        Rng(seed)
    }

    pub fn next(&mut self) -> u64 {
        self.0 = self.0.wrapping_add(0x9E37_79B9_7F4A_7C15);
        let mut z = self.0;
        z = (z ^ (z >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
        z = (z ^ (z >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
        z ^ (z >> 31)
    }

    /// Uniform in `[lo, hi]`.
    pub fn range(&mut self, lo: u64, hi: u64) -> u64 {
        lo + self.next() % (hi - lo + 1)
    }

    pub fn chance(&mut self, percent: u64) -> bool {
        self.next() % 100 < percent
    }

    pub fn value(&mut self) -> f32 {
        if self.chance(10) {
            return 0.0;
        }
        (self.range(0, 2_000_000) as f32 - 1_000_000.0) / 1000.0
    }
}

/// Values equal, or both NaN.
pub fn same(a: f32, b: f32) -> bool {
    a == b || (a.is_nan() && b.is_nan())
}

pub(crate) fn regions(rng: &mut Rng, chroms: &[(String, u32)]) -> Vec<(String, u32, u32)> {
    let mut out: Vec<(String, u32, u32)> = chroms.iter().map(|(n, l)| (n.clone(), 0, *l)).collect();
    for _ in 0..30 {
        let (name, len) = &chroms[rng.range(0, chroms.len() as u64 - 1) as usize];
        let start = rng.range(0, *len as u64 + 100) as u32;
        out.push((name.clone(), start, start + rng.range(0, 5_000) as u32));
    }
    out
}
