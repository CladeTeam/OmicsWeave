//! Deterministic input models; no production writer is used here.
use omicsweave_core::ContactUnit;
use std::collections::BTreeMap;

/// splitmix64, for deterministic fixtures from a proptest-chosen seed.
pub struct Rng(pub u64);

impl Rng {
    pub fn next(&mut self) -> u64 {
        self.0 = self.0.wrapping_add(0x9E37_79B9_7F4A_7C15);
        let mut z = self.0;
        z = (z ^ (z >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
        z = (z ^ (z >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
        z ^ (z >> 31)
    }

    pub fn range(&mut self, lo: u64, hi: u64) -> u64 {
        lo + self.next() % (hi - lo + 1)
    }

    pub fn chance(&mut self, percent: u64) -> bool {
        self.next() % 100 < percent
    }
}

/// One chromosome pair's map at one resolution. `chr1 <= chr2`; within a
/// chromosome, records keep `x <= y` (the stored half).
#[derive(Debug, Clone)]
pub struct Map {
    pub chr1: u32,
    pub chr2: u32,
    pub unit: ContactUnit,
    pub res: u32,
    pub records: Vec<(u32, u32, f32)>,
}

#[derive(Debug, Clone)]
pub struct Norm {
    pub kind: String,
    pub chrom: u32,
    pub unit: ContactUnit,
    pub res: u32,
    pub values: Vec<f64>,
}

/// An expected vector; `kind` "NONE" for raw counts.
#[derive(Debug, Clone)]
pub struct Expected {
    pub kind: String,
    pub unit: ContactUnit,
    pub res: u32,
    pub values: Vec<f64>,
    pub factors: Vec<(u32, f64)>,
}

#[derive(Debug, Clone)]
pub struct Spec {
    pub version: i32,
    pub chroms: Vec<(String, u64)>,
    pub bp_resolutions: Vec<u32>,
    pub frag_resolutions: Vec<u32>,
    pub sites: Vec<Vec<u32>>,
    pub maps: Vec<Map>,
    pub norms: Vec<Norm>,
    pub expected: Vec<Expected>,
    pub block_bins: u32,
    /// Integral values stored as short; otherwise float.
    pub short_values: bool,
    pub level: u32,
}

pub const NORMS: [&str; 2] = ["KR", "VC"];

impl Spec {
    pub fn resolutions(&self, unit: ContactUnit) -> &[u32] {
        match unit {
            ContactUnit::Bp => &self.bp_resolutions,
            ContactUnit::Frag => &self.frag_resolutions,
        }
    }

    /// A chromosome's span in `unit`.
    pub fn extent(&self, chrom: u32, unit: ContactUnit) -> u64 {
        match unit {
            ContactUnit::Bp => self.chroms[chrom as usize].1,
            ContactUnit::Frag => self.sites[chrom as usize].len() as u64,
        }
    }

    /// Bins of a chromosome (the last one partial).
    pub fn bins(&self, chrom: u32, unit: ContactUnit, res: u32) -> u32 {
        (self.extent(chrom, unit) / res as u64) as u32 + 1
    }

    pub fn units(&self) -> Vec<ContactUnit> {
        if self.frag_resolutions.is_empty() {
            vec![ContactUnit::Bp]
        } else {
            vec![ContactUnit::Bp, ContactUnit::Frag]
        }
    }

    pub fn norm(&self, kind: &str, chrom: u32, unit: ContactUnit, res: u32) -> Option<&[f64]> {
        self.norms
            .iter()
            .find(|n| n.kind == kind && n.chrom == chrom && n.unit == unit && n.res == res)
            .map(|n| n.values.as_slice())
    }

    pub fn expected(&self, kind: &str, unit: ContactUnit, res: u32) -> Option<&Expected> {
        self.expected
            .iter()
            .find(|e| e.kind == kind && e.unit == unit && e.res == res)
    }

    /// A vector value as stored (float from v9).
    pub fn stored(&self, v: f64) -> f64 {
        if self.version > 8 { v as f32 as f64 } else { v }
    }
}

pub fn random_spec(rng: &mut Rng, version: i32) -> Spec {
    let chroms: Vec<(String, u64)> = (0..rng.range(1, 3))
        .map(|i| (format!("chr{}", i + 1), rng.range(20_000, 200_000)))
        .collect();
    let (frag_resolutions, sites) = if rng.chance(50) {
        let mut sites = Vec::new();
        for (_, len) in &chroms {
            let mut s = Vec::new();
            for _ in 0..rng.range(30, 600) {
                s.push(rng.range(1, *len - 1) as u32);
            }
            s.sort_unstable();
            s.dedup();
            sites.push(s);
        }
        (vec![5, 1], sites)
    } else {
        (Vec::new(), Vec::new())
    };
    let mut spec = Spec {
        version,
        chroms,
        bp_resolutions: vec![5_000, 1_000],
        frag_resolutions,
        sites,
        maps: Vec::new(),
        norms: Vec::new(),
        expected: Vec::new(),
        block_bins: rng.range(4, 60) as u32,
        short_values: rng.chance(50),
        level: 0,
    };
    // As the Juicer tools write them: a pair present in the file has a
    // record at every resolution (possibly without blocks).
    let n = spec.chroms.len() as u32;
    let mut pairs = Vec::new();
    for c1 in 0..n {
        for c2 in c1..n {
            if c1 == c2 || rng.chance(70) {
                pairs.push((c1, c2));
            }
        }
    }
    for unit in spec.units() {
        for res in spec.resolutions(unit).to_vec() {
            for &(c1, c2) in &pairs {
                let (n1, n2) = (spec.bins(c1, unit, res) - 1, spec.bins(c2, unit, res) - 1);
                let mut cells = BTreeMap::new();
                for _ in 0..rng.range(0, 3_000) {
                    let (mut x, mut y) = (
                        rng.range(0, n1 as u64) as u32,
                        rng.range(0, n2 as u64) as u32,
                    );
                    if c1 == c2 {
                        // Contacts concentrate near the diagonal.
                        let reach = 1 << rng.range(0, 7);
                        y = (x + rng.range(0, reach) as u32).min(n2);
                        if y < x {
                            std::mem::swap(&mut x, &mut y);
                        }
                    }
                    let v = if spec.short_values {
                        rng.range(1, 500) as f32
                    } else {
                        rng.range(1, 100_000) as f32 / 37.0
                    };
                    cells.insert((x, y), v);
                }
                spec.maps.push(Map {
                    chr1: c1,
                    chr2: c2,
                    unit,
                    res,
                    records: cells.into_iter().map(|((x, y), v)| (x, y, v)).collect(),
                });
            }
            for kind in NORMS {
                for chrom in 0..n {
                    let values = (0..spec.bins(chrom, unit, res))
                        .map(|_| match rng.range(0, 30) {
                            0 => f64::NAN,
                            1 => 0.0,
                            _ => rng.range(1, 4_000) as f64 / 1000.0,
                        })
                        .collect();
                    spec.norms.push(Norm {
                        kind: kind.into(),
                        chrom,
                        unit,
                        res,
                        values,
                    });
                }
            }
            let longest = (0..n).map(|c| spec.bins(c, unit, res)).max().unwrap();
            for kind in ["NONE", NORMS[0], NORMS[1]] {
                // Shorter vectors than a chromosome too: distances past the
                // end read the last value.
                let len = rng.range(1, longest as u64 + 3);
                let values = (0..len)
                    .map(|_| rng.range(1, 90_000) as f64 / 7.0)
                    .collect();
                let mut factors = Vec::new();
                for c in 0..n {
                    if rng.chance(70) {
                        factors.push((c, rng.range(1, 3_000) as f64 / 1000.0));
                    }
                }
                spec.expected.push(Expected {
                    kind: kind.into(),
                    unit,
                    res,
                    values,
                    factors,
                });
            }
        }
    }
    spec
}
