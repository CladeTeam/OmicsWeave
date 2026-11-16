//! Deterministic coordinate-sorted SAM fixtures.

/// A small deterministic generator (splitmix64).
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

    /// Uniform in `lo..hi`.
    pub fn range(&mut self, lo: u64, hi: u64) -> u64 {
        lo + self.next() % (hi - lo).max(1)
    }

    pub fn chance(&mut self, percent: u64) -> bool {
        self.next() % 100 < percent
    }
}

pub struct SamSpec {
    pub refs: Vec<(&'static str, u64)>,
    pub per_ref: usize,
    /// Mean gap between read starts.
    pub gap: u64,
    pub long_cigars: bool,
}

fn seq(rng: &mut Rng, n: usize) -> String {
    let bases = b"ACGTN";
    (0..n)
        .map(|_| {
            let hi = if rng.chance(2) { 5 } else { 4 };
            bases[rng.range(0, hi) as usize] as char
        })
        .collect()
}

fn qual(rng: &mut Rng, n: usize) -> String {
    (0..n)
        .map(|_| (b'!' + rng.range(0, 41) as u8) as char)
        .collect()
}

/// A CIGAR of query length `l` mixing every operation.
fn cigar(rng: &mut Rng, l: usize) -> String {
    let mut out = String::new();
    let mut left = l;
    if rng.chance(20) && left > 10 {
        let s = rng.range(1, 6) as usize;
        out.push_str(&format!("{s}S"));
        left -= s;
    }
    if rng.chance(5) {
        out.insert_str(0, &format!("{}H", rng.range(1, 30)));
    }
    while left > 0 {
        let m = (rng.range(1, 60) as usize).min(left);
        let op = if rng.chance(30) {
            "="
        } else if rng.chance(10) {
            "X"
        } else {
            "M"
        };
        out.push_str(&format!("{m}{op}"));
        left -= m;
        if left == 0 {
            break;
        }
        match rng.range(0, 10) {
            0 => out.push_str(&format!("{}D", rng.range(1, 40))),
            1 => out.push_str(&format!("{}N", rng.range(50, 20_000))),
            2 if left > 2 => {
                let i = rng.range(1, 3) as usize;
                out.push_str(&format!("{i}I"));
                left -= i;
            }
            3 => out.push_str(&format!("{}P", rng.range(1, 4))),
            _ => {}
        }
    }
    out
}

fn rlen(cigar: &str) -> u64 {
    let mut n = 0u64;
    let mut num = 0u64;
    for c in cigar.chars() {
        if let Some(d) = c.to_digit(10) {
            num = num * 10 + d as u64;
        } else {
            if matches!(c, 'M' | 'D' | 'N' | '=' | 'X') {
                n += num;
            }
            num = 0;
        }
    }
    n
}

/// Random coordinate-sorted SAM text over `spec.refs`.
pub fn random_sam(rng: &mut Rng, spec: &SamSpec) -> Vec<u8> {
    let mut s = String::from("@HD\tVN:1.6\tSO:coordinate\n");
    for (name, len) in &spec.refs {
        s.push_str(&format!("@SQ\tSN:{name}\tLN:{len}\n"));
    }
    s.push_str("@RG\tID:rg1\tSM:x\n");
    let long_cigar = spec.long_cigars.then(|| "1M1I".repeat(33_000));
    let mut serial = 0;
    for (name, len) in &spec.refs {
        // Starts spread over the whole reference, some in dense piles.
        let mut starts = Vec::new();
        while starts.len() < spec.per_ref {
            let room = (*len / 4).min(40_000);
            let p = rng.range(1, len - room);
            starts.push(p);
            if rng.chance(10) {
                for _ in 0..rng.range(5, 40) {
                    starts.push(p + rng.range(0, spec.gap));
                }
            }
        }
        starts.sort_unstable();
        for pos in starts {
            serial += 1;
            let l = rng.range(20, 151) as usize;
            let tags = format!("RG:Z:rg1\tNM:i:{}\tXS:A:+", rng.range(0, 30));
            let (q, s_) = (
                if rng.chance(5) {
                    "*".to_string()
                } else {
                    qual(rng, l)
                },
                seq(rng, l),
            );
            if rng.chance(3) {
                // Unmapped, placed beside its mate.
                s.push_str(&format!(
                    "u{serial}\t{}\t{name}\t{pos}\t0\t*\t=\t{pos}\t0\t{s_}\t{q}\t{tags}\n",
                    4 | 1 | 0x80
                ));
                continue;
            }
            let c = if spec.long_cigars && rng.chance(1) {
                // More operations than a BAM CIGAR holds: htslib moves them
                // into a CG tag.
                let n = 33_000;
                let c = long_cigar.as_deref().unwrap();
                let seq_l = 2 * n;
                s.push_str(&format!(
                    "long{serial}\t0\t{name}\t{pos}\t60\t{c}\t*\t0\t0\t{}\t*\n",
                    seq(rng, seq_l)
                ));
                continue;
            } else {
                cigar(rng, l)
            };
            let span = rlen(&c).max(1);
            if pos + span >= *len {
                continue;
            }
            let paired = rng.chance(60);
            let reverse = rng.chance(50);
            let mut flag = if reverse { 0x10 } else { 0 };
            for (bit, p) in [(0x100, 3), (0x200, 2), (0x400, 8), (0x800, 2)] {
                if rng.chance(p) {
                    flag |= bit;
                }
            }
            let (rnext, pnext, tlen) = if paired {
                flag |= 1 | if rng.chance(85) { 2 } else { 0 };
                flag |= if rng.chance(50) { 0x40 } else { 0x80 };
                if rng.chance(5) {
                    flag |= 8;
                    ("=".to_string(), pos, 0i64)
                } else {
                    let d = rng.range(0, 600);
                    if rng.chance(50) {
                        let mate = pos + d;
                        ("=".to_string(), mate, (mate + 100 - pos) as i64)
                    } else {
                        let mate = pos.saturating_sub(d).max(1);
                        ("=".to_string(), mate, -((pos + span - mate) as i64))
                    }
                }
            } else {
                ("*".to_string(), 0, 0)
            };
            s.push_str(&format!(
                "r{serial}\t{flag}\t{name}\t{pos}\t{}\t{c}\t{rnext}\t{pnext}\t{tlen}\t{s_}\t{q}\t{tags}\n",
                rng.range(0, 61)
            ));
        }
    }
    for i in 0..rng.range(0, 30) {
        s.push_str(&format!("x{i}\t4\t*\t0\t0\t*\t*\t0\t0\tACGT\t*\n"));
    }
    s.into_bytes()
}
