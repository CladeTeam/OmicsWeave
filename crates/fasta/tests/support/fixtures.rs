//! Independent deterministic fixtures; byte writers do not call the implementation.

use omicsweave_fasta::{Fai, FaiBuilder, Result, TWOBIT_HEADER_LEN};

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

    /// Uniform coordinate in `[lo, hi]`.
    pub fn coord(&mut self, lo: u32, hi: u32) -> u32 {
        self.range(lo as u64, hi as u64) as u32
    }

    pub fn chance(&mut self, percent: u64) -> bool {
        self.next() % 100 < percent
    }

    pub fn pick<'a, T>(&mut self, items: &'a [T]) -> &'a T {
        &items[self.range(0, items.len() as u64 - 1) as usize]
    }
}

/// Random FASTA text: several records with varying line widths, `\n` or
/// `\r\n` endings, descriptions, blank lines, empty and duplicate sequences.
/// With `damage`, one edit likely makes it invalid.
pub fn random_fasta(rng: &mut Rng, damage: bool) -> Vec<u8> {
    const BASES: &[u8] = b"ACGTACGTACGTNacgtnRYKM";
    let eol: &[u8] = if rng.chance(20) { b"\r\n" } else { b"\n" };
    let mut out = Vec::new();
    let records = rng.range(1, 6);
    for i in 0..records {
        if i > 0 && rng.chance(15) {
            out.extend_from_slice(eol);
        }
        out.push(b'>');
        if rng.chance(10) {
            out.push(b' ');
        }
        let name = if i > 0 && rng.chance(10) {
            "dup".to_owned()
        } else {
            format!("seq{i}_{}", rng.next() % 1000)
        };
        out.extend_from_slice(name.as_bytes());
        if rng.chance(30) {
            out.extend_from_slice(if rng.chance(50) {
                b" description here"
            } else {
                b"\tdesc"
            });
        }
        out.extend_from_slice(eol);
        let len = if rng.chance(5) {
            0
        } else {
            rng.range(1, 3_000) as usize
        };
        let width = rng.range(1, 100) as usize;
        let bases: Vec<u8> = (0..len).map(|_| *rng.pick(BASES)).collect();
        let lines: Vec<&[u8]> = bases.chunks(width).collect();
        for (k, line) in lines.iter().enumerate() {
            out.extend_from_slice(line);
            let last = i + 1 == records && k + 1 == lines.len();
            if !last || rng.chance(90) {
                out.extend_from_slice(eol);
            }
        }
    }
    if damage && !out.is_empty() {
        let at = rng.range(0, out.len() as u64 - 1) as usize;
        match rng.range(0, 3) {
            0 => out.insert(at, b'A'),
            1 => out.insert(at, b'\n'),
            2 => out.truncate(at),
            _ => out
                .splice(at..at, b"\nxyz\n".iter().copied())
                .for_each(drop),
        }
    }
    out
}

/// Builds the `.fai` of `text`, fed in random chunks.
pub fn build_fai(text: &[u8], rng: &mut Rng) -> Result<Fai> {
    let mut builder = FaiBuilder::new();
    let mut at = 0;
    while at < text.len() {
        let n = (rng.range(1, 4_096) as usize).min(text.len() - at);
        builder.push(&text[at..at + n])?;
        at += n;
    }
    builder.finish()
}

/// One sequence for the 2bit test writer.
#[derive(Debug, Clone)]
pub struct TwoBitSeq {
    pub name: String,
    /// Uppercase ACGT.
    pub bases: Vec<u8>,
    pub n_blocks: Vec<(u32, u32)>,
    pub mask_blocks: Vec<(u32, u32)>,
}

fn random_blocks(rng: &mut Rng, len: u32, mean_gap: u64) -> Vec<(u32, u32)> {
    let mut blocks = Vec::new();
    let mut pos = rng.range(0, mean_gap) as u32;
    while pos < len {
        let size = (rng.range(1, 200) as u32).min(len - pos);
        blocks.push((pos, size));
        pos += size + rng.range(1, mean_gap) as u32;
    }
    blocks
}

pub fn random_twobit(rng: &mut Rng) -> Vec<TwoBitSeq> {
    (0..rng.range(1, 6))
        .map(|i| {
            let len = if rng.chance(5) {
                0
            } else {
                rng.range(1, 5_000) as u32
            };
            TwoBitSeq {
                name: format!("chr{i}_{}", rng.next() % 1000),
                bases: (0..len).map(|_| *rng.pick(b"ACGT")).collect(),
                n_blocks: if rng.chance(70) {
                    random_blocks(rng, len, 800)
                } else {
                    Vec::new()
                },
                mask_blocks: if rng.chance(70) {
                    random_blocks(rng, len, 300)
                } else {
                    Vec::new()
                },
            }
        })
        .collect()
}

/// Writes a 2bit file in the layout of `faToTwoBit`: header, directory,
/// then records in directory order.
pub fn write_twobit(seqs: &[TwoBitSeq], version: u32) -> Vec<u8> {
    let width = if version == 1 { 8 } else { 4 };
    let dir_len: usize = seqs.iter().map(|s| 1 + s.name.len() + width).sum();
    let records: Vec<Vec<u8>> = seqs
        .iter()
        .map(|s| {
            let mut r = Vec::new();
            r.extend_from_slice(&(s.bases.len() as u32).to_le_bytes());
            for blocks in [&s.n_blocks, &s.mask_blocks] {
                r.extend_from_slice(&(blocks.len() as u32).to_le_bytes());
                for b in blocks.iter() {
                    r.extend_from_slice(&b.0.to_le_bytes());
                }
                for b in blocks.iter() {
                    r.extend_from_slice(&b.1.to_le_bytes());
                }
            }
            r.extend_from_slice(&0u32.to_le_bytes());
            for chunk in s.bases.chunks(4) {
                let mut byte = 0u8;
                for k in 0..4 {
                    let code = match chunk.get(k) {
                        Some(b'C') => 1,
                        Some(b'A') => 2,
                        Some(b'G') => 3,
                        _ => 0,
                    };
                    byte |= code << (6 - 2 * k);
                }
                r.push(byte);
            }
            r
        })
        .collect();

    let mut out = Vec::new();
    for v in [0x1A41_2743u32, version, seqs.len() as u32, 0] {
        out.extend_from_slice(&v.to_le_bytes());
    }
    let mut offset = (TWOBIT_HEADER_LEN + dir_len) as u64;
    for (s, r) in seqs.iter().zip(&records) {
        out.push(s.name.len() as u8);
        out.extend_from_slice(s.name.as_bytes());
        if width == 8 {
            out.extend_from_slice(&offset.to_le_bytes());
        } else {
            out.extend_from_slice(&(offset as u32).to_le_bytes());
        }
        offset += r.len() as u64;
    }
    for r in &records {
        out.extend_from_slice(r);
    }
    out
}
