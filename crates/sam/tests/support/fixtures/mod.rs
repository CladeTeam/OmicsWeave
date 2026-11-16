mod corpus;
mod header;
mod record;

use super::rng::Rng;
use header::header;
use record::{cigar_code, cigar_text, record, word};

/// A random SAM text and, per data line, whether it names a reference the
/// header lacks (which this crate rejects and htslib does not).
pub struct Fixture {
    pub text: Vec<u8>,
    pub unknown: Vec<bool>,
}

/// A random SAM text: possibly no header, references with alternative
/// names, all record fields at their boundaries, and malformed fields.
pub fn random_sam(rng: &mut Rng) -> Fixture {
    let mut names = Vec::new();
    let header = if rng.chance(15) {
        Vec::new()
    } else {
        header(rng, &mut names)
    };
    let eol: &[u8] = if rng.chance(20) { b"\r\n" } else { b"\n" };
    let mut text = Vec::new();
    for l in &header {
        text.extend_from_slice(l);
        text.extend_from_slice(eol);
    }
    let mut unknown = Vec::new();
    let n = rng.range(if header.is_empty() { 1 } else { 0 }, 12);
    for i in 0..n {
        let (line, u) = record(rng, &names, header.is_empty() && i == 0);
        text.extend_from_slice(&line);
        text.extend_from_slice(eol);
        unknown.push(u);
    }
    Fixture { text, unknown }
}

/// SAM with records whose CIGARs exceed BAM's 65535 operations, given as
/// CIGAR text or as a CG tag behind a soft-clip placeholder.
pub fn long_cigar_sam(rng: &mut Rng) -> Fixture {
    let mut text = b"@HD\tVN:1.6\n@SQ\tSN:chr1\tLN:248956422\n".to_vec();
    let mut unknown = Vec::new();
    for _ in 0..rng.range(1, 3) {
        let n_ops = rng.range(65_530, 66_600);
        let mut ops: Vec<(u64, u8)> = (0..n_ops)
            .map(|i| {
                let op = if i % 2 == 0 {
                    b'M'
                } else {
                    *rng.pick(b"IDN=X")
                };
                (rng.range(1, 3), op)
            })
            .collect();
        if rng.chance(30) {
            ops.insert(0, (rng.range(1, 9), b'S'));
        }
        let qlen: u64 = ops
            .iter()
            .filter(|(_, op)| b"MIS=X".contains(op))
            .map(|(l, _)| l)
            .sum();
        let flag = *rng.pick(&[0, 16, 4, 99]);
        let seq = word(rng, b"ACGT", qlen, qlen);
        let qual = if rng.chance(50) {
            b"*".to_vec()
        } else {
            word(rng, b"#5?I", qlen, qlen)
        };
        let pos = rng.pick(&["1", "1000", "100000000"]).to_string();
        let mut line = format!("long{}\t{flag}\tchr1\t{pos}\t60\t", rng.range(0, 999)).into_bytes();
        let as_tag = rng.chance(40);
        if as_tag {
            line.extend_from_slice(format!("{qlen}S").as_bytes());
        } else {
            line.extend_from_slice(&cigar_text(&ops));
        }
        line.extend_from_slice(b"\t=\t1\t0\t");
        line.extend_from_slice(&seq);
        line.push(b'\t');
        line.extend_from_slice(&qual);
        line.extend_from_slice(b"\tNM:i:3");
        if as_tag {
            line.extend_from_slice(b"\tCG:B:I");
            for &(l, op) in &ops {
                line.extend_from_slice(format!(",{}", (l as u32) << 4 | cigar_code(op)).as_bytes());
            }
        }
        line.extend_from_slice(b"\tXZ:Z:after\n");
        text.extend_from_slice(&line);
        unknown.push(false);
    }
    Fixture { text, unknown }
}

pub fn cg_line(fake_count: usize, real_count: usize, before: usize, after: usize) -> Vec<u8> {
    let mut line = b"q\t0\tchr1\t1\t60\t8S".to_vec();
    for _ in 1..fake_count {
        line.extend_from_slice(b"0N");
    }
    line.extend_from_slice(b"\t*\t0\t0\tACGTACGT\tIIIIIIII\tAA:Z:");
    line.extend(std::iter::repeat_n(b'a', before));
    line.extend_from_slice(b"\tCG:B:I,128");
    for _ in 1..real_count {
        line.extend_from_slice(b",2");
    }
    line.extend_from_slice(b"\tCG:Z:shadow\tZZ:Z:");
    line.extend(std::iter::repeat_n(b'z', after));
    line
}
