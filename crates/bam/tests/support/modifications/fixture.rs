use crate::support::Rng;
use std::path::Path;

const LETTERS: &[u8] = b"mhfcaogebnN";

/// One random `MM` / `ML` / `MN` tag set for a read with `seq`, mostly
/// well formed: deltas within the read or just past it, ML of the right
/// length or one off, and the odd malformed entry.
fn random_tags(rng: &mut Rng, seq: &[u8], reverse: bool) -> String {
    let original: Vec<u8> = if reverse {
        seq.iter()
            .rev()
            .map(|&b| match b {
                b'A' => b'T',
                b'C' => b'G',
                b'G' => b'C',
                b'T' => b'A',
                x => x,
            })
            .collect()
    } else {
        seq.to_vec()
    };
    let mut mm = String::new();
    let mut n_ml = 0usize;
    for _ in 0..rng.range(0, 4) {
        let base = if rng.chance(2) {
            b'X'
        } else {
            b"ACGTNUCCA"[rng.range(0, 9) as usize]
        };
        mm.push(base as char);
        mm.push(if rng.chance(1) {
            '*'
        } else if rng.chance(25) {
            '-'
        } else {
            '+'
        });
        let n_codes = if rng.chance(15) {
            let chebi = [76792u64, 17596, 0, 27551][rng.range(0, 4) as usize];
            mm.push_str(&chebi.to_string());
            1
        } else {
            let n = rng.range(1, 4) as usize;
            for _ in 0..n {
                mm.push(LETTERS[rng.range(0, LETTERS.len() as u64) as usize] as char);
            }
            n
        };
        let marker = ["", ".", "?"][rng.range(0, 3) as usize];
        // Canonical bases of the read as sequenced.
        let canonical = match base {
            b'N' => original.len(),
            b'U' => original.iter().filter(|&&b| b == b'T').count(),
            b => original.iter().filter(|&&x| x == b).count(),
        };
        let mut left = canonical;
        let mut n_deltas = 0;
        let mut deltas = String::new();
        while left > 0 && rng.chance(80) {
            let d = rng.range(0, (left as u64).min(6)) as usize;
            deltas.push_str(&format!(",{d}"));
            n_deltas += 1;
            left -= d + 1;
        }
        if rng.chance(4) {
            deltas.push_str(&format!(",{}", rng.range(0, 4)));
            n_deltas += 1;
        }
        // htslib misreads a reverse read's ChEBI entry with a mark and no
        // deltas when there is no ML (it resumes at the mark); the decoder
        // reads it as the specification does, so it is left out here.
        if !(n_codes == 1 && mm.ends_with(|c: char| c.is_ascii_digit()) && n_deltas == 0) {
            mm.push_str(marker);
        }
        mm.push_str(&deltas);
        n_ml += n_deltas * n_codes;
        if !rng.chance(2) {
            mm.push(';');
        }
    }
    let legacy = rng.chance(10);
    let mut tags = format!("\t{}:Z:{mm}", if legacy { "Mm" } else { "MM" });
    if rng.chance(85) {
        let n = match rng.range(0, 20) {
            0 => n_ml + 1,
            1 if n_ml > 0 => n_ml - 1,
            _ => n_ml,
        };
        if n > 0 {
            let values: Vec<String> = (0..n).map(|_| rng.range(0, 256).to_string()).collect();
            let kind = if rng.chance(2) { 'S' } else { 'C' };
            tags.push_str(&format!(
                "\t{}:B:{kind},{}",
                if legacy { "Ml" } else { "ML" },
                values.join(",")
            ));
        }
    }
    if rng.chance(20) {
        let n = if rng.chance(20) {
            seq.len() + 1
        } else {
            seq.len()
        };
        tags.push_str(&format!("\tMN:i:{n}"));
    }
    tags
}

/// A CIGAR of `len` query bases with clips, insertions and deletions.
fn random_cigar(rng: &mut Rng, len: usize) -> String {
    if len == 0 {
        return "5M".into();
    }
    let mut ops = String::new();
    if rng.chance(10) {
        ops.push_str(&format!("{}H", rng.range(1, 20)));
    }
    let mut left = len;
    if left > 2 && rng.chance(30) {
        let s = rng.range(1, (left as u64 / 3).max(2)) as usize;
        ops.push_str(&format!("{s}S"));
        left -= s;
    }
    while left > 0 {
        let n = rng.range(1, (left as u64).min(30) + 1) as usize;
        let op = match rng.range(0, 10) {
            0 => 'I',
            1 => '=',
            2 => 'X',
            _ => 'M',
        };
        ops.push_str(&format!("{n}{op}"));
        left -= n;
        if left > 0 && rng.chance(20) {
            let d = if rng.chance(50) { 'D' } else { 'N' };
            ops.push_str(&format!("{}{d}", rng.range(1, 50)));
        }
    }
    if rng.chance(10) {
        ops.push_str(&format!("{}H", rng.range(1, 20)));
    }
    ops
}

/// A sorted SAM of `n` random reads with modification tags.
pub fn random_sam(seed: u64, n: usize) -> Vec<String> {
    let mut rng = Rng::new(seed);
    let mut lines = Vec::new();
    let mut pos = 1u64;
    for i in 0..n {
        pos += rng.range(0, 40);
        let len = if rng.chance(5) {
            0
        } else {
            rng.range(1, 150) as usize
        };
        let seq: Vec<u8> = (0..len)
            .map(|_| {
                if rng.chance(2) {
                    b"MRSWYKN"[rng.range(0, 7) as usize]
                } else {
                    b"ACGTACGTN"[rng.range(0, 9) as usize]
                }
            })
            .collect();
        let reverse = rng.chance(40);
        let unmapped = rng.chance(8);
        let flag = ((reverse as u16) * 0x10) | ((unmapped as u16) * 0x4);
        let cigar = if unmapped {
            "*".to_string()
        } else {
            random_cigar(&mut rng, len)
        };
        let (seq_text, qual) = if len == 0 {
            ("*".to_string(), "*".to_string())
        } else {
            (String::from_utf8(seq.clone()).unwrap(), "I".repeat(len))
        };
        let tags = if rng.chance(5) {
            String::new()
        } else {
            random_tags(&mut rng, &seq, reverse)
        };
        lines.push(format!(
            "r{i}\t{flag}\tchr1\t{pos}\t60\t{cigar}\t*\t0\t0\t{seq_text}\t{qual}{tags}"
        ));
    }
    lines
}

pub fn write_bam(dir: &Path, name: &str, lines: &[String]) -> std::path::PathBuf {
    let sam = dir.join(format!("{name}.sam"));
    let mut text = String::from("@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:chr1\tLN:1000000\n");
    for l in lines {
        text.push_str(l);
        text.push('\n');
    }
    std::fs::write(&sam, text).unwrap();
    let bam = dir.join(format!("{name}.bam"));
    crate::support::bam_of(&sam, &bam, 6);
    bam
}

/// The SAMtags §1.7 examples and the shapes of htslib's own base-mods
/// tests: interleaved codes, ChEBI codes, both strands of a duplex read,
/// explicit and implicit entries, reverse reads, `N` as the canonical
/// base, hard clips with `MN`.
pub fn spec_lines() -> Vec<String> {
    let s = "AGCTCTCCAGAGTCGNACGCCATYCGCGCGCCACCA";
    let rc: String = s
        .chars()
        .rev()
        .map(|b| match b {
            'A' => 'T',
            'C' => 'G',
            'G' => 'C',
            'T' => 'A',
            'Y' => 'R',
            x => x,
        })
        .collect();
    let q = "I".repeat(s.len());
    let n = s.len();
    let mut lines = Vec::new();
    let tagsets = [
        "\tMM:Z:C+m,2,2,1,4,1;\tML:B:C,102,128,153,179,204",
        "\tMM:Z:C+mh,2,2,0,0,4,1;\tML:B:C,77,159,103,133,128,108,154,82,179,57,204,25",
        "\tMM:Z:C+m,2,2,1,4,1;C+h,6,7;\tML:B:C,128,153,179,204,230,159,6",
        "\tMM:Z:C+m,2,2,1,4,1;C+76792,6,7;N+n,15,2;\tML:B:C,128,153,179,204,230,159,6,215,240",
        "\tMM:Z:C+m?,2,2,1;C+h.,6,7;\tML:B:C,128,153,179,159,6",
        "\tMM:Z:G-m,0,0,0;C+m,0,1,0;\tML:B:C,10,20,30,40,50,60",
        "\tMM:Z:C+m,2,2,1,4,1;",
        "\tMm:Z:C+m,2,2;\tMl:B:C,1,2,3",
        "\tMM:Z:A+a.,0,1;T-a,1;\tML:B:C,200,201,202",
        "\tMM:Z:C+m;",
        "\tMM:Z:;",
    ];
    for (i, tags) in tagsets.iter().enumerate() {
        let pos = 100 + i * 10;
        lines.push(format!(
            "f{i}\t0\tchr1\t{pos}\t60\t3S10M2I5M4D16M\t*\t0\t0\t{s}\t{q}{tags}"
        ));
        lines.push(format!(
            "r{i}\t16\tchr1\t{}\t60\t5H{n}M\t*\t0\t0\t{rc}\t{q}{tags}\tMN:i:{n}",
            pos + 1
        ));
    }
    lines
}
