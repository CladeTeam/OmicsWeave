use super::corpus::*;
use crate::support::rng::Rng;

pub(super) fn word(rng: &mut Rng, alphabet: &[u8], lo: u64, hi: u64) -> Vec<u8> {
    let n = rng.range(lo, hi);
    (0..n).map(|_| *rng.pick(alphabet)).collect()
}

fn number_field(rng: &mut Rng, plain: bool) -> String {
    match rng.range(0, 9) {
        0..=3 => rng.range(0, 3_000_000).to_string(),
        4..=7 => rng.pick(POSITIONS).to_string(),
        _ if plain => "100".into(),
        _ => rng
            .pick(&["", "+7", "2147483648", "4611686018427387904", "-1", "7x"])
            .to_string(),
    }
}

const QUERY_OPS: &[u8] = b"MIS=X";
const REF_OPS: &[u8] = b"DNPH";

/// A CIGAR consuming exactly `len` query bases, as (length, op) pairs.
fn cigar_ops(rng: &mut Rng, len: u64, max_ops: u64) -> Vec<(u64, u8)> {
    let mut ops = Vec::new();
    let mut left = len;
    let n = rng.range(1, max_ops);
    for i in 0..n {
        if left > 0 && (i == n - 1 || rng.chance(70)) {
            let l = if i == n - 1 { left } else { rng.range(1, left) };
            ops.push((l, *rng.pick(QUERY_OPS)));
            left -= l;
        } else {
            ops.push((rng.range(0, 60), *rng.pick(REF_OPS)));
        }
    }
    if left > 0 {
        ops.push((left, b'M'));
    }
    ops
}

pub(super) fn cigar_text(ops: &[(u64, u8)]) -> Vec<u8> {
    let mut s = Vec::new();
    for &(l, op) in ops {
        s.extend_from_slice(l.to_string().as_bytes());
        s.push(op);
    }
    s
}

pub(super) fn cigar_code(op: u8) -> u32 {
    b"MIDNSHP=XB".iter().position(|&c| c == op).unwrap() as u32
}

fn aux(rng: &mut Rng, plain: bool) -> Vec<u8> {
    let id = rng.pick(&["XA", "XB", "NM", "MD", "AS", "RG", "XS", "ZZ", "a1", "tp"]);
    let mut f = format!("{id}:").into_bytes();
    match rng.range(0, if plain { 9 } else { 11 }) {
        0 => {
            f.extend_from_slice(b"A:");
            f.push(*rng.pick(b"AcZ=!~0"));
        }
        1 | 2 => {
            let v = if rng.chance(70) {
                let v = rng.pick(INTS);
                if plain && matches!(*v, "-2147483649" | "4294967296") {
                    "-7".to_string()
                } else {
                    v.to_string()
                }
            } else {
                (rng.range(0, 10_000_000_000) as i64 - 5_000_000_000).to_string()
            };
            f.extend_from_slice(format!("{}:{v}", rng.pick(&["i", "I"])).as_bytes());
        }
        3 => f.extend_from_slice(format!("f:{}", rng.pick(FLOATS)).as_bytes()),
        4 | 5 => {
            f.extend_from_slice(b"Z:");
            f.extend(word(rng, b"ACGT acgt:;=,.-_0123456789", 0, 20));
            if !plain && rng.chance(10) {
                f.extend_from_slice("café".as_bytes());
            }
        }
        6 => {
            f.extend_from_slice(b"H:");
            let n = 2 * rng.range(0, 6);
            f.extend(word(rng, b"0123456789ABCDEFabcdef", n, n));
        }
        7 | 8 if rng.chance(15) => f = rng.pick(RETYPED).as_bytes().to_vec(),
        7 | 8 => {
            let sub = *rng.pick(b"cCsSiIf");
            f.extend_from_slice(format!("B:{}", sub as char).as_bytes());
            let n = rng.range(0, 8);
            for _ in 0..n {
                let v = if sub == b'f' {
                    let v = rng.pick(FLOATS).trim();
                    if plain && matches!(v, "5e" | "12abc") {
                        "2".to_string()
                    } else {
                        v.to_string()
                    }
                } else if plain || rng.chance(70) {
                    let (lo, hi): (i64, i64) = match sub {
                        b'c' => (-128, 127),
                        b'C' => (0, 255),
                        b's' => (-32768, 32767),
                        b'S' => (0, 65535),
                        b'i' => (-2147483648, 2147483647),
                        _ => (0, 4294967295),
                    };
                    let edges = [lo, hi, 0, lo + 1, hi - 1];
                    if rng.chance(50) {
                        rng.pick(&edges).to_string()
                    } else {
                        (lo + (rng.next() % (hi - lo + 1) as u64) as i64).to_string()
                    }
                } else {
                    rng.pick(INTS).to_string()
                };
                f.push(b',');
                f.extend_from_slice(v.as_bytes());
            }
            if !plain && rng.chance(5) {
                f.push(b',');
            }
        }
        9 => f.extend_from_slice(format!("d:{}", rng.pick(FLOATS)).as_bytes()),
        _ => f = rng.pick(BAD_AUX).as_bytes().to_vec(),
    }
    f
}

/// One data line. `plain` lines stay clear of deliberately malformed
/// fields, and so are what htslib's format detection recognises as SAM, as
/// the first line of a headerless file must be.
pub(super) fn record(rng: &mut Rng, names: &[String], plain: bool) -> (Vec<u8>, bool) {
    let plain = plain || rng.chance(60);
    let mut unknown = false;
    let mut f: Vec<Vec<u8>> = Vec::new();
    // QNAME
    f.push(if plain || rng.chance(95) {
        let mut q = word(
            rng,
            b"ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789",
            1,
            1,
        );
        q.extend(word(rng, b"ABCabc0123456789_:.-/#", 0, 24));
        q
    } else {
        rng.pick(&[String::new(), "q".repeat(254), "q".repeat(255)])
            .as_bytes()
            .to_vec()
    });
    // FLAG
    let flag = *rng.pick(&[
        0u32, 4, 16, 83, 99, 147, 163, 256, 2048, 0x900, 0xffff, 1, 8, 0x40,
    ]);
    let flag = if rng.chance(20) {
        rng.range(0, 65535) as u32
    } else {
        flag
    };
    f.push(
        match rng.range(0, if plain { 0 } else { 19 }) {
            0..=13 => flag.to_string(),
            14 => format!("0x{flag:x}"),
            15 => format!("0{flag:o}"),
            _ => rng
                .pick(&["", "65536", "0x10000", "1a", "-1", "0x", "08"])
                .to_string(),
        }
        .into_bytes(),
    );
    // RNAME
    let rname = if !plain && rng.chance(3) {
        unknown = true;
        "chrUnknown".to_string()
    } else if !names.is_empty() && rng.chance(85) {
        rng.pick(names).clone()
    } else {
        "*".to_string()
    };
    let placed = rname != "*";
    f.push(rname.into_bytes());
    // POS
    let pos = number_field(rng, plain);
    f.push(pos.clone().into_bytes());
    // MAPQ
    f.push(
        match rng.range(0, if plain { 3 } else { 5 }) {
            0..=3 => rng.range(0, 255).to_string(),
            _ => rng.pick(&["", "256", "60x"]).to_string(),
        }
        .into_bytes(),
    );
    // CIGAR, SEQ, QUAL
    let len = if rng.chance(10) { 0 } else { rng.range(1, 40) };
    let mut cg: Option<Vec<u8>> = None;
    let cigar = match rng.range(0, 99) {
        0..=14 => b"*".to_vec(),
        15..=79 => cigar_text(&cigar_ops(rng, len, 6)),
        80..=84 => {
            let wrong = len + rng.range(1, 3);
            cigar_text(&cigar_ops(rng, wrong, 4))
        }
        85..=91 if len > 0 && placed && pos != "0" => {
            // A CIGAR stored in a CG tag behind a soft-clip placeholder.
            let real = cigar_ops(rng, len, 8);
            let mut tag = format!("CG:B:{},", rng.pick(&["I", "I", "i"])).into_bytes();
            let codes: Vec<String> = real
                .iter()
                .map(|&(l, op)| ((l as u32) << 4 | cigar_code(op)).to_string())
                .collect();
            tag.extend_from_slice(codes.join(",").as_bytes());
            cg = Some(tag);
            if rng.chance(70) {
                format!("{len}S").into_bytes()
            } else {
                format!("{len}S{}N", rng.range(1, 100)).into_bytes()
            }
        }
        _ if plain => cigar_text(&cigar_ops(rng, len, 6)),
        _ => rng
            .pick(&[
                "10M5",
                "M",
                "0M",
                "268435455N",
                "268435456M",
                "10m",
                "4B",
                "+4M",
            ])
            .as_bytes()
            .to_vec(),
    };
    f.push(cigar);
    // RNEXT
    f.push(if !plain && rng.chance(2) {
        unknown = true;
        b"chrMate".to_vec()
    } else {
        match rng.range(0, 2) {
            0 => b"=".to_vec(),
            1 => b"*".to_vec(),
            _ if names.is_empty() => b"*".to_vec(),
            _ => rng.pick(names).clone().into_bytes(),
        }
    });
    // PNEXT, TLEN
    f.push(number_field(rng, plain).into_bytes());
    let tlens: &[&str] = &["0", "-1", "150", "-150", "+5", "2147483647", "-2147483648"];
    let wild: &[&str] = &[
        "2147483648",
        "-2147483649",
        "4611686018427387903",
        "-4611686018427387904",
        "4611686018427387904",
    ];
    let tlens = if plain || rng.chance(50) { tlens } else { wild };
    f.push(rng.pick(tlens).as_bytes().to_vec());
    let seq = if len == 0 && rng.chance(80) {
        b"*".to_vec()
    } else {
        word(rng, b"ACGTNacgtn=MRWSYKVHDB.*", len, len)
    };
    let l_seq = if seq == b"*" { 0 } else { seq.len() as u64 };
    f.push(seq);
    f.push(match rng.range(0, if plain { 89 } else { 99 }) {
        0..=79 => word(rng, &(33u8..=126).collect::<Vec<_>>(), l_seq, l_seq),
        80..=89 => b"*".to_vec(),
        90..=92 => word(rng, b"I", l_seq + 1, l_seq + 1),
        93 => vec![b' '; l_seq as usize],
        94 => vec![0xa0; l_seq as usize],
        95 => vec![0xa1; l_seq as usize],
        _ => Vec::new(),
    });
    for _ in 0..rng.range(0, 6) {
        f.push(aux(rng, plain));
    }
    if let Some(tag) = cg {
        let at = rng.range(11, f.len() as u64) as usize;
        f.insert(at, tag);
    }
    (f.join(&b'\t'), unknown)
}
