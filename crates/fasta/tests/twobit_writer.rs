//! 2bit written by `TwoBitWriter`, piece by piece: byte-identical to the
//! test writer's `faToTwoBit` layout and read back by lib2bit.

mod common;

use common::{Rng, TwoBitSeq, lib2bit, write_twobit};
use omicsweave_fasta::{TwoBitPiece, TwoBitRecord, TwoBitWriter};
use proptest::prelude::*;
use std::io::Cursor;

fn random_seq(rng: &mut Rng) -> Vec<u8> {
    let len = if rng.chance(5) {
        0
    } else {
        rng.range(1, 6_000) as usize
    };
    let mut seq = Vec::with_capacity(len);
    while seq.len() < len {
        let run = (rng.range(1, 300) as usize).min(len - seq.len());
        let alphabet: &[u8] = match rng.range(0, 6) {
            0 => b"NNNN",
            1 => b"nnnn",
            2 => b"acgt",
            3 => b"ACGTRYKMn-",
            _ => b"ACGT",
        };
        seq.extend((0..run).map(|_| *rng.pick(alphabet)));
    }
    seq
}

fn runs(seq: &[u8], hit: impl Fn(u8) -> bool) -> Vec<(u32, u32)> {
    let mut out = Vec::new();
    let mut i = 0;
    while i < seq.len() {
        if hit(seq[i]) {
            let s = i;
            while i < seq.len() && hit(seq[i]) {
                i += 1;
            }
            out.push((s as u32, (i - s) as u32));
        } else {
            i += 1;
        }
    }
    out
}

fn is_acgt(c: u8) -> bool {
    b"ACGTacgt".contains(&c)
}

proptest! {
    #![proptest_config(ProptestConfig { cases: 48, ..ProptestConfig::default() })]

    #[test]
    fn matches_the_reference_layout_and_lib2bit(seed in any::<u64>(), version in 0u32..=1) {
        let mut rng = Rng::new(seed);
        let seqs: Vec<(String, Vec<u8>)> = (0..rng.range(1, 6))
            .map(|i| (format!("s{i}"), random_seq(&mut rng)))
            .collect();
        let piece = 4 * rng.range(1, 400) as usize;

        let names = seqs.iter().map(|(n, _)| n.clone()).collect();
        let mut w = TwoBitWriter::create(Cursor::new(Vec::new()), names, version).unwrap();
        for (name, seq) in &seqs {
            let mut rec = TwoBitRecord::new();
            let mut packed = Vec::new();
            for (i, chunk) in seq.chunks(piece).enumerate() {
                let p = TwoBitPiece::encode((i * piece) as u64, chunk).unwrap();
                rec.push(&p).unwrap();
                packed.extend_from_slice(&p.packed);
            }
            w.add(name, &rec, &mut &packed[..]).unwrap();
        }
        let ours = w.finish().unwrap().into_inner();

        let reference: Vec<TwoBitSeq> = seqs
            .iter()
            .map(|(name, seq)| TwoBitSeq {
                name: name.clone(),
                bases: seq.iter().map(|&c| if is_acgt(c) { c.to_ascii_uppercase() } else { b'T' }).collect(),
                n_blocks: runs(seq, |c| !is_acgt(c)),
                mask_blocks: runs(seq, |c| c.is_ascii_lowercase()),
            })
            .collect();
        prop_assert_eq!(&ours, &write_twobit(&reference, version));

        let dir = tempfile::tempdir().unwrap();
        let path = dir.path().join("t.2bit");
        std::fs::write(&path, &ours).unwrap();
        for mask in [false, true] {
            let theirs = lib2bit::TwoBit::open(&path, mask).unwrap();
            for (name, seq) in &seqs {
                prop_assert_eq!(theirs.len(name), seq.len() as u32);
                if seq.is_empty() {
                    continue;
                }
                let got = theirs.sequence(name, 0, seq.len() as u32).unwrap();
                for (i, (&g, &c)) in got.iter().zip(seq).enumerate() {
                    if is_acgt(c) {
                        let want = if mask { c } else { c.to_ascii_uppercase() };
                        prop_assert_eq!(g, want, "{}:{}", name, i);
                    } else {
                        prop_assert_eq!(g.to_ascii_uppercase(), b'N', "{}:{}", name, i);
                    }
                }
            }
        }
    }
}
