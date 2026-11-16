//! Scalar sequence expectations with an explicit complement alphabet.
use omicsweave_core::{Span, Strand};

pub fn oriented(stored: &[u8], spans: &[Span], strand: Strand) -> Vec<u8> {
    let mut out: Vec<_> = spans
        .iter()
        .flat_map(|s| stored[s.start as usize..s.end as usize].iter().copied())
        .collect();
    if strand == Strand::Reverse {
        out.reverse();
        for byte in &mut out {
            let alphabet = b"ACGTUacgtuRYKMBVDHSWNrykmbvdhswn";
            let paired = b"TGCAAtgcaaYRMKVBHDSWNyrmkvbhdswn";
            if let Some(at) = alphabet.iter().position(|b| b == byte) {
                *byte = paired[at];
            }
        }
    }
    out
}

pub fn twobit_stored(bases: &[u8], mask: bool) -> Vec<u8> {
    bases
        .iter()
        .map(|&b| {
            if b"ACGTacgt".contains(&b) {
                if mask { b } else { b.to_ascii_uppercase() }
            } else {
                b'N'
            }
        })
        .collect()
}
