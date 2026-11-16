//! Per-base expectations computed directly from parsed feature columns.

use super::LEN;
use omicsweave_gtf::{Annotation, FeatureKind, IntervalSet, Strand};

/// Exonic bases from feature rows, without calling any derived view.
pub fn exon_bits(a: &Annotation, t: u32) -> Vec<bool> {
    let f = a.features();
    let mut bits = vec![false; LEN as usize];
    for i in a.transcripts().features(t) {
        let i = i as usize;
        if f.kind[i] == FeatureKind::Exon {
            bits[f.start[i] as usize..f.end[i] as usize].fill(true);
        }
    }
    bits
}

/// Maximal runs of `true`.
pub fn runs(bits: &[bool]) -> Vec<(u32, u32)> {
    let mut out = Vec::new();
    let mut i = 0;
    while i < bits.len() {
        if bits[i] {
            let s = i;
            while i < bits.len() && bits[i] {
                i += 1;
            }
            out.push((s as u32, i as u32));
        } else {
            i += 1;
        }
    }
    out
}

pub fn rows(set: &IntervalSet) -> Vec<(u32, u32, u32, Strand, u32)> {
    (0..set.len())
        .map(|i| {
            (
                set.seqid[i],
                set.start[i],
                set.end[i],
                set.strand[i],
                set.owner[i],
            )
        })
        .collect()
}

pub fn naive_introns(a: &Annotation, t: u32) -> Vec<(u32, u32)> {
    let r = runs(&exon_bits(a, t));
    r.windows(2).map(|w| (w[0].1, w[1].0)).collect()
}
