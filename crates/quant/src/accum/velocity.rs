//! RNA velocity: spliced, unspliced and ambiguous molecules per gene, as
//! layers of a gene matrix ([`crate::Counts::velocity`]).
//!
//! Rules (velocyto / STARsolo `Velocyto`, per read, then per molecule):
//! - A read belongs to a gene when its aligned blocks overlap the span of
//!   exactly one gene (introns included) on the strand the library gives;
//!   reads overlapping several genes are not counted (velocyto's
//!   multi-gene rule).
//! - Against each transcript of that gene whose span holds every block, the
//!   read is **exonic** when every block lies within one exon, and
//!   **intronic** otherwise (a block in an intron, or across an exon-intron
//!   boundary). Transcripts that do not hold the read are not consulted; a
//!   read no transcript holds is not counted.
//! - A read exonic against every transcript consulted is spliced, intronic
//!   against every one unspliced, and ambiguous when the transcripts
//!   disagree.
//! - A molecule (a UMI of a barcode and gene; without UMIs, a read, or a
//!   read pair by name) is spliced when all its reads are, unspliced when
//!   all are, and ambiguous otherwise.
//!
//! Each read contributes evidence bits ([`EXONIC`], [`INTRONIC`]) that the
//! matrix ORs per molecule and classifies when it is read.

use omicsweave_gtf::Annotation;

/// Evidence that a read is exonic against some transcript.
pub const EXONIC: u32 = 1;
/// Evidence that a read is intronic against some transcript.
pub const INTRONIC: u32 = 2;

/// Layer offset (from the `spliced` layer) of a molecule's evidence bits.
#[inline]
pub fn class(bits: u32) -> usize {
    match bits & (EXONIC | INTRONIC) {
        EXONIC => 0,
        INTRONIC => 1,
        _ => 2,
    }
}

/// The transcripts of each gene column: their spans and merged exons.
#[derive(Debug, Clone, Default)]
pub struct GeneModels {
    /// `transcripts[offsets[c]..offsets[c + 1]]` are column `c`'s.
    offsets: Vec<u32>,
    /// `(start, end, first exon, end of exons)`.
    transcripts: Vec<(u32, u32, u32, u32)>,
    exons: Vec<(u32, u32)>,
}

impl GeneModels {
    /// The models of `genes` (annotation gene rows, one per column, in
    /// column order).
    pub fn new(ann: &Annotation, genes: &[u32]) -> GeneModels {
        let mut m = GeneModels {
            offsets: vec![0],
            ..GeneModels::default()
        };
        let mut blocks = Vec::new();
        for &g in genes {
            for t in ann.genes().transcripts(g) {
                ann.exon_blocks(t, &mut blocks);
                if let (Some(first), Some(last)) = (blocks.first(), blocks.last()) {
                    let lo = m.exons.len() as u32;
                    m.exons.extend_from_slice(&blocks);
                    m.transcripts
                        .push((first.0, last.1, lo, m.exons.len() as u32));
                }
            }
            m.offsets.push(m.transcripts.len() as u32);
        }
        m
    }

    /// The evidence of a read with aligned `blocks` against the transcripts
    /// of column `col`: [`EXONIC`] and / or [`INTRONIC`]; 0 when no
    /// transcript holds it.
    pub fn classify(&self, col: u32, blocks: &[(u32, u32)]) -> u32 {
        let (Some(first), Some(last)) = (blocks.first(), blocks.last()) else {
            return 0;
        };
        let c = col as usize;
        let mut bits = 0;
        for &(ts, te, lo, hi) in
            &self.transcripts[self.offsets[c] as usize..self.offsets[c + 1] as usize]
        {
            if first.0 < ts || last.1 > te {
                continue;
            }
            let exons = &self.exons[lo as usize..hi as usize];
            let exonic = blocks.iter().all(|&(s, e)| {
                let i = exons.partition_point(|x| x.1 <= s);
                i < exons.len() && exons[i].0 <= s && e <= exons[i].1
            });
            bits |= if exonic { EXONIC } else { INTRONIC };
            // Both classes are terminal: additional transcripts cannot
            // change an ambiguous molecule's evidence.
            if bits == EXONIC | INTRONIC {
                break;
            }
        }
        bits
    }

    pub fn heap_bytes(&self) -> u64 {
        (4 * self.offsets.capacity() + 16 * self.transcripts.capacity() + 8 * self.exons.capacity())
            as u64
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn models() -> GeneModels {
        // One gene, two transcripts: exons [100,200) [300,400), and
        // [100,400) (the intron retained).
        GeneModels {
            offsets: vec![0, 2],
            transcripts: vec![(100, 400, 0, 2), (100, 400, 2, 3)],
            exons: vec![(100, 200), (300, 400), (100, 400)],
        }
    }

    #[test]
    fn reads_are_classified_against_every_transcript_holding_them() {
        let m = models();
        // Spliced across the intron: exonic in both.
        assert_eq!(m.classify(0, &[(150, 200), (300, 350)]), EXONIC);
        // In the intron: intronic in one, exonic in the other.
        assert_eq!(m.classify(0, &[(220, 260)]), EXONIC | INTRONIC);
        // Past the transcripts: held by none.
        assert_eq!(m.classify(0, &[(390, 420)]), 0);
        assert_eq!(class(EXONIC), 0);
        assert_eq!(class(INTRONIC), 1);
        assert_eq!(class(EXONIC | INTRONIC), 2);
    }
}
