//! The UCSC chain from reference to consensus (`-c`; consensus.c
//! `init_chain`, `push_chain_gap`, `print_chain`).

use std::fmt::Write;

/// Ungapped blocks and the gaps after them, on reference and consensus.
pub(crate) struct Chain {
    /// Per block: its length and the reference and consensus gaps after it.
    blocks: Vec<[i64; 3]>,
    ori_pos: i64,
    ref_last: i64,
    alt_last: i64,
    /// The region's reference length.
    length: i64,
}

impl Chain {
    pub fn new(ori_pos: i64, length: i64) -> Chain {
        Chain {
            blocks: Vec::new(),
            ori_pos,
            ref_last: ori_pos,
            alt_last: ori_pos,
            length,
        }
    }

    /// Starts another query while retaining space for its gap blocks.
    pub fn reset(&mut self, ori_pos: i64, length: i64) {
        self.blocks.clear();
        self.ori_pos = ori_pos;
        self.ref_last = ori_pos;
        self.alt_last = ori_pos;
        self.length = length;
    }

    /// Records a gap; one that starts where the last ended extends it.
    pub fn push_gap(&mut self, ref_start: i64, ref_len: i64, alt_start: i64, alt_len: i64) {
        match self.blocks.last_mut() {
            Some(last) if ref_start <= self.ref_last => {
                last[1] += ref_len;
                last[2] += alt_len;
            }
            _ => self
                .blocks
                .push([ref_start - self.ref_last, ref_len, alt_len]),
        }
        self.ref_last = ref_start + ref_len;
        self.alt_last = alt_start + alt_len;
    }

    /// The chain text of the region's sequence `chrom`, chain id 1.
    pub fn render(&self, chrom: &str, out: &mut String) {
        let ref_end = self.length + self.ori_pos;
        let tail = ref_end - self.ref_last;
        let alt_end = self.alt_last + tail;
        let score: i64 = self.blocks.iter().map(|b| b[0]).sum::<i64>() + tail;
        let ori = self.ori_pos;
        out.clear();
        let _ = writeln!(
            out,
            "chain {score} {chrom} {ref_end} + {ori} {ref_end} {chrom} {alt_end} + {ori} {alt_end} 1"
        );
        for [block, ref_gap, alt_gap] in &self.blocks {
            let _ = writeln!(out, "{block} {ref_gap} {alt_gap}");
        }
        let _ = write!(out, "{tail}\n\n");
    }
}
