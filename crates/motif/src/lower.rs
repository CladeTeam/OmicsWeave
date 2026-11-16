//! Step ⑥: an integer model becomes lookup blocks over packed k-mers.

use omicsweave_core::{Error, Result};

use crate::score::{FORBIDDEN, IntModel};

/// How symbols are packed for the scan: `bits` per symbol, `k` symbols per
/// k-mer code, so a block table has at most `2^(bits·k)` cells.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) struct Packing {
    pub n: usize,
    pub bits: u32,
    pub k: usize,
}

impl Packing {
    /// Cells one block table may hold.
    const TABLE_CELLS_LOG2: u32 = 10;

    pub fn new(alphabet_len: usize) -> Packing {
        let bits = (usize::BITS - (alphabet_len.max(2) - 1).leading_zeros()).max(1);
        let k = (Packing::TABLE_CELLS_LOG2 / bits).max(1) as usize;
        Packing {
            n: alphabet_len,
            bits,
            k,
        }
    }

    pub fn mask(&self) -> u32 {
        ((1u64 << (self.bits as usize * self.k)) - 1) as u32
    }
}

#[derive(Clone, Debug)]
pub(crate) enum Reader {
    /// The first `len` symbols of the k-mer at the block's offset:
    /// `kmer >> shift`.
    Kmer { shift: u32 },
    /// Symbols at arbitrary positions relative to the window, base `n`.
    Sparse { positions: Vec<i32> },
}

#[derive(Clone, Debug)]
pub(crate) struct Block {
    pub offset: i32,
    pub reader: Reader,
    pub table: Vec<i32>,
    pub max: i32,
}

/// Blocks in evaluation order with the bound used to prune.
#[derive(Clone, Debug)]
pub(crate) struct Lowered {
    pub blocks: Vec<Block>,
    /// `rest_max[j]`: the largest sum blocks `j..` can add.
    pub rest_max: Vec<i32>,
    /// Width of the window (positions that must be known).
    pub width: u32,
    /// Leftmost position read.
    pub lo: i32,
}

/// Lowers `model`. Every term lands in exactly one block, so a window's
/// block sum equals its integer score; a forbidden cell is stored as a value
/// negative enough that any window containing it sums below zero.
pub(crate) fn lower(model: &IntModel, packing: &Packing) -> Result<Lowered> {
    let k = packing.k as i32;
    let w = model.width as i32;
    let lo = model.lo();
    let n = packing.n;
    // (offset, length, terms) of k-mer blocks; sparse blocks per term.
    let mut kmer_blocks: Vec<(i32, i32, Vec<usize>)> = Vec::new();
    let mut sparse: Vec<usize> = Vec::new();
    if w - lo <= k {
        kmer_blocks.push((lo, w - lo, (0..model.terms.len()).collect()));
    } else {
        // Stab each term's interval of feasible offsets, greedily by its
        // right end: the fewest blocks that hold every term.
        let mut spans: Vec<(i32, i32, usize)> = Vec::new();
        for (i, t) in model.terms.iter().enumerate() {
            let a = t.positions[0];
            let b = *t.positions.last().unwrap();
            if b - a + 1 > k {
                sparse.push(i);
                continue;
            }
            spans.push(((b - k + 1).max(lo), a.min(w - k), i));
        }
        spans.sort_by_key(|&(_, upper, i)| (upper, i));
        for (lower, upper, i) in spans {
            // Chosen offsets increase with the sorted upper endpoints. The
            // first offset at least lower is the first greedy block that
            // contains this interval; all existing offsets are <= upper.
            let block = kmer_blocks.partition_point(|&(offset, _, _)| offset < lower);
            if block == kmer_blocks.len() {
                kmer_blocks.push((upper, k, vec![i]));
            } else {
                kmer_blocks[block].2.push(i);
            }
        }
    }
    let forbid = -(model.max as i64 + 1);
    let block_count = kmer_blocks.len() + sparse.len();
    if (block_count as i64 + 1) * (model.max as i64 + 1) > i32::MAX as i64 / 2 {
        return Err(Error::invalid(
            "motif scan",
            "integer scores too large for this many blocks",
        ));
    }
    let forbid = forbid as i32;
    let mut blocks = Vec::with_capacity(block_count);
    let bits = packing.bits as usize;
    let sym_mask = (1usize << bits) - 1;
    for (offset, len, members) in kmer_blocks {
        let len = len as usize;
        let size = 1usize << (bits * len);
        let mut table = vec![0i32; size];
        let mut digits = vec![0usize; len];
        'cells: for (idx, cell) in table.iter_mut().enumerate() {
            for (j, d) in digits.iter_mut().enumerate() {
                *d = (idx >> (bits * (len - 1 - j))) & sym_mask;
                if *d >= n {
                    continue 'cells;
                }
            }
            let mut sum = 0i64;
            for &ti in &members {
                let t = &model.terms[ti];
                let ci = t
                    .positions
                    .iter()
                    .fold(0usize, |acc, &p| acc * n + digits[(p - offset) as usize]);
                let v = t.table[ci];
                if v == FORBIDDEN {
                    sum = forbid as i64;
                    break;
                }
                sum += v as i64;
            }
            *cell = sum as i32;
        }
        let shift = (bits * (packing.k - len)) as u32;
        blocks.push(finish(offset, Reader::Kmer { shift }, table, forbid));
    }
    for ti in sparse {
        let t = &model.terms[ti];
        let offset = t.positions[0];
        let table = t
            .table
            .iter()
            .map(|&v| if v == FORBIDDEN { forbid } else { v })
            .collect();
        let positions = t.positions.iter().map(|p| p - offset).collect();
        blocks.push(finish(offset, Reader::Sparse { positions }, table, forbid));
    }
    // Most discriminating blocks first: they prune the most.
    let spread = |b: &Block| {
        let min = b
            .table
            .iter()
            .copied()
            .filter(|&v| v >= 0)
            .min()
            .unwrap_or(0);
        b.max - min
    };
    blocks.sort_by_cached_key(|b| std::cmp::Reverse(spread(b)));
    let mut rest_max = vec![0i32; blocks.len() + 1];
    for j in (0..blocks.len()).rev() {
        rest_max[j] = rest_max[j + 1] + blocks[j].max;
    }
    Ok(Lowered {
        blocks,
        rest_max,
        width: model.width,
        lo,
    })
}

fn finish(offset: i32, reader: Reader, table: Vec<i32>, forbid: i32) -> Block {
    let max = table
        .iter()
        .copied()
        .filter(|&v| v != forbid)
        .max()
        .unwrap_or(0);
    Block {
        offset,
        reader,
        table,
        max: max.max(0),
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::score::IntTerm;

    fn pwm(width: i32) -> IntModel {
        let terms = (0..width)
            .map(|p| IntTerm {
                positions: vec![p],
                table: vec![p, 1, 2, 3],
            })
            .collect::<Vec<_>>();
        let max = (0..width).map(|p| p.max(3)).sum();
        IntModel::for_tests(terms, max, 4)
    }

    #[test]
    fn packing_by_alphabet() {
        assert_eq!(Packing::new(4).k, 5);
        assert_eq!(Packing::new(20).k, 2);
        assert_eq!(Packing::new(255).k, 1);
    }

    #[test]
    fn pwm_blocks_cover_every_position_once() {
        let p = Packing::new(4);
        let l = lower(&pwm(12), &p).unwrap();
        let mut offsets: Vec<i32> = l.blocks.iter().map(|b| b.offset).collect();
        offsets.sort();
        assert_eq!(offsets, vec![0, 5, 7]);
        let l = lower(&pwm(3), &p).unwrap();
        assert_eq!(l.blocks.len(), 1);
        assert!(matches!(l.blocks[0].reader, Reader::Kmer { shift: 4 }));
    }

    #[test]
    fn long_pair_becomes_sparse() {
        let terms = vec![
            IntTerm {
                positions: vec![0, 9],
                table: (0..16).collect(),
            },
            IntTerm {
                positions: vec![4],
                table: vec![0, 1, 2, 3],
            },
        ];
        let l = lower(&IntModel::for_tests(terms, 18, 4), &Packing::new(4)).unwrap();
        assert_eq!(l.blocks.len(), 2);
        assert!(
            l.blocks
                .iter()
                .any(|b| matches!(&b.reader, Reader::Sparse { positions } if positions == &[0, 9]))
        );
    }
}
