//! BaMM (BaMMmotif2) inhomogeneous Markov models, `.ihbcp`: one model per
//! file, no header. Each position is a block of `K + 1` lines ended by a
//! blank line; line `k` holds the `A^(k+1)` conditional probabilities of
//! order `k`, the current symbol varying fastest and the oldest context
//! symbol slowest (as FFM's factor tables). Near the motif start the writer
//! tiles lower orders into full rows; FFM keeps the true order instead: at
//! position `j` the factor of order `min(j, K)`. The alphabet, told by the
//! width of the first line, is `ACGT` or BaMM's `ACGTMH`; the five-letter
//! `ACGTM` and `ACGTH` cannot be told apart and are refused. The ID is the
//! file name.

use std::path::Path;

use omicsweave_core::Result;
use omicsweave_ffm::{Alphabet, Factor, FactorKind, Motif, Support};

use super::{Entry, Parsed, bad, numbers, stem};

/// Blocks of numeric lines separated by blank lines.
fn blocks(text: &str) -> Option<Vec<Vec<Vec<f64>>>> {
    let mut out: Vec<Vec<Vec<f64>>> = Vec::new();
    let mut cur: Vec<Vec<f64>> = Vec::new();
    for line in text.lines() {
        let line = line.trim();
        if line.is_empty() {
            if !cur.is_empty() {
                out.push(std::mem::take(&mut cur));
            }
            continue;
        }
        cur.push(numbers(line).filter(|r| !r.is_empty())?);
    }
    if !cur.is_empty() {
        out.push(cur);
    }
    Some(out)
}

/// The order `K` and alphabet size of a block whose line `k` has
/// `A^(k+1)` values.
fn shape(block: &[Vec<f64>]) -> Option<(usize, usize)> {
    let a = block.first()?.len();
    if !(4..=6).contains(&a) {
        return None;
    }
    let mut want = a;
    for line in block {
        if line.len() != want {
            return None;
        }
        want = want.checked_mul(a)?;
    }
    Some((block.len() - 1, a))
}

/// Whether text is a BaMM model: every block has the same `A^(k+1)` shape.
pub(crate) fn is_bamm(text: &str, complete: bool) -> bool {
    let Some(mut bs) = blocks(text) else {
        return false;
    };
    // The last block of a sniffed prefix may be cut.
    if !complete && bs.len() > 1 && !text.ends_with("\n\n") {
        bs.pop();
    }
    let Some(first) = bs.first().and_then(|b| shape(b)) else {
        return false;
    };
    bs.iter().all(|b| shape(b) == Some(first))
}

fn alphabet(path: &Path, a: usize) -> Result<Alphabet> {
    match a {
        4 => Ok(Alphabet::dna()),
        6 => Ok(Alphabet::new("BaMM ACGTMH", b"ACGTMH").expect("valid")),
        _ => Err(bad(
            path,
            0,
            "a five-letter BaMM model is ACGTM or ACGTH, which the file does not say",
        )),
    }
}

pub(crate) fn parse(path: &Path, text: &str) -> Result<Parsed> {
    let bs = blocks(text).ok_or_else(|| bad(path, 0, "lines are not numbers"))?;
    let first = bs.first().and_then(|b| shape(b)).ok_or_else(|| {
        bad(
            path,
            0,
            "not a BaMM model: line k of a position needs A^(k+1) values",
        )
    })?;
    let (order, a) = first;
    let alphabet = alphabet(path, a)?;
    let width = bs.len() as u32;
    let mut factors = Vec::with_capacity(bs.len());
    for (j, mut block) in bs.into_iter().enumerate() {
        if shape(&block) != Some(first) {
            return Err(bad(
                path,
                0,
                format!("position {} has a different shape than position 1", j + 1),
            ));
        }
        let k = j.min(order);
        factors.push(Factor {
            kind: FactorKind::Conditional,
            positions: ((j - k) as u32..=j as u32).collect(),
            table: block.swap_remove(k),
        });
    }
    Ok(Parsed {
        alphabet,
        entries: vec![Entry {
            id: stem(path),
            name: String::new(),
            motif: Motif {
                width,
                support: Support::All,
                sites: None,
                factors,
                gaps: None,
            },
        }],
        skipped: Vec::new(),
    })
}
