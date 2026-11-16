//! Step ③ (contrast) and ④ (quantize): smoothed factors and the background
//! become additive log-odds terms, then integer terms.

use std::collections::HashMap;

use omicsweave_core::{Error, Result};
use omicsweave_ffm::FactorKind;

use crate::background::{Background, decode};
use crate::factor::Smoothed;
use crate::scoring::{Ambiguity, AmbiguityRule, LogScale, Quantizer, meme_round};

/// One additive term: a table over the global codes at a few window
/// positions (negative positions read the context before the window).
#[derive(Clone, Debug, PartialEq)]
pub struct Term {
    /// Strictly increasing positions.
    pub positions: Vec<i32>,
    /// `n^k` log-odds cells, last position fastest; `-∞` forbids.
    pub table: Vec<f64>,
}

/// A window score as a sum of terms: `S(x) = Σ_t T_t[x_{P_t}]`.
#[derive(Clone, Debug, PartialEq)]
pub struct ScoreModel {
    pub width: u32,
    pub terms: Vec<Term>,
    /// Per term, the cells of the ambiguity codes (see
    /// [`crate::Ambiguity`]); empty when the codes are not scored, or for a
    /// term over several positions.
    pub ambiguous: Vec<Vec<f64>>,
}

impl ScoreModel {
    /// The leftmost position read (≤ 0 when background terms reach into the
    /// context).
    pub fn lo(&self) -> i32 {
        self.terms
            .iter()
            .filter_map(|t| t.positions.first().copied())
            .min()
            .unwrap_or(0)
            .min(0)
    }
}

/// Builds the score model of one strand. `context` is how many symbols
/// before the window the background may condition on (the full order for
/// ordinary windows, fewer near a sequence start or an unknown symbol).
pub(crate) fn contrast(
    s: &Smoothed,
    background: &Background,
    context: usize,
    log: LogScale,
    ambiguity: &Ambiguity,
    weights: &[f64],
) -> ScoreModel {
    let n = background.alphabet_len();
    let m = background.order();
    let letters = background.letters();
    let mut used = vec![false; s.factors.len()];
    let mut single = vec![None; s.width as usize];
    for (i, f) in s.factors.iter().enumerate() {
        if f.kind != FactorKind::Potential && f.positions.len() == 1 {
            single[f.positions[0] as usize].get_or_insert(i);
        }
    }
    let mut terms: Vec<Term> = Vec::with_capacity(s.factors.len() + s.width as usize);
    for (p, &single_factor) in single.iter().enumerate() {
        if !s.covered[p] {
            continue;
        }
        let pi = p as i32;
        let len = m.min(p + context);
        if len == 0 {
            // A single-position factor at p absorbs the background as one
            // ratio, `log(q / b)`, the way MEME builds a PSSM.
            if let Some(i) = single_factor {
                used[i] = true;
                let table = s.factors[i]
                    .table
                    .iter()
                    .zip(letters)
                    .map(|(&q, &b)| log.log(q / b))
                    .collect();
                terms.push(Term {
                    positions: vec![pi],
                    table,
                });
                continue;
            }
        }
        let cond = background.conditional(len);
        terms.push(Term {
            positions: (pi - len as i32..=pi).collect(),
            table: cond.iter().map(|&b| -log.log(b)).collect(),
        });
    }
    for (f, _) in s.factors.iter().zip(&used).filter(|(_, u)| !**u) {
        terms.push(Term {
            positions: f.positions.clone(),
            table: f.table.iter().map(|&v| log.log(v)).collect(),
        });
    }
    let mut model = canonicalize(terms, n, s.width);
    if ambiguity.rule != AmbiguityRule::Skip {
        model.ambiguous = model
            .terms
            .iter()
            .map(|t| {
                if t.positions.len() != 1 {
                    return Vec::new();
                }
                ambiguity
                    .codes
                    .iter()
                    .map(|(_, set)| match ambiguity.rule {
                        AmbiguityRule::Sum => set.iter().map(|&c| t.table[c as usize]).sum(),
                        _ => {
                            let mut num = 0.0;
                            let mut den = 0.0;
                            for &c in set {
                                num += t.table[c as usize] * weights[c as usize];
                                den += weights[c as usize];
                            }
                            if den == 0.0 { 0.0 } else { num / den }
                        }
                    })
                    .collect()
            })
            .collect();
    }
    model
}

/// Merges terms over the same positions, folds a term into one whose
/// positions contain it, and orders terms by their last position.
fn canonicalize(terms: Vec<Term>, n: usize, width: u32) -> ScoreModel {
    let mut merged: Vec<Term> = Vec::with_capacity(terms.len());
    let mut by_positions = (terms.len() >= 64).then(|| HashMap::with_capacity(terms.len()));
    for t in terms {
        let same = match &by_positions {
            Some(index) => index.get(&t.positions).copied(),
            None => merged.iter().position(|u| u.positions == t.positions),
        };
        if let Some(i) = same {
            let u: &mut Term = &mut merged[i];
            for (a, b) in u.table.iter_mut().zip(&t.table) {
                *a += b;
            }
        } else {
            if let Some(index) = &mut by_positions {
                index.insert(t.positions.clone(), merged.len());
            }
            merged.push(t);
        }
    }
    // Fold smaller terms into a containing one, smallest first.
    let mut order: Vec<usize> = (0..merged.len()).collect();
    order.sort_by_key(|&i| merged[i].positions.len());
    let mut gone = vec![false; merged.len()];
    for &i in &order {
        let longer =
            order.partition_point(|&j| merged[j].positions.len() <= merged[i].positions.len());
        let host = order[longer..].iter().copied().find(|&j| {
            !gone[j]
                && merged[j].positions.len() > merged[i].positions.len()
                && merged[i]
                    .positions
                    .iter()
                    .all(|p| merged[j].positions.contains(p))
        });
        if let Some(j) = host {
            let small = merged[i].clone();
            broadcast_add(&mut merged[j], &small, n);
            gone[i] = true;
        }
    }
    let mut terms: Vec<Term> = merged
        .into_iter()
        .zip(gone)
        .filter(|(_, g)| !g)
        .map(|(t, _)| t)
        .collect();
    terms.sort_by_key(|t| (t.positions[t.positions.len() - 1], t.positions[0]));
    ScoreModel {
        width,
        terms,
        ambiguous: Vec::new(),
    }
}

/// Adds `small`'s table to `host`'s, `small.positions ⊂ host.positions`.
fn broadcast_add(host: &mut Term, small: &Term, n: usize) {
    let axes: Vec<usize> = small
        .positions
        .iter()
        .map(|p| host.positions.iter().position(|q| q == p).expect("subset"))
        .collect();
    let mut digits = vec![0usize; host.positions.len()];
    for (i, v) in host.table.iter_mut().enumerate() {
        decode(i, n, &mut digits);
        let j = axes.iter().fold(0usize, |acc, &a| acc * n + digits[a]);
        *v += small.table[j];
    }
}

/// The cell value of a forbidden combination.
pub const FORBIDDEN: i32 = i32::MIN;

/// A term with integer cells; [`FORBIDDEN`] marks a forbidden combination.
#[derive(Clone, Debug, PartialEq)]
pub struct IntTerm {
    pub positions: Vec<i32>,
    pub table: Vec<i32>,
}

/// A score model with integer cells: every allowed cell is ≥ 0, and an
/// integer score `S` stands for the real score `S / scale + offset`.
#[derive(Clone, Debug, PartialEq)]
pub struct IntModel {
    pub width: u32,
    pub terms: Vec<IntTerm>,
    /// Per term, the integer cells of the ambiguity codes (empty when not
    /// scored).
    pub ambiguous: Vec<Vec<i32>>,
    pub scale: f64,
    pub offset: f64,
    /// The largest integer score (the sum of each term's largest cell).
    pub max: i32,
    alphabet_len: usize,
}

impl IntModel {
    /// The real score of integer score `s`.
    pub fn raw(&self, s: i32) -> f64 {
        s as f64 / self.scale + self.offset
    }

    /// The leftmost position read.
    pub fn lo(&self) -> i32 {
        self.terms
            .iter()
            .filter_map(|t| t.positions.first().copied())
            .min()
            .unwrap_or(0)
            .min(0)
    }

    pub fn alphabet_len(&self) -> usize {
        self.alphabet_len
    }

    #[cfg(test)]
    pub(crate) fn for_tests(terms: Vec<IntTerm>, max: i32, alphabet_len: usize) -> IntModel {
        IntModel {
            width: terms
                .iter()
                .map(|t| *t.positions.last().unwrap() as u32 + 1)
                .max()
                .unwrap_or(0),
            terms,
            ambiguous: Vec::new(),
            scale: 1.0,
            offset: 0.0,
            max,
            alphabet_len,
        }
    }

    /// The integer score of a window whose symbols are alphabet codes or
    /// ambiguity codes (`ambiguity[i]` is the code index at `i`, `0xFF` for
    /// none). `None` when forbidden, or when an ambiguity code falls in a
    /// term that does not score it.
    pub fn score_ambiguous(&self, codes: &[u8], ambiguity: &[u8], at: usize) -> Option<i32> {
        let n = self.alphabet_len;
        let mut s = 0i64;
        for (ti, t) in self.terms.iter().enumerate() {
            let first = at.checked_add_signed(*t.positions.first()? as isize)?;
            let first_ambiguity = *ambiguity.get(first)?;
            let v = if t.positions.len() == 1 && first_ambiguity != 0xFF {
                *self.ambiguous.get(ti)?.get(first_ambiguity as usize)?
            } else {
                let mut idx = 0usize;
                for &p in &t.positions {
                    let i = at.checked_add_signed(p as isize)?;
                    if *ambiguity.get(i)? != 0xFF {
                        return None;
                    }
                    let c = *codes.get(i)? as usize;
                    if c >= n {
                        return None;
                    }
                    idx = idx.checked_mul(n)?.checked_add(c)?;
                }
                *t.table.get(idx)?
            };
            if v == FORBIDDEN {
                return None;
            }
            s = s.checked_add(v as i64)?;
        }
        s.try_into().ok()
    }

    /// The integer score of the window at `at` over global `codes` (no
    /// unknown symbols in the positions read), `None` when forbidden.
    pub fn score(&self, codes: &[u8], at: usize) -> Option<i32> {
        let n = self.alphabet_len;
        let mut s = 0i64;
        for t in &self.terms {
            if t.positions.is_empty() {
                return None;
            }
            let mut idx = 0usize;
            for &p in &t.positions {
                let i = at.checked_add_signed(p as isize)?;
                let c = *codes.get(i)? as usize;
                if c >= n {
                    return None;
                }
                idx = idx.checked_mul(n)?.checked_add(c)?;
            }
            let v = *t.table.get(idx)?;
            if v == FORBIDDEN {
                return None;
            }
            s = s.checked_add(v as i64)?;
        }
        s.try_into().ok()
    }
}

/// Integer cells from real ones. `scale` reuses another model's scale so
/// variants of one strand stay comparable. Ambiguity cells take part in the
/// range like any other cell.
pub(crate) fn quantize(
    model: &ScoreModel,
    quantizer: Quantizer,
    alphabet_len: usize,
    scale: Option<f64>,
) -> Result<IntModel> {
    let nint = |x: f64| -> i64 {
        if x >= 0.0 {
            (x + 0.5) as i64
        } else {
            (x - 0.5) as i64
        }
    };
    let cells = |i: usize| {
        let extra = model.ambiguous.get(i).map_or(&[][..], Vec::as_slice);
        finite(&model.terms[i]).chain(extra.iter().copied().filter(|v| v.is_finite()))
    };
    let span = || {
        let mut small = f64::INFINITY;
        let mut large = f64::NEG_INFINITY;
        for i in 0..model.terms.len() {
            for v in cells(i) {
                small = small.min(v);
                large = large.max(v);
            }
        }
        (small, large)
    };
    // The scale and each term's shift: cell = nint((x − shift) · scale).
    let (scale, shifts): (f64, Vec<f64>) = match quantizer {
        Quantizer::Fimo { range } => {
            let (mut small, large) = span();
            if !small.is_finite() {
                return Err(invalid("every window is forbidden"));
            }
            if large == small {
                small = large - 1.0;
            }
            let offset = small.floor();
            let scale = scale.unwrap_or_else(|| (range as f64 / (large - offset)).floor());
            if scale.is_nan() || scale < 1.0 {
                return Err(invalid(format!(
                    "scores span {} bits, more than the range {range} can scale",
                    large - offset
                )));
            }
            (scale, vec![offset; model.terms.len()])
        }
        Quantizer::Mast { range } => {
            let (small, large) = span();
            if !small.is_finite() {
                return Err(invalid("every window is forbidden"));
            }
            if large == small {
                return Err(invalid("the motif cannot be scaled: every score is equal"));
            }
            let scale = scale.unwrap_or_else(|| meme_round(range as f64 / (large - small), 14));
            (scale, vec![small; model.terms.len()])
        }
        Quantizer::Precision { states } => {
            let mut bounds = Vec::with_capacity(model.terms.len());
            for i in 0..model.terms.len() {
                let lo = cells(i).fold(f64::INFINITY, f64::min);
                let hi = cells(i).fold(f64::NEG_INFINITY, f64::max);
                if !lo.is_finite() {
                    return Err(invalid("every window is forbidden"));
                }
                bounds.push((lo, hi));
            }
            let spread: f64 = bounds.iter().map(|(lo, hi)| hi - lo).sum();
            let scale = scale.unwrap_or(if spread > 0.0 {
                (states - 1) as f64 / spread
            } else {
                1.0
            });
            (scale, bounds.iter().map(|(lo, _)| *lo).collect())
        }
    };
    if !(scale.is_finite() && scale > 0.0) {
        return Err(invalid("quantization scale is not finite and > 0"));
    }
    let convert = |v: f64, shift: f64| -> Result<i32> {
        if v.is_finite() {
            let scaled = (v - shift) * scale;
            if !scaled.is_finite() {
                return Err(invalid("quantized cell overflows"));
            }
            i32::try_from(nint(scaled)).map_err(|_| invalid("quantized cell overflows i32"))
        } else {
            Ok(FORBIDDEN)
        }
    };
    let mut terms = Vec::with_capacity(model.terms.len());
    let mut ambiguous = Vec::with_capacity(model.ambiguous.len());
    for (i, (t, &shift)) in model.terms.iter().zip(&shifts).enumerate() {
        terms.push(IntTerm {
            positions: t.positions.clone(),
            table: t
                .table
                .iter()
                .map(|&v| convert(v, shift))
                .collect::<Result<_>>()?,
        });
        if let Some(extra) = model.ambiguous.get(i) {
            ambiguous.push(
                extra
                    .iter()
                    .map(|&v| convert(v, shift))
                    .collect::<Result<Vec<i32>>>()?,
            );
        }
    }
    let offset = match quantizer {
        Quantizer::Precision { .. } => shifts.iter().sum(),
        _ => model.terms.len() as f64 * shifts.first().copied().unwrap_or(0.0),
    };
    let mut max = 0i64;
    for (i, t) in terms.iter().enumerate() {
        let allowed = || {
            t.table
                .iter()
                .chain(ambiguous.get(i).map_or(&[][..], Vec::as_slice))
                .copied()
                .filter(|&v| v != FORBIDDEN)
        };
        if allowed().min().is_some_and(|v| v < 0) {
            return Err(invalid("a quantized cell is negative"));
        }
        max += t
            .table
            .iter()
            .copied()
            .filter(|&v| v != FORBIDDEN)
            .max()
            .unwrap_or(0) as i64;
    }
    if max > 1 << 30 {
        return Err(invalid(format!("integer scores reach {max}, above 2^30")));
    }
    Ok(IntModel {
        width: model.width,
        terms,
        ambiguous,
        scale,
        offset,
        max: max as i32,
        alphabet_len,
    })
}

fn finite(t: &Term) -> impl Iterator<Item = f64> + '_ {
    t.table.iter().copied().filter(|v| v.is_finite())
}

fn invalid(detail: impl Into<String>) -> Error {
    Error::invalid("motif model", detail)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn public_scoring_rejects_missing_positions_codes_tables_and_overflow() {
        let mut model = IntModel::for_tests(
            vec![IntTerm {
                positions: vec![-1, 0],
                table: vec![0, 1, 2, 3],
            }],
            3,
            2,
        );
        assert_eq!(model.score(&[1, 0], 1), Some(2));
        assert_eq!(model.score(&[1, 0], 0), None);
        assert_eq!(model.score(&[1, 2], 1), None);
        assert_eq!(model.score(&[], usize::MAX), None);
        assert_eq!(model.score_ambiguous(&[1, 0], &[0xFF], 1), None);
        model.terms[0].table.clear();
        assert_eq!(model.score(&[1, 0], 1), None);
        model.terms[0].positions.clear();
        assert_eq!(model.lo(), 0);
        assert_eq!(model.score(&[1, 0], 1), None);
        assert_eq!(model.score_ambiguous(&[1, 0], &[0xFF; 2], 1), None);
        model.terms = vec![
            IntTerm {
                positions: vec![0],
                table: vec![i32::MAX; 2]
            };
            2
        ];
        assert_eq!(model.score(&[0], 0), None);
        assert_eq!(model.score_ambiguous(&[0], &[0xFF], 0), None);
    }

    #[test]
    fn canonicalized_terms_preserve_the_additive_model() {
        let terms = vec![
            Term {
                positions: vec![0],
                table: vec![1.0, 2.0],
            },
            Term {
                positions: vec![1],
                table: vec![3.0, 5.0],
            },
            Term {
                positions: vec![0],
                table: vec![7.0, 11.0],
            },
            Term {
                positions: vec![0, 3],
                table: vec![13.0, 17.0, 19.0, 23.0],
            },
        ];
        let compact = canonicalize(terms.clone(), 2, 4);
        for word in 0..16 {
            let codes: Vec<_> = (0..4).map(|i| (word >> i) & 1).collect();
            let score = |terms: &[Term]| {
                terms
                    .iter()
                    .map(|t| {
                        let index = t
                            .positions
                            .iter()
                            .fold(0usize, |a, &p| a * 2 + codes[p as usize]);
                        t.table[index]
                    })
                    .sum::<f64>()
            };
            assert_eq!(score(&terms), score(&compact.terms));
        }
    }

    #[test]
    fn quantization_rejects_unrepresentable_cells_and_scale() {
        let real = ScoreModel {
            width: 1,
            terms: vec![Term {
                positions: vec![0],
                table: vec![0.0, 4_294_967_297.0],
            }],
            ambiguous: Vec::new(),
        };
        assert!(quantize(&real, Quantizer::Precision { states: 128 }, 2, Some(1.0)).is_err());
        assert!(
            quantize(
                &real,
                Quantizer::Precision { states: 128 },
                2,
                Some(f64::INFINITY)
            )
            .is_err()
        );
    }
}
