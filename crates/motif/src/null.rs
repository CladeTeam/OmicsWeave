//! Step ⑤: the exact distribution of integer window scores under the
//! background, and conversions between p-values and scores.

use omicsweave_core::{Error, Result};

use crate::background::Background;
use crate::score::{FORBIDDEN, IntModel};
use crate::scoring::Rounding;

/// The distribution of a window's integer score when the window (and the
/// context it reads) is drawn from the background.
#[derive(Clone, Debug, PartialEq)]
pub struct NullDist {
    pmf: Vec<f64>,
    survival: Vec<f64>,
    forbidden: f64,
}

impl NullDist {
    /// `P(S = s)` for `s = 0 ..= max`.
    pub fn pmf(&self) -> &[f64] {
        &self.pmf
    }

    /// `P(S ≥ s)` for `s = 0 ..= max`.
    pub fn survival(&self) -> &[f64] {
        &self.survival
    }

    /// Probability that the window is forbidden (has no score).
    pub fn forbidden(&self) -> f64 {
        self.forbidden
    }

    /// `P(S ≥ s)`: the scorable probability below 0, and 0 above the largest
    /// score. Forbidden windows have no finite score.
    pub fn pvalue(&self, s: i32) -> f64 {
        if s <= 0 {
            return self.survival.first().copied().unwrap_or(0.0);
        }
        self.survival.get(s as usize).copied().unwrap_or(0.0)
    }

    /// The smallest score whose p-value, after `rounding`, is at most `p`;
    /// `None` when no score is that significant.
    pub fn cutoff(&self, p: f64, rounding: Rounding) -> Option<i32> {
        let ok = |s: usize| rounding.apply(self.survival[s]) <= p;
        let (mut lo, mut hi) = (0usize, self.survival.len());
        while lo < hi {
            let mid = lo + (hi - lo) / 2;
            if ok(mid) {
                hi = mid;
            } else {
                lo = mid + 1;
            }
        }
        (lo < self.survival.len()).then_some(lo as i32)
    }

    /// The largest score with positive probability.
    pub fn max_reachable(&self) -> Option<i32> {
        self.pmf.iter().rposition(|&p| p > 0.0).map(|s| s as i32)
    }

    fn from_pmf(pmf: Vec<f64>, forbidden: f64) -> NullDist {
        let mut survival = pmf.clone();
        for i in (0..survival.len().saturating_sub(1)).rev() {
            let p = pmf[i] + survival[i + 1];
            survival[i] = p.min(1.0);
        }
        NullDist {
            pmf,
            survival,
            forbidden,
        }
    }
}

/// The null distribution of `model` when the background chain starts
/// `context` symbols before the window (generating with the longest context
/// available, as the model's background terms do).
pub(crate) fn null_distribution(
    model: &IntModel,
    background: &Background,
    budget: u64,
) -> Result<NullDist> {
    let pwm = background.order() == 0 && model.terms.iter().all(|t| t.positions.len() == 1);
    if pwm {
        Ok(convolve(model, background.letters()))
    } else {
        banded(model, background, budget)
    }
}

/// The null distribution of single-position terms under the given letter
/// probabilities, used as they are.
pub(crate) fn null_with_letters(model: &IntModel, letters: &[f64]) -> Result<NullDist> {
    if model.terms.iter().any(|t| t.positions.len() != 1) || letters.len() != model.alphabet_len() {
        return Err(Error::invalid(
            "motif null distribution",
            "explicit null letters need single-position terms over the alphabet",
        ));
    }
    Ok(convolve(model, letters))
}

/// Position by position convolution for single-position terms under a
/// 0-order background, adding in MEME's order (letters outer, scores
/// ascending inner) so the table matches FIMO's bit for bit.
fn convolve(model: &IntModel, letters: &[f64]) -> NullDist {
    let size = model.max as usize + 1;
    let mut old = vec![0.0; size];
    let mut new = vec![0.0; size];
    old[0] = 1.0;
    let mut top = 0usize;
    let mut forbidden = 0.0;
    for t in &model.terms {
        let row_max = t
            .table
            .iter()
            .copied()
            .filter(|&s| s != FORBIDDEN)
            .max()
            .unwrap_or(0) as usize;
        new[..=top + row_max].fill(0.0);
        let mut mass = None;
        for (j, &s) in t.table.iter().enumerate() {
            if s == FORBIDDEN {
                let mass = *mass.get_or_insert_with(|| old[..=top].iter().sum::<f64>());
                forbidden += mass * letters[j];
                continue;
            }
            let s = s as usize;
            for (out, &o) in new[s..=s + top].iter_mut().zip(&old[..=top]) {
                if o != 0.0 {
                    *out += o * letters[j];
                }
            }
        }
        top += row_max;
        std::mem::swap(&mut old, &mut new);
    }
    old.truncate(top + 1);
    old.resize(size, 0.0);
    NullDist::from_pmf(old, forbidden)
}

/// Variable elimination along the footprint. The state at position `p` is
/// the symbols at earlier positions that are still read: by a term ending
/// at or after `p`, or as background context. A long-range pair of
/// positions therefore costs one extra symbol, not the whole span.
fn banded(model: &IntModel, background: &Background, budget: u64) -> Result<NullDist> {
    let n = model.alphabet_len();
    let m = background.order() as i32;
    let lo = model
        .terms
        .iter()
        .map(|t| t.positions[0])
        .min()
        .unwrap_or(0);
    let hi = model
        .terms
        .iter()
        .map(|t| t.positions[t.positions.len() - 1])
        .max()
        .unwrap_or(0);
    let stride = model.max as usize + 1;
    // The last position at which each position is still read by a term.
    let mut last_use = vec![i32::MIN; (hi - lo + 1) as usize];
    for t in &model.terms {
        let end = t.positions[t.positions.len() - 1];
        for &q in &t.positions {
            let u = &mut last_use[(q - lo) as usize];
            *u = (*u).max(end);
        }
    }
    let advance = |before: &mut Vec<i32>, p: i32| {
        let next = p + 1;
        before.retain(|&q| p < hi && (last_use[(q - lo) as usize] >= next || q >= next - m));
        if p < hi && (last_use[(p - lo) as usize] >= next || m > 0) {
            before.push(p);
        }
    };
    let mut cost = 0u128;
    let mut max_states = 1u128;
    let mut before = Vec::new();
    for p in lo..=hi {
        let states = (n as u128)
            .checked_pow(before.len() as u32)
            .unwrap_or(u128::MAX);
        cost = cost.saturating_add(states.saturating_mul(n as u128 * stride as u128));
        max_states = max_states.max(states);
        advance(&mut before, p);
    }
    if cost > budget as u128 || max_states.saturating_mul(stride as u128) > 1 << 27 {
        return Err(too_large(cost.min(u64::MAX as u128) as u64, budget));
    }
    let max_states = max_states as usize;
    let mut cur = vec![0.0; max_states * stride];
    let mut next = vec![0.0; max_states * stride];
    cur[0] = 1.0;
    let mut top = 0usize;
    let mut forbidden = 0.0;
    let mut digits = Vec::new();
    let mut after = Vec::new();
    let mut add = Vec::new();
    let mut prob = Vec::new();
    let mut to = Vec::new();
    let mut term_order: Vec<usize> = (0..model.terms.len()).collect();
    term_order.sort_unstable_by_key(|&i| *model.terms[i].positions.last().unwrap());
    let mut term_at = 0usize;
    for (h, p) in (lo..=hi).enumerate() {
        after.clear();
        after.extend_from_slice(&before);
        advance(&mut after, p);
        let mut slots = before.clone();
        slots.push(p);
        let slot = |q: i32| slots.iter().position(|&s| s == q).expect("live position");
        let k = slots.len();
        let combos = n.pow(k as u32);
        let ctx = (m as usize).min(h);
        let cond = background.conditional(ctx);
        let ctx_slots: Vec<usize> = (p - ctx as i32..p).map(slot).collect();
        let after_slots: Vec<usize> = after.iter().map(|&q| slot(q)).collect();
        let start = term_at;
        while term_at < term_order.len()
            && *model.terms[term_order[term_at]].positions.last().unwrap() == p
        {
            term_at += 1;
        }
        let ending: Vec<(&[i32], Vec<usize>)> = term_order[start..term_at]
            .iter()
            .map(|&i| &model.terms[i])
            .map(|t| {
                (
                    t.table.as_slice(),
                    t.positions.iter().map(|&q| slot(q)).collect(),
                )
            })
            .collect();
        // Per (state, symbol): score added, probability, next state.
        add.clear();
        add.resize(combos, 0i64);
        prob.clear();
        prob.resize(combos, 0.0);
        to.clear();
        to.resize(combos, 0usize);
        digits.resize(k, 0);
        for c in 0..combos {
            crate::background::decode(c, n, &mut digits);
            for (table, at) in &ending {
                if add[c] == i64::MIN {
                    break;
                }
                let idx = at.iter().fold(0usize, |a, &s| a * n + digits[s]);
                let v = table[idx];
                if v == FORBIDDEN {
                    add[c] = i64::MIN;
                } else {
                    add[c] += v as i64;
                }
            }
            let ci = ctx_slots.iter().fold(0usize, |a, &s| a * n + digits[s]);
            prob[c] = cond[ci * n + digits[k - 1]];
            to[c] = after_slots.iter().fold(0usize, |a, &s| a * n + digits[s]);
        }
        let add_max = add.iter().copied().max().unwrap_or(0).max(0) as usize;
        let states_in = combos / n;
        let states_out = n.pow(after.len() as u32);
        next[..states_out * stride]
            .iter_mut()
            .for_each(|v| *v = 0.0);
        for st in 0..states_in {
            let row = &cur[st * stride..st * stride + top + 1];
            let mass: f64 = row.iter().sum();
            if mass == 0.0 {
                continue;
            }
            for x in 0..n {
                let c = st * n + x;
                let pr = prob[c];
                if add[c] == i64::MIN {
                    forbidden += mass * pr;
                    continue;
                }
                let a = add[c] as usize;
                let base = to[c] * stride + a;
                for (o, &v) in next[base..base + top + 1].iter_mut().zip(row) {
                    if v != 0.0 {
                        *o += v * pr;
                    }
                }
            }
        }
        top += add_max;
        std::mem::swap(&mut cur, &mut next);
        std::mem::swap(&mut before, &mut after);
    }
    // No position is live after the last one: one state remains.
    let pmf = cur[..stride].to_vec();
    Ok(NullDist::from_pmf(pmf, forbidden))
}

fn too_large(cost: u64, budget: u64) -> Error {
    Error::invalid(
        "motif null distribution",
        format!("exact distribution needs ~{cost} steps, budget is {budget}"),
    )
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::score::IntTerm;

    fn model(terms: Vec<(Vec<i32>, Vec<i32>)>, n: usize) -> IntModel {
        let max = terms
            .iter()
            .map(|(_, t)| *t.iter().filter(|&&v| v != FORBIDDEN).max().unwrap())
            .sum();
        IntModel::for_tests(
            terms
                .into_iter()
                .map(|(positions, table)| IntTerm { positions, table })
                .collect(),
            max,
            n,
        )
    }

    #[test]
    fn convolution_matches_banded() {
        let m = model(
            vec![(vec![0], vec![0, 1, 2, 3]), (vec![1], vec![3, 0, 0, 1])],
            4,
        );
        let bg = Background::frequencies(&[0.1, 0.2, 0.3, 0.4]).unwrap();
        let a = convolve(&m, bg.letters());
        let b = banded(&m, &bg, u64::MAX).unwrap();
        for (x, y) in a.pmf().iter().zip(b.pmf()) {
            assert!((x - y).abs() < 1e-15);
        }
        assert!((a.pmf().iter().sum::<f64>() - 1.0).abs() < 1e-12);
        assert_eq!(a.survival()[0], 1.0);
    }

    #[test]
    fn pair_term_matches_enumeration() {
        // One term over positions 0 and 2 of a binary alphabet.
        let m = model(vec![(vec![0, 2], vec![0, 1, 2, FORBIDDEN])], 2);
        let bg = Background::frequencies(&[0.25, 0.75]).unwrap();
        let d = banded(&m, &bg, u64::MAX).unwrap();
        assert!((d.pmf()[0] - 0.0625).abs() < 1e-15);
        assert!((d.pmf()[1] - 0.1875).abs() < 1e-15);
        assert!((d.pmf()[2] - 0.1875).abs() < 1e-15);
        assert!((d.forbidden() - 0.5625).abs() < 1e-15);
        let distant = model(vec![(vec![0, 999], vec![0, 1, 2, FORBIDDEN])], 2);
        assert_eq!(banded(&distant, &bg, 1 << 20).unwrap(), d);
    }

    #[test]
    fn markov_background_restarts_at_the_footprint() {
        // P(x0 x1) under a first-order chain started at x0.
        let bg =
            Background::from_kmer_frequencies(2, vec![vec![1.0, 1.0], vec![3.0, 1.0, 1.0, 3.0]])
                .unwrap();
        let m = model(vec![(vec![0, 1], vec![0, 1, 1, 0])], 2);
        let d = banded(&m, &bg, u64::MAX).unwrap();
        // Same-symbol pairs score 0 with probability 0.5 · 0.75 · 2.
        assert!((d.pmf()[0] - 0.75).abs() < 1e-15);
        assert!((d.pmf()[1] - 0.25).abs() < 1e-15);
    }

    #[test]
    fn cutoff_respects_rounding() {
        let d = NullDist::from_pmf(vec![0.5, 0.25, 0.125, 0.125], 0.0);
        assert_eq!(d.survival(), &[1.0, 0.5, 0.25, 0.125]);
        assert_eq!(d.cutoff(0.25, Rounding::None), Some(2));
        assert_eq!(d.cutoff(0.1, Rounding::None), None);
        assert_eq!(d.pvalue(9), 0.0);
    }
}
