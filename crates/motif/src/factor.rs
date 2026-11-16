//! Step ① (estimate) and ② (orient): stored frequencies become smoothed
//! probabilities over global codes, oriented to a strand.

use omicsweave_core::{Error, Result};
use omicsweave_ffm::{FactorKind, MotifView};

use crate::background::decode;
use crate::scoring::Scoring;

/// A factor over global codes: smoothed probabilities for joint and
/// conditional factors, the stored weights for potentials.
#[derive(Clone, Debug)]
pub(crate) struct ProbFactor {
    pub kind: FactorKind,
    /// Strictly increasing positions in the motif window.
    pub positions: Vec<i32>,
    /// `n^k` cells, last position fastest.
    pub table: Vec<f64>,
}

/// A motif's factors on one strand.
#[derive(Clone, Debug)]
pub(crate) struct Smoothed {
    pub width: u32,
    pub factors: Vec<ProbFactor>,
    /// Positions the model generates (a joint factor's positions, a
    /// conditional factor's last one); each is contrasted with the
    /// background once.
    pub covered: Vec<bool>,
}

/// Normalizes, smooths and expands every factor of `view` to global codes.
/// `letters` is the background's letter distribution: the pseudocount's
/// prior mean and the distribution inside the rest category.
pub(crate) fn estimate(
    view: &MotifView<'_>,
    alphabet_len: usize,
    letters: &[f64],
    scoring: &Scoring,
) -> Result<Smoothed> {
    let n = alphabet_len;
    if view.alphabet_len() != n {
        return Err(invalid(format!(
            "motif is over {} symbols, alphabet has {n}",
            view.alphabet_len()
        )));
    }
    if view.gaps().is_some() {
        return Err(invalid(
            "a gapped profile is scored by alignment, not as a factor model",
        ));
    }
    let width = view.width();
    // Local code of each global code: its rank in the support, or the rest.
    let support = view.support();
    let r = view.arity();
    let mut local = vec![0usize; n];
    let mut local_bg = vec![0.0; r];
    match support {
        None => {
            for (g, l) in local.iter_mut().enumerate() {
                *l = g;
            }
            local_bg.copy_from_slice(letters);
        }
        Some(codes) => {
            let rest = codes.len();
            local.iter_mut().for_each(|l| *l = rest);
            for (i, &g) in codes.iter().enumerate() {
                local[g as usize] = i;
            }
            for (g, &l) in local.iter().enumerate() {
                local_bg[l] += letters[g];
            }
        }
    }
    let rest = support.map(<[u32]>::len);
    let sites = view.sites().unwrap_or(scoring.default_sites);
    let sites = if scoring.truncate_sites {
        sites.trunc()
    } else {
        sites
    };
    let alpha = scoring.pseudocount;
    let total = sites + alpha;
    let scaled_weights = total.is_infinite().then(|| {
        let scale = sites.max(alpha);
        (sites / scale, alpha / scale)
    });
    let smooth = |value: f64, prior: f64| match scaled_weights {
        Some((observed, pseudo)) => (value * observed + prior * pseudo) / (observed + pseudo),
        None => (value * sites + alpha * prior) / total,
    };

    let mut factors = Vec::with_capacity(view.factor_count());
    let mut covered = vec![false; width as usize];
    for f in view.factors() {
        let k = f.positions.len();
        let mut p = f.table.to_vec();
        match f.kind {
            FactorKind::Joint => {
                let size = p.len();
                normalize_groups(&mut p, size, scoring.normalize_tolerance)?
            }
            FactorKind::Conditional => normalize_groups(&mut p, r, scoring.normalize_tolerance)?,
            FactorKind::Potential => {}
        }
        if let (Some(prior), FactorKind::Joint, 1) = (&scoring.column_prior, f.kind, k) {
            p = prior.smooth(&p, sites, &local_bg);
            if p.len() != r
                || p.iter().any(|v| !v.is_finite() || *v < 0.0)
                || !p.iter().any(|&v| v > 0.0)
            {
                return Err(invalid("column prior returned an invalid frequency table"));
            }
        } else if alpha > 0.0 {
            match f.kind {
                FactorKind::Joint if k == 1 => {
                    for (v, &b) in p.iter_mut().zip(&local_bg) {
                        *v = smooth(*v, b);
                    }
                }
                FactorKind::Joint => {
                    // Shrink toward the product of the factor's own
                    // (smoothed) marginals: with little data the factor
                    // falls back to independent positions.
                    let mut marginals = vec![vec![0.0; r]; k];
                    let mut digits = vec![0usize; k];
                    for (i, &v) in p.iter().enumerate() {
                        decode(i, r, &mut digits);
                        for (m, &d) in marginals.iter_mut().zip(&digits) {
                            m[d] += v;
                        }
                    }
                    for m in &mut marginals {
                        for (v, &b) in m.iter_mut().zip(&local_bg) {
                            *v = smooth(*v, b);
                        }
                    }
                    for (i, v) in p.iter_mut().enumerate() {
                        decode(i, r, &mut digits);
                        let prior: f64 =
                            marginals.iter().zip(&digits).map(|(m, &d)| m[d]).product();
                        *v = smooth(*v, prior);
                    }
                }
                FactorKind::Conditional => {
                    for group in p.chunks_mut(r) {
                        for (v, &b) in group.iter_mut().zip(&local_bg) {
                            *v = smooth(*v, b);
                        }
                    }
                }
                FactorKind::Potential => {}
            }
        }
        let table = match rest {
            None => p,
            Some(rest) => expand(&p, f.kind, k, n, r, rest, &local, letters, &local_bg)?,
        };
        let positions: Vec<i32> = f.positions.iter().map(|&q| q as i32).collect();
        match f.kind {
            FactorKind::Joint => positions.iter().for_each(|&q| covered[q as usize] = true),
            FactorKind::Conditional => covered[positions[k - 1] as usize] = true,
            FactorKind::Potential => {}
        }
        factors.push(ProbFactor {
            kind: f.kind,
            positions,
            table,
        });
    }
    if factors.is_empty() {
        return Err(invalid("motif has no factors"));
    }
    Ok(Smoothed {
        width,
        factors,
        covered,
    })
}

/// Expands a table over local codes to global codes. A symbol of the rest
/// category takes the rest's probability times its share of the rest's
/// background mass, at every generated position.
#[allow(clippy::too_many_arguments)]
fn expand(
    p: &[f64],
    kind: FactorKind,
    k: usize,
    n: usize,
    r: usize,
    rest: usize,
    local: &[usize],
    letters: &[f64],
    local_bg: &[f64],
) -> Result<Vec<f64>> {
    let size = n
        .checked_pow(k as u32)
        .filter(|&size| size <= 1 << 24)
        .ok_or_else(|| invalid("expanded factor table exceeds 2^24 cells"))?;
    let mut out = vec![0.0; size];
    let mut digits = vec![0usize; k];
    for (g, o) in out.iter_mut().enumerate() {
        decode(g, n, &mut digits);
        let mut li = 0usize;
        let mut v = 1.0;
        for (j, &d) in digits.iter().enumerate() {
            let l = local[d];
            li = li * r + l;
            let generated = match kind {
                FactorKind::Joint => true,
                FactorKind::Conditional => j == k - 1,
                FactorKind::Potential => false,
            };
            if l == rest && generated {
                v *= letters[d] / local_bg[rest];
            }
        }
        *o = p[li] * v;
    }
    Ok(out)
}

fn normalize_groups(values: &mut [f64], group: usize, tolerance: f64) -> Result<()> {
    for g in values.chunks_mut(group) {
        let mut total = 0.0;
        for &value in g.iter() {
            total += value;
        }
        if total == f64::INFINITY {
            let max = g.iter().copied().fold(0.0, f64::max);
            let scaled: f64 = g.iter().map(|v| v / max).sum();
            for value in g {
                *value = (*value / max) / scaled;
            }
        } else {
            if !(total > 0.0 && total.is_finite()) {
                return Err(invalid("a normalization group has no mass"));
            }
            if 1.0 - tolerance > total || 1.0 + tolerance < total {
                for value in g {
                    *value /= total;
                }
            }
        }
    }
    Ok(())
}

/// The same factors on the reverse strand: position `p` becomes
/// `width − 1 − p` and every symbol its complement.
pub(crate) fn reverse(s: &Smoothed, alphabet_len: usize, complement: &[u8]) -> Smoothed {
    let n = alphabet_len;
    let w = s.width as i32;
    let factors = s
        .factors
        .iter()
        .map(|f| {
            let k = f.positions.len();
            let positions: Vec<i32> = f.positions.iter().rev().map(|&p| w - 1 - p).collect();
            let mut table = vec![0.0; f.table.len()];
            let mut digits = vec![0usize; k];
            for (i, t) in table.iter_mut().enumerate() {
                decode(i, n, &mut digits);
                // New axis j is old axis k − 1 − j, complemented.
                let old = digits
                    .iter()
                    .rev()
                    .fold(0usize, |acc, &d| acc * n + complement[d] as usize);
                *t = f.table[old];
            }
            ProbFactor {
                kind: f.kind,
                positions,
                table,
            }
        })
        .collect();
    let covered = s.covered.iter().rev().copied().collect();
    Smoothed {
        width: s.width,
        factors,
        covered,
    }
}

fn invalid(detail: impl Into<String>) -> Error {
    Error::invalid("motif model", detail)
}
