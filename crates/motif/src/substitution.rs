//! Substitution-matrix priors: motif columns smoothed toward letters that
//! are likely substitutes of the observed ones.
//!
//! A column `f` over `n` letters with background `b` gets the pseudocount
//! frequencies
//!
//! ```text
//! g_i = Σ_j (f_j / b_j) · q_ij
//! ```
//!
//! where `q_ij = b_i · b_j · exp(λ · s_ij)` are the target frequencies of an
//! integer log-odds substitution matrix `s` under the background, and `λ` is
//! the Karlin–Altschul parameter solving `Σ_ij b_i b_j exp(λ s_ij) = 1`. The
//! smoothed column is `(α·f + β·g) / (α + β)`.
//!
//! MCAST (MEME Suite 5.5.9) applies this to every motif column of its HMM
//! with the DNA PAM-1 transition/transversion matrix, `α = 20`, `β = 1`,
//! and the background it scores against; [`SubstitutionPrior::mcast_dna`]
//! is that configuration.

use std::collections::BTreeMap;
use std::fmt::Debug;

use omicsweave_core::{Error, Result};

use crate::scoring::ColumnPrior;

/// One step of the DNA mutation process, `Pr(i | j)` per 10 000 in the
/// order `ACGT`: identity 9900, transitions (A↔G, C↔T) 60, transversions 20.
const DNA_PAM1: [[f64; 4]; 4] = [
    [9900.0, 20.0, 60.0, 20.0],
    [20.0, 9900.0, 20.0, 60.0],
    [60.0, 20.0, 9900.0, 20.0],
    [20.0, 60.0, 20.0, 9900.0],
];

/// Letter frequencies the PAM process is defined under.
const DNA_PAM_FREQ: [f64; 4] = [0.25; 4];

/// Keeps `log` finite for impossible substitutions.
const EPSILON: f64 = 1e-200;

/// Bisection steps of the λ search.
const LAMBDA_STEPS: usize = 25;

fn nint(x: f64) -> i32 {
    if x >= 0.0 {
        (x + 0.5) as i32
    } else {
        (x - 0.5) as i32
    }
}

/// The DNA PAM log-odds matrix at evolutionary distance `dist ≥ 1`, row
/// major over `ACGT`.
///
/// The one-step conditional matrix is made symmetric in its joint form
/// (`joint_ij = (M_ij p_j + M_ji p_i) / 2`), raised to the power `dist`, and
/// each cell becomes `nint(c · ln((M_ij + 1e-200) / p_i))` with
/// `c = 2 / ln 2` below distance 170 and `3 / ln 2` from there on (half or
/// third bits, as BLAST's PAM tables).
pub fn pam_dna_scores(dist: u32) -> Result<Vec<i32>> {
    if dist == 0 {
        return Err(Error::invalid("PAM matrix", "distance must be at least 1"));
    }
    let n = 4;
    let p = DNA_PAM_FREQ;
    let mut one = [[0.0f64; 4]; 4];
    for i in 0..n {
        for j in 0..=i {
            let joint = (DNA_PAM1[i][j] * p[j] + DNA_PAM1[j][i] * p[i]) / 20000.0;
            one[i][j] = joint / p[j];
            if i != j {
                one[j][i] = joint / p[i];
            }
        }
    }
    let mut m = one;
    if dist > 256 {
        let mut power = one;
        m = [[0.0; 4]; 4];
        for (i, row) in m.iter_mut().enumerate() {
            row[i] = 1.0;
        }
        let mut remaining = dist;
        while remaining > 0 {
            if remaining & 1 != 0 {
                m = multiply(m, power);
            }
            remaining >>= 1;
            if remaining > 0 {
                power = multiply(power, power);
            }
        }
    } else {
        for _ in 1..dist {
            let mut product = [[0.0f64; 4]; 4];
            for (i, row) in product.iter_mut().enumerate() {
                for (j, cell) in row.iter_mut().enumerate() {
                    let mut s = 0.0;
                    for k in 0..n {
                        s += m[i][k] * one[k][j];
                    }
                    *cell = s;
                }
            }
            m = product;
        }
    }
    let factor = if dist < 170 {
        2.0 / 2f64.ln()
    } else {
        3.0 / 2f64.ln()
    };
    let mut out = Vec::with_capacity(n * n);
    for (i, row) in m.iter().enumerate() {
        for &v in row {
            out.push(nint(factor * ((v + EPSILON) / p[i]).ln()));
        }
    }
    Ok(out)
}

fn multiply(a: [[f64; 4]; 4], b: [[f64; 4]; 4]) -> [[f64; 4]; 4] {
    let mut out = [[0.0; 4]; 4];
    for (i, row) in out.iter_mut().enumerate() {
        for (j, cell) in row.iter_mut().enumerate() {
            for (k, right) in b.iter().enumerate() {
                *cell += a[i][k] * right[j];
            }
        }
    }
    out
}

/// The Karlin–Altschul λ of an integer score matrix (`n × n`, row major)
/// whose letters are drawn independently from `letters`.
///
/// The score distribution `P(s) = Σ_{s_ij = s} b_i b_j` is normalized to sum
/// 1; it must reach a positive score and have a negative mean. λ is the
/// positive root of `Σ_s P(s) e^{λ s} = 1`: an upper bound is found by
/// doubling from 1 until the sum reaches 1, then 25 bisection steps from
/// `[0, upper]` keep the lower end, the returned value.
pub fn karlin_lambda(scores: &[i32], letters: &[f64]) -> Result<f64> {
    let n = letters.len();
    if n.checked_mul(n) != Some(scores.len()) || n == 0 {
        return Err(invalid("score matrix and letters disagree in size"));
    }
    if letters.iter().any(|p| !p.is_finite() || *p < 0.0) {
        return Err(invalid(
            "letter probabilities must be finite and nonnegative",
        ));
    }
    let low = scores.iter().copied().fold(1, i32::min) as i64;
    let high = scores.iter().copied().fold(-1, i32::max) as i64;
    if low >= 0 {
        return Err(invalid("the lowest score must be negative"));
    }
    if high <= 0 {
        return Err(invalid("the highest score must be positive"));
    }
    let range = (high - low) as usize;
    if range > 4096 {
        return sparse_lambda(scores, letters, high as f64);
    }
    let mut pr = vec![0.0f64; range + 1];
    for i in 0..n {
        for j in 0..n {
            let s = scores[i * n + j] as i64;
            pr[(s - low) as usize] += letters[i] * letters[j];
        }
    }
    if !pr[(-low) as usize + 1..].iter().any(|&v| v != 0.0) {
        return Err(invalid("a positive score must be possible"));
    }
    if pr.iter().any(|&v| v < 0.0) {
        return Err(invalid("negative letter probabilities"));
    }
    let mut sum = 0.0;
    for &v in &pr {
        sum += v;
    }
    if !(sum.is_finite() && sum > 0.0) {
        return Err(invalid("letter probabilities have no finite positive mass"));
    }
    let p: Vec<f64> = pr.iter().map(|v| v / sum).collect();
    let mut mean = low as f64;
    for (i, &v) in p.iter().enumerate() {
        mean += i as f64 * v;
    }
    if mean >= 0.0 {
        return Err(invalid("the expected score must be negative"));
    }
    let moment = |lambda: f64| -> f64 {
        let beta = lambda.exp();
        let mut term = (lambda * (low - 1) as f64).exp();
        let mut s = 0.0;
        for &v in &p {
            term *= beta;
            s += v * term;
        }
        if s == 0.0 || s.is_nan() {
            p.iter()
                .enumerate()
                .filter(|(_, v)| **v > 0.0)
                .map(|(i, &v)| v * (lambda * (low + i as i64) as f64).exp())
                .sum()
        } else {
            s
        }
    };
    lambda_root(moment, 0.5)
}

fn sparse_lambda(scores: &[i32], letters: &[f64], high: f64) -> Result<f64> {
    let n = letters.len();
    let mut mass = BTreeMap::<i32, f64>::new();
    for (i, &a) in letters.iter().enumerate() {
        for (j, &b) in letters.iter().enumerate() {
            *mass.entry(scores[i * n + j]).or_default() += a * b;
        }
    }
    let total: f64 = mass.values().sum();
    if !(total.is_finite() && total > 0.0) {
        return Err(invalid("letter probabilities have no finite positive mass"));
    }
    for value in mass.values_mut() {
        *value /= total;
    }
    if !mass.iter().any(|(&score, &p)| score > 0 && p > 0.0) {
        return Err(invalid("a positive score must be possible"));
    }
    if mass
        .iter()
        .map(|(&score, &p)| score as f64 * p)
        .sum::<f64>()
        >= 0.0
    {
        return Err(invalid("the expected score must be negative"));
    }
    let moment = |lambda: f64| {
        mass.iter()
            .filter(|(_, p)| **p > 0.0)
            .map(|(&score, &p)| p * (lambda * score as f64).exp())
            .sum::<f64>()
    };
    lambda_root(moment, 0.5 / high)
}

fn lambda_root(moment: impl Fn(f64) -> f64, mut up: f64) -> Result<f64> {
    loop {
        up *= 2.0;
        if moment(up) >= 1.0 {
            break;
        }
    }
    let mut lambda = 0.0f64;
    for _ in 0..LAMBDA_STEPS {
        let mid = (lambda + up) / 2.0;
        if moment(mid) > 1.0 {
            up = mid;
        } else {
            lambda = mid;
        }
    }
    Ok(lambda)
}

/// Target frequencies `q_ij = (b_i · b_j) · exp(λ · s_ij)`, row major.
pub fn target_frequencies(scores: &[i32], letters: &[f64], lambda: f64) -> Vec<f64> {
    let n = letters.len();
    let mut q = Vec::with_capacity(n * n);
    for i in 0..n {
        for j in 0..n {
            q.push(letters[i] * letters[j] * (lambda * scores[i * n + j] as f64).exp());
        }
    }
    q
}

/// A column prior built from a substitution matrix (see the module docs).
#[derive(Clone, Debug, PartialEq)]
pub struct SubstitutionPrior {
    /// Weight of the observed column.
    pub alpha: f64,
    /// Weight of the pseudocount frequencies.
    pub beta: f64,
    /// Target frequencies `q_ij`, row major over the alphabet.
    pub target: Vec<f64>,
}

impl SubstitutionPrior {
    /// A prior from an integer score matrix under the background `letters`.
    pub fn new(
        scores: &[i32],
        letters: &[f64],
        alpha: f64,
        beta: f64,
    ) -> Result<SubstitutionPrior> {
        if !(alpha.is_finite()
            && beta.is_finite()
            && (alpha + beta).is_finite()
            && alpha >= 0.0
            && beta >= 0.0
            && alpha + beta > 0.0)
        {
            return Err(invalid("weights must be >= 0 with a positive sum"));
        }
        let lambda = karlin_lambda(scores, letters)?;
        Ok(SubstitutionPrior {
            alpha,
            beta,
            target: target_frequencies(scores, letters, lambda),
        })
    }

    /// MCAST's prior: DNA PAM-1, `α = 20`, `β = 1`, target frequencies under
    /// `letters` (the prepared background over `ACGT`).
    pub fn mcast_dna(letters: &[f64]) -> Result<SubstitutionPrior> {
        if letters.len() != 4 {
            return Err(invalid("the DNA prior needs 4 letter frequencies"));
        }
        SubstitutionPrior::new(&pam_dna_scores(1)?, letters, 20.0, 1.0)
    }

    /// The pseudocount frequencies `g` of a column under `letters`.
    pub fn pseudocounts(&self, column: &[f64], letters: &[f64]) -> Vec<f64> {
        let n = column.len();
        (0..n)
            .map(|i| {
                let mut g = 0.0;
                for j in 0..n {
                    g += (column[j] / letters[j]) * self.target[i * n + j];
                }
                g
            })
            .collect()
    }

    /// `(α·f + β·g) / (α + β)` per letter.
    pub fn smooth(&self, column: &[f64], letters: &[f64]) -> Vec<f64> {
        let g = self.pseudocounts(column, letters);
        column
            .iter()
            .zip(&g)
            .map(|(&f, &g)| (self.alpha * f + self.beta * g) / (self.alpha + self.beta))
            .collect()
    }
}

impl ColumnPrior for SubstitutionPrior {
    fn smooth(&self, column: &[f64], _sites: f64, letters: &[f64]) -> Vec<f64> {
        if column.len() * column.len() != self.target.len() || letters.len() != column.len() {
            // A column over a reduced support has no substitution matrix;
            // it is kept as observed.
            return column.to_vec();
        }
        SubstitutionPrior::smooth(self, column, letters)
    }
}

fn invalid(detail: &str) -> Error {
    Error::invalid("substitution prior", detail)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn pam1_dna_scores() {
        assert_eq!(
            pam_dna_scores(1).unwrap(),
            vec![
                4, -14, -11, -14, -14, 4, -14, -11, -11, -14, 4, -14, -14, -11, -14, 4
            ]
        );
        assert!(pam_dna_scores(0).is_err());
        let far = pam_dna_scores(50).unwrap();
        assert!(far[0] < 4 && far[2] > far[1]);
    }

    #[test]
    fn lambda_solves_the_moment_equation() {
        let s = pam_dna_scores(1).unwrap();
        let b = [0.3, 0.2, 0.2, 0.3];
        let l = karlin_lambda(&s, &b).unwrap();
        let q = target_frequencies(&s, &b, l);
        let total: f64 = q.iter().sum();
        assert!((total - 1.0).abs() < 1e-6, "{total}");
        assert!(karlin_lambda(&[1, 1, 1, 1], &[0.5, 0.5]).is_err());
        assert!(karlin_lambda(&[-1, -1, -1, -1], &[0.5, 0.5]).is_err());
    }

    #[test]
    fn mcast_lambda_under_symmetrized_nrdb() {
        // NRDB, averaged with the complement, raised by 5e-7, normalized.
        let raw = [0.281774f64, 0.222020, 0.228876, 0.267330];
        let at = (raw[0] + raw[3]) / 2.0;
        let cg = (raw[1] + raw[2]) / 2.0;
        let f = [at, cg, cg, at];
        let mut sum = 5e-7 * 4.0;
        for v in f {
            sum += v;
        }
        let mut b: Vec<f64> = f.iter().map(|v| (v + 5e-7) / sum).collect();
        let total: f64 = b.iter().sum();
        if total != 1.0 {
            b.iter_mut().for_each(|v| *v /= total);
        }
        let l = karlin_lambda(&pam_dna_scores(1).unwrap(), &b).unwrap();
        assert!((l - 0.34167200326919556).abs() < 1e-15, "{l}");
    }

    #[test]
    fn smoothing_mixes_toward_substitutes() {
        let b = [0.25; 4];
        let prior = SubstitutionPrior::mcast_dna(&b).unwrap();
        let col = [1.0, 0.0, 0.0, 0.0];
        let s = prior.smooth(&col, &b);
        assert!((s.iter().sum::<f64>() - 1.0).abs() < 1e-6);
        // A's transition partner G gets more than the transversions.
        assert!(s[2] > s[1] && s[2] > s[3] && s[0] > 0.9);
        assert_eq!(
            ColumnPrior::smooth(&prior, &[0.5, 0.5], 10.0, &[0.5, 0.5]),
            vec![0.5, 0.5]
        );
    }
    #[test]
    fn lambda_handles_wide_scores_and_rejects_invalid_frequencies() {
        let letters = [0.5, 0.5];
        let lambda = karlin_lambda(&[-1000, 1, 1, -1000], &letters).unwrap();
        assert!((lambda - 2f64.ln()).abs() < 1e-6);
        let lambda = karlin_lambda(
            &[-2_000_000_000, 1_000_000_000, 1_000_000_000, -2_000_000_000],
            &letters,
        )
        .unwrap();
        let expected = ((1.0 + 5f64.sqrt()) / 2.0).ln();
        assert!((lambda * 1e9 - expected).abs() < 1e-6);
        assert!(karlin_lambda(&[-1, 1, 1, -1], &[f64::NAN, 0.5]).is_err());
        assert!(karlin_lambda(&[-1, 1, 1, -1], &[0.0, 0.0]).is_err());
        assert!(karlin_lambda(&[-1, 1, 1, -1], &[-0.5, -0.5]).is_err());
        assert!(
            SubstitutionPrior::new(&pam_dna_scores(1).unwrap(), &[0.25; 4], f64::INFINITY, 1.0)
                .is_err()
        );
    }

    #[test]
    fn distant_pam_models_reach_the_stationary_distribution() {
        assert_eq!(pam_dna_scores(1_000_000).unwrap(), vec![0; 16]);
        assert_eq!(pam_dna_scores(u32::MAX).unwrap(), vec![0; 16]);
    }
}
