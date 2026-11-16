//! Statistics over window p-values, as replaceable interfaces.
//!
//! The scan reports exact single-window statistics. How they become
//! sequence-level or collection-level statistics (the best of many windows,
//! a combination of several motifs, a multiple-testing adjustment, a tail
//! fitted to sampled scores) is a choice, so each step is a trait with a few
//! implementations; closures of the right shape implement the traits too.
//! MEME Suite's formulas are provided as one implementation each.

use omicsweave_core::{Error, Result};

/// The p-value of the best of `n` independent tests that each have p-value
/// `p` under the null.
pub trait BestOf: Send + Sync {
    fn pvalue(&self, p: f64, n: f64) -> f64;
}

impl<F: Fn(f64, f64) -> f64 + Send + Sync> BestOf for F {
    fn pvalue(&self, p: f64, n: f64) -> f64 {
        self(p, n)
    }
}

/// `1 − (1 − p)^n`, computed without cancellation.
#[derive(Clone, Copy, Debug, Default)]
pub struct Sidak;

impl BestOf for Sidak {
    fn pvalue(&self, p: f64, n: f64) -> f64 {
        if n == 0.0 {
            return 0.0;
        }
        if p >= 1.0 {
            return 1.0;
        }
        -f64::exp_m1(n * f64::ln_1p(-p))
    }
}

/// MEME's form: `p·n` when that is at most 0.001, `1 − (1 − p)^n` otherwise.
#[derive(Clone, Copy, Debug, Default)]
pub struct MemeBestOf;

impl BestOf for MemeBestOf {
    fn pvalue(&self, p: f64, n: f64) -> f64 {
        let small = p * n;
        if small <= 0.001 {
            small
        } else {
            1.0 - (1.0 - p).powf(n)
        }
    }
}

/// One p-value from several independent ones.
pub trait Combine: Send + Sync {
    fn combine(&self, pvalues: &[f64]) -> f64;
}

impl<F: Fn(&[f64]) -> f64 + Send + Sync> Combine for F {
    fn combine(&self, pvalues: &[f64]) -> f64 {
        self(pvalues)
    }
}

/// The p-value of the product of `k` independent uniform p-values (Fisher's
/// method): `x · Σ_{i<k} (−ln x)^i / i!` for the observed product `x`. The
/// product is taken in the given order and the series summed as MEME's
/// `qfast` does.
#[derive(Clone, Copy, Debug, Default)]
pub struct Product;

impl Combine for Product {
    fn combine(&self, pvalues: &[f64]) -> f64 {
        let mut x = 1.0;
        for &p in pvalues {
            x *= p;
        }
        if x < f64::MIN_POSITIVE && !pvalues.contains(&0.0) {
            return log_product_tail(pvalues.len(), pvalues.iter().map(|p| p.ln()).sum());
        }
        product_tail(pvalues.len(), x)
    }
}

/// The product tail when a nonzero product cannot be represented as f64.
fn log_product_tail(n: usize, log_x: f64) -> f64 {
    let log_t = (-log_x).ln();
    let mut term = log_x;
    let mut sum = term;
    for i in 1..n {
        term += log_t - (i as f64).ln();
        let hi = sum.max(term);
        let lo = sum.min(term);
        sum = hi + (lo - hi).exp().ln_1p();
    }
    sum.exp().min(1.0)
}

/// `P(U_1 ⋯ U_n ≤ x)` for independent uniforms.
pub fn product_tail(n: usize, x: f64) -> f64 {
    if n == 0 {
        return 1.0;
    }
    if x == 0.0 {
        return 0.0;
    }
    let mlnx = -x.ln();
    let mut term = x;
    let mut phi = x;
    for i in 1..n {
        term *= mlnx / i as f64;
        phi += term;
    }
    phi
}

/// Smallest p-value times the number of p-values (Bonferroni on the minimum).
#[derive(Clone, Copy, Debug, Default)]
pub struct MinBonferroni;

impl Combine for MinBonferroni {
    fn combine(&self, pvalues: &[f64]) -> f64 {
        let min = pvalues.iter().copied().fold(1.0f64, f64::min);
        (min * pvalues.len() as f64).min(1.0)
    }
}

/// Adjusted values for a set of p-values (E-values, q-values, ...), in the
/// same order.
pub trait Adjust: Send + Sync {
    fn adjust(&self, pvalues: &[f64]) -> Vec<f64>;
}

impl<F: Fn(&[f64]) -> Vec<f64> + Send + Sync> Adjust for F {
    fn adjust(&self, pvalues: &[f64]) -> Vec<f64> {
        self(pvalues)
    }
}

/// Expected number of tests at least as significant: `p · tests`. `None`
/// counts the p-values given.
#[derive(Clone, Copy, Debug, Default)]
pub struct EValue {
    pub tests: Option<f64>,
}

impl Adjust for EValue {
    fn adjust(&self, pvalues: &[f64]) -> Vec<f64> {
        let n = self.tests.unwrap_or(pvalues.len() as f64);
        pvalues.iter().map(|p| p * n).collect()
    }
}

/// Benjamini–Hochberg q-values: `min_{j ≥ i} p_(j) · m / j`, at most 1.
#[derive(Clone, Copy, Debug, Default)]
pub struct BenjaminiHochberg {
    /// Number of tests when more were made than p-values given.
    pub tests: Option<f64>,
}

impl Adjust for BenjaminiHochberg {
    fn adjust(&self, pvalues: &[f64]) -> Vec<f64> {
        let m = self.tests.unwrap_or(pvalues.len() as f64);
        let mut order: Vec<usize> = (0..pvalues.len()).collect();
        order.sort_by(|&a, &b| pvalues[a].total_cmp(&pvalues[b]));
        let mut out = vec![0.0; pvalues.len()];
        let mut running = 1.0f64;
        for (rank, &i) in order.iter().enumerate().rev() {
            running = running.min(pvalues[i] * m / (rank + 1) as f64);
            out[i] = running;
        }
        out
    }
}

/// A right tail fitted to sampled scores: maps a score to the probability
/// that a null score is at least as large.
pub trait TailFit: Send + Sync {
    fn fit(&self, samples: &[f64]) -> Result<Box<dyn Fn(f64) -> f64 + Send + Sync>>;
}

/// The empirical tail: `(1 + #{samples ≥ x}) / (1 + n)`.
#[derive(Clone, Copy, Debug, Default)]
pub struct EmpiricalTail;

impl TailFit for EmpiricalTail {
    fn fit(&self, samples: &[f64]) -> Result<Box<dyn Fn(f64) -> f64 + Send + Sync>> {
        if samples.is_empty() {
            return Err(Error::invalid("tail fit", "no samples"));
        }
        let mut sorted = samples.to_vec();
        sorted.sort_by(f64::total_cmp);
        let n = sorted.len() as f64;
        Ok(Box::new(move |x| {
            let below = sorted.partition_point(|&s| s < x);
            (1.0 + (sorted.len() - below) as f64) / (1.0 + n)
        }))
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn product_tail_matches_closed_forms() {
        assert_eq!(product_tail(1, 0.3), 0.3);
        let x: f64 = 0.01;
        let two = x * (1.0 - x.ln());
        assert!((product_tail(2, x) - two).abs() < 1e-15);
        assert_eq!(Product.combine(&[0.5, 0.02]), product_tail(2, 0.01));
        assert_eq!(product_tail(0, 0.5), 1.0);
    }

    #[test]
    fn best_of_n() {
        assert_eq!(Sidak.pvalue(1.0, 0.0), 0.0);
        assert!((Sidak.pvalue(1e-6, 1000.0) - (1.0 - (1.0f64 - 1e-6).powi(1000))).abs() < 1e-12);
        assert_eq!(MemeBestOf.pvalue(1e-7, 100.0), 1e-7 * 100.0);
        assert!((MemeBestOf.pvalue(0.01, 10.0) - (1.0 - 0.99f64.powf(10.0))).abs() < 1e-15);
        let closure = |p: f64, n: f64| (p * n).min(1.0);
        assert_eq!(BestOf::pvalue(&closure, 0.1, 20.0), 1.0);
    }

    #[test]
    fn product_tail_preserves_probability_when_the_product_underflows() {
        let probabilities = vec![0.5; 2000];
        assert_eq!(probabilities.iter().product::<f64>(), 0.0);
        let p = Product.combine(&probabilities);
        assert!(p > 0.999 && p <= 1.0, "{p}");
        let mut with_zero = probabilities;
        with_zero[1000] = 0.0;
        assert_eq!(Product.combine(&with_zero), 0.0);
        let probabilities = vec![(-1.0f64).exp(); 744];
        let p = Product.combine(&probabilities);
        assert!((p - 0.495_124_652_497_848_3).abs() < 1e-10, "{p}");
        assert!(Product.combine(&vec![0.5; 1074]) > 0.999_999);
    }

    #[test]
    fn bh_is_monotone() {
        let q = BenjaminiHochberg::default().adjust(&[0.01, 0.04, 0.03, 0.5]);
        let t = 0.04 * 4.0 / 3.0;
        assert_eq!(q, vec![0.04, t, t, 0.5]);
    }
}
