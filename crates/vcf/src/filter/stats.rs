//! The statistics filter expressions call: htslib's `kt_fisher_exact` and
//! `kf_betai` (kfunc.c) and bcftools' `calc_binom_two_sided`, with the
//! same arithmetic in the same order.

unsafe extern "C" {
    /// The C library's log-gamma, which htslib's `lbinom` uses.
    fn lgamma_r(x: f64, sign: *mut std::ffi::c_int) -> f64;
}

fn c_lgamma(x: f64) -> f64 {
    let mut sign = 0;
    // SAFETY: libm writes only to the live, aligned local c_int. Its mutable
    // pointer is unaliased for the call and does not escape. The reentrant
    // entry point avoids shared signgam state between concurrent filters.
    unsafe { lgamma_r(x, &mut sign) }
}

const KF_GAMMA_EPS: f64 = 1e-14;
const KF_TINY: f64 = 1e-290;

// The constants are kfunc.c's, digit for digit.
#[allow(clippy::excessive_precision)]
fn kf_lgamma(z: f64) -> f64 {
    let mut x = 0.0;
    x += 0.1659470187408462e-06 / (z + 7.0);
    x += 0.9934937113930748e-05 / (z + 6.0);
    x -= 0.1385710331296526 / (z + 5.0);
    x += 12.50734324009056 / (z + 4.0);
    x -= 176.6150291498386 / (z + 3.0);
    x += 771.3234287757674 / (z + 2.0);
    x -= 1259.139216722289 / (z + 1.0);
    x += 676.5203681218835 / z;
    x += 0.9999999999995183;
    x.ln() - 5.58106146679532777 - z + (z - 0.5) * (z + 6.5).ln()
}

fn kf_betai_aux(a: f64, b: f64, x: f64) -> f64 {
    if x == 0.0 {
        return 0.0;
    }
    if x == 1.0 {
        return 1.0;
    }
    let mut f = 1.0;
    let mut c = f;
    let mut d = 0.0;
    for j in 1..200 {
        let m = f64::from(j >> 1);
        let aa = if j & 1 != 0 {
            -(a + m) * (a + b + m) * x / ((a + 2.0 * m) * (a + 2.0 * m + 1.0))
        } else {
            m * (b - m) * x / ((a + 2.0 * m - 1.0) * (a + 2.0 * m))
        };
        d = 1.0 + aa * d;
        if d < KF_TINY {
            d = KF_TINY;
        }
        c = 1.0 + aa / c;
        if c < KF_TINY {
            c = KF_TINY;
        }
        d = 1.0 / d;
        let dd = c * d;
        f *= dd;
        if (dd - 1.0).abs() < KF_GAMMA_EPS {
            break;
        }
    }
    (kf_lgamma(a + b) - kf_lgamma(a) - kf_lgamma(b) + a * x.ln() + b * (1.0 - x).ln()).exp() / a / f
}

/// Regularised incomplete beta function (`kf_betai`).
pub(super) fn kf_betai(a: f64, b: f64, x: f64) -> f64 {
    if x < (a + 1.0) / (a + b + 2.0) {
        kf_betai_aux(a, b, x)
    } else {
        1.0 - kf_betai_aux(b, a, 1.0 - x)
    }
}

/// bcftools' two-sided binomial test of `na` against `nb` successes.
pub(super) fn calc_binom_two_sided(na: i32, nb: i32, aprob: f64) -> f64 {
    if na == 0 && nb == 0 {
        return -1.0;
    }
    if na == nb {
        return 1.0;
    }
    let prob = if na > nb {
        2.0 * kf_betai(f64::from(na), f64::from(nb) + 1.0, aprob)
    } else {
        2.0 * kf_betai(f64::from(nb), f64::from(na) + 1.0, aprob)
    };
    if prob > 1.0 { 1.0 } else { prob }
}

// C `int` arithmetic below wraps on overflow, as x86-64 code does.

fn lbinom(n: i32, k: i32) -> f64 {
    if k == 0 || n == k {
        return 0.0;
    }
    c_lgamma(f64::from(n) + 1.0)
        - c_lgamma(f64::from(k) + 1.0)
        - c_lgamma(f64::from(n.wrapping_sub(k)) + 1.0)
}

fn hypergeo(n11: i32, n1_: i32, n_1: i32, n: i32) -> f64 {
    (lbinom(n1_, n11) + lbinom(n.wrapping_sub(n1_), n_1.wrapping_sub(n11)) - lbinom(n, n_1)).exp()
}

#[derive(Default)]
struct HgAcc {
    n11: i32,
    n1_: i32,
    n_1: i32,
    n: i32,
    p: f64,
}

fn hypergeo_acc(n11: i32, n1_: i32, n_1: i32, n: i32, aux: &mut HgAcc) -> f64 {
    if n1_ != 0 || n_1 != 0 || n != 0 {
        aux.n11 = n11;
        aux.n1_ = n1_;
        aux.n_1 = n_1;
        aux.n = n;
    } else {
        // n11 + n - n1_ - n_1
        let rest = n11
            .wrapping_add(aux.n)
            .wrapping_sub(aux.n1_)
            .wrapping_sub(aux.n_1);
        if n11 % 11 != 0 && rest != 0 {
            if n11 == aux.n11.wrapping_add(1) {
                aux.p *= f64::from(aux.n1_.wrapping_sub(aux.n11)) / f64::from(n11)
                    * f64::from(aux.n_1.wrapping_sub(aux.n11))
                    / f64::from(rest);
                aux.n11 = n11;
                return aux.p;
            }
            if n11 == aux.n11.wrapping_sub(1) {
                aux.p *= f64::from(aux.n11) / f64::from(aux.n1_.wrapping_sub(n11))
                    * f64::from(
                        aux.n11
                            .wrapping_add(aux.n)
                            .wrapping_sub(aux.n1_)
                            .wrapping_sub(aux.n_1),
                    )
                    / f64::from(aux.n_1.wrapping_sub(n11));
                aux.n11 = n11;
                return aux.p;
            }
        }
        aux.n11 = n11;
    }
    aux.p = hypergeo(aux.n11, aux.n1_, aux.n_1, aux.n);
    aux.p
}

/// Fisher's exact test of a 2x2 table: the two-sided p-value
/// (`kt_fisher_exact`'s `*two`). `None` for tables whose tails would take
/// longer to sum than any test can wait (counts in the hundreds of
/// millions).
pub(super) fn kt_fisher_two(n11: i32, n12: i32, n21: i32, n22: i32) -> Option<f64> {
    let n1_ = n11.wrapping_add(n12);
    let n_1 = n11.wrapping_add(n21);
    let n = n11.wrapping_add(n12).wrapping_add(n21).wrapping_add(n22);
    let max = if n_1 < n1_ { n_1 } else { n1_ };
    let mut min = n1_.wrapping_add(n_1).wrapping_sub(n);
    if min < 0 {
        min = 0;
    }
    if min == max {
        return Some(1.0);
    }
    if i64::from(max) - i64::from(min) > 50_000_000 {
        return None;
    }
    let mut aux = HgAcc::default();
    let q = hypergeo_acc(n11, n1_, n_1, n, &mut aux);
    if q == 0.0 {
        return Some(0.0);
    }
    let mut p = hypergeo_acc(min, 0, 0, 0, &mut aux);
    let mut left = 0.0;
    let mut i = min.wrapping_add(1);
    while p < 0.99999999 * q && i <= max {
        left += p;
        p = hypergeo_acc(i, 0, 0, 0, &mut aux);
        i = i.wrapping_add(1);
    }
    if p < 1.00000001 * q {
        left += p;
    }
    p = hypergeo_acc(max, 0, 0, 0, &mut aux);
    let mut right = 0.0;
    let mut j = max.wrapping_sub(1);
    while p < 0.99999999 * q && j >= 0 {
        right += p;
        p = hypergeo_acc(j, 0, 0, 0, &mut aux);
        j = j.wrapping_sub(1);
    }
    if p < 1.00000001 * q {
        right += p;
    }
    let two = left + right;
    Some(if two > 1.0 { 1.0 } else { two })
}
