//! Gapped profiles (FFM motifs with gap counts) scored by local alignment.
//!
//! A gapped profile is an alignment summarized per column: residue counts
//! (the motif's factors), deletion counts and the residues inserted after
//! each column ([`omicsweave_ffm::Gaps`]). It is scored the way GLAM2's
//! scanner (`glam2scan`, MEME 5.5.9) scores it:
//!
//! - **scores** ([`ProfileModel`]): column `i` with `m` residues, `d`
//!   deletions and `r` inserted residues after it (`a = m + d` aligned
//!   sequences) gets Beta posterior-mean rates
//!   `ds = (d + D) / (d + m + D + E)`, `ms = (m + E) / (d + m + D + E)`,
//!   `is = (r + I) / (r + a + I + J)`, `js = (a + J) / (r + a + I + J)` and
//!   residue probabilities `p` from a [`DirichletMixture`] posterior mean.
//!   Matching residue `x` scores `ln(js · ms · p_x / b_x)`, deleting the
//!   column `ln(ds · js)` and each residue inserted after it `ln(is)`; the
//!   last column has no `js` factor and admits no insertion. Symbols outside
//!   the alphabet (the wildcard) never match but may be inserted.
//! - **alignment**: local in the sequence, global in the profile — the best
//!   score of the columns `0 ..= i` ending at sequence position `j`, with
//!   free start anywhere and leading deletions allowed.
//! - **hits**: the best cell of the last row is traced back (preferring
//!   match, then deletion, then insertion on ties), its matched cells are
//!   forbidden, the affected cells are recomputed, and the next best is
//!   traced, until an alignment matches no residue. Hits may overlap; only
//!   identical (column, residue) pairs are excluded.
//! - **reverse strand**: the reverse complement is aligned with the same
//!   scores.
//!
//! [`ProfileProgram`] compiles profiles; its scan hands each sequence's hit
//! streams, in sequence order, to an [`AlignmentReducer`]. [`TopAlignments`]
//! keeps GLAM2's best `n` per motif and [`write_glam2_text`] prints them as
//! `glam2scan -t` does.

use std::cmp::Ordering;
use std::collections::BTreeMap;
use std::io::{self, Write};
use std::sync::atomic::{AtomicU64, AtomicUsize, Ordering as AtomicOrdering};
use std::sync::{Condvar, Mutex, mpsc};

use omicsweave_core::{Error, Result, Strand, Strings};
use omicsweave_ffm::{Alphabet, INVALID, MotifView};

use crate::{Background, ScanOptions, Strands};

/// Hits `glam2scan` reports by default (`-n`).
pub const GLAM2_HITS: usize = 25;

/// Width of the sequence-name field in GLAM2's text output.
const GLAM2_NAME_WIDTH: usize = 24;

const PROTEIN_SYMBOLS: &[u8] = b"ACDEFGHIKLMNPQRSTVWY";

/// Amino-acid abundances of Robinson & Robinson (1991), the background of
/// GLAM2 (and BLAST) for proteins, in [`PROTEIN_SYMBOLS`] order.
const ROBINSON_COUNTS: [f64; 20] = [
    35155.0, 8669.0, 24161.0, 28354.0, 17367.0, 33229.0, 9906.0, 23161.0, 25872.0, 40625.0,
    10101.0, 20212.0, 23435.0, 19208.0, 23105.0, 32070.0, 26311.0, 29012.0, 5990.0, 14488.0,
];

/// The published 20-component `recode3` Dirichlet mixture for amino acids
/// (columns in [`PROTEIN_SYMBOLS`] order), GLAM2's protein prior.
const RECODE3_WEIGHTS: [f64; 20] = [
    0.176513, 0.207622, 0.0669246, 0.0868259, 0.0593123, 0.0358616, 0.03427, 0.0428319, 0.047875,
    0.0466614, 0.0283695, 0.0301127, 0.0233828, 0.034662, 0.0270202, 0.0226822, 0.00898452,
    0.00716226, 0.00710292, 0.00582299,
];

const RECODE3_ALPHAS: [[f64; 20]; 20] = [
    [
        0.668781, 0.072626, 1.08101, 0.930936, 0.1235, 0.982387, 0.366701, 0.0586812, 0.849019,
        0.185073, 0.0648062, 0.979541, 0.519704, 0.555784, 0.636552, 1.04987, 0.60931, 0.129148,
        0.0676805, 0.201772,
    ],
    [
        0.880009, 0.197393, 0.303744, 0.497357, 0.657492, 0.391988, 0.348159, 0.936094, 0.527206,
        1.42844, 0.449302, 0.368825, 0.384576, 0.439775, 0.581516, 0.622624, 0.747064, 1.08007,
        0.235012, 0.606928,
    ],
    [
        0.153384, 0.0520756, 0.0073824, 0.0158439, 0.428964, 0.025533, 0.0185789, 0.845361,
        0.0282996, 2.42256, 0.424296, 0.0190716, 0.0313429, 0.0274578, 0.0252186, 0.028514,
        0.0519217, 0.522946, 0.0279653, 0.0664755,
    ],
    [
        0.794007, 0.0130659, 0.624236, 1.85769, 0.0290214, 0.115707, 0.123504, 0.22099, 1.52605,
        0.341371, 0.111114, 0.308302, 0.263545, 0.953727, 0.933444, 0.554741, 0.604551, 0.396451,
        0.00823516, 0.0420054,
    ],
    [
        0.740015, 0.187165, 0.0213261, 0.0456854, 0.118944, 0.0633687, 0.0170331, 1.06684,
        0.0380614, 0.733524, 0.138456, 0.0300644, 0.0718692, 0.0240143, 0.0301022, 0.0862989,
        0.367283, 1.70735, 0.0113856, 0.045079,
    ],
    [
        0.15978, 0.0261585, 0.0505181, 0.125524, 0.0350331, 0.102549, 0.157461, 0.0795041, 1.26261,
        0.189383, 0.0550608, 0.171028, 0.0844169, 0.290476, 1.44604, 0.129158, 0.138972, 0.0851144,
        0.0159134, 0.0637679,
    ],
    [
        0.308434, 0.0137217, 1.69731, 1.92422, 0.0361113, 0.162357, 0.07232, 0.0487895, 0.236135,
        0.0809074, 0.0286236, 0.213663, 0.181631, 0.320245, 0.104878, 0.218398, 0.141668,
        0.0747719, 0.0141705, 0.0453433,
    ],
    [
        0.00260287,
        9.99856e-06,
        0.00631292,
        0.00445502,
        0.00274753,
        1.03886e-05,
        1.02839e-05,
        0.000913052,
        0.0029241,
        0.00353485,
        0.00105128,
        0.00338172,
        1.04172e-05,
        0.00173574,
        0.00459583,
        0.00274255,
        0.00247625,
        0.00175366,
        1.02411e-05,
        0.00288489,
    ],
    [
        1.61043, 0.15522, 0.0378292, 0.0498243, 0.0406484, 0.529136, 0.0217524, 0.040597,
        0.0413396, 0.100193, 0.0509779, 0.0357917, 0.0931204, 0.0367156, 0.0330646, 0.529587,
        0.196607, 0.230878, 0.00909518, 0.0329275,
    ],
    [
        0.15525, 0.0136827, 0.0857138, 0.0508316, 0.0151451, 3.10555, 0.027169, 0.0140491,
        0.0654038, 0.0257501, 0.00901049, 0.127437, 0.0423873, 0.0345064, 0.0477247, 0.12452,
        0.0341196, 0.0230637, 0.00930115, 0.0187464,
    ],
    [
        0.225739, 0.0684326, 0.101072, 0.0813791, 0.0298832, 0.0915218, 0.0336807, 0.0833114,
        0.0931673, 0.0731542, 0.0419314, 0.230216, 0.087446, 0.0694702, 0.0751969, 1.13857,
        1.63158, 0.179083, 0.00912576, 0.0311963,
    ],
    [
        1.3431e-06,
        0.0166656,
        0.00743068,
        1.34592e-06,
        0.169076,
        0.00407061,
        0.00714122,
        1.98221,
        0.017522,
        0.816669,
        0.114773,
        0.00678027,
        0.0106392,
        0.0100244,
        0.0158968,
        0.00879658,
        0.0399043,
        1.81043,
        0.0150671,
        0.0517801,
    ],
    [
        0.063525, 0.0288391, 1.09265, 0.0959581, 0.00965196, 0.216914, 0.0730986, 0.0207325,
        0.08719, 0.0315107, 0.00960293, 0.752755, 0.059914, 0.0445321, 0.0312317, 0.287327,
        0.116896, 0.0249446, 0.00701663, 0.0305859,
    ],
    [
        0.294281, 0.019271, 0.12293, 0.162747, 0.0373667, 0.145029, 0.0412349, 0.0815261, 0.157594,
        0.151631, 0.021412, 0.0601581, 3.6966, 0.0809085, 0.101856, 0.23533, 0.135424, 0.140532,
        0.00900473, 0.0321389,
    ],
    [
        0.0844832, 0.0584945, 0.0411628, 0.045719, 0.847822, 0.0590839, 0.250253, 0.0675757,
        0.0562614, 0.168617, 0.0439737, 0.0794234, 0.028301, 0.0305672, 0.0598024, 0.0798202,
        0.0585385, 0.0858243, 0.227395, 1.30336,
    ],
    [
        0.0634034,
        0.0246167,
        1.3443e-06,
        0.00389272,
        1.12953,
        0.00796028,
        1.35032e-06,
        0.233395,
        1.34466e-06,
        0.541933,
        0.101309,
        1.36412e-06,
        0.027467,
        0.00704479,
        0.00802297,
        0.0248977,
        0.0276933,
        0.185467,
        0.183309,
        0.516892,
    ],
    [
        0.123696,
        0.0454619,
        0.0386434,
        0.351847,
        0.0560181,
        0.0439442,
        0.223229,
        0.01302,
        0.148699,
        0.19001,
        0.120964,
        0.098734,
        5.90055e-06,
        0.554971,
        0.219233,
        0.0453885,
        0.0564686,
        0.0614792,
        0.0410248,
        0.0800036,
    ],
    [
        0.0212037,
        3.18769,
        0.00745627,
        0.00382411,
        0.00691924,
        0.0126233,
        1.34375e-06,
        0.00724293,
        0.00522979,
        0.00785563,
        0.00489521,
        0.0105326,
        0.0136265,
        0.00505819,
        0.00677712,
        0.0251744,
        0.0235516,
        0.0371462,
        0.00187667,
        0.0038893,
    ],
    [
        0.0229376, 0.00427768, 0.00959934, 0.013608, 0.182277, 0.0227654, 0.0157344, 0.0226783,
        0.011561, 0.0803491, 0.0154283, 0.00899225, 0.00980608, 0.00600945, 0.0342359, 0.0216842,
        0.0189306, 0.0223176, 1.83914, 0.154565,
    ],
    [
        2.16602e-06,
        2.16245e-06,
        0.0198496,
        2.17942e-06,
        0.0246741,
        2.47051e-06,
        1.02563,
        0.0131152,
        2.16539e-06,
        0.00637704,
        2.1414e-06,
        0.0839371,
        0.0168135,
        0.0438887,
        0.0252951,
        0.0235533,
        0.0130626,
        0.00797507,
        2.16433e-06,
        0.0545531,
    ],
];

/// What GLAM2's built-in presets recognize an alphabet as.
#[derive(Clone, Copy, PartialEq, Eq)]
enum Kind {
    Nucleotide,
    Protein,
    Other,
}

fn kind(alphabet: &Alphabet) -> Kind {
    match alphabet.symbols() {
        b"ACGT" | b"ACGU" => Kind::Nucleotide,
        s if s == PROTEIN_SYMBOLS => Kind::Protein,
        _ => Kind::Other,
    }
}

fn invalid(what: &'static str, detail: impl Into<String>) -> Error {
    Error::invalid(what, detail)
}

/// DNA as GLAM2 reads it: [`Alphabet::dna`] with `U`/`u` read as `T`.
pub fn glam2_dna_alphabet() -> Alphabet {
    Alphabet::dna()
        .with_alias(b'U', b'T')
        .and_then(|a| a.with_alias(b'u', b'T'))
        .expect("valid preset")
}

/// GLAM2's own background: uniform for nucleotides, Robinson & Robinson
/// abundances for the 20 amino acids, uniform for any other alphabet.
pub fn glam2_background(alphabet: &Alphabet) -> Result<Background> {
    match kind(alphabet) {
        Kind::Protein => Background::frequencies(&ROBINSON_COUNTS),
        _ => Background::frequencies(&vec![1.0; alphabet.len()]),
    }
}

/// A Dirichlet mixture prior over residue probabilities: component `c` has
/// weight `weights[c]` and pseudocounts `alphas[c]` (one per symbol).
#[derive(Clone, Debug, PartialEq)]
pub struct DirichletMixture {
    pub weights: Vec<f64>,
    pub alphas: Vec<Vec<f64>>,
}

impl DirichletMixture {
    /// A checked mixture: at least one component, positive finite values,
    /// every component of the same dimension. Weights need not sum to 1.
    pub fn new(weights: Vec<f64>, alphas: Vec<Vec<f64>>) -> Result<DirichletMixture> {
        let m = DirichletMixture { weights, alphas };
        m.check(None)?;
        Ok(m)
    }

    /// One component with every pseudocount 1.
    pub fn uniform(alphabet_len: usize) -> DirichletMixture {
        DirichletMixture {
            weights: vec![1.0],
            alphas: vec![vec![1.0; alphabet_len]],
        }
    }

    /// GLAM2's nucleotide prior: one component, every pseudocount 0.4.
    pub fn glam2_dna() -> DirichletMixture {
        DirichletMixture {
            weights: vec![1.0],
            alphas: vec![vec![0.4; 4]],
        }
    }

    /// GLAM2's protein prior, the 20-component `recode3` mixture.
    pub fn glam2_protein() -> DirichletMixture {
        DirichletMixture {
            weights: RECODE3_WEIGHTS.to_vec(),
            alphas: RECODE3_ALPHAS.iter().map(|a| a.to_vec()).collect(),
        }
    }

    /// GLAM2's prior for an alphabet: [`Self::glam2_dna`] for nucleotides,
    /// [`Self::glam2_protein`] for the 20 amino acids, otherwise
    /// [`Self::uniform`].
    pub fn glam2(alphabet: &Alphabet) -> DirichletMixture {
        match kind(alphabet) {
            Kind::Nucleotide => DirichletMixture::glam2_dna(),
            Kind::Protein => DirichletMixture::glam2_protein(),
            Kind::Other => DirichletMixture::uniform(alphabet.len()),
        }
    }

    /// Parses a mixture file as `glam2scan -d` reads it: whitespace
    /// separated words; `Mixture=` is followed by a component weight,
    /// `Alpha=` by the sum of the pseudocounts (ignored) and then the
    /// pseudocounts, read while words are numbers. The i-th weight belongs
    /// to the i-th pseudocount list; other words are ignored.
    pub fn parse_glam2(text: &str, alphabet_len: usize) -> Result<DirichletMixture> {
        const WHAT: &str = "Dirichlet mixture file";
        let mut weights = Vec::new();
        let mut alphas: Vec<Vec<f64>> = Vec::new();
        let mut words = text.split_whitespace().peekable();
        let number = |w: &str| w.parse::<f64>().ok().filter(|x| !x.is_nan());
        while let Some(word) = words.next() {
            match word {
                "Mixture=" => {
                    let x = words
                        .next()
                        .and_then(number)
                        .ok_or_else(|| invalid(WHAT, "Mixture= is not followed by a number"))?;
                    if !(x > 0.0 && x.is_finite()) {
                        return Err(invalid(WHAT, format!("Mixture= value {x} is not > 0")));
                    }
                    weights.push(x);
                }
                "Alpha=" => {
                    words
                        .next()
                        .and_then(number)
                        .ok_or_else(|| invalid(WHAT, "Alpha= is not followed by the sum"))?;
                    let mut alpha = Vec::new();
                    while let Some(x) = words.peek().and_then(|w| number(w)) {
                        words.next();
                        if !(x > 0.0 && x.is_finite()) {
                            return Err(invalid(WHAT, format!("Alpha= value {x} is not > 0")));
                        }
                        alpha.push(x);
                    }
                    if alpha.is_empty() {
                        return Err(invalid(WHAT, "Alpha= has no pseudocounts"));
                    }
                    alphas.push(alpha);
                }
                _ => {}
            }
        }
        if weights.len() != alphas.len() {
            return Err(invalid(
                WHAT,
                format!(
                    "{} Mixture= and {} Alpha= entries",
                    weights.len(),
                    alphas.len()
                ),
            ));
        }
        let m = DirichletMixture { weights, alphas };
        m.check(Some(alphabet_len))?;
        Ok(m)
    }

    /// Symbols per component.
    pub fn dimension(&self) -> usize {
        self.alphas.first().map_or(0, Vec::len)
    }

    fn check(&self, alphabet_len: Option<usize>) -> Result<()> {
        const WHAT: &str = "Dirichlet mixture";
        if self.weights.is_empty() || self.weights.len() != self.alphas.len() {
            return Err(invalid(
                WHAT,
                format!(
                    "{} weights for {} components",
                    self.weights.len(),
                    self.alphas.len()
                ),
            ));
        }
        let dim = self.dimension();
        if dim == 0 || self.alphas.iter().any(|a| a.len() != dim) {
            return Err(invalid(WHAT, "components differ in dimension or are empty"));
        }
        if let Some(n) = alphabet_len
            && n != dim
        {
            return Err(invalid(WHAT, format!("{dim} symbols, alphabet has {n}")));
        }
        let positive = |v: &f64| v.is_finite() && *v > 0.0;
        if !self.weights.iter().all(positive) || !self.alphas.iter().flatten().all(positive) {
            return Err(invalid(
                WHAT,
                "weights and pseudocounts must be finite and > 0",
            ));
        }
        if self
            .alphas
            .iter()
            .any(|alpha| !alpha.iter().sum::<f64>().is_finite())
        {
            return Err(invalid(WHAT, "a component's pseudocount sum is not finite"));
        }
        Ok(())
    }

    /// Posterior mean residue probabilities after observing `counts`: each
    /// component's `(n_x + α_x) / (N + Σα)`, weighted by the component's
    /// posterior probability.
    ///
    /// Panics if the mixture has no components or their dimensions differ
    /// from `counts`.
    pub fn posterior_mean(&self, counts: &[u32]) -> Vec<f64> {
        assert!(
            !self.weights.is_empty()
                && self.weights.len() == self.alphas.len()
                && self.alphas.iter().all(|alpha| alpha.len() == counts.len()),
            "mixture components and counts must have matching dimensions"
        );
        let total: u64 = counts.iter().map(|&c| c as u64).sum();
        if self.weights.len() == 1 {
            let alpha = &self.alphas[0];
            let sum: f64 = alpha.iter().sum();
            // A single component has posterior probability one, so its
            // likelihood need not be evaluated at all.
            let part = 1.0 / (total as f64 + sum);
            return counts
                .iter()
                .zip(alpha)
                .map(|(&n, &a)| {
                    if part.is_finite() {
                        part * (n as f64 + a)
                    } else {
                        (n as f64 + a) / (total as f64 + sum)
                    }
                })
                .collect();
        }
        // Per component: log of its weight times the Dirichlet-multinomial
        // likelihood of the counts, up to a factor shared by all components.
        let mut log_post = Vec::with_capacity(self.weights.len());
        let mut sums = Vec::with_capacity(self.weights.len());
        for (&w, alpha) in self.weights.iter().zip(&self.alphas) {
            let mut sum = 0.0;
            for &a in alpha {
                sum += a;
            }
            let mut rising = 0.0;
            for (&a, &n) in alpha.iter().zip(counts) {
                rising += log_rising(a, n as u64);
            }
            log_post.push((w.ln() + rising) - log_rising(sum, total));
            sums.push(sum);
        }
        let mut max = log_post[0];
        for &x in &log_post {
            if x > max {
                max = x;
            }
        }
        let mut norm = 0.0;
        for &x in &log_post {
            norm += (x - max).exp();
        }
        let norm_log = norm.ln() + max;
        let mut p = vec![0.0; counts.len()];
        for ((&lp, &sum), alpha) in log_post.iter().zip(&sums).zip(&self.alphas) {
            let weight = if total > 256 {
                // Center before normalizing: adding a small log normalizer
                // to a large log likelihood would discard its low bits.
                (lp - max).exp() / norm
            } else {
                (lp - norm_log).exp()
            };
            let part = weight / (total as f64 + sum);
            for ((pj, &n), &a) in p.iter_mut().zip(counts).zip(alpha) {
                *pj += if part.is_finite() && (part != 0.0 || weight == 0.0) {
                    part * (n as f64 + a)
                } else {
                    weight * ((n as f64 + a) / (total as f64 + sum))
                };
            }
        }
        p
    }
}

/// `ln Γ(a + k) − ln Γ(a)` for positive finite a.
fn log_rising(mut a: f64, k: u64) -> f64 {
    if k <= 256 {
        let mut t = 0.0;
        for i in 0..k {
            t += (i as f64 + a).ln();
        }
        return t;
    }
    let k = k as f64;
    let mut shift = 0.0;
    // Above 32, the first omitted Stirling correction is below 2e-14.
    while a < 32.0 {
        shift += a.ln() - (a + k).ln();
        a += 1.0;
    }
    let correction = |x: f64| {
        let inv = x.recip();
        let square = inv * inv;
        inv * (1.0 / 12.0 + square * (-1.0 / 360.0 + square / 1260.0))
    };
    // This ratio form avoids subtracting two nearly equal log-gamma values
    // when the pseudocount is much larger than the observed count.
    shift + k * a.ln() + (a + k - 0.5) * (k / a).ln_1p() - k + correction(a + k) - correction(a)
}

/// How gap and residue counts become alignment scores.
#[derive(Clone, Debug, PartialEq)]
pub struct ProfileScoring {
    /// Beta pseudocount for a column being deleted (`glam2scan -D`).
    pub delete_pseudo: f64,
    /// Beta pseudocount for a column being present (`-E`).
    pub no_delete_pseudo: f64,
    /// Beta pseudocount for a residue inserted after a column (`-I`).
    pub insert_pseudo: f64,
    /// Beta pseudocount for no insertion after a column (`-J`).
    pub no_insert_pseudo: f64,
    /// Prior on each column's residue probabilities (`-d`).
    pub prior: DirichletMixture,
}

impl ProfileScoring {
    /// `glam2scan`'s defaults: D = 0.1, E = 2, I = 0.02, J = 1 and the
    /// alphabet's [`DirichletMixture::glam2`] prior.
    pub fn glam2(alphabet: &Alphabet) -> Result<ProfileScoring> {
        let scoring = ProfileScoring {
            delete_pseudo: 0.1,
            no_delete_pseudo: 2.0,
            insert_pseudo: 0.02,
            no_insert_pseudo: 1.0,
            prior: DirichletMixture::glam2(alphabet),
        };
        scoring.check(alphabet.len())?;
        Ok(scoring)
    }

    fn check(&self, alphabet_len: usize) -> Result<()> {
        for (name, v) in [
            ("delete", self.delete_pseudo),
            ("no-delete", self.no_delete_pseudo),
            ("insert", self.insert_pseudo),
            ("no-insert", self.no_insert_pseudo),
        ] {
            if !(v.is_finite() && v > 0.0) {
                return Err(invalid(
                    "profile scoring",
                    format!("{name} pseudocount {v} is not finite and > 0"),
                ));
            }
        }
        if !(self.delete_pseudo + self.no_delete_pseudo).is_finite()
            || !(self.insert_pseudo + self.no_insert_pseudo).is_finite()
        {
            return Err(invalid("profile scoring", "pseudocount sum is not finite"));
        }
        self.prior.check(Some(alphabet_len))
    }
}

/// One profile's alignment scores (natural log).
#[derive(Clone, Debug, PartialEq)]
pub struct ProfileModel {
    width: usize,
    /// Codes scored per column: the alphabet plus the wildcard.
    codes: usize,
    /// Per column, 256 entries so any code indexes without a bound check;
    /// codes past the wildcard are never looked up.
    matches: Vec<[f64; 256]>,
    deletes: Vec<f64>,
    inserts: Vec<f64>,
}

impl ProfileModel {
    /// Scores `view`'s counts against the letter frequencies `background`.
    pub fn build(
        view: &MotifView<'_>,
        background: &[f64],
        scoring: &ProfileScoring,
    ) -> Result<ProfileModel> {
        const WHAT: &str = "gapped profile";
        let n = background.len();
        scoring.check(n)?;
        let gaps = view.gaps().ok_or_else(|| {
            invalid(
                WHAT,
                "the motif has no gap counts; profile scoring needs a gapped profile",
            )
        })?;
        if view.alphabet_len() != n {
            return Err(invalid(
                WHAT,
                format!("{} symbols, alphabet has {n}", view.alphabet_len()),
            ));
        }
        if view.support().is_some() {
            return Err(invalid(
                WHAT,
                "residue counts must cover the whole alphabet",
            ));
        }
        let width = view.width() as usize;
        if view.factor_count() != width {
            return Err(invalid(
                WHAT,
                "one residue-count factor per column is required",
            ));
        }
        let count = |v: f64| -> Result<u32> {
            if v.fract() == 0.0 && (0.0..=u32::MAX as f64).contains(&v) {
                Ok(v as u32)
            } else {
                Err(invalid(WHAT, format!("count {v} is not a whole number")))
            }
        };
        let mut matches = Vec::with_capacity(width);
        let mut deletes = Vec::with_capacity(width);
        let mut inserts = Vec::with_capacity(width);
        let (d_sum, i_sum) = (
            scoring.delete_pseudo + scoring.no_delete_pseudo,
            scoring.insert_pseudo + scoring.no_insert_pseudo,
        );
        for (i, factor) in view.factors().enumerate() {
            let residues = factor
                .table
                .iter()
                .map(|&v| count(v))
                .collect::<Result<Vec<u32>>>()?;
            let m: u64 = residues.iter().map(|&c| c as u64).sum();
            let d = count(gaps.deletions()[i])? as u64;
            let r = match gaps.insertions().get(i) {
                Some(&v) => count(v)? as u64,
                None => 0,
            };
            let aligned = m + d;
            let ds = (d as f64 + scoring.delete_pseudo) / ((d + m) as f64 + d_sum);
            let ms = (m as f64 + scoring.no_delete_pseudo) / ((d + m) as f64 + d_sum);
            let is = (r as f64 + scoring.insert_pseudo) / ((r + aligned) as f64 + i_sum);
            let js = (aligned as f64 + scoring.no_insert_pseudo) / ((r + aligned) as f64 + i_sum);
            let p = scoring.prior.posterior_mean(&residues);
            let mut row = [f64::MIN; 256];
            if i + 1 < width {
                for ((r, &px), &b) in row.iter_mut().zip(&p).zip(background) {
                    *r = (js * ms * px / b).ln();
                }
                deletes.push((ds * js).ln());
                inserts.push(is.ln());
            } else {
                for ((r, &px), &b) in row.iter_mut().zip(&p).zip(background) {
                    *r = (ms * px / b).ln();
                }
                deletes.push(ds.ln());
                inserts.push(f64::MIN);
            }
            matches.push(row);
        }
        Ok(ProfileModel {
            width,
            codes: n + 1,
            matches,
            deletes,
            inserts,
        })
    }

    /// Columns.
    pub fn width(&self) -> usize {
        self.width
    }

    /// Scores of matching each code at `column`; the last entry, for the
    /// wildcard, is `f64::MIN`.
    pub fn match_scores(&self, column: usize) -> &[f64] {
        &self.matches[column][..self.codes]
    }

    /// Per column: the score of deleting it.
    pub fn delete_scores(&self) -> &[f64] {
        &self.deletes
    }

    /// Per column: the score of each residue inserted after it (`f64::MIN`
    /// after the last column).
    pub fn insert_scores(&self) -> &[f64] {
        &self.inserts
    }
}

/// Gapped profiles compiled for scanning under one alphabet, background and
/// strand choice.
#[derive(Clone, Debug)]
pub struct ProfileProgram {
    alphabet: Alphabet,
    strands: Strands,
    models: Vec<ProfileModel>,
}

impl ProfileProgram {
    /// Compiles `motifs`, which must be gapped profiles, against the 0-order
    /// letters of `background`.
    pub fn compile(
        motifs: &[MotifView<'_>],
        alphabet: &Alphabet,
        background: &Background,
        scoring: &ProfileScoring,
        strands: Strands,
    ) -> Result<ProfileProgram> {
        if background.alphabet_len() != alphabet.len() {
            return Err(invalid(
                "motif background",
                format!(
                    "{} symbols, alphabet has {}",
                    background.alphabet_len(),
                    alphabet.len()
                ),
            ));
        }
        if strands.reverse() && alphabet.complement().is_none() {
            return Err(invalid(
                "motif strands",
                format!("alphabet {:?} has no complement", alphabet.name()),
            ));
        }
        if motifs.len() > u32::MAX as usize {
            return Err(invalid("gapped profile", "more than 2^32 − 1 motifs"));
        }
        let models = motifs
            .iter()
            .enumerate()
            .map(|(i, view)| {
                ProfileModel::build(view, background.letters(), scoring).map_err(|e| match e {
                    Error::Invalid { what, detail } => Error::Invalid {
                        what,
                        detail: format!("motif {i}: {detail}"),
                    },
                    other => other,
                })
            })
            .collect::<Result<Vec<_>>>()?;
        Ok(ProfileProgram {
            alphabet: alphabet.clone(),
            strands,
            models,
        })
    }

    pub fn len(&self) -> usize {
        self.models.len()
    }

    pub fn is_empty(&self) -> bool {
        self.models.is_empty()
    }

    pub fn model(&self, motif: usize) -> &ProfileModel {
        &self.models[motif]
    }

    pub fn alphabet(&self) -> &Alphabet {
        &self.alphabet
    }

    pub fn strands(&self) -> Strands {
        self.strands
    }

    /// Scans `sequences` (ASCII bytes, one row per sequence) and returns
    /// what `reducer` keeps. Sequences are aligned in parallel and handed to
    /// the reducer in order, so the result does not depend on the thread
    /// count.
    pub fn scan<R: AlignmentReducer>(
        &self,
        sequences: Strings<'_>,
        mut reducer: R,
        options: &ScanOptions,
    ) -> Result<R::Output> {
        check_sequences(&sequences, self.models.iter().map(|m| m.width).max())?;
        reducer.start(self, sequences)?;
        let limit = reducer.per_strand_limit();
        let count = sequences.len();
        // Per motif, the reducer's latest floor. Floors never decrease, so a
        // worker reading a stale one only stops later than it could.
        let floors: Vec<AtomicU64> = (0..self.len())
            .map(|_| AtomicU64::new(f64::NEG_INFINITY.to_bits()))
            .collect();
        let floors = &floors;
        let raise = |reducer: &R| {
            for (m, f) in floors.iter().enumerate() {
                let x = reducer.floor(m as u32);
                if x > f64::from_bits(f.load(AtomicOrdering::Relaxed)) {
                    f.store(x.to_bits(), AtomicOrdering::Relaxed);
                }
            }
        };
        let threads = options
            .threads
            .unwrap_or_else(|| std::thread::available_parallelism().map_or(1, |n| n.get()))
            .clamp(1, count.max(1));
        if threads == 1 {
            let mut aligner = Aligner::default();
            for s in 0..count {
                let streams = aligner.sequence(self, sequences, s, limit, floors);
                reducer.sequence(self, s as u32, streams);
                raise(&reducer);
            }
            return Ok(reducer.finish(self));
        }
        let next = &AtomicUsize::new(0);
        // Sequences handed to the reducer; workers stay within a window
        // ahead of it so buffered results stay bounded.
        let gate = &(Mutex::new(0usize), Condvar::new());
        let ahead = threads * 4;
        std::thread::scope(|scope| {
            let (tx, rx) = mpsc::channel::<(usize, Vec<AlignmentStream>)>();
            for _ in 0..threads {
                let tx = tx.clone();
                scope.spawn(move || {
                    let mut aligner = Aligner::default();
                    loop {
                        let s = next.fetch_add(1, AtomicOrdering::Relaxed);
                        if s >= count {
                            break;
                        }
                        let mut done = gate.0.lock().expect("scan gate");
                        while s >= *done + ahead {
                            done = gate.1.wait(done).expect("scan gate");
                        }
                        drop(done);
                        let streams = aligner.sequence(self, sequences, s, limit, floors);
                        if tx.send((s, streams)).is_err() {
                            break;
                        }
                    }
                });
            }
            drop(tx);
            let mut pending = BTreeMap::new();
            let mut want = 0usize;
            for (s, streams) in rx {
                pending.insert(s, streams);
                while let Some(streams) = pending.remove(&want) {
                    reducer.sequence(self, want as u32, streams);
                    raise(&reducer);
                    want += 1;
                    *gate.0.lock().expect("scan gate") = want;
                    gate.1.notify_all();
                }
            }
        });
        Ok(reducer.finish(self))
    }

    /// Scans one sequence.
    pub fn scan_sequence<R: AlignmentReducer>(
        &self,
        sequence: &[u8],
        reducer: R,
        options: &ScanOptions,
    ) -> Result<R::Output> {
        let offsets = [0i64, sequence.len() as i64];
        self.scan(
            Strings {
                offsets: &offsets,
                data: sequence,
            },
            reducer,
            options,
        )
    }
}

fn check_sequences(sequences: &Strings<'_>, max_width: Option<usize>) -> Result<()> {
    let bad = |detail: String| Err(invalid("motif scan sequences", detail));
    if sequences.offsets.is_empty() {
        return bad("offsets are empty".into());
    }
    if sequences.len() > u32::MAX as usize {
        return bad("more than 2^32 − 1 sequences".into());
    }
    let mut prev = sequences.offsets[0];
    if prev < 0 {
        return bad("negative offset".into());
    }
    let mut longest = 0i64;
    for &o in &sequences.offsets[1..] {
        if o < prev {
            return bad("offsets decrease".into());
        }
        if o - prev >= u32::MAX as i64 {
            return bad("a sequence is longer than 2^32 − 2".into());
        }
        longest = longest.max(o - prev);
        prev = o;
    }
    if prev as u64 > sequences.data.len() as u64 {
        return bad("offsets run past the data".into());
    }
    let cells = (longest as usize + 1).checked_mul(max_width.unwrap_or(0) + 1);
    if cells.is_none_or(|c| c > isize::MAX as usize / 8) {
        return bad("the alignment matrix of the longest sequence is too large".into());
    }
    Ok(())
}

/// One aligned column of a hit.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct Column {
    /// Whether a residue fills the column (otherwise it is deleted).
    pub matched: bool,
    /// Along the scanned strand, from the hit's first residue: the matched
    /// residue, or for a deletion the residue it precedes.
    pub offset: u32,
}

/// One local alignment of a profile to a sequence strand.
#[derive(Clone, Debug, PartialEq)]
pub struct Alignment {
    pub sequence: u32,
    pub motif: u32,
    pub strand: Strand,
    /// Forward-strand half-open span of the aligned residues.
    pub start: u32,
    pub end: u32,
    /// Natural log; [`Self::bits`] in bits.
    pub score: f64,
    /// One per profile column, in column order. Residues of the span not
    /// covered by a matched column are insertions.
    pub columns: Vec<Column>,
}

impl Alignment {
    pub fn bits(&self) -> f64 {
        self.score / 2f64.ln()
    }

    /// Columns filled by a residue.
    pub fn matches(&self) -> usize {
        self.columns.iter().filter(|c| c.matched).count()
    }
}

/// The hits of one motif on one strand of a sequence, in the order the
/// search finds them (scores non-increasing).
#[derive(Clone, Debug, PartialEq)]
pub struct AlignmentStream {
    pub motif: u32,
    pub strand: Strand,
    pub alignments: Vec<Alignment>,
}

/// Consumes hit streams sequence by sequence.
pub trait AlignmentReducer {
    type Output;

    /// Called once before the scan.
    fn start(&mut self, program: &ProfileProgram, sequences: Strings<'_>) -> Result<()> {
        let _ = (program, sequences);
        Ok(())
    }

    /// The most hits each strand stream needs; the search stops there.
    fn per_strand_limit(&self) -> Option<usize>;

    /// A score at or below which no later hit of `motif` is needed; a
    /// stream ends before its first such hit. It must never decrease.
    fn floor(&self, motif: u32) -> f64 {
        let _ = motif;
        f64::NEG_INFINITY
    }

    /// One sequence's streams, motif by motif, forward before reverse.
    /// Sequences arrive in order.
    fn sequence(&mut self, program: &ProfileProgram, sequence: u32, streams: Vec<AlignmentStream>);

    fn finish(self, program: &ProfileProgram) -> Self::Output;
}

/// Every hit of every stream (up to `limit` per strand), in scan order.
#[derive(Clone, Debug, Default)]
pub struct AllAlignments {
    limit: Option<usize>,
    out: Vec<Alignment>,
}

impl AllAlignments {
    pub fn new(per_strand_limit: Option<usize>) -> AllAlignments {
        AllAlignments {
            limit: per_strand_limit,
            out: Vec::new(),
        }
    }
}

impl AlignmentReducer for AllAlignments {
    type Output = Vec<Alignment>;

    fn per_strand_limit(&self) -> Option<usize> {
        self.limit
    }

    fn sequence(&mut self, _: &ProfileProgram, _: u32, streams: Vec<AlignmentStream>) {
        for s in streams {
            self.out.extend(s.alignments);
        }
    }

    fn finish(self, _: &ProfileProgram) -> Vec<Alignment> {
        self.out
    }
}

/// The best `n` hits per motif as `glam2scan` keeps them: one bounded heap
/// fed stream by stream, each stream abandoned at its first hit that does
/// not beat the heap's worst; the survivors sorted by score, ties by
/// sequence name (byte order), remaining ties in heap order.
#[derive(Clone, Debug)]
pub struct TopAlignments<'a> {
    n: usize,
    names: Strings<'a>,
    heaps: Vec<Vec<Alignment>>,
}

impl<'a> TopAlignments<'a> {
    /// `names[i]` names sequence `i` (GLAM2 uses the first word of the
    /// FASTA title).
    pub fn glam2(n: usize, names: Strings<'a>) -> TopAlignments<'a> {
        TopAlignments {
            n,
            names,
            heaps: Vec::new(),
        }
    }

    /// GLAM2's rank: higher score first, then smaller name.
    fn rank(&self, a: &Alignment, b: &Alignment) -> Ordering {
        if a.score < b.score {
            Ordering::Greater
        } else if a.score > b.score {
            Ordering::Less
        } else {
            self.names
                .get(a.sequence as usize)
                .cmp(self.names.get(b.sequence as usize))
        }
    }

    /// Restores the heap (worst-ranked at the root) after appending at
    /// `n − 1`; positions are one-based.
    fn sift_up(&self, heap: &mut [Alignment], mut n: usize) {
        while n > 1 {
            let parent = n / 2;
            if self.rank(&heap[n - 1], &heap[parent - 1]) != Ordering::Greater {
                break;
            }
            heap.swap(n - 1, parent - 1);
            n = parent;
        }
    }

    /// Drops the root: the hole it leaves descends along the worse child to
    /// a leaf, which receives the last entry.
    fn pop_root(&self, heap: &mut Vec<Alignment>) {
        let n = heap.len();
        heap.swap(0, n - 1);
        let mut hole = 1;
        loop {
            let mut child = hole * 2;
            if child >= n {
                break;
            }
            if child + 1 < n && self.rank(&heap[child - 1], &heap[child]) == Ordering::Less {
                child += 1;
            }
            heap.swap(hole - 1, child - 1);
            hole = child;
        }
        self.sift_up(heap, hole);
        heap.pop();
    }
}

impl AlignmentReducer for TopAlignments<'_> {
    type Output = Vec<Vec<Alignment>>;

    fn start(&mut self, program: &ProfileProgram, sequences: Strings<'_>) -> Result<()> {
        if self.names.len() < sequences.len() {
            return Err(invalid(
                "motif scan sequences",
                format!(
                    "{} names for {} sequences",
                    self.names.len(),
                    sequences.len()
                ),
            ));
        }
        self.heaps = vec![Vec::new(); program.len()];
        Ok(())
    }

    // Truncating each stream at `n` hits loses nothing: a stream's scores
    // never increase (forbidding cells only lowers the matrix), so a later
    // hit of the stream never evicts an earlier one, and after `n` of its
    // hits the heap holds exactly those, none worse than hit `n + 1`, which
    // is therefore rejected and ends the stream.
    fn per_strand_limit(&self) -> Option<usize> {
        Some(self.n)
    }

    // A full heap rejects any hit not above its worst, whose score only
    // rises as hits are admitted.
    fn floor(&self, motif: u32) -> f64 {
        match self.heaps.get(motif as usize) {
            Some(heap) if self.n > 0 && heap.len() == self.n => heap[0].score,
            _ => f64::NEG_INFINITY,
        }
    }

    fn sequence(&mut self, _: &ProfileProgram, _: u32, streams: Vec<AlignmentStream>) {
        let mut heaps = std::mem::take(&mut self.heaps);
        for stream in streams {
            let heap = &mut heaps[stream.motif as usize];
            for hit in stream.alignments {
                if heap.len() == self.n {
                    if self.n == 0 || hit.score <= heap[0].score {
                        break;
                    }
                    self.pop_root(heap);
                }
                heap.push(hit);
                let len = heap.len();
                self.sift_up(heap, len);
            }
        }
        self.heaps = heaps;
    }

    fn finish(mut self, _: &ProfileProgram) -> Vec<Vec<Alignment>> {
        let mut heaps = std::mem::take(&mut self.heaps);
        for heap in &mut heaps {
            heap.sort_by(|a, b| self.rank(a, b));
        }
        heaps
    }
}

/// Reusable buffers for aligning one strand at a time.
#[derive(Default)]
struct Aligner {
    codes: Vec<u8>,
    reverse: Vec<u8>,
    /// `(width + 1) × (len + 1)` best scores, row `i` after `i` columns.
    dp: Vec<f64>,
    /// `width × len`: cells whose match an earlier hit already used.
    forbidden: Vec<bool>,
    /// Per column: matched, and its sequence position on the strand.
    trace: Vec<(bool, usize)>,
}

impl Aligner {
    fn sequence(
        &mut self,
        program: &ProfileProgram,
        sequences: Strings<'_>,
        s: usize,
        limit: Option<usize>,
        floors: &[AtomicU64],
    ) -> Vec<AlignmentStream> {
        let n = program.alphabet.len();
        let lut = program.alphabet.lut();
        let bytes = sequences.get(s);
        let wildcard = n as u8;
        self.codes.clear();
        self.codes
            .extend(bytes.iter().map(|&b| match lut[b as usize] {
                INVALID => wildcard,
                c => c,
            }));
        if program.strands.reverse() {
            let comp = program.alphabet.complement().expect("checked at compile");
            self.reverse.clear();
            self.reverse.extend(self.codes.iter().rev().map(|&c| {
                if c == wildcard {
                    wildcard
                } else {
                    comp[c as usize]
                }
            }));
        }
        let len = bytes.len() as u32;
        let mut out = Vec::new();
        for (mi, model) in program.models.iter().enumerate() {
            for strand in [Strand::Forward, Strand::Reverse] {
                let on = match strand {
                    Strand::Forward => program.strands.forward(),
                    _ => program.strands.reverse(),
                };
                if !on {
                    continue;
                }
                let codes = std::mem::take(match strand {
                    Strand::Forward => &mut self.codes,
                    _ => &mut self.reverse,
                });
                let mut alignments = Vec::new();
                let floor = f64::from_bits(floors[mi].load(AtomicOrdering::Relaxed));
                self.align(model, &codes, limit, floor, |score, start, end, columns| {
                    let (start, end) = match strand {
                        Strand::Forward => (start, end),
                        _ => (len - end, len - start),
                    };
                    alignments.push(Alignment {
                        sequence: s as u32,
                        motif: mi as u32,
                        strand,
                        start,
                        end,
                        score,
                        columns,
                    });
                });
                match strand {
                    Strand::Forward => self.codes = codes,
                    _ => self.reverse = codes,
                }
                out.push(AlignmentStream {
                    motif: mi as u32,
                    strand,
                    alignments,
                });
            }
        }
        out
    }

    /// The hit stream of `model` on one strand's codes: `emit(score, start,
    /// end, columns)` with strand coordinates, at most `limit` times and
    /// only while the score stays above `floor`.
    fn align(
        &mut self,
        model: &ProfileModel,
        seq: &[u8],
        limit: Option<usize>,
        floor: f64,
        mut emit: impl FnMut(f64, u32, u32, Vec<Column>),
    ) {
        let limit = limit.unwrap_or(usize::MAX);
        if limit == 0 {
            return;
        }
        let w = model.width;
        let l = seq.len();
        let cols = l + 1;
        let dp = &mut self.dp;
        dp.clear();
        dp.resize((w + 1) * cols, 0.0);
        for i in 0..w {
            dp[(i + 1) * cols] = model.deletes[i] + dp[i * cols];
        }
        for i in 0..w {
            let ms = &model.matches[i];
            let (ds, is) = (model.deletes[i], model.inserts[i]);
            let (above, below) = dp.split_at_mut((i + 1) * cols);
            let prev = &above[i * cols..];
            let row = &mut below[..cols];
            let (mut left, mut diag) = (row[0], prev[0]);
            for ((&c, &up), out) in seq.iter().zip(&prev[1..]).zip(&mut row[1..]) {
                let m = ms[c as usize] + diag;
                let d = ds + up;
                let mut v = left + is;
                if m > v {
                    v = m;
                }
                if d > v {
                    v = d;
                }
                *out = v;
                left = v;
                diag = up;
            }
        }
        let forbidden = &mut self.forbidden;
        forbidden.clear();
        forbidden.resize(w * l, false);
        let trace = &mut self.trace;
        trace.clear();
        trace.resize(w, (false, 0));
        let cell = |i: usize, j: usize| i * cols + j;
        let mut emitted = 0usize;
        loop {
            let last = &dp[cell(w, 0)..cell(w + 1, 0)];
            let mut end = 0usize;
            for k in 1..=l {
                if last[k] > last[end] {
                    end = k;
                }
            }
            let score = last[end];
            if score <= floor {
                return;
            }
            // Trace back from the best end: match, deletion, insertion on ties.
            let (mut i, mut j) = (w, end);
            let mut matched = 0usize;
            while i > 0 && j > 0 {
                let (ci, sj) = (i - 1, j - 1);
                let m = if forbidden[ci * l + sj] {
                    f64::MIN
                } else {
                    model.matches[ci][seq[sj] as usize] + dp[cell(ci, sj)]
                };
                let d = model.deletes[ci] + dp[cell(ci, j)];
                let ins = model.inserts[ci] + dp[cell(i, sj)];
                if m >= d && m >= ins {
                    trace[ci] = (true, sj);
                    forbidden[ci * l + sj] = true;
                    matched += 1;
                    i -= 1;
                    j -= 1;
                } else if d >= ins {
                    trace[ci] = (false, j);
                    i -= 1;
                } else {
                    j -= 1;
                }
            }
            for t in &mut trace[..i] {
                *t = (false, 0);
            }
            let start = j;
            if matched == 0 {
                return;
            }
            let columns = trace
                .iter()
                .map(|&(matched, pos)| Column {
                    matched,
                    offset: (pos - start) as u32,
                })
                .collect();
            emit(score, start as u32, end as u32, columns);
            emitted += 1;
            if emitted == limit {
                return;
            }
            // Recompute the cells the forbidden matches can lower, row by
            // row from the hit's start; a row stops at the first unchanged
            // cell past every change seen so far.
            let mut reach = start;
            for i in 0..w {
                let ms = &model.matches[i];
                let (ds, is) = (model.deletes[i], model.inserts[i]);
                for j in start..l {
                    let m = if forbidden[i * l + j] {
                        f64::MIN
                    } else {
                        ms[seq[j] as usize] + dp[cell(i, j)]
                    };
                    let d = ds + dp[cell(i, j + 1)];
                    let ins = is + dp[cell(i + 1, j)];
                    let mut v = m;
                    if d > v {
                        v = d;
                    }
                    if ins > v {
                        v = ins;
                    }
                    if v != dp[cell(i + 1, j + 1)] {
                        dp[cell(i + 1, j + 1)] = v;
                        reach = reach.max(j);
                    } else if j > reach {
                        break;
                    }
                }
            }
        }
    }
}

/// How hits print in GLAM2's text: sequence names and bytes, and the
/// character printed for each code.
pub struct Glam2Text<'a> {
    /// The scanned alphabet; reverse-strand residues print complemented.
    pub alphabet: &'a Alphabet,
    /// Per code, the character printed.
    pub symbols: Vec<u8>,
    /// Printed for an inserted byte outside the alphabet.
    pub wildcard: u8,
    pub sequence_names: Strings<'a>,
    pub sequences: Strings<'a>,
}

impl<'a> Glam2Text<'a> {
    /// GLAM2's characters: lowercase nucleotides with wildcard `N`;
    /// uppercase amino acids (and other alphabets' symbols) with wildcard
    /// `x`.
    pub fn new(
        alphabet: &'a Alphabet,
        sequence_names: Strings<'a>,
        sequences: Strings<'a>,
    ) -> Glam2Text<'a> {
        let (symbols, wildcard) = match kind(alphabet) {
            Kind::Nucleotide => (alphabet.symbols().to_ascii_lowercase(), b'N'),
            _ => (alphabet.symbols().to_vec(), b'x'),
        };
        Glam2Text {
            alphabet,
            symbols,
            wildcard,
            sequence_names,
            sequences,
        }
    }
}

/// `glam2scan`'s header: program, version and command line.
pub fn write_glam2_header<W: Write>(
    mut out: W,
    version: &str,
    command_line: &str,
) -> io::Result<()> {
    write!(out, "GLAM2SCAN\nVersion {version}\n\n{command_line}\n\n")
}

/// Writes `hits` as `glam2scan -t` prints them: a line marking columns
/// (`*`) and insertions (`.`), then the name (24 columns), the first
/// coordinate, the aligned residues (`.` for deletions), the last
/// coordinate (1-based, reversed on the `-` strand), the strand and the
/// score in bits (`%#.3g`), then an empty line.
pub fn write_glam2_text<W: Write>(
    mut out: W,
    hits: &[Alignment],
    text: &Glam2Text<'_>,
) -> io::Result<()> {
    let lut = text.alphabet.lut();
    let complement = text.alphabet.complement();
    let mut marks = Vec::new();
    let mut residues = Vec::new();
    for h in hits {
        let bytes = text.sequences.get(h.sequence as usize);
        let reverse = h.strand == Strand::Reverse;
        let residue = |k: u32| -> u8 {
            let b = if reverse {
                bytes[(h.end - 1 - k) as usize]
            } else {
                bytes[(h.start + k) as usize]
            };
            match (lut[b as usize], complement) {
                (INVALID, _) => text.wildcard,
                (c, Some(comp)) if reverse => text.symbols[comp[c as usize] as usize],
                (c, _) => text.symbols[c as usize],
            }
        };
        marks.clear();
        residues.clear();
        let mut at = 0u32;
        for c in &h.columns {
            while at < c.offset {
                marks.push(b'.');
                residues.push(residue(at));
                at += 1;
            }
            marks.push(b'*');
            if c.matched {
                residues.push(residue(at));
                at += 1;
            } else {
                residues.push(b'.');
            }
        }
        let (first, last) = if reverse {
            (h.end, h.start + 1)
        } else {
            (h.start + 1, h.end)
        };
        let first_text = first.to_string();
        let indent = GLAM2_NAME_WIDTH + first_text.len() + 2;
        out.write_all(&vec![b' '; indent])?;
        out.write_all(&marks)?;
        out.write_all(b"\n")?;
        let name = text.sequence_names.get(h.sequence as usize);
        let name = &name[..name.len().min(GLAM2_NAME_WIDTH)];
        out.write_all(name)?;
        out.write_all(&vec![b' '; GLAM2_NAME_WIDTH - name.len()])?;
        write!(out, " {first_text} ")?;
        out.write_all(&residues)?;
        writeln!(
            out,
            " {last} {} {}\n",
            if reverse { '-' } else { '+' },
            format_alt_g(h.bits(), 3)
        )?;
    }
    Ok(())
}

/// C's `%#.{precision}g`: `%g` keeping trailing zeros and the decimal point.
pub fn format_alt_g(x: f64, precision: usize) -> String {
    let p = precision.max(1);
    if x.is_nan() {
        return if x.is_sign_negative() { "-nan" } else { "nan" }.into();
    }
    if x.is_infinite() {
        return if x < 0.0 { "-inf" } else { "inf" }.into();
    }
    let sci = format!("{:.*e}", p - 1, x);
    let (mantissa, exp) = sci.split_once('e').expect("exponent");
    let exp: i32 = exp.parse().expect("exponent digits");
    let mut s = if exp < -4 || exp >= p as i32 {
        let mut m = mantissa.to_string();
        if !m.contains('.') {
            m.push('.');
        }
        let sign = if exp < 0 { '-' } else { '+' };
        format!("{m}e{sign}{:02}", exp.abs())
    } else {
        format!("{:.*}", (p as i32 - 1 - exp) as usize, x)
    };
    if !s.contains('.') {
        s.push('.');
    }
    s
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn rising_products_agree_with_direct_summation() {
        for a in [0.01, 0.4, 1.0, 1.6, 20.0, 1e6, 1e100] {
            for k in [0, 10, 256, 257, 4096, 65_536] {
                let expected: f64 = (0..k).map(|i| (a + i as f64).ln()).sum();
                let actual = log_rising(a, k);
                assert!(
                    (actual - expected).abs() <= expected.abs().max(1.0) * 8e-12,
                    "a={a}, k={k}: {actual} != {expected}"
                );
            }
        }
    }

    #[test]
    fn large_counts_have_finite_normalized_posterior_means() {
        let counts = [u32::MAX, u32::MAX, 0, 0];
        let single = DirichletMixture::glam2_dna().posterior_mean(&counts);
        let mixture =
            DirichletMixture::new(vec![0.25, 0.75], vec![vec![0.4; 4], vec![2.0; 4]]).unwrap();
        for mean in [single, mixture.posterior_mean(&counts)] {
            assert!(mean.iter().all(|p| p.is_finite() && *p > 0.0));
            assert!((mean.iter().sum::<f64>() - 1.0).abs() < 1e-12);
            assert_eq!(mean[0], mean[1]);
            assert_eq!(mean[2], mean[3]);
        }
        assert!(DirichletMixture::new(vec![1.0], vec![vec![f64::MAX; 4]]).is_err());
        let tiny = DirichletMixture::new(vec![1.0], vec![vec![1e-320; 2]]).unwrap();
        assert_eq!(tiny.posterior_mean(&[0, 0]), vec![0.5, 0.5]);
    }
}
