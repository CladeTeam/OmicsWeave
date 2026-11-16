//! Markov background models: the null hypothesis a window is scored and
//! tested against.

use omicsweave_core::{Error, Result, Strings};
use omicsweave_ffm::{Alphabet, INVALID};

/// NRDB letter frequencies MEME ships for DNA (`ACGT`).
const NRDB_DNA: [f64; 4] = [0.281774, 0.222020, 0.228876, 0.267330];

/// NRDB letter frequencies MEME ships for protein (`ACDEFGHIKLMNPQRSTVWY`).
const NRDB_PROTEIN: [f64; 20] = [
    0.073164, 0.018163, 0.051739, 0.062340, 0.040283, 0.069328, 0.022428, 0.056282, 0.058493,
    0.091712, 0.023067, 0.046077, 0.050674, 0.040755, 0.051897, 0.073802, 0.059411, 0.064362,
    0.013341, 0.032682,
];

const PROTEIN_SYMBOLS: &[u8] = b"ACDEFGHIKLMNPQRSTVWY";

/// Pseudocount MEME adds to every background letter before normalizing.
const MEME_BACKGROUND_EPSILON: f64 = 0.0000005;

/// An order-`m` Markov model over an alphabet's global codes.
///
/// It is defined by positive k-mer frequencies for `k = 1 ..= m + 1`, kept
/// as given. The probability of a symbol after a context of `c ≤ m` symbols
/// is the (c + 1)-mer frequency normalized over the last symbol. A sequence
/// is generated with the longest context available, so after its start (or
/// an unknown symbol) the chain restarts with lower orders.
#[derive(Clone, Debug, PartialEq)]
pub struct Background {
    alphabet_len: usize,
    /// `kmers[k - 1]`: `n^k` frequencies, first symbol most significant.
    kmers: Vec<Vec<f64>>,
    /// `conditional[c]`: `n^(c + 1)` probabilities, index `context * n + x`.
    conditional: Vec<Vec<f64>>,
}

impl Background {
    /// A model from positive k-mer frequencies, `kmers[k - 1]` holding the
    /// `n^k` frequencies of k-mers (first symbol most significant) for
    /// `k = 1 ..= order + 1`. Frequencies need not be normalized.
    pub fn from_kmer_frequencies(alphabet_len: usize, kmers: Vec<Vec<f64>>) -> Result<Background> {
        if alphabet_len == 0 || alphabet_len > 255 {
            return Err(invalid(format!("alphabet of {alphabet_len} symbols")));
        }
        if kmers.is_empty() {
            return Err(invalid("no frequencies"));
        }
        let mut need = 1usize;
        for (k, f) in kmers.iter().enumerate() {
            need = need
                .checked_mul(alphabet_len)
                .filter(|&n| n <= 1 << 24)
                .ok_or_else(|| invalid(format!("order {k} is too large")))?;
            if f.len() != need {
                return Err(invalid(format!(
                    "{}-mers need {need} frequencies, got {}",
                    k + 1,
                    f.len()
                )));
            }
            if let Some(v) = f.iter().find(|v| !(v.is_finite() && **v > 0.0)) {
                return Err(invalid(format!(
                    "{}-mer frequency {v} is not finite and > 0",
                    k + 1
                )));
            }
        }
        let conditional = kmers
            .iter()
            .map(|f| {
                let mut c = f.clone();
                for group in c.chunks_mut(alphabet_len) {
                    let sum: f64 = group.iter().sum();
                    if sum.is_finite() {
                        group.iter_mut().for_each(|v| *v /= sum);
                    } else {
                        let max = group.iter().copied().fold(0.0, f64::max);
                        let scaled: f64 = group.iter().map(|v| v / max).sum();
                        group.iter_mut().for_each(|v| *v = (*v / max) / scaled);
                    }
                }
                c
            })
            .collect();
        Ok(Background {
            alphabet_len,
            kmers,
            conditional,
        })
    }

    /// A 0-order model with the given (positive) letter frequencies.
    pub fn frequencies(frequencies: &[f64]) -> Result<Background> {
        Background::from_kmer_frequencies(frequencies.len(), vec![frequencies.to_vec()])
    }

    /// The 0-order uniform model.
    pub fn uniform(alphabet: &Alphabet) -> Background {
        let n = alphabet.len();
        Background::frequencies(&vec![1.0 / n as f64; n]).expect("valid uniform background")
    }

    /// The 0-order NRDB frequencies MEME uses for DNA (`ACGT`) and protein
    /// (the 20 amino acids in alphabetical order), normalized as MEME does.
    pub fn nrdb(alphabet: &Alphabet) -> Result<Background> {
        let mut f = match alphabet.symbols() {
            b"ACGT" => NRDB_DNA.to_vec(),
            s if s == PROTEIN_SYMBOLS => NRDB_PROTEIN.to_vec(),
            s => {
                return Err(invalid(format!(
                    "NRDB frequencies exist for ACGT and the 20 amino acids, not {:?}",
                    String::from_utf8_lossy(s)
                )));
            }
        };
        meme_normalize(&mut f, 0.0);
        Background::frequencies(&f)
    }

    /// A model from the text of a MEME background file: lines of a symbol
    /// chain and its frequency, `#` comments, chains of each length listed
    /// completely before longer ones. Symbols map through the alphabet's
    /// lookup table.
    pub fn from_meme(alphabet: &Alphabet, text: &str) -> Result<Background> {
        let n = alphabet.len();
        let mut kmers: Vec<Vec<f64>> = Vec::new();
        let mut seen: Vec<Vec<bool>> = Vec::new();
        for (line_no, line) in text.lines().enumerate() {
            let line = line.split('#').next().unwrap_or("");
            let mut fields = line.split_whitespace();
            let Some(chain) = fields.next() else { continue };
            let line_err = |detail: String| Error::Format {
                line: line_no as u64 + 1,
                detail,
            };
            let value: f64 = fields
                .next()
                .ok_or_else(|| line_err(format!("no frequency for {chain:?}")))?
                .parse()
                .map_err(|_| line_err(format!("bad frequency for {chain:?}")))?;
            if fields.next().is_some() {
                return Err(line_err(format!("junk after {chain:?}")));
            }
            let k = chain.len();
            if k > kmers.len() + 1 || k < kmers.len() {
                return Err(line_err(format!(
                    "{chain:?} is out of order: chains of each length must follow the shorter ones"
                )));
            }
            if k == kmers.len() + 1 {
                if let Some(prev) = seen.last()
                    && prev.iter().any(|s| !s)
                {
                    return Err(line_err(format!("{}-mers are incomplete", k - 1)));
                }
                let size = n
                    .checked_pow(k as u32)
                    .filter(|&s| s <= 1 << 24)
                    .ok_or_else(|| line_err(format!("{chain:?} is too long")))?;
                kmers.push(vec![0.0; size]);
                seen.push(vec![false; size]);
            }
            let mut index = 0usize;
            for &b in chain.as_bytes() {
                let code = alphabet.lut()[b as usize];
                if code == INVALID {
                    return Err(line_err(format!("{:?} is not in the alphabet", b as char)));
                }
                index = index * n + code as usize;
            }
            if seen[k - 1][index] {
                return Err(line_err(format!("{chain:?} is listed twice")));
            }
            seen[k - 1][index] = true;
            kmers[k - 1][index] = value;
        }
        if let Some(last) = seen.last()
            && last.iter().any(|s| !s)
        {
            return Err(invalid(format!("{}-mers are incomplete", seen.len())));
        }
        Background::from_kmer_frequencies(n, kmers)
    }

    /// Estimates an order-`order` model from sequences: k-mers over runs of
    /// alphabet symbols, each count plus `pseudocount` (> 0).
    pub fn estimate(
        alphabet: &Alphabet,
        sequences: Strings<'_>,
        order: usize,
        pseudocount: f64,
    ) -> Result<Background> {
        if !(pseudocount.is_finite() && pseudocount > 0.0) {
            return Err(invalid(format!("pseudocount {pseudocount} is not > 0")));
        }
        crate::scan::check_sequences(&sequences)?;
        if order > 32 {
            return Err(invalid("order above 32"));
        }
        let n = alphabet.len();
        let lut = alphabet.lut();
        let mut kmers: Vec<Vec<f64>> = Vec::with_capacity(order + 1);
        let mut size = 1usize;
        for _ in 0..=order {
            size = size
                .checked_mul(n)
                .filter(|&s| s <= 1 << 24)
                .ok_or_else(|| invalid(format!("order {order} is too large")))?;
            kmers.push(vec![pseudocount; size]);
        }
        for i in 0..sequences.len() {
            let mut run = 0usize;
            let mut window = 0usize;
            for &b in sequences.get(i) {
                let c = lut[b as usize];
                if c == INVALID {
                    run = 0;
                    window = 0;
                    continue;
                }
                window = (window * n + c as usize) % kmers[order].len();
                run += 1;
                for table in kmers.iter_mut().take(run.min(order + 1)) {
                    let size = table.len();
                    table[window % size] += 1.0;
                }
            }
        }
        Background::from_kmer_frequencies(n, kmers)
    }

    /// Averages every k-mer frequency with its reverse complement's.
    pub fn symmetrized(&self, complement: &[u8]) -> Result<Background> {
        let n = self.alphabet_len;
        if complement.len() != n {
            return Err(invalid("complement does not match the alphabet"));
        }
        if complement
            .iter()
            .enumerate()
            .any(|(i, &c)| complement.get(c as usize).copied() != Some(i as u8))
        {
            return Err(invalid(
                "complement must be an involution over the alphabet",
            ));
        }
        let kmers = self
            .kmers
            .iter()
            .enumerate()
            .map(|(k1, f)| {
                let k = k1 + 1;
                let mut out = f.clone();
                let mut digits = vec![0usize; k];
                for (i, o) in out.iter_mut().enumerate() {
                    decode(i, n, &mut digits);
                    let rc = digits
                        .iter()
                        .rev()
                        .fold(0usize, |acc, &d| acc * n + complement[d] as usize);
                    *o = finite_average(f[i], f[rc]);
                }
                out
            })
            .collect();
        Background::from_kmer_frequencies(n, kmers)
    }

    /// The same model cut to order `order` (no change if already lower).
    pub fn truncated(&self, order: usize) -> Background {
        let keep = order.saturating_add(1).min(self.kmers.len());
        Background {
            alphabet_len: self.alphabet_len,
            kmers: self.kmers[..keep].to_vec(),
            conditional: self.conditional[..keep].to_vec(),
        }
    }

    /// The 0-order letter frequencies after MEME's background preparation:
    /// averaged with the complement when `complement` is given, then each
    /// raised by 5e-7 and normalized.
    pub(crate) fn meme_prepared(&self, complement: Option<&[u8]>) -> Background {
        let mut f = self.kmers[0].clone();
        if let Some(comp) = complement {
            for i in 0..f.len() {
                let c = comp[i] as usize;
                if comp[c] as usize == i && i < c {
                    let avg = finite_average(f[i], f[c]);
                    f[i] = avg;
                    f[c] = avg;
                }
            }
        }
        let mut sum = MEME_BACKGROUND_EPSILON * f.len() as f64;
        for v in &f {
            sum += v;
        }
        if sum.is_finite() {
            for v in &mut f {
                *v = (*v + MEME_BACKGROUND_EPSILON) / sum;
            }
        } else {
            let max = f.iter().copied().fold(0.0, f64::max);
            let epsilon = MEME_BACKGROUND_EPSILON / max;
            let scaled: f64 = f.iter().map(|v| v / max + epsilon).sum();
            for v in &mut f {
                *v = (*v / max + epsilon) / scaled;
            }
        }
        meme_normalize(&mut f, 0.0);
        // MEME uses these letters as they are, without renormalizing.
        Background {
            alphabet_len: f.len(),
            kmers: vec![f.clone()],
            conditional: vec![f],
        }
    }

    /// Generates sequences of the given lengths from the chain (restarted at
    /// each sequence), as the canonical symbols of `alphabet`, in the
    /// `(offsets, data)` layout of [`Strings`]. Deterministic for a seed.
    ///
    /// # Panics
    ///
    /// Panics if the alphabet size differs from the background, or the
    /// combined lengths cannot fit the output's address space and i64 offsets.
    pub fn sample(&self, alphabet: &Alphabet, lengths: &[usize], seed: u64) -> (Vec<i64>, Vec<u8>) {
        let n = self.alphabet_len;
        let m = self.order();
        assert_eq!(
            alphabet.len(),
            n,
            "sample alphabet and background sizes differ"
        );
        let total = lengths
            .iter()
            .try_fold(0usize, |sum, &len| sum.checked_add(len))
            .filter(|&len| len <= isize::MAX as usize && len as u128 <= i64::MAX as u128)
            .expect("sample lengths exceed the output address space or i64 offsets");
        let symbols = alphabet.symbols();
        let mut state = seed ^ 0x9E37_79B9_7F4A_7C15;
        let mut next = move || {
            // splitmix64
            state = state.wrapping_add(0x9E37_79B9_7F4A_7C15);
            let mut z = state;
            z = (z ^ (z >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
            z = (z ^ (z >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
            ((z ^ (z >> 31)) >> 11) as f64 / (1u64 << 53) as f64
        };
        let mut offsets = Vec::with_capacity(lengths.len() + 1);
        let mut data = Vec::with_capacity(total);
        offsets.push(0i64);
        let context_states = self.conditional[m].len() / n;
        for &len in lengths {
            let mut ctx = 0usize;
            for i in 0..len {
                let c = m.min(i);
                let probs = &self.conditional[c][ctx * n..(ctx + 1) * n];
                let u = next();
                let mut acc = 0.0;
                let mut x = n - 1;
                for (j, &p) in probs.iter().enumerate() {
                    acc += p;
                    if u < acc {
                        x = j;
                        break;
                    }
                }
                // The validated longest conditional table has n^(m+1)
                // cells, bounding the multiply before the rolling remainder.
                ctx = (ctx * n + x) % context_states;
                data.push(symbols[x]);
            }
            offsets.push(data.len() as i64);
        }
        (offsets, data)
    }

    pub fn alphabet_len(&self) -> usize {
        self.alphabet_len
    }

    /// The Markov order `m`.
    pub fn order(&self) -> usize {
        self.kmers.len() - 1
    }

    /// Frequencies of k-mers as given, `k = 1 ..= order + 1`.
    pub fn kmer_frequencies(&self, k: usize) -> Option<&[f64]> {
        self.kmers.get(k.checked_sub(1)?).map(Vec::as_slice)
    }

    /// Normalized letter frequencies.
    pub fn letters(&self) -> &[f64] {
        &self.conditional[0]
    }

    /// `b(x | context)` for a context of `c ≤ order` symbols: `n^(c + 1)`
    /// probabilities, index `context * n + x`.
    pub fn conditional(&self, c: usize) -> &[f64] {
        &self.conditional[c]
    }
}

/// MEME's `normalize_subarray`: divides by the sum unless it is within
/// `tolerance` of 1.
pub(crate) fn meme_normalize(values: &mut [f64], tolerance: f64) -> f64 {
    let mut total = 0.0;
    for v in values.iter() {
        total += v;
    }
    if (1.0 - tolerance > total) || (1.0 + tolerance < total) {
        for v in values.iter_mut() {
            *v /= total;
        }
    }
    total
}

/// Writes the base-`n` digits of `index` into `out`, most significant first.
pub(crate) fn decode(mut index: usize, n: usize, out: &mut [usize]) {
    for d in out.iter_mut().rev() {
        *d = index % n;
        index /= n;
    }
}

fn finite_average(a: f64, b: f64) -> f64 {
    let sum = a + b;
    if sum.is_finite() {
        sum / 2.0
    } else {
        a / 2.0 + b / 2.0
    }
}

fn invalid(detail: impl Into<String>) -> Error {
    Error::invalid("motif background", detail)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn conditionals_normalize_per_context() {
        let bg =
            Background::from_kmer_frequencies(2, vec![vec![1.0, 3.0], vec![1.0, 1.0, 2.0, 6.0]])
                .unwrap();
        assert_eq!(bg.order(), 1);
        assert_eq!(bg.letters(), &[0.25, 0.75]);
        assert_eq!(bg.conditional(1), &[0.5, 0.5, 0.25, 0.75]);
    }

    #[test]
    fn large_finite_frequencies_normalize_and_metadata_bounds_are_checked() {
        let bg = Background::frequencies(&[1e308; 4]).unwrap();
        assert_eq!(bg.letters(), &[0.25; 4]);
        assert_eq!(bg.symmetrized(&[3, 2, 1, 0]).unwrap().letters(), &[0.25; 4]);
        assert_eq!(bg.meme_prepared(Some(&[3, 2, 1, 0])).letters(), &[0.25; 4]);
        assert_eq!(bg.meme_prepared(None).letters(), &[0.25; 4]);
        assert!(bg.symmetrized(&[255; 4]).is_err());
        assert!(bg.symmetrized(&[1, 2, 3, 0]).is_err());
        assert_eq!(bg.truncated(usize::MAX), bg);
        assert!(
            Background::estimate(
                &Alphabet::dna(),
                Strings {
                    offsets: &[0],
                    data: &[],
                },
                usize::MAX,
                1.0
            )
            .is_err()
        );
    }

    #[test]
    fn meme_file_parses_by_order() {
        let text = "# order 0\nA 0.3\nC 0.2\nG 0.2\nT 0.3\n";
        let bg = Background::from_meme(&Alphabet::dna(), text).unwrap();
        assert_eq!(bg.order(), 0);
        assert_eq!(bg.kmer_frequencies(1).unwrap(), &[0.3, 0.2, 0.2, 0.3]);
        assert!(Background::from_meme(&Alphabet::dna(), "A 0.5\nC 0.5\n").is_err());
    }

    #[test]
    fn estimate_counts_runs() {
        let data = b"ACNGT";
        let offsets = [0i64, 5];
        let bg = Background::estimate(
            &Alphabet::dna(),
            Strings {
                offsets: &offsets,
                data,
            },
            1,
            1.0,
        )
        .unwrap();
        assert_eq!(bg.kmer_frequencies(1).unwrap(), &[2.0, 2.0, 2.0, 2.0]);
        let di = bg.kmer_frequencies(2).unwrap();
        assert_eq!(di[1], 2.0); // AC
        assert_eq!(di[2 * 4 + 3], 2.0); // GT
        assert_eq!(di[4 + 2], 1.0); // C and G are separated by N
    }

    #[test]
    fn samples_follow_the_chain() {
        let bg = Background::from_kmer_frequencies(
            4,
            vec![vec![1.0; 4], {
                // After A always C; otherwise uniform.
                let mut f = vec![1.0; 16];
                f[0] = 1e-9;
                f[2] = 1e-9;
                f[3] = 1e-9;
                f[1] = 1.0;
                f
            }],
        )
        .unwrap();
        let (offsets, data) = bg.sample(&Alphabet::dna(), &[5000, 0, 10], 7);
        assert_eq!(offsets, vec![0, 5000, 5000, 5010]);
        for w in data[..5000].windows(2) {
            if w[0] == b'A' {
                assert_eq!(w[1], b'C');
            }
        }
        assert_eq!(
            bg.sample(&Alphabet::dna(), &[50], 7),
            bg.sample(&Alphabet::dna(), &[50], 7)
        );
    }

    #[test]
    fn symmetrized_averages_reverse_complements() {
        let bg = Background::frequencies(&[0.5, 0.125, 0.25, 0.125]).unwrap();
        let s = bg.symmetrized(&[3, 2, 1, 0]).unwrap();
        assert_eq!(
            s.kmer_frequencies(1).unwrap(),
            &[0.3125, 0.1875, 0.1875, 0.3125]
        );
    }
    #[test]
    fn sampling_rolls_two_symbol_context_and_restarts_at_each_sequence() {
        let alphabet = Alphabet::dna();
        let mut pairs = vec![1.0; 16];
        pairs[1] = 1e300;
        let mut triples = vec![1.0; 64];
        for context in 0..16 {
            let next = match context {
                1 => 3,
                7 => 2,
                14 => 0,
                8 => 1,
                _ => (context % 4 + 1) % 4,
            };
            triples[context * 4 + next] = 1e300;
        }
        let background =
            Background::from_kmer_frequencies(4, vec![vec![1e300, 1.0, 1.0, 1.0], pairs, triples])
                .unwrap();
        let (offsets, data) = background.sample(&alphabet, &[0, 12, 4], 7);
        assert_eq!(offsets, [0, 0, 12, 16]);
        assert_eq!(data, b"ACTGACTGACTGACTG");
    }

    #[test]
    #[should_panic(expected = "sample lengths exceed")]
    fn sampling_checks_the_total_before_allocating() {
        let alphabet = Alphabet::dna();
        Background::uniform(&alphabet).sample(&alphabet, &[usize::MAX, 1], 0);
    }

    #[test]
    #[should_panic(expected = "sample alphabet and background sizes differ")]
    fn sampling_checks_the_alphabet_before_allocating() {
        Background::frequencies(&[0.5, 0.5])
            .unwrap()
            .sample(&Alphabet::dna(), &[], 0);
    }
}
