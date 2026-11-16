//! Scoring policies: how stored frequencies become window scores and how
//! scores and p-values are reported.

use std::fmt::Debug;
use std::sync::Arc;

use omicsweave_core::{Error, Result};
use omicsweave_ffm::Alphabet;

use crate::background::Background;

/// MEME's truncated `log2(e)`; its bits are `ln(x)` times this.
#[allow(clippy::approx_constant)]
const MEME_LOG2_E: f64 = 1.44269504;

/// MEME's truncated `ln(10)`, used by its rounding macro.
#[allow(clippy::approx_constant)]
const MEME_LN_10: f64 = 2.30258509299405;

/// Which strands of the sequence a motif is tested on.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub enum Strands {
    Forward,
    Reverse,
    Both,
}

impl Strands {
    pub fn forward(self) -> bool {
        matches!(self, Strands::Forward | Strands::Both)
    }

    pub fn reverse(self) -> bool {
        matches!(self, Strands::Reverse | Strands::Both)
    }
}

/// The logarithm scores are expressed in.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub enum LogScale {
    /// Bits: `log2(x)`.
    Bits,
    /// Bits as MEME computes them, `ln(x) · 1.44269504`.
    MemeBits,
    /// Nats: `ln(x)`.
    Natural,
    /// Hundredths of a bit as MEME computes them, `ln(x) / ln(2) · 100`.
    Centibits,
}

impl LogScale {
    /// The logarithm of `x`; `-∞` for `x ≤ 0`.
    pub fn log(self, x: f64) -> f64 {
        if x <= 0.0 {
            return f64::NEG_INFINITY;
        }
        match self {
            LogScale::Bits => x.log2(),
            LogScale::MemeBits => x.ln() * MEME_LOG2_E,
            LogScale::Natural => x.ln(),
            LogScale::Centibits => (x.ln() / 2f64.ln()) * 100.0,
        }
    }
}

/// How real-valued scores become the integers the scan and the null
/// distribution work with.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub enum Quantizer {
    /// MEME's `scale_score_matrix`: one offset `floor(min)` for every term
    /// and an integer scale `floor(range / (max − offset))`, so each term's
    /// cells lie in `0 ..= range`.
    Fimo { range: u32 },
    /// Each term shifted by its own minimum, all scaled together so the best
    /// window scores `states − 1`.
    Precision { states: u32 },
    /// MAST's `scale_lo`: offset `min` (not floored) and the real scale
    /// `range / (max − min)` kept to 14 significant digits, so each term's
    /// cells lie in `0 ..= range`.
    Mast { range: u32 },
}

/// Rounding applied to reported scores and p-values (and so to p-value
/// thresholds).
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub enum Rounding {
    None,
    /// MEME's `RND(x, d)`: `d` significant digits, ties to even.
    Significant(u32),
}

impl Rounding {
    pub fn apply(self, x: f64) -> f64 {
        match self {
            Rounding::None => x,
            Rounding::Significant(d) => meme_round(x, d),
        }
    }
}

/// MEME's `RND` macro.
pub(crate) fn meme_round(x: f64, digits: u32) -> f64 {
    if !x.is_finite() {
        return x;
    }
    let d = digits as f64;
    if x > 0.0 {
        let exponent = (d - 1.0 - x.ln() / MEME_LN_10).ceil();
        let z = 10f64.powf(exponent);
        if !z.is_finite() || z == 0.0 {
            return round_decimal_step(x, exponent);
        }
        (z * x).round_ties_even() / z
    } else if x < 0.0 {
        let exponent = (d - 1.0 - (-x).ln() / MEME_LN_10).ceil();
        let z = 10f64.powf(exponent);
        if !z.is_finite() || z == 0.0 {
            return round_decimal_step(x, exponent);
        }
        -(z * -x).round_ties_even() / z
    } else {
        0.0
    }
}

fn round_decimal_step(x: f64, exponent: f64) -> f64 {
    let step = 10f64.powf(-exponent);
    if step == 0.0 || !step.is_finite() {
        return x;
    }
    (x / step).round_ties_even() * step
}

/// How the given background is prepared before use.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub enum BackgroundRule {
    /// Used as given.
    AsGiven,
    /// As MEME prepares it: cut to 0-order; averaged with the complement
    /// when both strands are scanned; every letter raised by 5e-7 and the
    /// whole normalized.
    Meme,
}

/// How a window containing an ambiguity code is scored.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub enum AmbiguityRule {
    /// Not scored, like any window with an unknown symbol (FIMO).
    Skip,
    /// The mean of the member symbols' scores, weighted by the null letter
    /// probabilities (MAST).
    Mean,
    /// The sum of the member symbols' scores (MCAST).
    Sum,
}

/// Ambiguity codes: bytes outside the alphabet that stand for a set of its
/// symbols. Only single-position terms score them; a window where an
/// ambiguity code falls in a wider term is not scored.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Ambiguity {
    pub rule: AmbiguityRule,
    /// Each code: its byte and the global codes it stands for, ascending.
    pub codes: Vec<(u8, Vec<u8>)>,
    /// A code (by its byte) that every other byte outside the alphabet
    /// stands for, as tools that read any unknown character as `N` do.
    pub unknown: Option<u8>,
}

impl Ambiguity {
    /// No ambiguity codes.
    pub fn skip() -> Ambiguity {
        Ambiguity {
            rule: AmbiguityRule::Skip,
            codes: Vec::new(),
            unknown: None,
        }
    }

    /// The IUPAC nucleotide codes `R Y S W K M B D H V N` (either case) over
    /// an alphabet whose symbols are `ACGT` or `ACGU`. Bytes the alphabet
    /// already maps are left to it.
    pub fn iupac(alphabet: &Alphabet, rule: AmbiguityRule) -> Result<Ambiguity> {
        let s = alphabet.symbols();
        if s != b"ACGT" && s != b"ACGU" {
            return Err(Error::invalid(
                "motif ambiguity",
                "IUPAC codes need the symbols ACGT or ACGU",
            ));
        }
        let sets: [(u8, &[u8]); 11] = [
            (b'R', &[0, 2]),
            (b'Y', &[1, 3]),
            (b'S', &[1, 2]),
            (b'W', &[0, 3]),
            (b'K', &[2, 3]),
            (b'M', &[0, 1]),
            (b'B', &[1, 2, 3]),
            (b'D', &[0, 2, 3]),
            (b'H', &[0, 1, 3]),
            (b'V', &[0, 1, 2]),
            (b'N', &[0, 1, 2, 3]),
        ];
        let mut codes = Vec::new();
        for (byte, set) in sets {
            for b in [byte, byte.to_ascii_lowercase()] {
                if alphabet.lut()[b as usize] == omicsweave_ffm::INVALID {
                    codes.push((b, set.to_vec()));
                }
            }
        }
        Ok(Ambiguity {
            rule,
            codes,
            unknown: None,
        })
    }

    /// Adds `byte` as another spelling of the set of code `of` (e.g. `X` for
    /// `N`).
    pub fn with_alias(mut self, byte: u8, of: u8) -> Result<Ambiguity> {
        let set = self
            .codes
            .iter()
            .find(|(b, _)| *b == of)
            .map(|(_, s)| s.clone())
            .ok_or_else(|| {
                Error::invalid("motif ambiguity", format!("{:?} is not a code", of as char))
            })?;
        if self.codes.iter().any(|(b, _)| *b == byte) {
            return Err(Error::invalid(
                "motif ambiguity",
                format!("{:?} is already a code", byte as char),
            ));
        }
        self.codes.push((byte, set));
        Ok(self)
    }

    /// Reads every other byte outside the alphabet as code `of`.
    pub fn with_unknown(mut self, of: u8) -> Result<Ambiguity> {
        if !self.codes.iter().any(|(b, _)| *b == of) {
            return Err(Error::invalid(
                "motif ambiguity",
                format!("{:?} is not a code", of as char),
            ));
        }
        self.unknown = Some(of);
        Ok(self)
    }

    pub(crate) fn check(&self, alphabet: &Alphabet) -> Result<()> {
        if let Some(u) = self.unknown
            && !self.codes.iter().any(|(b, _)| *b == u)
        {
            return Err(Error::invalid(
                "motif ambiguity",
                format!(
                    "unknown bytes stand for {:?}, which is not a code",
                    u as char
                ),
            ));
        }
        if self.codes.len() > 254 {
            return Err(Error::invalid("motif ambiguity", "more than 254 codes"));
        }
        let mut seen = [false; 256];
        for (b, set) in &self.codes {
            if std::mem::replace(&mut seen[*b as usize], true) {
                return Err(Error::invalid(
                    "motif ambiguity",
                    format!("{:?} is listed more than once", *b as char),
                ));
            }
            if alphabet.lut()[*b as usize] != omicsweave_ffm::INVALID {
                return Err(Error::invalid(
                    "motif ambiguity",
                    format!("{:?} is an alphabet symbol", *b as char),
                ));
            }
            if set.is_empty()
                || set.windows(2).any(|w| w[0] >= w[1])
                || set.iter().any(|&c| c as usize >= alphabet.len())
            {
                return Err(Error::invalid(
                    "motif ambiguity",
                    format!(
                        "the set of {:?} is not ascending alphabet codes",
                        *b as char
                    ),
                ));
            }
        }
        Ok(())
    }
}

/// A prior for single-position columns that replaces the pseudocount: it
/// maps a column's normalized frequencies (over the motif's local codes) to
/// smoothed probabilities.
pub trait ColumnPrior: Send + Sync + Debug {
    /// `column`: normalized frequencies; `sites`: observations behind them
    /// (after the policy's `default_sites` / `truncate_sites`); `letters`:
    /// the background's letter distribution over the same codes.
    fn smooth(&self, column: &[f64], sites: f64, letters: &[f64]) -> Vec<f64>;
}

/// A scoring policy. [`Scoring::fimo`] reproduces FIMO; the default is the
/// exact policy described in `docs/motif.md`.
#[derive(Clone, Debug)]
pub struct Scoring {
    /// Total prior mass α of the pseudocount; 0 disables smoothing.
    pub pseudocount: f64,
    /// Observations assumed for a motif without `sites`.
    pub default_sites: f64,
    /// Use the integer part of `sites`.
    pub truncate_sites: bool,
    /// Normalization groups whose sum is within this of 1 are used as is.
    pub normalize_tolerance: f64,
    pub log: LogScale,
    pub quantizer: Quantizer,
    pub rounding: Rounding,
    pub background: BackgroundRule,
    /// Most elementary steps one null distribution may take; larger models
    /// fail to compile instead of being approximated.
    pub null_budget: u64,
    /// Letter probabilities the null distribution is computed with, used as
    /// given (not normalized); `None` uses the prepared background. Only for
    /// single-position terms under a 0-order background.
    pub null_letters: Option<Vec<f64>>,
    /// The reverse strand reuses the forward strand's null distribution, as
    /// tools that build one p-value table per motif do. Only valid when the
    /// reverse model is a permutation of the forward one (single-position
    /// terms, a background symmetric under complement).
    pub shared_reverse_null: bool,
    /// Replaces the pseudocount for single-position joint factors.
    pub column_prior: Option<Arc<dyn ColumnPrior>>,
    /// Ambiguity codes and how windows containing them are scored.
    pub ambiguity: Ambiguity,
}

impl Default for Scoring {
    fn default() -> Scoring {
        Scoring::exact()
    }
}

impl Scoring {
    /// FIMO 5.5.9: pseudocount 0.1, 20 sites by default (integer part),
    /// rows within 1e-5 of 1 kept, MEME bits, `Fimo { range: 1000 }`, ten
    /// significant digits, MEME's background preparation.
    pub fn fimo() -> Scoring {
        Scoring {
            pseudocount: 0.1,
            default_sites: 20.0,
            truncate_sites: true,
            normalize_tolerance: 0.00001,
            log: LogScale::MemeBits,
            quantizer: Quantizer::Fimo { range: 1000 },
            rounding: Rounding::Significant(10),
            background: BackgroundRule::Meme,
            null_budget: 1 << 33,
            null_letters: None,
            shared_reverse_null: false,
            column_prior: None,
            ambiguity: Ambiguity::skip(),
        }
    }

    /// Exact normalization, bits, `Precision { states: 1 << 16 }`, no
    /// rounding, background used as given.
    pub fn exact() -> Scoring {
        Scoring {
            pseudocount: 0.1,
            default_sites: 20.0,
            truncate_sites: false,
            normalize_tolerance: 0.0,
            log: LogScale::Bits,
            quantizer: Quantizer::Precision { states: 1 << 16 },
            rounding: Rounding::None,
            background: BackgroundRule::AsGiven,
            null_budget: 1 << 33,
            null_letters: None,
            shared_reverse_null: false,
            column_prior: None,
            ambiguity: Ambiguity::skip(),
        }
    }

    pub(crate) fn check(&self) -> Result<()> {
        let bad = |detail: String| Err(Error::invalid("motif scoring", detail));
        if !(self.pseudocount.is_finite() && self.pseudocount >= 0.0) {
            return bad(format!(
                "pseudocount {} is not finite and >= 0",
                self.pseudocount
            ));
        }
        if !(self.default_sites.is_finite() && self.default_sites > 0.0) {
            return bad(format!("default sites {} is not > 0", self.default_sites));
        }
        if !(self.normalize_tolerance.is_finite() && self.normalize_tolerance >= 0.0) {
            return bad("normalize tolerance is not >= 0".into());
        }
        match self.quantizer {
            Quantizer::Fimo { range } if range == 0 || range > 1 << 20 => {
                return bad(format!("Fimo range {range} is not in 1..=2^20"));
            }
            Quantizer::Precision { states } if !(2..=1 << 24).contains(&states) => {
                return bad(format!("precision states {states} is not in 2..=2^24"));
            }
            Quantizer::Mast { range } if range == 0 || range > 1 << 20 => {
                return bad(format!("MAST range {range} is not in 1..=2^20"));
            }
            _ => {}
        }
        if let Some(l) = &self.null_letters
            && l.iter().any(|v| !(v.is_finite() && *v >= 0.0))
        {
            return bad("null letters must be finite and >= 0".into());
        }
        Ok(())
    }

    /// The background a scan uses, after this policy's preparation.
    pub fn prepare_background(
        &self,
        background: &Background,
        alphabet: &Alphabet,
        strands: Strands,
    ) -> Result<Background> {
        if background.alphabet_len() != alphabet.len() {
            return Err(Error::invalid(
                "motif background",
                format!(
                    "{} symbols, alphabet has {}",
                    background.alphabet_len(),
                    alphabet.len()
                ),
            ));
        }
        Ok(match self.background {
            BackgroundRule::AsGiven => background.clone(),
            BackgroundRule::Meme => {
                let complement = if strands == Strands::Both {
                    alphabet.complement()
                } else {
                    None
                };
                background.meme_prepared(complement)
            }
        })
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn meme_rounding_keeps_ten_digits() {
        let r = Rounding::Significant(10);
        assert_eq!(r.apply(1.234567890123e-5), 1.23456789e-5);
        assert_eq!(r.apply(-12.3456789012345), -12.3456789);
        assert_eq!(r.apply(0.0), 0.0);
    }

    #[test]
    fn significant_rounding_handles_subnormals_and_extreme_values() {
        let rounding = Rounding::Significant(6);
        for magnitude in [1e-308, 3.456789e-310, 1e-320, f64::from_bits(1)] {
            for sign in [1.0, -1.0] {
                let value = sign * magnitude;
                let rounded = rounding.apply(value);
                assert!(rounded.is_finite());
                assert_eq!(rounded.is_sign_negative(), value.is_sign_negative());
                let tolerance = (magnitude * 1e-5).max(f64::from_bits(1));
                assert!((rounded - value).abs() <= tolerance);
            }
        }
        assert_eq!(rounding.apply(f64::from_bits(1)), f64::from_bits(1));
        assert!(!rounding.apply(f64::MAX).is_nan());
        assert!(!rounding.apply(-f64::MAX).is_nan());
        assert_eq!(rounding.apply(f64::INFINITY), f64::INFINITY);
        assert!(rounding.apply(f64::NAN).is_nan());
        assert_eq!(meme_round(0.12345, u32::MAX), 0.12345);
    }

    #[test]
    fn log_scales() {
        assert_eq!(LogScale::Bits.log(8.0), 3.0);
        assert_eq!(LogScale::Bits.log(0.0), f64::NEG_INFINITY);
        assert!((LogScale::MemeBits.log(8.0) - 3.0).abs() < 1e-8);
    }
}
