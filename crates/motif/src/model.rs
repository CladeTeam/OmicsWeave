//! Compiled motifs: steps ①–⑤ for every strand, with score and p-value
//! conversions. Usable on their own for threshold calibration.

use omicsweave_core::{Error, Result, Strand};
use omicsweave_ffm::{Alphabet, MotifView};

use crate::background::Background;
use crate::factor::{estimate, reverse};
use crate::null::{NullDist, null_distribution, null_with_letters};
use crate::score::{IntModel, ScoreModel, contrast, quantize};
use crate::scoring::{Rounding, Scoring, Strands};

/// A significance threshold on windows.
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum Threshold {
    /// Windows whose reported p-value is at most this.
    PValue(f64),
    /// Windows whose reported score is at least this.
    Score(f64),
}

/// One strand of a motif for windows with a given amount of background
/// context: its real and integer score models and null distribution.
#[derive(Clone, Debug)]
pub struct Variant {
    real: ScoreModel,
    int: IntModel,
    null: NullDist,
    rounding: Rounding,
}

impl Variant {
    pub fn score_model(&self) -> &ScoreModel {
        &self.real
    }

    pub fn int_model(&self) -> &IntModel {
        &self.int
    }

    pub fn null(&self) -> &NullDist {
        &self.null
    }

    /// The reported score of integer score `s`.
    pub fn score(&self, s: i32) -> f64 {
        self.rounding.apply(self.int.raw(s))
    }

    /// The reported p-value of integer score `s`.
    pub fn pvalue(&self, s: i32) -> f64 {
        self.rounding.apply(self.null.pvalue(s))
    }

    /// Whether integer score `s` passes `threshold`.
    pub fn meets(&self, s: i32, threshold: Threshold) -> bool {
        match threshold {
            Threshold::PValue(p) => self.pvalue(s) <= p,
            Threshold::Score(x) => self.score(s) >= x,
        }
    }

    /// The smallest integer score passing `threshold`; `None` when none does.
    pub fn cutoff(&self, threshold: Threshold) -> Option<i32> {
        match threshold {
            Threshold::PValue(p) => self.null.cutoff(p, self.rounding),
            Threshold::Score(x) => {
                let (mut lo, mut hi) = (0i32, self.int.max + 1);
                while lo < hi {
                    let mid = lo + (hi - lo) / 2;
                    if self.score(mid) >= x {
                        hi = mid;
                    } else {
                        lo = mid + 1;
                    }
                }
                (lo <= self.int.max).then_some(lo)
            }
        }
    }

    /// The smallest p-value any window reaches.
    pub fn p_min(&self) -> f64 {
        self.null.max_reachable().map_or(1.0, |s| self.pvalue(s))
    }

    /// Background probability that a window passes `threshold`.
    pub fn alpha_eff(&self, threshold: Threshold) -> f64 {
        self.cutoff(threshold).map_or(0.0, |s| self.null.pvalue(s))
    }
}

/// A motif on one strand. `variants()[c]` scores windows with `c` symbols of
/// background context available; the last one is the ordinary case.
#[derive(Clone, Debug)]
pub struct StrandModel {
    strand: Strand,
    variants: Vec<Variant>,
}

impl StrandModel {
    pub fn strand(&self) -> Strand {
        self.strand
    }

    /// The model for windows with the full background context.
    pub fn full(&self) -> &Variant {
        self.variants.last().expect("at least one variant")
    }

    /// The model for windows with `context` symbols of context.
    pub fn variant(&self, context: usize) -> &Variant {
        &self.variants[context.min(self.variants.len() - 1)]
    }

    pub fn variants(&self) -> &[Variant] {
        &self.variants
    }
}

/// A motif compiled under a background and a scoring policy.
#[derive(Clone, Debug)]
pub struct MotifModel {
    width: u32,
    strands: Vec<StrandModel>,
}

impl MotifModel {
    /// Compiles one motif. The background is prepared by `scoring` first.
    pub fn compile(
        view: &MotifView<'_>,
        alphabet: &Alphabet,
        background: &Background,
        scoring: &Scoring,
        strands: Strands,
    ) -> Result<MotifModel> {
        scoring.check()?;
        scoring.ambiguity.check(alphabet)?;
        check_strands(alphabet, strands)?;
        let prepared = scoring.prepare_background(background, alphabet, strands)?;
        MotifModel::build(view, alphabet, &prepared, scoring, strands)
    }

    /// Compiles with an already prepared background.
    pub(crate) fn build(
        view: &MotifView<'_>,
        alphabet: &Alphabet,
        background: &Background,
        scoring: &Scoring,
        strands: Strands,
    ) -> Result<MotifModel> {
        let n = alphabet.len();
        if u64::from(view.width()) + background.order() as u64 > i32::MAX as u64 {
            return Err(Error::invalid(
                "motif model",
                "width plus context exceeds 2^31 - 1",
            ));
        }
        let forward = estimate(view, n, background.letters(), scoring)?;
        let mut out: Vec<StrandModel> = Vec::with_capacity(2);
        let null_of = |int: &IntModel| -> Result<NullDist> {
            match &scoring.null_letters {
                Some(letters) if background.order() == 0 => null_with_letters(int, letters),
                Some(_) => Err(Error::invalid(
                    "motif null distribution",
                    "explicit null letters need a 0-order background",
                )),
                None => null_distribution(int, background, scoring.null_budget),
            }
        };
        let mut orient = |s: &crate::factor::Smoothed, strand: Strand| -> Result<()> {
            let m = background.order();
            let weights = scoring
                .null_letters
                .as_deref()
                .unwrap_or(background.letters());
            let full_real = contrast(s, background, m, scoring.log, &scoring.ambiguity, weights);
            let full_int = quantize(&full_real, scoring.quantizer, n, None)?;
            let mut variants = Vec::with_capacity(m + 1);
            for c in 0..m {
                let real = contrast(s, background, c, scoring.log, &scoring.ambiguity, weights);
                let int = quantize(&real, scoring.quantizer, n, Some(full_int.scale))?;
                let null = null_of(&int)?;
                variants.push(Variant {
                    real,
                    int,
                    null,
                    rounding: scoring.rounding,
                });
            }
            let shared = (strand == Strand::Reverse && scoring.shared_reverse_null)
                .then(|| out.first().map(|f| f.full().null.clone()))
                .flatten();
            let null = match shared {
                Some(null) => null,
                None => null_of(&full_int)?,
            };
            variants.push(Variant {
                real: full_real,
                int: full_int,
                null,
                rounding: scoring.rounding,
            });
            out.push(StrandModel { strand, variants });
            Ok(())
        };
        if strands.forward() {
            orient(&forward, Strand::Forward)?;
        }
        if strands.reverse() {
            let complement = alphabet
                .complement()
                .ok_or_else(|| Error::invalid("motif strands", "alphabet has no complement"))?;
            let rev = reverse(&forward, n, complement);
            orient(&rev, Strand::Reverse)?;
        }
        Ok(MotifModel {
            width: view.width(),
            strands: out,
        })
    }

    pub fn width(&self) -> u32 {
        self.width
    }

    /// One model per scanned strand, forward first.
    pub fn strands(&self) -> &[StrandModel] {
        &self.strands
    }

    pub fn strand(&self, strand: Strand) -> Option<&StrandModel> {
        self.strands.iter().find(|s| s.strand == strand)
    }
}

pub(crate) fn check_strands(alphabet: &Alphabet, strands: Strands) -> Result<()> {
    if strands.reverse() && alphabet.complement().is_none() {
        return Err(Error::invalid(
            "motif strands",
            format!("alphabet {:?} has no complement", alphabet.name()),
        ));
    }
    Ok(())
}
