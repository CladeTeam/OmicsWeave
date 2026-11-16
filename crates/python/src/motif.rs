//! Motif matching. Policies, statistics and reductions come in as the
//! package's dataclasses, read attribute by attribute; native steps are
//! told apart by their `kind`, and any other callable is called back with
//! the GIL (its first exception is raised once the native call returns).
//! Presets and defaults are this crate's, handed to Python as plain values.

use std::sync::{Arc, Mutex};

use numpy::{IntoPyArray, PyArray1, PyArray2, PyArrayMethods, PyReadonlyArray1};
use omicsweave_core::{Error, Result, Strand, Strings};
use omicsweave_ffm::Alphabet;
use omicsweave_motif::mast::{
    self, MastHit, MastOptions, MastReport, MastSequence, MotifBest, StrandMode,
};
use omicsweave_motif::mcast::{
    self, BackgroundNull, Cluster, ClusterNull, ClusterOptions, ClusterSearch, ClusterTail,
    Clusters, GcBinnedExponential, MemeQValue, MemeSynthetic, OutputFilter, Report, ScoredCluster,
    Statistics, Unconditional,
};
use omicsweave_motif::profile::{
    self, Alignment, AllAlignments, Column, DirichletMixture, Glam2Text,
    ProfileProgram as Profiles, ProfileScoring, TopAlignments,
};
use omicsweave_motif::stats::{
    Adjust, BenjaminiHochberg, BestOf, Combine, EValue, EmpiricalTail, MemeBestOf, MinBonferroni,
    Product, Sidak, TailFit,
};
use omicsweave_motif::substitution::SubstitutionPrior;
use omicsweave_motif::{
    Ambiguity, AmbiguityRule, Background, BackgroundRule, BestEntry, ColumnPrior, CountMatrix,
    Demand, Dense, DenseScores, HitEntry, Hits, LogScale, Program, Quantizer, Reducer, Rounding,
    ScanInfo, ScanOptions, Scoring, SequenceBest, SequenceChunk, SequenceCounts, Strands,
    Threshold, Windows, fimo,
};
use pyo3::exceptions::{PyIndexError, PyTypeError, PyValueError};
use pyo3::prelude::*;
use pyo3::types::PyDict;

use crate::error::OrRaise;
use crate::ffm::{self as conv, attr, byte_of, bytes_of, floats, motif_set, views};

type Floats<'py> = Bound<'py, PyArray1<f64>>;
type StringArrays<'py> = (PyReadonlyArray1<'py, i64>, PyReadonlyArray1<'py, u8>);
type StringColumns<'py> = (Bound<'py, PyArray1<i64>>, Bound<'py, PyArray1<u8>>);

fn invalid(detail: impl Into<String>) -> PyErr {
    PyValueError::new_err(detail.into())
}

fn slice<'a, T: numpy::Element>(a: &'a PyReadonlyArray1<'_, T>, what: &str) -> PyResult<&'a [T]> {
    crate::shared::slice(a, what)
}

fn strings<'a>(s: &'a StringArrays<'_>, what: &str) -> PyResult<Strings<'a>> {
    Ok(Strings {
        offsets: slice(&s.0, what)?,
        data: slice(&s.1, what)?,
    })
}

fn strands(s: &str) -> PyResult<Strands> {
    Ok(match s {
        "both" => Strands::Both,
        "+" => Strands::Forward,
        "-" => Strands::Reverse,
        _ => return Err(invalid(format!("strands must be both, + or -, not {s:?}"))),
    })
}

fn strand(s: &str) -> PyResult<Strand> {
    Ok(match s {
        "+" => Strand::Forward,
        "-" => Strand::Reverse,
        _ => return Err(invalid(format!("strand must be + or -, not {s:?}"))),
    })
}

fn strand_code(s: Strand) -> i8 {
    match s {
        Strand::Forward => 1,
        Strand::Reverse => -1,
        _ => 0,
    }
}

fn strand_of(code: i8) -> PyResult<Strand> {
    match code {
        1 => Ok(Strand::Forward),
        -1 => Ok(Strand::Reverse),
        _ => Err(invalid(format!("strand codes are 1 and -1, not {code}"))),
    }
}

fn kind(o: &Bound<'_, PyAny>) -> Option<String> {
    o.getattr("kind").ok().and_then(|k| k.extract().ok())
}

fn column<'py, T: numpy::Element + Copy>(o: &Bound<'py, PyAny>, name: &str) -> PyResult<Vec<T>> {
    let a: PyReadonlyArray1<'py, T> = o.getattr(name)?.extract().map_err(|_| {
        PyTypeError::new_err(format!(
            "column {name} must be a numpy array of {}",
            std::any::type_name::<T>()
        ))
    })?;
    Ok(a.as_array().iter().copied().collect())
}

// ---------------------------------------------------------------------------
// Python callbacks

/// The first exception raised by a callback during one native call.
#[derive(Clone, Default)]
struct Errors(Arc<Mutex<Option<PyErr>>>);

impl Errors {
    fn set(&self, e: PyErr) {
        let mut slot = self.0.lock().expect("callback errors lock");
        if slot.is_none() {
            *slot = Some(e);
        }
    }

    fn check(&self) -> PyResult<()> {
        match self.0.lock().expect("callback errors lock").take() {
            Some(e) => Err(e),
            None => Ok(()),
        }
    }
}

/// A Python callable standing in for a native step.
struct Callback {
    f: Py<PyAny>,
    errors: Errors,
}

impl std::fmt::Debug for Callback {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.write_str("Callback")
    }
}

impl Callback {
    fn new(f: &Bound<'_, PyAny>, errors: &Errors, what: &str) -> PyResult<Callback> {
        if !f.is_callable() {
            return Err(PyTypeError::new_err(format!(
                "{what} must be one of the package's steps or a callable"
            )));
        }
        Ok(Callback {
            f: f.clone().unbind(),
            errors: errors.clone(),
        })
    }

    fn call<T>(
        &self,
        fallback: T,
        f: impl FnOnce(Python<'_>, &Bound<'_, PyAny>) -> PyResult<T>,
    ) -> T {
        Python::attach(|py| f(py, self.f.bind(py))).unwrap_or_else(|e| {
            self.errors.set(e);
            fallback
        })
    }
}

fn vec_of(o: &Bound<'_, PyAny>) -> PyResult<Vec<f64>> {
    let np = o.py().import("numpy")?;
    let a = np.call_method1("asarray", (o, "float64"))?;
    let a: PyReadonlyArray1<'_, f64> = a.extract()?;
    Ok(a.as_array().iter().copied().collect())
}

impl BestOf for Callback {
    fn pvalue(&self, p: f64, n: f64) -> f64 {
        self.call(f64::NAN, |_, f| f.call1((p, n))?.extract())
    }
}

impl Combine for Callback {
    fn combine(&self, pvalues: &[f64]) -> f64 {
        self.call(f64::NAN, |py, f| {
            f.call1((pvalues.to_vec().into_pyarray(py),))?.extract()
        })
    }
}

impl Adjust for Callback {
    fn adjust(&self, pvalues: &[f64]) -> Vec<f64> {
        self.call(vec![f64::NAN; pvalues.len()], |py, f| {
            let out = vec_of(&f.call1((pvalues.to_vec().into_pyarray(py),))?)?;
            if out.len() != pvalues.len() {
                return Err(invalid(format!(
                    "an adjustment returned {} values for {} p-values",
                    out.len(),
                    pvalues.len()
                )));
            }
            Ok(out)
        })
    }
}

impl TailFit for Callback {
    fn fit(&self, samples: &[f64]) -> Result<Box<dyn Fn(f64) -> f64 + Send + Sync>> {
        let tail = self.call(None, |py, f| {
            Ok(Some(
                f.call1((samples.to_vec().into_pyarray(py),))?.unbind(),
            ))
        });
        let Some(tail) = tail else {
            return Err(Error::invalid("tail fit", "the Python fit raised"));
        };
        let tail = Callback {
            f: tail,
            errors: self.errors.clone(),
        };
        Ok(Box::new(move |x| {
            tail.call(f64::NAN, |_, g| g.call1((x,))?.extract())
        }))
    }
}

impl ColumnPrior for Callback {
    fn smooth(&self, column: &[f64], sites: f64, letters: &[f64]) -> Vec<f64> {
        self.call(column.to_vec(), |py, f| {
            let out = vec_of(&f.call1((
                column.to_vec().into_pyarray(py),
                sites,
                letters.to_vec().into_pyarray(py),
            ))?)?;
            if out.len() != column.len() {
                return Err(invalid(format!(
                    "a column prior returned {} values for {}",
                    out.len(),
                    column.len()
                )));
            }
            Ok(out)
        })
    }
}

/// A boxed tail fit, so [`Unconditional`] can hold any.
struct AnyTail(Box<dyn TailFit>);

impl TailFit for AnyTail {
    fn fit(&self, samples: &[f64]) -> Result<Box<dyn Fn(f64) -> f64 + Send + Sync>> {
        self.0.fit(samples)
    }
}

// ---------------------------------------------------------------------------
// Statistics steps

fn best_of(o: &Bound<'_, PyAny>, errors: &Errors) -> PyResult<Arc<dyn BestOf>> {
    Ok(match kind(o).as_deref() {
        Some("sidak") => Arc::new(Sidak),
        Some("meme") => Arc::new(MemeBestOf),
        _ => Arc::new(Callback::new(o, errors, "best_of")?),
    })
}

fn combine(o: &Bound<'_, PyAny>, errors: &Errors) -> PyResult<Arc<dyn Combine>> {
    Ok(match kind(o).as_deref() {
        Some("product") => Arc::new(Product),
        Some("min-bonferroni") => Arc::new(MinBonferroni),
        _ => Arc::new(Callback::new(o, errors, "combine")?),
    })
}

fn adjust(o: &Bound<'_, PyAny>, errors: &Errors) -> PyResult<Box<dyn Adjust>> {
    Ok(match kind(o).as_deref() {
        Some("evalue") => Box::new(EValue {
            tests: attr(o, "tests")?,
        }),
        Some("benjamini-hochberg") => Box::new(BenjaminiHochberg {
            tests: attr(o, "tests")?,
        }),
        Some("meme-qvalue") => Box::new(MemeQValue {
            bootstraps: attr(o, "bootstraps")?,
            samples: attr(o, "samples")?,
            lambdas: attr(o, "lambdas")?,
            max_lambda: attr(o, "max_lambda")?,
            min_values: attr(o, "min_values")?,
            seed: attr(o, "seed")?,
        }),
        _ => Box::new(Callback::new(o, errors, "adjust")?),
    })
}

fn tail_fit(o: &Bound<'_, PyAny>, errors: &Errors) -> PyResult<Box<dyn TailFit>> {
    Ok(match kind(o).as_deref() {
        Some("empirical") => Box::new(EmpiricalTail),
        _ => Box::new(Callback::new(o, errors, "a tail fit")?),
    })
}

fn cluster_null(o: &Bound<'_, PyAny>, alphabet: &Alphabet) -> PyResult<Box<dyn ClusterNull>> {
    Ok(match kind(o).as_deref() {
        Some("meme-synthetic") => Box::new(MemeSynthetic {
            seed: attr(o, "seed")?,
            round_length: attr(o, "round_length")?,
            min_rounds: attr(o, "min_rounds")?,
            max_rounds: attr(o, "max_rounds")?,
            give_up_after: attr(o, "give_up_after")?,
            want: attr(o, "want")?,
            max_samples: attr(o, "max_samples")?,
        }),
        Some("background") => Box::new(BackgroundNull {
            background: background(alphabet.len(), &o.getattr("background")?)?,
            length: attr(o, "length")?,
            sequences: attr(o, "sequences")?,
            seed: attr(o, "seed")?,
        }),
        _ => {
            return Err(PyTypeError::new_err(
                "a cluster null is MemeSynthetic or BackgroundNull",
            ));
        }
    })
}

fn cluster_tail(o: &Bound<'_, PyAny>, errors: &Errors) -> PyResult<Box<dyn ClusterTail>> {
    Ok(match kind(o).as_deref() {
        Some("gc-binned-exponential") => Box::new(GcBinnedExponential {
            min_scores: attr(o, "min_scores")?,
            max_ratio: attr(o, "max_ratio")?,
            max_gap: attr(o, "max_gap")?,
        }),
        Some("unconditional") => Box::new(Unconditional(AnyTail(tail_fit(
            &o.getattr("fit")?,
            errors,
        )?))),
        _ => {
            return Err(PyTypeError::new_err(
                "a cluster tail is GcBinnedExponential or Unconditional",
            ));
        }
    })
}

fn statistics(o: &Bound<'_, PyAny>, alphabet: &Alphabet, errors: &Errors) -> PyResult<Statistics> {
    Ok(Statistics {
        null: cluster_null(&o.getattr("null")?, alphabet)?,
        tail: cluster_tail(&o.getattr("tail")?, errors)?,
        evalue: adjust(&o.getattr("evalue")?, errors)?,
        qvalue: adjust(&o.getattr("qvalue")?, errors)?,
    })
}

#[pyfunction]
pub fn stats_best_of<'py>(
    py: Python<'py>,
    step: &Bound<'py, PyAny>,
    p: PyReadonlyArray1<'py, f64>,
    n: PyReadonlyArray1<'py, f64>,
) -> PyResult<Floats<'py>> {
    let errors = Errors::default();
    let f = best_of(step, &errors)?;
    let (p, n) = (slice(&p, "p")?, slice(&n, "n")?);
    let out: Vec<f64> = p.iter().zip(n).map(|(&p, &n)| f.pvalue(p, n)).collect();
    errors.check()?;
    Ok(out.into_pyarray(py))
}

#[pyfunction]
pub fn stats_combine(step: &Bound<'_, PyAny>, pvalues: PyReadonlyArray1<'_, f64>) -> PyResult<f64> {
    let errors = Errors::default();
    let out = combine(step, &errors)?.combine(slice(&pvalues, "pvalues")?);
    errors.check()?;
    Ok(out)
}

#[pyfunction]
pub fn stats_adjust<'py>(
    py: Python<'py>,
    step: &Bound<'py, PyAny>,
    pvalues: PyReadonlyArray1<'py, f64>,
) -> PyResult<Floats<'py>> {
    let errors = Errors::default();
    let out = adjust(step, &errors)?.adjust(slice(&pvalues, "pvalues")?);
    errors.check()?;
    Ok(out.into_pyarray(py))
}

/// Each native step's and reduction's defaults, by group.
#[pyfunction]
pub fn motif_defaults(py: Python<'_>) -> PyResult<Bound<'_, PyDict>> {
    let out = PyDict::new(py);
    let group = |pairs: &[(&str, f64)]| -> PyResult<Bound<'_, PyDict>> {
        let d = PyDict::new(py);
        for (k, v) in pairs {
            d.set_item(k, v)?;
        }
        Ok(d)
    };
    let int = |pairs: &[(&str, u64)]| -> PyResult<Bound<'_, PyDict>> {
        let d = PyDict::new(py);
        for (k, v) in pairs {
            d.set_item(k, v)?;
        }
        Ok(d)
    };
    let m = MastOptions::new(StrandMode::Combine);
    out.set_item(
        "mast",
        group(&[
            ("hit_threshold", m.hit_threshold),
            ("evalue_threshold", m.evalue_threshold),
        ])?,
    )?;
    let c = ClusterOptions::default();
    let clusters = group(&[("pthresh", c.pthresh), ("egcost", c.egcost)])?;
    clusters.set_item("max_gap", c.max_gap)?;
    clusters.set_item("overlap", c.overlap)?;
    out.set_item("clusters", clusters)?;
    let q = MemeQValue::default();
    let qvalue = int(&[
        ("bootstraps", q.bootstraps as u64),
        ("samples", q.samples as u64),
        ("lambdas", q.lambdas as u64),
        ("min_values", q.min_values as u64),
        ("seed", q.seed as u64),
    ])?;
    qvalue.set_item("max_lambda", q.max_lambda)?;
    out.set_item("meme_qvalue", qvalue)?;
    let s = MemeSynthetic::mcast(0);
    let synthetic = int(&[
        ("round_length", s.round_length as u64),
        ("min_rounds", s.min_rounds as u64),
        ("max_rounds", s.max_rounds as u64),
        ("give_up_after", s.give_up_after as u64),
        ("max_samples", s.max_samples as u64),
    ])?;
    synthetic.set_item("want", s.want)?;
    out.set_item("meme_synthetic", synthetic)?;
    let g = GcBinnedExponential::mcast(c.max_gap);
    let gc = int(&[("min_scores", g.min_scores as u64)])?;
    gc.set_item("max_ratio", g.max_ratio)?;
    out.set_item("gc_binned_exponential", gc)?;
    out.set_item("glam2scan", int(&[("top", profile::GLAM2_HITS as u64)])?)?;
    let f = OutputFilter::default();
    out.set_item(
        "mcast_output",
        group(&[
            ("evalue", f.evalue),
            ("pvalue", f.pvalue),
            ("qvalue", f.qvalue),
        ])?,
    )?;
    Ok(out)
}

// ---------------------------------------------------------------------------
// Backgrounds

/// A `Background` dataclass (its `kmers`) over `n` symbols.
fn background(n: usize, o: &Bound<'_, PyAny>) -> PyResult<Background> {
    kmers_background(n, attr(o, "kmers")?)
}

fn kmers_background(n: usize, kmers: Vec<PyReadonlyArray1<'_, f64>>) -> PyResult<Background> {
    let kmers = kmers.iter().map(|k| k.as_array().to_vec()).collect();
    Background::from_kmer_frequencies(n, kmers).or_raise()
}

fn kmer_arrays<'py>(py: Python<'py>, b: &Background) -> Vec<Floats<'py>> {
    (1..=b.order() + 1)
        .filter_map(|k| b.kmer_frequencies(k))
        .map(|f| f.to_vec().into_pyarray(py))
        .collect()
}

#[pyfunction]
pub fn background_nrdb<'py>(
    py: Python<'py>,
    alphabet: &Bound<'py, PyAny>,
) -> PyResult<Vec<Floats<'py>>> {
    let a = conv::alphabet(alphabet)?;
    Ok(kmer_arrays(py, &Background::nrdb(&a).or_raise()?))
}

#[pyfunction]
pub fn background_glam2<'py>(
    py: Python<'py>,
    alphabet: &Bound<'py, PyAny>,
) -> PyResult<Vec<Floats<'py>>> {
    let a = conv::alphabet(alphabet)?;
    Ok(kmer_arrays(py, &profile::glam2_background(&a).or_raise()?))
}

#[pyfunction]
pub fn background_meme<'py>(
    py: Python<'py>,
    alphabet: &Bound<'py, PyAny>,
    text: &str,
) -> PyResult<Vec<Floats<'py>>> {
    let a = conv::alphabet(alphabet)?;
    Ok(kmer_arrays(
        py,
        &Background::from_meme(&a, text).or_raise()?,
    ))
}

#[pyfunction]
pub fn background_estimate<'py>(
    py: Python<'py>,
    alphabet: &Bound<'py, PyAny>,
    offsets: PyReadonlyArray1<'py, i64>,
    data: PyReadonlyArray1<'py, u8>,
    order: usize,
    pseudocount: f64,
) -> PyResult<Vec<Floats<'py>>> {
    let a = conv::alphabet(alphabet)?;
    let s = Strings {
        offsets: slice(&offsets, "offsets")?,
        data: slice(&data, "data")?,
    };
    let b = py
        .detach(|| Background::estimate(&a, s, order, pseudocount))
        .or_raise()?;
    Ok(kmer_arrays(py, &b))
}

#[pyfunction]
pub fn background_symmetrized<'py>(
    py: Python<'py>,
    kmers: Vec<PyReadonlyArray1<'py, f64>>,
    alphabet: &Bound<'py, PyAny>,
) -> PyResult<Vec<Floats<'py>>> {
    let a = conv::alphabet(alphabet)?;
    let complement = a
        .complement()
        .ok_or_else(|| invalid(format!("alphabet {:?} has no complement", a.name())))?;
    let b = kmers_background(a.len(), kmers)?;
    Ok(kmer_arrays(py, &b.symmetrized(complement).or_raise()?))
}

#[pyfunction]
pub fn background_sample<'py>(
    py: Python<'py>,
    kmers: Vec<PyReadonlyArray1<'py, f64>>,
    alphabet: &Bound<'py, PyAny>,
    lengths: PyReadonlyArray1<'py, u64>,
    seed: u64,
) -> PyResult<StringColumns<'py>> {
    let a = conv::alphabet(alphabet)?;
    let b = kmers_background(a.len(), kmers)?;
    let lengths: Vec<usize> = slice(&lengths, "lengths")?
        .iter()
        .map(|&l| usize::try_from(l).map_err(|_| invalid("sequence length exceeds address space")))
        .collect::<PyResult<_>>()?;
    lengths
        .iter()
        .try_fold(0usize, |total, &n| total.checked_add(n))
        .filter(|&total| total <= isize::MAX as usize && total <= i64::MAX as usize)
        .ok_or_else(|| invalid("sequence lengths exceed the offset or address space"))?;
    let (offsets, data) = py.detach(|| b.sample(&a, &lengths, seed));
    Ok((offsets.into_pyarray(py), data.into_pyarray(py)))
}

// ---------------------------------------------------------------------------
// Scoring policies

fn log_scale(s: &str) -> PyResult<LogScale> {
    Ok(match s {
        "bits" => LogScale::Bits,
        "meme-bits" => LogScale::MemeBits,
        "natural" => LogScale::Natural,
        "centibits" => LogScale::Centibits,
        _ => {
            return Err(invalid(format!(
                "log must be bits, meme-bits, natural or centibits, not {s:?}"
            )));
        }
    })
}

fn log_name(l: LogScale) -> &'static str {
    match l {
        LogScale::Bits => "bits",
        LogScale::MemeBits => "meme-bits",
        LogScale::Natural => "natural",
        LogScale::Centibits => "centibits",
    }
}

fn rule(s: &str) -> PyResult<AmbiguityRule> {
    Ok(match s {
        "skip" => AmbiguityRule::Skip,
        "mean" => AmbiguityRule::Mean,
        "sum" => AmbiguityRule::Sum,
        _ => {
            return Err(invalid(format!(
                "ambiguity rule must be skip, mean or sum, not {s:?}"
            )));
        }
    })
}

fn rule_name(r: AmbiguityRule) -> &'static str {
    match r {
        AmbiguityRule::Skip => "skip",
        AmbiguityRule::Mean => "mean",
        AmbiguityRule::Sum => "sum",
    }
}

fn ambiguity(o: &Bound<'_, PyAny>, alphabet: &Alphabet) -> PyResult<Ambiguity> {
    let mut codes = Vec::new();
    for item in o.getattr("codes")?.call_method0("items")?.try_iter()? {
        let (code, members): (String, String) = item?.extract()?;
        let mut set = bytes_of(&members, "ambiguity members")?
            .into_iter()
            .map(|b| {
                alphabet
                    .symbols()
                    .iter()
                    .position(|&s| s == b)
                    .map(|c| c as u8)
                    .ok_or_else(|| invalid(format!("{:?} is not an alphabet symbol", b as char)))
            })
            .collect::<PyResult<Vec<u8>>>()?;
        set.sort_unstable();
        set.dedup();
        codes.push((byte_of(&code, "an ambiguity code")?, set));
    }
    let unknown = attr::<Option<String>>(o, "unknown")?
        .map(|u| byte_of(&u, "the unknown code"))
        .transpose()?;
    Ok(Ambiguity {
        rule: rule(&attr::<String>(o, "rule")?)?,
        codes,
        unknown,
    })
}

type AmbiguityFields = (&'static str, Vec<(String, String)>, Option<String>);

fn ambiguity_fields(a: &Ambiguity, alphabet: Option<&Alphabet>) -> AmbiguityFields {
    let codes = a
        .codes
        .iter()
        .map(|(b, set)| {
            let members = match alphabet {
                Some(al) => set
                    .iter()
                    .map(|&c| al.symbols()[c as usize] as char)
                    .collect(),
                None => String::new(),
            };
            ((*b as char).to_string(), members)
        })
        .collect();
    (
        rule_name(a.rule),
        codes,
        a.unknown.map(|u| (u as char).to_string()),
    )
}

fn scoring(o: &Bound<'_, PyAny>, alphabet: &Alphabet, errors: &Errors) -> PyResult<Scoring> {
    let q = o.getattr("quantizer")?;
    let size: u32 = attr(&q, "size")?;
    let quantizer = match attr::<String>(&q, "kind")?.as_str() {
        "fimo" => Quantizer::Fimo { range: size },
        "precision" => Quantizer::Precision { states: size },
        "mast" => Quantizer::Mast { range: size },
        k => {
            return Err(invalid(format!(
                "quantizer must be fimo, precision or mast, not {k:?}"
            )));
        }
    };
    let rounding = match attr::<Option<u32>>(o, "rounding")? {
        Some(d) => Rounding::Significant(d),
        None => Rounding::None,
    };
    let background = match attr::<String>(o, "background")?.as_str() {
        "as-given" => BackgroundRule::AsGiven,
        "meme" => BackgroundRule::Meme,
        b => {
            return Err(invalid(format!(
                "background must be as-given or meme, not {b:?}"
            )));
        }
    };
    let letters = o.getattr("null_letters")?;
    let prior = o.getattr("column_prior")?;
    let column_prior: Option<Arc<dyn ColumnPrior>> = if prior.is_none() {
        None
    } else if kind(&prior).as_deref() == Some("substitution") {
        Some(Arc::new(SubstitutionPrior {
            alpha: attr(&prior, "alpha")?,
            beta: attr(&prior, "beta")?,
            target: floats(&prior.getattr("target")?.call_method0("ravel")?)?,
        }))
    } else {
        Some(Arc::new(Callback::new(&prior, errors, "column_prior")?))
    };
    Ok(Scoring {
        pseudocount: attr(o, "pseudocount")?,
        default_sites: attr(o, "default_sites")?,
        truncate_sites: attr(o, "truncate_sites")?,
        normalize_tolerance: attr(o, "normalize_tolerance")?,
        log: log_scale(&attr::<String>(o, "log")?)?,
        quantizer,
        rounding,
        background,
        null_budget: attr(o, "null_budget")?,
        null_letters: if letters.is_none() {
            None
        } else {
            Some(floats(&letters)?)
        },
        shared_reverse_null: attr(o, "shared_reverse_null")?,
        column_prior,
        ambiguity: ambiguity(&o.getattr("ambiguity")?, alphabet)?,
    })
}

fn scoring_dict<'py>(
    py: Python<'py>,
    s: &Scoring,
    alphabet: Option<&Alphabet>,
    prior: Option<&SubstitutionPrior>,
) -> PyResult<Bound<'py, PyDict>> {
    let d = PyDict::new(py);
    d.set_item("pseudocount", s.pseudocount)?;
    d.set_item("default_sites", s.default_sites)?;
    d.set_item("truncate_sites", s.truncate_sites)?;
    d.set_item("normalize_tolerance", s.normalize_tolerance)?;
    d.set_item("log", log_name(s.log))?;
    let q = match s.quantizer {
        Quantizer::Fimo { range } => ("fimo", range),
        Quantizer::Precision { states } => ("precision", states),
        Quantizer::Mast { range } => ("mast", range),
    };
    d.set_item("quantizer", q)?;
    d.set_item(
        "rounding",
        match s.rounding {
            Rounding::Significant(n) => Some(n),
            Rounding::None => None,
        },
    )?;
    d.set_item(
        "background",
        match s.background {
            BackgroundRule::AsGiven => "as-given",
            BackgroundRule::Meme => "meme",
        },
    )?;
    d.set_item("null_budget", s.null_budget)?;
    d.set_item(
        "null_letters",
        s.null_letters.as_ref().map(|l| l.clone().into_pyarray(py)),
    )?;
    d.set_item("shared_reverse_null", s.shared_reverse_null)?;
    d.set_item("column_prior", prior.map(|p| prior_fields(py, p)))?;
    d.set_item("ambiguity", ambiguity_fields(&s.ambiguity, alphabet))?;
    Ok(d)
}

fn prior_fields<'py>(py: Python<'py>, p: &SubstitutionPrior) -> (f64, f64, Bound<'py, PyAny>) {
    let n = (p.target.len() as f64).sqrt() as usize;
    let target = p.target.clone().into_pyarray(py);
    let target = target
        .reshape([n, n])
        .map(Bound::into_any)
        .unwrap_or_else(|_| p.target.clone().into_pyarray(py).into_any());
    (p.alpha, p.beta, target)
}

#[pyfunction]
pub fn scoring_preset<'py>(py: Python<'py>, name: &str) -> PyResult<Bound<'py, PyDict>> {
    let s = match name {
        "fimo" => Scoring::fimo(),
        "exact" => Scoring::exact(),
        _ => return Err(invalid(format!("no scoring preset {name:?}"))),
    };
    scoring_dict(py, &s, None, None)
}

fn strand_mode(s: &str) -> PyResult<StrandMode> {
    Ok(match s {
        "combine" => StrandMode::Combine,
        "separate" => StrandMode::Separate,
        "forward" => StrandMode::Forward,
        _ => {
            return Err(invalid(format!(
                "mode must be combine, separate or forward, not {s:?}"
            )));
        }
    })
}

#[pyfunction]
pub fn scoring_mast<'py>(
    py: Python<'py>,
    alphabet: &Bound<'py, PyAny>,
    letters: Option<PyReadonlyArray1<'py, f64>>,
    mode: &str,
) -> PyResult<Bound<'py, PyDict>> {
    let a = conv::alphabet(alphabet)?;
    let given = letters.as_ref().map(|l| slice(l, "letters")).transpose()?;
    let s = mast::scoring(&a, given, strand_mode(mode)?).or_raise()?;
    scoring_dict(py, &s, Some(&a), None)
}

#[pyfunction]
pub fn scoring_mcast<'py>(
    py: Python<'py>,
    kmers: Vec<PyReadonlyArray1<'py, f64>>,
    alphabet: &Bound<'py, PyAny>,
) -> PyResult<Bound<'py, PyDict>> {
    let a = conv::alphabet(alphabet)?;
    let b = kmers_background(a.len(), kmers)?;
    let s = mcast::scoring(&b, &a).or_raise()?;
    let prepared = s.prepare_background(&b, &a, Strands::Both).or_raise()?;
    let prior = SubstitutionPrior::mcast_dna(prepared.letters()).or_raise()?;
    scoring_dict(py, &s, Some(&a), Some(&prior))
}

#[pyfunction]
pub fn ambiguity_iupac(alphabet: &Bound<'_, PyAny>, rule_: &str) -> PyResult<AmbiguityFields> {
    let a = conv::alphabet(alphabet)?;
    let amb = Ambiguity::iupac(&a, rule(rule_)?).or_raise()?;
    Ok(ambiguity_fields(&amb, Some(&a)))
}

#[pyfunction]
pub fn substitution_mcast<'py>(
    py: Python<'py>,
    letters: PyReadonlyArray1<'py, f64>,
) -> PyResult<(f64, f64, Bound<'py, PyAny>)> {
    let p = SubstitutionPrior::mcast_dna(slice(&letters, "letters")?).or_raise()?;
    Ok(prior_fields(py, &p))
}

/// MEME's significant-digit rounding of `values`.
#[pyfunction]
pub fn round_significant<'py>(
    py: Python<'py>,
    values: PyReadonlyArray1<'py, f64>,
    digits: u32,
) -> PyResult<Floats<'py>> {
    let r = Rounding::Significant(digits);
    let out: Vec<f64> = slice(&values, "values")?
        .iter()
        .map(|&x| r.apply(x))
        .collect();
    Ok(out.into_pyarray(py))
}

// ---------------------------------------------------------------------------
// Reductions

fn threshold(o: &Bound<'_, PyAny>) -> PyResult<Threshold> {
    let (k, v): (String, f64) = attr(o, "threshold")?;
    Ok(match k.as_str() {
        "pvalue" => Threshold::PValue(v),
        _ => Threshold::Score(v),
    })
}

/// One reduction of a scan, as its native reducer.
enum AnyReducer {
    Hits(Hits),
    Best(SequenceBest),
    Counts(SequenceCounts),
    Scores(Dense),
    Mast(mast::Mast),
    Clusters(Clusters),
}

enum AnyOutput {
    Hits(Vec<HitEntry>),
    Best(Vec<BestEntry>),
    Counts(CountMatrix),
    Scores(DenseScores),
    Mast(MastReport),
    Clusters(Vec<Cluster>),
}

macro_rules! each {
    ($self:expr, $r:ident => $e:expr) => {
        match $self {
            AnyReducer::Hits($r) => $e,
            AnyReducer::Best($r) => $e,
            AnyReducer::Counts($r) => $e,
            AnyReducer::Scores($r) => $e,
            AnyReducer::Mast($r) => $e,
            AnyReducer::Clusters($r) => $e,
        }
    };
}

impl Reducer for AnyReducer {
    type Output = AnyOutput;

    fn needs(&self, out: &mut Vec<Demand>) {
        each!(self, r => r.needs(out))
    }

    fn prepare(&mut self, program: &Program, info: &ScanInfo<'_>) -> Result<()> {
        each!(self, r => r.prepare(program, info))
    }

    fn fork(&self) -> Self {
        match self {
            AnyReducer::Hits(r) => AnyReducer::Hits(r.fork()),
            AnyReducer::Best(r) => AnyReducer::Best(r.fork()),
            AnyReducer::Counts(r) => AnyReducer::Counts(r.fork()),
            AnyReducer::Scores(r) => AnyReducer::Scores(r.fork()),
            AnyReducer::Mast(r) => AnyReducer::Mast(r.fork()),
            AnyReducer::Clusters(r) => AnyReducer::Clusters(r.fork()),
        }
    }

    fn windows(&mut self, batch: &Windows<'_>) {
        each!(self, r => r.windows(batch))
    }

    fn sequence(&mut self, program: &Program, chunk: &SequenceChunk<'_>) {
        each!(self, r => r.sequence(program, chunk))
    }

    fn merge(&mut self, other: Self) {
        match (self, other) {
            (AnyReducer::Hits(a), AnyReducer::Hits(b)) => a.merge(b),
            (AnyReducer::Best(a), AnyReducer::Best(b)) => a.merge(b),
            (AnyReducer::Counts(a), AnyReducer::Counts(b)) => a.merge(b),
            (AnyReducer::Scores(a), AnyReducer::Scores(b)) => a.merge(b),
            (AnyReducer::Mast(a), AnyReducer::Mast(b)) => a.merge(b),
            (AnyReducer::Clusters(a), AnyReducer::Clusters(b)) => a.merge(b),
            _ => unreachable!("forks keep their variant"),
        }
    }

    fn finish(self, program: &Program, info: &ScanInfo<'_>) -> AnyOutput {
        match self {
            AnyReducer::Hits(r) => AnyOutput::Hits(r.finish(program, info)),
            AnyReducer::Best(r) => AnyOutput::Best(r.finish(program, info)),
            AnyReducer::Counts(r) => AnyOutput::Counts(r.finish(program, info)),
            AnyReducer::Scores(r) => AnyOutput::Scores(r.finish(program, info)),
            AnyReducer::Mast(r) => AnyOutput::Mast(r.finish(program, info)),
            AnyReducer::Clusters(r) => AnyOutput::Clusters(r.finish(program, info)),
        }
    }
}

/// A reduction: its reducer and, for clusters, what scores them afterwards.
struct Spec {
    reducer: AnyReducer,
    clusters: Option<(ClusterOptions, Option<Statistics>)>,
}

fn spec(o: &Bound<'_, PyAny>, program: &Program, errors: &Errors) -> PyResult<Spec> {
    let part = || -> PyResult<Option<usize>> {
        let p: Option<usize> = attr(o, "part")?;
        if let Some(p) = p
            && p >= program.parts().len()
        {
            return Err(PyIndexError::new_err(format!(
                "part {p} of a program with {} parts",
                program.parts().len()
            )));
        }
        Ok(p)
    };
    let kind = kind(o).ok_or_else(|| PyTypeError::new_err("expected a reduction"))?;
    let mut clusters = None;
    let reducer = match kind.as_str() {
        "hits" => {
            let r = Hits::new(threshold(o)?);
            AnyReducer::Hits(match part()? {
                Some(p) => r.in_part(p),
                None => r,
            })
        }
        "best" => {
            let r = SequenceBest::new();
            AnyReducer::Best(match part()? {
                Some(p) => r.in_part(p),
                None => r,
            })
        }
        "counts" => {
            let r = SequenceCounts::new(threshold(o)?);
            AnyReducer::Counts(match part()? {
                Some(p) => r.in_part(p),
                None => r,
            })
        }
        "scores" => AnyReducer::Scores(Dense::new()),
        "mast" => {
            let options = MastOptions {
                mode: strand_mode(&attr::<String>(o, "mode")?)?,
                hit_threshold: attr(o, "hit_threshold")?,
                evalue_threshold: attr(o, "evalue_threshold")?,
                best_of: best_of(&o.getattr("best_of")?, errors)?,
                combine: combine(&o.getattr("combine")?, errors)?,
                adjust: Arc::from(adjust(&o.getattr("adjust")?, errors)?),
            };
            let r = mast::Mast::new(options);
            AnyReducer::Mast(match part()? {
                Some(p) => r.in_part(p),
                None => r,
            })
        }
        "clusters" => {
            let options = ClusterOptions {
                pthresh: attr(o, "pthresh")?,
                max_gap: attr(o, "max_gap")?,
                egcost: attr(o, "egcost")?,
                threshold: attr(o, "threshold")?,
                segment: attr(o, "segment")?,
                overlap: attr(o, "overlap")?,
                part: part()?,
            };
            let stats = o.getattr("statistics")?;
            let stats = if stats.is_none() {
                None
            } else {
                Some(statistics(&stats, program.alphabet(), errors)?)
            };
            clusters = Some((options.clone(), stats));
            AnyReducer::Clusters(Clusters::new(options))
        }
        k => {
            return Err(PyTypeError::new_err(format!(
                "unknown reduction kind {k:?}"
            )));
        }
    };
    Ok(Spec { reducer, clusters })
}

fn hit_columns<'py>(py: Python<'py>, d: &Bound<'py, PyDict>, hits: &[HitEntry]) -> PyResult<()> {
    d.set_item(
        "sequence",
        hits.iter()
            .map(|h| h.sequence)
            .collect::<Vec<_>>()
            .into_pyarray(py),
    )?;
    d.set_item(
        "motif",
        hits.iter()
            .map(|h| h.motif)
            .collect::<Vec<_>>()
            .into_pyarray(py),
    )?;
    d.set_item(
        "start",
        hits.iter()
            .map(|h| h.start)
            .collect::<Vec<_>>()
            .into_pyarray(py),
    )?;
    d.set_item(
        "end",
        hits.iter()
            .map(|h| h.end)
            .collect::<Vec<_>>()
            .into_pyarray(py),
    )?;
    d.set_item(
        "strand",
        hits.iter()
            .map(|h| strand_code(h.strand))
            .collect::<Vec<_>>()
            .into_pyarray(py),
    )?;
    d.set_item(
        "score",
        hits.iter()
            .map(|h| h.score)
            .collect::<Vec<_>>()
            .into_pyarray(py),
    )?;
    d.set_item(
        "pvalue",
        hits.iter()
            .map(|h| h.pvalue)
            .collect::<Vec<_>>()
            .into_pyarray(py),
    )?;
    Ok(())
}

fn mast_dict<'py>(py: Python<'py>, d: &Bound<'py, PyDict>, r: &MastReport) -> PyResult<()> {
    let rows = &r.sequences;
    d.set_item(
        "sequence",
        rows.iter()
            .map(|s| s.sequence)
            .collect::<Vec<_>>()
            .into_pyarray(py),
    )?;
    d.set_item(
        "strand",
        rows.iter()
            .map(|s| s.strand.map_or(0, strand_code))
            .collect::<Vec<_>>()
            .into_pyarray(py),
    )?;
    d.set_item(
        "pvalue",
        rows.iter()
            .map(|s| s.pvalue)
            .collect::<Vec<_>>()
            .into_pyarray(py),
    )?;
    d.set_item(
        "evalue",
        rows.iter()
            .map(|s| s.evalue)
            .collect::<Vec<_>>()
            .into_pyarray(py),
    )?;
    d.set_item(
        "ranked",
        r.ranked
            .iter()
            .map(|&i| i as u64)
            .collect::<Vec<_>>()
            .into_pyarray(py),
    )?;
    let mut indptr = vec![0u64];
    let best: Vec<&MotifBest> = rows
        .iter()
        .flat_map(|s| {
            indptr.push(indptr.last().copied().unwrap_or(0) + s.best.len() as u64);
            s.best.iter()
        })
        .collect();
    let b = PyDict::new(py);
    b.set_item("indptr", indptr.into_pyarray(py))?;
    b.set_item(
        "motif",
        best.iter()
            .map(|m| m.motif)
            .collect::<Vec<_>>()
            .into_pyarray(py),
    )?;
    b.set_item(
        "score",
        best.iter()
            .map(|m| m.score)
            .collect::<Vec<_>>()
            .into_pyarray(py),
    )?;
    b.set_item(
        "pvalue",
        best.iter()
            .map(|m| m.pvalue)
            .collect::<Vec<_>>()
            .into_pyarray(py),
    )?;
    b.set_item(
        "sequence_pvalue",
        best.iter()
            .map(|m| m.sequence_pvalue)
            .collect::<Vec<_>>()
            .into_pyarray(py),
    )?;
    d.set_item("best", b)?;
    let mut indptr = vec![0u64];
    let hits: Vec<&MastHit> = rows
        .iter()
        .flat_map(|s| {
            indptr.push(indptr.last().copied().unwrap_or(0) + s.hits.len() as u64);
            s.hits.iter()
        })
        .collect();
    let h = PyDict::new(py);
    h.set_item("indptr", indptr.into_pyarray(py))?;
    h.set_item(
        "motif",
        hits.iter()
            .map(|m| m.motif)
            .collect::<Vec<_>>()
            .into_pyarray(py),
    )?;
    h.set_item(
        "start",
        hits.iter()
            .map(|m| m.start)
            .collect::<Vec<_>>()
            .into_pyarray(py),
    )?;
    h.set_item(
        "end",
        hits.iter()
            .map(|m| m.end)
            .collect::<Vec<_>>()
            .into_pyarray(py),
    )?;
    h.set_item(
        "strand",
        hits.iter()
            .map(|m| strand_code(m.strand))
            .collect::<Vec<_>>()
            .into_pyarray(py),
    )?;
    h.set_item(
        "score",
        hits.iter()
            .map(|m| m.score)
            .collect::<Vec<_>>()
            .into_pyarray(py),
    )?;
    h.set_item(
        "pvalue",
        hits.iter()
            .map(|m| m.pvalue)
            .collect::<Vec<_>>()
            .into_pyarray(py),
    )?;
    d.set_item("hits", h)?;
    Ok(())
}

fn cluster_dict<'py>(
    py: Python<'py>,
    d: &Bound<'py, PyDict>,
    r: &Report,
    threshold: f64,
) -> PyResult<()> {
    let c = &r.clusters;
    macro_rules! col {
        ($name:literal, $f:expr) => {
            d.set_item($name, c.iter().map($f).collect::<Vec<_>>().into_pyarray(py))?;
        };
    }
    col!("id", |s: &ScoredCluster| s.id);
    col!("sequence", |s: &ScoredCluster| s.cluster.sequence);
    col!("start", |s: &ScoredCluster| s.cluster.start);
    col!("end", |s: &ScoredCluster| s.cluster.end);
    col!("score", |s: &ScoredCluster| s.cluster.score);
    col!("excess", |s: &ScoredCluster| s.cluster.excess);
    col!("hits", |s: &ScoredCluster| s.cluster.hits);
    col!("gc", |s: &ScoredCluster| s.cluster.gc);
    col!("segment_start", |s: &ScoredCluster| s.cluster.segment_start);
    col!("segment_len", |s: &ScoredCluster| s.cluster.segment_len);
    col!("pvalue", |s: &ScoredCluster| s.pvalue);
    col!("evalue", |s: &ScoredCluster| s.evalue);
    col!("qvalue", |s: &ScoredCluster| s.qvalue);
    d.set_item("statistics", r.statistics)?;
    d.set_item("threshold", threshold)?;
    Ok(())
}

/// A finished reduction, with the clusters' report when it is one.
enum Done {
    Plain(AnyOutput),
    Clusters(Report, f64),
}

fn output_dict<'py>(
    py: Python<'py>,
    out: &Done,
    sequences: usize,
    motifs: usize,
) -> PyResult<Bound<'py, PyDict>> {
    let d = PyDict::new(py);
    match out {
        Done::Plain(AnyOutput::Hits(h)) => {
            d.set_item("kind", "hits")?;
            hit_columns(py, &d, h)?;
        }
        Done::Plain(AnyOutput::Best(b)) => {
            d.set_item("kind", "best")?;
            d.set_item(
                "sequence",
                b.iter()
                    .map(|h| h.sequence)
                    .collect::<Vec<_>>()
                    .into_pyarray(py),
            )?;
            d.set_item(
                "motif",
                b.iter()
                    .map(|h| h.motif)
                    .collect::<Vec<_>>()
                    .into_pyarray(py),
            )?;
            d.set_item(
                "start",
                b.iter()
                    .map(|h| h.start)
                    .collect::<Vec<_>>()
                    .into_pyarray(py),
            )?;
            d.set_item(
                "strand",
                b.iter()
                    .map(|h| strand_code(h.strand))
                    .collect::<Vec<_>>()
                    .into_pyarray(py),
            )?;
            d.set_item(
                "score",
                b.iter()
                    .map(|h| h.score)
                    .collect::<Vec<_>>()
                    .into_pyarray(py),
            )?;
            d.set_item(
                "pvalue",
                b.iter()
                    .map(|h| h.pvalue)
                    .collect::<Vec<_>>()
                    .into_pyarray(py),
            )?;
        }
        Done::Plain(AnyOutput::Counts(c)) => {
            d.set_item("kind", "counts")?;
            d.set_item("indptr", c.indptr.clone().into_pyarray(py))?;
            d.set_item("motifs", c.motifs.clone().into_pyarray(py))?;
            d.set_item("counts", c.counts.clone().into_pyarray(py))?;
            d.set_item("shape", (sequences, motifs))?;
        }
        Done::Plain(AnyOutput::Scores(s)) => {
            d.set_item("kind", "scores")?;
            d.set_item("offsets", s.offsets.clone().into_pyarray(py))?;
            d.set_item("values", s.values.clone().into_pyarray(py))?;
            d.set_item("motifs", s.motifs)?;
            d.set_item("strands", s.strands)?;
        }
        Done::Plain(AnyOutput::Mast(r)) => {
            d.set_item("kind", "mast")?;
            mast_dict(py, &d, r)?;
        }
        Done::Plain(AnyOutput::Clusters(_)) => unreachable!("clusters finish with their report"),
        Done::Clusters(r, threshold) => {
            d.set_item("kind", "clusters")?;
            cluster_dict(py, &d, r, *threshold)?;
        }
    }
    Ok(d)
}

// ---------------------------------------------------------------------------
// Programs

/// Compiled motifs.
#[pyclass(frozen, module = "omicsweave._native")]
pub struct MotifProgram {
    program: Program,
}

type TermFields<'py> = (Bound<'py, PyArray1<i32>>, Bound<'py, PyAny>);

#[pymethods]
impl MotifProgram {
    #[new]
    #[pyo3(signature = (source, alphabet, kmers, scoring, strands, threads))]
    fn new(
        py: Python<'_>,
        source: &Bound<'_, PyAny>,
        alphabet: &Bound<'_, PyAny>,
        kmers: Vec<PyReadonlyArray1<'_, f64>>,
        scoring: &Bound<'_, PyAny>,
        strands: &str,
        threads: Option<usize>,
    ) -> PyResult<MotifProgram> {
        let alphabet = conv::alphabet(alphabet)?;
        let set = motif_set(source, &alphabet)?;
        let background = kmers_background(alphabet.len(), kmers)?;
        let errors = Errors::default();
        let scoring = self::scoring(scoring, &alphabet, &errors)?;
        let strands = self::strands(strands)?;
        let program = py.detach(|| {
            Program::compile(
                &views(&set),
                &alphabet,
                &background,
                &scoring,
                strands,
                threads,
            )
        });
        errors.check()?;
        Ok(MotifProgram {
            program: program.or_raise()?,
        })
    }

    #[staticmethod]
    fn combine(programs: Vec<Bound<'_, MotifProgram>>) -> PyResult<MotifProgram> {
        let all: Vec<Program> = programs.iter().map(|p| p.get().program.clone()).collect();
        Ok(MotifProgram {
            program: Program::combine(all).or_raise()?,
        })
    }

    fn __len__(&self) -> usize {
        self.program.len()
    }

    fn widths<'py>(&self, py: Python<'py>) -> Bound<'py, PyArray1<u32>> {
        (0..self.program.len())
            .map(|i| self.program.motif(i).width())
            .collect::<Vec<_>>()
            .into_pyarray(py)
    }

    /// Motif `motif` on `strand`, for windows with the full background
    /// context: its real and integer terms, quantization and null.
    fn model<'py>(
        &self,
        py: Python<'py>,
        motif: usize,
        strand: &str,
    ) -> PyResult<Bound<'py, PyDict>> {
        if motif >= self.program.len() {
            return Err(PyIndexError::new_err("motif index out of range"));
        }
        let s = self::strand(strand)?;
        let m = self.program.motif(motif).strand(s).ok_or_else(|| {
            invalid(format!(
                "motif {motif} was not compiled for strand {strand}"
            ))
        })?;
        let v = m.full();
        let real = v.score_model();
        let int = v.int_model();
        let n = self.program.alphabet().len();
        let shaped = |table: Bound<'py, PyAny>, order: usize| -> PyResult<Bound<'py, PyAny>> {
            table.call_method1("reshape", (vec![n; order],))
        };
        let d = PyDict::new(py);
        d.set_item("width", real.width)?;
        let terms = real
            .terms
            .iter()
            .map(|t| -> PyResult<TermFields<'py>> {
                Ok((
                    t.positions.clone().into_pyarray(py),
                    shaped(
                        t.table.clone().into_pyarray(py).into_any(),
                        t.positions.len(),
                    )?,
                ))
            })
            .collect::<PyResult<Vec<_>>>()?;
        d.set_item("terms", terms)?;
        let int_terms = int
            .terms
            .iter()
            .map(|t| -> PyResult<TermFields<'py>> {
                Ok((
                    t.positions.clone().into_pyarray(py),
                    shaped(
                        t.table.clone().into_pyarray(py).into_any(),
                        t.positions.len(),
                    )?,
                ))
            })
            .collect::<PyResult<Vec<_>>>()?;
        d.set_item("int_terms", int_terms)?;
        d.set_item("scale", int.scale)?;
        d.set_item("offset", int.offset)?;
        d.set_item("max", int.max)?;
        d.set_item("pmf", v.null().pmf().to_vec().into_pyarray(py))?;
        d.set_item("survival", v.null().survival().to_vec().into_pyarray(py))?;
        d.set_item("forbidden", v.null().forbidden())?;
        Ok(d)
    }

    #[pyo3(signature = (offsets, data, reductions, threads, tile_len))]
    fn scan<'py>(
        &self,
        py: Python<'py>,
        offsets: PyReadonlyArray1<'py, i64>,
        data: PyReadonlyArray1<'py, u8>,
        reductions: Vec<Bound<'py, PyAny>>,
        threads: Option<usize>,
        tile_len: Option<usize>,
    ) -> PyResult<Vec<Bound<'py, PyDict>>> {
        let errors = Errors::default();
        let mut reducers = Vec::with_capacity(reductions.len());
        let mut after = Vec::with_capacity(reductions.len());
        for r in &reductions {
            let s = spec(r, &self.program, &errors)?;
            reducers.push(s.reducer);
            after.push(s.clusters);
        }
        let options = ScanOptions {
            threads,
            tile_len: tile_len.unwrap_or(ScanOptions::default().tile_len),
        };
        let sequences = Strings {
            offsets: slice(&offsets, "offsets")?,
            data: slice(&data, "data")?,
        };
        let program = &self.program;
        let done = py.detach(|| -> Result<Vec<Done>> {
            let outputs = program.scan(sequences, reducers, &options)?;
            outputs
                .into_iter()
                .zip(after)
                .map(|(out, after)| match (out, after) {
                    (AnyOutput::Clusters(found), Some((cluster_options, stats))) => {
                        let search = ClusterSearch {
                            program,
                            options: &cluster_options,
                            scan: &options,
                        };
                        let report = search.report(&found, stats.as_ref())?;
                        Ok(Done::Clusters(report, search.threshold()))
                    }
                    (out, _) => Ok(Done::Plain(out)),
                })
                .collect()
        });
        errors.check()?;
        let done = done.or_raise()?;
        let n = sequences.len();
        done.iter()
            .map(|d| output_dict(py, d, n, self.program.len()))
            .collect()
    }
}

// ---------------------------------------------------------------------------
// Gapped profiles

fn mixture_fields<'py>(
    py: Python<'py>,
    m: &DirichletMixture,
) -> PyResult<(Floats<'py>, Bound<'py, PyArray2<f64>>)> {
    Ok((
        m.weights.clone().into_pyarray(py),
        PyArray2::from_vec2(py, &m.alphas).map_err(|e| invalid(e.to_string()))?,
    ))
}

#[pyfunction]
pub fn dirichlet_glam2<'py>(
    py: Python<'py>,
    alphabet: &Bound<'py, PyAny>,
) -> PyResult<(Floats<'py>, Bound<'py, PyArray2<f64>>)> {
    mixture_fields(py, &DirichletMixture::glam2(&conv::alphabet(alphabet)?))
}

#[pyfunction]
pub fn dirichlet_uniform(
    py: Python<'_>,
    n: usize,
) -> PyResult<(Floats<'_>, Bound<'_, PyArray2<f64>>)> {
    mixture_fields(py, &DirichletMixture::uniform(n))
}

#[pyfunction]
pub fn dirichlet_parse_glam2<'py>(
    py: Python<'py>,
    text: &str,
    n: usize,
) -> PyResult<(Floats<'py>, Bound<'py, PyArray2<f64>>)> {
    mixture_fields(py, &DirichletMixture::parse_glam2(text, n).or_raise()?)
}

type MixtureFields<'py> = (Floats<'py>, Bound<'py, PyArray2<f64>>);

#[pyfunction]
pub fn profile_scoring_glam2<'py>(
    py: Python<'py>,
    alphabet: &Bound<'py, PyAny>,
) -> PyResult<(f64, f64, f64, f64, MixtureFields<'py>)> {
    let s = ProfileScoring::glam2(&conv::alphabet(alphabet)?).or_raise()?;
    Ok((
        s.delete_pseudo,
        s.no_delete_pseudo,
        s.insert_pseudo,
        s.no_insert_pseudo,
        mixture_fields(py, &s.prior)?,
    ))
}

fn profile_scoring(o: &Bound<'_, PyAny>) -> PyResult<ProfileScoring> {
    let prior = o.getattr("prior")?;
    let alphas: Vec<Vec<f64>> = prior.getattr("alphas")?.call_method0("tolist")?.extract()?;
    Ok(ProfileScoring {
        delete_pseudo: attr(o, "delete_pseudo")?,
        no_delete_pseudo: attr(o, "no_delete_pseudo")?,
        insert_pseudo: attr(o, "insert_pseudo")?,
        no_insert_pseudo: attr(o, "no_insert_pseudo")?,
        prior: DirichletMixture::new(floats(&prior.getattr("weights")?)?, alphas).or_raise()?,
    })
}

fn alignment_dict<'py>(py: Python<'py>, hits: &[Alignment]) -> PyResult<Bound<'py, PyDict>> {
    let d = PyDict::new(py);
    d.set_item(
        "sequence",
        hits.iter()
            .map(|h| h.sequence)
            .collect::<Vec<_>>()
            .into_pyarray(py),
    )?;
    d.set_item(
        "motif",
        hits.iter()
            .map(|h| h.motif)
            .collect::<Vec<_>>()
            .into_pyarray(py),
    )?;
    d.set_item(
        "strand",
        hits.iter()
            .map(|h| strand_code(h.strand))
            .collect::<Vec<_>>()
            .into_pyarray(py),
    )?;
    d.set_item(
        "start",
        hits.iter()
            .map(|h| h.start)
            .collect::<Vec<_>>()
            .into_pyarray(py),
    )?;
    d.set_item(
        "end",
        hits.iter()
            .map(|h| h.end)
            .collect::<Vec<_>>()
            .into_pyarray(py),
    )?;
    d.set_item(
        "score",
        hits.iter()
            .map(|h| h.score)
            .collect::<Vec<_>>()
            .into_pyarray(py),
    )?;
    let mut indptr = vec![0u64];
    let mut matched = Vec::new();
    let mut offset = Vec::new();
    for h in hits {
        for c in &h.columns {
            matched.push(c.matched);
            offset.push(c.offset);
        }
        indptr.push(matched.len() as u64);
    }
    d.set_item("column_indptr", indptr.into_pyarray(py))?;
    d.set_item("matched", matched.into_pyarray(py))?;
    d.set_item("offset", offset.into_pyarray(py))?;
    Ok(d)
}

fn alignments(o: &Bound<'_, PyAny>) -> PyResult<Vec<Alignment>> {
    let sequence: Vec<u32> = column(o, "sequence")?;
    let motif: Vec<u32> = column(o, "motif")?;
    let strand: Vec<i8> = column(o, "strand")?;
    let start: Vec<u32> = column(o, "start")?;
    let end: Vec<u32> = column(o, "end")?;
    let score: Vec<f64> = column(o, "score")?;
    let indptr: Vec<u64> = column(o, "column_indptr")?;
    let matched: Vec<bool> = column(o, "matched")?;
    let offset: Vec<u32> = column(o, "offset")?;
    let n = score.len();
    if [
        sequence.len(),
        motif.len(),
        strand.len(),
        start.len(),
        end.len(),
    ]
    .iter()
    .any(|&l| l != n)
        || indptr.len() != n + 1
        || indptr.last().copied() != Some(matched.len() as u64)
        || offset.len() != matched.len()
    {
        return Err(invalid("alignment columns differ in length"));
    }
    (0..n)
        .map(|i| {
            let (a, b) = (indptr[i] as usize, indptr[i + 1] as usize);
            if a > b {
                return Err(invalid("column_indptr decreases"));
            }
            Ok(Alignment {
                sequence: sequence[i],
                motif: motif[i],
                strand: strand_of(strand[i])?,
                start: start[i],
                end: end[i],
                score: score[i],
                columns: (a..b)
                    .map(|k| Column {
                        matched: matched[k],
                        offset: offset[k],
                    })
                    .collect(),
            })
        })
        .collect()
}

/// Compiled gapped profiles.
#[pyclass(frozen, module = "omicsweave._native")]
pub struct ProfileProgram {
    program: Profiles,
}

#[pymethods]
impl ProfileProgram {
    #[new]
    fn new(
        py: Python<'_>,
        source: &Bound<'_, PyAny>,
        alphabet: &Bound<'_, PyAny>,
        kmers: Vec<PyReadonlyArray1<'_, f64>>,
        scoring: &Bound<'_, PyAny>,
        strands: &str,
    ) -> PyResult<ProfileProgram> {
        let alphabet = conv::alphabet(alphabet)?;
        let set = motif_set(source, &alphabet)?;
        let background = kmers_background(alphabet.len(), kmers)?;
        let scoring = profile_scoring(scoring)?;
        let strands = self::strands(strands)?;
        let program = py
            .detach(|| Profiles::compile(&views(&set), &alphabet, &background, &scoring, strands))
            .or_raise()?;
        Ok(ProfileProgram { program })
    }

    fn __len__(&self) -> usize {
        self.program.len()
    }

    #[allow(clippy::too_many_arguments)]
    #[pyo3(signature = (offsets, data, top, names, limit, threads))]
    fn scan<'py>(
        &self,
        py: Python<'py>,
        offsets: PyReadonlyArray1<'py, i64>,
        data: PyReadonlyArray1<'py, u8>,
        top: Option<usize>,
        names: Option<StringArrays<'py>>,
        limit: Option<usize>,
        threads: Option<usize>,
    ) -> PyResult<Bound<'py, PyDict>> {
        let sequences = Strings {
            offsets: slice(&offsets, "offsets")?,
            data: slice(&data, "data")?,
        };
        let options = ScanOptions {
            threads,
            ..ScanOptions::default()
        };
        let program = &self.program;
        let hits = match top {
            Some(n) => {
                let names = names
                    .as_ref()
                    .ok_or_else(|| invalid("top alignments need sequence names"))?;
                let names = strings(names, "names")?;
                if names.len() != sequences.len() {
                    return Err(invalid("one name per sequence"));
                }
                py.detach(|| program.scan(sequences, TopAlignments::glam2(n, names), &options))
                    .or_raise()?
                    .into_iter()
                    .flatten()
                    .collect::<Vec<_>>()
            }
            None => py
                .detach(|| program.scan(sequences, AllAlignments::new(limit), &options))
                .or_raise()?,
        };
        alignment_dict(py, &hits)
    }
}

// ---------------------------------------------------------------------------
// Text outputs

fn text(bytes: Vec<u8>) -> String {
    conv::text_of(&bytes)
}

fn io(e: std::io::Error) -> PyErr {
    crate::error::to_py(&Error::from(e))
}

fn hit_entries(o: &Bound<'_, PyAny>) -> PyResult<Vec<HitEntry>> {
    let sequence: Vec<u32> = column(o, "sequence")?;
    let motif: Vec<u32> = column(o, "motif")?;
    let start: Vec<u32> = column(o, "start")?;
    let end: Vec<u32> = column(o, "end")?;
    let strand: Vec<i8> = column(o, "strand")?;
    let score: Vec<f64> = column(o, "score")?;
    let pvalue: Vec<f64> = column(o, "pvalue")?;
    let n = score.len();
    if [
        sequence.len(),
        motif.len(),
        start.len(),
        end.len(),
        strand.len(),
        pvalue.len(),
    ]
    .iter()
    .any(|&l| l != n)
    {
        return Err(invalid("hit columns differ in length"));
    }
    (0..n)
        .map(|i| {
            Ok(HitEntry {
                score: score[i],
                pvalue: pvalue[i],
                sequence: sequence[i],
                motif: motif[i],
                start: start[i],
                end: end[i],
                strand: strand_of(strand[i])?,
            })
        })
        .collect()
}

fn check_rows(
    hits_seq: impl Iterator<Item = u32>,
    motif: impl Iterator<Item = u32>,
    sequences: usize,
    motifs: usize,
) -> PyResult<()> {
    if hits_seq.into_iter().any(|s| s as usize >= sequences) {
        return Err(invalid("a sequence number is out of range"));
    }
    if motif.into_iter().any(|m| m as usize >= motifs) {
        return Err(invalid("a motif number is out of range"));
    }
    Ok(())
}

#[pyfunction]
#[allow(clippy::too_many_arguments)]
pub fn text_fimo(
    hits: &Bound<'_, PyAny>,
    alphabet: &Bound<'_, PyAny>,
    motif_ids: Vec<String>,
    alt_ids: Vec<String>,
    names: StringArrays<'_>,
    sequences: StringArrays<'_>,
    origins: Option<PyReadonlyArray1<'_, u64>>,
) -> PyResult<String> {
    let alphabet = conv::alphabet(alphabet)?;
    let hits = hit_entries(hits)?;
    let names = strings(&names, "names")?;
    let sequences = strings(&sequences, "sequences")?;
    let origins = origins.as_ref().map(|o| slice(o, "origins")).transpose()?;
    check_rows(
        hits.iter().map(|h| h.sequence),
        hits.iter().map(|h| h.motif),
        sequences.len().min(names.len()),
        motif_ids.len(),
    )?;
    if hits
        .iter()
        .any(|h| h.end as usize > sequences.get(h.sequence as usize).len())
    {
        return Err(invalid("a hit runs past its sequence"));
    }
    if origins.is_some_and(|o| o.len() < sequences.len()) {
        return Err(invalid("one origin per sequence"));
    }
    let ids: Vec<&str> = motif_ids.iter().map(String::as_str).collect();
    let alts: Vec<&str> = alt_ids.iter().map(String::as_str).collect();
    let mut out = Vec::new();
    fimo::write_text(
        &mut out,
        &hits,
        &fimo::FimoText {
            alphabet: &alphabet,
            motif_ids: &ids,
            motif_alt_ids: &alts,
            sequence_names: names,
            sequences,
            origins,
        },
    )
    .map_err(io)?;
    Ok(text(out))
}

fn mast_report(o: &Bound<'_, PyAny>) -> PyResult<MastReport> {
    let sequence: Vec<u32> = column(o, "sequence")?;
    let strand: Vec<i8> = column(o, "strand")?;
    let pvalue: Vec<f64> = column(o, "pvalue")?;
    let evalue: Vec<f64> = column(o, "evalue")?;
    let ranked: Vec<u64> = column(o, "ranked")?;
    let b = o.getattr("best")?;
    let b_ptr: Vec<u64> = column(&b, "indptr")?;
    let b_motif: Vec<u32> = column(&b, "motif")?;
    let b_score: Vec<f64> = column(&b, "score")?;
    let b_p: Vec<f64> = column(&b, "pvalue")?;
    let b_sp: Vec<f64> = column(&b, "sequence_pvalue")?;
    let h = o.getattr("hits")?;
    let h_ptr: Vec<u64> = column(&h, "indptr")?;
    let h_motif: Vec<u32> = column(&h, "motif")?;
    let h_start: Vec<u32> = column(&h, "start")?;
    let h_end: Vec<u32> = column(&h, "end")?;
    let h_strand: Vec<i8> = column(&h, "strand")?;
    let h_score: Vec<f64> = column(&h, "score")?;
    let h_p: Vec<f64> = column(&h, "pvalue")?;
    let n = sequence.len();
    let csr_ok = |ptr: &[u64], len: usize| {
        ptr.len() == n + 1
            && ptr.windows(2).all(|w| w[0] <= w[1])
            && ptr.last() == Some(&(len as u64))
    };
    if [strand.len(), pvalue.len(), evalue.len()]
        .iter()
        .any(|&l| l != n)
        || !csr_ok(&b_ptr, b_motif.len())
        || !csr_ok(&h_ptr, h_motif.len())
        || [b_score.len(), b_p.len(), b_sp.len()]
            .iter()
            .any(|&l| l != b_motif.len())
        || [
            h_start.len(),
            h_end.len(),
            h_strand.len(),
            h_score.len(),
            h_p.len(),
        ]
        .iter()
        .any(|&l| l != h_motif.len())
        || ranked.iter().any(|&r| r as usize >= n)
    {
        return Err(invalid("MAST report columns are inconsistent"));
    }
    let sequences = (0..n)
        .map(|i| {
            let best = (b_ptr[i] as usize..b_ptr[i + 1] as usize)
                .map(|k| MotifBest {
                    motif: b_motif[k],
                    score: b_score[k],
                    pvalue: b_p[k],
                    sequence_pvalue: b_sp[k],
                })
                .collect();
            let hits = (h_ptr[i] as usize..h_ptr[i + 1] as usize)
                .map(|k| {
                    Ok(MastHit {
                        motif: h_motif[k],
                        start: h_start[k],
                        end: h_end[k],
                        strand: strand_of(h_strand[k])?,
                        score: h_score[k],
                        pvalue: h_p[k],
                    })
                })
                .collect::<PyResult<Vec<_>>>()?;
            Ok(MastSequence {
                sequence: sequence[i],
                strand: if strand[i] == 0 {
                    None
                } else {
                    Some(strand_of(strand[i])?)
                },
                pvalue: pvalue[i],
                evalue: evalue[i],
                best,
                hits,
            })
        })
        .collect::<PyResult<Vec<_>>>()?;
    Ok(MastReport {
        sequences,
        ranked: ranked.into_iter().map(|r| r as usize).collect(),
    })
}

#[pyfunction]
pub fn text_mast_hit_list(
    report: &Bound<'_, PyAny>,
    motif_ids: Vec<String>,
    alt_ids: Vec<String>,
    numbers: Vec<u32>,
    names: StringArrays<'_>,
    stranded: bool,
) -> PyResult<String> {
    let report = mast_report(report)?;
    let names = strings(&names, "names")?;
    if motif_ids.len() != alt_ids.len() || motif_ids.len() != numbers.len() {
        return Err(invalid("one ID, alternate ID and number per motif"));
    }
    check_rows(
        report.sequences.iter().map(|s| s.sequence),
        report
            .sequences
            .iter()
            .flat_map(|s| s.hits.iter().map(|h| h.motif)),
        names.len(),
        motif_ids.len(),
    )?;
    let ids: Vec<&str> = motif_ids.iter().map(String::as_str).collect();
    let alts: Vec<&str> = alt_ids.iter().map(String::as_str).collect();
    let mut out = Vec::new();
    mast::write_hit_list(
        &mut out,
        &report,
        &mast::MastText {
            sequence_names: names,
            motif_ids: &ids,
            motif_alt_ids: &alts,
            motif_numbers: &numbers,
            stranded,
        },
    )
    .map_err(io)?;
    Ok(text(out))
}

fn cluster_report(o: &Bound<'_, PyAny>) -> PyResult<Report> {
    let id: Vec<u32> = column(o, "id")?;
    let sequence: Vec<u32> = column(o, "sequence")?;
    let start: Vec<u32> = column(o, "start")?;
    let end: Vec<u32> = column(o, "end")?;
    let score: Vec<f64> = column(o, "score")?;
    let excess: Vec<f64> = column(o, "excess")?;
    let hits: Vec<u32> = column(o, "hits")?;
    let gc: Vec<f64> = column(o, "gc")?;
    let segment_start: Vec<u32> = column(o, "segment_start")?;
    let segment_len: Vec<u32> = column(o, "segment_len")?;
    let pvalue: Vec<f64> = column(o, "pvalue")?;
    let evalue: Vec<f64> = column(o, "evalue")?;
    let qvalue: Vec<f64> = column(o, "qvalue")?;
    let n = id.len();
    if [
        sequence.len(),
        start.len(),
        end.len(),
        score.len(),
        excess.len(),
        hits.len(),
        gc.len(),
        segment_start.len(),
        segment_len.len(),
        pvalue.len(),
        evalue.len(),
        qvalue.len(),
    ]
    .iter()
    .any(|&l| l != n)
    {
        return Err(invalid("cluster columns differ in length"));
    }
    Ok(Report {
        clusters: (0..n)
            .map(|i| ScoredCluster {
                cluster: Cluster {
                    sequence: sequence[i],
                    start: start[i],
                    end: end[i],
                    score: score[i],
                    excess: excess[i],
                    hits: hits[i],
                    gc: gc[i],
                    segment_start: segment_start[i],
                    segment_len: segment_len[i],
                },
                id: id[i],
                pvalue: pvalue[i],
                evalue: evalue[i],
                qvalue: qvalue[i],
            })
            .collect(),
        statistics: attr(o, "statistics")?,
    })
}

#[pyfunction]
pub fn text_mcast(
    clusters: &Bound<'_, PyAny>,
    names: StringArrays<'_>,
    sequences: StringArrays<'_>,
    origins: Option<PyReadonlyArray1<'_, u64>>,
    filter: (Option<f64>, Option<f64>, Option<f64>),
) -> PyResult<String> {
    let report = cluster_report(clusters)?;
    let names = strings(&names, "names")?;
    let sequences = strings(&sequences, "sequences")?;
    let origins = origins.as_ref().map(|o| slice(o, "origins")).transpose()?;
    let n = sequences.len().min(names.len());
    if report.clusters.iter().any(|c| {
        c.cluster.sequence as usize >= n
            || c.cluster.end as usize > sequences.get(c.cluster.sequence as usize).len()
    }) {
        return Err(invalid("a cluster lies outside the sequences"));
    }
    if origins.is_some_and(|o| o.len() < sequences.len()) {
        return Err(invalid("one origin per sequence"));
    }
    let d = OutputFilter::default();
    let mut out = Vec::new();
    mcast::write_text(
        &mut out,
        &report,
        &mcast::McastText {
            sequence_names: names,
            sequences,
            origins,
            filter: OutputFilter {
                evalue: filter.0.unwrap_or(d.evalue),
                pvalue: filter.1.unwrap_or(d.pvalue),
                qvalue: filter.2.unwrap_or(d.qvalue),
            },
        },
    )
    .map_err(io)?;
    Ok(text(out))
}

#[pyfunction]
pub fn text_glam2(
    alignments: &Bound<'_, PyAny>,
    alphabet: &Bound<'_, PyAny>,
    names: StringArrays<'_>,
    sequences: StringArrays<'_>,
    header: Option<(String, String)>,
) -> PyResult<String> {
    let alphabet = conv::alphabet(alphabet)?;
    let hits = self::alignments(alignments)?;
    let names = strings(&names, "names")?;
    let sequences = strings(&sequences, "sequences")?;
    let n = sequences.len().min(names.len());
    if hits.iter().any(|h| {
        h.sequence as usize >= n || h.end as usize > sequences.get(h.sequence as usize).len()
    }) {
        return Err(invalid("an alignment lies outside the sequences"));
    }
    let mut out = Vec::new();
    if let Some((version, command)) = header {
        profile::write_glam2_header(&mut out, &version, &command).map_err(io)?;
    }
    profile::write_glam2_text(
        &mut out,
        &hits,
        &Glam2Text::new(&alphabet, names, sequences),
    )
    .map_err(io)?;
    Ok(text(out))
}

/// Registers the motif classes and functions.
pub fn register(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_class::<conv::MotifCollection>()?;
    m.add_class::<conv::FfmWriter>()?;
    m.add_class::<MotifProgram>()?;
    m.add_class::<ProfileProgram>()?;
    m.add_function(wrap_pyfunction!(conv::alphabet_preset, m)?)?;
    for f in [
        wrap_pyfunction!(motif_defaults, m)?,
        wrap_pyfunction!(stats_best_of, m)?,
        wrap_pyfunction!(stats_combine, m)?,
        wrap_pyfunction!(stats_adjust, m)?,
        wrap_pyfunction!(background_nrdb, m)?,
        wrap_pyfunction!(background_glam2, m)?,
        wrap_pyfunction!(background_meme, m)?,
        wrap_pyfunction!(background_estimate, m)?,
        wrap_pyfunction!(background_symmetrized, m)?,
        wrap_pyfunction!(background_sample, m)?,
        wrap_pyfunction!(scoring_preset, m)?,
        wrap_pyfunction!(scoring_mast, m)?,
        wrap_pyfunction!(scoring_mcast, m)?,
        wrap_pyfunction!(ambiguity_iupac, m)?,
        wrap_pyfunction!(substitution_mcast, m)?,
        wrap_pyfunction!(round_significant, m)?,
        wrap_pyfunction!(dirichlet_glam2, m)?,
        wrap_pyfunction!(dirichlet_uniform, m)?,
        wrap_pyfunction!(dirichlet_parse_glam2, m)?,
        wrap_pyfunction!(profile_scoring_glam2, m)?,
        wrap_pyfunction!(text_fimo, m)?,
        wrap_pyfunction!(text_mast_hit_list, m)?,
        wrap_pyfunction!(text_mcast, m)?,
        wrap_pyfunction!(text_glam2, m)?,
    ] {
        m.add_function(f)?;
    }
    Ok(())
}
