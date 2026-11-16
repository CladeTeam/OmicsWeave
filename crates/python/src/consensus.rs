//! Consensus sequences (`bcftools consensus`) of a BCF, read by region
//! through its CSI, over a reference.

use std::path::PathBuf;
use std::sync::Arc;

use numpy::{IntoPyArray, PyArray1};
use omicsweave_consensus::{
    Alleles, Consensus, ConsensusOptions, Mark, Mask, MaskWith, Overlap, RecordFilter,
};
use omicsweave_core::{Error, FileSource, Region, Result, Strand};
use omicsweave_vcf::BcfReader;
use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;
use pyo3::pybacked::PyBackedStr;

use crate::error::{Done, Fail, OrRaise};
use crate::sequence::SequenceReader;
use crate::shared::{Shared, path_str};
use crate::variant::open_index;

/// One query's choice of alleles, as the Python selector classes encode it:
/// (`all`, [], None), (`samples`, names, None) or
/// (`haplotype`, [sample], haplotype).
type Selector = (PyBackedStr, Vec<PyBackedStr>, Option<PyBackedStr>);

struct Builder {
    consensus: Consensus,
}

impl Builder {
    fn build(
        &mut self,
        reference: &Shared<crate::sequence::Sequences>,
        region: Region<'_>,
        selector: &Selector,
        load: bool,
    ) -> Done<usize> {
        let names: Vec<&str>;
        let mut local_names = [""; 4];
        let alleles = match (&*selector.0, selector.1.as_slice(), &selector.2) {
            ("all", [], None) => Alleles::AllAlt,
            ("samples", samples, None) => {
                if samples.len() <= local_names.len() {
                    for (out, name) in local_names.iter_mut().zip(samples) {
                        *out = name;
                    }
                    Alleles::Samples(&local_names[..samples.len()])
                } else {
                    names = samples.iter().map(|s| &**s).collect();
                    Alleles::Samples(&names)
                }
            }
            ("haplotype", [sample], Some(haplotype)) => Alleles::Haplotype { sample, haplotype },
            _ => {
                return Err(Fail::Py(PyValueError::new_err(format!(
                    "bad allele selector {selector:?}"
                ))));
            }
        };
        if load {
            reference.ensure_open()?;
            let mut guard = reference.lock();
            let reader = guard.as_mut().ok_or(Fail::Closed)?;
            self.consensus.load(reader, region)?;
        }
        Ok(self.consensus.query(alleles)?)
    }
}

fn strand(s: &str) -> PyResult<Strand> {
    match s.as_bytes() {
        [b] => Strand::from_byte(*b),
        _ => None,
    }
    .ok_or_else(|| PyValueError::new_err(format!("strand {s:?} is not one of + - . ?")))
}

fn byte(what: &str, c: Option<char>) -> PyResult<Option<u8>> {
    c.map(|c| {
        u8::try_from(c)
            .ok()
            .filter(u8::is_ascii)
            .ok_or_else(|| PyValueError::new_err(format!("{what} must be one ASCII character")))
    })
    .transpose()
}

fn mark(what: &str, m: Option<&str>) -> PyResult<Option<Mark>> {
    Ok(match m {
        None => None,
        Some("upper") => Some(Mark::Upper),
        Some("lower") => Some(Mark::Lower),
        Some(s) => {
            let mut chars = s.chars();
            match (chars.next(), chars.next()) {
                (Some(c), None) => Some(Mark::Char(byte(what, Some(c))?.unwrap())),
                _ => {
                    return Err(PyValueError::new_err(format!(
                        "{what} must be 'upper', 'lower' or one character"
                    )));
                }
            }
        }
    })
}

#[pyclass(frozen, module = "omicsweave._native")]
pub struct ConsensusBuilder {
    inner: Shared<Builder>,
    reference: Py<SequenceReader>,
    samples: Vec<String>,
    index: PathBuf,
}

#[pymethods]
impl ConsensusBuilder {
    #[new]
    #[allow(clippy::too_many_arguments)]
    #[pyo3(signature = (bcf, index, reference, iupac_codes, missing, absent, mark_del, mark_ins, mark_snv, mask, mask_with, regions_overlap, chain, include, exclude))]
    fn new(
        py: Python<'_>,
        bcf: PathBuf,
        index: Option<PathBuf>,
        reference: Py<SequenceReader>,
        iupac_codes: bool,
        missing: Option<char>,
        absent: Option<char>,
        mark_del: Option<char>,
        mark_ins: Option<&str>,
        mark_snv: Option<&str>,
        mask: Option<PathBuf>,
        mask_with: &str,
        regions_overlap: &str,
        chain: bool,
        include: Option<String>,
        exclude: Option<String>,
    ) -> PyResult<ConsensusBuilder> {
        let mask_with = match mask_with {
            "upper" => MaskWith::Upper,
            "lower" => MaskWith::Lower,
            s => match byte("mask_with", s.chars().next())? {
                Some(c) if s.len() == 1 => MaskWith::Char(c),
                _ => {
                    return Err(PyValueError::new_err(
                        "mask_with must be 'upper', 'lower' or one character",
                    ));
                }
            },
        };
        let regions_overlap = match regions_overlap {
            "start" => Overlap::Start,
            "record" => Overlap::Record,
            "variant" => Overlap::Variant,
            other => {
                return Err(PyValueError::new_err(format!(
                    "regions_overlap must be 'start', 'record' or 'variant', not {other:?}"
                )));
            }
        };
        let mut options = ConsensusOptions {
            iupac_codes,
            missing: byte("missing", missing)?,
            absent: byte("absent", absent)?,
            mark_del: byte("mark_del", mark_del)?,
            mark_ins: mark("mark_ins", mark_ins)?,
            mark_snv: mark("mark_snv", mark_snv)?,
            mask: None,
            regions_overlap,
            chain,
        };
        let built = py.detach(|| -> Result<(Consensus, PathBuf)> {
            let filter = match (&include, &exclude) {
                (Some(_), Some(_)) => {
                    return Err(Error::invalid(
                        "consensus",
                        "both include and exclude given",
                    ));
                }
                (Some(e), None) => Some(RecordFilter::Include(e)),
                (None, Some(e)) => Some(RecordFilter::Exclude(e)),
                (None, None) => None,
            };
            if let Some(p) = &mask {
                options.mask = Some(Arc::new(Mask::parse(&std::fs::read(p)?, mask_with)?));
            }
            let (csi, at) = open_index(&bcf, index.as_deref())?;
            let variants = BcfReader::open(FileSource::open(&bcf)?, csi)?;
            Ok((Consensus::new(variants, filter, options)?, at))
        });
        let (consensus, index) = built.or_raise()?;
        Ok(ConsensusBuilder {
            samples: consensus.samples().to_vec(),
            inner: Shared::new(Builder { consensus }),
            reference,
            index,
        })
    }

    fn samples(&self) -> Vec<String> {
        self.samples.clone()
    }

    #[getter]
    fn index_path(&self) -> String {
        path_str(&self.index)
    }

    fn clone_context(&self, py: Python<'_>, reference: Py<SequenceReader>) -> PyResult<Self> {
        reference.get().inner.ensure_process()?;
        let consensus = self.inner.run(py, |b| Ok(b.consensus.try_clone()?))?;
        Ok(Self {
            inner: Shared::new(Builder { consensus }),
            reference,
            samples: self.samples.clone(),
            index: self.index.clone(),
        })
    }

    /// The consensus of `chrom:[start, end)` with `selector`, 5' to 3' of
    /// `strand`: its bases and its chain (with `chain`).
    #[allow(clippy::too_many_arguments)]
    fn sequence<'py>(
        &self,
        py: Python<'py>,
        chrom: &str,
        start: u32,
        end: u32,
        strand: &str,
        selector: Selector,
    ) -> PyResult<(Bound<'py, PyArray1<u8>>, Option<String>)> {
        let strand = self::strand(strand)?;
        let reference = &self.reference.get().inner;
        let (bases, chain) = self.inner.run(py, |b| -> Done<_> {
            let n = b.build(reference, Region::new(chrom, start, end), &selector, true)?;
            let mut bases = vec![0u8; n];
            b.consensus.write(strand, &mut bases)?;
            Ok((bases, b.consensus.chain().map(str::to_owned)))
        })?;
        Ok((bases.into_pyarray(py), chain))
    }

    /// The consensus of `chrom:[start, end)` with each selector, 5' to 3'
    /// of `strand`, packed: sequence `i` is `bases[offsets[i]..offsets[i + 1]]`.
    #[allow(clippy::type_complexity, clippy::too_many_arguments)]
    fn sequences<'py>(
        &self,
        py: Python<'py>,
        chrom: &str,
        start: u32,
        end: u32,
        strand: &str,
        selectors: Vec<Selector>,
    ) -> PyResult<(
        Bound<'py, PyArray1<i64>>,
        Bound<'py, PyArray1<u8>>,
        Vec<Option<String>>,
    )> {
        let strand = self::strand(strand)?;
        let reference = &self.reference.get().inner;
        let (offsets, bases, chains) = self.inner.run(py, |b| -> Done<_> {
            let region = Region::new(chrom, start, end);
            let mut offsets = Vec::with_capacity(selectors.len() + 1);
            offsets.push(0i64);
            let mut bases = Vec::new();
            let mut chains = Vec::with_capacity(selectors.len());
            for (index, selector) in selectors.iter().enumerate() {
                // Only selectors in this call share the loaded region.
                let n = b.build(reference, region, selector, index == 0)?;
                let at = bases.len();
                bases.resize(at + n, 0u8);
                b.consensus.write(strand, &mut bases[at..])?;
                offsets.push(bases.len() as i64);
                chains.push(b.consensus.chain().map(str::to_owned));
            }
            Ok((offsets, bases, chains))
        })?;
        Ok((offsets.into_pyarray(py), bases.into_pyarray(py), chains))
    }

    fn close(&self, py: Python<'_>) {
        drop(self.inner.take(py));
    }
}
