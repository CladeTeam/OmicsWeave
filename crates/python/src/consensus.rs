//! Consensus sequences (`bcftools consensus`) of a VCF over a reference.

use std::fs::File;
use std::io::BufReader;
use std::path::PathBuf;
use std::sync::Arc;

use numpy::{IntoPyArray, PyArray1};
use omicsweave_consensus::{
    Alleles, Consensus, ConsensusOptions, Mark, Mask, MaskWith, Overlap, RecordFilter, VcfStore,
};
use omicsweave_core::{Error, Region, Result};
use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;

use crate::error::{Done, Fail, OrRaise};
use crate::sequence::SequenceReader;
use crate::shared::Shared;

/// One query's choice of alleles, as the Python selector classes encode it:
/// (`all`, [], None), (`samples`, names, None) or
/// (`haplotype`, [sample], haplotype).
type Selector = (String, Vec<String>, Option<String>);

struct Builder {
    consensus: Consensus,
    /// The region `consensus` holds loaded.
    loaded: Option<(String, u32, u32)>,
}

impl Builder {
    fn build(
        &mut self,
        reference: &Shared<crate::sequence::Sequences>,
        region: Region<'_>,
        selector: &Selector,
    ) -> Done<usize> {
        let cached = self.loaded.as_ref().is_some_and(|(chrom, start, end)| {
            chrom == region.chrom && *start == region.start && *end == region.end
        });
        if !cached {
            self.loaded = None;
            let mut guard = reference.lock();
            let reader = guard.as_mut().ok_or(Fail::Closed)?;
            self.consensus.load(reader, region)?;
            self.loaded = Some((region.chrom.to_owned(), region.start, region.end));
        }
        let names: Vec<&str>;
        let alleles = match (selector.0.as_str(), selector.1.as_slice(), &selector.2) {
            ("all", [], None) => Alleles::AllAlt,
            ("samples", samples, None) => {
                names = samples.iter().map(String::as_str).collect();
                Alleles::Samples(&names)
            }
            ("haplotype", [sample], Some(haplotype)) => Alleles::Haplotype { sample, haplotype },
            _ => {
                return Err(Fail::Py(PyValueError::new_err(format!(
                    "bad allele selector {selector:?}"
                ))));
            }
        };
        Ok(self.consensus.query(alleles)?)
    }
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
    records: usize,
    skipped: usize,
}

#[pymethods]
impl ConsensusBuilder {
    #[new]
    #[allow(clippy::too_many_arguments)]
    #[pyo3(signature = (vcf, reference, iupac_codes, missing, absent, mark_del, mark_ins, mark_snv, mask, mask_with, regions_overlap, chain, include, exclude))]
    fn new(
        py: Python<'_>,
        vcf: PathBuf,
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
        let built = py.detach(|| -> Result<(VcfStore, Option<Mask>)> {
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
            let store = VcfStore::read(BufReader::new(File::open(&vcf)?), filter)?;
            let mask = match &mask {
                Some(p) => Some(Mask::parse(&std::fs::read(p)?, mask_with)?),
                None => None,
            };
            Ok((store, mask))
        });
        let (store, mask) = built.or_raise()?;
        options.mask = mask.map(Arc::new);
        let samples = store.samples().to_vec();
        let (records, skipped) = (store.n_records(), store.skipped_records());
        Ok(ConsensusBuilder {
            inner: Shared::new(Builder {
                consensus: Consensus::new(Arc::new(store), options),
                loaded: None,
            }),
            reference,
            samples,
            records,
            skipped,
        })
    }

    fn samples(&self) -> Vec<String> {
        self.samples.clone()
    }

    /// (records kept, malformed lines skipped).
    fn records(&self) -> (usize, usize) {
        (self.records, self.skipped)
    }

    /// The consensus of `chrom:[start, end)` with `selector`: its bases and
    /// its chain (with `chain`).
    fn sequence<'py>(
        &self,
        py: Python<'py>,
        chrom: &str,
        start: u32,
        end: u32,
        selector: Selector,
    ) -> PyResult<(Bound<'py, PyArray1<u8>>, Option<String>)> {
        let reference = &self.reference.get().inner;
        let (bases, chain) = self.inner.run(py, |b| -> Done<_> {
            let n = b.build(reference, Region::new(chrom, start, end), &selector)?;
            let mut bases = vec![0u8; n];
            b.consensus.write(&mut bases)?;
            Ok((bases, b.consensus.chain().map(str::to_owned)))
        })?;
        Ok((bases.into_pyarray(py), chain))
    }

    /// The consensus of `chrom:[start, end)` with each selector, packed:
    /// sequence `i` is `bases[offsets[i]..offsets[i + 1]]`.
    #[allow(clippy::type_complexity)]
    fn sequences<'py>(
        &self,
        py: Python<'py>,
        chrom: &str,
        start: u32,
        end: u32,
        selectors: Vec<Selector>,
    ) -> PyResult<(
        Bound<'py, PyArray1<i64>>,
        Bound<'py, PyArray1<u8>>,
        Vec<Option<String>>,
    )> {
        let reference = &self.reference.get().inner;
        let (offsets, bases, chains) = self.inner.run(py, |b| -> Done<_> {
            let region = Region::new(chrom, start, end);
            let mut offsets = Vec::with_capacity(selectors.len() + 1);
            offsets.push(0i64);
            let mut bases = Vec::new();
            let mut chains = Vec::with_capacity(selectors.len());
            for selector in &selectors {
                let n = b.build(reference, region, selector)?;
                let at = bases.len();
                bases.resize(at + n, 0u8);
                b.consensus.write(&mut bases[at..])?;
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
