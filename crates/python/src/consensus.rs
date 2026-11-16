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

use crate::cache::{self, PyBlockCache};
use crate::error::{Done, Fail, OrRaise};
use crate::sequence::SequenceReader;
use crate::shared::{Shared, path_str};
use crate::variant::open_index;

/// One query's choice of alleles, read from a Python selector class
/// (`types.AllAlt`, `Samples`, `Haplotype`, `Record`, `Variant`) by its
/// `kind` and its attributes.
enum Selection {
    AllAlt,
    Samples(Vec<PyBackedStr>),
    Haplotype {
        sample: PyBackedStr,
        haplotype: PyBackedStr,
    },
    Record {
        index: usize,
        alt: usize,
    },
    Variant {
        start: i64,
        reference: PyBackedStr,
        alt: PyBackedStr,
    },
}

impl Selection {
    fn read(o: &Bound<'_, PyAny>) -> PyResult<Selection> {
        fn attr<'py, T: FromPyObjectOwned<'py>>(o: &Bound<'py, PyAny>, name: &str) -> PyResult<T> {
            o.getattr(name)?.extract::<T>().map_err(Into::into)
        }
        let kind: PyBackedStr = attr(o, "kind")?;
        Ok(match &*kind {
            "all" => Selection::AllAlt,
            "samples" => Selection::Samples(attr(o, "names")?),
            "haplotype" => Selection::Haplotype {
                sample: attr(o, "sample")?,
                haplotype: attr(o, "haplotype")?,
            },
            "record" => Selection::Record {
                index: attr(o, "index")?,
                alt: attr(o, "alt")?,
            },
            "variant" => Selection::Variant {
                start: attr(o, "start")?,
                reference: attr(o, "ref")?,
                alt: attr(o, "alt")?,
            },
            other => {
                return Err(PyValueError::new_err(format!(
                    "unknown allele selection {other:?}"
                )));
            }
        })
    }
}

struct Builder {
    consensus: Consensus,
}

impl Builder {
    fn load(
        &mut self,
        reference: &Shared<crate::sequence::Sequences>,
        region: Region<'_>,
    ) -> Done<()> {
        reference.ensure_open()?;
        let mut guard = reference.lock();
        let reader = guard.as_mut().ok_or(Fail::Closed)?;
        Ok(self.consensus.load(reader, region)?)
    }

    fn build(&mut self, selection: &Selection) -> Done<usize> {
        let names: Vec<&str>;
        let mut local_names = [""; 4];
        let alleles = match selection {
            Selection::AllAlt => Alleles::AllAlt,
            Selection::Samples(samples) => {
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
            Selection::Haplotype { sample, haplotype } => Alleles::Haplotype { sample, haplotype },
            Selection::Record { index, alt } => Alleles::Record {
                index: *index,
                alt: *alt,
            },
            Selection::Variant {
                start,
                reference,
                alt,
            } => Alleles::Variant {
                start: *start,
                reference: reference.as_bytes(),
                alt: alt.as_bytes(),
            },
        };
        Ok(self.consensus.query(alleles)?)
    }
}

/// Packs strings as Arrow offsets and bytes.
fn pack<'s>(strings: impl Iterator<Item = &'s [u8]>) -> (Vec<i64>, Vec<u8>) {
    let mut offsets = vec![0i64];
    let mut data = Vec::new();
    for s in strings {
        data.extend_from_slice(s);
        offsets.push(data.len() as i64);
    }
    (offsets, data)
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
    /// The consensus of the BCF at `bcf` (through its CSI) over
    /// `reference`, shaped by `options` (a `query.ConsensusOptions`, read
    /// attribute by attribute); `cache` holds the BCF's decoded blocks.
    #[new]
    #[pyo3(signature = (bcf, index, reference, options, cache=None))]
    fn new(
        py: Python<'_>,
        bcf: PathBuf,
        index: Option<PathBuf>,
        reference: Py<SequenceReader>,
        options: &Bound<'_, PyAny>,
        cache: Option<PyRef<'_, PyBlockCache>>,
    ) -> PyResult<ConsensusBuilder> {
        let cache = cache::handle(cache)?;
        fn attr<'py, T: FromPyObjectOwned<'py>>(o: &Bound<'py, PyAny>, name: &str) -> PyResult<T> {
            o.getattr(name)?.extract::<T>().map_err(Into::into)
        }
        let o = options;
        let iupac_codes: bool = attr(o, "iupac_codes")?;
        let missing: Option<char> = attr(o, "missing")?;
        let absent: Option<char> = attr(o, "absent")?;
        let mark_del: Option<char> = attr(o, "mark_del")?;
        let mark_ins: Option<String> = attr(o, "mark_ins")?;
        let mark_snv: Option<String> = attr(o, "mark_snv")?;
        let mask: Option<PathBuf> = attr(o, "mask")?;
        let mask_with: String = attr(o, "mask_with")?;
        let regions_overlap: String = attr(o, "regions_overlap")?;
        let chain: bool = attr(o, "chain")?;
        let include: Option<String> = attr(o, "include")?;
        let exclude: Option<String> = attr(o, "exclude")?;
        let mask_with = match mask_with.as_str() {
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
        let regions_overlap = match regions_overlap.as_str() {
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
            mark_ins: mark("mark_ins", mark_ins.as_deref())?,
            mark_snv: mark("mark_snv", mark_snv.as_deref())?,
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
            let mut variants = BcfReader::open(FileSource::open(&bcf)?, csi)?;
            variants.set_cache(cache);
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

    fn set_cache(&self, py: Python<'_>, cache: Option<PyRef<'_, PyBlockCache>>) -> PyResult<()> {
        let cache = cache::handle(cache)?;
        self.inner.run(py, |b| {
            b.consensus.variants_mut().set_cache(cache);
            Ok(())
        })
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
        selector: &Bound<'py, PyAny>,
    ) -> PyResult<(Bound<'py, PyArray1<u8>>, Option<String>)> {
        let strand = self::strand(strand)?;
        let selection = Selection::read(selector)?;
        let reference = &self.reference.get().inner;
        let (bases, chain) = self.inner.run(py, |b| -> Done<_> {
            b.load(reference, Region::new(chrom, start, end))?;
            let n = b.build(&selection)?;
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
        selectors: Vec<Bound<'py, PyAny>>,
    ) -> PyResult<(
        Bound<'py, PyArray1<i64>>,
        Bound<'py, PyArray1<u8>>,
        Vec<Option<String>>,
    )> {
        let strand = self::strand(strand)?;
        let selections = selectors
            .iter()
            .map(Selection::read)
            .collect::<PyResult<Vec<_>>>()?;
        let reference = &self.reference.get().inner;
        let (offsets, bases, chains) = self.inner.run(py, |b| -> Done<_> {
            let mut offsets = Vec::with_capacity(selections.len() + 1);
            offsets.push(0i64);
            let mut bases = Vec::new();
            let mut chains = Vec::with_capacity(selections.len());
            // Only selections in this call share the loaded region.
            if !selections.is_empty() {
                b.load(reference, Region::new(chrom, start, end))?;
            }
            for selection in &selections {
                let n = b.build(selection)?;
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

    /// The records `chrom:[start, end)` takes, in the order a `Record`
    /// selection counts them: starts, ends (start plus the reference span),
    /// REF and the comma-joined ALTs as Arrow strings, and whether the mask
    /// leaves each out.
    #[allow(clippy::type_complexity)]
    fn records<'py>(
        &self,
        py: Python<'py>,
        chrom: &str,
        start: u32,
        end: u32,
    ) -> PyResult<(
        Bound<'py, PyArray1<i64>>,
        Bound<'py, PyArray1<i64>>,
        (Bound<'py, PyArray1<i64>>, Bound<'py, PyArray1<u8>>),
        (Bound<'py, PyArray1<i64>>, Bound<'py, PyArray1<u8>>),
        Bound<'py, PyArray1<bool>>,
    )> {
        let reference = &self.reference.get().inner;
        let (starts, ends, refs, alts, masked) = self.inner.run(py, |b| -> Done<_> {
            b.load(reference, Region::new(chrom, start, end))?;
            let records: Vec<_> = b.consensus.records()?.collect();
            let starts = records.iter().map(|r| r.start()).collect::<Vec<_>>();
            let ends = records
                .iter()
                .map(|r| r.start() + r.len())
                .collect::<Vec<_>>();
            let refs = pack(records.iter().map(|r| r.allele(0)));
            let joined: Vec<Vec<u8>> = records
                .iter()
                .map(|r| {
                    (1..r.n_alleles())
                        .map(|i| r.allele(i))
                        .collect::<Vec<_>>()
                        .join(&b","[..])
                })
                .collect();
            let alts = pack(joined.iter().map(Vec::as_slice));
            let masked = records.iter().map(|r| r.masked).collect::<Vec<_>>();
            Ok((starts, ends, refs, alts, masked))
        })?;
        Ok((
            starts.into_pyarray(py),
            ends.into_pyarray(py),
            (refs.0.into_pyarray(py), refs.1.into_pyarray(py)),
            (alts.0.into_pyarray(py), alts.1.into_pyarray(py)),
            masked.into_pyarray(py),
        ))
    }

    fn close(&self, py: Python<'_>) {
        drop(self.inner.take(py));
    }
}

/// Whether `text` is a haplotype `bcftools consensus -H` takes.
#[pyfunction]
pub fn is_haplotype(text: &str) -> bool {
    omicsweave_consensus::is_haplotype(text)
}
