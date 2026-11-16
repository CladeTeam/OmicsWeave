//! Variants (BCF + CSI): the records of a region as columns, and their
//! genotype matrix.

use std::collections::HashMap;
use std::path::{Path, PathBuf};
use std::sync::Arc;

use numpy::{IntoPyArray, PyArray1};
use omicsweave_bgzf::{BinIndex, locate_index};
use omicsweave_core::{Error, FileSource, Region, Result};
use omicsweave_vcf::{BcfReader, Column, RecordFilter, RecordQuery, Strings, Values};
use pyo3::IntoPyObjectExt;
use pyo3::prelude::*;
use pyo3::types::{PyDict, PyList};

use crate::cache::{self, PyBlockCache};
use crate::error::{Done, OrRaise};
use crate::shared::{Shared, path_str, threads_or_available};

/// The CSI at `given`, or else `<bcf>.csi` as htslib finds it; it must
/// exist (indexes are built explicitly, never on open).
pub(crate) fn open_index(bcf: &Path, given: Option<&Path>) -> Result<(BinIndex, PathBuf)> {
    let path = match given {
        Some(p) => p.to_path_buf(),
        None => locate_index(bcf).ok_or_else(|| Error::Io {
            kind: std::io::ErrorKind::NotFound,
            detail: format!(
                "no index of {} found (.csi beside it); build it with Bcf.build_index()",
                bcf.display()
            ),
        })?,
    };
    Ok((BinIndex::read_file(&path)?, path))
}

type PyStrings<'py> = (Bound<'py, PyArray1<i64>>, Bound<'py, PyArray1<u8>>);

fn py_strings(py: Python<'_>, s: Strings) -> PyStrings<'_> {
    (s.offsets.into_pyarray(py), s.data.into_pyarray(py))
}

/// (start, end, id, ref, alt, ploidy, gt, phased).
type Matrix<'py> = (
    Bound<'py, PyArray1<i64>>,
    Bound<'py, PyArray1<i64>>,
    PyStrings<'py>,
    PyStrings<'py>,
    PyStrings<'py>,
    usize,
    Bound<'py, PyArray1<i8>>,
    Bound<'py, PyArray1<bool>>,
);

/// `(name, kind, width, values)`: kind `flag`, `int`, `float` or `str`
/// (values then `(offsets, data)`).
fn py_column<'py>(py: Python<'py>, c: Column) -> PyResult<Bound<'py, PyAny>> {
    let (kind, width, values) = match c.values {
        Values::Flag(v) => ("flag", 1, v.into_pyarray(py).into_any()),
        Values::Int { width, values } => ("int", width, values.into_pyarray(py).into_any()),
        Values::Float { width, values } => ("float", width, values.into_pyarray(py).into_any()),
        Values::Str(s) => ("str", 1, py_strings(py, s).into_bound_py_any(py)?),
    };
    (c.name, kind, width, values).into_bound_py_any(py)
}

/// Builds the `.csi` of a sorted BGZF BCF beside it, as `bcftools index`
/// does; returns its path.
#[pyfunction]
#[pyo3(signature = (path, min_shift=omicsweave_vcf::CSI_MIN_SHIFT, threads=None))]
pub fn build_bcf_index(
    py: Python<'_>,
    path: PathBuf,
    min_shift: u32,
    threads: Option<usize>,
) -> PyResult<String> {
    let threads = threads_or_available(threads)?;
    py.detach(|| omicsweave_vcf::build_index_file(&path, min_shift, threads))
        .map(|at| path_str(&at))
        .or_raise()
}

#[pyclass(frozen, module = "omicsweave._native")]
pub struct VariantReader {
    inner: Shared<BcfReader<FileSource>>,
    chroms: Arc<[(String, Option<u64>)]>,
    samples: Arc<[String]>,
    all_sample_ids: Arc<[usize]>,
    index: PathBuf,
}

#[pymethods]
impl VariantReader {
    #[new]
    #[pyo3(signature = (path, index=None, cache=None))]
    fn new(
        py: Python<'_>,
        path: PathBuf,
        index: Option<PathBuf>,
        cache: Option<PyRef<'_, PyBlockCache>>,
    ) -> PyResult<VariantReader> {
        let cache = cache::handle(cache)?;
        py.detach(|| -> Result<VariantReader> {
            let (index, at) = open_index(&path, index.as_deref())?;
            let mut reader = BcfReader::open(FileSource::open(&path)?, index)?;
            reader.set_cache(cache);
            let header = reader.header();
            let mut first_sample = HashMap::with_capacity(header.samples.len());
            let all_sample_ids = header
                .samples
                .iter()
                .enumerate()
                .map(|(i, name)| *first_sample.entry(name).or_insert(i))
                .collect();
            Ok(VariantReader {
                chroms: header
                    .contigs
                    .iter()
                    .filter(|c| !c.name.is_empty())
                    .map(|c| (c.name.clone(), c.len))
                    .collect(),
                samples: header.samples.clone().into(),
                all_sample_ids,
                index: at,
                inner: Shared::new(reader),
            })
        })
        .or_raise()
    }

    /// (name, length or None) of every contig the header declares.
    fn chroms(&self) -> Vec<(String, Option<u64>)> {
        self.chroms.to_vec()
    }

    fn samples(&self) -> Vec<String> {
        self.samples.to_vec()
    }

    #[getter]
    fn index_path(&self) -> String {
        path_str(&self.index)
    }

    fn set_cache(&self, py: Python<'_>, cache: Option<PyRef<'_, PyBlockCache>>) -> PyResult<()> {
        let cache = cache::handle(cache)?;
        self.inner.run(py, |r| {
            r.set_cache(cache);
            Ok(())
        })
    }

    fn clone_context(&self, py: Python<'_>) -> PyResult<Self> {
        let reader = self.inner.run(py, |r| Ok(r.try_clone()?))?;
        Ok(Self {
            inner: Shared::new(reader),
            chroms: self.chroms.clone(),
            samples: self.samples.clone(),
            all_sample_ids: Arc::clone(&self.all_sample_ids),
            index: self.index.clone(),
        })
    }

    /// The genotype matrix of `samples` (header indexes) over the records
    /// overlapping `chrom:[start, end)`.
    #[pyo3(signature = (chrom, start, end, samples=None))]
    fn genotypes<'py>(
        &self,
        py: Python<'py>,
        chrom: &str,
        start: u32,
        end: u32,
        samples: Option<Vec<usize>>,
    ) -> PyResult<Matrix<'py>> {
        let samples = samples.as_deref().unwrap_or(&self.all_sample_ids);
        let m = self.inner.run(py, |r| -> Done<_> {
            Ok(r.genotypes(Region::new(chrom, start, end), samples)?)
        })?;
        Ok((
            m.start.into_pyarray(py),
            m.end.into_pyarray(py),
            py_strings(py, m.id),
            py_strings(py, m.reference),
            py_strings(py, m.alt),
            m.ploidy,
            m.gt.into_pyarray(py),
            m.phased.into_pyarray(py),
        ))
    }

    /// The records overlapping `chrom:[start, end)` that `include` keeps
    /// (or `exclude` leaves out; bcftools expressions), with the genotype
    /// matrix of `samples`, QUAL, FILTER and the `info` and `format`
    /// fields: a dict of columns, fields as `(name, kind, width, values)`.
    #[pyo3(signature = (chrom, start, end, samples=None, info=Vec::new(), format=Vec::new(),
                        include=None, exclude=None))]
    #[allow(clippy::too_many_arguments)]
    fn records<'py>(
        &self,
        py: Python<'py>,
        chrom: &str,
        start: u32,
        end: u32,
        samples: Option<Vec<usize>>,
        info: Vec<String>,
        format: Vec<String>,
        include: Option<String>,
        exclude: Option<String>,
    ) -> PyResult<Bound<'py, PyDict>> {
        let samples = samples.as_deref().unwrap_or(&self.all_sample_ids);
        let query = RecordQuery { info, format };
        let t = self.inner.run(py, |r| -> Done<_> {
            let mut filter = match (&include, &exclude) {
                (Some(_), Some(_)) => {
                    return Err(
                        Error::invalid("variants", "give include or exclude, not both").into(),
                    );
                }
                (Some(e), None) => Some(RecordFilter::include(r.header(), e)?),
                (None, Some(e)) => Some(RecordFilter::exclude(r.header(), e)?),
                (None, None) => None,
            };
            Ok(r.records(
                Region::new(chrom, start, end),
                samples,
                &query,
                filter.as_mut(),
            )?)
        })?;
        let m = t.genotypes;
        let d = PyDict::new(py);
        d.set_item("start", m.start.into_pyarray(py))?;
        d.set_item("end", m.end.into_pyarray(py))?;
        d.set_item("id", py_strings(py, m.id))?;
        d.set_item("ref", py_strings(py, m.reference))?;
        d.set_item("alt", py_strings(py, m.alt))?;
        d.set_item("ploidy", m.ploidy)?;
        d.set_item("gt", m.gt.into_pyarray(py))?;
        d.set_item("phased", m.phased.into_pyarray(py))?;
        d.set_item("qual", t.qual.into_pyarray(py))?;
        d.set_item("filters", py_strings(py, t.filters))?;
        let columns = |cs: Vec<Column>| -> PyResult<Bound<'py, PyList>> {
            let items = cs
                .into_iter()
                .map(|c| py_column(py, c))
                .collect::<PyResult<Vec<_>>>()?;
            PyList::new(py, items)
        };
        d.set_item("info", columns(t.info)?)?;
        d.set_item("format", columns(t.format)?)?;
        Ok(d)
    }

    fn close(&self, py: Python<'_>) {
        drop(self.inner.take(py));
    }
}
