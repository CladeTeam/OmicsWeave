//! Variants (BCF + CSI): the genotype matrix of a region.

use std::collections::HashMap;
use std::path::{Path, PathBuf};
use std::sync::Arc;

use numpy::{IntoPyArray, PyArray1};
use omicsweave_bgzf::{BinIndex, locate_index};
use omicsweave_core::{Error, FileSource, Region, Result};
use omicsweave_vcf::{BcfReader, Strings};
use pyo3::prelude::*;
use pyo3::types::PyDict;

use crate::cache::{self, CacheKey, CacheNamespace, PyBlockCache};
use crate::error::{Done, OrRaise};
use crate::shared::{Shared, path_str, seek_stats};

/// The CSI at `given`, or else `<bcf>.csi` as htslib finds it; it must
/// exist (indexes are built explicitly, never on open).
pub(crate) fn open_index(bcf: &Path, given: Option<&Path>) -> Result<(BinIndex, PathBuf)> {
    let path = match given {
        Some(p) => p.to_path_buf(),
        None => locate_index(bcf).ok_or_else(|| Error::Io {
            kind: std::io::ErrorKind::NotFound,
            detail: format!(
                "no index of {} found (.csi beside it); write it with `ow convert --index`",
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

    fn cache_namespace(&self, py: Python<'_>) -> PyResult<CacheNamespace> {
        self.inner
            .run(py, |r| Ok(CacheNamespace::new(r.cache_namespace())))
    }

    fn describe_cache(
        &self,
        py: Python<'_>,
        chrom: &str,
        start: u32,
        end: u32,
    ) -> PyResult<Vec<CacheKey>> {
        self.inner.run(py, |r| {
            let mut keys = Vec::new();
            r.describe_cache(Region::new(chrom, start, end), &mut keys)?;
            Ok(keys.into_iter().map(CacheKey::new).collect())
        })
    }

    fn prefetch_cache(&self, py: Python<'_>, chrom: &str, start: u32, end: u32) -> PyResult<()> {
        self.inner.run(
            py,
            |r| Ok(r.prefetch_cache(Region::new(chrom, start, end))?),
        )
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

    fn stream_stats<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyDict>> {
        let stats = self.inner.run(py, |r| Ok(r.stream_stats()))?;
        seek_stats(py, stats)
    }

    fn reset_stream_stats(&self, py: Python<'_>) -> PyResult<()> {
        self.inner.run(py, |r| {
            r.reset_stats();
            Ok(())
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

    fn close(&self, py: Python<'_>) {
        drop(self.inner.take(py));
    }
}
