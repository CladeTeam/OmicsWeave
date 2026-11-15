//! Signal tracks (bigWig): reading values and intervals, and writing.

use std::fs::File;
use std::io::BufWriter;
use std::path::PathBuf;

use numpy::{
    IntoPyArray, PyArray1, PyReadonlyArray1, PyReadwriteArray1, PyReadwriteArray2,
    PyUntypedArrayMethods,
};
use omicsweave_bbi::{BbiReader, Chrom, Kind, WriteOptions};
use omicsweave_core::{
    BinSource, BinSpec, Error, FileSource, Region, Result, Stat, TrackEntry, TrackRead,
};
use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;

use crate::cache::{self, CacheKey, CacheNamespace, PyBlockCache};
use crate::error::{Done, OrRaise};
use crate::shared::{LocusKey, Shared, path_str, save_sidecar, slice, with_locus};

/// Opens a BBI file, which must be of `kind`.
pub fn open_bbi(source: FileSource, kind: Kind) -> Result<BbiReader<FileSource>> {
    let reader = BbiReader::open(source)?;
    if reader.header().kind != kind {
        return Err(Error::invalid(
            "BBI file",
            match kind {
                Kind::BigWig => "a bigBed file, not a bigWig",
                Kind::BigBed => "a bigWig file, not a bigBed",
            },
        ));
    }
    Ok(reader)
}

pub fn bbi_chroms(reader: &BbiReader<FileSource>) -> Vec<(String, u32)> {
    reader
        .chroms()
        .iter()
        .map(|(_, c)| (c.name.clone(), c.len))
        .collect()
}

struct Track {
    reader: BbiReader<FileSource>,
    entries: Vec<TrackEntry>,
}

#[pyclass(frozen, module = "omicsweave._native")]
pub struct TrackReader {
    inner: Shared<Track>,
    chroms: Vec<(String, u32)>,
    /// (bases covered, min, max, sum, sum of squares).
    summary: Option<(u64, f64, f64, f64, f64)>,
}

#[pymethods]
impl TrackReader {
    #[new]
    #[pyo3(signature = (path, cache=None))]
    fn new(
        py: Python<'_>,
        path: PathBuf,
        cache: Option<PyRef<'_, PyBlockCache>>,
    ) -> PyResult<TrackReader> {
        let cache = cache::handle(cache)?;
        py.detach(|| -> Result<TrackReader> {
            let mut reader = open_bbi(FileSource::open(&path)?, Kind::BigWig)?;
            reader.set_cache(cache)?;
            let summary = reader
                .total_summary()?
                .map(|s| (s.bases_covered, s.min, s.max, s.sum, s.sum_squares));
            Ok(TrackReader {
                chroms: bbi_chroms(&reader),
                summary,
                inner: Shared::new(Track {
                    reader,
                    entries: Vec::new(),
                }),
            })
        })
        .or_raise()
    }

    fn chroms(&self) -> Vec<(String, u32)> {
        self.chroms.clone()
    }

    fn summary(&self) -> Option<(u64, f64, f64, f64, f64)> {
        self.summary
    }

    fn set_cache(&self, py: Python<'_>, cache: Option<PyRef<'_, PyBlockCache>>) -> PyResult<()> {
        let cache = cache::handle(cache)?;
        self.inner.run(py, |t| Ok(t.reader.set_cache(cache)?))
    }

    fn cache_namespace(&self, py: Python<'_>) -> PyResult<CacheNamespace> {
        self.inner
            .run(py, |t| Ok(CacheNamespace::new(t.reader.cache_namespace())))
    }

    fn describe_cache(&self, py: Python<'_>, key: LocusKey) -> PyResult<Vec<CacheKey>> {
        self.inner.run(py, |t| {
            let mut keys = Vec::new();
            with_locus(&key, |locus| t.reader.describe_cache(locus, &mut keys))?;
            Ok(keys.into_iter().map(CacheKey::new).collect())
        })
    }

    fn prefetch_cache(&self, py: Python<'_>, key: LocusKey) -> PyResult<()> {
        self.inner.run(py, |t| {
            Ok(with_locus(&key, |locus| t.reader.prefetch_cache(locus))?)
        })
    }

    fn clone_context(&self, py: Python<'_>) -> PyResult<Self> {
        let reader = self.inner.run(py, |t| Ok(t.reader.try_clone()?))?;
        Ok(Self {
            inner: Shared::new(Track {
                reader,
                entries: Vec::new(),
            }),
            chroms: self.chroms.clone(),
            summary: self.summary,
        })
    }

    /// Writes one value per position of `key` into `out`, 5' to 3'.
    fn values(
        &self,
        py: Python<'_>,
        key: LocusKey,
        mut out: PyReadwriteArray1<'_, f32>,
    ) -> PyResult<()> {
        let out = out
            .as_slice_mut()
            .map_err(|_| PyValueError::new_err("out must be a contiguous float32 array"))?;
        self.inner.run(py, |t| -> Done<()> {
            Ok(with_locus(&key, |locus| t.reader.values(locus, out))?)
        })
    }

    /// Writes the `stats` of the bins of `key` into `out`, shaped
    /// `(len(stats), bins)`, C-contiguous.
    #[pyo3(signature = (key, resolution, phase, stats, exact, out))]
    #[allow(clippy::too_many_arguments)]
    fn bins(
        &self,
        py: Python<'_>,
        key: LocusKey,
        resolution: u32,
        phase: u32,
        stats: Vec<String>,
        exact: bool,
        mut out: PyReadwriteArray2<'_, f32>,
    ) -> PyResult<()> {
        let stats = stats
            .iter()
            .map(|s| {
                Stat::parse(s).ok_or_else(|| {
                    PyValueError::new_err(format!(
                        "unknown statistic {s:?}; one of mean, min, max, sum, std, coverage"
                    ))
                })
            })
            .collect::<PyResult<Vec<_>>>()?;
        if out.shape()[0] != stats.len() {
            return Err(PyValueError::new_err(format!(
                "out has {} rows, need one per statistic ({})",
                out.shape()[0],
                stats.len()
            )));
        }
        if !out.is_c_contiguous() {
            return Err(PyValueError::new_err(
                "out must be a C-contiguous float32 array",
            ));
        }
        let out = out
            .as_slice_mut()
            .map_err(|_| PyValueError::new_err("out must be a C-contiguous float32 array"))?;
        let spec = BinSpec {
            resolution,
            phase,
            stats: &stats,
            source: if exact {
                BinSource::Exact
            } else {
                BinSource::Zoom
            },
        };
        self.inner.run(py, |t| -> Done<()> {
            Ok(with_locus(&key, |locus| t.reader.bins(locus, &spec, out))?)
        })
    }

    /// Zoom reductions in file order, coarsening.
    fn zoom_reductions(&self, py: Python<'_>) -> PyResult<Vec<u32>> {
        self.inner.run(py, |t| {
            Ok(t.reader.zooms().iter().map(|z| z.reduction).collect())
        })
    }

    /// The intervals with data in `chrom:[start, end)`, clipped to it, as
    /// (start, end, value) columns.
    #[allow(clippy::type_complexity)]
    fn intervals<'py>(
        &self,
        py: Python<'py>,
        chrom: &str,
        start: u32,
        end: u32,
    ) -> PyResult<(
        Bound<'py, PyArray1<u32>>,
        Bound<'py, PyArray1<u32>>,
        Bound<'py, PyArray1<f32>>,
    )> {
        let (s, e, v) = self.inner.run(py, |t| -> Done<_> {
            let n = t.reader.query_entries(Region::new(chrom, start, end))?;
            t.entries.resize(n, TrackEntry::default());
            t.reader.write_entries(&mut t.entries)?;
            let s: Vec<_> = t.entries.iter().map(|x| x.start).collect();
            let e: Vec<_> = t.entries.iter().map(|x| x.end).collect();
            let v: Vec<_> = t.entries.iter().map(|x| x.value).collect();
            Ok((s, e, v))
        })?;
        Ok((s.into_pyarray(py), e.into_pyarray(py), v.into_pyarray(py)))
    }

    fn close(&self, py: Python<'_>) {
        drop(self.inner.take(py));
    }
}

/// The options of a BBI writer, as the Python options class holds them.
pub fn write_options(
    items_per_slot: Option<u32>,
    block_size: u32,
    compress: bool,
    level: i32,
    zoom_levels: usize,
    zoom_reductions: Option<Vec<u32>>,
) -> WriteOptions {
    WriteOptions {
        items_per_slot,
        block_size,
        compress,
        level,
        zoom_levels,
        zoom_reductions,
    }
}

pub fn bbi_chrom_list(chroms: Vec<(String, u32)>) -> Vec<Chrom> {
    chroms
        .into_iter()
        .map(|(name, len)| Chrom { name, len })
        .collect()
}

#[pyclass(frozen, module = "omicsweave._native")]
pub struct BigWigWriter {
    inner: Shared<omicsweave_bbi::BigWigWriter<BufWriter<File>>>,
    path: PathBuf,
}

#[pymethods]
impl BigWigWriter {
    #[new]
    #[allow(clippy::too_many_arguments)]
    #[pyo3(signature = (path, chroms, items_per_slot, block_size, compress, level, zoom_levels, zoom_reductions))]
    fn new(
        py: Python<'_>,
        path: PathBuf,
        chroms: Vec<(String, u32)>,
        items_per_slot: Option<u32>,
        block_size: u32,
        compress: bool,
        level: i32,
        zoom_levels: usize,
        zoom_reductions: Option<Vec<u32>>,
    ) -> PyResult<BigWigWriter> {
        let options = write_options(
            items_per_slot,
            block_size,
            compress,
            level,
            zoom_levels,
            zoom_reductions,
        );
        py.detach(|| -> Result<BigWigWriter> {
            let out = BufWriter::new(File::create(&path)?);
            let writer =
                omicsweave_bbi::BigWigWriter::create(out, bbi_chrom_list(chroms), options)?;
            Ok(BigWigWriter {
                inner: Shared::new_writer(writer),
                path,
            })
        })
        .or_raise()
    }

    /// Adds items of chromosome `chrom` (its index in the header).
    fn add(
        &self,
        py: Python<'_>,
        chrom: u32,
        start: PyReadonlyArray1<'_, u32>,
        end: PyReadonlyArray1<'_, u32>,
        value: PyReadonlyArray1<'_, f32>,
    ) -> PyResult<()> {
        let (s, e, v) = (
            slice(&start, "start")?,
            slice(&end, "end")?,
            slice(&value, "value")?,
        );
        if s.len() != e.len() || s.len() != v.len() {
            return Err(PyValueError::new_err(
                "start, end and value differ in length",
            ));
        }
        self.inner.run(py, |w| -> Done<()> {
            for i in 0..s.len() {
                w.add(chrom, s[i], e[i], v[i])?;
            }
            Ok(())
        })
    }

    /// Ends the file; with `sidecar`, saves its sidecar beside it. Returns
    /// the sidecar's path.
    fn close(&self, py: Python<'_>, sidecar: bool) -> PyResult<Option<String>> {
        let Some(writer) = self.inner.take(py) else {
            return Ok(None);
        };
        py.detach(|| -> Result<Option<String>> {
            let (out, made) = writer.finish()?;
            out.into_inner().map_err(|e| e.into_error())?.sync_all()?;
            if !sidecar {
                return Ok(None);
            }
            let at = save_sidecar(&self.path, |id| made.to_bytes(id))?;
            Ok(Some(path_str(&at)))
        })
        .or_raise()
    }
}
