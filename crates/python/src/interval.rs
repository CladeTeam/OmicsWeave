//! Interval features (bigBed): reading by region as columns, and writing
//! bigBed.

use std::fs::File;
use std::io::BufWriter;
use std::path::PathBuf;

use numpy::{IntoPyArray, PyArray1, PyReadonlyArray1};
use omicsweave_bbi::{BbiReader, BedSchema, Kind};
use omicsweave_core::{
    FileSource, IntervalEntry, IntervalRead, IntervalSize, IntervalsOut, Region, Result,
    StringColumn,
};
use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;

use crate::cache::{self, PyBlockCache};
use crate::error::{Done, OrRaise};
use crate::shared::{Shared, Workers, path_str, save_sidecar, slice};
use crate::track::{bbi_chrom_list, bbi_chroms, open_bbi, write_options};

struct Intervals {
    reader: Box<BbiReader<FileSource>>,
    entries: Vec<IntervalEntry>,
}

type Columns = (Vec<u32>, Vec<u32>, Vec<i64>, Vec<u8>);

impl Intervals {
    fn query(&mut self, region: Region<'_>) -> Result<Columns> {
        let source: &mut dyn IntervalRead = self.reader.as_mut();
        let size: IntervalSize = source.query(region)?;
        self.entries.resize(size.entries, IntervalEntry::default());
        let mut offsets = vec![0i64; size.entries + 1];
        let mut data = vec![0u8; size.annotation_bytes];
        source.write(IntervalsOut {
            entries: &mut self.entries,
            annotations: StringColumn {
                offsets: &mut offsets,
                data: &mut data,
            },
        })?;
        let (starts, ends) = self.entries.iter().map(|e| (e.start, e.end)).unzip();
        Ok((starts, ends, offsets, data))
    }
}

#[pyclass(frozen, module = "omicsweave._native")]
pub struct IntervalReader {
    inner: Shared<Intervals>,
    chroms: Vec<(String, u32)>,
    auto_sql: Option<String>,
    field_count: (u16, u16),
}

#[pymethods]
impl IntervalReader {
    #[new]
    #[pyo3(signature = (path, cache=None))]
    fn new(
        py: Python<'_>,
        path: PathBuf,
        cache: Option<PyRef<'_, PyBlockCache>>,
    ) -> PyResult<IntervalReader> {
        let cache = cache::handle(cache)?;
        py.detach(|| -> Result<IntervalReader> {
            let mut reader = open_bbi(FileSource::open(&path)?, Kind::BigBed)?;
            reader.set_cache(cache)?;
            let header = reader.header();
            Ok(IntervalReader {
                chroms: bbi_chroms(&reader),
                auto_sql: reader.auto_sql()?,
                field_count: (header.field_count, header.defined_field_count),
                inner: Shared::new(Intervals {
                    reader: Box::new(reader),
                    entries: Vec::new(),
                }),
            })
        })
        .or_raise()
    }

    fn chroms(&self) -> Vec<(String, u32)> {
        self.chroms.clone()
    }

    fn set_cache(&self, py: Python<'_>, cache: Option<PyRef<'_, PyBlockCache>>) -> PyResult<()> {
        let cache = cache::handle(cache)?;
        self.inner.run(py, |t| Ok(t.reader.set_cache(cache)?))
    }

    fn clone_context(&self, py: Python<'_>) -> PyResult<Self> {
        let reader = self
            .inner
            .run(py, |t| Ok(Box::new(t.reader.try_clone()?)))?;
        Ok(Self {
            inner: Shared::new(Intervals {
                reader,
                entries: Vec::new(),
            }),
            chroms: self.chroms.clone(),
            auto_sql: self.auto_sql.clone(),
            field_count: self.field_count,
        })
    }

    #[getter]
    fn auto_sql(&self) -> Option<String> {
        self.auto_sql.clone()
    }

    /// (fields, standard fields).
    #[getter]
    fn field_count(&self) -> (u16, u16) {
        self.field_count
    }

    /// The features overlapping `chrom:[start, end)` as columns: starts,
    /// ends, and the annotations as Arrow string offsets and bytes.
    #[allow(clippy::type_complexity)]
    fn query<'py>(
        &self,
        py: Python<'py>,
        chrom: &str,
        start: u32,
        end: u32,
    ) -> PyResult<(
        Bound<'py, PyArray1<u32>>,
        Bound<'py, PyArray1<u32>>,
        Bound<'py, PyArray1<i64>>,
        Bound<'py, PyArray1<u8>>,
    )> {
        let (s, e, o, d) = self.inner.run(py, |t| -> Done<Columns> {
            Ok(t.query(Region::new(chrom, start, end))?)
        })?;
        Ok((
            s.into_pyarray(py),
            e.into_pyarray(py),
            o.into_pyarray(py),
            d.into_pyarray(py),
        ))
    }

    fn close(&self, py: Python<'_>) {
        drop(self.inner.take(py));
    }
}

#[pyclass(frozen, module = "omicsweave._native")]
pub struct BigBedWriter {
    inner: Shared<omicsweave_bbi::BigBedWriter<BufWriter<File>>>,
    path: PathBuf,
    /// The pool the writer's jobs run on, alive as long as it is.
    _workers: Workers,
}

#[pymethods]
impl BigBedWriter {
    #[new]
    #[allow(clippy::too_many_arguments)]
    #[pyo3(signature = (path, chroms, field_count, defined_field_count, auto_sql, items_per_slot, block_size, compress, level, zoom_levels, zoom_reductions, threads=1))]
    fn new(
        py: Python<'_>,
        path: PathBuf,
        chroms: Vec<(String, u32)>,
        field_count: u16,
        defined_field_count: u16,
        auto_sql: Option<String>,
        items_per_slot: Option<u32>,
        block_size: u32,
        compress: bool,
        level: i32,
        zoom_levels: usize,
        zoom_reductions: Option<Vec<u32>>,
        threads: usize,
    ) -> PyResult<BigBedWriter> {
        let workers = Workers::new(threads)?;
        let options = write_options(
            items_per_slot,
            block_size,
            compress,
            level,
            zoom_levels,
            zoom_reductions,
        );
        let schema = BedSchema {
            field_count,
            defined_field_count,
            auto_sql,
        };
        py.detach(|| -> Result<BigBedWriter> {
            let out = BufWriter::new(File::create(&path)?);
            let writer = omicsweave_bbi::BigBedWriter::create_with(
                out,
                bbi_chrom_list(chroms),
                schema,
                options,
                workers.exec(),
            )?;
            Ok(BigBedWriter {
                inner: Shared::new_writer(writer),
                path,
                _workers: workers,
            })
        })
        .or_raise()
    }

    /// Adds features of chromosome `chrom` (its index in the header); the
    /// annotations come as Arrow string offsets and bytes.
    fn add(
        &self,
        py: Python<'_>,
        chrom: u32,
        start: PyReadonlyArray1<'_, u32>,
        end: PyReadonlyArray1<'_, u32>,
        offsets: PyReadonlyArray1<'_, i64>,
        data: PyReadonlyArray1<'_, u8>,
    ) -> PyResult<()> {
        let (s, e) = (slice(&start, "start")?, slice(&end, "end")?);
        let (o, d) = (slice(&offsets, "offsets")?, slice(&data, "data")?);
        if s.len() != e.len() || o.len() != s.len() + 1 {
            return Err(PyValueError::new_err(
                "start, end and annotations differ in length",
            ));
        }
        for w in o.windows(2) {
            let (a, b) = (usize::try_from(w[0]), usize::try_from(w[1]));
            match (a, b) {
                (Ok(a), Ok(b)) if a <= b && b <= d.len() => {}
                _ => return Err(PyValueError::new_err("annotation offsets are invalid")),
            }
        }
        self.inner.run(py, |w| -> Done<()> {
            // Validate the entire batch before adding any feature. Each offset
            // fits usize and bounds a slice of data; no per-feature ranges
            // need to be retained while the writer consumes the columns.
            for ((&start, &end), offsets) in s.iter().zip(e).zip(o.windows(2)) {
                w.add(
                    chrom,
                    start,
                    end,
                    &d[offsets[0] as usize..offsets[1] as usize],
                )?;
            }
            Ok(())
        })
    }

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
