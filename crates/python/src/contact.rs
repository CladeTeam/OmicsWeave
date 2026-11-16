//! Contact maps (`.hic`): reading tiles, and writing maps or contacts.

use std::fs::File;
use std::io::BufWriter;
use std::path::PathBuf;

use numpy::{IntoPyArray, PyArray1, PyReadonlyArray1, PyReadwriteArray2, PyUntypedArrayMethods};
use omicsweave_core::{
    ContactEntry, ContactQuery, ContactRead, ContactUnit, FileSource, Region, Result,
};
use omicsweave_hic::{
    BlockBins, BlockLayout, BuildOptions, Chrom, HicHeader, HicReader, HicSidecar, Record,
};
use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;

use crate::cache::{self, PyBlockCache};
use crate::error::{Done, OrRaise};
use crate::shared::{Shared, Workers, path_str, save_sidecar, slice};

fn unit(name: &str) -> PyResult<ContactUnit> {
    ContactUnit::parse(name)
        .ok_or_else(|| PyValueError::new_err(format!("unit must be 'BP' or 'FRAG', not {name:?}")))
}

/// A tile key as the Python query class passes it.
type Tile = (String, u32, u32, String, u32, u32, u32, String);

fn query(tile: &Tile) -> PyResult<ContactQuery<'_>> {
    let (rc, rs, re, cc, cs, ce, resolution, u) = tile;
    Ok(ContactQuery::new(
        Region::new(rc, *rs, *re),
        Region::new(cc, *cs, *ce),
        *resolution,
    )
    .with_unit(unit(u)?))
}

struct Contacts {
    reader: HicReader<FileSource>,
    entries: Vec<ContactEntry>,
}

#[pyclass(frozen, module = "omicsweave._native")]
pub struct ContactReader {
    inner: Shared<Contacts>,
    version: i32,
    genome: String,
    attributes: Vec<(String, String)>,
    chroms: Vec<(String, u64)>,
    bp_resolutions: Vec<u32>,
    frag_resolutions: Vec<u32>,
}

#[pymethods]
impl ContactReader {
    #[new]
    #[pyo3(signature = (path, cache=None))]
    fn new(
        py: Python<'_>,
        path: PathBuf,
        cache: Option<PyRef<'_, PyBlockCache>>,
    ) -> PyResult<ContactReader> {
        let cache = cache::handle(cache)?;
        py.detach(|| -> Result<ContactReader> {
            let mut reader = HicReader::open(FileSource::open(&path)?)?;
            reader.set_cache(cache);
            let h = reader.header();
            Ok(ContactReader {
                version: h.version,
                genome: h.genome.clone(),
                attributes: h.attributes.clone(),
                chroms: h.chroms.iter().map(|c| (c.name.clone(), c.len)).collect(),
                bp_resolutions: h.bp_resolutions.clone(),
                frag_resolutions: h.frag_resolutions.clone(),
                inner: Shared::new(Contacts {
                    reader,
                    entries: Vec::new(),
                }),
            })
        })
        .or_raise()
    }

    #[getter]
    fn version(&self) -> i32 {
        self.version
    }

    #[getter]
    fn genome(&self) -> String {
        self.genome.clone()
    }

    fn attributes(&self) -> Vec<(String, String)> {
        self.attributes.clone()
    }

    fn set_cache(&self, py: Python<'_>, cache: Option<PyRef<'_, PyBlockCache>>) -> PyResult<()> {
        let cache = cache::handle(cache)?;
        self.inner.run(py, |c| {
            c.reader.set_cache(cache);
            Ok(())
        })
    }

    /// The chromosome pairs, by name, that have a contact map, in
    /// chromosome order.
    fn maps(&self, py: Python<'_>) -> PyResult<Vec<(String, String)>> {
        let chroms = &self.chroms;
        self.inner.run(py, |c| {
            let mut pairs: Vec<(u32, u32)> = c.reader.footer().matrices.keys().copied().collect();
            pairs.sort_unstable();
            Ok(pairs
                .into_iter()
                .filter_map(|(a, b)| {
                    Some((
                        chroms.get(a as usize)?.0.clone(),
                        chroms.get(b as usize)?.0.clone(),
                    ))
                })
                .collect())
        })
    }

    fn clone_context(&self, py: Python<'_>) -> PyResult<Self> {
        let reader = self.inner.run(py, |c| Ok(c.reader.try_clone_context()?))?;
        Ok(Self {
            inner: Shared::new(Contacts {
                reader,
                entries: Vec::new(),
            }),
            version: self.version,
            genome: self.genome.clone(),
            attributes: self.attributes.clone(),
            chroms: self.chroms.clone(),
            bp_resolutions: self.bp_resolutions.clone(),
            frag_resolutions: self.frag_resolutions.clone(),
        })
    }

    fn chroms(&self) -> Vec<(String, u64)> {
        self.chroms.clone()
    }

    fn resolutions(&self, unit_name: &str) -> PyResult<Vec<u32>> {
        Ok(match unit(unit_name)? {
            ContactUnit::Bp => self.bp_resolutions.clone(),
            ContactUnit::Frag => self.frag_resolutions.clone(),
        })
    }

    fn fragment_count(&self, py: Python<'_>, chrom: &str) -> PyResult<u32> {
        self.inner
            .run(py, |c| -> Done<u32> { Ok(c.reader.fragment_count(chrom)?) })
    }

    /// The restriction sites of `chrom`.
    fn sites<'py>(&self, py: Python<'py>, chrom: &str) -> PyResult<Bound<'py, PyArray1<u32>>> {
        let sites = self.inner.run(py, |c| -> Done<Vec<u32>> {
            let index = c.reader.chrom_index(chrom)?;
            Ok(c.reader.sites(index)?.to_vec())
        })?;
        Ok(sites.into_pyarray(py))
    }

    /// Writes the tile row-major into `out`; rows (columns) run from the
    /// high coordinate down when `flip_rows` (`flip_cols`), for keys on the
    /// reverse strand.
    #[pyo3(signature = (tile, out, flip_rows=false, flip_cols=false))]
    fn matrix(
        &self,
        py: Python<'_>,
        tile: Tile,
        mut out: PyReadwriteArray2<'_, f32>,
        flip_rows: bool,
        flip_cols: bool,
    ) -> PyResult<()> {
        let q = query(&tile)?;
        q.check().or_raise()?;
        let (rows, cols) = q.shape();
        if out.shape() != [rows, cols] {
            return Err(PyValueError::new_err(format!(
                "out shape must be ({rows}, {cols}), not {:?}",
                out.shape()
            )));
        }
        // A contiguous NumPy array may also have column-major strides;
        // the reader writes cells in row-major order.
        if !out.is_c_contiguous() {
            return Err(PyValueError::new_err(
                "out must be a C-contiguous float32 array",
            ));
        }
        let out = out
            .as_slice_mut()
            .map_err(|_| PyValueError::new_err("out must be a C-contiguous float32 array"))?;
        self.inner.run(py, |c| -> Done<()> {
            c.reader.matrix(q, out)?;
            omicsweave_core::orient::flip_tile(out, rows, cols, flip_rows, flip_cols);
            Ok(())
        })
    }

    /// The tile's nonzero cells as (row, col, value) columns.
    #[allow(clippy::type_complexity)]
    fn cells<'py>(
        &self,
        py: Python<'py>,
        tile: Tile,
    ) -> PyResult<(
        Bound<'py, PyArray1<u32>>,
        Bound<'py, PyArray1<u32>>,
        Bound<'py, PyArray1<f32>>,
    )> {
        let q = query(&tile)?;
        let (r, k, v) = self.inner.run(py, |c| -> Done<_> {
            let n = c.reader.query_entries(q)?;
            c.entries.resize(n, ContactEntry::default());
            c.reader.write_entries(&mut c.entries)?;
            let r: Vec<_> = c.entries.iter().map(|x| x.row).collect();
            let k: Vec<_> = c.entries.iter().map(|x| x.col).collect();
            let v: Vec<_> = c.entries.iter().map(|x| x.value).collect();
            Ok((r, k, v))
        })?;
        Ok((r.into_pyarray(py), k.into_pyarray(py), v.into_pyarray(py)))
    }

    fn close(&self, py: Python<'_>) {
        drop(self.inner.take(py));
    }
}

/// A header as the Python header class passes it.
#[allow(clippy::too_many_arguments)]
fn header(
    version: i32,
    genome: String,
    attributes: Vec<(String, String)>,
    chroms: Vec<(String, u64)>,
    bp_resolutions: Vec<u32>,
    frag_resolutions: Vec<u32>,
    sites: Vec<Vec<u32>>,
) -> HicHeader {
    let mut h = HicHeader::new(
        &genome,
        chroms
            .into_iter()
            .map(|(name, len)| Chrom { name, len })
            .collect(),
        bp_resolutions,
    );
    h.version = version;
    h.attributes = attributes;
    h.frag_resolutions = frag_resolutions;
    h.sites = sites;
    h
}

fn columns<'a>(
    a: &'a PyReadonlyArray1<'_, u32>,
    b: &'a PyReadonlyArray1<'_, u32>,
    v: &'a PyReadonlyArray1<'_, f32>,
) -> PyResult<(&'a [u32], &'a [u32], &'a [f32])> {
    let (a, b, v) = (slice(a, "x")?, slice(b, "y")?, slice(v, "value")?);
    if a.len() != b.len() || a.len() != v.len() {
        return Err(PyValueError::new_err("columns differ in length"));
    }
    Ok((a, b, v))
}

fn close_with_sidecar(
    path: &std::path::Path,
    made: HicSidecar,
    sidecar: bool,
) -> Result<Option<String>> {
    if !sidecar {
        return Ok(None);
    }
    let at = save_sidecar(path, |id| made.to_bytes(id))?;
    Ok(Some(path_str(&at)))
}

#[pyclass(frozen, module = "omicsweave._native")]
pub struct HicWriter {
    inner: Shared<omicsweave_hic::HicWriter<BufWriter<File>>>,
    path: PathBuf,
    /// The pool the writer's jobs run on, alive as long as it is.
    _workers: Workers,
}

/// Bins per block side, within and between chromosomes (`None`: the
/// Juicer grid); `None` as a whole is the default layout.
type Blocks = Option<(Option<u32>, Option<u32>)>;

fn layout(blocks: Blocks) -> BlockLayout {
    let bins = |b: Option<u32>| b.map_or(BlockBins::Juicer, BlockBins::Fixed);
    match blocks {
        None => BlockLayout::default(),
        Some((intra, inter)) => BlockLayout {
            intra: bins(intra),
            inter: bins(inter),
        },
    }
}

#[pymethods]
impl HicWriter {
    #[new]
    #[allow(clippy::too_many_arguments)]
    #[pyo3(signature = (path, version, genome, attributes, chroms, bp_resolutions, frag_resolutions, sites, level, blocks=None, threads=1))]
    fn new(
        py: Python<'_>,
        path: PathBuf,
        version: i32,
        genome: String,
        attributes: Vec<(String, String)>,
        chroms: Vec<(String, u64)>,
        bp_resolutions: Vec<u32>,
        frag_resolutions: Vec<u32>,
        sites: Vec<Vec<u32>>,
        level: u32,
        blocks: Blocks,
        threads: usize,
    ) -> PyResult<HicWriter> {
        let workers = Workers::new(threads)?;
        let h = header(
            version,
            genome,
            attributes,
            chroms,
            bp_resolutions,
            frag_resolutions,
            sites,
        );
        py.detach(|| -> Result<HicWriter> {
            let out = BufWriter::new(File::create(&path)?);
            let writer = omicsweave_hic::HicWriter::create_exec(
                out,
                h,
                level,
                layout(blocks),
                workers.exec(),
            )?;
            Ok(HicWriter {
                inner: Shared::new_writer(writer),
                path,
                _workers: workers,
            })
        })
        .or_raise()
    }

    /// Writes the map of chromosomes `chr1` x `chr2` (header indexes) at one
    /// resolution, from bin columns.
    #[allow(clippy::too_many_arguments)]
    fn write_matrix(
        &self,
        py: Python<'_>,
        chr1: u32,
        chr2: u32,
        unit_name: &str,
        bin_size: u32,
        x: PyReadonlyArray1<'_, u32>,
        y: PyReadonlyArray1<'_, u32>,
        value: PyReadonlyArray1<'_, f32>,
    ) -> PyResult<()> {
        let u = unit(unit_name)?;
        let (x, y, v) = columns(&x, &y, &value)?;
        self.inner.run(py, |w| -> Done<()> {
            let records: Vec<Record> = (0..x.len())
                .map(|i| Record {
                    x: x[i],
                    y: y[i],
                    value: v[i],
                })
                .collect();
            Ok(w.write_matrix(chr1, chr2, u, bin_size, &records)?)
        })
    }

    fn close(&self, py: Python<'_>, sidecar: bool) -> PyResult<Option<String>> {
        let Some(writer) = self.inner.take(py) else {
            return Ok(None);
        };
        py.detach(|| -> Result<Option<String>> {
            let (out, made) = writer.finish_with_sidecar()?;
            out.into_inner().map_err(|e| e.into_error())?.sync_all()?;
            close_with_sidecar(&self.path, made, sidecar)
        })
        .or_raise()
    }
}

#[pyclass(frozen, module = "omicsweave._native")]
pub struct HicBuilder {
    inner: Shared<omicsweave_hic::HicBuilder>,
    path: PathBuf,
    options: BuildOptions,
    workers: Workers,
}

#[pymethods]
impl HicBuilder {
    #[new]
    #[allow(clippy::too_many_arguments)]
    #[pyo3(signature = (path, version, genome, attributes, chroms, bp_resolutions, frag_resolutions, sites, level, blocks=None, threads=1))]
    fn new(
        path: PathBuf,
        version: i32,
        genome: String,
        attributes: Vec<(String, String)>,
        chroms: Vec<(String, u64)>,
        bp_resolutions: Vec<u32>,
        frag_resolutions: Vec<u32>,
        sites: Vec<Vec<u32>>,
        level: u32,
        blocks: Blocks,
        threads: usize,
    ) -> PyResult<HicBuilder> {
        let h = header(
            version,
            genome,
            attributes,
            chroms,
            bp_resolutions,
            frag_resolutions,
            sites,
        );
        Ok(HicBuilder {
            inner: Shared::new(omicsweave_hic::HicBuilder::new(h).or_raise()?),
            path,
            options: BuildOptions {
                level,
                blocks: layout(blocks),
            },
            workers: Workers::new(threads)?,
        })
    }

    /// Adds contacts: `count[i]` between base `pos1[i]` of `chr1[i]` and
    /// base `pos2[i]` of `chr2[i]` (header indexes).
    fn add(
        &self,
        py: Python<'_>,
        chr1: PyReadonlyArray1<'_, u32>,
        pos1: PyReadonlyArray1<'_, u32>,
        chr2: PyReadonlyArray1<'_, u32>,
        pos2: PyReadonlyArray1<'_, u32>,
        count: PyReadonlyArray1<'_, f32>,
    ) -> PyResult<()> {
        let (c1, p1, n) = columns(&chr1, &pos1, &count)?;
        let (c2, p2) = (slice(&chr2, "chr2")?, slice(&pos2, "pos2")?);
        if c2.len() != c1.len() || p2.len() != c1.len() {
            return Err(PyValueError::new_err("columns differ in length"));
        }
        self.inner.run(py, |b| -> Done<()> {
            for i in 0..c1.len() {
                b.add(c1[i], p1[i], c2[i], p2[i], n[i])?;
            }
            Ok(())
        })
    }

    /// Bins every contact at every resolution and writes the file.
    fn close(&self, py: Python<'_>, sidecar: bool) -> PyResult<Option<String>> {
        let Some(builder) = self.inner.take(py) else {
            return Ok(None);
        };
        py.detach(|| -> Result<Option<String>> {
            let out = BufWriter::new(File::create(&self.path)?);
            let (out, made) = builder.finish_exec(out, &self.options, self.workers.exec())?;
            out.into_inner().map_err(|e| e.into_error())?.sync_all()?;
            close_with_sidecar(&self.path, made, sidecar)
        })
        .or_raise()
    }
}
