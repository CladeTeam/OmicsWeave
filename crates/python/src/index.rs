//! Building companion files explicitly, and telling formats apart.

use std::path::PathBuf;

use numpy::IntoPyArray;
use omicsweave_bbi::{BbiSidecar, BlockInfo, LevelBlocks};
use omicsweave_convert::Format;
use omicsweave_core::sidecar::sidecar_path;
use omicsweave_core::{Error, Result};
use omicsweave_fasta::{build_index_files, fai_path, gzi_path};
use omicsweave_hic::{BlockSize, HicSidecar, MatrixBlocks};
use pyo3::IntoPyObjectExt;
use pyo3::prelude::*;
use pyo3::types::PyDict;

use crate::error::OrRaise;
use crate::shared::path_str;

/// The name of the format the file at `path` holds, from its content.
#[pyfunction]
pub fn detect(py: Python<'_>, path: PathBuf) -> PyResult<&'static str> {
    py.detach(|| Format::detect(&path).map(Format::name))
        .or_raise()
}

/// Builds the `.fai` (and, for bgzip FASTA, the `.gzi`) beside `path`;
/// returns their paths.
#[pyfunction]
pub fn build_fasta_index(py: Python<'_>, path: PathBuf) -> PyResult<Vec<String>> {
    py.detach(|| -> Result<Vec<String>> {
        let (_, gzi) = build_index_files(&path)?;
        let mut saved = vec![path_str(&fai_path(&path))];
        if gzi.is_some() {
            saved.push(path_str(&gzi_path(&path)));
        }
        Ok(saved)
    })
    .or_raise()
}

/// Builds the sidecar of a bigWig, bigBed or `.hic` file beside it;
/// returns its path.
#[pyfunction]
pub fn build_sidecar(py: Python<'_>, path: PathBuf) -> PyResult<String> {
    py.detach(|| -> Result<String> {
        match Format::detect(&path)? {
            Format::BigWig | Format::BigBed => drop(BbiSidecar::build_file(&path)?),
            Format::Hic => drop(HicSidecar::build_file(&path)?),
            _ => {
                return Err(Error::invalid(
                    "sidecar",
                    "only bigWig, bigBed and .hic files have sidecars",
                ));
            }
        }
        Ok(path_str(&sidecar_path(&path)))
    })
    .or_raise()
}

type Columns<'py> = Bound<'py, PyDict>;

fn bbi_level<'py>(py: Python<'py>, level: &LevelBlocks) -> PyResult<(u32, Columns<'py>)> {
    let b = &level.blocks;
    let d = PyDict::new(py);
    let col = |f: fn(&BlockInfo) -> u64| b.iter().map(f).collect::<Vec<u64>>();
    d.set_item("start_chrom", col(|x| x.start.0 as u64).into_pyarray(py))?;
    d.set_item("start", col(|x| x.start.1 as u64).into_pyarray(py))?;
    d.set_item("end_chrom", col(|x| x.end.0 as u64).into_pyarray(py))?;
    d.set_item("end", col(|x| x.end.1 as u64).into_pyarray(py))?;
    d.set_item("offset", col(|x| x.offset).into_pyarray(py))?;
    d.set_item("size", col(|x| x.size as u64).into_pyarray(py))?;
    d.set_item("inflated", col(|x| x.inflated as u64).into_pyarray(py))?;
    d.set_item("items", col(|x| x.items as u64).into_pyarray(py))?;
    Ok((level.reduction, d))
}

fn hic_matrix<'py>(py: Python<'py>, m: &MatrixBlocks) -> PyResult<Bound<'py, PyAny>> {
    let b = &m.blocks;
    let d = PyDict::new(py);
    let col = |f: fn(&BlockSize) -> u64| b.iter().map(f).collect::<Vec<u64>>();
    d.set_item("number", col(|x| x.number as u64).into_pyarray(py))?;
    d.set_item("position", col(|x| x.position).into_pyarray(py))?;
    d.set_item("size", col(|x| x.size as u64).into_pyarray(py))?;
    d.set_item("inflated", col(|x| x.inflated as u64).into_pyarray(py))?;
    d.set_item("records", col(|x| x.records as u64).into_pyarray(py))?;
    (
        m.chr1,
        m.chr2,
        m.unit.as_str(),
        m.bin_size,
        m.block_bin_count,
        m.block_column_count,
        d,
    )
        .into_bound_py_any(py)
}

/// The sidecar beside `path` (a bigWig, bigBed or `.hic` file), checked
/// against the file: `("bbi", uncompress buffer size, [(reduction,
/// columns)])` or `("hic", version, [(chr1, chr2, unit, bin size, block
/// bins, block columns, columns)])`.
#[pyfunction]
pub fn read_sidecar(py: Python<'_>, path: PathBuf) -> PyResult<Bound<'_, PyAny>> {
    enum Loaded {
        Bbi(BbiSidecar),
        Hic(HicSidecar),
    }
    let loaded = py
        .detach(|| -> Result<Loaded> {
            Ok(match Format::detect(&path)? {
                Format::BigWig | Format::BigBed => Loaded::Bbi(BbiSidecar::load(&path)?),
                Format::Hic => Loaded::Hic(HicSidecar::load(&path)?),
                _ => {
                    return Err(Error::invalid(
                        "sidecar",
                        "only bigWig, bigBed and .hic files have sidecars",
                    ));
                }
            })
        })
        .or_raise()?;
    match loaded {
        Loaded::Bbi(s) => {
            let levels = s
                .levels
                .iter()
                .map(|l| bbi_level(py, l))
                .collect::<PyResult<Vec<_>>>()?;
            ("bbi", s.uncompress_buf_size, levels).into_bound_py_any(py)
        }
        Loaded::Hic(s) => {
            let maps = s
                .matrices
                .iter()
                .map(|m| hic_matrix(py, m))
                .collect::<PyResult<Vec<_>>>()?;
            ("hic", s.version, maps).into_bound_py_any(py)
        }
    }
}
