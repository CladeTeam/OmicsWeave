//! Building companion files explicitly, and telling formats apart.

use std::path::PathBuf;

use omicsweave_bbi::BbiSidecar;
use omicsweave_convert::Format;
use omicsweave_core::sidecar::sidecar_path;
use omicsweave_core::{Error, Result};
use omicsweave_fasta::{build_index_files, fai_path, gzi_path};
use omicsweave_hic::HicSidecar;
use pyo3::prelude::*;

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
