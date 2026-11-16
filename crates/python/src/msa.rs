//! Thin bindings for native bigMaf projection, with atomic query/write calls.

use std::path::PathBuf;

use numpy::{PyArray2, PyArrayMethods, PyUntypedArrayMethods};
use omicsweave_core::{FileSource, Region, Result, Strand};
use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;

use crate::cache::{self, PyBlockCache};
use crate::error::{Done, OrRaise};
use crate::shared::Shared;

#[pyclass(frozen, module = "omicsweave._native")]
pub struct MsaReader {
    inner: Shared<omicsweave_msa::MsaReader<FileSource>>,
    chroms: Vec<(String, u32)>,
}

#[pymethods]
impl MsaReader {
    #[new]
    #[pyo3(signature = (path, species=None, fill=b'-', cache=None))]
    fn new(
        py: Python<'_>,
        path: PathBuf,
        species: Option<Vec<String>>,
        fill: u8,
        cache: Option<PyRef<'_, PyBlockCache>>,
    ) -> PyResult<Self> {
        let cache = cache::handle(cache)?;
        py.detach(|| -> Result<Self> {
            let mut reader =
                omicsweave_msa::MsaReader::open(FileSource::open(path)?, species, fill)?;
            reader.set_cache(cache)?;
            let chroms = reader
                .chroms()
                .iter()
                .map(|(_, c)| (c.name.clone(), c.len))
                .collect();
            Ok(Self {
                inner: Shared::new(reader),
                chroms,
            })
        })
        .or_raise()
    }

    fn chroms(&self) -> Vec<(String, u32)> {
        self.chroms.clone()
    }

    fn set_cache(&self, py: Python<'_>, cache: Option<PyRef<'_, PyBlockCache>>) -> PyResult<()> {
        let cache = cache::handle(cache)?;
        self.inner.run(py, |r| Ok(r.set_cache(cache)?))
    }

    fn clone_context(&self, py: Python<'_>) -> PyResult<Self> {
        let reader = self.inner.run(py, |r| Ok(r.try_clone()?))?;
        Ok(Self {
            inner: Shared::new(reader),
            chroms: self.chroms.clone(),
        })
    }

    fn fetch<'py>(
        &self,
        py: Python<'py>,
        chrom: &str,
        start: u32,
        end: u32,
        reverse: bool,
    ) -> PyResult<(Vec<String>, Bound<'py, PyArray2<u8>>)> {
        let strand = if reverse {
            Strand::Reverse
        } else {
            Strand::Forward
        };
        let (names, array) = self.inner.run(py, |r| -> Done<_> {
            let shape = r.query(Region::new(chrom, start, end))?;
            let array = Python::attach(|py| -> Done<Py<PyArray2<u8>>> {
                // The NumPy C API allocates an owning array without invoking
                // a replaceable Python function while the reader is locked.
                let array = PyArray2::<u8>::zeros(py, [shape.species, shape.len], false);
                {
                    let mut borrowed = array.readwrite();
                    let out = borrowed
                        .as_slice_mut()
                        .map_err(|_| PyValueError::new_err("MSA array must be contiguous"))?;
                    py.detach(|| r.write(strand, out))?;
                }
                Ok(array.unbind())
            })?;
            Ok((r.species().to_vec(), array))
        })?;
        Ok((names, array.into_bound(py)))
    }

    fn fetch_into(
        &self,
        py: Python<'_>,
        chrom: &str,
        start: u32,
        end: u32,
        reverse: bool,
        out: &Bound<'_, PyArray2<u8>>,
    ) -> PyResult<()> {
        let shape = [out.shape()[0], out.shape()[1]];
        let mut out = out
            .try_readwrite()
            .map_err(|error| PyValueError::new_err(format!("out cannot be written: {error}")))?;
        let out = out
            .as_slice_mut()
            .map_err(|_| PyValueError::new_err("out must be a contiguous uint8 array"))?;
        self.inner.run(py, |r| -> Done<()> {
            let size = r.query(Region::new(chrom, start, end))?;
            if shape != [size.species, size.len] {
                return Err(PyValueError::new_err(format!(
                    "out must have shape ({}, {}), got ({}, {})",
                    size.species, size.len, shape[0], shape[1],
                ))
                .into());
            }
            Ok(r.write(
                if reverse {
                    Strand::Reverse
                } else {
                    Strand::Forward
                },
                out,
            )?)
        })
    }

    fn close(&self, py: Python<'_>) {
        drop(self.inner.take(py));
    }
}
