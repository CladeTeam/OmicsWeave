//! State shared between Python threads, and small helpers.

use std::path::{Path, PathBuf};
use std::sync::{Mutex, MutexGuard};

use numpy::{Element, PyReadonlyArray1};
use omicsweave_core::sidecar::{self, SourceIdentity};
use omicsweave_core::{ReadAt, Result};
use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;

use crate::error::{Done, Fail};

/// A reader or writer, `None` once closed. The lock is taken only with the
/// GIL released.
pub struct Shared<T>(Mutex<Option<T>>);

impl<T: Send> Shared<T> {
    pub fn new(value: T) -> Shared<T> {
        Shared(Mutex::new(Some(value)))
    }

    pub fn lock(&self) -> MutexGuard<'_, Option<T>> {
        self.0
            .lock()
            .unwrap_or_else(|poisoned| poisoned.into_inner())
    }

    /// Runs `f` on the value with the GIL released.
    pub fn run<R: Send>(
        &self,
        py: Python<'_>,
        f: impl FnOnce(&mut T) -> Done<R> + Send,
    ) -> PyResult<R> {
        Ok(py.detach(|| match self.lock().as_mut() {
            Some(value) => f(value),
            None => Err(Fail::Closed),
        })?)
    }

    /// Takes the value out, closing it.
    pub fn take(&self, py: Python<'_>) -> Option<T> {
        py.detach(|| self.lock().take())
    }
}

/// The contiguous data of an input array.
pub fn slice<'a, T: Element>(array: &'a PyReadonlyArray1<'_, T>, what: &str) -> PyResult<&'a [T]> {
    array
        .as_slice()
        .map_err(|_| PyValueError::new_err(format!("{what} must be a contiguous array")))
}

/// Saves the sidecar of the file just written at `path`.
pub fn save_sidecar(path: &Path, bytes: impl FnOnce(SourceIdentity) -> Vec<u8>) -> Result<PathBuf> {
    let identity = SourceIdentity::of(&omicsweave_core::FileSource::open(path)?)?;
    let at = sidecar::sidecar_path(path);
    sidecar::write_file(&at, &bytes(identity))?;
    Ok(at)
}

/// The first `n` bytes of a source (fewer if it is shorter).
pub fn head(source: &dyn ReadAt, n: usize) -> Result<Vec<u8>> {
    let mut buf = vec![0; n.min(source.len() as usize)];
    source.read_exact_at(&mut buf, 0)?;
    Ok(buf)
}

pub fn path_str(path: &Path) -> String {
    path.to_string_lossy().into_owned()
}
