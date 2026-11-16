//! A process-local store of decoded blocks that readers share. Which
//! blocks a query needs and how they are scheduled belong to the native
//! layer (loaders); Python attaches a cache, sizes it and reads its
//! counters.

use std::sync::Arc;

use crate::shared::process_id;
use omicsweave_core::cache::BlockCache;
use pyo3::exceptions::PyRuntimeError;
use pyo3::prelude::*;
use pyo3::types::PyDict;

fn check_process(process: u32) -> PyResult<()> {
    if process != process_id() {
        return Err(PyRuntimeError::new_err(
            "cache objects cannot be reused after fork; create a cache in this process",
        ));
    }
    Ok(())
}

#[pyclass(frozen, name = "BlockCache", module = "omicsweave._native")]
pub struct PyBlockCache {
    pub inner: Arc<BlockCache>,
    process: u32,
}

impl PyBlockCache {
    pub fn handle(&self) -> PyResult<Arc<BlockCache>> {
        check_process(self.process)?;
        Ok(self.inner.clone())
    }
}

pub fn handle(cache: Option<PyRef<'_, PyBlockCache>>) -> PyResult<Option<Arc<BlockCache>>> {
    cache.map(|c| c.handle()).transpose()
}

#[pymethods]
impl PyBlockCache {
    #[new]
    fn new(max_bytes: usize, max_entries: usize) -> Self {
        Self {
            inner: Arc::new(BlockCache::new(max_bytes, max_entries)),
            process: std::process::id(),
        }
    }

    fn stats<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyDict>> {
        let cache = self.handle()?;
        let s = py.detach(|| cache.stats());
        let d = PyDict::new(py);
        d.set_item("hits", s.hits)?;
        d.set_item("misses", s.misses)?;
        d.set_item("loads", s.loads)?;
        d.set_item("waits", s.waits)?;
        d.set_item("admission_rejections", s.admission_rejections)?;
        d.set_item("evictions", s.evictions)?;
        d.set_item("resident_bytes", s.resident_bytes)?;
        d.set_item("resident_entries", s.resident_entries)?;
        d.set_item("live_bytes", s.live_bytes)?;
        d.set_item("live_entries", s.live_entries)?;
        d.set_item("in_flight", s.in_flight)?;
        d.set_item("max_bytes", s.max_bytes)?;
        d.set_item("max_entries", s.max_entries)?;
        Ok(d)
    }

    /// Drops every entry; blocks leased by a running query stay valid.
    fn clear(&self, py: Python<'_>) -> PyResult<()> {
        let cache = self.handle()?;
        py.detach(|| cache.clear());
        Ok(())
    }

    fn resize(&self, py: Python<'_>, max_bytes: usize, max_entries: usize) -> PyResult<()> {
        let cache = self.handle()?;
        py.detach(|| cache.resize(max_bytes, max_entries));
        Ok(())
    }
}
