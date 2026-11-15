//! Explicit process-local control of caller-owned decoded blocks.

use std::sync::Arc;

use crate::shared::process_id;
use omicsweave_core::cache::{BlockCache, BlockKey, SourceToken};
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

#[pyclass(frozen, eq, hash, skip_from_py_object, module = "omicsweave.cache")]
#[derive(Clone, PartialEq, Eq, Hash)]
pub struct CacheNamespace {
    pub token: SourceToken,
    process: u32,
}

impl CacheNamespace {
    pub fn new(token: SourceToken) -> Self {
        Self {
            token,
            process: process_id(),
        }
    }
}

#[pymethods]
impl CacheNamespace {
    fn __repr__(&self) -> String {
        format!("CacheNamespace({})", self.token.get())
    }
}

#[pyclass(
    frozen,
    eq,
    hash,
    skip_from_py_object,
    name = "BlockKey",
    module = "omicsweave.cache"
)]
#[derive(Clone, PartialEq, Eq, Hash)]
pub struct CacheKey {
    key: BlockKey,
    process: u32,
}

impl CacheKey {
    pub fn new(key: BlockKey) -> Self {
        Self {
            key,
            process: process_id(),
        }
    }
}

#[pymethods]
impl CacheKey {
    #[getter]
    fn source(&self) -> PyResult<CacheNamespace> {
        check_process(self.process)?;
        Ok(CacheNamespace::new(self.key.source))
    }

    fn __repr__(&self) -> String {
        format!(
            "CacheKey(source={}, kind={}, stream={:?}, block={}, variant={})",
            self.key.source.get(),
            self.key.kind,
            self.key.stream,
            self.key.block,
            self.key.variant
        )
    }
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

    fn keys(&self, py: Python<'_>) -> PyResult<Vec<CacheKey>> {
        let cache = self.handle()?;
        Ok(py
            .detach(|| cache.keys())
            .into_iter()
            .map(CacheKey::new)
            .collect())
    }

    fn evict(&self, py: Python<'_>, key: &CacheKey) -> PyResult<bool> {
        check_process(key.process)?;
        let cache = self.handle()?;
        Ok(py.detach(|| cache.evict(key.key)))
    }

    #[pyo3(signature = (source=None))]
    fn clear(&self, py: Python<'_>, source: Option<&CacheNamespace>) -> PyResult<()> {
        let cache = self.handle()?;
        let source = source
            .map(|s| {
                check_process(s.process)?;
                Ok::<_, PyErr>(s.token)
            })
            .transpose()?;
        py.detach(|| match source {
            Some(s) => cache.clear_source(s),
            None => cache.clear(),
        });
        Ok(())
    }

    fn resize(&self, py: Python<'_>, max_bytes: usize, max_entries: usize) -> PyResult<()> {
        let cache = self.handle()?;
        py.detach(|| cache.resize(max_bytes, max_entries));
        Ok(())
    }
}
