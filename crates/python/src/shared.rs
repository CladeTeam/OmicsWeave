//! State shared between Python threads, and small helpers.

use std::path::{Path, PathBuf};
#[cfg(unix)]
use std::sync::OnceLock;
use std::sync::atomic::{AtomicBool, AtomicU32, Ordering};
use std::sync::{Mutex, MutexGuard};

use numpy::{Element, PyReadonlyArray1};
use omicsweave_bgzf::SeekStats;
use omicsweave_core::sidecar::{self, SourceIdentity};
use omicsweave_core::{Error, Locus, ReadAt, Result, Span, Strand};
use pyo3::exceptions::{PyRuntimeError, PyValueError};
use pyo3::prelude::*;
use pyo3::pybacked::PyBackedStr;
use pyo3::types::PyDict;

use crate::error::{Done, Fail};

static PROCESS_ID: AtomicU32 = AtomicU32::new(0);

#[cfg(unix)]
unsafe extern "C" {
    fn getpid() -> std::ffi::c_int;
    fn pthread_atfork(
        prepare: Option<unsafe extern "C" fn()>,
        parent: Option<unsafe extern "C" fn()>,
        child: Option<unsafe extern "C" fn()>,
    ) -> std::ffi::c_int;
}

#[cfg(unix)]
unsafe extern "C" fn native_after_fork_child() {
    // SAFETY: POSIX getpid is async-signal-safe and cannot fail. The
    // lock-free AtomicU32 store allocates nothing and cannot unwind.
    // No Python API or inherited synchronization primitive is touched.
    let process = unsafe { getpid() };
    PROCESS_ID.store(process as u32, Ordering::Relaxed);
}

/// Registers the process marker before Python's child callbacks can run.
pub fn initialize_process_guard() -> PyResult<()> {
    #[cfg(unix)]
    {
        static REGISTERED: OnceLock<std::ffi::c_int> = OnceLock::new();
        let error = *REGISTERED.get_or_init(|| {
            // SAFETY: the child callback is a static C-ABI function which
            // only performs async-signal-safe work, takes no arguments and
            // does not unwind. The libc registration retains its address.
            unsafe { pthread_atfork(None, None, Some(native_after_fork_child)) }
        });
        if error != 0 {
            return Err(PyRuntimeError::new_err(format!(
                "could not register the fork process guard (error {error})"
            )));
        }
    }
    PROCESS_ID.store(std::process::id(), Ordering::Relaxed);
    Ok(())
}

/// The native at-fork hook updates this before Python child callbacks run.
/// The hot path needs only an atomic load, without a getpid syscall per query.
pub fn process_id() -> u32 {
    let process = PROCESS_ID.load(Ordering::Relaxed);
    if process == 0 {
        let process = std::process::id();
        PROCESS_ID.store(process, Ordering::Relaxed);
        process
    } else {
        process
    }
}

#[pyfunction]
pub fn _after_fork_child() {
    PROCESS_ID.store(std::process::id(), Ordering::Relaxed);
}

/// A reader or writer, `None` once closed. The lock is taken only with the
/// GIL released.
pub struct Shared<T> {
    value: Mutex<Option<T>>,
    process: u32,
    abandon_inherited_writer: bool,
    closed: AtomicBool,
}

impl<T: Send> Shared<T> {
    pub fn new(value: T) -> Shared<T> {
        Shared {
            value: Mutex::new(Some(value)),
            process: std::process::id(),
            abandon_inherited_writer: false,
            closed: AtomicBool::new(false),
        }
    }

    /// A sink's destructor can flush buffered bytes. A forked copy must never
    /// write them to its parent's file or wait for vanished writer threads.
    pub fn new_writer(value: T) -> Shared<T> {
        Shared {
            value: Mutex::new(Some(value)),
            process: std::process::id(),
            abandon_inherited_writer: true,
            closed: AtomicBool::new(false),
        }
    }

    pub fn ensure_process(&self) -> Done<()> {
        if self.process != process_id() {
            return Err(Fail::Py(PyRuntimeError::new_err(
                "reader/writer cannot be reused after fork; reopen it in this process",
            )));
        }
        Ok(())
    }

    /// Rejects invalid execution state before a caller allocates an output.
    pub fn ensure_open(&self) -> Done<()> {
        self.ensure_process()?;
        if self.closed.load(Ordering::Acquire) {
            return Err(Fail::Closed);
        }
        Ok(())
    }

    pub fn lock(&self) -> MutexGuard<'_, Option<T>> {
        self.value
            .lock()
            .unwrap_or_else(|poisoned| poisoned.into_inner())
    }

    /// Runs `f` on the value with the GIL released.
    pub fn run<R: Send>(
        &self,
        py: Python<'_>,
        f: impl FnOnce(&mut T) -> Done<R> + Send,
    ) -> PyResult<R> {
        self.ensure_process()?;
        Ok(py.detach(|| match self.lock().as_mut() {
            Some(value) => f(value),
            None => Err(Fail::Closed),
        })?)
    }

    /// Takes the value out, closing it.
    pub fn take(&self, py: Python<'_>) -> Option<T> {
        if self.process != process_id() {
            return None;
        }
        py.detach(|| {
            let value = self.lock().take();
            self.closed.store(true, Ordering::Release);
            value
        })
    }
}

impl<T> Drop for Shared<T> {
    fn drop(&mut self) {
        if self.abandon_inherited_writer && self.process != std::process::id() {
            // Destruction provides exclusive access, so get_mut cannot wait
            // on an inherited lock. Forgetting just this foreign sink also
            // suppresses BufWriter flushing and worker-pool joins; the child
            // OS process releases its retained file descriptors at exit.
            let value = self
                .value
                .get_mut()
                .unwrap_or_else(|p| p.into_inner())
                .take();
            if let Some(writer) = value {
                std::mem::forget(writer);
            }
        }
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

/// Physical BGZF work, separate from cache lookup counters and output size.
pub fn seek_stats(py: Python<'_>, stats: SeekStats) -> PyResult<Bound<'_, PyDict>> {
    let out = PyDict::new(py);
    out.set_item("read_calls", stats.read_calls)?;
    out.set_item("read_bytes", stats.read_bytes)?;
    out.set_item("inflate_calls", stats.inflate_calls)?;
    out.set_item("compacted_bytes", stats.compacted_bytes)?;
    out.set_item("discarded_bytes", stats.discarded_bytes)?;
    Ok(out)
}

/// A dense-query key as the Python layer passes it: chromosome, ascending
/// `(start, end)` spans and the strand as `+`, `-`, `.` or `?`.
pub type LocusKey = (PyBackedStr, Vec<(u32, u32)>, PyBackedStr);

/// Runs `f` on the [`Locus`] of `key`.
pub fn with_locus<R>(key: &LocusKey, f: impl FnOnce(&Locus<'_>) -> Result<R>) -> Result<R> {
    let (chrom, spans, strand) = key;
    let strand = match strand.as_bytes() {
        [b] => Strand::from_byte(*b),
        _ => None,
    }
    .ok_or_else(|| Error::invalid("strand", format!("{strand:?} is not one of + - . ?")))?;
    // Common loci fit on the stack; the borrowed Python strings stay owned
    // by key while the decoding callback runs with the GIL released.
    if spans.len() <= 4 {
        let mut local = [Span::new(0, 0); 4];
        for (out, &(s, e)) in local.iter_mut().zip(spans) {
            *out = Span::new(s, e);
        }
        f(&Locus::new(chrom, &local[..spans.len()], strand)?)
    } else {
        let spans: Vec<Span> = spans.iter().map(|&(s, e)| Span::new(s, e)).collect();
        f(&Locus::new(chrom, &spans, strand)?)
    }
}
