//! [`omicsweave_core::Error`] as the exceptions of `omicsweave.errors`,
//! one class per variant with its fields as attributes.

use omicsweave_core::Error;
use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;
use pyo3::types::PyDict;

/// Why work done with the GIL released failed.
pub enum Fail {
    Core(Error),
    /// The reader or writer was closed.
    Closed,
    Py(PyErr),
}

impl From<Error> for Fail {
    fn from(e: Error) -> Fail {
        Fail::Core(e)
    }
}

impl From<std::io::Error> for Fail {
    fn from(e: std::io::Error) -> Fail {
        Fail::Core(e.into())
    }
}

impl From<PyErr> for Fail {
    fn from(e: PyErr) -> Fail {
        Fail::Py(e)
    }
}

impl From<Fail> for PyErr {
    fn from(f: Fail) -> PyErr {
        match f {
            Fail::Core(e) => to_py(&e),
            Fail::Closed => PyValueError::new_err("I/O operation on closed file"),
            Fail::Py(e) => e,
        }
    }
}

pub type Done<T> = Result<T, Fail>;

/// Results of the format crates, raised as `omicsweave.errors`.
pub trait OrRaise<T> {
    fn or_raise(self) -> PyResult<T>;
}

impl<T> OrRaise<T> for omicsweave_core::Result<T> {
    fn or_raise(self) -> PyResult<T> {
        self.map_err(|e| to_py(&e))
    }
}

pub fn to_py(e: &Error) -> PyErr {
    Python::attach(|py| build(py, e).unwrap_or_else(|err| err))
}

fn build(py: Python<'_>, e: &Error) -> PyResult<PyErr> {
    // Errors are a cold path. Python's module registry already retains this
    // module; an extra native lazy lock could be inherited mid-initialization.
    let module = py.import("omicsweave.errors")?;
    let kw = PyDict::new(py);
    let class = match e {
        Error::Truncated { what, need, have } => {
            kw.set_item("what", what)?;
            kw.set_item("need", need)?;
            kw.set_item("have", have)?;
            "TruncatedError"
        }
        Error::BadMagic { what, found } => {
            kw.set_item("what", what)?;
            kw.set_item("found", found)?;
            "BadMagicError"
        }
        Error::BigEndian => "BigEndianError",
        Error::Unsupported(feature) => {
            kw.set_item("feature", feature)?;
            "UnsupportedError"
        }
        Error::OutOfRange {
            what,
            offset,
            lo,
            hi,
        } => {
            kw.set_item("what", what)?;
            kw.set_item("offset", offset)?;
            kw.set_item("lo", lo)?;
            kw.set_item("hi", hi)?;
            "OutOfRangeError"
        }
        Error::Invalid { what, detail } => {
            kw.set_item("what", what)?;
            kw.set_item("detail", detail)?;
            "InvalidError"
        }
        Error::TooDeep(depth) => {
            kw.set_item("depth", depth)?;
            "TooDeepError"
        }
        Error::Format { line, detail } => {
            kw.set_item("line", line)?;
            kw.set_item("detail", detail)?;
            "FormatError"
        }
        Error::Inflate(detail) => {
            kw.set_item("detail", detail)?;
            "InflateError"
        }
        Error::UnknownChrom(chrom) => {
            kw.set_item("chrom", chrom)?;
            "UnknownChromError"
        }
        Error::InvalidRegion { start, end } => {
            kw.set_item("start", start)?;
            kw.set_item("end", end)?;
            "InvalidRegionError"
        }
        Error::RegionOutOfRange { start, end, len } => {
            kw.set_item("start", start)?;
            kw.set_item("end", end)?;
            kw.set_item("length", len)?;
            "RegionOutOfRangeError"
        }
        Error::OutputSize { what, have, need } => {
            kw.set_item("what", what)?;
            kw.set_item("have", have)?;
            kw.set_item("need", need)?;
            "OutputSizeError"
        }
        Error::NotQueried => "NotQueriedError",
        Error::StaleSidecar => "StaleSidecarError",
        Error::Io { kind, detail } => {
            kw.set_item("kind", format!("{kind:?}"))?;
            kw.set_item("detail", detail)?;
            "IoError"
        }
    };
    let value = module.getattr(class)?.call((e.to_string(),), Some(&kw))?;
    Ok(PyErr::from_value(value))
}
