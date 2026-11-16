//! FFM motif collections: a resident collection read by number or ID, a
//! writer, and the conversions between the package's `Alphabet` / `Motif`
//! dataclasses and the format crate's values. Dataclasses are read
//! attribute by attribute; a missing or mistyped attribute is an error.

use std::fs::File;
use std::io::BufWriter;
use std::path::PathBuf;
use std::sync::{Arc, Mutex};

use numpy::{IntoPyArray, PyArray1, PyReadonlyArray1};
use omicsweave_core::{FileSource, Result};
use omicsweave_ffm::{
    Alphabet, FactorKind, FfmWriter as Writer, Gaps, INVALID, Key, KeyType, Motif, MotifSet,
    MotifSetBuilder, MotifView, Schema, Support, WriterOptions,
};
use pyo3::exceptions::{PyIndexError, PyTypeError, PyValueError};
use pyo3::prelude::*;
use pyo3::types::{PyBytes, PyDict, PyList};

use crate::error::OrRaise;
use crate::shared::Shared;

pub fn attr<'py, T: FromPyObjectOwned<'py>>(o: &Bound<'py, PyAny>, name: &str) -> PyResult<T> {
    o.getattr(name)?.extract::<T>().map_err(Into::into)
}

pub fn ffm<T>(r: omicsweave_ffm::Result<T>) -> PyResult<T> {
    r.or_raise()
}

/// A one-character string as the byte it names (code points below 256).
pub fn byte_of(s: &str, what: &str) -> PyResult<u8> {
    let mut chars = s.chars();
    match (chars.next(), chars.next()) {
        (Some(c), None) if (c as u32) < 256 => Ok(c as u32 as u8),
        _ => Err(PyValueError::new_err(format!(
            "{what} must be one character below U+0100, not {s:?}"
        ))),
    }
}

/// A string as bytes, one per character (code points below 256).
pub fn bytes_of(s: &str, what: &str) -> PyResult<Vec<u8>> {
    s.chars()
        .map(|c| {
            u8::try_from(c as u32).map_err(|_| {
                PyValueError::new_err(format!("{what} has a character above U+00FF: {c:?}"))
            })
        })
        .collect()
}

/// Bytes as a string, one character per byte.
pub fn text_of(bytes: &[u8]) -> String {
    bytes.iter().map(|&b| b as char).collect()
}

// ---------------------------------------------------------------------------
// Alphabets

/// An `omicsweave.types.Alphabet`.
pub fn alphabet(o: &Bound<'_, PyAny>) -> PyResult<Alphabet> {
    let name: String = attr(o, "name")?;
    let symbols = bytes_of(&attr::<String>(o, "symbols")?, "alphabet symbols")?;
    let mut a = ffm(Alphabet::new(name, &symbols))?;
    for item in o.getattr("aliases")?.call_method0("items")?.try_iter()? {
        let (byte, symbol): (String, String) = item?.extract()?;
        a = ffm(a.with_alias(
            byte_of(&byte, "an alias")?,
            byte_of(&symbol, "an alias target")?,
        ))?;
    }
    if let Some(c) = attr::<Option<String>>(o, "complement")? {
        a = ffm(a.with_complement(&bytes_of(&c, "alphabet complement")?))?;
    }
    Ok(a)
}

pub type AlphabetFields<'py> = (String, String, Option<String>, Bound<'py, PyDict>);

/// `(name, symbols, complement, aliases)`: the fields of the Python alphabet.
pub fn alphabet_fields<'py>(py: Python<'py>, a: &Alphabet) -> PyResult<AlphabetFields<'py>> {
    let symbols = a.symbols();
    let complement = a
        .complement()
        .map(|c| c.iter().map(|&k| symbols[k as usize] as char).collect());
    let aliases = PyDict::new(py);
    for (byte, &code) in a.lut().iter().enumerate() {
        if code != INVALID && !symbols.contains(&(byte as u8)) {
            aliases.set_item(
                (byte as u8 as char).to_string(),
                (symbols[code as usize] as char).to_string(),
            )?;
        }
    }
    Ok((a.name().to_owned(), text_of(symbols), complement, aliases))
}

#[pyfunction]
pub fn alphabet_preset<'py>(py: Python<'py>, name: &str) -> PyResult<AlphabetFields<'py>> {
    let a = match name {
        "dna" => Alphabet::dna(),
        "rna" => Alphabet::rna(),
        "protein" => Alphabet::protein(),
        _ => {
            return Err(PyValueError::new_err(format!(
                "no alphabet preset {name:?}"
            )));
        }
    };
    alphabet_fields(py, &a)
}

// ---------------------------------------------------------------------------
// Motifs

fn factor_kind(s: &str) -> PyResult<FactorKind> {
    Ok(match s {
        "joint" => FactorKind::Joint,
        "conditional" => FactorKind::Conditional,
        "potential" => FactorKind::Potential,
        _ => {
            return Err(PyValueError::new_err(format!(
                "factor kind must be joint, conditional or potential, not {s:?}"
            )));
        }
    })
}

pub fn kind_name(k: FactorKind) -> &'static str {
    match k {
        FactorKind::Joint => "joint",
        FactorKind::Conditional => "conditional",
        FactorKind::Potential => "potential",
    }
}

pub fn floats(o: &Bound<'_, PyAny>) -> PyResult<Vec<f64>> {
    let a: PyReadonlyArray1<'_, f64> = o.extract()?;
    Ok(a.as_array().iter().copied().collect())
}

/// An `omicsweave.types.Motif`'s model.
pub fn motif(o: &Bound<'_, PyAny>) -> PyResult<Motif> {
    let support = match attr::<Option<Vec<u8>>>(o, "support")? {
        Some(codes) => Support::Codes(codes),
        None => Support::All,
    };
    let mut factors = Vec::new();
    for f in o.getattr("factors")?.try_iter()? {
        let f = f?;
        factors.push(omicsweave_ffm::Factor {
            kind: factor_kind(&attr::<String>(&f, "kind")?)?,
            positions: attr(&f, "positions")?,
            table: floats(&f.getattr("table")?)?,
        });
    }
    let gaps = o.getattr("gaps")?;
    let gaps = if gaps.is_none() {
        None
    } else {
        Some(Gaps {
            deletions: floats(&gaps.getattr("deletions")?)?,
            insertions: floats(&gaps.getattr("insertions")?)?,
        })
    };
    Ok(Motif {
        width: attr(o, "width")?,
        support,
        sites: attr(o, "sites")?,
        factors,
        gaps,
    })
}

type FactorFields<'py> = (
    &'static str,
    Bound<'py, PyArray1<u32>>,
    Bound<'py, PyArray1<f64>>,
);
type GapFields<'py> = (Bound<'py, PyArray1<f64>>, Bound<'py, PyArray1<f64>>);
type MotifFields<'py> = (
    u32,
    Option<Bound<'py, PyBytes>>,
    Option<f64>,
    Vec<FactorFields<'py>>,
    Option<GapFields<'py>>,
    Vec<String>,
);

/// Motifs ready to compile: a loaded collection, or the package's `Motif`
/// values packed into one (their IDs are not needed and need not be unique).
pub fn motif_set(source: &Bound<'_, PyAny>, alphabet: &Alphabet) -> PyResult<Arc<MotifSet>> {
    if let Ok(c) = source.cast::<MotifCollection>() {
        let set = c.get().set()?;
        if set.alphabet() != alphabet {
            return Err(PyValueError::new_err(
                "alphabet differs from the collection's",
            ));
        }
        return Ok(set);
    }
    let list = source
        .cast::<PyList>()
        .map_err(|_| PyTypeError::new_err("motifs must be a MotifCollection or a list of Motif"))?;
    let mut builder = MotifSetBuilder::new(alphabet.clone(), Schema::default());
    for (i, m) in list.iter().enumerate() {
        let model = motif(&m)?;
        builder
            .push(&i.to_string(), &model, &[])
            .map_err(|e| PyValueError::new_err(format!("motif {i}: {e}")))?;
    }
    Ok(Arc::new(builder.finish()))
}

/// Every motif of `set`, in order.
pub fn views(set: &MotifSet) -> Vec<MotifView<'_>> {
    (0..set.len()).filter_map(|i| set.motif(i)).collect()
}

// ---------------------------------------------------------------------------
// The resident collection

/// An FFM collection read whole and checked on open.
#[pyclass(frozen, module = "omicsweave._native")]
pub struct MotifCollection {
    set: Mutex<Option<Arc<MotifSet>>>,
}

impl MotifCollection {
    pub fn set(&self) -> PyResult<Arc<MotifSet>> {
        self.set
            .lock()
            .expect("motif collection lock")
            .clone()
            .ok_or_else(|| PyValueError::new_err("I/O operation on closed file"))
    }
}

#[pymethods]
impl MotifCollection {
    #[new]
    fn new(py: Python<'_>, path: PathBuf) -> PyResult<MotifCollection> {
        let set = py
            .detach(|| -> Result<MotifSet> { MotifSet::load(FileSource::open(&path)?) })
            .or_raise()?;
        Ok(MotifCollection {
            set: Mutex::new(Some(Arc::new(set))),
        })
    }

    fn alphabet<'py>(&self, py: Python<'py>) -> PyResult<AlphabetFields<'py>> {
        alphabet_fields(py, self.set()?.alphabet())
    }

    fn schema(&self) -> PyResult<Vec<(String, &'static str)>> {
        Ok(self
            .set()?
            .schema()
            .keys()
            .iter()
            .map(|k| {
                let kind = match k.kind {
                    KeyType::String => "string",
                    KeyType::Integer => "integer",
                    KeyType::Float => "float",
                };
                (k.name.clone(), kind)
            })
            .collect())
    }

    fn ids(&self) -> PyResult<Vec<String>> {
        let set = self.set()?;
        Ok((0..set.len())
            .map(|i| set.id(i).unwrap_or_default().to_owned())
            .collect())
    }

    fn __len__(&self) -> PyResult<usize> {
        Ok(self.set()?.len())
    }

    fn widths<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyArray1<u32>>> {
        let set = self.set()?;
        let w: Vec<u32> = views(&set).iter().map(MotifView::width).collect();
        Ok(w.into_pyarray(py))
    }

    /// `(width, support, sites, factors, gaps, annotation values)`.
    fn motif<'py>(&self, py: Python<'py>, i: usize) -> PyResult<MotifFields<'py>> {
        let set = self.set()?;
        let view = set
            .motif(i)
            .ok_or_else(|| PyIndexError::new_err("motif index out of range"))?;
        let support = view.support().map(|codes| {
            let bytes: Vec<u8> = codes.iter().map(|&c| c as u8).collect();
            PyBytes::new(py, &bytes)
        });
        let factors = view
            .factors()
            .map(|f| {
                (
                    kind_name(f.kind),
                    f.positions.to_vec().into_pyarray(py),
                    f.table.to_vec().into_pyarray(py),
                )
            })
            .collect();
        let gaps = view.gaps().map(|g| {
            (
                g.deletions().to_vec().into_pyarray(py),
                g.insertions().to_vec().into_pyarray(py),
            )
        });
        let values = set
            .annotation(i)
            .map(|a| a.values().map(str::to_owned).collect())
            .unwrap_or_default();
        Ok((view.width(), support, view.sites(), factors, gaps, values))
    }

    fn close(&self) {
        self.set.lock().expect("motif collection lock").take();
    }
}

// ---------------------------------------------------------------------------
// The writer

type Sink = Writer<BufWriter<File>>;

/// Writes one FFM collection; the file is written on `close`.
#[pyclass(module = "omicsweave._native")]
pub struct FfmWriter {
    inner: Shared<Sink>,
    keys: usize,
}

#[pymethods]
impl FfmWriter {
    #[new]
    fn new(
        path: PathBuf,
        alphabet: &Bound<'_, PyAny>,
        schema: Vec<(String, String)>,
        block_records: usize,
        block_values: usize,
    ) -> PyResult<FfmWriter> {
        let alphabet = self::alphabet(alphabet)?;
        let keys = schema
            .iter()
            .map(|(name, kind)| {
                let kind = match kind.as_str() {
                    "string" => KeyType::String,
                    "integer" => KeyType::Integer,
                    "float" => KeyType::Float,
                    _ => {
                        return Err(PyValueError::new_err(format!(
                            "annotation type must be string, integer or float, not {kind:?}"
                        )));
                    }
                };
                Ok(Key::new(name.clone(), kind))
            })
            .collect::<PyResult<Vec<_>>>()?;
        let n = keys.len();
        let schema = ffm(Schema::new(keys))?;
        let options = WriterOptions {
            block_records,
            block_values,
            ..WriterOptions::default()
        };
        let file = File::create(&path)
            .map_err(omicsweave_core::Error::from)
            .or_raise()?;
        let writer = ffm(Writer::new(BufWriter::new(file), alphabet, schema, options))?;
        Ok(FfmWriter {
            inner: Shared::new_writer(writer),
            keys: n,
        })
    }

    fn add(&self, py: Python<'_>, motif: &Bound<'_, PyAny>, values: Vec<String>) -> PyResult<()> {
        self.inner.ensure_open()?;
        if values.len() != self.keys {
            return Err(PyValueError::new_err("one annotation value per schema key"));
        }
        let id: String = attr(motif, "id")?;
        let model = self::motif(motif)?;
        let values: Vec<&str> = values.iter().map(String::as_str).collect();
        self.inner
            .run(py, |writer| Ok(writer.push(&id, &model, &values)?))
    }

    fn __len__(&self, py: Python<'_>) -> PyResult<usize> {
        self.inner.ensure_process()?;
        Ok(py.detach(|| self.inner.lock().as_ref().map_or(0, Sink::len)))
    }

    fn close(&self, py: Python<'_>) -> PyResult<()> {
        let Some(writer) = self.inner.take(py) else {
            return Ok(());
        };
        py.detach(|| -> Result<()> {
            let out = writer.finish()?;
            out.into_inner()
                .map_err(|e| omicsweave_core::Error::from(e.into_error()))?
                .sync_all()?;
            Ok(())
        })
        .or_raise()
    }
}
