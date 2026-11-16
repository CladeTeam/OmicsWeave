//! Sequences: `.fa` (+ `.fai`), `.fa.gz` (+ `.fai`, `.gzi`) and `.2bit`;
//! and writing FASTA with its indexes.

use std::collections::HashMap;
use std::fs::File;
use std::io::BufWriter;
use std::path::{Path, PathBuf};
use std::sync::Arc;

use numpy::{
    Element, PY_ARRAY_API, PyArray1, PyArrayDescrMethods, PyArrayMethods, PyReadwriteArray1,
};
use omicsweave_core::cache::{BlockCache, BlockKey, SourceToken};
use omicsweave_core::{Error, FileSource, Locus, Result, SequenceRead};
use omicsweave_fasta::{Fai, FastaReader, Gzi, TWOBIT_MAGIC, TwoBitReader, fai_path, gzi_path};
use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;

use crate::cache::{self, CacheKey, CacheNamespace, PyBlockCache};
use crate::error::{Done, OrRaise};
use crate::shared::{LocusKey, Shared, head, path_str, with_locus};

pub enum Sequences {
    Fasta(FastaReader<FileSource>),
    TwoBit(TwoBitReader<FileSource>),
}

impl Sequences {
    fn set_cache(&mut self, cache: Option<Arc<BlockCache>>) -> Result<()> {
        match self {
            Self::Fasta(r) => r.set_cache(cache),
            Self::TwoBit(r) => r.set_cache(cache),
        }
    }

    fn namespace(&self) -> SourceToken {
        match self {
            Self::Fasta(r) => r.cache_namespace(),
            Self::TwoBit(r) => r.cache_namespace(),
        }
    }

    fn describe(&self, locus: &Locus<'_>, keys: &mut Vec<BlockKey>) -> Result<()> {
        match self {
            Self::Fasta(r) => r.describe_cache(locus, keys),
            Self::TwoBit(r) => r.describe_cache(locus, keys),
        }
    }

    fn prefetch(&mut self, locus: &Locus<'_>) -> Result<()> {
        match self {
            Self::Fasta(r) => r.prefetch_cache(locus),
            Self::TwoBit(r) => r.prefetch_cache(locus),
        }
    }

    pub(crate) fn try_clone(&self) -> Result<Self> {
        match self {
            Self::Fasta(r) => Ok(Self::Fasta(r.try_clone()?)),
            Self::TwoBit(r) => Ok(Self::TwoBit(r.try_clone()?)),
        }
    }
}

impl SequenceRead for Sequences {
    fn chrom_len(&self, chrom: &str) -> Result<u32> {
        match self {
            Sequences::Fasta(r) => r.chrom_len(chrom),
            Sequences::TwoBit(r) => r.chrom_len(chrom),
        }
    }

    fn sequence(&mut self, locus: &Locus<'_>, out: &mut [u8]) -> Result<()> {
        match self {
            Sequences::Fasta(r) => r.sequence(locus, out),
            Sequences::TwoBit(r) => r.sequence(locus, out),
        }
    }
}

/// What a sequence file holds, from its first bytes.
#[derive(Clone, Copy, PartialEq, Eq)]
enum SequenceKind {
    Fasta,
    Bgzf,
    TwoBit,
}

impl SequenceKind {
    fn of(head: &[u8]) -> SequenceKind {
        if head.len() >= 4 {
            let magic = u32::from_le_bytes(head[..4].try_into().unwrap());
            if magic == TWOBIT_MAGIC || magic.swap_bytes() == TWOBIT_MAGIC {
                return SequenceKind::TwoBit;
            }
        }
        if head.starts_with(&[0x1f, 0x8b]) {
            SequenceKind::Bgzf
        } else {
            SequenceKind::Fasta
        }
    }
}

/// The index at `given`, or else beside the FASTA file at `beside`; it
/// must exist (indexes are built explicitly, never on open).
fn index<T>(given: Option<&Path>, beside: PathBuf, parse: fn(&[u8]) -> Result<T>) -> Result<T> {
    let path = given.map_or(beside, Path::to_path_buf);
    match std::fs::read(&path) {
        Ok(bytes) => parse(&bytes),
        Err(e) if e.kind() == std::io::ErrorKind::NotFound => Err(Error::Io {
            kind: e.kind(),
            detail: format!(
                "index {} not found; build it with omicsweave.build_index",
                path.display()
            ),
        }),
        Err(e) => Err(e.into()),
    }
}

/// Opens the sequence file at `path` by its content (FASTA, bgzip FASTA,
/// 2bit) with the indexes beside it.
pub(crate) fn open_sequences(path: &Path) -> Result<Sequences> {
    let source = FileSource::open(path)?;
    Ok(match SequenceKind::of(&head(&source, 4)?) {
        SequenceKind::TwoBit => Sequences::TwoBit(TwoBitReader::open(source, false)?),
        kind => Sequences::Fasta(open_fasta(
            path,
            source,
            kind == SequenceKind::Bgzf,
            None,
            None,
        )?),
    })
}

/// Opens `.fa` with its `.fai`, or `.fa.gz` with its `.fai` and `.gzi`.
fn open_fasta(
    path: &Path,
    source: FileSource,
    bgzf: bool,
    fai: Option<&Path>,
    gzi: Option<&Path>,
) -> Result<FastaReader<FileSource>> {
    let found_fai = index(fai, fai_path(path), Fai::parse)?;
    let found_gzi = if bgzf {
        Some(index(gzi, gzi_path(path), Gzi::parse)?)
    } else {
        None
    };
    Ok(FastaReader::new(source, found_fai, found_gzi))
}

fn zeros(py: Python<'_>, length: usize) -> PyResult<Bound<'_, PyArray1<u8>>> {
    let mut dimension = isize::try_from(length)
        .map_err(|_| PyValueError::new_err("sequence output exceeds NumPy dimension limit"))?;
    // SAFETY: the one-element dimension array remains live for this call,
    // its nonnegative value fits npy_intp, and the dtype is exactly uint8.
    // PyArray_Zeros steals the owned dtype reference and returns an owned,
    // fully initialized ndarray or NULL with a Python allocation exception.
    // The fallible owned-pointer conversion preserves that exception.
    let array = unsafe {
        let pointer = PY_ARRAY_API.PyArray_Zeros(
            py,
            1,
            &mut dimension,
            u8::get_dtype(py).into_dtype_ptr(),
            0,
        );
        Bound::from_owned_ptr_or_err(py, pointer)?
    };
    Ok(array.cast_into::<PyArray1<u8>>()?)
}

#[pyclass(frozen, module = "omicsweave._native")]
pub struct SequenceReader {
    pub inner: Shared<Sequences>,
    chroms: Arc<[(Arc<str>, u32)]>,
    lengths: Arc<HashMap<Arc<str>, u32>>,
}

impl SequenceReader {
    fn checked_length(&self, key: &LocusKey) -> PyResult<usize> {
        self.inner.ensure_open()?;
        with_locus(key, |locus| {
            let length = self
                .lengths
                .get(locus.chrom)
                .ok_or_else(|| Error::UnknownChrom(locus.chrom.to_owned()))?;
            locus.check_within(*length)?;
            Ok(locus.len() as usize)
        })
        .or_raise()
    }
}

#[pymethods]
impl SequenceReader {
    #[new]
    #[pyo3(signature = (path, fai=None, gzi=None, soft_mask=false, cache=None))]
    fn new(
        py: Python<'_>,
        path: PathBuf,
        fai: Option<PathBuf>,
        gzi: Option<PathBuf>,
        soft_mask: bool,
        cache: Option<PyRef<'_, PyBlockCache>>,
    ) -> PyResult<SequenceReader> {
        let cache = cache::handle(cache)?;
        py.detach(|| -> Result<SequenceReader> {
            let source = FileSource::open(&path)?;
            let kind = SequenceKind::of(&head(&source, 4)?);
            let mut reader = match kind {
                SequenceKind::TwoBit => Sequences::TwoBit(TwoBitReader::open(source, soft_mask)?),
                _ => Sequences::Fasta(open_fasta(
                    &path,
                    source,
                    kind == SequenceKind::Bgzf,
                    fai.as_deref(),
                    gzi.as_deref(),
                )?),
            };
            if cache.is_some() {
                reader.set_cache(cache)?;
            }
            let chroms: Vec<(Arc<str>, u32)> = match &reader {
                Sequences::Fasta(r) => r
                    .fai()
                    .entries()
                    .iter()
                    .map(|e| (Arc::from(e.name.as_str()), e.len))
                    .collect(),
                Sequences::TwoBit(r) => r
                    .names()
                    .map(|n| Ok((Arc::from(n), r.chrom_len(n)?)))
                    .collect::<Result<_>>()?,
            };
            let lengths = Arc::new(
                chroms
                    .iter()
                    .map(|(name, len)| (Arc::clone(name), *len))
                    .collect(),
            );
            Ok(SequenceReader {
                inner: Shared::new(reader),
                chroms: chroms.into(),
                lengths,
            })
        })
        .or_raise()
    }

    /// (name, length) of every sequence, in file order.
    fn chroms(&self) -> Vec<(String, u32)> {
        self.chroms
            .iter()
            .map(|(name, len)| (name.to_string(), *len))
            .collect()
    }

    fn set_cache(&self, py: Python<'_>, cache: Option<PyRef<'_, PyBlockCache>>) -> PyResult<()> {
        let cache = cache::handle(cache)?;
        self.inner.run(py, |r| Ok(r.set_cache(cache)?))
    }

    fn cache_namespace(&self, py: Python<'_>) -> PyResult<CacheNamespace> {
        self.inner
            .run(py, |r| Ok(CacheNamespace::new(r.namespace())))
    }

    fn describe_cache(&self, py: Python<'_>, key: LocusKey) -> PyResult<Vec<CacheKey>> {
        self.inner.run(py, |r| {
            let mut keys = Vec::new();
            with_locus(&key, |locus| r.describe(locus, &mut keys))?;
            Ok(keys.into_iter().map(CacheKey::new).collect())
        })
    }

    fn prefetch_cache(&self, py: Python<'_>, key: LocusKey) -> PyResult<()> {
        self.inner
            .run(py, |r| Ok(with_locus(&key, |locus| r.prefetch(locus))?))
    }

    fn clone_context(&self, py: Python<'_>) -> PyResult<Self> {
        let r = self.inner.run(py, |r| Ok(r.try_clone()?))?;
        Ok(Self {
            inner: Shared::new(r),
            chroms: self.chroms.clone(),
            lengths: self.lengths.clone(),
        })
    }

    /// Writes the bases of `key` into `out`, 5' to 3' of its strand.
    fn sequence(
        &self,
        py: Python<'_>,
        key: LocusKey,
        mut out: PyReadwriteArray1<'_, u8>,
    ) -> PyResult<()> {
        let out = out
            .as_slice_mut()
            .map_err(|_| PyValueError::new_err("out must be a contiguous uint8 array"))?;
        self.inner.run(py, |r| -> Done<()> {
            Ok(with_locus(&key, |locus| r.sequence(locus, out))?)
        })
    }

    /// The bases of `key` in new NumPy storage.
    fn fetch<'py>(&self, py: Python<'py>, key: LocusKey) -> PyResult<Bound<'py, PyArray1<u8>>> {
        let length = self.checked_length(&key)?;
        // Initialized, privately allocated storage permits a safe mutable
        // slice and preserves NumPy ownership/resizing of default outputs.
        let out = zeros(py, length)?;
        {
            let mut borrow = out.readwrite();
            let bytes = borrow
                .as_slice_mut()
                .expect("new NumPy array is contiguous");
            self.inner.run(py, |r| {
                Ok(with_locus(&key, |locus| r.sequence(locus, bytes))?)
            })?;
        }
        Ok(out)
    }

    /// Validates a key before the Python caller allocates a large output.
    fn sequence_length(&self, key: LocusKey) -> PyResult<usize> {
        self.checked_length(&key)
    }

    fn close(&self, py: Python<'_>) {
        drop(self.inner.take(py));
    }
}

/// FASTA being written, with its `.fai` (and `.gzi`) built as it goes.
#[pyclass(frozen, module = "omicsweave._native")]
pub struct FastaWriter {
    inner: Shared<omicsweave_fasta::FastaWriter<BufWriter<File>>>,
    path: PathBuf,
}

#[pymethods]
impl FastaWriter {
    /// `level` compresses as BGZF; `None` writes plain text.
    #[new]
    #[pyo3(signature = (path, line_width=60, level=None))]
    fn new(
        py: Python<'_>,
        path: PathBuf,
        line_width: usize,
        level: Option<i32>,
    ) -> PyResult<FastaWriter> {
        py.detach(|| -> Result<FastaWriter> {
            let out = BufWriter::new(File::create(&path)?);
            let writer = match level {
                Some(level) => omicsweave_fasta::FastaWriter::bgzf(out, line_width, level)?,
                None => omicsweave_fasta::FastaWriter::new(out, line_width),
            };
            Ok(FastaWriter {
                inner: Shared::new_writer(writer),
                path,
            })
        })
        .or_raise()
    }

    fn add(&self, py: Python<'_>, name: &str, seq: &[u8]) -> PyResult<()> {
        self.inner
            .run(py, |w| -> Done<()> { Ok(w.add(name, seq)?) })
    }

    /// Ends the file and saves its indexes beside it; returns their paths.
    fn close(&self, py: Python<'_>) -> PyResult<Vec<String>> {
        let Some(writer) = self.inner.take(py) else {
            return Ok(Vec::new());
        };
        py.detach(|| -> Result<Vec<String>> {
            let (out, fai, gzi) = writer.finish()?;
            out.into_inner().map_err(|e| e.into_error())?.sync_all()?;
            fai.save(&self.path)?;
            let mut saved = vec![path_str(&fai_path(&self.path))];
            if let Some(gzi) = gzi {
                gzi.save(&self.path)?;
                saved.push(path_str(&gzi_path(&self.path)));
            }
            Ok(saved)
        })
        .or_raise()
    }
}
