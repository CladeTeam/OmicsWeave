//! Sequences: `.fa` (+ `.fai`), `.fa.gz` (+ `.fai`, `.gzi`) and `.2bit`;
//! and writing FASTA with its indexes.

use std::fs::File;
use std::io::BufWriter;
use std::path::{Path, PathBuf};

use numpy::PyReadwriteArray1;
use omicsweave_core::{Error, FileSource, Region, Result, SequenceRead};
use omicsweave_fasta::{Fai, FastaReader, Gzi, TWOBIT_MAGIC, TwoBitReader, fai_path, gzi_path};
use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;

use crate::error::{Done, OrRaise};
use crate::shared::{Shared, head, path_str};

pub enum Sequences {
    Fasta(FastaReader<FileSource>),
    TwoBit(TwoBitReader<FileSource>),
}

impl SequenceRead for Sequences {
    fn chrom_len(&self, chrom: &str) -> Result<u32> {
        match self {
            Sequences::Fasta(r) => r.chrom_len(chrom),
            Sequences::TwoBit(r) => r.chrom_len(chrom),
        }
    }

    fn sequence(&mut self, region: Region<'_>, out: &mut [u8]) -> Result<()> {
        match self {
            Sequences::Fasta(r) => r.sequence(region, out),
            Sequences::TwoBit(r) => r.sequence(region, out),
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

#[pyclass(frozen, module = "omicsweave._native")]
pub struct SequenceReader {
    pub inner: Shared<Sequences>,
    chroms: Vec<(String, u32)>,
}

#[pymethods]
impl SequenceReader {
    #[new]
    #[pyo3(signature = (path, fai=None, gzi=None, soft_mask=false))]
    fn new(
        py: Python<'_>,
        path: PathBuf,
        fai: Option<PathBuf>,
        gzi: Option<PathBuf>,
        soft_mask: bool,
    ) -> PyResult<SequenceReader> {
        py.detach(|| -> Result<SequenceReader> {
            let source = FileSource::open(&path)?;
            let kind = SequenceKind::of(&head(&source, 4)?);
            let reader = match kind {
                SequenceKind::TwoBit => Sequences::TwoBit(TwoBitReader::open(source, soft_mask)?),
                _ => Sequences::Fasta(open_fasta(
                    &path,
                    source,
                    kind == SequenceKind::Bgzf,
                    fai.as_deref(),
                    gzi.as_deref(),
                )?),
            };
            let chroms = match &reader {
                Sequences::Fasta(r) => r
                    .fai()
                    .entries()
                    .iter()
                    .map(|e| (e.name.clone(), e.len))
                    .collect(),
                Sequences::TwoBit(r) => r
                    .names()
                    .map(|n| Ok((n.to_owned(), r.chrom_len(n)?)))
                    .collect::<Result<_>>()?,
            };
            Ok(SequenceReader {
                inner: Shared::new(reader),
                chroms,
            })
        })
        .or_raise()
    }

    /// (name, length) of every sequence, in file order.
    fn chroms(&self) -> Vec<(String, u32)> {
        self.chroms.clone()
    }

    /// Writes the bases of `chrom:[start, end)` into `out`.
    fn sequence(
        &self,
        py: Python<'_>,
        chrom: &str,
        start: u32,
        end: u32,
        mut out: PyReadwriteArray1<'_, u8>,
    ) -> PyResult<()> {
        let out = out
            .as_slice_mut()
            .map_err(|_| PyValueError::new_err("out must be a contiguous uint8 array"))?;
        self.inner.run(py, |r| -> Done<()> {
            Ok(r.sequence(Region::new(chrom, start, end), out)?)
        })
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
                inner: Shared::new(writer),
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
