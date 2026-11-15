//! Alignments (BAM, CRAM): records by region as columns, and building the
//! BAI / CSI / CRAI index.

use std::path::{Path, PathBuf};
use std::sync::Arc;

use numpy::{IntoPyArray, PyArray1};
use omicsweave_bam::{BamReader, BinIndex, IndexChoice, build_index_file, locate_index};
use omicsweave_core::{
    AlignmentEntry, AlignmentFields, AlignmentRead, AlignmentsOut, CigarColumn, Error, FileSource,
    Region, Result, StringColumn,
};
use omicsweave_cram::{CRAM_MAGIC, CramHeader, CramIndex, CramOptions, CramReader};
use pyo3::prelude::*;
use pyo3::types::PyBytes;
use pyo3::types::PyDict;

use crate::cache::{self, CacheKey, CacheNamespace, PyBlockCache};
use crate::error::{Done, OrRaise};
use crate::sequence::{SequenceReader, Sequences, open_sequences};
use crate::shared::{Shared, head, path_str, seek_stats};

/// The index at `given`, or else the one htslib finds beside `bam`; it
/// must exist (indexes are built explicitly, never on open).
fn open_index(bam: &Path, given: Option<&Path>) -> Result<(BinIndex, PathBuf)> {
    let path = match given {
        Some(p) => p.to_path_buf(),
        None => locate_index(bam).ok_or_else(|| Error::Io {
            kind: std::io::ErrorKind::NotFound,
            detail: format!(
                "no index of {} found (.csi or .bai beside it); build it with Bam.build_index()",
                bam.display()
            ),
        })?,
    };
    Ok((BinIndex::read_file(&path)?, path))
}

fn open_bam(path: &Path, index: Option<&Path>) -> Result<(BamReader<FileSource>, PathBuf)> {
    let (index, at) = open_index(path, index)?;
    Ok((BamReader::open(FileSource::open(path)?, index)?, at))
}

fn chroms_of(reader: &BamReader<FileSource>) -> Vec<(String, u32)> {
    reader
        .header()
        .references()
        .iter()
        .map(|r| (r.name.clone(), r.len))
        .collect()
}

type Cram = CramReader<FileSource, Sequences>;

/// The `.crai` at `given`, or else the one htslib finds beside `cram`.
fn open_cram(
    path: &Path,
    index: Option<&Path>,
    reference: Option<Sequences>,
) -> Result<(Cram, PathBuf)> {
    let source = FileSource::open(path)?;
    let header = CramHeader::read(&source)?;
    let at = match index {
        Some(p) => p.to_path_buf(),
        None => omicsweave_cram::locate_index(path).ok_or_else(|| Error::Io {
            kind: std::io::ErrorKind::NotFound,
            detail: format!(
                "no index of {} found (.crai beside it); build it with Cram.build_index()",
                path.display()
            ),
        })?,
    };
    let index = CramIndex::read_file(&at, header.references().len())?;
    let options = CramOptions {
        // htslib names reads stored without a name after the file.
        name_prefix: path
            .file_name()
            .map(|n| n.to_string_lossy().as_bytes().to_vec())
            .unwrap_or_default(),
        ..CramOptions::default()
    };
    Ok((CramReader::open(source, index, reference, options)?, at))
}

fn no_cache() -> Error {
    Error::Unsupported("block caches of CRAM files")
}

/// A BAM or a CRAM reader.
#[allow(clippy::large_enum_variant)]
enum Reads {
    Bam(BamReader<FileSource>),
    Cram(Box<Cram>),
}

impl Reads {
    fn read(&mut self) -> &mut dyn AlignmentRead {
        match self {
            Reads::Bam(r) => r,
            Reads::Cram(r) => &mut **r,
        }
    }

    fn try_clone(&self) -> Result<Reads> {
        Ok(match self {
            Reads::Bam(r) => Reads::Bam(r.try_clone()?),
            Reads::Cram(r) => {
                let reference = match r.reference() {
                    Some(s) => Some(s.try_clone()?),
                    None => None,
                };
                Reads::Cram(Box::new(r.with_parts(r.source().try_clone()?, reference)))
            }
        })
    }
}

type Strings = Option<(Vec<i64>, Vec<u8>)>;

/// One query's output, column by column.
struct Columns {
    start: Vec<i64>,
    end: Vec<i64>,
    flag: Vec<u16>,
    mapq: Vec<u8>,
    mate_tid: Vec<i32>,
    mate_start: Vec<i32>,
    tlen: Vec<i32>,
    names: Strings,
    seqs: Strings,
    quals: Strings,
    cigars: Option<(Vec<i64>, Vec<u32>)>,
    aux: Strings,
}

struct Alignments {
    reader: Reads,
    entries: Vec<AlignmentEntry>,
}

fn strings(wanted: bool, rows: usize, bytes: usize) -> Strings {
    wanted.then(|| (vec![0i64; rows + 1], vec![0u8; bytes]))
}

fn column<'a>(s: &'a mut Strings) -> Option<StringColumn<'a>> {
    s.as_mut()
        .map(|(offsets, data)| StringColumn { offsets, data })
}

impl Alignments {
    fn query(&mut self, region: Region<'_>, fields: AlignmentFields) -> Result<Columns> {
        let reader = self.reader.read();
        let size = reader.query(region, fields)?;
        let n = size.records;
        self.entries.resize(n, AlignmentEntry::default());
        let mut names = strings(fields.contains(AlignmentFields::NAME), n, size.name_bytes);
        let mut seqs = strings(fields.contains(AlignmentFields::SEQ), n, size.seq_bytes);
        let mut quals = strings(fields.contains(AlignmentFields::QUAL), n, size.seq_bytes);
        let mut aux = strings(fields.contains(AlignmentFields::AUX), n, size.aux_bytes);
        let mut cigars = fields
            .contains(AlignmentFields::CIGAR)
            .then(|| (vec![0i64; n + 1], vec![0u32; size.cigar_ops]));
        reader.write(AlignmentsOut {
            entries: &mut self.entries,
            names: column(&mut names),
            seqs: column(&mut seqs),
            quals: column(&mut quals),
            cigars: cigars
                .as_mut()
                .map(|(offsets, ops)| CigarColumn { offsets, ops }),
            aux: column(&mut aux),
        })?;
        let e = &self.entries;
        Ok(Columns {
            start: e.iter().map(|x| x.start).collect(),
            end: e.iter().map(|x| x.end).collect(),
            flag: e.iter().map(|x| x.flag).collect(),
            mapq: e.iter().map(|x| x.mapq).collect(),
            mate_tid: e.iter().map(|x| x.mate_tid).collect(),
            mate_start: e.iter().map(|x| x.mate_start).collect(),
            tlen: e.iter().map(|x| x.tlen).collect(),
            names,
            seqs,
            quals,
            cigars,
            aux,
        })
    }
}

type PyStrings<'py> = Option<(Bound<'py, PyArray1<i64>>, Bound<'py, PyArray1<u8>>)>;

fn py_strings(py: Python<'_>, s: Strings) -> PyStrings<'_> {
    s.map(|(o, d)| (o.into_pyarray(py), d.into_pyarray(py)))
}

#[pyclass(frozen, module = "omicsweave._native")]
pub struct AlignmentReader {
    inner: Shared<Alignments>,
    chroms: Arc<[(String, u32)]>,
    text: Arc<[u8]>,
    index: PathBuf,
}

#[pymethods]
impl AlignmentReader {
    /// A BAM, or a CRAM with `reference` (a native `SequenceReader`, or
    /// the path of a FASTA / bgzip FASTA / 2bit file) for its slices that
    /// need one.
    #[new]
    #[pyo3(signature = (path, index=None, cache=None, reference=None))]
    fn new(
        py: Python<'_>,
        path: PathBuf,
        index: Option<PathBuf>,
        cache: Option<PyRef<'_, PyBlockCache>>,
        reference: Option<Bound<'_, PyAny>>,
    ) -> PyResult<AlignmentReader> {
        let cache = cache::handle(cache)?;
        enum Given {
            Reader(Sequences),
            Path(PathBuf),
        }
        let given = match reference {
            None => None,
            Some(r) => match r.extract::<PyRef<'_, SequenceReader>>() {
                Ok(seq) => Some(Given::Reader(seq.inner.run(py, |s| Ok(s.try_clone()?))?)),
                Err(_) => Some(Given::Path(r.extract::<PathBuf>()?)),
            },
        };
        py.detach(|| -> Result<AlignmentReader> {
            let is_cram = head(&FileSource::open(&path)?, 4)? == CRAM_MAGIC;
            if !is_cram {
                if given.is_some() {
                    return Err(Error::invalid(
                        "AlignmentReader",
                        "a reference is given for a BAM; only CRAM files take one",
                    ));
                }
                let (mut reader, index) = open_bam(&path, index.as_deref())?;
                reader.set_cache(cache);
                return Ok(AlignmentReader {
                    chroms: chroms_of(&reader).into(),
                    text: Arc::from(reader.header().text()),
                    index,
                    inner: Shared::new(Alignments {
                        reader: Reads::Bam(reader),
                        entries: Vec::new(),
                    }),
                });
            }
            if cache.is_some() {
                return Err(no_cache());
            }
            let reference = match given {
                None => None,
                Some(Given::Reader(s)) => Some(s),
                Some(Given::Path(p)) => Some(open_sequences(&p)?),
            };
            let (reader, index) = open_cram(&path, index.as_deref(), reference)?;
            let chroms: Vec<(String, u32)> = reader
                .header()
                .references()
                .iter()
                .map(|r| (r.name.clone(), r.len))
                .collect();
            Ok(AlignmentReader {
                chroms: chroms.into(),
                text: Arc::from(reader.header().text()),
                index,
                inner: Shared::new(Alignments {
                    reader: Reads::Cram(Box::new(reader)),
                    entries: Vec::new(),
                }),
            })
        })
        .or_raise()
    }

    /// (name, length) of every reference, in reference-id order.
    fn chroms(&self) -> Vec<(String, u32)> {
        self.chroms.to_vec()
    }

    /// The SAM header text as stored.
    fn header_text<'py>(&self, py: Python<'py>) -> Bound<'py, PyBytes> {
        PyBytes::new(py, &self.text)
    }

    /// The index file read.
    #[getter]
    fn index_path(&self) -> String {
        path_str(&self.index)
    }

    fn set_cache(&self, py: Python<'_>, cache: Option<PyRef<'_, PyBlockCache>>) -> PyResult<()> {
        let cache = cache::handle(cache)?;
        self.inner.run(py, |a| match &mut a.reader {
            Reads::Bam(r) => {
                r.set_cache(cache);
                Ok(())
            }
            Reads::Cram(_) if cache.is_none() => Ok(()),
            Reads::Cram(_) => Err(no_cache().into()),
        })
    }

    fn cache_namespace(&self, py: Python<'_>) -> PyResult<CacheNamespace> {
        self.inner.run(py, |a| match &a.reader {
            Reads::Bam(r) => Ok(CacheNamespace::new(r.cache_namespace())),
            Reads::Cram(_) => Err(no_cache().into()),
        })
    }

    fn describe_cache(
        &self,
        py: Python<'_>,
        chrom: &str,
        start: u32,
        end: u32,
    ) -> PyResult<Vec<CacheKey>> {
        self.inner.run(py, |a| {
            let mut keys = Vec::new();
            let Reads::Bam(r) = &mut a.reader else {
                return Err(no_cache().into());
            };
            r.describe_cache(Region::new(chrom, start, end), &mut keys)?;
            Ok(keys.into_iter().map(CacheKey::new).collect())
        })
    }

    fn prefetch_cache(&self, py: Python<'_>, chrom: &str, start: u32, end: u32) -> PyResult<()> {
        self.inner.run(py, |a| match &mut a.reader {
            Reads::Bam(r) => Ok(r.prefetch_cache(Region::new(chrom, start, end))?),
            Reads::Cram(_) => Err(no_cache().into()),
        })
    }

    fn clone_context(&self, py: Python<'_>) -> PyResult<Self> {
        let reader = self.inner.run(py, |a| Ok(a.reader.try_clone()?))?;
        Ok(Self {
            inner: Shared::new(Alignments {
                reader,
                entries: Vec::new(),
            }),
            chroms: self.chroms.clone(),
            text: self.text.clone(),
            index: self.index.clone(),
        })
    }

    fn stream_stats<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyDict>> {
        enum Stats {
            Bam(omicsweave_bgzf::SeekStats),
            Cram(omicsweave_cram::CramStats),
        }
        let stats = self.inner.run(py, |a| {
            Ok(match &a.reader {
                Reads::Bam(r) => Stats::Bam(r.stats()),
                Reads::Cram(r) => Stats::Cram(r.stats()),
            })
        })?;
        match stats {
            Stats::Bam(s) => seek_stats(py, s),
            Stats::Cram(s) => {
                let d = PyDict::new(py);
                d.set_item("containers", s.containers)?;
                d.set_item("slices_decoded", s.slices_decoded)?;
                d.set_item("records_decoded", s.records_decoded)?;
                d.set_item("bytes_read", s.bytes_read)?;
                Ok(d)
            }
        }
    }

    fn reset_stream_stats(&self, py: Python<'_>) -> PyResult<()> {
        self.inner.run(py, |a| {
            match &mut a.reader {
                Reads::Bam(r) => r.reset_stats(),
                Reads::Cram(r) => r.reset_stats(),
            }
            Ok(())
        })
    }

    /// The records overlapping `chrom:[start, end)`: the fixed fields as
    /// columns, then (offsets, data) of names, sequences, qualities, CIGARs
    /// and aux fields, each `None` unless `fields` asks for it.
    #[allow(clippy::type_complexity)]
    fn query<'py>(
        &self,
        py: Python<'py>,
        chrom: &str,
        start: u32,
        end: u32,
        fields: u8,
    ) -> PyResult<(
        Bound<'py, PyArray1<i64>>,
        Bound<'py, PyArray1<i64>>,
        Bound<'py, PyArray1<u16>>,
        Bound<'py, PyArray1<u8>>,
        Bound<'py, PyArray1<i32>>,
        Bound<'py, PyArray1<i32>>,
        Bound<'py, PyArray1<i32>>,
        PyStrings<'py>,
        PyStrings<'py>,
        PyStrings<'py>,
        Option<(Bound<'py, PyArray1<i64>>, Bound<'py, PyArray1<u32>>)>,
        PyStrings<'py>,
    )> {
        let fields = AlignmentFields::from_bits(fields).or_raise()?;
        let c = self.inner.run(py, |a| -> Done<Columns> {
            Ok(a.query(Region::new(chrom, start, end), fields)?)
        })?;
        Ok((
            c.start.into_pyarray(py),
            c.end.into_pyarray(py),
            c.flag.into_pyarray(py),
            c.mapq.into_pyarray(py),
            c.mate_tid.into_pyarray(py),
            c.mate_start.into_pyarray(py),
            c.tlen.into_pyarray(py),
            py_strings(py, c.names),
            py_strings(py, c.seqs),
            py_strings(py, c.quals),
            c.cigars
                .map(|(o, d)| (o.into_pyarray(py), d.into_pyarray(py))),
            py_strings(py, c.aux),
        ))
    }

    fn close(&self, py: Python<'_>) {
        drop(self.inner.take(py));
    }
}

/// Builds the BAI or CSI index of a sorted BAM beside it; returns its
/// kind (`bai` or `csi`) and path.
#[pyfunction]
#[pyo3(signature = (path, kind, min_shift, threads))]
pub fn build_bam_index(
    py: Python<'_>,
    path: PathBuf,
    kind: &str,
    min_shift: u32,
    threads: usize,
) -> PyResult<(&'static str, String)> {
    py.detach(|| -> Result<(&'static str, String)> {
        let choice = IndexChoice::from_name(kind, min_shift)?;
        let (kind, at) = build_index_file(&path, choice, threads)?;
        Ok((kind.name(), path_str(&at)))
    })
    .or_raise()
}

/// Builds the `.crai` of a CRAM beside it, as `samtools index` does;
/// returns its path.
#[pyfunction]
pub fn build_cram_index(py: Python<'_>, path: PathBuf) -> PyResult<String> {
    py.detach(|| -> Result<String> { Ok(path_str(&omicsweave_cram::build_index_file(&path)?)) })
        .or_raise()
}
