//! Alignments (BAM, CRAM): records by region as columns, and building the
//! BAI / CSI / CRAI index.

use std::path::{Path, PathBuf};
use std::sync::Arc;

use numpy::{IntoPyArray, PyArray1};
use omicsweave_bam::{BamReader, BinIndex, IndexChoice, build_index_file, locate_index};
use omicsweave_core::{
    AlignmentEntry, AlignmentFields, AlignmentRead, AlignmentsOut, CigarColumn, Error, FileSource,
    ModificationColumn, ModificationEntry, Region, Result, StringColumn,
};
use omicsweave_cram::{CRAM_MAGIC, CramHeader, CramIndex, CramOptions, CramReader};
use pyo3::prelude::*;
use pyo3::types::PyBytes;

use crate::cache::{self, PyBlockCache};
use crate::error::{Done, OrRaise};
use crate::sequence::{SequenceReader, Sequences};
use crate::shared::{Shared, head, path_str, threads_or_available};

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
    decode: &Decode,
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
        verify_checksums: decode.verify_checksums,
        decode_md: decode.decode_md,
        // htslib names reads stored without a name after the file.
        name_prefix: match &decode.name_prefix {
            Some(p) => p.as_bytes().to_vec(),
            None => path
                .file_name()
                .map(|n| n.to_string_lossy().as_bytes().to_vec())
                .unwrap_or_default(),
        },
    };
    Ok((CramReader::open(source, index, reference, options)?, at))
}

/// How CRAM records are decoded.
struct Decode {
    verify_checksums: bool,
    decode_md: bool,
    name_prefix: Option<String>,
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
            Reads::Cram(r) => Reads::Cram(Box::new(
                r.try_clone_with(FileSource::try_clone, Sequences::try_clone)?,
            )),
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
    modifications: Option<Modifications>,
}

/// Base modifications of every record, column by column.
struct Modifications {
    offsets: Vec<i64>,
    ref_pos: Vec<i64>,
    qpos: Vec<i32>,
    code: Vec<i32>,
    canonical: Vec<u8>,
    probability: Vec<u8>,
    flags: Vec<u8>,
}

impl Modifications {
    fn of(offsets: Vec<i64>, rows: &[ModificationEntry]) -> Modifications {
        let n = rows.len();
        // Separate vectorized passes are faster for small batches; a single
        // pass avoids repeatedly streaming large record arrays through cache.
        if n < 16_384 {
            return Modifications {
                offsets,
                ref_pos: rows.iter().map(|m| m.ref_pos).collect(),
                qpos: rows.iter().map(|m| m.qpos).collect(),
                code: rows.iter().map(|m| m.code).collect(),
                canonical: rows.iter().map(|m| m.canonical).collect(),
                probability: rows.iter().map(|m| m.probability).collect(),
                flags: rows.iter().map(|m| m.flags).collect(),
            };
        }
        let mut columns = Modifications {
            offsets,
            ref_pos: Vec::with_capacity(n),
            qpos: Vec::with_capacity(n),
            code: Vec::with_capacity(n),
            canonical: Vec::with_capacity(n),
            probability: Vec::with_capacity(n),
            flags: Vec::with_capacity(n),
        };
        for (i, m) in rows.iter().enumerate() {
            // SAFETY: every vector owns a distinct allocation with capacity
            // at least n, and i < rows.len() == n. No vector is reallocated
            // during the loop. Typed writes preserve each column's element
            // type and initialize slot i once; lengths stay zero until all
            // slots of every column are written.
            unsafe {
                columns.ref_pos.as_mut_ptr().add(i).write(m.ref_pos);
                columns.qpos.as_mut_ptr().add(i).write(m.qpos);
                columns.code.as_mut_ptr().add(i).write(m.code);
                columns.canonical.as_mut_ptr().add(i).write(m.canonical);
                columns.probability.as_mut_ptr().add(i).write(m.probability);
                columns.flags.as_mut_ptr().add(i).write(m.flags);
            }
        }
        // SAFETY: the loop initialized exactly n primitive values in each
        // allocation, within its capacity, and retained no pointers to them.
        unsafe {
            columns.ref_pos.set_len(n);
            columns.qpos.set_len(n);
            columns.code.set_len(n);
            columns.canonical.set_len(n);
            columns.probability.set_len(n);
            columns.flags.set_len(n);
        }
        columns
    }
}

struct Alignments {
    reader: Reads,
    entries: Vec<AlignmentEntry>,
    mods: Vec<ModificationEntry>,
}

fn strings(wanted: bool, rows: usize, bytes: usize) -> Strings {
    wanted.then(|| (vec![0i64; rows + 1], vec![0u8; bytes]))
}

fn column<'a>(s: &'a mut Strings) -> Option<StringColumn<'a>> {
    s.as_mut()
        .map(|(offsets, data)| StringColumn { offsets, data })
}

impl Alignments {
    fn query(
        &mut self,
        region: Region<'_>,
        fields: AlignmentFields,
        modification_codes: &[i32],
    ) -> Result<Columns> {
        let reader = self.reader.read();
        let size = reader.query_filtered(region, fields, modification_codes)?;
        let n = size.records;
        self.entries.resize(n, AlignmentEntry::default());
        let mut names = strings(fields.contains(AlignmentFields::NAME), n, size.name_bytes);
        let mut seqs = strings(fields.contains(AlignmentFields::SEQ), n, size.seq_bytes);
        let mut quals = strings(fields.contains(AlignmentFields::QUAL), n, size.seq_bytes);
        let mut aux = strings(fields.contains(AlignmentFields::AUX), n, size.aux_bytes);
        let mut cigars = fields
            .contains(AlignmentFields::CIGAR)
            .then(|| (vec![0i64; n + 1], vec![0u32; size.cigar_ops]));
        let want_mods = fields.contains(AlignmentFields::MODIFICATIONS);
        let mut mod_offsets = want_mods.then(|| vec![0i64; n + 1]);
        self.mods
            .resize(size.modifications, ModificationEntry::default());
        reader.write(AlignmentsOut {
            entries: &mut self.entries,
            names: column(&mut names),
            seqs: column(&mut seqs),
            quals: column(&mut quals),
            cigars: cigars
                .as_mut()
                .map(|(offsets, ops)| CigarColumn { offsets, ops }),
            aux: column(&mut aux),
            modifications: mod_offsets.as_mut().map(|offsets| ModificationColumn {
                offsets,
                entries: &mut self.mods,
            }),
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
            modifications: mod_offsets.map(|o| Modifications::of(o, &self.mods)),
        })
    }
}

type PyStrings<'py> = Option<(Bound<'py, PyArray1<i64>>, Bound<'py, PyArray1<u8>>)>;

type PyModifications<'py> = (
    Bound<'py, PyArray1<i64>>,
    Bound<'py, PyArray1<i64>>,
    Bound<'py, PyArray1<i32>>,
    Bound<'py, PyArray1<i32>>,
    Bound<'py, PyArray1<u8>>,
    Bound<'py, PyArray1<u8>>,
    Bound<'py, PyArray1<u8>>,
);

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
    /// A BAM, or a CRAM with `reference` (a native `SequenceReader`, whose
    /// open file and cache the reader shares) for its slices that need one.
    /// `verify_checksums`, `decode_md` and `name_prefix` shape CRAM
    /// decoding (htslib's defaults: checks on, MD / NM generated, unnamed
    /// reads named after the file).
    #[new]
    #[pyo3(signature = (path, index=None, cache=None, reference=None, *,
                        verify_checksums=true, decode_md=true, name_prefix=None))]
    #[allow(clippy::too_many_arguments)]
    fn new(
        py: Python<'_>,
        path: PathBuf,
        index: Option<PathBuf>,
        cache: Option<PyRef<'_, PyBlockCache>>,
        reference: Option<PyRef<'_, SequenceReader>>,
        verify_checksums: bool,
        decode_md: bool,
        name_prefix: Option<String>,
    ) -> PyResult<AlignmentReader> {
        let cache = cache::handle(cache)?;
        let given = match reference {
            None => None,
            Some(seq) => Some(seq.inner.run(py, |s| Ok(s.try_clone()?))?),
        };
        let decode = Decode {
            verify_checksums,
            decode_md,
            name_prefix,
        };
        let cram_options = !verify_checksums || !decode_md || decode.name_prefix.is_some();
        py.detach(|| -> Result<AlignmentReader> {
            let is_cram = head(&FileSource::open(&path)?, 4)? == CRAM_MAGIC;
            if !is_cram {
                if given.is_some() {
                    return Err(Error::invalid(
                        "AlignmentReader",
                        "a reference is given for a BAM; only CRAM files take one",
                    ));
                }
                if cram_options {
                    return Err(Error::invalid(
                        "AlignmentReader",
                        "verify_checksums, decode_md and name_prefix are CRAM decoding options",
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
                        mods: Vec::new(),
                    }),
                });
            }
            let (mut reader, index) = open_cram(&path, index.as_deref(), given, &decode)?;
            reader.set_cache(cache);
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
                    mods: Vec::new(),
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
            Reads::Cram(r) => {
                r.set_cache(cache);
                Ok(())
            }
        })
    }

    fn clone_context(&self, py: Python<'_>) -> PyResult<Self> {
        let reader = self.inner.run(py, |a| Ok(a.reader.try_clone()?))?;
        Ok(Self {
            inner: Shared::new(Alignments {
                reader,
                entries: Vec::new(),
                mods: Vec::new(),
            }),
            chroms: self.chroms.clone(),
            text: self.text.clone(),
            index: self.index.clone(),
        })
    }

    /// The records overlapping `chrom:[start, end)`: the fixed fields as
    /// a tuple of columns, then (offsets, data) of names, sequences, qualities, CIGARs
    /// and aux fields, and (offsets, ref_pos, qpos, code, canonical,
    /// probability, flags) of base modifications, each `None` unless
    /// `fields` asks for it. Modifications keep the codes in
    /// `modification_codes` (letters as ASCII, ChEBI numbers negated), or
    /// all when it is empty.
    #[allow(clippy::type_complexity)]
    #[pyo3(signature = (chrom, start, end, fields, modification_codes=Vec::new()))]
    fn query<'py>(
        &self,
        py: Python<'py>,
        chrom: &str,
        start: u32,
        end: u32,
        fields: u8,
        modification_codes: Vec<i32>,
    ) -> PyResult<(
        (
            Bound<'py, PyArray1<i64>>,
            Bound<'py, PyArray1<i64>>,
            Bound<'py, PyArray1<u16>>,
            Bound<'py, PyArray1<u8>>,
            Bound<'py, PyArray1<i32>>,
            Bound<'py, PyArray1<i32>>,
            Bound<'py, PyArray1<i32>>,
        ),
        PyStrings<'py>,
        PyStrings<'py>,
        PyStrings<'py>,
        Option<(Bound<'py, PyArray1<i64>>, Bound<'py, PyArray1<u32>>)>,
        PyStrings<'py>,
        Option<PyModifications<'py>>,
    )> {
        let fields = AlignmentFields::from_bits(fields).or_raise()?;
        let c = self.inner.run(py, |a| -> Done<Columns> {
            Ok(a.query(Region::new(chrom, start, end), fields, &modification_codes)?)
        })?;
        Ok((
            (
                c.start.into_pyarray(py),
                c.end.into_pyarray(py),
                c.flag.into_pyarray(py),
                c.mapq.into_pyarray(py),
                c.mate_tid.into_pyarray(py),
                c.mate_start.into_pyarray(py),
                c.tlen.into_pyarray(py),
            ),
            py_strings(py, c.names),
            py_strings(py, c.seqs),
            py_strings(py, c.quals),
            c.cigars
                .map(|(o, d)| (o.into_pyarray(py), d.into_pyarray(py))),
            py_strings(py, c.aux),
            c.modifications.map(|m| {
                (
                    m.offsets.into_pyarray(py),
                    m.ref_pos.into_pyarray(py),
                    m.qpos.into_pyarray(py),
                    m.code.into_pyarray(py),
                    m.canonical.into_pyarray(py),
                    m.probability.into_pyarray(py),
                    m.flags.into_pyarray(py),
                )
            }),
        ))
    }

    fn close(&self, py: Python<'_>) {
        drop(self.inner.take(py));
    }
}

/// Builds the BAI or CSI index of a sorted BAM beside it; returns its
/// kind (`bai` or `csi`) and path.
#[pyfunction]
#[pyo3(signature = (path, kind, min_shift, threads=None))]
pub fn build_bam_index(
    py: Python<'_>,
    path: PathBuf,
    kind: &str,
    min_shift: u32,
    threads: Option<usize>,
) -> PyResult<(&'static str, String)> {
    let threads = threads_or_available(threads)?;
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

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn modification_columns_preserve_values_and_record_offsets() {
        for n in [0, 7, 16_383, 16_384, 16_385, 200_000] {
            let rows: Vec<_> = (0..n)
                .map(|i| ModificationEntry {
                    ref_pos: i as i64 * 17 - 1,
                    qpos: (i % 137) as i32,
                    code: if i % 2 == 0 { b'm' as i32 } else { -17596 },
                    canonical: b"ACGTN"[i % 5],
                    probability: (i * 73) as u8,
                    flags: (i % 8) as u8,
                    _pad: [0; 5],
                })
                .collect();
            let offsets = vec![0, n as i64 / 3, n as i64 / 3, n as i64];
            let columns = Modifications::of(offsets.clone(), &rows);
            assert_eq!(columns.offsets, offsets);
            assert_eq!(
                columns.ref_pos,
                rows.iter().map(|m| m.ref_pos).collect::<Vec<_>>()
            );
            assert_eq!(
                columns.qpos,
                rows.iter().map(|m| m.qpos).collect::<Vec<_>>()
            );
            assert_eq!(
                columns.code,
                rows.iter().map(|m| m.code).collect::<Vec<_>>()
            );
            assert_eq!(
                columns.canonical,
                rows.iter().map(|m| m.canonical).collect::<Vec<_>>()
            );
            assert_eq!(
                columns.probability,
                rows.iter().map(|m| m.probability).collect::<Vec<_>>()
            );
            assert_eq!(
                columns.flags,
                rows.iter().map(|m| m.flags).collect::<Vec<_>>()
            );
        }
    }
}
