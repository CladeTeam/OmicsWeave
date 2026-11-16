//! Gene annotation (GTF / GFF3): columnar tables, lookups, overlap queries
//! and derived interval views.
//!
//! The annotation is immutable once read, so it is shared without a lock;
//! every call that walks it releases the GIL.

use std::fs::File;
use std::io::BufWriter;
use std::path::PathBuf;
use std::sync::Arc;

use numpy::{IntoPyArray, PyArray1, PyReadonlyArray1};
use omicsweave_core::{Error, Result};
use omicsweave_gtf as gtf;
use omicsweave_gtf::{
    Aliases, AttributePolicy, Dialect, Flank, GtfLayout, Level, NONE, ReadOptions, Row, Select,
    Selection, Strand, StringPool, Tags,
};
use pyo3::exceptions::{PyRuntimeError, PyValueError};
use pyo3::prelude::*;
use pyo3::types::PyDict;

use crate::error::OrRaise;
use crate::shared::{process_id, slice};

fn check_process(process: u32) -> PyResult<()> {
    if process != process_id() {
        return Err(PyRuntimeError::new_err(
            "annotation query state cannot be reused after fork; reopen or rebuild it in this process",
        ));
    }
    Ok(())
}

type PyStrings<'py> = (Bound<'py, PyArray1<i64>>, Bound<'py, PyArray1<u8>>);

/// Per-query rows as (offsets, rows).
type BatchHits<'py> = (Bound<'py, PyArray1<u64>>, Bound<'py, PyArray1<u32>>);

/// (gene rows, transcript rows).
type Rows<'py> = (Bound<'py, PyArray1<u32>>, Bound<'py, PyArray1<u32>>);

/// The strings `ids` of `pool` as Arrow offsets and bytes (`NONE` is `""`).
fn strings(pool: &StringPool, ids: impl Iterator<Item = u32>) -> (Vec<i64>, Vec<u8>) {
    let mut offsets = vec![0i64];
    let mut data = Vec::new();
    for id in ids {
        if id != NONE {
            data.extend_from_slice(pool.bytes(id));
        }
        offsets.push(data.len() as i64);
    }
    (offsets, data)
}

fn py_strings(py: Python<'_>, s: (Vec<i64>, Vec<u8>)) -> PyStrings<'_> {
    (s.0.into_pyarray(py), s.1.into_pyarray(py))
}

fn strands(s: &[Strand]) -> Vec<i8> {
    s.iter().map(|x| x.sign()).collect()
}

fn parse_strand(v: i8) -> Strand {
    match v {
        1 => Strand::Forward,
        -1 => Strand::Reverse,
        _ => Strand::Unstranded,
    }
}

fn level(s: &str) -> PyResult<Level> {
    Level::parse(s).ok_or_else(|| {
        PyValueError::new_err(format!("level must be gene, transcript or exon, not {s:?}"))
    })
}

fn aliases(pairs: Option<Vec<(String, String)>>) -> Option<Aliases> {
    pairs.map(Aliases::from_pairs)
}

/// A dict of named numpy columns.
struct Columns<'py> {
    dict: Bound<'py, PyDict>,
}

impl<'py> Columns<'py> {
    fn new(py: Python<'py>) -> Columns<'py> {
        Columns {
            dict: PyDict::new(py),
        }
    }

    fn array<T: numpy::Element>(self, name: &str, v: Vec<T>) -> PyResult<Self> {
        let py = self.dict.py();
        self.dict.set_item(name, v.into_pyarray(py))?;
        Ok(self)
    }

    fn strings(self, name: &str, s: (Vec<i64>, Vec<u8>)) -> PyResult<Self> {
        let py = self.dict.py();
        self.dict.set_item(name, py_strings(py, s))?;
        Ok(self)
    }
}

#[pyclass(frozen, name = "Annotation", module = "omicsweave._native")]
pub struct GeneAnnotation {
    inner: Arc<gtf::Annotation>,
    process: u32,
}

/// An interval set: derived from an annotation or built from columns.
#[pyclass(frozen, module = "omicsweave._native")]
pub struct IntervalSet {
    inner: Arc<gtf::IntervalSet>,
    process: u32,
}

impl IntervalSet {
    fn wrap(set: gtf::IntervalSet) -> IntervalSet {
        IntervalSet {
            inner: Arc::new(set),
            process: std::process::id(),
        }
    }
}

#[pymethods]
impl IntervalSet {
    /// A set from columns; `seqid` numbers `chroms`, `strand` is `1`, `-1`
    /// or `0`.
    #[staticmethod]
    #[pyo3(signature = (chroms, seqid, start, end, strand, owner=None, names=None, owner_level=None))]
    #[allow(clippy::too_many_arguments)]
    fn from_columns(
        py: Python<'_>,
        chroms: Vec<String>,
        seqid: PyReadonlyArray1<'_, u32>,
        start: PyReadonlyArray1<'_, u32>,
        end: PyReadonlyArray1<'_, u32>,
        strand: PyReadonlyArray1<'_, i8>,
        owner: Option<PyReadonlyArray1<'_, u32>>,
        names: Option<Vec<String>>,
        owner_level: Option<&str>,
    ) -> PyResult<IntervalSet> {
        let seqid = slice(&seqid, "seqid")?.to_vec();
        let start = slice(&start, "start")?.to_vec();
        let end = slice(&end, "end")?.to_vec();
        let strand: Vec<Strand> = slice(&strand, "strand")?
            .iter()
            .map(|&s| parse_strand(s))
            .collect();
        let owner = owner
            .map(|o| slice(&o, "owner").map(<[u32]>::to_vec))
            .transpose()?;
        let owner_level = owner_level.map(level).transpose()?;
        py.detach(|| {
            let chroms: Vec<&str> = chroms.iter().map(String::as_str).collect();
            let names: Option<Vec<&str>> = names
                .as_ref()
                .map(|v| v.iter().map(String::as_str).collect());
            gtf::IntervalSet::from_columns(
                &chroms,
                seqid,
                start,
                end,
                strand,
                owner,
                names.as_deref(),
                owner_level,
            )
            .map(IntervalSet::wrap)
        })
        .or_raise()
    }

    fn __len__(&self) -> usize {
        self.inner.len()
    }

    fn chroms(&self) -> Vec<String> {
        (0..self.inner.seqid_count() as u32)
            .map(|s| self.inner.seqid_name(s).to_string())
            .collect()
    }

    #[getter]
    fn owner_level(&self) -> Option<&'static str> {
        self.inner.owner_level().map(Level::name)
    }

    /// seqid, start, end, strand, owner, and name (or `None` when no row
    /// has one).
    fn columns<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyDict>> {
        let s = Arc::clone(&self.inner);
        let (named, names) = py.detach(|| {
            let named = s.name.iter().any(|&n| n != NONE);
            (
                named,
                named.then(|| strings(s.pool(), s.name.iter().copied())),
            )
        });
        let c = Columns::new(py)
            .array("seqid", s.seqid.clone())?
            .array("start", s.start.clone())?
            .array("end", s.end.clone())?
            .array("strand", strands(&s.strand))?
            .array("owner", s.owner.clone())?;
        let c = match names {
            Some(n) if named => c.strings("name", n)?,
            _ => {
                c.dict.set_item("name", py.None())?;
                c
            }
        };
        Ok(c.dict)
    }

    /// Rows overlapping `chrom:[start, end)`, in start order.
    fn overlaps<'py>(
        &self,
        py: Python<'py>,
        chrom: &str,
        start: u32,
        end: u32,
    ) -> PyResult<Bound<'py, PyArray1<u32>>> {
        check_process(self.process)?;
        let s = Arc::clone(&self.inner);
        let hits = py.detach(|| {
            let mut out = Vec::new();
            if let Some(seqid) = s.seqid_of(chrom) {
                s.overlaps(seqid, start, end, &mut out);
            }
            out
        });
        Ok(hits.into_pyarray(py))
    }

    /// Overlaps of many regions: (offsets, rows), region `i`'s rows being
    /// `rows[offsets[i]:offsets[i + 1]]`.
    fn overlaps_batch<'py>(
        &self,
        py: Python<'py>,
        chroms: Vec<String>,
        starts: PyReadonlyArray1<'_, u32>,
        ends: PyReadonlyArray1<'_, u32>,
    ) -> PyResult<BatchHits<'py>> {
        check_process(self.process)?;
        let starts = slice(&starts, "starts")?;
        let ends = slice(&ends, "ends")?;
        if chroms.len() != starts.len() || ends.len() != starts.len() {
            return Err(PyValueError::new_err(
                "chroms, starts and ends differ in length",
            ));
        }
        let s = Arc::clone(&self.inner);
        let (off, ids) = py.detach(|| {
            let queries = chroms
                .iter()
                .zip(starts)
                .zip(ends)
                .map(|((c, &a), &b)| (s.seqid_of(c).unwrap_or(u32::MAX), a, b));
            let (mut off, mut ids) = (Vec::new(), Vec::new());
            s.index().overlaps_batch(queries, &mut off, &mut ids);
            (off, ids)
        });
        Ok((off.into_pyarray(py), ids.into_pyarray(py)))
    }
}

fn selection(
    genes: Option<PyReadonlyArray1<'_, u32>>,
    transcripts: Option<PyReadonlyArray1<'_, u32>>,
) -> PyResult<Option<Selection>> {
    match (genes, transcripts) {
        (Some(g), Some(t)) => Ok(Some(Selection {
            genes: slice(&g, "genes")?.to_vec(),
            transcripts: slice(&t, "transcripts")?.to_vec(),
        })),
        (None, None) => Ok(None),
        _ => Err(PyValueError::new_err(
            "a selection needs genes and transcripts",
        )),
    }
}

fn check_selection(ann: &gtf::Annotation, sel: &Option<Selection>) -> Result<()> {
    if let Some(s) = sel {
        let ok = |v: &[u32], n: usize| {
            v.windows(2).all(|w| w[0] < w[1]) && v.iter().all(|&x| (x as usize) < n)
        };
        if !ok(&s.genes, ann.genes().len()) || !ok(&s.transcripts, ann.transcripts().len()) {
            return Err(Error::invalid(
                "selection",
                "rows must be ascending, distinct and in range",
            ));
        }
    }
    Ok(())
}

impl GeneAnnotation {
    /// The annotation, shared (for quantification).
    pub(crate) fn shared(&self) -> PyResult<Arc<gtf::Annotation>> {
        check_process(self.process)?;
        Ok(Arc::clone(&self.inner))
    }

    fn with<R: Send>(
        &self,
        py: Python<'_>,
        f: impl FnOnce(&gtf::Annotation) -> Result<R> + Send,
    ) -> PyResult<R> {
        check_process(self.process)?;
        let a = Arc::clone(&self.inner);
        py.detach(|| f(&a)).or_raise()
    }

    /// Positions in `exons()` of exon feature rows.
    fn exon_positions(ann: &gtf::Annotation, rows: Vec<u32>) -> Vec<u32> {
        let exons = ann.exons();
        rows.into_iter()
            .map(|f| exons.binary_search(&f).map_or(NONE, |p| p as u32))
            .collect()
    }
}

#[pymethods]
impl GeneAnnotation {
    /// Reads a GTF or GFF3 file (plain, gzip or BGZF). `attributes` is
    /// `"all"`, `"none"`, or a list of keys to keep.
    #[new]
    #[pyo3(signature = (path, attributes=None, keys=None, aliases=None, strict=false, threads=0, dialect=None))]
    #[allow(clippy::too_many_arguments)]
    fn new(
        py: Python<'_>,
        path: PathBuf,
        attributes: Option<&str>,
        keys: Option<Vec<String>>,
        aliases: Option<Vec<(String, String)>>,
        strict: bool,
        threads: usize,
        dialect: Option<&str>,
    ) -> PyResult<GeneAnnotation> {
        let attributes = match (attributes.unwrap_or("all"), keys) {
            ("all", None) => AttributePolicy::All,
            ("none", None) => AttributePolicy::None,
            ("only", Some(k)) => AttributePolicy::Only(k),
            (other, _) => {
                return Err(PyValueError::new_err(format!(
                    "attributes must be all, none or a list of keys, not {other:?}"
                )));
            }
        };
        let dialect = match dialect {
            None => None,
            Some("gtf") => Some(Dialect::Gtf),
            Some("gff3") => Some(Dialect::Gff3),
            Some(other) => {
                return Err(PyValueError::new_err(format!(
                    "dialect must be gtf or gff3, not {other:?}"
                )));
            }
        };
        let opts = ReadOptions {
            dialect,
            attributes,
            aliases: self::aliases(aliases),
            strict,
            threads,
            ..ReadOptions::default()
        };
        let ann = py
            .detach(|| gtf::Annotation::read(&path, &opts))
            .or_raise()?;
        Ok(GeneAnnotation {
            inner: Arc::new(ann),
            process: std::process::id(),
        })
    }

    #[getter]
    fn dialect(&self) -> &'static str {
        match self.inner.dialect() {
            Dialect::Gtf => "gtf",
            Dialect::Gff3 => "gff3",
        }
    }

    fn chroms(&self) -> Vec<String> {
        self.inner.seqid_names().map(str::to_string).collect()
    }

    fn heap_bytes(&self) -> usize {
        self.inner.heap_bytes()
    }

    fn genes<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyDict>> {
        let a = Arc::clone(&self.inner);
        let g = a.genes();
        let (id, name, biotype, source, ftype) = py.detach(|| {
            let p = a.pool();
            (
                strings(p, g.id.iter().copied()),
                strings(p, g.name.iter().copied()),
                strings(p, g.biotype.iter().copied()),
                strings(p, g.source.iter().copied()),
                strings(p, g.feature_type.iter().copied()),
            )
        });
        let off = &g.transcript_offsets;
        Ok(Columns::new(py)
            .array("seqid", g.seqid.clone())?
            .array("start", g.start.clone())?
            .array("end", g.end.clone())?
            .array("strand", strands(&g.strand))?
            .strings("id", id)?
            .strings("name", name)?
            .strings("biotype", biotype)?
            .strings("source", source)?
            .strings("type", ftype)?
            .array("transcript_start", off[..off.len() - 1].to_vec())?
            .array("transcript_end", off[1..].to_vec())?
            .array(
                "synthesized",
                g.record.iter().map(|&r| r == NONE).collect::<Vec<bool>>(),
            )?
            .dict)
    }

    fn transcripts<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyDict>> {
        check_process(self.process)?;
        let a = Arc::clone(&self.inner);
        let t = a.transcripts();
        let (id, name, biotype, source, ftype, tags, lengths, cds) = py.detach(|| {
            let p = a.pool();
            (
                strings(p, t.id.iter().copied()),
                strings(p, t.name.iter().copied()),
                strings(p, t.biotype.iter().copied()),
                strings(p, t.source.iter().copied()),
                strings(p, t.feature_type.iter().copied()),
                strings(p, t.tag_values.iter().copied()),
                a.transcript_lengths().to_vec(),
                a.cds_lengths().to_vec(),
            )
        });
        let off = &t.feature_offsets;
        Ok(Columns::new(py)
            .array("seqid", t.seqid.clone())?
            .array("start", t.start.clone())?
            .array("end", t.end.clone())?
            .array("strand", strands(&t.strand))?
            .array("gene", t.gene.clone())?
            .strings("id", id)?
            .strings("name", name)?
            .strings("biotype", biotype)?
            .strings("source", source)?
            .strings("type", ftype)?
            .array(
                "tag_bits",
                t.tag_bits.iter().map(|b| b.0).collect::<Vec<u32>>(),
            )?
            .strings("tag_values", tags)?
            .array(
                "tag_offsets",
                t.tag_offsets
                    .iter()
                    .map(|&o| o as i64)
                    .collect::<Vec<i64>>(),
            )?
            .array("feature_start", off[..off.len() - 1].to_vec())?
            .array("feature_end", off[1..].to_vec())?
            .array("length", lengths)?
            .array("cds_length", cds)?
            .array(
                "synthesized",
                t.record.iter().map(|&r| r == NONE).collect::<Vec<bool>>(),
            )?
            .dict)
    }

    /// Feature columns, of every row or of the exon rows only; `feature`
    /// holds the row numbers.
    #[pyo3(signature = (exons_only=false))]
    fn features<'py>(&self, py: Python<'py>, exons_only: bool) -> PyResult<Bound<'py, PyDict>> {
        check_process(self.process)?;
        let a = Arc::clone(&self.inner);
        let f = a.features();
        let (rows, ftype, seqid) = py.detach(|| {
            let rows: Vec<u32> = if exons_only {
                a.exons().to_vec()
            } else {
                (0..f.len() as u32).collect()
            };
            let ftype = strings(a.pool(), rows.iter().map(|&r| f.feature_type[r as usize]));
            let seqid: Vec<u32> = rows.iter().map(|&r| a.feature_seqid(r)).collect();
            (rows, ftype, seqid)
        });
        let pick = |v: &[u32]| rows.iter().map(|&r| v[r as usize]).collect::<Vec<u32>>();
        let strand: Vec<i8> = rows.iter().map(|&r| f.strand[r as usize].sign()).collect();
        let kind: Vec<u8> = rows.iter().map(|&r| f.kind[r as usize] as u8).collect();
        let phase: Vec<u8> = rows.iter().map(|&r| f.phase[r as usize]).collect();
        Ok(Columns::new(py)
            .array("transcript", pick(&f.transcript))?
            .array("seqid", seqid)?
            .array("kind", kind)?
            .array("start", pick(&f.start))?
            .array("end", pick(&f.end))?
            .array("strand", strand)?
            .array("phase", phase)?
            .strings("type", ftype)?
            .array("feature", rows)?
            .dict)
    }

    /// Feature rows of exons, ascending.
    fn exons<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyArray1<u32>>> {
        check_process(self.process)?;
        let a = Arc::clone(&self.inner);
        Ok(py.detach(|| a.exons().to_vec()).into_pyarray(py))
    }

    /// Counts, and problems as (kind names, lines, ids).
    fn report<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyDict>> {
        let a = &self.inner;
        let r = a.report();
        let d = PyDict::new(py);
        d.set_item("records", r.records)?;
        d.set_item("comments", r.comments)?;
        d.set_item("genes", r.genes)?;
        d.set_item("transcripts", r.transcripts)?;
        d.set_item("features", r.features)?;
        d.set_item("others", r.others)?;
        d.set_item("synthesized_genes", r.synthesized_genes)?;
        d.set_item("synthesized_transcripts", r.synthesized_transcripts)?;
        let kinds: Vec<&str> = r.problems.iter().map(|p| p.kind.name()).collect();
        let lines: Vec<u64> = r.problems.iter().map(|p| p.line).collect();
        let ids: Vec<String> = r.problems.iter().map(|p| a.str(p.id).to_string()).collect();
        d.set_item("problem_kinds", kinds)?;
        d.set_item("problem_lines", lines.into_pyarray(py))?;
        d.set_item("problem_ids", ids)?;
        Ok(d)
    }

    /// Gene rows by id (`by="id"`) or name.
    fn gene_rows<'py>(
        &self,
        py: Python<'py>,
        key: &str,
        by: &str,
    ) -> PyResult<Bound<'py, PyArray1<u32>>> {
        let rows = self.with(py, |a| {
            Ok(match by {
                "name" => a.genes_by_name(key),
                _ => a.genes_by_id(key),
            }
            .to_vec())
        })?;
        Ok(rows.into_pyarray(py))
    }

    fn transcript_rows<'py>(
        &self,
        py: Python<'py>,
        key: &str,
        by: &str,
    ) -> PyResult<Bound<'py, PyArray1<u32>>> {
        let rows = self.with(py, |a| {
            Ok(match by {
                "name" => a.transcripts_by_name(key),
                _ => a.transcripts_by_id(key),
            }
            .to_vec())
        })?;
        Ok(rows.into_pyarray(py))
    }

    /// The kept attributes of a gene, transcript or feature row.
    fn attributes(&self, level: &str, row: u32) -> PyResult<Vec<(String, String)>> {
        let a = &self.inner;
        let (r, n) = match level {
            "gene" => (Row::Gene(row), a.genes().len()),
            "transcript" => (Row::Transcript(row), a.transcripts().len()),
            "feature" => (Row::Feature(row), a.features().len()),
            _ => return Err(PyValueError::new_err(format!("bad level {level:?}"))),
        };
        if row as usize >= n {
            return Err(PyValueError::new_err(format!(
                "{level} row {row} out of range"
            )));
        }
        Ok(a.attributes(r)
            .map(|(k, v)| (k.to_string(), v.to_string()))
            .collect())
    }

    /// Rows at `level` overlapping `chrom:[start, end)`; exon rows are
    /// positions in `exons()`.
    fn overlaps<'py>(
        &self,
        py: Python<'py>,
        level: &str,
        chrom: &str,
        start: u32,
        end: u32,
    ) -> PyResult<Bound<'py, PyArray1<u32>>> {
        let lv = self::level(level)?;
        let hits = self.with(py, |a| {
            let mut out = Vec::new();
            a.overlaps(
                lv,
                omicsweave_core::Region::new(chrom, start, end),
                &mut out,
            );
            Ok(if lv == Level::Exon {
                GeneAnnotation::exon_positions(a, out)
            } else {
                out
            })
        })?;
        Ok(hits.into_pyarray(py))
    }

    /// Overlaps of many regions at `level`: (offsets, rows).
    fn overlaps_batch<'py>(
        &self,
        py: Python<'py>,
        level: &str,
        chroms: Vec<String>,
        starts: PyReadonlyArray1<'_, u32>,
        ends: PyReadonlyArray1<'_, u32>,
    ) -> PyResult<BatchHits<'py>> {
        let lv = self::level(level)?;
        let starts = slice(&starts, "starts")?;
        let ends = slice(&ends, "ends")?;
        if chroms.len() != starts.len() || ends.len() != starts.len() {
            return Err(PyValueError::new_err(
                "chroms, starts and ends differ in length",
            ));
        }
        let (off, ids) = self.with(py, |a| {
            let queries = chroms
                .iter()
                .zip(starts)
                .zip(ends)
                .map(|((c, &s), &e)| (a.seqid(c).unwrap_or(u32::MAX), s, e));
            let (mut off, mut ids) = (Vec::new(), Vec::new());
            a.interval_index(lv)
                .overlaps_batch(queries, &mut off, &mut ids);
            if lv == Level::Exon {
                ids = GeneAnnotation::exon_positions(a, ids);
            }
            Ok((off, ids))
        })?;
        Ok((off.into_pyarray(py), ids.into_pyarray(py)))
    }

    /// The rows at `level` nearest to `chrom:[start, end)` and their
    /// distance (`-1` when the chromosome has none).
    fn nearest<'py>(
        &self,
        py: Python<'py>,
        level: &str,
        chrom: &str,
        start: u32,
        end: u32,
    ) -> PyResult<(Bound<'py, PyArray1<u32>>, i64)> {
        let lv = self::level(level)?;
        let (hits, d) = self.with(py, |a| {
            let mut out = Vec::new();
            let d = a
                .seqid(chrom)
                .and_then(|s| a.interval_index(lv).nearest(s, start, end, &mut out));
            let out = if lv == Level::Exon {
                GeneAnnotation::exon_positions(a, out)
            } else {
                out
            };
            Ok((out, d.map_or(-1, |d| d as i64)))
        })?;
        Ok((hits.into_pyarray(py), d))
    }

    /// (gene rows, transcript rows) passing the filters.
    #[pyo3(signature = (gene_biotypes=None, transcript_biotypes=None, tags=None, tag_values=None, genes=None, chroms=None))]
    #[allow(clippy::too_many_arguments)]
    fn select<'py>(
        &self,
        py: Python<'py>,
        gene_biotypes: Option<Vec<String>>,
        transcript_biotypes: Option<Vec<String>>,
        tags: Option<Vec<String>>,
        tag_values: Option<Vec<String>>,
        genes: Option<Vec<String>>,
        chroms: Option<Vec<String>>,
    ) -> PyResult<Rows<'py>> {
        let mut f = Select::new();
        if let Some(v) = gene_biotypes {
            f = f.gene_biotypes(v);
        }
        if let Some(v) = transcript_biotypes {
            f = f.transcript_biotypes(v);
        }
        if let Some(v) = tags {
            let mut bits = Tags::default();
            for t in &v {
                bits |= Tags::parse(t)
                    .ok_or_else(|| PyValueError::new_err(format!("unknown tag {t:?}")))?;
            }
            f = f.tags(bits);
        }
        if let Some(v) = tag_values {
            f = f.tag_values(v);
        }
        if let Some(v) = genes {
            f = f.genes(v);
        }
        if let Some(v) = chroms {
            f = f.seqids(v);
        }
        let s = self.with(py, |a| Ok(a.select(&f)))?;
        Ok((s.genes.into_pyarray(py), s.transcripts.into_pyarray(py)))
    }

    #[pyo3(signature = (genes=None, transcripts=None))]
    fn introns(
        &self,
        py: Python<'_>,
        genes: Option<PyReadonlyArray1<'_, u32>>,
        transcripts: Option<PyReadonlyArray1<'_, u32>>,
    ) -> PyResult<IntervalSet> {
        let sel = selection(genes, transcripts)?;
        self.with(py, |a| {
            check_selection(a, &sel)?;
            Ok(IntervalSet::wrap(match &sel {
                Some(s) => a.derive().selection(s).introns(),
                None => a.introns().clone(),
            }))
        })
    }

    #[pyo3(signature = (genes=None, transcripts=None))]
    fn exonic_regions(
        &self,
        py: Python<'_>,
        genes: Option<PyReadonlyArray1<'_, u32>>,
        transcripts: Option<PyReadonlyArray1<'_, u32>>,
    ) -> PyResult<IntervalSet> {
        let sel = selection(genes, transcripts)?;
        self.with(py, |a| {
            check_selection(a, &sel)?;
            Ok(IntervalSet::wrap(match &sel {
                Some(s) => a.derive().selection(s).exonic_regions(),
                None => a.exonic_regions().clone(),
            }))
        })
    }

    /// (set, gene offsets, genes, transcript offsets, transcripts, donor,
    /// acceptor).
    #[pyo3(signature = (genes=None, transcripts=None))]
    #[allow(clippy::type_complexity)]
    fn junctions<'py>(
        &self,
        py: Python<'py>,
        genes: Option<PyReadonlyArray1<'_, u32>>,
        transcripts: Option<PyReadonlyArray1<'_, u32>>,
    ) -> PyResult<(
        IntervalSet,
        Bound<'py, PyArray1<u32>>,
        Bound<'py, PyArray1<u32>>,
        Bound<'py, PyArray1<u32>>,
        Bound<'py, PyArray1<u32>>,
        Bound<'py, PyArray1<u32>>,
        Bound<'py, PyArray1<u32>>,
    )> {
        let sel = selection(genes, transcripts)?;
        let (j, donor, acceptor) = self.with(py, |a| {
            check_selection(a, &sel)?;
            let j = match &sel {
                Some(s) => a.derive().selection(s).junctions(),
                None => a.junctions().clone(),
            };
            let n = j.len() as u32;
            let donor = (0..n).map(|i| j.donor(i)).collect::<Vec<u32>>();
            let acceptor = (0..n).map(|i| j.acceptor(i)).collect::<Vec<u32>>();
            Ok((j, donor, acceptor))
        })?;
        Ok((
            IntervalSet::wrap(j.set),
            j.genes.offsets.into_pyarray(py),
            j.genes.values.into_pyarray(py),
            j.transcripts.offsets.into_pyarray(py),
            j.transcripts.values.into_pyarray(py),
            donor.into_pyarray(py),
            acceptor.into_pyarray(py),
        ))
    }

    /// (set, part, transcript offsets, transcripts, gene offsets, genes).
    #[pyo3(signature = (aggregate=false, genes=None, transcripts=None))]
    #[allow(clippy::type_complexity)]
    fn exon_bins<'py>(
        &self,
        py: Python<'py>,
        aggregate: bool,
        genes: Option<PyReadonlyArray1<'_, u32>>,
        transcripts: Option<PyReadonlyArray1<'_, u32>>,
    ) -> PyResult<(
        IntervalSet,
        Bound<'py, PyArray1<u32>>,
        Bound<'py, PyArray1<u32>>,
        Bound<'py, PyArray1<u32>>,
        Bound<'py, PyArray1<u32>>,
        Bound<'py, PyArray1<u32>>,
    )> {
        let sel = selection(genes, transcripts)?;
        let b = self.with(py, |a| {
            check_selection(a, &sel)?;
            Ok(match &sel {
                Some(s) => a.derive().selection(s).exon_bins(aggregate),
                None => a.exon_bins(aggregate).clone(),
            })
        })?;
        Ok((
            IntervalSet::wrap(b.set),
            b.part.into_pyarray(py),
            b.transcripts.offsets.into_pyarray(py),
            b.transcripts.values.into_pyarray(py),
            b.genes.offsets.into_pyarray(py),
            b.genes.values.into_pyarray(py),
        ))
    }

    /// Windows per gene or transcript: `kind` is `promoter` (around the
    /// TSS), `tes` (around the TES) or `body` (the body extended);
    /// `chrom_lengths` clips them.
    #[pyo3(signature = (kind, level, upstream, downstream, chrom_lengths=None, genes=None, transcripts=None))]
    #[allow(clippy::too_many_arguments)]
    fn windows(
        &self,
        py: Python<'_>,
        kind: &str,
        level: &str,
        upstream: u32,
        downstream: u32,
        chrom_lengths: Option<Vec<(String, u32)>>,
        genes: Option<PyReadonlyArray1<'_, u32>>,
        transcripts: Option<PyReadonlyArray1<'_, u32>>,
    ) -> PyResult<IntervalSet> {
        let lv = self::level(level)?;
        if lv == Level::Exon {
            return Err(PyValueError::new_err("windows are per gene or transcript"));
        }
        if !matches!(kind, "promoter" | "tes" | "body") {
            return Err(PyValueError::new_err(format!("bad window kind {kind:?}")));
        }
        let sel = selection(genes, transcripts)?;
        let flank = Flank::new(upstream, downstream);
        self.with(py, |a| {
            check_selection(a, &sel)?;
            let lengths = chrom_lengths.map(|c| {
                let pairs: Vec<(&str, u32)> = c.iter().map(|(n, l)| (n.as_str(), *l)).collect();
                a.seqid_lengths(&pairs, None)
            });
            let mut d = a.derive();
            if let Some(s) = &sel {
                d = d.selection(s);
            }
            if let Some(l) = &lengths {
                d = d.clip(l);
            }
            Ok(IntervalSet::wrap(match kind {
                "promoter" => d.promoters(lv, flank),
                "tes" => d.tes_windows(lv, flank),
                _ => d.gene_bodies(lv, flank),
            }))
        })
    }

    /// For each name, the annotation seqid it refers to (`-1` for none).
    #[pyo3(signature = (names, aliases=None))]
    fn match_chroms<'py>(
        &self,
        py: Python<'py>,
        names: Vec<String>,
        aliases: Option<Vec<(String, String)>>,
    ) -> Bound<'py, PyArray1<i64>> {
        let aliases = self::aliases(aliases);
        let names: Vec<&str> = names.iter().map(String::as_str).collect();
        self.inner
            .match_seqids(&names, aliases.as_ref())
            .into_iter()
            .map(|s| s.map_or(-1, |s| s as i64))
            .collect::<Vec<i64>>()
            .into_pyarray(py)
    }

    /// The spans of `part` (`exons`, `cds`, `utr5`, `utr3`) of every
    /// transcript in `rows`, as CSR columns: `offsets` (`len(rows) + 1`),
    /// `start`, `end`; plus each transcript's seqid and strand.
    #[allow(clippy::type_complexity)]
    fn transcript_parts<'py>(
        &self,
        py: Python<'py>,
        rows: PyReadonlyArray1<'_, u32>,
        part: &str,
    ) -> PyResult<(
        Bound<'py, PyArray1<u32>>,
        Bound<'py, PyArray1<i8>>,
        Bound<'py, PyArray1<i64>>,
        Bound<'py, PyArray1<u32>>,
        Bound<'py, PyArray1<u32>>,
    )> {
        let part = gtf::TranscriptPart::parse(part).ok_or_else(|| {
            PyValueError::new_err(format!(
                "part must be 'exons', 'cds', 'utr5' or 'utr3', not {part:?}"
            ))
        })?;
        let rows = slice(&rows, "rows")?.to_vec();
        let (seqid, strand, offsets, starts, ends) = self.with(py, |a| {
            let t = a.transcripts();
            if let Some(&bad) = rows.iter().find(|&&r| r as usize >= t.len()) {
                return Err(Error::invalid(
                    "transcript rows",
                    format!("row {bad} of {} transcripts", t.len()),
                ));
            }
            let mut blocks = Vec::new();
            let mut offsets = Vec::with_capacity(rows.len() + 1);
            let (mut starts, mut ends) = (Vec::new(), Vec::new());
            offsets.push(0i64);
            for &r in &rows {
                a.part_blocks(r, part, &mut blocks);
                starts.extend(blocks.iter().map(|b| b.0));
                ends.extend(blocks.iter().map(|b| b.1));
                offsets.push(starts.len() as i64);
            }
            let seqid: Vec<u32> = rows.iter().map(|&r| t.seqid[r as usize]).collect();
            let strand: Vec<Strand> = rows.iter().map(|&r| t.strand[r as usize]).collect();
            Ok((seqid, strands(&strand), offsets, starts, ends))
        })?;
        Ok((
            seqid.into_pyarray(py),
            strand.into_pyarray(py),
            offsets.into_pyarray(py),
            starts.into_pyarray(py),
            ends.into_pyarray(py),
        ))
    }

    /// Writes GTF; `layout` is `auto`, `original` or `model`.
    fn write_gtf(&self, py: Python<'_>, path: PathBuf, layout: &str) -> PyResult<()> {
        let layout = match layout {
            "auto" => GtfLayout::Auto,
            "original" => GtfLayout::Original,
            "model" => GtfLayout::Model,
            other => return Err(PyValueError::new_err(format!("bad layout {other:?}"))),
        };
        self.with(py, |a| {
            let mut w = BufWriter::new(File::create(&path)?);
            a.write_gtf(&mut w, layout)?;
            std::io::Write::flush(&mut w)?;
            Ok(())
        })
    }

    /// Writes transcripts (all, or the rows given) as BED12.
    #[pyo3(signature = (path, transcripts=None))]
    fn write_bed12(
        &self,
        py: Python<'_>,
        path: PathBuf,
        transcripts: Option<PyReadonlyArray1<'_, u32>>,
    ) -> PyResult<()> {
        let rows = transcripts
            .map(|t| slice(&t, "transcripts").map(<[u32]>::to_vec))
            .transpose()?;
        self.with(py, |a| {
            if let Some(r) = &rows
                && r.iter().any(|&t| t as usize >= a.transcripts().len())
            {
                return Err(Error::invalid("write_bed12", "transcript row out of range"));
            }
            let mut w = BufWriter::new(File::create(&path)?);
            a.write_bed12(&mut w, rows.as_deref())?;
            std::io::Write::flush(&mut w)?;
            Ok(())
        })
    }
}
