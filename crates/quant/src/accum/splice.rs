//! Splice junctions: every intron the reads cross, with STAR's `SJ.out.tab`
//! statistics, and optionally a junction × row count matrix.
//!
//! Definitions (STAR `SJ.out.tab`):
//! - A **junction** is the reference interval of an `N` operation.
//! - A read pair crossing a junction counts once, even when both mates
//!   cross it (with `pairs`; mates are paired as for read-pair counting,
//!   and a mate whose pair cannot be completed counts alone); every
//!   alignment of a multi-mapper counts.
//! - **Unique** / **multi**: crossing reads (pairs) with `NH` 1 / `NH` > 1.
//! - **Overhang** of a read at a junction: the shorter of the two aligned
//!   segments next to it, a segment being a run of `M` / `=` / `X` not
//!   interrupted by an insertion, deletion, clip or intron; the table
//!   keeps the largest of any read.
//! - **Motif** from the reference's first two and last two intron bases:
//!   1 `GT/AG`, 2 `CT/AC`, 3 `GC/AG`, 4 `CT/GC`, 5 `AT/AC`, 6 `GT/AT`,
//!   0 any other (and without a reference). **Strand**: `+` for odd
//!   motifs, `-` for even ones; for non-canonical junctions the
//!   annotation's strand when annotated, else undefined.
//! - **Annotated**: an intron of a transcript of the task's annotation.
//!
//! The matrix has a column per junction of the table, in table order, and
//! counts reads (pairs) or UMIs per row key that also pass its own filter.

use std::io::Read;
use std::path::PathBuf;
use std::sync::Arc;

use omicsweave_convert::stage::{Sorter, Spill, read_or_end};
use omicsweave_core::{Error, Result};
use omicsweave_gtf::{Annotation, Strand};
use smallvec::SmallVec;

use super::counts::{
    Dedup, Entry, Fold, GlobalRows, LocalRows, Reduce, RowKey, Rows, collapse as collapse_cells,
};
use super::coverage::Tally;
use super::{Accumulator, BatchInfo, Delivery, MemoryNeed, Needs, PlanInfo, RunCtx};
use crate::model::{BarcodeTable, Cells, CountMatrix, FeatureTable, Junction, JunctionTable};
use crate::prelude::{Filter, MateKey, Mates, Pairs, ReadView, expects_mate, umi_key};
use crate::reference::Reference;
use crate::sink::json::Json;

/// The junction × row matrix of a [`Junctions`] product.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct JunctionMatrix {
    pub rows: RowKey,
    /// Only these row keys, one per line, rows in this order.
    pub whitelist: Option<PathBuf>,
    pub umi_tag: [u8; 2],
    /// [`Dedup::None`] (reads or pairs) or [`Dedup::Umi`].
    pub dedup: Dedup,
    /// Which crossing reads the matrix counts (beyond the table's filter).
    pub filter: Filter,
}

impl JunctionMatrix {
    /// STARsolo `SJ`: UMIs per `CB`, unique alignments, `0x704` excluded.
    pub fn single_cell() -> JunctionMatrix {
        JunctionMatrix {
            rows: RowKey::Tag(*b"CB"),
            whitelist: None,
            umi_tag: *b"UB",
            dedup: Dedup::Umi,
            filter: Filter {
                unique: true,
                ..Filter::new()
            },
        }
    }

    /// Reads (pairs) per input file, unique alignments.
    pub fn bulk() -> JunctionMatrix {
        JunctionMatrix {
            rows: RowKey::File,
            whitelist: None,
            umi_tag: *b"UB",
            dedup: Dedup::None,
            filter: Filter {
                require: 0,
                exclude: omicsweave_bam::flag::UNMAPPED,
                min_mapq: 0,
                unique: true,
            },
        }
    }
}

/// Splice junctions as `<name>/junctions.tab` (the columns of STAR's
/// `SJ.out.tab`), and with `matrix` a 10x directory over them in `<name>/`.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Junctions {
    pub name: String,
    pub filter: Filter,
    /// Count a pair crossing a junction once.
    pub pairs: bool,
    pub matrix: Option<JunctionMatrix>,
}

impl Default for Junctions {
    /// Every mapped alignment that is not QC-failed or supplementary
    /// (STAR's chimeric segments), duplicates and secondary alignments
    /// included; pairs counted once; no matrix.
    fn default() -> Junctions {
        use omicsweave_bam::flag;
        Junctions {
            name: "junctions".into(),
            filter: Filter {
                require: 0,
                exclude: flag::UNMAPPED | flag::QC_FAIL | flag::SUPPLEMENTARY,
                min_mapq: 0,
                unique: false,
            },
            pairs: true,
            matrix: None,
        }
    }
}

/// One junction crossing: the table's columns before folding.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) struct JEntry {
    tid: u32,
    start: u32,
    end: u32,
    unique: u32,
    multi: u32,
    overhang: u32,
}

impl JEntry {
    #[inline]
    fn key(&self) -> (u32, u32, u32) {
        (self.tid, self.start, self.end)
    }

    #[inline]
    fn add(&mut self, o: &JEntry) {
        self.unique += o.unique;
        self.multi += o.multi;
        self.overhang = self.overhang.max(o.overhang);
    }
}

fn u32_at(b: &[u8], i: usize) -> u32 {
    u32::from_le_bytes(b[i..i + 4].try_into().unwrap())
}

impl Spill for JEntry {
    type Key = (u32, u32, u32);
    fn key(&self) -> Self::Key {
        JEntry::key(self)
    }
    fn mem_size(&self) -> usize {
        0
    }
    fn encode(&self, out: &mut Vec<u8>) {
        for v in [
            self.tid,
            self.start,
            self.end,
            self.unique,
            self.multi,
            self.overhang,
        ] {
            out.extend_from_slice(&v.to_le_bytes());
        }
    }
    fn decode(input: &mut impl Read) -> Result<Option<JEntry>> {
        let mut b = [0u8; 24];
        if !read_or_end(input, &mut b)? {
            return Ok(None);
        }
        Ok(Some(JEntry {
            tid: u32_at(&b, 0),
            start: u32_at(&b, 4),
            end: u32_at(&b, 8),
            unique: u32_at(&b, 12),
            multi: u32_at(&b, 16),
            overhang: u32_at(&b, 20),
        }))
    }
}

/// One matrix observation of a junction, keyed by the junction first so
/// that columns can be numbered once every junction is known.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) struct MEntry {
    tid: u32,
    start: u32,
    end: u32,
    row: u32,
    sub: u64,
    value: u32,
}

impl MEntry {
    #[inline]
    fn key(&self) -> (u32, u32, u32, u32, u64) {
        (self.tid, self.start, self.end, self.row, self.sub)
    }
}

impl Spill for MEntry {
    type Key = (u32, u32, u32, u32, u64);
    fn key(&self) -> Self::Key {
        MEntry::key(self)
    }
    fn mem_size(&self) -> usize {
        0
    }
    fn encode(&self, out: &mut Vec<u8>) {
        for v in [self.tid, self.start, self.end, self.row] {
            out.extend_from_slice(&v.to_le_bytes());
        }
        out.extend_from_slice(&self.sub.to_le_bytes());
        out.extend_from_slice(&self.value.to_le_bytes());
    }
    fn decode(input: &mut impl Read) -> Result<Option<MEntry>> {
        let mut b = [0u8; 28];
        if !read_or_end(input, &mut b)? {
            return Ok(None);
        }
        Ok(Some(MEntry {
            tid: u32_at(&b, 0),
            start: u32_at(&b, 4),
            end: u32_at(&b, 8),
            row: u32_at(&b, 12),
            sub: u64::from_le_bytes(b[16..24].try_into().unwrap()),
            value: u32_at(&b, 24),
        }))
    }
}

/// The junctions of a read: `(start, end, overhang)`, in CIGAR order.
pub fn read_junctions(r: &ReadView<'_>) -> SmallVec<[(u32, u32, u32); 2]> {
    use omicsweave_bam::op;
    let mut segments: SmallVec<[u32; 8]> = SmallVec::new();
    // (start, end, index of the segment before it).
    let mut found: SmallVec<[(u32, u32, usize); 2]> = SmallVec::new();
    let mut run = 0u32;
    let mut at = r.pos();
    for (o, len) in r.rec.cigar() {
        match o {
            op::MATCH | op::EQUAL | op::DIFF => {
                run += len;
                at = at.saturating_add(len);
            }
            op::REF_SKIP => {
                segments.push(run);
                found.push((at, at.saturating_add(len), segments.len() - 1));
                run = 0;
                at = at.saturating_add(len);
            }
            op::DEL => {
                segments.push(run);
                run = 0;
                at = at.saturating_add(len);
            }
            op::INS | op::SOFT_CLIP => {
                segments.push(run);
                run = 0;
            }
            _ => {}
        }
    }
    segments.push(run);
    found
        .into_iter()
        .map(|(s, e, i)| (s, e, segments[i].min(segments[i + 1])))
        .collect()
}

/// The motif of an intron from its first and last two bases.
pub fn motif(donor: &[u8], acceptor: &[u8]) -> u8 {
    match (donor, acceptor) {
        (b"GT", b"AG") => 1,
        (b"CT", b"AC") => 2,
        (b"GC", b"AG") => 3,
        (b"CT", b"GC") => 4,
        (b"AT", b"AC") => 5,
        (b"GT", b"AT") => 6,
        _ => 0,
    }
}

/// A read's crossings waiting for its mate.
#[derive(Debug, Clone)]
pub(crate) struct JHalf {
    /// `(tid, start, end, overhang)`.
    junctions: SmallVec<[(u32, u32, u32, u32); 2]>,
    multi: bool,
    /// Matrix row and UMI, when the read counts in the matrix.
    row: Option<(u32, u64)>,
}

impl JHalf {
    fn join(mut self, o: JHalf) -> JHalf {
        for j in o.junctions {
            match self
                .junctions
                .iter_mut()
                .find(|x| (x.0, x.1, x.2) == (j.0, j.1, j.2))
            {
                Some(x) => x.3 = x.3.max(j.3),
                None => self.junctions.push(j),
            }
        }
        self.multi |= o.multi;
        self.row = self.row.or(o.row);
        self
    }
}

/// A [`Junctions`] product resolved against a task.
pub struct JunctionStage {
    pub(crate) spec: Junctions,
    pub(crate) rows: Option<Rows>,
    pub(crate) reduce: Reduce,
    /// The annotation and, by tid, its seqid.
    pub(crate) annotation: Option<(Arc<Annotation>, Vec<Option<u32>>)>,
    pub(crate) reference: Option<Arc<Reference>>,
    pub(crate) table_bytes: u64,
}

pub struct JunctionBatch {
    table: Vec<JEntry>,
    matrix: Vec<MEntry>,
    local: LocalRows,
    pairs: Pairs<JHalf>,
    dangling: Vec<(MateKey, JHalf)>,
    tally: Tally,
    spliced: u64,
}

pub struct JunctionTotal<'p> {
    table: Sorter<'p, JEntry>,
    matrix: Sorter<'p, MEntry>,
    rows: GlobalRows,
    mates: Mates<JHalf>,
    tally: Tally,
    spliced: u64,
    scratch: (Vec<JEntry>, Vec<MEntry>),
}

impl JunctionStage {
    fn emit(&self, h: &JHalf, table: &mut Vec<JEntry>, matrix: &mut Vec<MEntry>) {
        for &(tid, start, end, overhang) in &h.junctions {
            table.push(JEntry {
                tid,
                start,
                end,
                unique: !h.multi as u32,
                multi: h.multi as u32,
                overhang,
            });
            if let Some((row, sub)) = h.row {
                matrix.push(MEntry {
                    tid,
                    start,
                    end,
                    row,
                    sub,
                    value: 1,
                });
            }
        }
    }

    /// The matrix row and UMI of `r`, when it counts in the matrix.
    fn row_of(&self, r: &ReadView<'_>, local: &mut LocalRows) -> Option<(u32, u64)> {
        let m = self.spec.matrix.as_ref()?;
        if m.filter.reject(r).is_some() {
            return None;
        }
        let row = self.rows.as_ref()?.of(r, local)?;
        let sub = match m.dedup {
            Dedup::Umi => umi_key(r.tag_bytes(m.umi_tag)?),
            _ => 0,
        };
        Some((row, sub))
    }

    fn annotate(
        &self,
        j: &mut Junction,
        lookup: &mut Option<crate::reference::Lookup<'_>>,
    ) -> Result<()> {
        let mut ann_strand = 0;
        if let Some((ann, seqids)) = &self.annotation
            && let Some(Some(seqid)) = seqids.get(j.tid as usize)
        {
            let js = ann.junctions();
            if let Some(row) = js.find_any(*seqid, j.start, j.end) {
                j.annotated = true;
                ann_strand = match js.set.strand[row as usize] {
                    Strand::Forward => 1,
                    Strand::Reverse => 2,
                    _ => 0,
                };
            }
        }
        if let Some(l) = lookup
            && j.end - j.start >= 4
            && let (Some(d), Some(a)) = (
                l.bases(j.tid, j.start, j.start + 2)?,
                l.bases(j.tid, j.end - 2, j.end)?,
            )
        {
            j.motif = motif(&d, &a);
        }
        j.strand = match j.motif {
            0 => ann_strand,
            m if m % 2 == 1 => 1,
            _ => 2,
        };
        Ok(())
    }
}

/// Sorts and folds a batch's table entries.
fn collapse_table(v: &mut Vec<JEntry>) {
    v.sort_unstable_by_key(JEntry::key);
    let mut w = 0;
    for i in 0..v.len() {
        if w > 0 && v[w - 1].key() == v[i].key() {
            let x = v[i];
            v[w - 1].add(&x);
        } else {
            v[w] = v[i];
            w += 1;
        }
    }
    v.truncate(w);
}

fn collapse_matrix(v: &mut Vec<MEntry>, reduce: Reduce) {
    v.sort_unstable_by_key(MEntry::key);
    let mut w = 0;
    for i in 0..v.len() {
        if w > 0 && v[w - 1].key() == v[i].key() {
            v[w - 1].value = reduce.combine(v[w - 1].value, v[i].value);
        } else {
            v[w] = v[i];
            w += 1;
        }
    }
    v.truncate(w);
}

/// The table, the matrix (with no junction at all, none), the filter's
/// tally, spliced records and spilled bytes.
type Products = (JunctionTable, Option<CountMatrix>, Tally, u64, u64);

impl Accumulator for JunctionStage {
    type Batch = JunctionBatch;
    type Total<'p> = JunctionTotal<'p>;
    type Product = Products;

    fn name(&self) -> &str {
        &self.spec.name
    }

    fn needs(&self) -> Needs {
        let mut n = Needs::CIGAR | Needs::TAGS;
        if self.spec.pairs {
            n |= Needs::MATES;
        }
        n
    }

    fn estimate(&self, info: &PlanInfo) -> MemoryNeed {
        MemoryNeed {
            fixed: self.table_bytes,
            per_batch: info.batch_bytes / 16,
            min_buffer: 2 << 20,
            weight: 2,
        }
    }

    fn batch(&self, _info: BatchInfo) -> JunctionBatch {
        JunctionBatch {
            table: Vec::new(),
            matrix: Vec::new(),
            local: LocalRows::default(),
            pairs: Pairs::default(),
            dangling: Vec::new(),
            tally: Tally::default(),
            spliced: 0,
        }
    }

    fn observe(&self, b: &mut JunctionBatch, r: &ReadView<'_>) -> Result<()> {
        b.tally.seen += 1;
        if let Some(why) = self.spec.filter.reject(r) {
            b.tally.reject(why);
            return Ok(());
        }
        let found = read_junctions(r);
        if found.is_empty() {
            return Ok(());
        }
        b.spliced += 1;
        let tid = r.tid();
        let half = JHalf {
            junctions: found.iter().map(|&(s, e, o)| (tid, s, e, o)).collect(),
            multi: r.nh() > 1,
            row: self.row_of(r, &mut b.local),
        };
        if self.spec.pairs && expects_mate(r) {
            if let Some((x, y)) = b.pairs.offer(r, half) {
                self.emit(&x.join(y), &mut b.table, &mut b.matrix);
            }
        } else {
            self.emit(&half, &mut b.table, &mut b.matrix);
        }
        Ok(())
    }

    fn seal(&self, b: &mut JunctionBatch) -> Result<()> {
        collapse_table(&mut b.table);
        collapse_matrix(&mut b.matrix, self.reduce);
        b.dangling = std::mem::take(&mut b.pairs).into_dangling();
        Ok(())
    }

    fn total<'p>(&'p self, rc: &RunCtx<'p>) -> Result<JunctionTotal<'p>> {
        let share = if self.spec.matrix.is_some() {
            rc.buffer / 2
        } else {
            rc.buffer
        };
        Ok(JunctionTotal {
            table: Sorter::new(rc.pool, share, &rc.tmp),
            matrix: Sorter::new(rc.pool, share, &rc.tmp),
            rows: GlobalRows::default(),
            mates: Mates::new(std::mem::size_of::<JHalf>()),
            tally: Tally::default(),
            spliced: 0,
            scratch: (Vec::new(), Vec::new()),
        })
    }

    fn merge<'p>(
        &'p self,
        t: &mut JunctionTotal<'p>,
        b: JunctionBatch,
        frontier: u64,
        rc: &RunCtx<'p>,
    ) -> Result<()> {
        t.tally.add(&b.tally);
        t.spliced += b.spliced;
        let JunctionBatch {
            table,
            mut matrix,
            local,
            dangling,
            ..
        } = b;
        let map = match &self.rows {
            Some(rows) if rows.local() => {
                let map = t.rows.remap(local);
                for e in &mut matrix {
                    e.row = map[e.row as usize];
                }
                Some(map)
            }
            _ => None,
        };
        let (out_t, out_m) = &mut t.scratch;
        out_t.clear();
        out_m.clear();
        for (key, mut half) in dangling {
            if let (Some(m), Some((row, _))) = (&map, &mut half.row) {
                *row = m[*row as usize];
            }
            if let Some((x, y)) = t.mates.offer(key, half) {
                self.emit(&x.join(y), out_t, out_m);
            }
        }
        t.mates.expire(frontier, |h| self.emit(&h, out_t, out_m));
        if t.rows.bytes() + t.mates.bytes() > rc.buffer.max(64 << 20) {
            return Err(Error::invalid(
                "budget",
                format!(
                    "{}: {} row keys and {} mates waiting for their pair outgrow the memory \
                     budget; give a whitelist or raise the memory",
                    self.spec.name,
                    t.rows.names.len(),
                    t.mates.len()
                ),
            ));
        }
        let mut table = table;
        table.append(out_t);
        matrix.append(out_m);
        t.table.extend(table)?;
        t.matrix.extend(matrix)
    }

    fn end_input<'p>(
        &'p self,
        t: &mut JunctionTotal<'p>,
        _input: u32,
        _rc: &RunCtx<'p>,
    ) -> Result<()> {
        let (mut table, mut matrix) = (Vec::new(), Vec::new());
        t.mates.drain(|h| self.emit(&h, &mut table, &mut matrix));
        t.table.extend(table)?;
        t.matrix.extend(matrix)
    }

    fn finish<'p>(&'p self, t: JunctionTotal<'p>, rc: &RunCtx<'p>) -> Result<Products> {
        let spilled = t.table.spilled() + t.matrix.spilled();
        let mut merged = t.table.finish()?;
        let mut folded: Vec<JEntry> = Vec::new();
        while let Some(e) = merged.next_record()? {
            match folded.last_mut() {
                Some(last) if last.key() == e.key() => last.add(&e),
                _ => folded.push(e),
            }
        }
        drop(merged);
        let mut lookup = match &self.reference {
            Some(r) => Some(r.lookup(0)?),
            None => None,
        };
        let mut junctions = Vec::with_capacity(folded.len());
        for e in folded {
            let mut j = Junction {
                tid: e.tid,
                start: e.start,
                end: e.end,
                strand: 0,
                motif: 0,
                annotated: false,
                unique: e.unique,
                multi: e.multi,
                max_overhang: e.overhang,
            };
            self.annotate(&mut j, &mut lookup)?;
            junctions.push(j);
        }
        let chroms = rc.inputs[0].refs.clone();
        let matrix = match (&self.spec.matrix, &self.rows) {
            (Some(_), Some(rows)) if !junctions.is_empty() => {
                let mut features = FeatureTable::new("Splice Junction");
                for j in &junctions {
                    let chrom = &chroms[j.tid as usize].0;
                    let id = format!("{chrom}:{}-{}", j.start, j.end);
                    let strand = match j.strand {
                        1 => '+',
                        2 => '-',
                        _ => '.',
                    };
                    let len = (j.end - j.start) as u64;
                    features.push(&id, &id, chrom, j.start, j.end, strand, len);
                }
                let mut cells = t.matrix.finish()?;
                let mut sorter = Sorter::new(rc.pool, rc.buffer / 2, &rc.tmp);
                let mut chunk: Vec<Entry> = Vec::new();
                let mut col = 0usize;
                let mut held: Option<MEntry> = None;
                let fold = Fold::plain(self.reduce);
                let flush = |e: MEntry, col: &mut usize, chunk: &mut Vec<Entry>| -> Result<()> {
                    let key = (e.tid, e.start, e.end);
                    while *col < junctions.len()
                        && (
                            junctions[*col].tid,
                            junctions[*col].start,
                            junctions[*col].end,
                        ) < key
                    {
                        *col += 1;
                    }
                    chunk.push(Entry {
                        row: e.row,
                        col: *col as u32,
                        layer: 0,
                        sub: e.sub,
                        value: e.value,
                    });
                    Ok(())
                };
                while let Some(e) = cells.next_record()? {
                    match &mut held {
                        Some(h) if h.key() == e.key() => {
                            h.value = self.reduce.combine(h.value, e.value)
                        }
                        _ => {
                            if let Some(h) = held.replace(e) {
                                flush(h, &mut col, &mut chunk)?;
                                if chunk.len() >= 1 << 16 {
                                    collapse_cells(&mut chunk, fold);
                                    sorter.extend(std::mem::take(&mut chunk))?;
                                }
                            }
                        }
                    }
                }
                if let Some(h) = held {
                    flush(h, &mut col, &mut chunk)?;
                }
                collapse_cells(&mut chunk, fold);
                sorter.extend(chunk)?;
                let names = match rows.fixed_names() {
                    Some(n) => n.to_vec(),
                    None => t.rows.names,
                };
                Some(CountMatrix {
                    rows: BarcodeTable { names },
                    features: Arc::new(features),
                    layers: vec!["counts".into()],
                    named_layers: false,
                    cells: Cells {
                        merge: Some(sorter.finish()?),
                        fold,
                        head: None,
                    },
                })
            }
            _ => None,
        };
        Ok((
            JunctionTable { chroms, junctions },
            matrix,
            t.tally,
            t.spliced,
            spilled,
        ))
    }

    fn deliver(
        &self,
        (table, matrix, tally, spliced, spilled): Products,
        d: &mut Delivery<'_>,
    ) -> Result<Json> {
        let mut files = Vec::new();
        let rel = format!("{}/junctions.tab", self.spec.name);
        let path = d.path(&rel)?;
        crate::sink::junctions::write(&table, &path, d)?;
        files.push(Json::from(rel));
        let n = table.junctions.len() as u64;
        let annotated = table.junctions.iter().filter(|j| j.annotated).count() as u64;
        let canonical = table.junctions.iter().filter(|j| j.motif > 0).count() as u64;
        let mut fields = vec![
            ("product", Json::from("junctions")),
            ("junctions", Json::from(n)),
            ("annotated", Json::from(annotated)),
            ("canonical", Json::from(canonical)),
            ("spliced_records", Json::from(spliced)),
        ];
        match matrix {
            Some(m) => {
                let rows = m.rows.len() as u64;
                let written = crate::sink::matrix::write(m, &self.spec.name, false, d)?;
                fields.push(("matrix_rows", Json::from(rows)));
                fields.push(("matrix_nonzero", Json::from(written.nonzero)));
                files.extend(written.files.into_iter().map(Json::from));
            }
            None if self.spec.matrix.is_some() => fields.push(("matrix_rows", Json::from(0u64))),
            None => {}
        }
        fields.push(("spilled_bytes", Json::from(spilled)));
        fields.extend(tally.json());
        fields.push(("files", Json::Arr(files)));
        Ok(Json::obj(fields))
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn motifs_follow_star() {
        assert_eq!(motif(b"GT", b"AG"), 1);
        assert_eq!(motif(b"CT", b"AC"), 2);
        assert_eq!(motif(b"GC", b"AG"), 3);
        assert_eq!(motif(b"CT", b"GC"), 4);
        assert_eq!(motif(b"AT", b"AC"), 5);
        assert_eq!(motif(b"GT", b"AT"), 6);
        assert_eq!(motif(b"GG", b"AG"), 0);
    }
}
