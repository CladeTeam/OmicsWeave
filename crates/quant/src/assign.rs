//! Feature assignment: feature spaces (the columns of a count matrix) over
//! [`omicsweave_gtf::IntervalSet`]s, and the per-batch sweep cursors that
//! assign each read once.
//!
//! A [`FeatureSpace`] is an interval set whose rows map to columns: the
//! merged exons of a gene all map to the gene, a peak or a genome bin to
//! itself. Spaces are registered once per task ([`Spaces`]) and shared by
//! every accumulator that counts over them; within a batch,
//! [`ReadView::assign`](crate::ReadView::assign) runs a query against a
//! space once per record and caches the result.

use std::cell::RefCell;
use std::path::Path;
use std::sync::Arc;

use omicsweave_bbi::BbiReader;
use omicsweave_core::{
    Error, FileSource, IntervalEntry, IntervalRead, IntervalsOut, Region, Result, StringColumn,
};
use omicsweave_gtf::{
    Aliases, Annotation, Flank, IntervalSet, IntervalSetBuilder, Level, NONE, Selection, Strand,
    SweepCursor, match_names,
};

use crate::model::FeatureTable;
use crate::prelude::{Hits, ReadView};

/// A space registered with a task.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub struct SpaceId(pub usize);

/// What part of a read is matched against features.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum Query {
    /// The aligned blocks (`M` / `=` / `X`).
    Blocks,
    /// The aligned blocks with deletions joined.
    BlocksWithDeletions,
    /// The whole alignment span, introns included.
    Span,
    /// The pair's fragment ([`ReadView::fragment`]).
    Fragment,
    /// The Tn5-corrected fragment ([`ReadView::atac_fragment`]).
    AtacFragment,
    /// The first Tn5 insertion of the corrected fragment (its start).
    CutStart,
    /// The second Tn5 insertion (its end - 1).
    CutEnd,
    FivePrime,
    ThreePrime,
}

impl Query {
    pub(crate) fn code(self) -> u8 {
        self as u8
    }
}

/// The columns of a count matrix as intervals.
#[derive(Debug, Clone)]
pub struct FeatureSpace {
    /// One row per column.
    pub features: Arc<FeatureTable>,
    /// The intervals; `column[i]` is the column of row `i`.
    pub set: IntervalSet,
    pub column: Vec<u32>,
    /// What chromosome names of the set match: the annotation's seqids,
    /// or the set's own names.
    names: Vec<String>,
}

fn strand_char(s: Strand) -> char {
    match s {
        Strand::Forward => '+',
        Strand::Reverse => '-',
        _ => '.',
    }
}

impl FeatureSpace {
    /// The rows of an annotation-derived set whose owner is a gene, one
    /// column per gene of `genes` (in that order); `lengths` gives each
    /// gene's feature length.
    fn of_genes(
        ann: &Annotation,
        set: IntervalSet,
        genes: &[u32],
        kind: &str,
        length: impl Fn(u32) -> u64,
    ) -> FeatureSpace {
        let g = ann.genes();
        let mut col_of = vec![NONE; g.len()];
        for (c, &gene) in genes.iter().enumerate() {
            col_of[gene as usize] = c as u32;
        }
        let keep: Vec<u32> = (0..set.len() as u32)
            .filter(|&r| {
                let o = set.owner[r as usize];
                o != NONE && col_of[o as usize] != NONE
            })
            .collect();
        let set = set.select(&keep);
        let column = set.owner.iter().map(|&o| col_of[o as usize]).collect();
        let mut t = FeatureTable::new(kind);
        for &gene in genes {
            let i = gene as usize;
            let name = match ann.str(g.name[i]) {
                "" => ann.str(g.id[i]),
                n => n,
            };
            t.push(
                ann.str(g.id[i]),
                name,
                ann.seqid_name(g.seqid[i]),
                g.start[i],
                g.end[i],
                strand_char(g.strand[i]),
                length(gene),
            );
        }
        FeatureSpace {
            features: Arc::new(t),
            set,
            column,
            names: ann.seqid_names().map(str::to_owned).collect(),
        }
    }

    /// Genes by their merged exons (featureCounts' meta-features, STARsolo
    /// `Gene`); a gene's length is its exonic length.
    pub fn gene_exons(ann: &Annotation, selection: Option<&Selection>) -> FeatureSpace {
        let genes = gene_rows(ann, selection);
        let mut d = ann.derive();
        if let Some(s) = selection {
            d = d.selection(s);
        }
        let set = d.exonic_regions();
        let mut len = vec![0u64; ann.genes().len()];
        for i in 0..set.len() {
            len[set.owner[i] as usize] += (set.end[i] - set.start[i]) as u64;
        }
        FeatureSpace::of_genes(ann, set, &genes, "Gene Expression", |g| len[g as usize])
    }

    /// Genes by their whole span, introns included (STARsolo `GeneFull`).
    pub fn gene_spans(ann: &Annotation, selection: Option<&Selection>) -> FeatureSpace {
        let genes = gene_rows(ann, selection);
        let mut d = ann.derive();
        if let Some(s) = selection {
            d = d.selection(s);
        }
        let set = d.gene_bodies(Level::Gene, Flank::new(0, 0));
        let g = ann.genes();
        FeatureSpace::of_genes(ann, set, &genes, "Gene Expression", |i| {
            (g.end[i as usize] - g.start[i as usize]) as u64
        })
    }

    /// One window per gene: around its TSS (`promoter`), or its body
    /// extended upstream and downstream (`body`, gene activity), clipped to
    /// `lengths` (by annotation seqid).
    pub fn gene_windows(
        ann: &Annotation,
        selection: Option<&Selection>,
        lengths: &[Option<u32>],
        promoter: bool,
        flank: Flank,
    ) -> FeatureSpace {
        let genes = gene_rows(ann, selection);
        let mut d = ann.derive().clip(lengths);
        if let Some(s) = selection {
            d = d.selection(s);
        }
        let set = if promoter {
            d.promoters(Level::Gene, flank)
        } else {
            d.gene_bodies(Level::Gene, flank)
        };
        let mut len = vec![0u64; ann.genes().len()];
        for i in 0..set.len() {
            len[set.owner[i] as usize] = (set.end[i] - set.start[i]) as u64;
        }
        let kind = if promoter {
            "Promoter"
        } else {
            "Gene Activity"
        };
        let mut s = FeatureSpace::of_genes(ann, set, &genes, kind, |g| len[g as usize]);
        // A window's coordinates describe the feature better than the gene's.
        let t = Arc::make_mut(&mut s.features);
        for r in 0..s.set.len() {
            let c = s.column[r] as usize;
            t.start[c] = s.set.start[r];
            t.end[c] = s.set.end[r];
        }
        s
    }

    /// Fixed tiles of `size` bases over every chromosome (the last one
    /// shorter), named `chrom:start-end`.
    pub fn bins(chroms: &[(String, u32)], size: u32) -> Result<FeatureSpace> {
        if size == 0 {
            return Err(Error::invalid("bins", "the bin size must be at least 1"));
        }
        let mut b = IntervalSetBuilder::new().chroms(chroms.iter().map(|c| c.0.as_str()))?;
        let mut t = FeatureTable::new("Bins");
        for (name, len) in chroms {
            let mut s = 0u32;
            while s < *len {
                let e = s.saturating_add(size).min(*len);
                b.push(name, s, e, Strand::Unstranded, NONE, None)?;
                let id = format!("{name}:{s}-{e}");
                t.push(&id, &id, name, s, e, '.', (e - s) as u64);
                s = e;
            }
        }
        let set = b.build();
        Ok(FeatureSpace {
            column: (0..set.len() as u32).collect(),
            features: Arc::new(t),
            names: chroms.iter().map(|c| c.0.clone()).collect(),
            set,
        })
    }

    /// Every feature of a bigBed file, one column each in file order (by
    /// chromosome id, then start), named `chrom:start-end`; the strand is
    /// the sixth BED field when there is one.
    pub fn big_bed(path: &Path, kind: &str) -> Result<FeatureSpace> {
        let mut r = BbiReader::open(FileSource::open(path)?)?;
        let chroms: Vec<(String, u32)> = r
            .chroms()
            .iter()
            .map(|(_, c)| (c.name.clone(), c.len))
            .collect();
        let mut b = IntervalSetBuilder::new().chroms(chroms.iter().map(|c| c.0.as_str()))?;
        let mut t = FeatureTable::new(kind);
        for (name, len) in &chroms {
            let size = r.query(
                Region {
                    chrom: name,
                    start: 0,
                    end: *len,
                },
                false,
            )?;
            let mut entries = vec![IntervalEntry::default(); size.entries];
            let mut offsets = vec![0i64; size.entries + 1];
            let mut data = vec![0u8; size.annotation_bytes];
            r.write(IntervalsOut {
                entries: &mut entries,
                annotations: StringColumn {
                    offsets: &mut offsets,
                    data: &mut data,
                },
            })?;
            for e in &entries {
                let a = e.annotation as usize;
                let rest = &data[offsets[a] as usize..offsets[a + 1] as usize];
                let strand = rest
                    .split(|&c| c == b'\t')
                    .nth(2)
                    .and_then(|f| match f {
                        b"+" => Some(Strand::Forward),
                        b"-" => Some(Strand::Reverse),
                        _ => None,
                    })
                    .unwrap_or(Strand::Unstranded);
                b.push(name, e.start, e.end, strand, NONE, None)?;
                let id = format!("{name}:{}-{}", e.start, e.end);
                t.push(
                    &id,
                    &id,
                    name,
                    e.start,
                    e.end,
                    strand_char(strand),
                    (e.end - e.start) as u64,
                );
            }
        }
        let set = b.build();
        Ok(FeatureSpace {
            column: (0..set.len() as u32).collect(),
            features: Arc::new(t),
            names: chroms.into_iter().map(|c| c.0).collect(),
            set,
        })
    }

    /// Windows around the TSS of every gene: `flank` bases upstream and
    /// downstream, each its own column (for per-barcode TSS counts).
    pub fn tss_windows(ann: &Annotation, lengths: &[Option<u32>], flank: Flank) -> FeatureSpace {
        let set = ann.derive().clip(lengths).promoters(Level::Gene, flank);
        let mut t = FeatureTable::new("TSS");
        for r in 0..set.len() {
            let g = set.owner[r] as usize;
            t.push(
                ann.str(ann.genes().id[g]),
                ann.str(ann.genes().id[g]),
                ann.seqid_name(set.seqid[r]),
                set.start[r],
                set.end[r],
                strand_char(set.strand[r]),
                (set.end[r] - set.start[r]) as u64,
            );
        }
        FeatureSpace {
            column: (0..set.len() as u32).collect(),
            features: Arc::new(t),
            names: ann.seqid_names().map(str::to_owned).collect(),
            set,
        }
    }

    /// Every row of `set` as its own column (regions built by the caller:
    /// introns, exon bins, peaks), with id `chrom:start-end` and the row's
    /// name, or the id when it has none. Invalid public columns are rejected;
    /// the overlap index is built from the owned coordinates on first query.
    pub fn of_intervals(mut set: IntervalSet, kind: &str) -> Result<FeatureSpace> {
        let n = set.len();
        if n > u32::MAX as usize || set.seqid_count() > u32::MAX as usize {
            return Err(Error::invalid(
                "feature space",
                "too many intervals or chromosomes",
            ));
        }
        if [
            set.seqid.len(),
            set.end.len(),
            set.strand.len(),
            set.owner.len(),
            set.name.len(),
        ]
        .iter()
        .any(|&len| len != n)
        {
            return Err(Error::invalid(
                "feature space",
                "interval columns differ in length",
            ));
        }
        if set.seqid.iter().any(|&s| s as usize >= set.seqid_count()) {
            return Err(Error::invalid(
                "feature space",
                "seqid outside the chromosome list",
            ));
        }
        if set
            .name
            .iter()
            .any(|&name| name != NONE && name as usize >= set.pool().len())
        {
            return Err(Error::invalid(
                "feature space",
                "name outside the string pool",
            ));
        }
        if let Some((&start, &end)) = set.start.iter().zip(&set.end).find(|(s, e)| s > e) {
            return Err(Error::InvalidRegion { start, end });
        }
        // Public coordinate columns may have been edited after an index was built.
        set.clear_index();
        let mut t = FeatureTable::new(kind);
        for r in 0..set.len() {
            let chrom = set.seqid_name(set.seqid[r]);
            let (start, end) = (set.start[r], set.end[r]);
            let id = format!("{chrom}:{start}-{end}");
            let name = set.name(r).filter(|n| !n.is_empty()).unwrap_or(id.as_str());
            t.push(
                &id,
                name,
                chrom,
                start,
                end,
                strand_char(set.strand[r]),
                (end - start) as u64,
            );
        }
        Ok(FeatureSpace {
            column: (0..set.len() as u32).collect(),
            features: Arc::new(t),
            names: (0..set.seqid_count() as u32)
                .map(|i| set.seqid_name(i).to_owned())
                .collect(),
            set,
        })
    }

    /// A space over the rows of any interval set: `column[i]` is row `i`'s
    /// column of `features`; `chroms` names the set's seqids.
    pub fn from_set(
        set: IntervalSet,
        column: Vec<u32>,
        features: FeatureTable,
        chroms: Vec<String>,
    ) -> Result<FeatureSpace> {
        if column.len() != set.len() || column.iter().any(|&c| c as usize >= features.len()) {
            return Err(Error::invalid(
                "feature space",
                "every row needs a column of the feature table",
            ));
        }
        Ok(FeatureSpace {
            features: Arc::new(features),
            set,
            column,
            names: chroms,
        })
    }

    /// For each reference of a BAM, the set's seqid of the same chromosome
    /// (names matched as [`match_names`] matches them).
    pub fn match_references(&self, refs: &[&str], aliases: Option<&Aliases>) -> Vec<Option<u32>> {
        let names: Vec<&str> = self.names.iter().map(String::as_str).collect();
        match_names(&names, refs, aliases)
    }

    /// Bytes held: features, intervals and their index.
    pub fn heap_bytes(&self) -> u64 {
        (self.set.heap_bytes() + 4 * self.column.capacity()) as u64 + self.features.heap_bytes()
    }
}

/// The gene rows a selection keeps (all genes without one), ascending.
fn gene_rows(ann: &Annotation, selection: Option<&Selection>) -> Vec<u32> {
    match selection {
        Some(s) => s.genes.clone(),
        None => (0..ann.genes().len() as u32).collect(),
    }
}

/// The feature spaces of a task, with each input's reference → seqid
/// mapping.
#[derive(Default)]
pub struct Spaces {
    pub(crate) spaces: Vec<(String, Arc<FeatureSpace>)>,
    /// `maps[space][input][tid]`.
    pub(crate) maps: Vec<Vec<Vec<Option<u32>>>>,
}

impl Spaces {
    /// Registers `space` under `key` (a description of how it was built),
    /// or returns the space already registered under it.
    pub(crate) fn register(
        &mut self,
        key: String,
        build: impl FnOnce() -> Result<FeatureSpace>,
    ) -> Result<SpaceId> {
        if let Some(i) = self.spaces.iter().position(|(k, _)| *k == key) {
            return Ok(SpaceId(i));
        }
        self.spaces.push((key, Arc::new(build()?)));
        Ok(SpaceId(self.spaces.len() - 1))
    }

    pub fn get(&self, id: SpaceId) -> &Arc<FeatureSpace> {
        &self.spaces[id.0].1
    }

    /// Matches every space against the references of every input.
    pub(crate) fn match_inputs(&mut self, refs: &[Vec<String>], aliases: Option<&Aliases>) {
        self.maps = self
            .spaces
            .iter()
            .map(|(_, s)| {
                refs.iter()
                    .map(|r| {
                        let names: Vec<&str> = r.iter().map(String::as_str).collect();
                        s.match_references(&names, aliases)
                    })
                    .collect()
            })
            .collect();
    }

    pub fn heap_bytes(&self) -> u64 {
        self.spaces.iter().map(|(_, s)| s.heap_bytes()).sum()
    }
}

/// One space's sweep over a batch.
struct Sweep<'a> {
    space: &'a FeatureSpace,
    map: &'a [Option<u32>],
    tid: Option<u32>,
    cursor: Option<SweepCursor<'a>>,
}

/// The sweep cursors of one batch, one per space.
pub struct BatchSpaces<'a> {
    sweeps: Vec<RefCell<Sweep<'a>>>,
}

impl<'a> BatchSpaces<'a> {
    pub(crate) fn new(spaces: &'a Spaces, input: usize) -> BatchSpaces<'a> {
        BatchSpaces {
            sweeps: spaces
                .spaces
                .iter()
                .enumerate()
                .map(|(i, (_, s))| {
                    RefCell::new(Sweep {
                        space: s,
                        map: &spaces.maps[i][input],
                        tid: None,
                        cursor: None,
                    })
                })
                .collect(),
        }
    }

    pub(crate) fn assign(
        &self,
        id: SpaceId,
        r: &ReadView<'_>,
        query: Query,
        strand: Option<bool>,
    ) -> Hits {
        self.sweep(id, r, strand, |visit| match query {
            Query::Blocks => r.blocks(false).iter().for_each(|&(s, e)| visit(s, e)),
            Query::BlocksWithDeletions => r.blocks(true).iter().for_each(|&(s, e)| visit(s, e)),
            Query::Span => visit(r.pos(), r.end().max(r.pos() + 1)),
            Query::Fragment => {
                if let Some(f) = r.fragment(0) {
                    visit(f.start, f.end)
                }
            }
            Query::AtacFragment => {
                if let Some(f) = r.atac_fragment(0) {
                    visit(f.start, f.end)
                }
            }
            Query::CutStart => {
                if let Some(f) = r.atac_fragment(0) {
                    visit(f.start, f.start + 1)
                }
            }
            Query::CutEnd => {
                if let Some(f) = r.atac_fragment(0) {
                    visit(f.end - 1, f.end)
                }
            }
            Query::FivePrime => {
                let p = r.five_prime();
                visit(p, p + 1)
            }
            Query::ThreePrime => {
                let p = r.three_prime();
                visit(p, p + 1)
            }
        })
    }

    /// The columns `[start, end)` overlaps, for an interval at or after the
    /// read's position.
    pub(crate) fn assign_at(&self, id: SpaceId, r: &ReadView<'_>, start: u32, end: u32) -> Hits {
        self.sweep(id, r, None, |visit| visit(start, end))
    }

    /// Moves the space's cursor to `r` and collects the columns of every
    /// interval `intervals` visits, ascending and distinct.
    fn sweep(
        &self,
        id: SpaceId,
        r: &ReadView<'_>,
        strand: Option<bool>,
        intervals: impl FnOnce(&mut dyn FnMut(u32, u32)),
    ) -> Hits {
        let mut hits = Hits::new();
        let mut sweep = self.sweeps[id.0].borrow_mut();
        let sw = &mut *sweep;
        let tid = r.tid();
        let Some(seqid) = sw.map.get(tid as usize).copied().flatten() else {
            return hits;
        };
        if sw.tid != Some(tid) {
            match &mut sw.cursor {
                Some(c) => c.seek(seqid),
                None => sw.cursor = Some(sw.space.set.cursor(seqid)),
            }
            sw.tid = Some(tid);
        }
        let cursor = sw.cursor.as_mut().expect("seeked above");
        cursor.advance(r.pos());
        let set = &sw.space.set;
        let column = &sw.space.column;
        let mut visit = |s: u32, e: u32| {
            cursor.for_each_overlap(s, e, |row| {
                let keep = match strand {
                    None => true,
                    Some(rev) => match set.strand[row as usize] {
                        Strand::Forward => !rev,
                        Strand::Reverse => rev,
                        _ => true,
                    },
                };
                if keep {
                    hits.push(column[row as usize]);
                }
            })
        };
        intervals(&mut visit);
        hits.sort_unstable();
        hits.dedup();
        hits
    }
}
