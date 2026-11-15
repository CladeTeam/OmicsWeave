//! [`Annotation`]: the tables, lookups and lazily built indexes.

use std::sync::{Arc, OnceLock};

use omicsweave_core::Region;

use crate::alias::{Aliases, match_names};
use crate::derived::{ExonBins, Junctions};
use crate::index::{IntervalIndex, SweepCursor};
use crate::model::{
    AttributePolicy, Attributes, Dialect, FeatureKind, Features, Genes, Others, Row, RowRef,
    Transcripts,
};
use crate::report::Report;
use crate::set::{IntervalSet, Level};
use crate::strings::{IdIndex, IdMap, NONE, StringPool};

/// A gene annotation (GTF or GFF3) as columnar tables.
///
/// Genes → transcripts → features form a hierarchy of contiguous row
/// ranges ([`Genes::transcripts`], [`Transcripts::features`]); every string
/// is interned once in [`Annotation::pool`]. Indexes, lookups and derived
/// views are built on first use and cached; the annotation is immutable
/// and `Sync`, so threads share it by reference.
#[derive(Debug)]
pub struct Annotation {
    pub(crate) pool: Arc<StringPool>,
    pub(crate) dialect: Dialect,
    pub(crate) policy: AttributePolicy,
    pub(crate) seqid_names: Arc<[u32]>,
    pub(crate) seqid_of: IdMap<u32, u32>,
    pub(crate) genes: Genes,
    pub(crate) transcripts: Transcripts,
    pub(crate) features: Features,
    pub(crate) others: Others,
    pub(crate) attributes: Attributes,
    pub(crate) rows: Vec<RowRef>,
    /// (records before the line, text) of every `#` line.
    pub(crate) comments: Vec<(u32, u32)>,
    pub(crate) report: Report,
    pub(crate) cache: Cache,
}

#[derive(Debug, Default)]
pub(crate) struct Cache {
    pub gene_ids: OnceLock<IdIndex>,
    pub gene_names: OnceLock<IdIndex>,
    pub transcript_ids: OnceLock<IdIndex>,
    pub transcript_names: OnceLock<IdIndex>,
    pub exons: OnceLock<Vec<u32>>,
    pub gene_index: OnceLock<IntervalIndex>,
    pub transcript_index: OnceLock<IntervalIndex>,
    pub exon_index: OnceLock<IntervalIndex>,
    pub introns: OnceLock<IntervalSet>,
    pub exonic: OnceLock<IntervalSet>,
    pub junctions: OnceLock<Junctions>,
    pub bins: [OnceLock<ExonBins>; 2],
    pub transcript_lengths: OnceLock<Vec<u32>>,
    pub cds_lengths: OnceLock<Vec<u32>>,
}

impl Annotation {
    /// Every string of the annotation.
    pub fn pool(&self) -> &StringPool {
        &self.pool
    }

    /// String `id` of the pool (`""` for [`NONE`]).
    #[inline]
    pub fn str(&self, id: u32) -> &str {
        self.pool.str(id)
    }

    pub fn dialect(&self) -> Dialect {
        self.dialect
    }

    pub fn attribute_policy(&self) -> &AttributePolicy {
        &self.policy
    }

    pub fn genes(&self) -> &Genes {
        &self.genes
    }

    pub fn transcripts(&self) -> &Transcripts {
        &self.transcripts
    }

    pub fn features(&self) -> &Features {
        &self.features
    }

    pub fn others(&self) -> &Others {
        &self.others
    }

    /// The flat attribute table, by file record.
    pub fn attribute_table(&self) -> &Attributes {
        &self.attributes
    }

    pub fn report(&self) -> &Report {
        &self.report
    }

    /// Chromosomes, numbered in order of first appearance.
    pub fn seqid_count(&self) -> usize {
        self.seqid_names.len()
    }

    pub fn seqid_name(&self, seqid: u32) -> &str {
        self.pool.str(self.seqid_names[seqid as usize])
    }

    pub fn seqid_names(&self) -> impl ExactSizeIterator<Item = &str> + '_ {
        self.seqid_names.iter().map(|&n| self.pool.str(n))
    }

    /// The seqid called `name`.
    pub fn seqid(&self, name: &str) -> Option<u32> {
        let id = self.pool.get(name.as_bytes())?;
        self.seqid_of.get(&id).copied()
    }

    /// For each of `targets` (a BAM header's reference names, say), the
    /// annotation seqid it refers to: same name, `aliases`, or the
    /// built-in UCSC ↔ Ensembl spellings (`chr1` ↔ `1`, `chrM` ↔ `MT`).
    pub fn match_seqids(&self, targets: &[&str], aliases: Option<&Aliases>) -> Vec<Option<u32>> {
        let names: Vec<&str> = self.seqid_names().collect();
        match_names(&names, targets, aliases)
    }

    /// Lengths by annotation seqid from `(name, length)` pairs of another
    /// file, matched as by [`Annotation::match_seqids`]; for clipping
    /// derived windows.
    pub fn seqid_lengths(
        &self,
        chroms: &[(&str, u32)],
        aliases: Option<&Aliases>,
    ) -> Vec<Option<u32>> {
        let names: Vec<&str> = chroms.iter().map(|c| c.0).collect();
        let mut out = vec![None; self.seqid_count()];
        for (i, s) in self.match_seqids(&names, aliases).into_iter().enumerate() {
            if let Some(s) = s {
                out[s as usize].get_or_insert(chroms[i].1);
            }
        }
        out
    }

    /// File records in file order, as the table rows they became.
    pub fn rows(&self) -> impl ExactSizeIterator<Item = Row> + '_ {
        self.rows.iter().map(|r| r.get())
    }

    /// The file record behind a row ([`NONE`] when synthesized).
    pub fn record_of(&self, row: Row) -> u32 {
        match row {
            Row::Gene(i) => self.genes.record[i as usize],
            Row::Transcript(i) => self.transcripts.record[i as usize],
            Row::Feature(i) => self.features.record[i as usize],
            Row::Other(i) => self.others.record[i as usize],
        }
    }

    /// The attributes of a row as (key, value), in file order, as kept by
    /// the attribute policy. GFF3 values are as written (not
    /// percent-decoded; multiple values stay comma-joined).
    pub fn attributes(&self, row: Row) -> impl Iterator<Item = (&str, &str)> + '_ {
        self.attributes
            .of(self.record_of(row))
            .map(|(k, v)| (self.pool.str(k), self.pool.str(v)))
    }

    /// The first value of attribute `key` of a row.
    pub fn attribute(&self, row: Row, key: &str) -> Option<&str> {
        let k = self.pool.get(key.as_bytes())?;
        self.attributes
            .of(self.record_of(row))
            .find(|&(kk, _)| kk == k)
            .map(|(_, v)| self.pool.str(v))
    }

    /// The seqid of feature row `f`.
    #[inline]
    pub fn feature_seqid(&self, f: u32) -> u32 {
        self.transcripts.seqid[self.features.transcript[f as usize] as usize]
    }

    /// The gene of feature row `f`.
    #[inline]
    pub fn feature_gene(&self, f: u32) -> u32 {
        self.transcripts.gene[self.features.transcript[f as usize] as usize]
    }

    fn lookup<'a>(&'a self, cell: &'a OnceLock<IdIndex>, keys: &[u32], s: &str) -> &'a [u32] {
        let Some(id) = self.pool.get(s.as_bytes()) else {
            return &[];
        };
        cell.get_or_init(|| IdIndex::build(keys)).rows(id)
    }

    /// Genes with id `id`, ascending (an id may repeat on several seqids).
    pub fn genes_by_id(&self, id: &str) -> &[u32] {
        self.lookup(&self.cache.gene_ids, &self.genes.id, id)
    }

    /// The first gene with id `id`.
    pub fn gene(&self, id: &str) -> Option<u32> {
        self.genes_by_id(id).first().copied()
    }

    pub fn genes_by_name(&self, name: &str) -> &[u32] {
        self.lookup(&self.cache.gene_names, &self.genes.name, name)
    }

    pub fn transcripts_by_id(&self, id: &str) -> &[u32] {
        self.lookup(&self.cache.transcript_ids, &self.transcripts.id, id)
    }

    pub fn transcript(&self, id: &str) -> Option<u32> {
        self.transcripts_by_id(id).first().copied()
    }

    pub fn transcripts_by_name(&self, name: &str) -> &[u32] {
        self.lookup(&self.cache.transcript_names, &self.transcripts.name, name)
    }

    /// Feature rows that are exons, ascending.
    pub fn exons(&self) -> &[u32] {
        self.cache.exons.get_or_init(|| {
            (0..self.features.len() as u32)
                .filter(|&f| self.features.kind[f as usize] == FeatureKind::Exon)
                .collect()
        })
    }

    /// The overlap index of genes, transcripts or exons; ids are gene
    /// rows, transcript rows, or feature rows.
    pub fn interval_index(&self, level: Level) -> &IntervalIndex {
        let n = self.seqid_count();
        match level {
            Level::Gene => self.cache.gene_index.get_or_init(|| {
                IntervalIndex::build(n, &self.genes.seqid, &self.genes.start, &self.genes.end)
            }),
            Level::Transcript => self.cache.transcript_index.get_or_init(|| {
                IntervalIndex::build(
                    n,
                    &self.transcripts.seqid,
                    &self.transcripts.start,
                    &self.transcripts.end,
                )
            }),
            Level::Exon => self.cache.exon_index.get_or_init(|| {
                IntervalIndex::build_with(
                    n,
                    self.exons().iter().map(|&f| {
                        let i = f as usize;
                        (
                            self.feature_seqid(f),
                            self.features.start[i],
                            self.features.end[i],
                            f,
                        )
                    }),
                )
            }),
        }
    }

    /// Appends the ids at `level` overlapping `region`; an unknown
    /// chromosome has none.
    pub fn overlaps(&self, level: Level, region: Region<'_>, out: &mut Vec<u32>) {
        if let Some(s) = self.seqid(region.chrom) {
            self.interval_index(level)
                .overlaps(s, region.start, region.end, out);
        }
    }

    /// A sweep cursor over `seqid` at `level`.
    pub fn cursor(&self, level: Level, seqid: u32) -> SweepCursor<'_> {
        self.interval_index(level).cursor(seqid)
    }

    /// An empty interval set sharing this annotation's seqids and strings.
    pub(crate) fn new_set(&self, owner: Option<Level>) -> IntervalSet {
        IntervalSet::empty(Arc::clone(&self.pool), Arc::clone(&self.seqid_names), owner)
    }

    /// Heap bytes of the tables, strings and whatever caches are built.
    pub fn heap_bytes(&self) -> usize {
        let g = &self.genes;
        let t = &self.transcripts;
        let f = &self.features;
        let o = &self.others;
        let a = &self.attributes;
        let v4 = |v: &Vec<u32>| 4 * v.capacity();
        let genes = [
            &g.seqid,
            &g.start,
            &g.end,
            &g.id,
            &g.name,
            &g.biotype,
            &g.source,
            &g.feature_type,
            &g.score,
            &g.record,
            &g.transcript_offsets,
        ]
        .iter()
        .map(|v| v4(v))
        .sum::<usize>()
            + g.strand.capacity()
            + g.phase.capacity();
        let transcripts = [
            &t.seqid,
            &t.start,
            &t.end,
            &t.gene,
            &t.id,
            &t.name,
            &t.biotype,
            &t.source,
            &t.feature_type,
            &t.score,
            &t.record,
            &t.tag_values,
            &t.tag_offsets,
            &t.feature_offsets,
        ]
        .iter()
        .map(|v| v4(v))
        .sum::<usize>()
            + t.strand.capacity()
            + t.phase.capacity()
            + 4 * t.tag_bits.capacity();
        let features = [
            &f.transcript,
            &f.start,
            &f.end,
            &f.feature_type,
            &f.source,
            &f.score,
            &f.record,
        ]
        .iter()
        .map(|v| v4(v))
        .sum::<usize>()
            + f.kind.capacity()
            + f.strand.capacity()
            + f.phase.capacity();
        let others = [
            &o.seqid,
            &o.start,
            &o.end,
            &o.feature_type,
            &o.source,
            &o.score,
            &o.id,
            &o.record,
        ]
        .iter()
        .map(|v| v4(v))
        .sum::<usize>()
            + o.strand.capacity()
            + o.phase.capacity();
        let attrs = a.heap_bytes();
        let c = &self.cache;
        let ids = [
            &c.gene_ids,
            &c.gene_names,
            &c.transcript_ids,
            &c.transcript_names,
        ]
        .iter()
        .filter_map(|x| x.get())
        .map(IdIndex::heap_bytes)
        .sum::<usize>();
        let indexes = [&c.gene_index, &c.transcript_index, &c.exon_index]
            .iter()
            .filter_map(|x| x.get())
            .map(IntervalIndex::heap_bytes)
            .sum::<usize>();
        let sets = [&c.introns, &c.exonic]
            .iter()
            .filter_map(|x| x.get())
            .map(IntervalSet::heap_bytes)
            .sum::<usize>();
        genes
            + transcripts
            + features
            + others
            + attrs
            + self.pool.heap_bytes()
            + 4 * self.rows.capacity()
            + 8 * self.comments.capacity()
            + ids
            + indexes
            + sets
            + c.exons.get().map_or(0, |v| 4 * v.capacity())
            + c.junctions.get().map_or(0, Junctions::heap_bytes)
            + c.bins
                .iter()
                .filter_map(OnceLock::get)
                .map(ExonBins::heap_bytes)
                .sum::<usize>()
            + [&c.transcript_lengths, &c.cds_lengths]
                .iter()
                .filter_map(|c| c.get())
                .map(|v| 4 * v.capacity())
                .sum::<usize>()
    }

    /// Display text of a pool id that may be [`NONE`].
    pub(crate) fn text(&self, id: u32) -> &str {
        if id == NONE { "" } else { self.pool.str(id) }
    }
}
