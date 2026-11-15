//! Views computed from the tables: exon blocks, introns, exonic regions,
//! junctions, disjoint exon bins, TSS / TES / promoter / gene-body windows,
//! transcript and CDS lengths.
//!
//! Every interval view is an [`IntervalSet`] whose owner column names the
//! gene or transcript each row comes from, so it is indexed and swept
//! exactly like genes and exons. The views of the whole annotation that
//! take no parameters are cached on the [`Annotation`]; [`Derive`] builds
//! them for a [`Selection`], clipped to chromosome lengths, or with
//! windows.

use crate::annotation::Annotation;
use crate::model::{FeatureKind, Strand};
use crate::select::Selection;
use crate::set::{Csr, IntervalSet, Level};
use crate::strings::{IdMap, NONE, id_map};

/// Bases upstream and downstream of a site, in the direction of
/// transcription (for `-` strand rows, upstream is the higher coordinate;
/// rows without a known strand are read as `+`).
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub struct Flank {
    pub upstream: u32,
    pub downstream: u32,
}

impl Flank {
    pub fn new(upstream: u32, downstream: u32) -> Flank {
        Flank {
            upstream,
            downstream,
        }
    }
}

/// The annotated splice junctions: introns of all transcripts, each
/// distinct (seqid, start, end, strand) once.
///
/// `set` rows are sorted by (seqid, start, end, strand); `set.owner` is the
/// first gene, `genes` and `transcripts` list every gene and transcript
/// with the junction. Intron coordinates are 0-based, half-open: the
/// donor site is the first intron base on the transcript strand, the
/// acceptor the last.
#[derive(Debug, Clone)]
pub struct Junctions {
    pub set: IntervalSet,
    pub genes: Csr,
    pub transcripts: Csr,
    lookup: IdMap<(u64, u64), u32>,
}

impl Junctions {
    pub(crate) fn heap_bytes(&self) -> usize {
        self.set.heap_bytes()
            + self.genes.heap_bytes()
            + self.transcripts.heap_bytes()
            + 24 * self.lookup.capacity()
    }

    pub fn len(&self) -> usize {
        self.set.len()
    }

    pub fn is_empty(&self) -> bool {
        self.set.is_empty()
    }

    /// The junction with exactly these coordinates and strand.
    pub fn find(&self, seqid: u32, start: u32, end: u32, strand: Strand) -> Option<u32> {
        self.lookup.get(&jkey(seqid, start, end, strand)).copied()
    }

    /// The junction with these coordinates on any strand (`+`, then `-`,
    /// then unstranded), for reads of an unstranded library.
    pub fn find_any(&self, seqid: u32, start: u32, end: u32) -> Option<u32> {
        [
            Strand::Forward,
            Strand::Reverse,
            Strand::Unstranded,
            Strand::Unknown,
        ]
        .into_iter()
        .find_map(|s| self.find(seqid, start, end, s))
    }

    /// The donor position of junction `i` (first intron base, 5' side).
    pub fn donor(&self, i: u32) -> u32 {
        let i = i as usize;
        if self.set.strand[i].is_reverse() {
            self.set.end[i] - 1
        } else {
            self.set.start[i]
        }
    }

    /// The acceptor position of junction `i` (last intron base, 3' side).
    pub fn acceptor(&self, i: u32) -> u32 {
        let i = i as usize;
        if self.set.strand[i].is_reverse() {
            self.set.start[i]
        } else {
            self.set.end[i] - 1
        }
    }
}

#[inline]
fn jkey(seqid: u32, start: u32, end: u32, strand: Strand) -> (u64, u64) {
    (
        ((seqid as u64) << 8) | strand as u64,
        ((start as u64) << 32) | end as u64,
    )
}

/// Disjoint exon bins, as DEXSeq's exonic parts: within a gene (or an
/// aggregate of genes whose exons overlap on the same strand), the exonic
/// bases cut wherever the set of transcripts covering them changes.
///
/// `set.owner` is the gene (the first of an aggregate) and `set.name` its
/// id; `part` numbers the bins of a gene from 1 in genomic order;
/// `transcripts` lists the transcripts covering each bin, and `genes` the
/// genes of each bin's aggregate (one gene unless aggregated).
#[derive(Debug, Clone)]
pub struct ExonBins {
    pub set: IntervalSet,
    pub part: Vec<u32>,
    pub transcripts: Csr,
    pub genes: Csr,
}

impl ExonBins {
    pub(crate) fn heap_bytes(&self) -> usize {
        self.set.heap_bytes()
            + 4 * self.part.capacity()
            + self.transcripts.heap_bytes()
            + self.genes.heap_bytes()
    }
}

/// A part of a transcript, read as the spans of a spliced key.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum TranscriptPart {
    /// The exon blocks ([`Annotation::exon_blocks`]): the mature
    /// transcript.
    Exons,
    /// The coding blocks, stop codon included: CDS rows merged with the
    /// codon rows, since GENCODE and Ensembl GTF list `stop_codon` outside
    /// `CDS` while GFF3 usually includes it. A complete ORF is then a
    /// multiple of three long.
    Cds,
    /// The exonic bases 5' of the coding blocks.
    Utr5,
    /// The exonic bases 3' of the coding blocks (after the stop codon).
    Utr3,
}

impl TranscriptPart {
    pub fn parse(s: &str) -> Option<TranscriptPart> {
        Some(match s {
            "exons" => TranscriptPart::Exons,
            "cds" => TranscriptPart::Cds,
            "utr5" => TranscriptPart::Utr5,
            "utr3" => TranscriptPart::Utr3,
            _ => return None,
        })
    }
}

impl Annotation {
    /// A builder of derived views over all of the annotation.
    pub fn derive(&self) -> Derive<'_> {
        Derive {
            ann: self,
            selection: None,
            lengths: None,
        }
    }

    /// The merged exon blocks of transcript `t`, ascending, into `out`
    /// (cleared first): its exon rows, or, for a transcript without exon
    /// rows, its CDS, UTR and codon rows; overlapping or abutting pieces
    /// are merged.
    pub fn exon_blocks(&self, t: u32, out: &mut Vec<(u32, u32)>) {
        out.clear();
        let f = &self.features;
        let rows = self.transcripts.features(t);
        let has_exon = rows
            .clone()
            .any(|i| f.kind[i as usize] == FeatureKind::Exon);
        for i in rows {
            let i = i as usize;
            let k = f.kind[i];
            if (has_exon && k != FeatureKind::Exon) || (!has_exon && !k.is_exonic()) {
                continue;
            }
            push_merged(out, f.start[i], f.end[i]);
        }
    }

    /// The blocks of `part` of transcript `t`, ascending and merged, into
    /// `out` (cleared first); see [`TranscriptPart`]. A non-coding
    /// transcript has no coding part and no UTRs.
    pub fn part_blocks(&self, t: u32, part: TranscriptPart, out: &mut Vec<(u32, u32)>) {
        if part == TranscriptPart::Exons {
            return self.exon_blocks(t, out);
        }
        out.clear();
        let f = &self.features;
        for i in self.transcripts.features(t) {
            let i = i as usize;
            if f.kind[i].is_coding() {
                push_merged(out, f.start[i], f.end[i]);
            }
        }
        if part == TranscriptPart::Cds {
            return;
        }
        let (Some(first), Some(last)) = (out.first(), out.last()) else {
            return;
        };
        let (lo, hi) = (first.0, last.1);
        let low_side =
            (part == TranscriptPart::Utr5) != self.transcripts.strand[t as usize].is_reverse();
        let mut exons = Vec::new();
        self.exon_blocks(t, &mut exons);
        out.clear();
        for (s, e) in exons {
            let (s, e) = if low_side {
                (s, e.min(lo))
            } else {
                (s.max(hi), e)
            };
            if s < e {
                out.push((s, e));
            }
        }
    }

    /// Introns of every transcript (owner: transcript), cached.
    pub fn introns(&self) -> &IntervalSet {
        self.cache.introns.get_or_init(|| self.derive().introns())
    }

    /// Merged exonic regions of every gene (owner: gene), cached.
    pub fn exonic_regions(&self) -> &IntervalSet {
        self.cache
            .exonic
            .get_or_init(|| self.derive().exonic_regions())
    }

    /// The annotated junctions, cached.
    pub fn junctions(&self) -> &Junctions {
        self.cache
            .junctions
            .get_or_init(|| self.derive().junctions())
    }

    /// Disjoint exon bins per gene, or per aggregate of overlapping genes,
    /// cached.
    pub fn exon_bins(&self, aggregate: bool) -> &ExonBins {
        self.cache.bins[aggregate as usize].get_or_init(|| self.derive().exon_bins(aggregate))
    }

    /// Spliced length of every transcript (sum of its exon blocks), cached.
    pub fn transcript_lengths(&self) -> &[u32] {
        self.cache.transcript_lengths.get_or_init(|| {
            let mut blocks = Vec::new();
            (0..self.transcripts.len() as u32)
                .map(|t| {
                    self.exon_blocks(t, &mut blocks);
                    blocks.iter().map(|&(s, e)| e - s).sum()
                })
                .collect()
        })
    }

    /// Coding length of every transcript (merged CDS rows; stop codons
    /// outside CDS rows are not counted), cached.
    pub fn cds_lengths(&self) -> &[u32] {
        self.cache.cds_lengths.get_or_init(|| {
            let f = &self.features;
            let mut blocks = Vec::new();
            (0..self.transcripts.len() as u32)
                .map(|t| {
                    blocks.clear();
                    for i in self.transcripts.features(t) {
                        let i = i as usize;
                        if f.kind[i] == FeatureKind::Cds {
                            push_merged(&mut blocks, f.start[i], f.end[i]);
                        }
                    }
                    blocks.iter().map(|&(s, e)| e - s).sum()
                })
                .collect()
        })
    }
}

/// Appends `[s, e)` to sorted blocks, merging with the last on overlap or
/// contact.
#[inline]
fn push_merged(out: &mut Vec<(u32, u32)>, s: u32, e: u32) {
    if let Some(last) = out.last_mut()
        && s <= last.1
    {
        last.1 = last.1.max(e);
        return;
    }
    out.push((s, e));
}

/// Derived views of an annotation, optionally restricted to a
/// [`Selection`] and clipped to chromosome lengths. Built on each call.
#[derive(Clone, Copy)]
pub struct Derive<'a> {
    ann: &'a Annotation,
    selection: Option<&'a Selection>,
    lengths: Option<&'a [Option<u32>]>,
}

impl<'a> Derive<'a> {
    /// Only the selected genes and transcripts.
    pub fn selection(mut self, selection: &'a Selection) -> Self {
        self.selection = Some(selection);
        self
    }

    /// Clips windows to these lengths, by annotation seqid
    /// ([`Annotation::seqid_lengths`]); windows are always clipped at 0.
    pub fn clip(mut self, lengths: &'a [Option<u32>]) -> Self {
        self.lengths = Some(lengths);
        self
    }

    fn transcripts(&self) -> Box<dyn Iterator<Item = u32> + 'a> {
        match self.selection {
            Some(s) => Box::new(s.transcripts.iter().copied()),
            None => Box::new(0..self.ann.transcripts.len() as u32),
        }
    }

    fn genes(&self) -> Box<dyn Iterator<Item = u32> + 'a> {
        match self.selection {
            Some(s) => Box::new(s.genes.iter().copied()),
            None => Box::new(0..self.ann.genes.len() as u32),
        }
    }

    /// The transcripts of gene `g` that are in scope.
    fn gene_transcripts(&self, g: u32) -> impl Iterator<Item = u32> + 'a {
        let range = self.ann.genes.transcripts(g);
        let sel = self.selection;
        range.filter(move |t| sel.is_none_or(|s| s.transcripts.binary_search(t).is_ok()))
    }

    /// Introns between the exon blocks of each transcript (owner:
    /// transcript; name: transcript id).
    pub fn introns(&self) -> IntervalSet {
        let a = self.ann;
        let tx = &a.transcripts;
        let mut set = a.new_set(Some(Level::Transcript));
        let mut blocks = Vec::new();
        for t in self.transcripts() {
            a.exon_blocks(t, &mut blocks);
            let i = t as usize;
            for w in blocks.windows(2) {
                set.push(tx.seqid[i], w[0].1, w[1].0, tx.strand[i], t, tx.id[i]);
            }
        }
        set
    }

    /// The union of the exon blocks of each gene's transcripts (owner:
    /// gene; name: gene id).
    pub fn exonic_regions(&self) -> IntervalSet {
        let a = self.ann;
        let g = &a.genes;
        let mut set = a.new_set(Some(Level::Gene));
        let mut all = Vec::new();
        let mut blocks = Vec::new();
        let mut merged = Vec::new();
        for gene in self.genes() {
            all.clear();
            for t in self.gene_transcripts(gene) {
                a.exon_blocks(t, &mut blocks);
                all.extend_from_slice(&blocks);
            }
            all.sort_unstable();
            merged.clear();
            for &(s, e) in &all {
                push_merged(&mut merged, s, e);
            }
            let i = gene as usize;
            for &(s, e) in &merged {
                set.push(g.seqid[i], s, e, g.strand[i], gene, g.id[i]);
            }
        }
        set
    }

    /// The distinct introns with their genes and transcripts.
    pub fn junctions(&self) -> Junctions {
        let a = self.ann;
        let tx = &a.transcripts;
        let mut all: Vec<(u32, u32, u32, Strand, u32, u32)> = Vec::new();
        let mut blocks = Vec::new();
        for t in self.transcripts() {
            a.exon_blocks(t, &mut blocks);
            let i = t as usize;
            for w in blocks.windows(2) {
                all.push((tx.seqid[i], w[0].1, w[1].0, tx.strand[i], tx.gene[i], t));
            }
        }
        all.sort_unstable();
        let mut set = a.new_set(Some(Level::Gene));
        let mut genes = Csr::new();
        let mut transcripts = Csr::new();
        let mut lookup = id_map(all.len());
        let mut i = 0;
        while i < all.len() {
            let (s, st, en, strand, gene, _) = all[i];
            let mut j = i;
            let row = set.len() as u32;
            while j < all.len() && (all[j].0, all[j].1, all[j].2, all[j].3) == (s, st, en, strand) {
                j += 1;
            }
            let group = &all[i..j];
            let mut gs: Vec<u32> = group.iter().map(|x| x.4).collect();
            gs.dedup();
            genes.push_row(gs);
            transcripts.push_row(group.iter().map(|x| x.5));
            set.push(s, st, en, strand, gene, NONE);
            lookup.insert(jkey(s, st, en, strand), row);
            i = j;
        }
        Junctions {
            set,
            genes,
            transcripts,
            lookup,
        }
    }

    /// Genes in scope grouped into aggregates: genes whose exonic regions
    /// overlap on the same seqid and strand, transitively. Each aggregate
    /// is sorted, and aggregates are ordered by their first gene.
    fn aggregates(&self, exonic: &IntervalSet) -> Vec<Vec<u32>> {
        let mut parent: IdMap<u32, u32> = id_map(exonic.len());
        fn find(p: &mut IdMap<u32, u32>, x: u32) -> u32 {
            let mut r = x;
            while let Some(&q) = p.get(&r) {
                if q == r {
                    break;
                }
                r = q;
            }
            let mut y = x;
            while y != r {
                let q = p[&y];
                p.insert(y, r);
                y = q;
            }
            r
        }
        let mut rows: Vec<u32> = (0..exonic.len() as u32).collect();
        rows.sort_unstable_by_key(|&r| {
            let r = r as usize;
            (
                exonic.seqid[r],
                exonic.strand[r],
                exonic.start[r],
                exonic.end[r],
            )
        });
        let mut cluster: Option<(u32, Strand, u32, u32)> = None;
        for &r in &rows {
            let r = r as usize;
            let (s, strand, st, en, g) = (
                exonic.seqid[r],
                exonic.strand[r],
                exonic.start[r],
                exonic.end[r],
                exonic.owner[r],
            );
            parent.entry(g).or_insert(g);
            match &mut cluster {
                Some((cs, cstrand, cend, rep)) if *cs == s && *cstrand == strand && st < *cend => {
                    *cend = (*cend).max(en);
                    let a = find(&mut parent, g);
                    let b = find(&mut parent, *rep);
                    if a != b {
                        let (lo, hi) = (a.min(b), a.max(b));
                        parent.insert(hi, lo);
                    }
                }
                _ => cluster = Some((s, strand, en, g)),
            }
        }
        let mut groups: IdMap<u32, Vec<u32>> = id_map(parent.len());
        for g in self.genes() {
            let root = if parent.contains_key(&g) {
                find(&mut parent, g)
            } else {
                g
            };
            groups.entry(root).or_default().push(g);
        }
        let mut out: Vec<Vec<u32>> = groups.into_values().collect();
        for v in &mut out {
            v.sort_unstable();
        }
        out.sort_unstable_by_key(|v| v[0]);
        out
    }

    /// Disjoint exon bins per gene, or per aggregate of genes with
    /// overlapping exons on the same strand (`aggregate`).
    pub fn exon_bins(&self, aggregate: bool) -> ExonBins {
        let a = self.ann;
        let g = &a.genes;
        let units: Vec<Vec<u32>> = if aggregate {
            let exonic = self.exonic_regions();
            self.aggregates(&exonic)
        } else {
            self.genes().map(|x| vec![x]).collect()
        };
        let mut out = ExonBins {
            set: a.new_set(Some(Level::Gene)),
            part: Vec::new(),
            transcripts: Csr::new(),
            genes: Csr::new(),
        };
        let mut events: Vec<(u32, bool, u32)> = Vec::new();
        let mut blocks = Vec::new();
        let mut active: Vec<u32> = Vec::new();
        for unit in units {
            let first = unit[0] as usize;
            events.clear();
            for &gene in &unit {
                for t in self.gene_transcripts(gene) {
                    a.exon_blocks(t, &mut blocks);
                    for &(s, e) in &blocks {
                        events.push((s, true, t));
                        events.push((e, false, t));
                    }
                }
            }
            events.sort_unstable();
            active.clear();
            let mut part = 0;
            let mut prev = 0;
            let mut i = 0;
            let start_row = out.set.len();
            while i < events.len() {
                let p = events[i].0;
                if !active.is_empty() && prev < p {
                    let last = out.set.len();
                    let merge = last > start_row
                        && out.set.end[last - 1] == prev
                        && out.transcripts.get(last as u32 - 1) == active.as_slice();
                    if merge {
                        out.set.end[last - 1] = p;
                    } else {
                        part += 1;
                        out.set.push(
                            g.seqid[first],
                            prev,
                            p,
                            g.strand[first],
                            first as u32,
                            g.id[first],
                        );
                        out.part.push(part);
                        out.transcripts.push_row(active.iter().copied());
                        out.genes.push_row(unit.iter().copied());
                    }
                }
                while i < events.len() && events[i].0 == p {
                    let (_, add, t) = events[i];
                    match active.binary_search(&t) {
                        Ok(at) if !add => {
                            active.remove(at);
                        }
                        Err(at) if add => active.insert(at, t),
                        _ => {}
                    }
                    i += 1;
                }
                prev = p;
            }
        }
        out
    }

    fn clip_to(&self, seqid: u32, s: u32, e: u32) -> (u32, u32) {
        let len = self
            .lengths
            .and_then(|l| l.get(seqid as usize).copied().flatten())
            .unwrap_or(u32::MAX);
        let e = e.min(len);
        (s.min(e), e)
    }

    /// Rows at `level` (gene or transcript) as (row, seqid, start, end,
    /// strand, id).
    fn level_rows(&self, level: Level) -> Vec<(u32, u32, u32, u32, Strand, u32)> {
        let a = self.ann;
        match level {
            Level::Gene => {
                let g = &a.genes;
                self.genes()
                    .map(|r| {
                        let i = r as usize;
                        (r, g.seqid[i], g.start[i], g.end[i], g.strand[i], g.id[i])
                    })
                    .collect()
            }
            _ => {
                let t = &a.transcripts;
                self.transcripts()
                    .map(|r| {
                        let i = r as usize;
                        (r, t.seqid[i], t.start[i], t.end[i], t.strand[i], t.id[i])
                    })
                    .collect()
            }
        }
    }

    fn windows(&self, level: Level, f: impl Fn(u32, u32, Strand) -> (u32, u32)) -> IntervalSet {
        let level = if level == Level::Exon {
            Level::Transcript
        } else {
            level
        };
        let mut set = self.ann.new_set(Some(level));
        for (row, seqid, s, e, strand, id) in self.level_rows(level) {
            let (ws, we) = f(s, e, strand);
            let (ws, we) = self.clip_to(seqid, ws, we);
            set.push(seqid, ws, we, strand, row, id);
        }
        set
    }

    /// Windows around each TSS: `flank.upstream` bases before it and
    /// `flank.downstream` bases from it on (`GenomicRanges::promoters`).
    /// `level` is genes or transcripts (`Exon` means transcripts).
    pub fn promoters(&self, level: Level, flank: Flank) -> IntervalSet {
        self.windows(level, |s, e, strand| {
            if strand.is_reverse() {
                (
                    e.saturating_sub(flank.downstream),
                    e.saturating_add(flank.upstream),
                )
            } else {
                (
                    s.saturating_sub(flank.upstream),
                    s.saturating_add(flank.downstream),
                )
            }
        })
    }

    /// The first base of each gene or transcript, as a 1-base interval.
    pub fn tss(&self, level: Level) -> IntervalSet {
        self.promoters(level, Flank::new(0, 1))
    }

    /// Windows around each TES (the last transcribed base's end):
    /// `flank.upstream` bases inside the body and `flank.downstream` bases
    /// past it.
    pub fn tes_windows(&self, level: Level, flank: Flank) -> IntervalSet {
        self.windows(level, |s, e, strand| {
            if strand.is_reverse() {
                (
                    s.saturating_sub(flank.downstream),
                    s.saturating_add(flank.upstream),
                )
            } else {
                (
                    e.saturating_sub(flank.upstream),
                    e.saturating_add(flank.downstream),
                )
            }
        })
    }

    /// The last base of each gene or transcript, as a 1-base interval.
    pub fn tes(&self, level: Level) -> IntervalSet {
        self.tes_windows(level, Flank::new(1, 0))
    }

    /// Each body extended by `flank.upstream` bases before its TSS and
    /// `flank.downstream` bases past its TES: the "gene activity" regions
    /// of ArchR (5 kb upstream) and Signac (2 kb upstream).
    pub fn gene_bodies(&self, level: Level, flank: Flank) -> IntervalSet {
        self.windows(level, |s, e, strand| {
            if strand.is_reverse() {
                (
                    s.saturating_sub(flank.downstream),
                    e.saturating_add(flank.upstream),
                )
            } else {
                (
                    s.saturating_sub(flank.upstream),
                    e.saturating_add(flank.downstream),
                )
            }
        })
    }
}
