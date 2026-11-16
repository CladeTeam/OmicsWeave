//! Single-cell ATAC: the fragments file and per-barcode metrics; the
//! matrices (peaks, tiles, gene activity) are [`Counts`] products the
//! [`Atac`] preset adds.
//!
//! Definitions (10x Genomics `fragments.tsv`, ArchR / Signac counting):
//! - A **fragment** is a read pair passing the filter, both mates mapped
//!   to the same reference, taken from its leftmost mate: from that mate's
//!   5' end + 4 to the other mate's 5' end - 5 (`pos + 4` to
//!   `pos + |TLEN| - 5`, 0-based half-open), the Tn5 shift.
//! - Its **cut sites** (Tn5 insertions) are its first base and its last
//!   base (`start`, `end - 1`).
//! - **Duplicates**: fragments of one barcode with the same reference,
//!   start and end are one fragment; the fifth column of the file counts
//!   the read pairs behind it.
//! - **Per-barcode metrics**: read pairs, fragments, fragments with a cut
//!   site in a peak and the cut sites in peaks, and cut sites within
//!   ±50 bp of a gene TSS (101 bp) and in the two 100 bp flanks 1900-2000
//!   bp away; ArchR's TSS enrichment of a barcode is
//!   `(tss_cuts / 101) / (tss_flank_cuts / 200)`.

use std::path::PathBuf;

use omicsweave_core::Result;
use omicsweave_gtf::Flank;

use super::counts::{Counts, FeatureSource, GlobalRows, LocalRows, Rows};
use super::coverage::Tally;
use super::{Accumulator, BatchInfo, Delivery, MemoryNeed, Needs, PlanInfo, RunCtx};
use crate::assign::{Query, SpaceId};
use crate::model::{BarcodeMetrics, BarcodeTable, Frag, FragmentTable};
use crate::prelude::{Filter, ReadView};
use crate::sink::json::Json;
use crate::spool::{Fixed, Spool, Spooled, sort_spooled};
use omicsweave_convert::stage::{Spill, read_or_end};
use std::io::Read;

/// The single-cell ATAC preset: fragments, metrics and the matrices asked
/// for, all over the barcodes of the same reads.
#[derive(Debug, Clone)]
pub struct Atac {
    /// The output directory (`atac`).
    pub name: String,
    pub barcode_tag: [u8; 2],
    /// Only these barcodes, rows in this order.
    pub whitelist: Option<PathBuf>,
    pub filter: Filter,
    /// Fragments longer than this are dropped (0: no limit).
    pub max_fragment: u32,
    /// Write `fragments.tsv.gz` with its tabix index.
    pub fragments: bool,
    /// Peaks (bigBed): a peak matrix and peak metrics.
    pub peaks: Option<PathBuf>,
    /// A genome tile matrix of this bin size.
    pub tiles: Option<u32>,
    /// A gene activity matrix of fragments over gene bodies with these
    /// flanks (needs an annotation).
    pub gene_activity: Option<Flank>,
}

impl Default for Atac {
    /// Proper pairs with MAPQ ≥ 30, primary, not QC-failed (duplicates are
    /// kept and collapsed), 5 kb tiles, and gene activity over each gene
    /// body with 2 kb upstream (Signac) when there is an annotation.
    fn default() -> Atac {
        use omicsweave_bam::flag;
        Atac {
            name: "atac".into(),
            barcode_tag: *b"CB",
            whitelist: None,
            filter: Filter {
                require: flag::PAIRED | flag::PROPER_PAIR,
                exclude: flag::UNMAPPED
                    | flag::MATE_UNMAPPED
                    | flag::SECONDARY
                    | flag::QC_FAIL
                    | flag::SUPPLEMENTARY,
                min_mapq: 30,
                unique: false,
            },
            max_fragment: 0,
            fragments: true,
            peaks: None,
            tiles: Some(5000),
            gene_activity: Some(Flank::new(2000, 0)),
        }
    }
}

impl Atac {
    /// The count products of the preset.
    pub fn matrices(&self, has_annotation: bool) -> Vec<Counts> {
        let mut out = Vec::new();
        let mk = |sub: &str, features: FeatureSource| {
            let mut c = Counts::atac(&format!("{}/{sub}", self.name), features, self.filter);
            c.rows = super::counts::RowKey::Tag(self.barcode_tag);
            c.whitelist = self.whitelist.clone();
            c.max_fragment = self.max_fragment;
            c
        };
        if let Some(p) = &self.peaks {
            out.push(mk("peaks", FeatureSource::BigBed(p.clone())));
        }
        if let Some(b) = self.tiles {
            out.push(mk("tiles", FeatureSource::Bins(b)));
        }
        if let (Some(f), true) = (self.gene_activity, has_annotation) {
            let mut c = mk("gene_activity", FeatureSource::GeneActivity(f));
            c.unit = super::counts::CountUnit::Fragment { tn5: true };
            c.dedup = super::counts::Dedup::Position;
            out.push(c);
        }
        out
    }
}

/// A fragment of a batch, before duplicates are collapsed.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
struct BFrag {
    /// `tid << 32 | start`.
    key: u64,
    end: u32,
    row: u32,
    count: u32,
    peak_cuts: u8,
    tss_cuts: u8,
    flank_cuts: u8,
}

impl BFrag {
    fn order(&self) -> (u64, u32, u32) {
        (self.key, self.end, self.row)
    }
}

impl Fixed for BFrag {
    const SIZE: usize = 23;
    fn put(&self, out: &mut [u8]) {
        out[..8].copy_from_slice(&self.key.to_le_bytes());
        for (at, value) in [(8, self.end), (12, self.row), (16, self.count)] {
            out[at..at + 4].copy_from_slice(&value.to_le_bytes());
        }
        out[20] = self.peak_cuts;
        out[21] = self.tss_cuts;
        out[22] = self.flank_cuts;
    }
    fn get(b: &[u8]) -> Self {
        let u = |i| u32::from_le_bytes(b[i..i + 4].try_into().unwrap());
        Self {
            key: u64::from_le_bytes(b[..8].try_into().unwrap()),
            end: u(8),
            row: u(12),
            count: u(16),
            peak_cuts: b[20],
            tss_cuts: b[21],
            flank_cuts: b[22],
        }
    }
}

impl Spill for BFrag {
    type Key = (u64, u32, u32);
    fn key(&self) -> Self::Key {
        self.order()
    }
    fn mem_size(&self) -> usize {
        0
    }
    fn encode(&self, out: &mut Vec<u8>) {
        let mut b = [0u8; 23];
        self.put(&mut b);
        out.extend_from_slice(&b);
    }
    fn decode(input: &mut impl Read) -> Result<Option<Self>> {
        let mut b = [0u8; 23];
        Ok(read_or_end(input, &mut b)?.then(|| Self::get(&b)))
    }
}

/// Sorts and folds duplicates.
fn collapse(v: &mut Vec<BFrag>) {
    v.sort_unstable_by_key(BFrag::order);
    let mut w = 0;
    for i in 0..v.len() {
        if w > 0 && v[w - 1].order() == v[i].order() {
            v[w - 1].count += v[i].count;
        } else {
            v[w] = v[i];
            w += 1;
        }
    }
    v.truncate(w);
}

/// The fragments stage of [`Atac`], resolved against a task.
pub struct FragmentStage {
    pub(crate) spec: Atac,
    pub(crate) rows: Rows,
    pub(crate) peaks: Option<SpaceId>,
    /// TSS centres and flanks.
    pub(crate) tss: Option<(SpaceId, SpaceId)>,
    pub(crate) table_bytes: u64,
}

pub struct FragmentBatch {
    frags: Vec<BFrag>,
    local: LocalRows,
    tally: Tally,
    no_row: u64,
}

pub struct FragmentTotal {
    held: Vec<BFrag>,
    scratch: Vec<BFrag>,
    rows: GlobalRows,
    spool: Spool<Frag>,
    metrics: Vec<[u64; 6]>,
    tally: Tally,
    no_row: u64,
    deferred: Option<Spool<BFrag>>,
}

const METRICS: [&str; 6] = [
    "read_pairs",
    "fragments",
    "peak_fragments",
    "peak_cuts",
    "tss_cuts",
    "tss_flank_cuts",
];

impl FragmentTotal {
    fn emit(&mut self, f: &BFrag) -> Result<()> {
        let row = f.row as usize;
        if self.metrics.len() <= row {
            self.metrics.resize(row + 1, [0; 6]);
        }
        let m = &mut self.metrics[row];
        m[0] += f.count as u64;
        m[1] += 1;
        m[2] += (f.peak_cuts > 0) as u64;
        m[3] += f.peak_cuts as u64;
        m[4] += f.tss_cuts as u64;
        m[5] += f.flank_cuts as u64;
        self.spool.push(Frag {
            tid: (f.key >> 32) as u32,
            start: f.key as u32,
            end: f.end,
            row: f.row,
            count: f.count,
        })
    }

    /// Merges `next` (sorted, collapsed) with the held fragments and emits
    /// those starting before `limit`.
    fn advance(&mut self, next: &[BFrag], limit: u64, rc: &RunCtx<'_>) -> Result<()> {
        let retained = self.held.len() - self.held.partition_point(|f| f.key < limit) + next.len()
            - next.partition_point(|f| f.key < limit);
        if self.deferred.is_none()
            && retained.saturating_mul(size_of::<BFrag>()) as u64 > (rc.buffer / 4).max(64 << 10)
        {
            self.spool.stream_to_disk()?;
            let mut deferred = Spool::new((rc.buffer / 8).max(8 << 10), &rc.tmp);
            for f in std::mem::take(&mut self.held) {
                deferred.push(f)?;
            }
            self.scratch = Vec::new();
            self.deferred = Some(deferred);
        }
        if let Some(deferred) = &mut self.deferred {
            for &f in next {
                deferred.push(f)?;
            }
            return Ok(());
        }
        let held = std::mem::take(&mut self.held);
        let mut all = std::mem::take(&mut self.scratch);
        all.clear();
        let (mut i, mut j) = (0, 0);
        while i < held.len() || j < next.len() {
            let take_held =
                j == next.len() || (i < held.len() && held[i].order() <= next[j].order());
            let f = if take_held {
                i += 1;
                held[i - 1]
            } else {
                j += 1;
                next[j - 1]
            };
            match all.last_mut() {
                Some(last) if last.order() == f.order() => last.count += f.count,
                _ => all.push(f),
            }
        }
        let cut = all.partition_point(|f| f.key < limit);
        for f in &all[..cut] {
            self.emit(f)?;
        }
        let mut rest = held;
        rest.clear();
        rest.extend_from_slice(&all[cut..]);
        self.held = rest;
        self.scratch = all;
        Ok(())
    }
}

fn cuts(r: &ReadView<'_>, space: SpaceId) -> u8 {
    (!r.assign(space, Query::CutStart, None).is_empty()) as u8
        + (!r.assign(space, Query::CutEnd, None).is_empty()) as u8
}

impl Accumulator for FragmentStage {
    type Batch = FragmentBatch;
    type Total<'p> = FragmentTotal;
    type Product = (FragmentTable, BarcodeMetrics, Tally, u64);

    fn name(&self) -> &str {
        &self.spec.name
    }

    fn needs(&self) -> Needs {
        Needs::CIGAR | Needs::TAGS | Needs::ASSIGN
    }

    fn estimate(&self, info: &PlanInfo) -> MemoryNeed {
        MemoryNeed {
            fixed: self.table_bytes
                + self.rows.fixed_names().map_or(0, |names| {
                    (names.len() as u64).saturating_mul(2 * size_of::<[u64; 6]>() as u64)
                }),
            per_batch: info.batch_bytes / 8,
            min_buffer: 2 << 20,
            weight: 1,
        }
    }

    fn batch(&self, _info: BatchInfo) -> FragmentBatch {
        FragmentBatch {
            frags: Vec::new(),
            local: LocalRows::default(),
            tally: Tally::default(),
            no_row: 0,
        }
    }

    fn observe(&self, b: &mut FragmentBatch, r: &ReadView<'_>) -> Result<()> {
        b.tally.seen += 1;
        if let Some(why) = self.spec.filter.reject(r) {
            b.tally.reject(why);
            return Ok(());
        }
        let Some(f) = r.atac_fragment(self.spec.max_fragment) else {
            return Ok(());
        };
        let Some(row) = self.rows_of(r, &mut b.local) else {
            b.no_row += 1;
            return Ok(());
        };
        let (tss_cuts, flank_cuts) = match self.tss {
            Some((c, fl)) => (cuts(r, c), cuts(r, fl)),
            None => (0, 0),
        };
        b.frags.push(BFrag {
            key: (r.tid() as u64) << 32 | f.start as u64,
            end: f.end,
            row,
            count: 1,
            peak_cuts: self.peaks.map_or(0, |p| cuts(r, p)),
            tss_cuts,
            flank_cuts,
        });
        Ok(())
    }

    fn seal(&self, b: &mut FragmentBatch) -> Result<()> {
        collapse(&mut b.frags);
        Ok(())
    }

    fn total<'p>(&'p self, rc: &RunCtx<'p>) -> Result<FragmentTotal> {
        Ok(FragmentTotal {
            held: Vec::new(),
            scratch: Vec::new(),
            rows: GlobalRows::default(),
            spool: Spool::new(rc.buffer, &rc.tmp),
            deferred: None,
            metrics: Vec::new(),
            tally: Tally::default(),
            no_row: 0,
        })
    }

    fn merge<'p>(
        &'p self,
        t: &mut FragmentTotal,
        b: FragmentBatch,
        frontier: u64,
        rc: &RunCtx<'p>,
    ) -> Result<()> {
        t.tally.add(&b.tally);
        t.no_row += b.no_row;
        let FragmentBatch {
            mut frags, local, ..
        } = b;
        if matches!(self.rows, Rows::Tag(_)) {
            let map = t.rows.remap(local);
            for f in &mut frags {
                f.row = map[f.row as usize];
            }
            if t.rows.bytes() > rc.buffer.max(64 << 20) {
                return Err(omicsweave_core::Error::invalid(
                    "budget",
                    "ATAC row keys outgrow the memory budget; give a whitelist or raise the memory",
                ));
            }
            // Renumbering reorders only fragments of equal position.
            let mut i = 0;
            while i < frags.len() {
                let mut j = i + 1;
                while j < frags.len()
                    && (frags[j].key, frags[j].end) == (frags[i].key, frags[i].end)
                {
                    j += 1;
                }
                if j - i > 1 {
                    frags[i..j].sort_unstable_by_key(|f| f.row);
                }
                i = j;
            }
        }
        // Later records start at or after the frontier's position, so their
        // fragments start at least 4 bases further.
        let limit = frontier.saturating_sub(1).saturating_add(4);
        t.advance(&frags, limit, rc)
    }

    fn end_input<'p>(&'p self, t: &mut FragmentTotal, _input: u32, rc: &RunCtx<'p>) -> Result<()> {
        if let Some(deferred) = t.deferred.take() {
            let mut pending: Option<BFrag> = None;
            let spilled = sort_spooled(
                deferred,
                rc.pool,
                (rc.buffer / 2).max(1),
                &rc.tmp,
                |batch| {
                    for &f in batch {
                        if let Some(last) = &mut pending
                            && last.order() == f.order()
                        {
                            last.count += f.count;
                            continue;
                        }
                        if let Some(last) = pending.replace(f) {
                            t.emit(&last)?;
                        }
                    }
                    Ok(())
                },
            )?;
            if let Some(last) = pending {
                t.emit(&last)?;
            }
            t.spool.account_spilled(spilled);
        }
        t.advance(&[], u64::MAX, rc)
    }

    fn finish<'p>(
        &'p self,
        t: FragmentTotal,
        rc: &RunCtx<'p>,
    ) -> Result<(FragmentTable, BarcodeMetrics, Tally, u64)> {
        let names = match self.rows.fixed_names() {
            Some(n) => n.to_vec(),
            None => t.rows.names,
        };
        let mut metrics = t.metrics;
        metrics.resize(names.len(), [0; 6]);
        let columns = (0..6)
            .map(|c| {
                (
                    METRICS[c].to_string(),
                    metrics.iter().map(|m| m[c]).collect(),
                )
            })
            .collect();
        let barcodes = BarcodeTable { names };
        let fragments: Spooled<Frag> = t.spool.finish()?;
        Ok((
            FragmentTable {
                chroms: rc.inputs[0].refs.clone(),
                barcodes: barcodes.clone(),
                fragments,
            },
            BarcodeMetrics { barcodes, columns },
            t.tally,
            t.no_row,
        ))
    }

    fn deliver(
        &self,
        (table, metrics, tally, no_row): (FragmentTable, BarcodeMetrics, Tally, u64),
        d: &mut Delivery<'_>,
    ) -> Result<Json> {
        let mut files = Vec::new();
        let n = table.fragments.len();
        let spilled = table.fragments.spilled();
        if self.spec.fragments {
            let rel = format!("{}/fragments.tsv.gz", self.spec.name);
            let path = d.path(&rel)?;
            crate::sink::fragments::write(table, &path, d)?;
            files.push(Json::from(rel.clone()));
            files.push(Json::from(format!("{rel}.tbi")));
        }
        let rel = format!("{}/barcode_metrics.csv", self.spec.name);
        let path = d.path(&rel)?;
        crate::sink::table::write_metrics(&metrics, &path, d)?;
        files.push(Json::from(rel));
        let mut fields = vec![
            ("product", Json::from("atac-fragments")),
            ("fragments", Json::from(n)),
            ("barcodes", Json::from(metrics.barcodes.len() as u64)),
            ("unassigned_no_barcode", Json::from(no_row)),
            ("spilled_bytes", Json::from(spilled)),
        ];
        fields.extend(tally.json());
        fields.push(("files", Json::Arr(files)));
        Ok(Json::obj(fields))
    }
}

impl FragmentStage {
    fn rows_of(&self, r: &ReadView<'_>, local: &mut LocalRows) -> Option<u32> {
        match &self.rows {
            Rows::Tag(t) => r.tag_bytes(*t).map(|v| local.id(v)),
            Rows::Listed(t, ids, _) => ids.get(r.tag_bytes(*t)?).copied(),
            _ => None,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::source::Input;
    use omicsweave_core::exec::Pool;
    use std::collections::BTreeMap;
    use std::sync::Arc;

    #[test]
    fn same_start_fragment_piles_use_bounded_deferred_sorting() {
        let dir = tempfile::tempdir().unwrap();
        let pool = Pool::new(2);
        let inputs = [Input {
            path: dir.path().join("in.bam"),
            refs: vec![("chr1".into(), 200_000)],
            read_groups: vec![],
            len: 0,
        }];
        let small = RunCtx {
            pool: &pool,
            buffer: 1 << 20,
            tmp: dir.path().to_owned(),
            inputs: &inputs,
        };
        let large = RunCtx {
            pool: &pool,
            buffer: 128 << 20,
            tmp: dir.path().to_owned(),
            inputs: &inputs,
        };
        let stage = FragmentStage {
            spec: Atac::default(),
            rows: Rows::File(Arc::new(
                (0..4).map(|i| format!("B{i}").into_bytes()).collect(),
            )),
            peaks: None,
            tss: None,
            table_bytes: 0,
        };
        assert_eq!(stage.spec.max_fragment, 0);
        let mut a = stage.total(&small).unwrap();
        let mut b = stage.total(&large).unwrap();
        let mut expected: BTreeMap<(u64, u32, u32), BFrag> = BTreeMap::new();
        for batch in 0..16u32 {
            let mut values = Vec::new();
            for i in 0..2000u32 {
                let offset = 1 + (i * 391 + batch * 127) % 20_000;
                let f = BFrag {
                    key: 100_004,
                    end: 100_004 + offset,
                    row: offset % 4,
                    count: 1 + i % 2,
                    peak_cuts: (offset % 3) as u8,
                    tss_cuts: (offset % 2) as u8,
                    flank_cuts: 0,
                };
                expected
                    .entry(f.order())
                    .and_modify(|e| e.count += f.count)
                    .or_insert(f);
                values.push(f);
            }
            if batch == 0 {
                let f = BFrag {
                    key: 10,
                    end: 20,
                    row: 0,
                    count: 1,
                    peak_cuts: 0,
                    tss_cuts: 0,
                    flank_cuts: 0,
                };
                expected.insert(f.order(), f);
                values.push(f);
            }
            collapse(&mut values);
            a.advance(&values, 100_004, &small).unwrap();
            b.advance(&values, 100_004, &large).unwrap();
        }
        assert!(a.deferred.is_some());
        assert_eq!(a.held.capacity(), 0);
        assert_eq!(a.scratch.capacity(), 0);
        assert!(b.deferred.is_none());
        stage.end_input(&mut a, 0, &small).unwrap();
        stage.end_input(&mut b, 0, &large).unwrap();
        let want: Vec<_> = expected
            .values()
            .map(|f| (f.key, f.end, f.row, f.count))
            .collect();
        let mut pair_counts = vec![0u64; 4];
        for f in expected.values() {
            pair_counts[f.row as usize] += f.count as u64;
        }
        for (state, rc) in [(a, &small), (b, &large)] {
            let (table, metrics, _, _) = stage.finish(state, rc).unwrap();
            assert_eq!(metrics.columns[0].1, pair_counts);
            let mut got = Vec::new();
            table
                .fragments
                .drain(|f| {
                    got.push((
                        ((f.tid as u64) << 32) | f.start as u64,
                        f.end,
                        f.row,
                        f.count,
                    ));
                    Ok(())
                })
                .unwrap();
            assert_eq!(got, want);
        }
    }
}
