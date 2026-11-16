//! From raw records to the hierarchy: grouping (GTF by ids, GFF3 by
//! `ID` / `Parent`), synthesizing missing genes and transcripts, checks,
//! and sorting into the final tables.

use std::sync::Arc;

use omicsweave_core::{Error, Result};

use crate::alias::Aliases;
use crate::annotation::{Annotation, Cache};
use crate::model::{
    AttributePolicy, Attributes, Dialect, FeatureKind, Features, Genes, Others, PHASE_NONE, Row,
    RowRef, Strand, Tags, Transcripts,
};
use crate::parse::{GENE_BIOTYPE, GENE_ID, GENE_NAME, ID, Raw, TRANSCRIPT_ID, TX_BIOTYPE, TX_NAME};
use crate::report::{ProblemKind, Report};
use crate::strings::{IdMap, NONE, id_map};

struct GeneTmp {
    record: u32,
    seqid: u32,
    id: u32,
}

struct TxTmp {
    record: u32,
    seqid: u32,
    id: u32,
    gene: u32,
}

/// Genes, transcripts and (transcript, record) feature rows found by the
/// dialect's grouping.
#[derive(Default)]
struct Groups {
    genes: Vec<GeneTmp>,
    txs: Vec<TxTmp>,
    feats: Vec<(u32, u32)>,
}

#[inline]
fn key(seqid: u32, id: u32) -> u64 {
    ((seqid as u64) << 32) | id as u64
}

fn gene_like(t: &[u8]) -> bool {
    t == b"gene" || t == b"pseudogene" || t.ends_with(b"_gene")
}

pub(crate) struct BuildOptions<'a> {
    pub dialect: Dialect,
    pub policy: AttributePolicy,
    pub aliases: Option<&'a Aliases>,
    pub strict: bool,
}

pub(crate) fn build(mut raw: Raw, opts: BuildOptions<'_>) -> Result<Annotation> {
    let mut report = Report {
        records: raw.len(),
        comments: raw.comments.len(),
        ..Report::default()
    };
    let (seqid_names, seqid_of) = number_seqids(&mut raw, opts.aliases)?;
    let mut kinds: IdMap<u32, Option<FeatureKind>> = id_map(64);
    let kind_of: Vec<Option<FeatureKind>> = raw
        .feature_type
        .iter()
        .map(|&t| {
            *kinds
                .entry(t)
                .or_insert_with(|| FeatureKind::of_type(raw.pool.bytes(t)))
        })
        .collect();
    let groups = match opts.dialect {
        Dialect::Gtf => group_gtf(&raw, &mut report),
        Dialect::Gff3 => group_gff3(&raw, &kind_of, &mut report),
    };
    let mut ann = assemble(raw, groups, &kind_of, report, opts.dialect, opts.policy)?;
    ann.seqid_names = Arc::from(seqid_names);
    ann.seqid_of = seqid_of;
    ann.report.finish();
    if opts.strict
        && let Some(p) = ann.report.problems.first()
    {
        let id = ann.text(p.id);
        return Err(Error::Format {
            line: p.line,
            detail: if id.is_empty() {
                p.kind.describe().to_string()
            } else {
                format!("{} ({id})", p.kind.describe())
            },
        });
    }
    Ok(ann)
}

/// Replaces seqid names by numbers in first-appearance order (after
/// renaming through `aliases`).
fn number_seqids(raw: &mut Raw, aliases: Option<&Aliases>) -> Result<(Vec<u32>, IdMap<u32, u32>)> {
    let mut names: Vec<u32> = Vec::new();
    let mut of: IdMap<u32, u32> = id_map(64);
    let mut seen: IdMap<u32, u32> = id_map(64);
    for r in 0..raw.len() {
        let sid = raw.seqid[r];
        let s = match seen.get(&sid) {
            Some(&s) => s,
            None => {
                let name = match aliases {
                    Some(a) => {
                        let renamed = a.rename(raw.pool.str(sid)).to_string();
                        raw.pool.intern(renamed.as_bytes())?
                    }
                    None => sid,
                };
                let n = names.len() as u32;
                let s = *of.entry(name).or_insert(n);
                if s == n {
                    names.push(name);
                }
                seen.insert(sid, s);
                s
            }
        };
        raw.seqid[r] = s;
    }
    Ok((names, of))
}

fn group_gtf(raw: &Raw, report: &mut Report) -> Groups {
    let sid_gene = raw.pool.get(b"gene");
    let sid_tx = raw.pool.get(b"transcript");
    let mut g = Groups::default();
    let mut gene_map: IdMap<u64, u32> = id_map(1024);
    let mut tx_map: IdMap<u64, u32> = id_map(4096);
    for r in 0..raw.len() {
        let s = &raw.std[r];
        let seqid = raw.seqid[r];
        let t = Some(raw.feature_type[r]);
        let line = raw.line[r];
        let rec = r as u32;
        if s[GENE_ID] == NONE {
            report.push(ProblemKind::MissingGeneId, line, s[TRANSCRIPT_ID]);
            continue;
        }
        let mut gene_of = |g: &mut Groups| {
            let n = g.genes.len() as u32;
            let gene = *gene_map.entry(key(seqid, s[GENE_ID])).or_insert(n);
            if gene == n {
                g.genes.push(GeneTmp {
                    record: NONE,
                    seqid,
                    id: s[GENE_ID],
                });
            }
            gene
        };
        if t == sid_gene {
            let gene = gene_of(&mut g);
            let gt = &mut g.genes[gene as usize];
            if gt.record != NONE {
                report.push(ProblemKind::DuplicateId, line, s[GENE_ID]);
            } else {
                gt.record = rec;
            }
            continue;
        }
        if s[TRANSCRIPT_ID] == NONE {
            if t == sid_tx {
                report.push(ProblemKind::MissingTranscriptId, line, s[GENE_ID]);
            }
            continue;
        }
        let n = g.txs.len() as u32;
        let tx = *tx_map.entry(key(seqid, s[TRANSCRIPT_ID])).or_insert(n);
        if tx == n {
            let gene = gene_of(&mut g);
            g.txs.push(TxTmp {
                record: NONE,
                seqid,
                id: s[TRANSCRIPT_ID],
                gene,
            });
        } else if g.genes[g.txs[tx as usize].gene as usize].id != s[GENE_ID] {
            report.push(ProblemKind::GeneConflict, line, s[TRANSCRIPT_ID]);
        }
        if t == sid_tx {
            let tt = &mut g.txs[tx as usize];
            if tt.record != NONE {
                report.push(ProblemKind::DuplicateId, line, s[TRANSCRIPT_ID]);
            } else {
                tt.record = rec;
            }
        } else {
            g.feats.push((tx, rec));
        }
    }
    g
}

const E_CHILD: u8 = 1;
const NE_CHILD: u8 = 2;
const TX_CHILD: u8 = 4;
const STRONG_TX_CHILD: u8 = 8;
const TX: u8 = 16;
const STRONG: u8 = 32;
const GENE: u8 = 64;

fn group_gff3(raw: &Raw, kind_of: &[Option<FeatureKind>], report: &mut Report) -> Groups {
    let n = raw.len();
    let is_e = |r: usize| kind_of[r].is_some();
    let mut first: IdMap<u32, u32> = id_map(n / 2);
    let mut dup = vec![false; n];
    #[allow(clippy::needless_range_loop)]
    for r in 0..n {
        let id = raw.std[r][ID];
        if id == NONE {
            continue;
        }
        match first.get(&id) {
            Some(&f) => {
                // Lines sharing an ID are one discontinuous feature (RefSeq
                // CDS); only a repeated gene or transcript is a problem.
                if !is_e(r) || !is_e(f as usize) {
                    dup[r] = true;
                    report.push(ProblemKind::DuplicateId, raw.line[r], id);
                }
            }
            None => {
                first.insert(id, r as u32);
            }
        }
    }
    let mut parent = vec![NONE; n];
    let mut flags = vec![0u8; n];
    let mut parents_precede_children = true;
    for r in 0..n {
        if dup[r] {
            continue;
        }
        let source = if raw.std[r][ID] == NONE {
            r as u32
        } else {
            first[&raw.std[r][ID]]
        };
        for &p in raw.parents(r) {
            match first.get(&p) {
                Some(&pr) => {
                    parents_precede_children &= pr < source;
                    if parent[r] == NONE {
                        parent[r] = pr;
                    }
                    flags[pr as usize] |= if is_e(r) { E_CHILD } else { NE_CHILD };
                }
                None => report.push(ProblemKind::UnknownParent, raw.line[r], p),
            }
        }
    }
    if !parents_precede_children {
        check_parent_cycles(raw, &first, &dup, report);
    }
    // A transcript is the parent of exon-level features, or a leaf below
    // another non-feature record (a transcript without exons).
    for r in 0..n {
        if is_e(r) || dup[r] {
            continue;
        }
        let p = parent[r];
        if flags[r] & E_CHILD != 0 {
            flags[r] |= TX | STRONG;
        } else if p != NONE && !is_e(p as usize) && flags[r] & NE_CHILD == 0 {
            flags[r] |= TX;
        }
    }
    for r in 0..n {
        let p = parent[r];
        if flags[r] & TX != 0 && p != NONE {
            flags[p as usize] |= TX_CHILD;
            if flags[r] & STRONG != 0 {
                flags[p as usize] |= STRONG_TX_CHILD;
            }
        }
    }
    // A record with exons of its own and transcripts with exons below it
    // is a gene; its direct exons do not fit the hierarchy.
    for f in flags.iter_mut() {
        if *f & STRONG != 0 && *f & STRONG_TX_CHILD != 0 {
            *f &= !(TX | STRONG);
        }
    }
    for r in 0..n {
        if is_e(r) || dup[r] || flags[r] & TX != 0 {
            continue;
        }
        if flags[r] & TX_CHILD != 0
            || (parent[r] == NONE && gene_like(raw.pool.bytes(raw.feature_type[r])))
        {
            flags[r] |= GENE;
        }
    }
    let gene_id = |r: usize| {
        let s = &raw.std[r];
        if s[GENE_ID] != NONE {
            s[GENE_ID]
        } else {
            s[ID]
        }
    };
    let tx_id = |r: usize| {
        let s = &raw.std[r];
        if s[TRANSCRIPT_ID] != NONE {
            s[TRANSCRIPT_ID]
        } else {
            s[ID]
        }
    };
    let mut g = Groups::default();
    let mut gene_of = vec![NONE; n];
    let mut tx_of = vec![NONE; n];
    for r in 0..n {
        if flags[r] & GENE != 0 {
            gene_of[r] = g.genes.len() as u32;
            g.genes.push(GeneTmp {
                record: r as u32,
                seqid: raw.seqid[r],
                id: gene_id(r),
            });
        }
    }
    let mut synthesized: IdMap<u64, u32> = id_map(64);
    for r in 0..n {
        if flags[r] & TX == 0 {
            continue;
        }
        let seqid = raw.seqid[r];
        let line = raw.line[r];
        if raw.parents(r).len() > 1 {
            report.push(ProblemKind::MultipleParents, line, tx_id(r));
        }
        let p = parent[r];
        if p == NONE {
            if gene_like(raw.pool.bytes(raw.feature_type[r])) {
                let gene = g.genes.len() as u32;
                gene_of[r] = gene;
                g.genes.push(GeneTmp {
                    record: r as u32,
                    seqid,
                    id: gene_id(r),
                });
                tx_of[r] = g.txs.len() as u32;
                g.txs.push(TxTmp {
                    record: NONE,
                    seqid,
                    id: gene_id(r),
                    gene,
                });
            } else {
                let gid = gene_id(r);
                let next = g.genes.len() as u32;
                let gene = *synthesized.entry(key(seqid, gid)).or_insert(next);
                if gene == next {
                    g.genes.push(GeneTmp {
                        record: NONE,
                        seqid,
                        id: gid,
                    });
                }
                tx_of[r] = g.txs.len() as u32;
                g.txs.push(TxTmp {
                    record: r as u32,
                    seqid,
                    id: tx_id(r),
                    gene,
                });
            }
        } else if gene_of[p as usize] != NONE {
            if raw.seqid[p as usize] != seqid {
                report.push(ProblemKind::SeqidConflict, line, tx_id(r));
                continue;
            }
            tx_of[r] = g.txs.len() as u32;
            g.txs.push(TxTmp {
                record: r as u32,
                seqid,
                id: tx_id(r),
                gene: gene_of[p as usize],
            });
        } else {
            report.push(ProblemKind::BadNesting, line, tx_id(r));
        }
    }
    #[allow(clippy::needless_range_loop)]
    for r in 0..n {
        if !is_e(r) || dup[r] {
            continue;
        }
        let line = raw.line[r];
        let parents = raw.parents(r);
        if parents.is_empty() {
            report.push(ProblemKind::MissingParent, line, raw.std[r][ID]);
            continue;
        }
        for &pid in parents {
            let Some(&p) = first.get(&pid) else {
                continue;
            };
            let t = tx_of[p as usize];
            if t == NONE {
                report.push(ProblemKind::BadNesting, line, pid);
            } else if raw.seqid[p as usize] != raw.seqid[r] {
                report.push(ProblemKind::SeqidConflict, line, pid);
            } else {
                g.feats.push((t, r as u32));
            }
        }
    }
    g
}

/// A parent-first logical ID graph is acyclic by construction. For other orders,
/// iterative DFS checks every parent edge, including secondary parents,
/// without consuming call-stack depth on malformed hierarchies.
fn check_parent_cycles(raw: &Raw, first: &IdMap<u32, u32>, dup: &[bool], report: &mut Report) {
    let mut state = vec![0u8; raw.len()];
    // Discontinuous feature records share one logical ID and contribute
    // all of their edges. Link those rows without allocating per-ID Vecs.
    let mut next_record = vec![NONE; raw.len()];
    for r in 0..raw.len() {
        if dup[r] {
            continue;
        }
        let node = first.get(&raw.std[r][ID]).copied().unwrap_or(r as u32) as usize;
        if node != r {
            next_record[r] = next_record[node];
            next_record[node] = r as u32;
        }
    }
    let mut stack = Vec::new();
    for root in 0..raw.len() {
        if dup[root] || state[root] != 0 {
            continue;
        }
        // Non-first records are traversed through their logical node.
        if first
            .get(&raw.std[root][ID])
            .is_some_and(|&r| r as usize != root)
        {
            continue;
        }
        state[root] = 1;
        stack.push((root, root, 0));
        while let Some((node, r, edge)) = stack.last_mut() {
            let parents = raw.parents(*r);
            if *edge == parents.len() {
                if next_record[*r] != NONE {
                    *r = next_record[*r] as usize;
                    *edge = 0;
                    continue;
                }
                state[*node] = 2;
                stack.pop();
                continue;
            }
            let pid = parents[*edge];
            *edge += 1;
            let Some(&p) = first.get(&pid) else {
                continue;
            };
            let p = p as usize;
            match state[p] {
                0 => {
                    state[p] = 1;
                    stack.push((p, p, 0));
                }
                1 => report.push(ProblemKind::BadNesting, raw.line[*r], pid),
                _ => {}
            }
        }
    }
}

/// Counting sort of `items` by `group(item)` into CSR offsets.
fn bucket<T: Copy>(items: &[T], groups: usize, group: impl Fn(&T) -> u32) -> (Vec<u32>, Vec<T>) {
    let mut offsets = vec![0u32; groups + 1];
    for it in items {
        offsets[group(it) as usize + 1] += 1;
    }
    for i in 0..groups {
        offsets[i + 1] += offsets[i];
    }
    let mut at: Vec<u32> = offsets[..groups].to_vec();
    // Every destination is assigned exactly once by the disjoint bucket
    // cursors. A copied, initialized buffer avoids an Option per row and
    // a final unwrapping pass.
    let mut slots = items.to_vec();
    for it in items {
        let g = group(it) as usize;
        slots[at[g] as usize] = *it;
        at[g] += 1;
    }
    (offsets, slots)
}

fn assemble(
    mut raw: Raw,
    groups: Groups,
    kind_of: &[Option<FeatureKind>],
    mut report: Report,
    dialect: Dialect,
    policy: AttributePolicy,
) -> Result<Annotation> {
    let Groups {
        genes: gt,
        txs: tt,
        feats,
    } = groups;
    // Multi-parent GFF3 records can expand to more feature rows than
    // input records. RowRef reserves the upper two bits for its table.
    if feats.len() > crate::model::ROW_LIMIT as usize {
        return Err(Error::Unsupported("more than 2^30 annotation feature rows"));
    }
    let kind = |rec: u32| kind_of[rec as usize].unwrap_or(FeatureKind::Other);
    // Features of each transcript, by (start, end, kind, record).
    let (f_off, mut f_rows) = bucket(&feats, tt.len(), |&(t, _)| t);
    for t in 0..tt.len() {
        let rows = &mut f_rows[f_off[t] as usize..f_off[t + 1] as usize];
        rows.sort_unstable_by_key(|&(_, r)| {
            (raw.start[r as usize], raw.end[r as usize], kind(r), r)
        });
    }
    // Transcript spans and strands; synthesized ones from their features.
    let mut t_span = Vec::with_capacity(tt.len());
    for (t, tx) in tt.iter().enumerate() {
        let rows = &f_rows[f_off[t] as usize..f_off[t + 1] as usize];
        t_span.push(if tx.record != NONE {
            let r = tx.record as usize;
            (raw.start[r], raw.end[r], raw.strand[r])
        } else {
            let s = rows
                .iter()
                .map(|&(_, r)| raw.start[r as usize])
                .min()
                .unwrap_or(0);
            let e = rows
                .iter()
                .map(|&(_, r)| raw.end[r as usize])
                .max()
                .unwrap_or(0);
            let strand = rows
                .iter()
                .min_by_key(|&&(_, r)| r)
                .map_or(Strand::Unstranded, |&(_, r)| raw.strand[r as usize]);
            (s, e, strand)
        });
    }
    let tx_items: Vec<u32> = (0..tt.len() as u32).collect();
    let (t_off, mut t_rows) = bucket(&tx_items, gt.len(), |&t| tt[t as usize].gene);
    for g in 0..gt.len() {
        let rows = &mut t_rows[t_off[g] as usize..t_off[g + 1] as usize];
        rows.sort_unstable_by_key(|&t| (t_span[t as usize].0, t_span[t as usize].1, t));
    }
    let mut g_span = Vec::with_capacity(gt.len());
    for (g, gene) in gt.iter().enumerate() {
        let rows = &t_rows[t_off[g] as usize..t_off[g + 1] as usize];
        g_span.push(if gene.record != NONE {
            let r = gene.record as usize;
            (raw.start[r], raw.end[r], raw.strand[r])
        } else {
            let s = rows
                .iter()
                .map(|&t| t_span[t as usize].0)
                .min()
                .unwrap_or(0);
            let e = rows
                .iter()
                .map(|&t| t_span[t as usize].1)
                .max()
                .unwrap_or(0);
            let strand = rows
                .first()
                .map_or(Strand::Unstranded, |&t| t_span[t as usize].2);
            (s, e, strand)
        });
    }
    let mut gene_order: Vec<u32> = (0..gt.len() as u32).collect();
    gene_order.sort_unstable_by_key(|&g| {
        let (s, e, _) = g_span[g as usize];
        (gt[g as usize].seqid, s, e, g)
    });

    let sid_gene = raw.pool.intern(b"gene")?;
    let sid_tx = raw.pool.intern(b"transcript")?;
    let n_rec = raw.len();
    let mut rows = vec![RowRef(u32::MAX); n_rec];
    let mut genes = Genes {
        transcript_offsets: vec![0],
        ..Genes::default()
    };
    let mut txs = Transcripts {
        tag_offsets: vec![0],
        feature_offsets: vec![0],
        ..Transcripts::default()
    };
    let mut fs = Features::default();
    for &g in &gene_order {
        let gene = &gt[g as usize];
        let gi = genes.len() as u32;
        let (gs, ge, gstrand) = g_span[g as usize];
        let t_list = &t_rows[t_off[g as usize] as usize..t_off[g as usize + 1] as usize];
        let mut g_src = NONE;
        let mut flagged_strand = false;
        let mut flagged_outside = false;
        for &t in t_list {
            let tx = &tt[t as usize];
            let ti = txs.len() as u32;
            let (ts, te, tstrand) = t_span[t as usize];
            let fl = &f_rows[f_off[t as usize] as usize..f_off[t as usize + 1] as usize];
            let t_line = if tx.record != NONE {
                raw.line[tx.record as usize]
            } else {
                fl.first().map_or(0, |&(_, r)| raw.line[r as usize])
            };
            if tstrand != gstrand && !flagged_strand {
                report.push(ProblemKind::StrandConflict, t_line, tx.id);
                flagged_strand = true;
            }
            if (ts < gs || te > ge) && !flagged_outside {
                report.push(ProblemKind::OutsideParent, t_line, tx.id);
                flagged_outside = true;
            }
            let mut f_src = NONE;
            let mut f_strand = false;
            let mut f_outside = false;
            for &(_, r) in fl {
                let ri = r as usize;
                let fi = fs.len() as u32;
                if raw.strand[ri] != tstrand && !f_strand {
                    report.push(ProblemKind::StrandConflict, raw.line[ri], tx.id);
                    f_strand = true;
                }
                if (raw.start[ri] < ts || raw.end[ri] > te) && !f_outside {
                    report.push(ProblemKind::OutsideParent, raw.line[ri], tx.id);
                    f_outside = true;
                }
                if f_src == NONE {
                    f_src = raw.source[ri];
                }
                fs.transcript.push(ti);
                fs.kind.push(kind(r));
                fs.start.push(raw.start[ri]);
                fs.end.push(raw.end[ri]);
                fs.strand.push(raw.strand[ri]);
                fs.phase.push(raw.phase[ri]);
                fs.feature_type.push(raw.feature_type[ri]);
                fs.source.push(raw.source[ri]);
                fs.score.push(raw.score[ri]);
                fs.record.push(r);
                if rows[ri].0 == u32::MAX {
                    rows[ri] = RowRef::new(Row::Feature(fi));
                }
            }
            let recs = || {
                std::iter::once(tx.record)
                    .filter(|&r| r != NONE)
                    .chain(fl.iter().map(|&(_, r)| r))
            };
            let pick = |slot: usize| {
                recs()
                    .map(|r| raw.std[r as usize][slot])
                    .find(|&v| v != NONE)
                    .unwrap_or(NONE)
            };
            let tag_rec = recs().find(|&r| !raw.tags(r as usize).is_empty());
            let mut bits = Tags::default();
            if let Some(r) = tag_rec {
                for &v in raw.tags(r as usize) {
                    bits |= Tags::of_name(raw.pool.bytes(v));
                    txs.tag_values.push(v);
                }
            }
            txs.tag_offsets.push(
                u32::try_from(txs.tag_values.len())
                    .map_err(|_| Error::Unsupported("more than 2^32 transcript tag values"))?,
            );
            txs.tag_bits.push(bits);
            let (src, ftype, score, phase) = if tx.record != NONE {
                let r = tx.record as usize;
                rows[r] = RowRef::new(Row::Transcript(ti));
                (
                    raw.source[r],
                    raw.feature_type[r],
                    raw.score[r],
                    raw.phase[r],
                )
            } else {
                report.synthesized_transcripts += 1;
                (f_src, sid_tx, NONE, PHASE_NONE)
            };
            if g_src == NONE {
                g_src = src;
            }
            txs.seqid.push(tx.seqid);
            txs.start.push(ts);
            txs.end.push(te);
            txs.strand.push(tstrand);
            txs.gene.push(gi);
            txs.id.push(tx.id);
            txs.name.push(pick(TX_NAME));
            txs.biotype.push(pick(TX_BIOTYPE));
            txs.source.push(src);
            txs.feature_type.push(ftype);
            txs.score.push(score);
            txs.phase.push(phase);
            txs.record.push(tx.record);
            txs.feature_offsets.push(fs.len() as u32);
        }
        // Gene fields: the gene record first, then its transcripts' and
        // their features' records in order.
        let pick = |slot: usize| -> u32 {
            if gene.record != NONE {
                let v = raw.std[gene.record as usize][slot];
                if v != NONE {
                    return v;
                }
            }
            for &t in t_list {
                let tx = &tt[t as usize];
                if tx.record != NONE {
                    let v = raw.std[tx.record as usize][slot];
                    if v != NONE {
                        return v;
                    }
                }
                for &(_, r) in &f_rows[f_off[t as usize] as usize..f_off[t as usize + 1] as usize] {
                    let v = raw.std[r as usize][slot];
                    if v != NONE {
                        return v;
                    }
                }
            }
            NONE
        };
        let (src, ftype, score, phase) = if gene.record != NONE {
            let r = gene.record as usize;
            rows[r] = RowRef::new(Row::Gene(gi));
            (
                raw.source[r],
                raw.feature_type[r],
                raw.score[r],
                raw.phase[r],
            )
        } else {
            report.synthesized_genes += 1;
            (g_src, sid_gene, NONE, PHASE_NONE)
        };
        genes.seqid.push(gene.seqid);
        genes.start.push(gs);
        genes.end.push(ge);
        genes.strand.push(gstrand);
        genes.id.push(gene.id);
        genes.name.push(pick(GENE_NAME));
        genes.biotype.push(pick(GENE_BIOTYPE));
        genes.source.push(src);
        genes.feature_type.push(ftype);
        genes.score.push(score);
        genes.phase.push(phase);
        genes.record.push(gene.record);
        genes.transcript_offsets.push(txs.len() as u32);
    }
    if fs.len() >= crate::model::ROW_LIMIT as usize || txs.len() >= crate::model::ROW_LIMIT as usize
    {
        return Err(Error::Unsupported("more than 2^30 annotation rows"));
    }
    let mut others = Others::default();
    for (r, row) in rows.iter_mut().enumerate() {
        if row.0 != u32::MAX {
            continue;
        }
        *row = RowRef::new(Row::Other(others.len() as u32));
        let s = &raw.std[r];
        others.seqid.push(raw.seqid[r]);
        others.start.push(raw.start[r]);
        others.end.push(raw.end[r]);
        others.strand.push(raw.strand[r]);
        others.phase.push(raw.phase[r]);
        others.feature_type.push(raw.feature_type[r]);
        others.source.push(raw.source[r]);
        others.score.push(raw.score[r]);
        others
            .id
            .push(if s[ID] != NONE { s[ID] } else { s[GENE_ID] });
        others.record.push(r as u32);
    }
    report.genes = genes.len();
    report.transcripts = txs.len();
    report.features = fs.len();
    report.others = others.len();
    let attributes = Attributes {
        pairs: std::mem::take(&mut raw.pairs),
        ids: std::mem::take(&mut raw.attrs),
        offsets: std::mem::take(&mut raw.attr_offsets),
    };
    let comments = std::mem::take(&mut raw.comments);
    let mut pool = std::mem::take(&mut raw.pool);
    drop(raw);
    pool.shrink();
    Ok(Annotation {
        pool: Arc::new(pool),
        dialect,
        policy,
        seqid_names: Arc::from(Vec::new()),
        seqid_of: id_map(0),
        genes,
        transcripts: txs,
        features: fs,
        others,
        attributes,
        rows,
        comments,
        report,
        cache: Cache::default(),
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::parse::{Keep, parse_batch};

    #[test]
    fn cycle_check_unifies_discontinuous_feature_ids() {
        let raw = parse_batch(
            b"c\t.\tCDS\t1\t2\t.\t+\t0\tID=a\nc\t.\tCDS\t3\t4\t.\t+\t0\tID=b;Parent=a\nc\t.\tCDS\t5\t6\t.\t+\t0\tID=a;Parent=b\n",
            1, Dialect::Gff3, &Keep::All,
        ).unwrap();
        let mut first = id_map(raw.len());
        for r in 0..raw.len() {
            first.entry(raw.std[r][ID]).or_insert(r as u32);
        }
        let mut report = Report::default();
        check_parent_cycles(&raw, &first, &[false; 3], &mut report);
        assert_eq!(report.count(ProblemKind::BadNesting), 1);
        assert_eq!(report.problems[0].line, 2);
        assert_eq!(raw.pool.str(report.problems[0].id), "a");
        let ann = crate::Annotation::parse(
            b"##gff-version 3\nc\t.\tCDS\t1\t2\t.\t+\t0\tID=a\nc\t.\tCDS\t3\t4\t.\t+\t0\tID=b;Parent=a\nc\t.\tCDS\t5\t6\t.\t+\t0\tID=a;Parent=b\n",
            &crate::ReadOptions::default(),
        ).unwrap();
        assert_eq!(ann.report().count(ProblemKind::DuplicateId), 0);
        assert!(ann.report().count(ProblemKind::BadNesting) > 0);
    }
}
