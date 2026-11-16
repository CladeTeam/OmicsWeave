//! Transcript regions, bins and stranded windows against base-set oracles.
use crate::support::{
    LEN,
    base_sets::{exon_bits, naive_introns, rows, runs},
    random, read,
};
use omicsweave_gtf::{Annotation, Flank, Level, Strand};
use std::collections::{BTreeMap, BTreeSet};
type JunctionKey = (u32, u32, u32, Strand);

#[test]
fn introns_exonic_regions_and_lengths() {
    for seed in 0..8 {
        let a = random(seed);
        let tx = a.transcripts();
        let mut want = Vec::new();
        for t in 0..tx.len() as u32 {
            let i = t as usize;
            for (s, e) in naive_introns(&a, t) {
                want.push((tx.seqid[i], s, e, tx.strand[i], t));
            }
            let len: u32 = exon_bits(&a, t).iter().filter(|&&b| b).count() as u32;
            assert_eq!(a.transcript_lengths()[i], len);
        }
        assert_eq!(rows(a.introns()), want);

        let g = a.genes();
        let mut want = Vec::new();
        for gene in 0..g.len() as u32 {
            let mut bits = vec![false; LEN as usize];
            for t in g.transcripts(gene) {
                for (b, e) in bits.iter_mut().zip(exon_bits(&a, t)) {
                    *b |= e;
                }
            }
            let i = gene as usize;
            for (s, e) in runs(&bits) {
                want.push((g.seqid[i], s, e, g.strand[i], gene));
            }
        }
        assert_eq!(rows(a.exonic_regions()), want);
    }
}

#[test]
fn junctions_collect_genes_and_transcripts() {
    for seed in 0..8 {
        let a = random(seed);
        let tx = a.transcripts();
        let mut want: BTreeMap<JunctionKey, (BTreeSet<u32>, Vec<u32>)> = BTreeMap::new();
        for t in 0..tx.len() as u32 {
            let i = t as usize;
            for (s, e) in naive_introns(&a, t) {
                let v = want.entry((tx.seqid[i], s, e, tx.strand[i])).or_default();
                v.0.insert(tx.gene[i]);
                v.1.push(t);
            }
        }
        let j = a.junctions();
        assert_eq!(j.len(), want.len());
        for (k, ((s, st, en, strand), (genes, txs))) in want.iter().enumerate() {
            let k = k as u32;
            let i = k as usize;
            assert_eq!(
                (
                    j.set.seqid[i],
                    j.set.start[i],
                    j.set.end[i],
                    j.set.strand[i]
                ),
                (*s, *st, *en, *strand)
            );
            assert_eq!(j.genes.get(k), genes.iter().copied().collect::<Vec<_>>());
            assert_eq!(j.transcripts.get(k), txs.as_slice());
            assert_eq!(j.find(*s, *st, *en, *strand), Some(k));
            assert_eq!(
                j.find_any(*s, *st, *en).map(|x| j.set.start[x as usize]),
                Some(*st)
            );
            let (d, ac) = (j.donor(k), j.acceptor(k));
            if *strand == Strand::Reverse {
                assert_eq!((d, ac), (en - 1, *st));
            } else {
                assert_eq!((d, ac), (*st, en - 1));
            }
        }
        assert_eq!(j.find(0, 1, 2, Strand::Forward), None);
    }
}

/// Bins of a unit of genes: runs of bases with the same non-empty set of
/// covering transcripts.
fn naive_bins(a: &Annotation, genes: &[u32]) -> Vec<(u32, u32, Vec<u32>)> {
    let mut cover: Vec<Vec<u32>> = vec![Vec::new(); LEN as usize];
    let mut txs: Vec<u32> = genes
        .iter()
        .flat_map(|&g| a.genes().transcripts(g))
        .collect();
    txs.sort_unstable();
    for &t in &txs {
        for (b, e) in cover.iter_mut().zip(exon_bits(a, t)) {
            if e {
                b.push(t);
            }
        }
    }
    let mut out: Vec<(u32, u32, Vec<u32>)> = Vec::new();
    for (p, set) in cover.into_iter().enumerate() {
        if set.is_empty() {
            continue;
        }
        let p = p as u32;
        if let Some(last) = out.last_mut()
            && last.1 == p
            && last.2 == set
        {
            last.1 += 1;
            continue;
        }
        out.push((p, p + 1, set));
    }
    out
}

#[test]
fn transcript_parts_merge_codons_and_select_utrs_by_strand() {
    use omicsweave_gtf::TranscriptPart;

    for strand in ['+', '-', '.', '?'] {
        for explicit_exons in [false, true] {
            let mut text = String::new();
            use std::fmt::Write;
            let mut pieces = vec![
                ("UTR", 1, 5, '.'),
                ("CDS", 6, 10, '2'),
                ("CDS", 21, 24, '0'),
                ("stop_codon", 25, 27, '0'),
                ("UTR", 28, 30, '.'),
            ];
            if explicit_exons {
                pieces.extend([("exon", 1, 10, '.'), ("exon", 21, 30, '.')]);
            }
            for (kind, start, end, phase) in pieces.into_iter().rev() {
                writeln!(text, "c\t.\t{kind}\t{start}\t{end}\t.\t{strand}\t{phase}\tgene_id \"g\"; transcript_id \"t\";").unwrap();
            }
            let ann = read(&text);
            let mut out = vec![(99, 100)];
            for (part, expected) in [
                (TranscriptPart::Exons, vec![(0, 10), (20, 30)]),
                (TranscriptPart::Cds, vec![(5, 10), (20, 27)]),
                (
                    TranscriptPart::Utr5,
                    if strand == '-' {
                        vec![(27, 30)]
                    } else {
                        vec![(0, 5)]
                    },
                ),
                (
                    TranscriptPart::Utr3,
                    if strand == '-' {
                        vec![(0, 5)]
                    } else {
                        vec![(27, 30)]
                    },
                ),
            ] {
                ann.part_blocks(0, part, &mut out);
                assert_eq!(out, expected);
            }
            assert_eq!(ann.cds_lengths(), [9]);
            assert_eq!(ann.transcript_lengths(), [20]);
        }
    }
    let noncoding = read("c\t.\texon\t1\t10\t.\t+\t.\tgene_id \"g\"; transcript_id \"t\";\n");
    for part in [
        TranscriptPart::Cds,
        TranscriptPart::Utr5,
        TranscriptPart::Utr3,
    ] {
        let mut out = vec![(99, 100)];
        noncoding.part_blocks(0, part, &mut out);
        assert!(out.is_empty());
    }
}

#[test]
fn shared_exon_boundaries_and_sparse_transcripts_match_base_sets() {
    use omicsweave_gtf::Selection;
    use std::fmt::Write;

    let mut text = String::new();
    for i in 0..128 {
        // Shared starts/ends exercise bulk updates; the additional pieces
        // interleave removals and additions across transcript identities.
        for (s, e) in [(1, 20), (41 + i % 8, 60 + i % 8), (81, 100)] {
            writeln!(
                text,
                "c\t.\texon\t{s}\t{e}\t.\t+\t.\tgene_id \"g\"; transcript_id \"t{i}\";"
            )
            .unwrap();
        }
    }
    let ann = read(&text);
    let bins = ann.exon_bins(false);
    let transcript_bits: Vec<_> = (0..ann.transcripts().len() as u32)
        .map(|t| exon_bits(&ann, t))
        .collect();
    let want = naive_bins(&ann, &[0]);
    let got: Vec<_> = (0..bins.set.len())
        .map(|i| {
            (
                bins.set.start[i],
                bins.set.end[i],
                bins.transcripts.get(i as u32).to_vec(),
            )
        })
        .collect();
    assert_eq!(got, want);
    for transcripts in [vec![], vec![0], vec![0, 0, 63, 127], (0..128).collect()] {
        let selection = Selection {
            genes: vec![0],
            transcripts: transcripts.clone(),
        };
        let sparse = ann.derive().selection(&selection).exon_bins(false);
        for (p, _) in transcript_bits[0].iter().take(100).enumerate() {
            let mut want: Vec<_> = transcripts
                .iter()
                .copied()
                .filter(|&t| transcript_bits[t as usize][p])
                .collect();
            want.dedup();
            let row = (0..sparse.set.len())
                .find(|&i| sparse.set.start[i] <= p as u32 && (p as u32) < sparse.set.end[i]);
            assert_eq!(
                row.map_or(&[][..], |i| sparse.transcripts.get(i as u32)),
                want
            );
        }
    }
}

#[test]
fn exon_bins_per_gene_and_aggregate() {
    for seed in 0..6 {
        let a = random(seed);
        let g = a.genes();
        let bins = a.exon_bins(false);
        let mut want = Vec::new();
        for gene in 0..g.len() as u32 {
            for (n, (s, e, txs)) in naive_bins(&a, &[gene]).into_iter().enumerate() {
                want.push((gene, n as u32 + 1, s, e, txs));
            }
        }
        let got: Vec<_> = (0..bins.set.len())
            .map(|i| {
                (
                    bins.set.owner[i],
                    bins.part[i],
                    bins.set.start[i],
                    bins.set.end[i],
                    bins.transcripts.get(i as u32).to_vec(),
                )
            })
            .collect();
        assert_eq!(got, want);

        // Aggregates: genes sharing exonic bases on the same seqid and
        // strand, transitively.
        let n = g.len();
        let mut parent: Vec<usize> = (0..n).collect();
        fn find(p: &mut [usize], x: usize) -> usize {
            if p[x] != x {
                let r = find(p, p[x]);
                p[x] = r;
            }
            p[x]
        }
        let exonic = a.exonic_regions();
        for i in 0..exonic.len() {
            for j in 0..i {
                if exonic.seqid[i] == exonic.seqid[j]
                    && exonic.strand[i] == exonic.strand[j]
                    && exonic.start[i] < exonic.end[j]
                    && exonic.start[j] < exonic.end[i]
                {
                    let (x, y) = (
                        find(&mut parent, exonic.owner[i] as usize),
                        find(&mut parent, exonic.owner[j] as usize),
                    );
                    parent[x.max(y)] = x.min(y);
                }
            }
        }
        let mut groups: BTreeMap<usize, Vec<u32>> = BTreeMap::new();
        for gene in 0..n {
            groups
                .entry(find(&mut parent, gene))
                .or_default()
                .push(gene as u32);
        }
        let agg = a.exon_bins(true);
        let mut want = Vec::new();
        for members in groups.values() {
            for (k, (s, e, txs)) in naive_bins(&a, members).into_iter().enumerate() {
                want.push((members[0], k as u32 + 1, s, e, txs, members.clone()));
            }
        }
        let got: Vec<_> = (0..agg.set.len())
            .map(|i| {
                (
                    agg.set.owner[i],
                    agg.part[i],
                    agg.set.start[i],
                    agg.set.end[i],
                    agg.transcripts.get(i as u32).to_vec(),
                    agg.genes.get(i as u32).to_vec(),
                )
            })
            .collect();
        assert_eq!(got, want);
    }
}

#[test]
fn windows_follow_strand_and_clip() {
    let a = random(3);
    let g = a.genes();
    let lengths: Vec<Option<u32>> = (0..a.seqid_count())
        .map(|s| (s != 1).then_some(LEN - 100))
        .collect();
    let d = a.derive().clip(&lengths);
    let flank = Flank::new(250, 40);
    let promoters = d.promoters(Level::Gene, flank);
    let bodies = d.gene_bodies(Level::Gene, Flank::new(500, 0));
    let tes = d.tes_windows(Level::Gene, Flank::new(30, 70));
    let tss = a.derive().tss(Level::Gene);
    assert_eq!(promoters.len(), g.len());
    for i in 0..g.len() {
        let (s, e) = (g.start[i] as i64, g.end[i] as i64);
        let len = lengths[g.seqid[i] as usize].map_or(i64::MAX, |l| l as i64);
        let clip = |a: i64, b: i64| {
            let b = b.min(len).max(0);
            (a.max(0).min(b) as u32, b as u32)
        };
        let rev = g.strand[i] == Strand::Reverse;
        // 1-based GenomicRanges::promoters, then back to half-open.
        let (tss1, want_p) = if rev {
            (e, clip(e - 40, e + 250))
        } else {
            (s + 1, clip(s - 250, s + 40))
        };
        assert_eq!((promoters.start[i], promoters.end[i]), want_p, "gene {i}");
        assert_eq!(promoters.owner[i], i as u32);
        assert_eq!(a.str(promoters.name[i]), a.str(g.id[i]));
        assert_eq!((tss.start[i] as i64 + 1, tss.end[i] as i64), (tss1, tss1));
        let want_b = if rev {
            clip(s, e + 500)
        } else {
            clip(s - 500, e)
        };
        assert_eq!((bodies.start[i], bodies.end[i]), want_b);
        let want_t = if rev {
            clip(s - 70, s + 30)
        } else {
            clip(e - 30, e + 70)
        };
        assert_eq!((tes.start[i], tes.end[i]), want_t);
    }
    // Windows index like any other set.
    let mut hits = Vec::new();
    promoters.overlaps(
        promoters.seqid[0],
        promoters.start[0],
        promoters.end[0].max(promoters.start[0] + 1),
        &mut hits,
    );
    assert!(hits.contains(&0) || promoters.start[0] == promoters.end[0]);
    let t = a.derive().promoters(Level::Transcript, flank);
    assert_eq!(t.len(), a.transcripts().len());
    assert_eq!(t.owner_level(), Some(Level::Transcript));
}
