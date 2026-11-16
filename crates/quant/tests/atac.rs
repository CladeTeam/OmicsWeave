//! Single-cell ATAC against the 10x / ArchR definitions computed on the
//! generated reads: fragments (Tn5 shift, duplicates collapsed per
//! barcode), cut-site matrices over peaks and tiles, gene activity, and the
//! tabix index htslib builds for the same file.

mod common;
mod support;

use std::collections::{BTreeMap, BTreeSet};
use std::sync::Arc;

use common::{Read, Rng, Spec, big_bed, genes, gtf, gunzip, read_10x, synth, write_bam};
use omicsweave_gtf::{Annotation, ReadOptions};
use omicsweave_quant::{Atac, Budget, Filter, Product, Task, quantify};

fn passes(r: &Read, f: &Filter) -> bool {
    r.flag & f.exclude == 0 && r.flag & f.require == f.require && r.mapq >= f.min_mapq
}

/// A fragment: reference, start, end, barcode.
type FragmentKey = (usize, u32, u32, String);

/// Fragment → read pairs.
fn naive_fragments(
    reads: &[Read],
    refs: &[(String, u32)],
    f: &Filter,
) -> BTreeMap<FragmentKey, u32> {
    let mut out = BTreeMap::new();
    for r in reads {
        if !passes(r, f) || !r.leader() {
            continue;
        }
        let Some(cb) = r.tag("CB") else { continue };
        let s = r.pos + 4;
        let e = r.pos + r.tlen.unsigned_abs() as u32;
        if e < 5 || e - 5 <= s {
            continue;
        }
        let _ = refs;
        *out.entry((r.tid, s, e - 5, cb.to_owned())).or_insert(0) += 1;
    }
    out
}

#[test]
fn fragments_matrices_and_metrics_match_definitions() {
    let dir = tempfile::tempdir().unwrap();
    let spec = Spec {
        fragments: 6000,
        repeats: 25,
        ..Spec::default()
    };
    let mut rng = Rng::new(21);
    let reads = synth(&mut rng, &spec);
    let gs = genes(&mut rng, &spec.refs, 10);
    let bam = write_bam(dir.path(), "atac", &spec.refs, &reads, 21);
    let gtf_path = dir.path().join("g.gtf");
    std::fs::write(&gtf_path, gtf(&spec.refs, &gs)).unwrap();
    let ann = Arc::new(Annotation::read(&gtf_path, &ReadOptions::default()).unwrap());
    // Peaks: sorted, some adjacent, none overlapping.
    let mut peaks = Vec::new();
    for (t, (_, len)) in spec.refs.iter().enumerate() {
        let mut at = rng.range(0, 500) as u32;
        while at + 600 < *len {
            let w = rng.range(100, 500) as u32;
            peaks.push((t, at, at + w));
            at += w + rng.range(0, 3000) as u32;
        }
    }
    let bb = dir.path().join("peaks.bb");
    big_bed(&bb, &spec.refs, &peaks);

    let atac = Atac {
        peaks: Some(bb),
        tiles: Some(2000),
        ..Atac::default()
    };
    let out = dir.path().join("out");
    let mut task = Task::new(vec![bam.clone()], &out);
    task.annotation = Some(ann);
    task.products.push(Product::Atac(atac.clone()));
    quantify(task.clone(), Budget::new(4, 512 << 20), None).unwrap();

    // Fragments: every line, in order.
    let want = naive_fragments(&reads, &spec.refs, &atac.filter);
    let text = gunzip(&out.join("atac/fragments.tsv.gz"));
    let mut got = BTreeMap::new();
    let mut last = (0usize, 0u32);
    for l in text.lines() {
        let f: Vec<&str> = l.split('\t').collect();
        let t = spec.refs.iter().position(|r| r.0 == f[0]).unwrap();
        let s: u32 = f[1].parse().unwrap();
        assert!((t, s) >= last, "fragments out of order");
        last = (t, s);
        let prev = got.insert(
            (t, s, f[2].parse().unwrap(), f[3].to_owned()),
            f[4].parse().unwrap(),
        );
        assert!(prev.is_none(), "duplicate fragment line {l}");
    }
    assert_eq!(got, want);
    assert!(want.values().any(|&c| c > 1), "the data has duplicates");

    // Barcodes in order of first fragment.
    let mut order = Vec::new();
    let mut seen = BTreeSet::new();
    for r in &reads {
        if passes(r, &atac.filter) && r.leader() {
            let e = r.pos + r.tlen.unsigned_abs() as u32;
            if e >= 5
                && e - 5 > r.pos + 4
                && let Some(cb) = r.tag("CB")
                && seen.insert(cb.to_owned())
            {
                order.push(cb.to_owned());
            }
        }
    }
    for m in ["peaks", "tiles", "gene_activity"] {
        let b: Vec<String> = gunzip(&out.join(format!("atac/{m}/barcodes.tsv.gz")))
            .lines()
            .map(str::to_owned)
            .collect();
        assert_eq!(b, order, "{m} barcodes");
    }

    // Cut sites per (barcode, feature), each distinct fragment once.
    let site_matrix = |feature: &dyn Fn(usize, u32) -> Option<String>| {
        let mut m: BTreeMap<(String, String), u32> = BTreeMap::new();
        for (t, s, e, cb) in want.keys() {
            for p in [*s, *e - 1] {
                if let Some(id) = feature(*t, p) {
                    *m.entry((cb.clone(), id)).or_insert(0) += 1;
                }
            }
        }
        m
    };
    let peak_of = |t: usize, p: u32| {
        peaks
            .iter()
            .find(|&&(pt, s, e)| pt == t && s <= p && p < e)
            .map(|&(_, s, e)| format!("{}:{s}-{e}", spec.refs[t].0))
    };
    assert_eq!(read_10x(&out.join("atac/peaks")), site_matrix(&peak_of));
    let tile_of = |t: usize, p: u32| {
        let len = spec.refs[t].1;
        let b = p / 2000;
        Some(format!(
            "{}:{}-{}",
            spec.refs[t].0,
            b * 2000,
            ((b + 1) * 2000).min(len)
        ))
    };
    assert_eq!(read_10x(&out.join("atac/tiles")), site_matrix(&tile_of));

    // Gene activity: distinct fragments overlapping the body + 2 kb
    // upstream of each gene.
    let mut ga: BTreeMap<(String, String), u32> = BTreeMap::new();
    for (t, s, e, cb) in want.keys() {
        for g in &gs {
            if g.tid != *t {
                continue;
            }
            let (ws, we) = if g.strand == '-' {
                (g.start, (g.end + 2000).min(spec.refs[*t].1))
            } else {
                (g.start.saturating_sub(2000), g.end)
            };
            if *s < we && ws < *e {
                *ga.entry((cb.clone(), g.id.clone())).or_insert(0) += 1;
            }
        }
    }
    assert_eq!(read_10x(&out.join("atac/gene_activity")), ga);

    // Metrics: read pairs, fragments, and peak fragments per barcode.
    let metrics = std::fs::read_to_string(out.join("atac/barcode_metrics.csv")).unwrap();
    let mut lines = metrics.lines();
    assert_eq!(
        lines.next().unwrap(),
        "barcode,read_pairs,fragments,peak_fragments,peak_cuts,tss_cuts,tss_flank_cuts"
    );
    for l in lines {
        let f: Vec<&str> = l.split(',').collect();
        let mine: Vec<(&FragmentKey, &u32)> = want.iter().filter(|(k, _)| k.3 == f[0]).collect();
        let pairs: u32 = mine.iter().map(|(_, c)| **c).sum();
        let in_peak = mine
            .iter()
            .filter(|(k, _)| peak_of(k.0, k.1).is_some() || peak_of(k.0, k.2 - 1).is_some())
            .count();
        assert_eq!(f[1].parse::<u32>().unwrap(), pairs, "{}", f[0]);
        assert_eq!(f[2].parse::<usize>().unwrap(), mine.len(), "{}", f[0]);
        assert_eq!(f[3].parse::<usize>().unwrap(), in_peak, "{}", f[0]);
    }

    support::assert_execution_invariant(
        dir.path(),
        |target| {
            let mut t = task.clone();
            t.out_dir = target.to_owned();
            t
        },
        8 << 10,
        1,
    );

    support::assert_storage_invariant(
        &out,
        |target| {
            let mut t = task.clone();
            t.out_dir = target.to_owned();
            t
        },
        &[
            ("atac/peaks", "matrix"),
            ("atac/tiles", "matrix"),
            ("atac/gene_activity", "matrix"),
        ],
    );

    // The tabix index is the one htslib builds from the same file.
    common::assert_tabix_as_htslib(&out.join("atac/fragments.tsv.gz"), dir.path());
}
