//! Coverage against `samtools depth` (aligned bases, with and without
//! deletions, MAPQ and flag filters) and against per-base reference
//! counts for the other units and strand splits.

mod common;

use std::path::Path;

use common::{Read, Rng, Spec, bigwig_values, reverse_for, samtools_depth, synth, write_bam};
use omicsweave_quant::{
    Budget, Coverage, CoverageUnit, Filter, Product, Scale, StrandSplit, Strandedness, Task,
    quantify,
};

fn run(bam: &Path, out: &Path, cov: Coverage, threads: usize, memory: u64) {
    let mut task = Task::new(vec![bam.to_owned()], out);
    task.products.push(Product::Coverage(cov));
    quantify(task, Budget::new(threads, memory), None).unwrap();
}

/// Per-base counts of what `unit` makes of the reads passing `filter`.
fn naive(
    reads: &[Read],
    refs: &[(String, u32)],
    unit: CoverageUnit,
    filter: &Filter,
    strand: Option<(bool, u8)>,
) -> Vec<Vec<u32>> {
    let mut d: Vec<Vec<u32>> = refs.iter().map(|r| vec![0; r.1 as usize]).collect();
    for r in reads {
        if r.flag & filter.exclude != 0
            || r.flag & filter.require != filter.require
            || r.mapq < filter.min_mapq
            || (filter.unique && r.nh() > 1)
        {
            continue;
        }
        if let Some((rev, lib)) = strand
            && reverse_for(r, lib) != rev
        {
            continue;
        }
        let ivs: Vec<(u32, u32)> = match unit {
            CoverageUnit::Bases { deletions } => r.blocks(deletions),
            CoverageUnit::Span => vec![(r.pos, r.end())],
            CoverageUnit::FivePrime => {
                let p = if r.reverse() { r.end() - 1 } else { r.pos };
                vec![(p, p + 1)]
            }
            CoverageUnit::ThreePrime => {
                let p = if r.reverse() { r.pos } else { r.end() - 1 };
                vec![(p, p + 1)]
            }
            CoverageUnit::Fragment { max_len } => {
                if r.flag & 1 == 0 {
                    vec![(r.pos, r.end())]
                } else if r.leader() && (max_len == 0 || r.tlen.unsigned_abs() as u32 <= max_len) {
                    vec![(r.pos, r.pos + r.tlen.unsigned_abs() as u32)]
                } else {
                    vec![]
                }
            }
        };
        for (s, e) in ivs {
            for p in s..e {
                d[r.tid][p as usize] += 1;
            }
        }
    }
    d
}

fn assert_same(ours: &[Vec<f32>], want: &[Vec<u32>], what: &str) {
    for (t, (a, b)) in ours.iter().zip(want).enumerate() {
        for (p, (x, y)) in a.iter().zip(b).enumerate() {
            assert_eq!(*x, *y as f32, "{what}: reference {t} base {p}");
        }
    }
}

#[test]
fn aligned_bases_match_samtools_depth() {
    let dir = tempfile::tempdir().unwrap();
    for seed in 0..3u64 {
        let spec = Spec {
            paired: seed != 1,
            ..Spec::default()
        };
        let reads = synth(&mut Rng::new(seed), &spec);
        let bam = write_bam(dir.path(), &format!("s{seed}"), &spec.refs, &reads, seed);
        for (deletions, mapq) in [(false, 0u8), (true, 0), (false, 20)] {
            let mut extra = vec![];
            if deletions {
                extra.push("-J".to_string());
            }
            let q = mapq.to_string();
            if mapq > 0 {
                extra.push("-Q".into());
                extra.push(q.clone());
            }
            let extra: Vec<&str> = extra.iter().map(String::as_str).collect();
            let want = samtools_depth(&bam, &extra, &spec.refs);
            let out = dir.path().join(format!("o{seed}{deletions}{mapq}"));
            let cov = Coverage {
                unit: CoverageUnit::Bases { deletions },
                filter: Filter {
                    min_mapq: mapq,
                    ..Filter::new()
                },
                ..Coverage::default()
            };
            run(&bam, &out, cov.clone(), 3, 256 << 20);
            let ours = bigwig_values(&out.join("coverage.bw"), &spec.refs);
            assert_same(
                &ours,
                &want,
                &format!("seed {seed} -J {deletions} -Q {mapq}"),
            );
            // samtools depth agrees with the reference counts too.
            let mine = naive(&reads, &spec.refs, cov.unit, &cov.filter, None);
            assert_eq!(mine, want, "reference counts vs samtools, seed {seed}");
        }
    }
}

#[test]
fn units_and_strands_match_reference_counts() {
    let dir = tempfile::tempdir().unwrap();
    let spec = Spec::default();
    let reads = synth(&mut Rng::new(7), &spec);
    let bam = write_bam(dir.path(), "u", &spec.refs, &reads, 7);
    let units = [
        CoverageUnit::Span,
        CoverageUnit::FivePrime,
        CoverageUnit::ThreePrime,
        CoverageUnit::Fragment { max_len: 0 },
        CoverageUnit::Fragment { max_len: 400 },
    ];
    for (i, unit) in units.into_iter().enumerate() {
        for split in [None, Some(1u8), Some(2u8)] {
            let filter = Filter {
                unique: i % 2 == 0,
                min_mapq: 5,
                ..Filter::new()
            };
            let out = dir.path().join(format!("u{i}{split:?}"));
            let strands = match split {
                None => StrandSplit::Merged,
                Some(l) => StrandSplit::Split(Strandedness::from_code(l).unwrap()),
            };
            run(
                &bam,
                &out,
                Coverage {
                    unit,
                    strands,
                    filter,
                    ..Coverage::default()
                },
                4,
                128 << 20,
            );
            match split {
                None => {
                    let ours = bigwig_values(&out.join("coverage.bw"), &spec.refs);
                    let want = naive(&reads, &spec.refs, unit, &filter, None);
                    assert_same(&ours, &want, &format!("{unit:?}"));
                }
                Some(lib) => {
                    for (rev, file) in [
                        (false, "coverage.forward.bw"),
                        (true, "coverage.reverse.bw"),
                    ] {
                        let ours = bigwig_values(&out.join(file), &spec.refs);
                        let want = naive(&reads, &spec.refs, unit, &filter, Some((rev, lib)));
                        assert_same(&ours, &want, &format!("{unit:?} lib {lib} {file}"));
                    }
                }
            }
        }
    }
}

#[test]
fn scaling_and_spilling_do_not_change_runs() {
    let dir = tempfile::tempdir().unwrap();
    let spec = Spec {
        fragments: 6000,
        ..Spec::default()
    };
    let reads = synth(&mut Rng::new(11), &spec);
    let bam = write_bam(dir.path(), "x", &spec.refs, &reads, 11);
    let base = Coverage {
        sidecar: true,
        ..Coverage::default()
    };
    // One thread and a roomy budget, against many threads and a spool so
    // small the runs spill.
    run(&bam, &dir.path().join("a"), base.clone(), 1, 1 << 30);
    let mut task = Task::new(vec![bam.clone()], dir.path().join("b"));
    task.products.push(Product::Coverage(base.clone()));
    let mut p = omicsweave_quant::plan(task, Budget::new(8, 256 << 20)).unwrap();
    p.stages[0].1 = 16 << 10;
    omicsweave_quant::run(&p, None).unwrap();
    let summary = std::fs::read_to_string(dir.path().join("b/summary.json")).unwrap();
    assert!(
        summary
            .lines()
            .any(|l| l.contains("spilled_bytes") && !l.contains(": 0,")),
        "the small spool should spill"
    );
    let a = std::fs::read(dir.path().join("a/coverage.bw")).unwrap();
    let b = std::fs::read(dir.path().join("b/coverage.bw")).unwrap();
    assert_eq!(a, b);
    assert!(dir.path().join("a/coverage.bw.sidecar").exists());
    let summary = std::fs::read_to_string(dir.path().join("a/summary.json")).unwrap();
    assert!(summary.contains("\"counted\""));
    let counted = reads
        .iter()
        .filter(|r| r.flag & Filter::DEFAULT_EXCLUDE == 0)
        .count() as f64;
    run(
        &bam,
        &dir.path().join("c"),
        Coverage {
            scale: Scale::Cpm,
            ..base
        },
        2,
        1 << 30,
    );
    let raw = bigwig_values(&dir.path().join("a/coverage.bw"), &spec.refs);
    let cpm = bigwig_values(&dir.path().join("c/coverage.bw"), &spec.refs);
    for (x, y) in raw.iter().flatten().zip(cpm.iter().flatten()) {
        let want = (*x as f64 * 1e6 / counted) as f32;
        assert!((y - want).abs() <= want.abs() * 1e-6, "{y} vs {want}");
    }
}
