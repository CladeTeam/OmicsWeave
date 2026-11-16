use crate::support::mcast::{Rng, program};
use omicsweave_core::Strings;
use omicsweave_ffm::Motif;
use omicsweave_motif::mcast::{self, ClusterOptions, cluster_threshold};
use omicsweave_motif::{Hits, LogScale, Program, ScanOptions};
use std::collections::BTreeMap;

fn sharp(consensus: &[u8]) -> Motif {
    let rows: Vec<Vec<f64>> = consensus
        .iter()
        .map(|&c| {
            let k = b"ACGT".iter().position(|&x| x == c).unwrap();
            (0..4).map(|j| if j == k { 0.91 } else { 0.03 }).collect()
        })
        .collect();
    Motif::pfm(&rows, Some(20.0))
}

fn filler(rng: &mut Rng, n: usize) -> Vec<u8> {
    // Low-complexity filler without strong matches to the test motifs.
    (0..n).map(|_| b"ACGT"[rng.below(4)]).collect()
}

const SITE: &[u8] = b"TTGACGCATG";

/// Copies of `SITE` separated by `gaps` in poly-A, which holds no hit.
fn planted(gaps: &[usize]) -> Vec<u8> {
    let mut s = vec![b'A'; 200];
    s.extend_from_slice(SITE);
    for &g in gaps {
        s.extend(std::iter::repeat_n(b'A', g));
        s.extend_from_slice(SITE);
    }
    s.extend(std::iter::repeat_n(b'A', 200));
    s
}

fn site_scores(p: &Program, seq: &[u8], pthresh: f64) -> BTreeMap<u32, f64> {
    let hits = p
        .scan_sequence(seq, Hits::pvalue(pthresh), &ScanOptions::default())
        .unwrap();
    let log_t = LogScale::MemeBits.log(pthresh);
    let mut best = BTreeMap::new();
    for h in hits {
        let s = log_t - LogScale::MemeBits.log(h.pvalue);
        let e = best.entry(h.start).or_insert(s);
        if s > *e {
            *e = s;
        }
    }
    best
}

#[test]
fn cluster_score_is_hits_minus_gaps() {
    let p = program(&[sharp(SITE)]);
    let opts = ClusterOptions::default();
    let d = cluster_threshold(opts.pthresh, opts.max_gap, 2, 1.0);
    let ext = d / opts.max_gap as f64;
    let seq = planted(&[7, 23]);
    let sites = site_scores(&p, &seq, opts.pthresh);
    let found = mcast::search(
        &p,
        Strings {
            offsets: &[0, seq.len() as i64],
            data: &seq,
        },
        &opts,
        None,
        &ScanOptions::default(),
    )
    .unwrap();
    let c = found
        .clusters
        .iter()
        .find(|c| c.cluster.start == 200)
        .expect("planted cluster");
    let expected = sites[&200] + sites[&217] + sites[&250] - (7.0 + 23.0) * ext;
    assert_eq!(c.cluster.end, 260);
    assert_eq!(c.cluster.hits, 3);
    assert!(
        (c.cluster.score - expected).abs() < 1e-9,
        "{} vs {expected}",
        c.cluster.score
    );
}

#[test]
fn gaps_of_max_gap_split_and_adjacent_sites_join() {
    let p = program(&[sharp(SITE)]);
    let opts = ClusterOptions::default();
    let d = cluster_threshold(opts.pthresh, opts.max_gap, 2, 1.0);
    let scan = ScanOptions::default();
    let run = |seq: &[u8]| {
        mcast::search(
            &p,
            Strings {
                offsets: &[0, seq.len() as i64],
                data: seq,
            },
            &opts,
            None,
            &scan,
        )
        .unwrap()
        .clusters
    };
    let split = planted(&[opts.max_gap as usize]);
    let found = run(&split);
    let at: Vec<u32> = found.iter().map(|c| c.cluster.start).collect();
    assert!(at.contains(&200) && at.contains(&260), "{at:?}");
    let joined = planted(&[opts.max_gap as usize - 1]);
    let found = run(&joined);
    assert!(
        found
            .iter()
            .any(|c| c.cluster.start == 200 && c.cluster.end == 269)
    );
    let adjacent = planted(&[0]);
    let sites = site_scores(&p, &adjacent, opts.pthresh);
    let found = run(&adjacent);
    let c = found.iter().find(|c| c.cluster.start == 200).unwrap();
    assert_eq!(c.cluster.end, 220);
    let expected = sites[&200] + sites[&210] - d;
    assert!((c.cluster.score - expected).abs() < 1e-9);
}

#[test]
fn results_do_not_depend_on_tiles() {
    let p = program(&[sharp(SITE), sharp(b"GGGCATTACCA")]);
    let mut rng = Rng(77);
    let mut data = Vec::new();
    let mut offsets = vec![0i64];
    for _ in 0..5 {
        let len = 3000 + rng.below(5000);
        let mut s = filler(&mut rng, len);
        for _ in 0..12 {
            let at = rng.below(s.len() - 20);
            s[at..at + 10].copy_from_slice(SITE);
        }
        data.extend(s);
        offsets.push(data.len() as i64);
    }
    let seqs = Strings {
        offsets: &offsets,
        data: &data,
    };
    let opts = ClusterOptions {
        segment: Some(4000),
        ..ClusterOptions::default()
    };
    let a = mcast::search(&p, seqs, &opts, None, &ScanOptions::default()).unwrap();
    let small = ScanOptions {
        tile_len: 1024,
        ..ScanOptions::default()
    };
    let b = mcast::search(&p, seqs, &opts, None, &small).unwrap();
    assert!(!a.clusters.is_empty());
    let key =
        |r: &mcast::Report| -> Vec<_> { r.clusters.iter().map(|c| (c.id, c.cluster)).collect() };
    assert_eq!(key(&a), key(&b));
}

#[test]
fn invalid_cluster_parameters_fail_and_the_program_can_be_reused() {
    let p = program(&[sharp(SITE)]);
    let seq = planted(&[7, 23]);
    let sequences = Strings {
        offsets: &[0, seq.len() as i64],
        data: &seq,
    };
    let defaults = ClusterOptions::default();
    for options in [
        ClusterOptions {
            pthresh: 0.0,
            ..defaults.clone()
        },
        ClusterOptions {
            pthresh: f64::NAN,
            ..defaults.clone()
        },
        ClusterOptions {
            pthresh: 1.1,
            ..defaults.clone()
        },
        ClusterOptions {
            max_gap: 0,
            ..defaults.clone()
        },
        ClusterOptions {
            egcost: f64::INFINITY,
            ..defaults.clone()
        },
        ClusterOptions {
            threshold: Some(-1.0),
            ..defaults.clone()
        },
        ClusterOptions {
            segment: Some(defaults.overlap),
            ..defaults.clone()
        },
    ] {
        let error =
            mcast::search(&p, sequences, &options, None, &ScanOptions::default()).unwrap_err();
        assert!(error.to_string().contains("MCAST options"), "{error}");
    }
    assert!(
        !mcast::search(&p, sequences, &defaults, None, &ScanOptions::default())
            .unwrap()
            .clusters
            .is_empty()
    );
}
