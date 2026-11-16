//! MCAST 5.5.9 interoperability; OMICSWEAVE_MCAST_STATS=1 also compares p/E/q values
//! using 1e8 random letters per case on each side.
use crate::mcast_support::{Rng, program};
use omicsweave_core::Strings;
use omicsweave_ffm::Motif;
use omicsweave_motif::ScanOptions;
use omicsweave_motif::mcast::{
    self, ClusterOptions, McastText, OutputFilter, Statistics, write_text,
};
use std::collections::BTreeMap;
use std::fmt::Write as _;
use std::process::Command;

fn revcomp(s: &[u8]) -> Vec<u8> {
    s.iter()
        .rev()
        .map(|&c| match c {
            b'A' => b'T',
            b'C' => b'G',
            b'G' => b'C',
            _ => b'A',
        })
        .collect()
}

struct Case {
    meme: String,
    motifs: Vec<Motif>,
    fasta: String,
    names: Vec<String>,
    data: Vec<u8>,
    offsets: Vec<i64>,
}

fn case(seed: u64, max_gap: usize, long: usize) -> Case {
    let mut rng = Rng(seed);
    let mut meme = String::from(
        "MEME version 4\n\nALPHABET= ACGT\n\nstrands: + -\n\nBackground letter frequencies\nA 0.25 C 0.25 G 0.25 T 0.25\n\n",
    );
    let mut motifs = Vec::new();
    let mut consensus = Vec::new();
    let count = 1 + rng.below(3);
    for i in 0..count {
        let width = 6 + rng.below(9);
        let mut rows = Vec::new();
        let mut text_rows = String::new();
        let mut cons = Vec::new();
        for _ in 0..width {
            let top = rng.below(4);
            let mut raw: Vec<f64> = (0..4)
                .map(|j| {
                    if j == top {
                        3.0 + 6.0 * rng.unit()
                    } else if rng.below(4) == 0 {
                        0.0
                    } else {
                        rng.unit()
                    }
                })
                .collect();
            if rng.below(7) == 0 {
                raw[(top + 2) % 4] = raw[top];
            }
            let sum: f64 = raw.iter().sum();
            let row: Vec<f64> = raw
                .iter()
                .map(|v| format!("{:.6}", v / sum).parse().unwrap())
                .collect();
            for v in &row {
                write!(text_rows, " {v:.6}").unwrap();
            }
            text_rows.push('\n');
            cons.push(b"ACGT"[top]);
            rows.push(row);
        }
        write!(
            meme,
            "MOTIF M{i}\nletter-probability matrix: alength= 4 w= {width} nsites= 20 E= 0\n{text_rows}\n"
        )
        .unwrap();
        motifs.push(Motif::pfm(&rows, Some(20.0)));
        consensus.push(cons);
    }
    let mut fasta = String::new();
    let mut names = Vec::new();
    let mut data = Vec::new();
    let mut offsets = vec![0i64];
    let nseq = 6 + rng.below(4);
    for si in 0..nseq {
        let target = if si == 0 && long > 0 {
            long
        } else {
            300 + rng.below(3000)
        };
        let mut s: Vec<u8> = Vec::with_capacity(target + 200);
        let mut gaps = vec![0usize, 1, max_gap - 1, max_gap, max_gap + 1, 2, 5];
        while s.len() < target {
            // A cluster: several sites with chosen gaps.
            s.extend((0..20 + rng.below(400)).map(|_| b"ACGT"[rng.below(4)]));
            let sites = 1 + rng.below(4);
            for k in 0..sites {
                if k > 0 {
                    let g = if gaps.is_empty() || rng.below(3) == 0 {
                        rng.below(max_gap + 10)
                    } else {
                        gaps.swap_remove(rng.below(gaps.len()))
                    };
                    s.extend((0..g).map(|_| b"ACGT"[rng.below(4)]));
                }
                let m = rng.below(consensus.len());
                let mut site: Vec<u8> = consensus[m]
                    .iter()
                    .map(|&c| {
                        if rng.below(12) == 0 {
                            b"ACGT"[rng.below(4)]
                        } else {
                            c
                        }
                    })
                    .collect();
                if rng.below(2) == 0 {
                    site = revcomp(&site);
                }
                s.extend(site);
            }
        }
        // Lower case, ambiguity codes and N runs.
        for b in s.iter_mut() {
            match rng.below(400) {
                0 => *b = b'N',
                1 => *b = b"RYSWKMBDHV"[rng.below(10)],
                2..=40 => *b = b.to_ascii_lowercase(),
                _ => {}
            }
        }
        let name = format!("seq{si}");
        writeln!(fasta, ">{name}").unwrap();
        for chunk in s.chunks(60) {
            fasta.push_str(std::str::from_utf8(chunk).unwrap());
            fasta.push('\n');
        }
        data.extend(&s);
        offsets.push(data.len() as i64);
        names.push(name);
    }
    Case {
        meme,
        motifs,
        fasta,
        names,
        data,
        offsets,
    }
}

/// Rows of an MCAST table keyed by (sequence, start, stop).
fn rows(text: &str) -> BTreeMap<(String, u64, u64), Vec<String>> {
    text.lines()
        .filter(|l| l.starts_with("cluster-"))
        .map(|l| {
            let f: Vec<String> = l.split('\t').map(str::to_string).collect();
            (
                (f[1].clone(), f[2].parse().unwrap(), f[3].parse().unwrap()),
                f,
            )
        })
        .collect()
}

#[test]
#[ignore = "requires MEME 5.5.9 mcast and OMICSWEAVE_MCAST; run explicitly on a worker"]
fn matches_mcast() {
    let bin = std::env::var_os("OMICSWEAVE_MCAST")
        .expect("set OMICSWEAVE_MCAST to the MEME 5.5.9 executable before running ignored interoperability tests");
    let with_stats = std::env::var_os("OMICSWEAVE_MCAST_STATS").is_some();
    let dir = tempfile::tempdir().unwrap();
    let runs: [(u64, f64, u32, usize); 6] = [
        (1, 0.0005, 50, 0),
        (2, 0.0005, 30, 0),
        (3, 0.001, 50, 0),
        (4, 0.0002, 40, 0),
        (5, 0.0005, 50, 420_000),
        (6, 0.002, 20, 0),
    ];
    let mut compared = 0usize;
    let mut stats_equal = 0usize;
    let mut stats_total = 0usize;
    for (seed, pthresh, max_gap, long) in runs {
        let c = case(seed, max_gap as usize, long);
        let meme = dir.path().join(format!("m{seed}.meme"));
        let fasta = dir.path().join(format!("s{seed}.fa"));
        let out_dir = dir.path().join(format!("out{seed}"));
        std::fs::write(&meme, &c.meme).unwrap();
        std::fs::write(&fasta, &c.fasta).unwrap();
        std::fs::create_dir_all(&out_dir).unwrap();
        let out = Command::new(&bin)
            .args(["--oc"])
            .arg(&out_dir)
            .args(["--verbosity", "1", "--no-pgc", "--output-pthresh", "1"])
            .args(["--motif-pthresh", &pthresh.to_string()])
            .args(["--max-gap", &max_gap.to_string()])
            .arg(&meme)
            .arg(&fasta)
            .output()
            .unwrap();
        assert!(
            out.status.success(),
            "mcast failed: {}",
            String::from_utf8_lossy(&out.stderr)
        );
        let theirs = std::fs::read_to_string(out_dir.join("mcast.tsv")).unwrap();

        let p = program(&c.motifs);
        let options = ClusterOptions {
            pthresh,
            max_gap,
            ..ClusterOptions::default()
        };
        let stats = Statistics::mcast(0, max_gap);
        let seqs = Strings {
            offsets: &c.offsets,
            data: &c.data,
        };
        let report = mcast::search(
            &p,
            seqs,
            &options,
            with_stats.then_some(&stats),
            &ScanOptions::default(),
        )
        .unwrap();
        let mut name_data = Vec::new();
        let mut name_offsets = vec![0i64];
        for n in &c.names {
            name_data.extend_from_slice(n.as_bytes());
            name_offsets.push(name_data.len() as i64);
        }
        let mut ours = Vec::new();
        write_text(
            &mut ours,
            &report,
            &McastText {
                sequence_names: Strings {
                    offsets: &name_offsets,
                    data: &name_data,
                },
                sequences: seqs,
                origins: None,
                filter: OutputFilter {
                    evalue: f64::MAX,
                    pvalue: 1.0,
                    qvalue: 1.0,
                },
            },
        )
        .unwrap();
        let ours = String::from_utf8(ours).unwrap();
        let (a, b) = (rows(&ours), rows(&theirs));
        let ka: Vec<_> = a.keys().collect();
        let kb: Vec<_> = b.keys().collect();
        assert_eq!(ka, kb, "case {seed}: cluster coordinates differ");
        for (k, ra) in &a {
            let rb = &b[k];
            assert_eq!(ra[0], rb[0], "case {seed} {k:?}: cluster id");
            assert_eq!(ra[4], rb[4], "case {seed} {k:?}: score");
            assert_eq!(ra[8], rb[8], "case {seed} {k:?}: matched sequence");
            if with_stats {
                stats_total += 1;
                if ra[5..8] == rb[5..8] {
                    stats_equal += 1;
                } else {
                    eprintln!(
                        "case {seed} {k:?}: ours {:?} mcast {:?}",
                        &ra[5..8],
                        &rb[5..8]
                    );
                }
            }
        }
        compared += a.len();
    }
    eprintln!("compared {compared} clusters");
    if with_stats {
        eprintln!("p/E/q agree on {stats_equal} of {stats_total} clusters");
        assert_eq!(stats_equal, stats_total);
    }
    assert!(compared > 0);
}
