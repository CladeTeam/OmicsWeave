//! Explicit MEME 5.5.9 interoperability comparison.

use std::fmt::Write as _;
use std::process::Command;

use omicsweave_core::Strings;
use omicsweave_ffm::{Alphabet, FfmReader, FfmWriter, Motif, MotifRef, Schema, WriterOptions};
use omicsweave_motif::fimo::{FimoText, write_text};
use omicsweave_motif::{Background, Hits, Program, ScanOptions, Scoring, Strands};

use crate::rng::Rng;

struct Case {
    meme: String,
    motifs: Vec<Motif>,
    ids: Vec<String>,
    alts: Vec<String>,
    fasta: String,
    names: Vec<String>,
    data: Vec<u8>,
    offsets: Vec<i64>,
    bfile: String,
    bg: Vec<f64>,
}

fn case(seed: u64) -> Case {
    let mut rng = Rng(seed);
    let mut meme = String::from(
        "MEME version 4\n\nALPHABET= ACGT\n\nstrands: + -\n\nBackground letter frequencies\nA 0.25 C 0.25 G 0.25 T 0.25\n\n",
    );
    let mut motifs = Vec::new();
    let mut ids = Vec::new();
    let mut alts = Vec::new();
    for i in 0..12 {
        let width = 1 + rng.below(20);
        let mut rows = Vec::new();
        let mut text_rows = String::new();
        for _ in 0..width {
            let mut raw: Vec<f64> = (0..4)
                .map(|_| {
                    let x = rng.unit();
                    if rng.below(5) == 0 { 0.0 } else { x * x }
                })
                .collect();
            if raw.iter().all(|&v| v == 0.0) {
                raw[rng.below(4)] = 1.0;
            }
            let sum: f64 = raw.iter().sum();
            // One motif carries rows MEME has to renormalize.
            let scale = if i == 3 { 1.00003 } else { 1.0 };
            let row: Vec<f64> = raw
                .iter()
                .map(|v| {
                    format!("{:.6}", (v / sum * scale).min(1.0))
                        .parse()
                        .unwrap()
                })
                .collect();
            for v in &row {
                write!(text_rows, " {v:.6}").unwrap();
            }
            text_rows.push('\n');
            rows.push(row);
        }
        let sites = match i % 4 {
            0 => None,
            1 => Some(12.7),
            _ => Some((5 + rng.below(60)) as f64),
        };
        let id = format!("M{i}");
        let alt = format!("ALT{i}");
        write!(
            meme,
            "MOTIF {id} {alt}\nletter-probability matrix: alength= 4 w= {width}"
        )
        .unwrap();
        if let Some(s) = sites {
            write!(meme, " nsites= {s}").unwrap();
        }
        meme.push_str(" E= 0\n");
        meme.push_str(&text_rows);
        meme.push('\n');
        motifs.push(Motif::pfm(&rows, sites));
        ids.push(id);
        alts.push(alt);
    }
    let mut fasta = String::new();
    let mut names = Vec::new();
    let mut data = Vec::new();
    let mut offsets = vec![0i64];
    for i in 0..25 {
        let len = match i {
            0 => 0,
            1 => 2,
            _ => rng.below(2500),
        };
        let name = format!("seq{i}");
        writeln!(fasta, ">{name} description").unwrap();
        let start = data.len();
        for _ in 0..len {
            data.push(match rng.below(200) {
                0..=2 => b'N',
                3 => b'R',
                4 => b'y',
                5..=24 => b"acgt"[rng.below(4)],
                _ => b"ACGT"[rng.below(4)],
            });
        }
        for chunk in data[start..].chunks(60) {
            fasta.push_str(std::str::from_utf8(chunk).unwrap());
            fasta.push('\n');
        }
        names.push(name);
        offsets.push(data.len() as i64);
    }
    let bg: Vec<f64> = (0..4).map(|_| 0.05 + 0.4 * rng.unit()).collect();
    let mut bfile = String::from("# 0-order\n");
    for (s, v) in "ACGT".chars().zip(&bg) {
        writeln!(bfile, "{s} {v:.6}").unwrap();
    }
    Case {
        meme,
        motifs,
        ids,
        alts,
        fasta,
        names,
        data,
        offsets,
        bfile,
        bg,
    }
}

fn ours(c: &Case, background: &Background, strands: Strands, thresh: f64) -> String {
    let alphabet = Alphabet::dna();
    let mut w = FfmWriter::new(
        Vec::new(),
        alphabet.clone(),
        Schema::default(),
        WriterOptions::default(),
    )
    .unwrap();
    for (id, m) in c.ids.iter().zip(&c.motifs) {
        w.push(id, m, &[]).unwrap();
    }
    let reader = FfmReader::open(w.finish().unwrap()).unwrap();
    let refs: Vec<MotifRef> = (0..reader.len())
        .map(|i| reader.motif(i).unwrap())
        .collect();
    let views: Vec<_> = refs.iter().map(|r| r.view()).collect();
    let program = Program::compile(
        &views,
        &alphabet,
        background,
        &Scoring::fimo(),
        strands,
        None,
    )
    .unwrap();
    let seqs = Strings {
        offsets: &c.offsets,
        data: &c.data,
    };
    let hits = program
        .scan(seqs, Hits::pvalue(thresh), &ScanOptions::default())
        .unwrap();
    let mut name_data = Vec::new();
    let mut name_offsets = vec![0i64];
    for n in &c.names {
        name_data.extend_from_slice(n.as_bytes());
        name_offsets.push(name_data.len() as i64);
    }
    let ids: Vec<&str> = c.ids.iter().map(String::as_str).collect();
    let alts: Vec<&str> = c.alts.iter().map(String::as_str).collect();
    let mut out = Vec::new();
    write_text(
        &mut out,
        &hits,
        &FimoText {
            alphabet: &alphabet,
            motif_ids: &ids,
            motif_alt_ids: &alts,
            sequence_names: Strings {
                offsets: &name_offsets,
                data: &name_data,
            },
            sequences: seqs,
            origins: None,
        },
    )
    .unwrap();
    String::from_utf8(out).unwrap()
}

#[test]
#[ignore = "requires MEME 5.5.9 fimo and OMICSWEAVE_FIMO; run explicitly on a worker"]
fn matches_fimo_text_output() {
    let fimo = std::env::var_os("OMICSWEAVE_FIMO")
        .expect("set OMICSWEAVE_FIMO to the MEME 5.5.9 executable before running ignored interoperability tests");
    let dir = tempfile::tempdir().unwrap();
    let alphabet = Alphabet::dna();
    for seed in [3u64, 17, 101] {
        let c = case(seed);
        let meme = dir.path().join("motifs.meme");
        let fasta = dir.path().join("seqs.fa");
        let bfile = dir.path().join("bg.txt");
        std::fs::write(&meme, &c.meme).unwrap();
        std::fs::write(&fasta, &c.fasta).unwrap();
        std::fs::write(&bfile, &c.bfile).unwrap();
        let file_bg = Background::from_meme(&alphabet, &c.bfile).unwrap();
        assert_eq!(file_bg.kmer_frequencies(1).unwrap().len(), c.bg.len());
        let runs: Vec<(Vec<String>, Background, Strands, f64)> = vec![
            (
                vec!["--bfile".into(), bfile.display().to_string()],
                file_bg,
                Strands::Both,
                1e-3,
            ),
            (
                vec!["--bfile".into(), "--uniform--".into(), "--norc".into()],
                Background::uniform(&alphabet),
                Strands::Forward,
                1e-2,
            ),
            (
                vec!["--bfile".into(), "--nrdb--".into()],
                Background::nrdb(&alphabet).unwrap(),
                Strands::Both,
                1e-4,
            ),
        ];
        for (args, bg, strands, thresh) in runs {
            let out = Command::new(&fimo)
                .args(["--text", "--no-qvalue", "--verbosity", "1"])
                .args(&args)
                .args(["--thresh", &thresh.to_string()])
                .arg(&meme)
                .arg(&fasta)
                .output()
                .unwrap();
            assert!(
                out.status.success(),
                "{}",
                String::from_utf8_lossy(&out.stderr)
            );
            let theirs = String::from_utf8(out.stdout).unwrap();
            let mine = ours(&c, &bg, strands, thresh);
            let t: Vec<&str> = theirs.lines().collect();
            let m: Vec<&str> = mine.lines().collect();
            for (i, (a, b)) in t.iter().zip(&m).enumerate() {
                assert_eq!(a, b, "seed {seed} {args:?}: line {i}");
            }
            assert_eq!(t.len(), m.len(), "seed {seed} {args:?}: line count");
            assert!(t.len() > 1, "seed {seed} {args:?}: no hits to compare");
        }
    }
}
