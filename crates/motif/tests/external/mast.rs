//! Explicit MEME 5.5.9 interoperability comparison.

use std::fmt::Write as _;
use std::process::Command;

use omicsweave_core::Strings;
use omicsweave_ffm::{Alphabet, FfmReader, FfmWriter, Motif, MotifRef, Schema, WriterOptions};
use omicsweave_motif::mast::{
    Mast, MastOptions, MastReport, MastText, StrandMode, compile, correlations, format_scores,
    write_hit_list,
};
use omicsweave_motif::{Background, ScanOptions};

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
    motif_bg: Vec<f64>,
    bfile: String,
    file_bg: Vec<f64>,
}

fn case(seed: u64, ambiguous: bool) -> Case {
    let mut rng = Rng(seed);
    let motif_bg: Vec<f64> = {
        let raw: Vec<f64> = (0..4).map(|_| 0.15 + 0.2 * rng.unit()).collect();
        let sum: f64 = raw.iter().sum();
        raw.iter()
            .map(|v| format!("{:.3}", v / sum).parse().unwrap())
            .collect()
    };
    let mut meme = String::from("MEME version 4\n\nALPHABET= ACGT\n\nstrands: + -\n\n");
    meme.push_str("Background letter frequencies\n");
    for (s, v) in "ACGT".chars().zip(&motif_bg) {
        write!(meme, "{s} {v:.3} ").unwrap();
    }
    meme.push_str("\n\n");
    let mut motifs = Vec::new();
    let mut ids = Vec::new();
    let mut alts = Vec::new();
    for i in 0..6 {
        let width = 4 + rng.below(12);
        let mut rows = Vec::new();
        let mut text = String::new();
        for _ in 0..width {
            let mut raw: Vec<f64> = (0..4)
                .map(|_| {
                    let x = rng.unit();
                    if rng.below(6) == 0 { 0.0 } else { x * x * x }
                })
                .collect();
            if raw.iter().all(|&v| v == 0.0) {
                raw[rng.below(4)] = 1.0;
            }
            let sum: f64 = raw.iter().sum();
            let row: Vec<f64> = raw
                .iter()
                .map(|v| format!("{:.6}", v / sum).parse().unwrap())
                .collect();
            for v in &row {
                write!(text, " {v:.6}").unwrap();
            }
            text.push('\n');
            rows.push(row);
        }
        let sites = match i % 3 {
            0 => None,
            1 => Some(17.6),
            _ => Some((5 + rng.below(40)) as f64),
        };
        let id = format!("M{i}");
        let alt = if i == 2 {
            String::new()
        } else {
            format!("alt{i}")
        };
        write!(
            meme,
            "MOTIF {id} {alt}\nletter-probability matrix: alength= 4 w= {width}"
        )
        .unwrap();
        if let Some(s) = sites {
            write!(meme, " nsites= {s}").unwrap();
        }
        meme.push_str(" E= 0\n");
        meme.push_str(&text);
        meme.push('\n');
        motifs.push(Motif::pfm(&rows, sites));
        ids.push(id);
        alts.push(alt);
    }
    let mut fasta = String::new();
    let mut names = Vec::new();
    let mut data = Vec::new();
    let mut offsets = vec![0i64];
    for i in 0..30 {
        let len = match i {
            0 => 3,
            _ => 20 + rng.below(800),
        };
        let name = format!("s{i}");
        writeln!(fasta, ">{name}").unwrap();
        let start = data.len();
        for _ in 0..len {
            let b = match rng.below(300) {
                0 if ambiguous => b'N',
                1 if ambiguous => b"RYKMSWBDHV"[rng.below(10)],
                2..=40 => b"acgt"[rng.below(4)],
                _ => b"ACGT"[rng.below(4)],
            };
            data.push(b);
        }
        // Plant a motif's consensus so some sequences carry real signal.
        if i % 3 == 1 && len > 40 {
            let m = &motifs[i % motifs.len()];
            let at = start + rng.below(len - m.width as usize);
            for (j, f) in m.factors.iter().enumerate() {
                let best = (0..4)
                    .max_by(|&a, &b| f.table[a].total_cmp(&f.table[b]))
                    .unwrap();
                data[at + j] = b"ACGT"[best];
            }
        }
        for chunk in data[start..].chunks(60) {
            fasta.push_str(std::str::from_utf8(chunk).unwrap());
            fasta.push('\n');
        }
        names.push(name);
        offsets.push(data.len() as i64);
    }
    let file_bg: Vec<f64> = (0..4).map(|_| 0.05 + 0.3 * rng.unit()).collect();
    let mut bfile = String::from("# 0-order\n");
    for (s, v) in "ACGT".chars().zip(&file_bg) {
        writeln!(bfile, "{s} {v:.6}").unwrap();
    }
    let file_bg = file_bg
        .iter()
        .map(|v| format!("{v:.6}").parse().unwrap())
        .collect();
    Case {
        meme,
        motifs,
        ids,
        alts,
        fasta,
        names,
        data,
        offsets,
        motif_bg,
        bfile,
        file_bg,
    }
}

fn ours(
    c: &Case,
    bg: &Background,
    null_given: Option<&[f64]>,
    mode: StrandMode,
) -> (MastReport, String, Vec<Vec<f64>>) {
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
    let program = compile(&views, &alphabet, bg, null_given, mode, None).unwrap();
    let seqs = Strings {
        offsets: &c.offsets,
        data: &c.data,
    };
    let report = program
        .scan(
            seqs,
            Mast::new(MastOptions::new(mode)),
            &ScanOptions::default(),
        )
        .unwrap();
    let mut name_data = Vec::new();
    let mut name_offsets = vec![0i64];
    for n in &c.names {
        name_data.extend_from_slice(n.as_bytes());
        name_offsets.push(name_data.len() as i64);
    }
    let ids: Vec<&str> = c.ids.iter().map(String::as_str).collect();
    let alts: Vec<&str> = c.alts.iter().map(String::as_str).collect();
    let numbers: Vec<u32> = (1..=ids.len() as u32).collect();
    let mut out = Vec::new();
    write_hit_list(
        &mut out,
        &report,
        &MastText {
            sequence_names: Strings {
                offsets: &name_offsets,
                data: &name_data,
            },
            motif_ids: &ids,
            motif_alt_ids: &alts,
            motif_numbers: &numbers,
            stranded: true,
        },
    )
    .unwrap();
    let corr = correlations(&program, &(0..program.len()).collect::<Vec<_>>());
    (report, String::from_utf8(out).unwrap(), corr)
}

/// (name, strand, combined p-value, E-value) of each `<score>` in mast.xml.
fn xml_scores(xml: &str) -> Vec<(String, String, String, String)> {
    let attr = |line: &str, key: &str| -> String {
        let k = format!("{key}=\"");
        let i = line.find(&k).expect("attribute") + k.len();
        line[i..].split('"').next().unwrap().to_string()
    };
    let mut out = Vec::new();
    let mut name = String::new();
    for line in xml.lines() {
        let t = line.trim();
        if t.starts_with("<sequence ") {
            name = attr(t, "name");
        } else if t.starts_with("<score ") {
            out.push((
                name.clone(),
                attr(t, "strand"),
                attr(t, "combined_pvalue"),
                attr(t, "evalue"),
            ));
        }
    }
    out
}

type Run<'a> = (Vec<String>, &'a Background, Option<&'a [f64]>, StrandMode);

#[test]
#[ignore = "requires MEME 5.5.9 mast and OMICSWEAVE_MAST; run explicitly on a worker"]
fn matches_mast() {
    let mast = std::env::var_os("OMICSWEAVE_MAST")
        .expect("set OMICSWEAVE_MAST to the MEME 5.5.9 executable before running ignored interoperability tests");
    let dir = tempfile::tempdir().unwrap();
    for (seed, ambiguous) in [(5u64, false), (23, false), (41, true), (77, true)] {
        let c = case(seed, ambiguous);
        let meme = dir.path().join("motifs.meme");
        let fasta = dir.path().join("seqs.fa");
        let bfile = dir.path().join("bg.txt");
        std::fs::write(&meme, &c.meme).unwrap();
        std::fs::write(&fasta, &c.fasta).unwrap();
        std::fs::write(&bfile, &c.bfile).unwrap();
        let motif_bg = Background::frequencies(&c.motif_bg).unwrap();
        let file_bg = Background::frequencies(&c.file_bg).unwrap();
        let runs: Vec<Run<'_>> = vec![
            (vec![], &motif_bg, None, StrandMode::Combine),
            (vec!["-sep".into()], &motif_bg, None, StrandMode::Separate),
            (vec!["-norc".into()], &motif_bg, None, StrandMode::Forward),
            (
                vec!["-bfile".into(), bfile.display().to_string()],
                &file_bg,
                Some(&c.file_bg),
                StrandMode::Combine,
            ),
        ];
        for (args, bg, null_given, mode) in runs {
            let label = format!("seed {seed} {args:?}");
            let (report, hit_list, corr) = ours(&c, bg, null_given, mode);

            let out = Command::new(&mast)
                .args(["-nostatus", "-hit_list"])
                .args(&args)
                .arg(&meme)
                .arg(&fasta)
                .output()
                .unwrap();
            assert!(
                out.status.success(),
                "{label}: {}",
                String::from_utf8_lossy(&out.stderr)
            );
            let theirs: Vec<String> = String::from_utf8(out.stdout)
                .unwrap()
                .lines()
                .filter(|l| !l.starts_with('#'))
                .map(str::to_string)
                .collect();
            let mine: Vec<&str> = hit_list.lines().collect();
            for (i, (a, b)) in theirs.iter().zip(&mine).enumerate() {
                assert_eq!(a, b, "{label}: hit list line {i}");
            }
            assert_eq!(theirs.len(), mine.len(), "{label}: hit list lines");

            let outdir = dir.path().join("out");
            let out = Command::new(&mast)
                .args(["-nostatus", "-oc"])
                .arg(&outdir)
                .args(&args)
                .arg(&meme)
                .arg(&fasta)
                .output()
                .unwrap();
            assert!(
                out.status.success(),
                "{label}: {}",
                String::from_utf8_lossy(&out.stderr)
            );
            let xml = std::fs::read_to_string(outdir.join("mast.xml")).unwrap();
            for line in xml
                .lines()
                .filter(|l| l.trim().starts_with("<correlation "))
            {
                let get = |k: &str| {
                    let key = format!("{k}=\"");
                    let i = line.find(&key).unwrap() + key.len();
                    line[i..].split('"').next().unwrap().to_string()
                };
                let a: usize = get("idx_a").parse().unwrap();
                let b: usize = get("idx_b").parse().unwrap();
                assert_eq!(
                    get("value"),
                    format!("{:.2}", corr[b][a]),
                    "{label}: correlation"
                );
            }
            let theirs = xml_scores(&xml);
            let mine: Vec<(String, String, String, String)> = report
                .ranked
                .iter()
                .map(|&i| {
                    let s = &report.sequences[i];
                    let (p, e) = format_scores(s);
                    let strand = match s.strand {
                        Some(omicsweave_core::Strand::Reverse) => "reverse",
                        Some(_) => "forward",
                        None if mode == StrandMode::Forward => "forward",
                        None => "both",
                    };
                    (
                        c.names[s.sequence as usize].clone(),
                        strand.to_string(),
                        p,
                        e,
                    )
                })
                .collect();
            if mode == StrandMode::Separate {
                // mast.xml lists each sequence once, at its better strand,
                // with both strands' scores; compare as sets.
                let mut a = theirs.clone();
                let mut b = mine.clone();
                a.sort();
                b.sort();
                assert_eq!(a, b, "{label}: scores");
            } else {
                assert_eq!(theirs, mine, "{label}: scores");
            }
            assert!(!theirs.is_empty(), "{label}: nothing ranked");
        }
    }
}
