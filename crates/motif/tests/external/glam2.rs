use crate::glam2_support::*;
use omicsweave_ffm::Alphabet;
use omicsweave_motif::Strands;
use omicsweave_motif::profile::{
    DirichletMixture, GLAM2_HITS, ProfileScoring, glam2_dna_alphabet, write_glam2_header,
};
use std::process::Command;

struct Run {
    n: usize,
    strands: Strands,
    scoring: ProfileScoring,
}

fn run(args: &str, alphabet: &Alphabet, mixture: Option<&str>) -> Run {
    let mut scoring = ProfileScoring::glam2(alphabet).unwrap();
    let mut n = GLAM2_HITS;
    let mut strands = Strands::Forward;
    let words: Vec<&str> = args.split_whitespace().collect();
    let mut i = 0;
    while i < words.len() {
        match words[i] {
            "-2" => strands = Strands::Both,
            "-n" => n = words[i + 1].parse().unwrap(),
            "-D" => scoring.delete_pseudo = words[i + 1].parse().unwrap(),
            "-E" => scoring.no_delete_pseudo = words[i + 1].parse().unwrap(),
            "-I" => scoring.insert_pseudo = words[i + 1].parse().unwrap(),
            "-J" => scoring.no_insert_pseudo = words[i + 1].parse().unwrap(),
            "-d" => {
                scoring.prior =
                    DirichletMixture::parse_glam2(mixture.unwrap(), alphabet.len()).unwrap()
            }
            other => panic!("unexpected {other}"),
        }
        i += if words[i] == "-2" { 1 } else { 2 };
    }
    Run {
        n,
        strands,
        scoring,
    }
}

#[test]
#[ignore = "requires MEME 5.5.9 glam2scan and OMICSWEAVE_GLAM2SCAN; run explicitly on a worker"]
fn matches_glam2scan_text_output() {
    let binary = std::env::var_os("OMICSWEAVE_GLAM2SCAN")
        .expect("set OMICSWEAVE_GLAM2SCAN to the MEME 5.5.9 executable before running ignored interoperability tests");
    let binary = binary.into_string().unwrap();
    let dir = tempfile::tempdir().unwrap();
    let mixture =
        "Mixture= 0.7\nAlpha= 1.2 0.3 0.3 0.3 0.3\nMixture= 0.3\nAlpha= 9 5.0 1.0 2.0 1.0\n";
    let mix_path = dir.path().join("mix.txt");
    std::fs::write(&mix_path, mixture).unwrap();
    let mix_arg = mix_path.display().to_string();
    let dna = glam2_dna_alphabet();
    let protein = Alphabet::protein();
    let mut compared = 0;
    for seed in [5u64, 23, 71, 1009] {
        let mut rng = Rng(seed);
        let cases: Vec<(&str, &Alphabet, String, String, Vec<String>)> = vec![
            (
                "n",
                &dna,
                random_motif(&mut rng, b"acgtAC", 10, 14),
                random_fasta(&mut rng, b"acgtACGTNnu-", 45, 140),
                vec![
                    "-2".into(),
                    "-n 50".into(),
                    "-2 -n 300 -D 0.5 -E 3 -I 0.3 -J 2".into(),
                    format!("-2 -n 150 -d {mix_arg}"),
                    "-2 -n 0".into(),
                ],
            ),
            (
                "n",
                &dna,
                random_motif(&mut rng, b"acgt", 3, 40),
                random_fasta(&mut rng, b"acgtN", 30, 60),
                vec!["-2 -n 1".into(), "-2 -n 400".into()],
            ),
            (
                "p",
                &protein,
                random_motif(&mut rng, b"ACDEFGHIKLMNPQRSTVWYacd", 8, 25),
                random_fasta(&mut rng, b"ACDEFGHIKLMNPQRSTVWYXBZ*acdef", 30, 90),
                vec!["".into(), "-n 200".into(), "-n 40 -D 0.05 -I 0.1".into()],
            ),
        ];
        for (alph_arg, alphabet, motif_text, fasta_text, runs) in &cases {
            let motif_path = dir.path().join("motif.glam2");
            let fasta_path = dir.path().join("seqs.fa");
            std::fs::write(&motif_path, motif_text).unwrap();
            std::fs::write(&fasta_path, fasta_text).unwrap();
            let motif = parse_glam2_motif(motif_text, alphabet);
            let fasta = Fasta::parse(fasta_text);
            for args in runs {
                let r = run(args, alphabet, Some(mixture));
                let mut argv = vec!["-t".to_string()];
                argv.extend(args.split_whitespace().map(String::from));
                argv.push(alph_arg.to_string());
                argv.push(motif_path.display().to_string());
                argv.push(fasta_path.display().to_string());
                let out = Command::new(&binary).args(&argv).output().unwrap();
                assert!(
                    out.status.success(),
                    "glam2scan {argv:?} failed: {}",
                    String::from_utf8_lossy(&out.stderr)
                );
                let theirs = String::from_utf8(out.stdout).unwrap();
                let version = theirs
                    .lines()
                    .nth(1)
                    .and_then(|l| l.strip_prefix("Version "))
                    .unwrap();
                let p = program(
                    std::slice::from_ref(&motif),
                    alphabet,
                    &r.scoring,
                    r.strands,
                )
                .unwrap();
                let mut ours = Vec::new();
                write_glam2_header(&mut ours, version, &format!("{binary} {}", argv.join(" ")))
                    .unwrap();
                let mut ours = String::from_utf8(ours).unwrap();
                ours.push_str(&hits_text(&p, &fasta, r.n, None));
                assert_eq!(
                    ours, theirs,
                    "seed {seed}, glam2scan {argv:?}\nmotif:\n{motif_text}"
                );
                compared += 1;
            }
        }
    }
    eprintln!("{compared} glam2scan runs agree byte for byte");
}
