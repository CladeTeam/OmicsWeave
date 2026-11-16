use crate::support::scan_fixtures::{Rng, dinucleotide, markov, pwm, reader, sequences};
use crate::support::scan_oracle::Collect;
use omicsweave_core::Strings;
use omicsweave_ffm::{Alphabet, Motif, MotifRef};
use omicsweave_motif::{Background, Hits, Program, ScanOptions, Scoring, Strands};

#[test]
fn combined_parts_scan_like_separate_programs() {
    let mut rng = Rng(61);
    let alphabet = Alphabet::dna();
    let a_motifs: Vec<Motif> = (0..3).map(|i| pwm(&mut rng, 5 + i)).collect();
    let b_motifs = vec![markov(&mut rng, 7), dinucleotide(&mut rng, 6)];
    let compile = |motifs: &[Motif], bg: &Background, scoring: &Scoring, strands| {
        let r = reader(motifs);
        let refs: Vec<MotifRef> = (0..r.len()).map(|i| r.motif(i).unwrap()).collect();
        let views: Vec<_> = refs.iter().map(|m| m.view()).collect();
        Program::compile(&views, &alphabet, bg, scoring, strands, Some(2)).unwrap()
    };
    let a = compile(
        &a_motifs,
        &Background::nrdb(&alphabet).unwrap(),
        &Scoring::fimo(),
        Strands::Both,
    );
    let markov_bg = Background::from_kmer_frequencies(
        4,
        vec![
            (0..4).map(|_| 0.5 + rng.unit()).collect(),
            (0..16).map(|_| 0.5 + rng.unit()).collect(),
        ],
    )
    .unwrap();
    let b = compile(&b_motifs, &markov_bg, &Scoring::exact(), Strands::Forward);
    let combined = Program::combine([a.clone(), b.clone()]).unwrap();
    assert_eq!(combined.parts().len(), 2);
    assert_eq!(combined.parts()[1].motifs(), 3..5);
    let (data, offsets) = sequences(&mut rng, 8);
    let seqs = Strings {
        offsets: &offsets,
        data: &data,
    };
    let opt = ScanOptions {
        threads: Some(3),
        tile_len: 2000,
    };
    let all = combined.scan(seqs, Collect::default(), &opt).unwrap();
    let mut want = a.scan(seqs, Collect::default(), &opt).unwrap();
    want.extend(
        b.scan(seqs, Collect::default(), &opt)
            .unwrap()
            .into_iter()
            .map(|mut w| {
                w.0 += 3;
                w
            }),
    );
    want.sort();
    assert_eq!(all, want);

    // A reducer limited to one part sees only that part's motifs.
    let hits = combined
        .scan(seqs, Hits::pvalue(1e-2).in_part(1), &opt)
        .unwrap();
    let separate = b.scan(seqs, Hits::pvalue(1e-2), &opt).unwrap();
    assert_eq!(hits.len(), separate.len());
    assert!(hits.iter().all(|h| h.motif >= 3));
}

#[test]
fn incompatible_ambiguity_dictionaries_cannot_be_combined() {
    use omicsweave_motif::{Ambiguity, AmbiguityRule};
    let alphabet = Alphabet::dna();
    let reader = reader(&[Motif::pfm(&[vec![4.0, 3.0, 2.0, 1.0]], None)]);
    let motif = reader.motif(0).unwrap();
    let make = |members| {
        let mut scoring = Scoring::exact();
        scoring.ambiguity = Ambiguity {
            rule: AmbiguityRule::Mean,
            codes: vec![(b'X', members)],
            unknown: None,
        };
        Program::compile(
            &[motif.view()],
            &alphabet,
            &Background::uniform(&alphabet),
            &scoring,
            Strands::Forward,
            Some(1),
        )
        .unwrap()
    };
    let error = Program::combine([make(vec![0]), make(vec![1])]).unwrap_err();
    assert!(error.to_string().contains("share"), "{error}");
    assert!(Program::combine([]).is_err());
}
