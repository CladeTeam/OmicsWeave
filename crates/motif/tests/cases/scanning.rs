use crate::support::scan_fixtures::{Rng, motif_zoo, pwm};
use crate::support::scan_oracle::check;
use omicsweave_ffm::{Alphabet, Motif};
use omicsweave_motif::{Background, Scoring};

#[test]
fn zero_order_background_exact_scoring() {
    let mut rng = Rng(7);
    let motifs = motif_zoo(&mut rng);
    let bg = Background::frequencies(&[0.3, 0.2, 0.25, 0.25]).unwrap();
    check(&motifs, &bg, &Scoring::exact(), 11);
}

#[test]
fn ambiguity_codes_score_by_mean_and_sum() {
    let mut rng = Rng(71);
    let motifs = motif_zoo(&mut rng);
    let alphabet = Alphabet::dna();
    let bg = Background::frequencies(&[0.3, 0.2, 0.2, 0.3]).unwrap();
    for rule in [
        omicsweave_motif::AmbiguityRule::Mean,
        omicsweave_motif::AmbiguityRule::Sum,
    ] {
        let mut scoring = Scoring::exact();
        scoring.ambiguity = omicsweave_motif::Ambiguity::iupac(&alphabet, rule).unwrap();
        check(&motifs, &bg, &scoring, 73);
        // Every other unknown byte read as N.
        scoring.ambiguity = scoring.ambiguity.with_unknown(b'N').unwrap();
        check(&motifs, &bg, &scoring, 79);
    }
}

#[test]
fn fimo_scoring_of_pwms() {
    let mut rng = Rng(19);
    let motifs: Vec<Motif> = (0..6).map(|i| pwm(&mut rng, 3 + i * 3)).collect();
    let bg = Background::nrdb(&Alphabet::dna()).unwrap();
    check(&motifs, &bg, &Scoring::fimo(), 23);
}

#[test]
fn markov_backgrounds() {
    let mut rng = Rng(31);
    let motifs = motif_zoo(&mut rng);
    for order in [1usize, 2] {
        let kmers = (1..=order + 1)
            .map(|k| {
                (0..4usize.pow(k as u32))
                    .map(|_| 0.5 + rng.unit())
                    .collect()
            })
            .collect();
        let bg = Background::from_kmer_frequencies(4, kmers).unwrap();
        check(&motifs, &bg, &Scoring::exact(), 37 + order as u64);
    }
}

/// Sums a real score model over a footprint whose position `lo` is index 0.

#[test]
fn malformed_sequence_views_fail_before_scanning_and_program_remains_usable() {
    use crate::support::resident::resident;
    use omicsweave_core::Strings;
    use omicsweave_motif::{Hits, Program, ScanOptions, Strands};
    let alphabet = Alphabet::dna();
    let reader = resident(&Motif::pfm(&[vec![4.0, 3.0, 2.0, 1.0]], None));
    let motif = reader.motif(0).unwrap();
    let program = Program::compile(
        &[motif.view()],
        &alphabet,
        &Background::uniform(&alphabet),
        &Scoring::exact(),
        Strands::Forward,
        Some(1),
    )
    .unwrap();
    for offsets in [&[][..], &[-1, 1], &[0, 2, 1], &[0, 5]] {
        let error = program
            .scan(
                Strings {
                    offsets,
                    data: b"ACGT",
                },
                Hits::new(omicsweave_motif::Threshold::Score(f64::NEG_INFINITY)),
                &ScanOptions::default(),
            )
            .unwrap_err();
        assert!(
            error.to_string().contains("motif scan sequences"),
            "offsets={offsets:?}: {error}"
        );
    }
    let hits = program
        .scan_sequence(
            b"ACGT",
            Hits::new(omicsweave_motif::Threshold::Score(f64::NEG_INFINITY)),
            &ScanOptions {
                threads: Some(0),
                tile_len: 0,
            },
        )
        .unwrap();
    assert_eq!(hits.len(), 4);
    assert_eq!(
        hits.iter().map(|h| h.start).collect::<Vec<_>>(),
        vec![0, 1, 2, 3]
    );
    assert!(
        program
            .scan_sequence(b"", Hits::pvalue(1.0), &ScanOptions::default())
            .unwrap()
            .is_empty()
    );
}
