use crate::support::resident::resident;
use omicsweave_ffm::{Alphabet, Factor, FactorKind, Motif, Support};
use omicsweave_motif::{
    Ambiguity, AmbiguityRule, Background, ColumnPrior, Hits, Program, Quantizer, ScanOptions,
    Scoring, Strands, Threshold,
};
use std::sync::Arc;

#[test]
fn ambiguous_scores_above_the_ordinary_maximum_still_pass_thresholds() {
    let reader = resident(&Motif::pfm(&[vec![20.0, 20.0, 1.0, 1.0]], None));
    let motif = reader.motif(0).unwrap();
    let mut scoring = Scoring::exact();
    scoring.quantizer = Quantizer::Precision { states: 128 };
    scoring.ambiguity = Ambiguity {
        rule: AmbiguityRule::Sum,
        codes: vec![(b'X', vec![0, 1])],
        unknown: None,
    };
    let program = Program::compile(
        &[motif.view()],
        &Alphabet::dna(),
        &Background::uniform(&Alphabet::dna()),
        &scoring,
        Strands::Forward,
        Some(1),
    )
    .unwrap();
    let max = program.motif(0).strands()[0]
        .full()
        .int_model()
        .raw(program.motif(0).strands()[0].full().int_model().max);
    let above = program
        .scan_sequence(
            b"X",
            Hits::new(Threshold::Score(max + 0.1)),
            &ScanOptions::default(),
        )
        .unwrap();
    assert_eq!(above.len(), 1);
    assert!(above[0].score > max);
}

#[test]
fn sparse_motifs_with_unrepresentable_width_are_rejected_before_expansion() {
    let reader = resident(&Motif {
        width: u32::MAX,
        support: Support::All,
        sites: None,
        factors: vec![Factor {
            kind: FactorKind::Potential,
            positions: vec![0],
            table: vec![1.0; 4],
        }],
        gaps: None,
    });
    let motif = reader.motif(0).unwrap();
    assert!(
        Program::compile(
            &[motif.view()],
            &Alphabet::dna(),
            &Background::uniform(&Alphabet::dna()),
            &Scoring::exact(),
            Strands::Forward,
            Some(1),
        )
        .is_err()
    );
}

#[test]
fn reduced_support_factors_reject_an_oversized_global_expansion() {
    let reader = resident(&Motif {
        width: 13,
        support: Support::Codes(vec![0]),
        sites: None,
        factors: vec![Factor {
            kind: FactorKind::Joint,
            positions: (0..13).collect(),
            table: vec![1.0; 1 << 13],
        }],
        gaps: None,
    });
    let motif = reader.motif(0).unwrap();
    let result = Program::compile(
        &[motif.view()],
        &Alphabet::dna(),
        &Background::uniform(&Alphabet::dna()),
        &Scoring::exact(),
        Strands::Forward,
        Some(1),
    );
    assert!(
        result
            .unwrap_err()
            .to_string()
            .contains("expanded factor table")
    );
}

#[derive(Debug)]
struct FixedPrior(Vec<f64>);

impl ColumnPrior for FixedPrior {
    fn smooth(&self, _column: &[f64], _sites: f64, _letters: &[f64]) -> Vec<f64> {
        self.0.clone()
    }
}

#[test]
fn invalid_column_prior_outputs_are_compile_errors() {
    let reader = resident(&Motif::pfm(&[vec![4.0, 3.0, 2.0, 1.0]], None));
    let motif = reader.motif(0).unwrap();
    for values in [
        vec![],
        vec![0.25; 3],
        vec![f64::NAN; 4],
        vec![-1.0; 4],
        vec![0.0; 4],
    ] {
        let mut scoring = Scoring::exact();
        scoring.column_prior = Some(Arc::new(FixedPrior(values)));
        let result = Program::compile(
            &[motif.view()],
            &Alphabet::dna(),
            &Background::uniform(&Alphabet::dna()),
            &scoring,
            Strands::Forward,
            Some(1),
        );
        assert!(result.unwrap_err().to_string().contains("column prior"));
    }
}

#[test]
fn ambiguity_bytes_must_identify_one_code() {
    let reader = resident(&Motif::pfm(&[vec![4.0, 3.0, 2.0, 1.0]], None));
    let motif = reader.motif(0).unwrap();
    let mut scoring = Scoring::exact();
    scoring.ambiguity = Ambiguity {
        rule: AmbiguityRule::Mean,
        codes: vec![(b'X', vec![0]), (b'X', vec![1])],
        unknown: Some(b'X'),
    };
    assert!(
        Program::compile(
            &[motif.view()],
            &Alphabet::dna(),
            &Background::uniform(&Alphabet::dna()),
            &scoring,
            Strands::Forward,
            Some(1),
        )
        .is_err()
    );
}

#[test]
fn finite_large_counts_and_smoothing_weights_keep_the_same_model() {
    let alphabet = Alphabet::dna();
    let background = Background::uniform(&alphabet);
    let huge = resident(&Motif::pfm(&[vec![1e308, 5e307, 2.5e307, 2.5e307]], None));
    let ordinary = resident(&Motif::pfm(&[vec![4.0, 2.0, 1.0, 1.0]], None));
    let huge = huge.motif(0).unwrap();
    let ordinary = ordinary.motif(0).unwrap();
    let mut scoring = Scoring::exact();
    scoring.pseudocount = 0.0;
    let compile = |view, policy: &Scoring| {
        Program::compile(
            &[view],
            &alphabet,
            &background,
            policy,
            Strands::Forward,
            Some(1),
        )
        .unwrap()
    };
    let a = compile(huge.view(), &scoring);
    let b = compile(ordinary.view(), &scoring);
    assert_eq!(
        a.motif(0).strands()[0].full().score_model(),
        b.motif(0).strands()[0].full().score_model()
    );
    scoring.default_sites = f64::MAX;
    scoring.pseudocount = f64::MAX;
    let a = compile(ordinary.view(), &scoring);
    scoring.default_sites = 1.0;
    scoring.pseudocount = 1.0;
    let b = compile(ordinary.view(), &scoring);
    assert_eq!(
        a.motif(0).strands()[0].full().score_model(),
        b.motif(0).strands()[0].full().score_model()
    );
}

#[test]
fn single_column_scores_are_log_odds_of_counts_and_background() {
    let reader = resident(&Motif::pfm(&[vec![4.0, 2.0, 1.0, 1.0]], None));
    let motif = reader.motif(0).unwrap();
    let alphabet = Alphabet::dna();
    let mut scoring = Scoring::exact();
    scoring.pseudocount = 0.0;
    let program = Program::compile(
        &[motif.view()],
        &alphabet,
        &Background::frequencies(&[0.1, 0.2, 0.3, 0.4]).unwrap(),
        &scoring,
        Strands::Forward,
        Some(1),
    )
    .unwrap();
    let hits = program
        .scan_sequence(
            b"ACGT",
            Hits::new(omicsweave_motif::Threshold::Score(f64::NEG_INFINITY)),
            &ScanOptions::default(),
        )
        .unwrap();
    assert_eq!(hits.len(), 4);
    let model = program.motif(0).strands()[0].full().int_model();
    for (hit, (probability, background)) in
        hits.iter()
            .zip([(0.5_f64, 0.1_f64), (0.25, 0.2), (0.125, 0.3), (0.125, 0.4)])
    {
        let expected = (probability / background).log2();
        assert!(
            (hit.score - expected).abs() <= 0.5 / model.scale + 1e-12,
            "base={}, expected={expected}, hit={hit:?}",
            hit.start
        );
    }
}
