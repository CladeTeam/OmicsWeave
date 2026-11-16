use crate::support::scan_fixtures::{Rng, dinucleotide, markov, pwm, reader, subset};
use omicsweave_ffm::Alphabet;
use omicsweave_motif::{Background, Scoring, Strands, Threshold};

fn real_score(model: &omicsweave_motif::ScoreModel, codes: &[u8], at: usize) -> f64 {
    model
        .terms
        .iter()
        .map(|t| {
            let idx = t.positions.iter().fold(0usize, |acc, &p| {
                acc * 4 + codes[(at as i64 + p as i64) as usize] as usize
            });
            t.table[idx]
        })
        .sum()
}

#[test]
fn null_distributions_match_enumeration() {
    let mut rng = Rng(43);
    let motifs = vec![
        pwm(&mut rng, 4),
        dinucleotide(&mut rng, 4),
        markov(&mut rng, 4),
        subset(&mut rng),
    ];
    let r = reader(&motifs);
    let alphabet = Alphabet::dna();
    for bg in [
        Background::frequencies(&[0.1, 0.4, 0.3, 0.2]).unwrap(),
        Background::from_kmer_frequencies(
            4,
            vec![
                (0..4).map(|_| 0.5 + rng.unit()).collect(),
                (0..16).map(|_| 0.5 + rng.unit()).collect(),
            ],
        )
        .unwrap(),
    ] {
        for i in 0..r.len() {
            let mref = r.motif(i).unwrap();
            let model = omicsweave_motif::MotifModel::compile(
                &mref.view(),
                &alphabet,
                &bg,
                &Scoring::exact(),
                Strands::Both,
            )
            .unwrap();
            for sm in model.strands() {
                for (c, v) in sm.variants().iter().enumerate() {
                    let int = v.int_model();
                    let lo = int.lo();
                    let len = (int.width as i32 - lo) as usize;
                    let mut pmf = vec![0.0; int.max as usize + 1];
                    let mut forbidden = 0.0;
                    let mut codes = vec![0u8; len];
                    for idx in 0..4usize.pow(len as u32) {
                        let mut x = idx;
                        for cd in codes.iter_mut().rev() {
                            *cd = (x % 4) as u8;
                            x /= 4;
                        }
                        // Chain restarted at the footprint's first position.
                        let mut pr = 1.0;
                        for h in 0..len {
                            let ctx = bg.order().min(h);
                            let cidx = codes[h - ctx..h]
                                .iter()
                                .fold(0usize, |a, &d| a * 4 + d as usize);
                            pr *= bg.conditional(ctx)[cidx * 4 + codes[h] as usize];
                        }
                        let at = (-lo) as usize;
                        match int.score(&codes, at) {
                            Some(s) => {
                                pmf[s as usize] += pr;
                                // Quantization stays within its rounding bound.
                                let real = real_score(v.score_model(), &codes, at);
                                let bound = int.terms.len() as f64 / (2.0 * int.scale) + 1e-9;
                                assert!((int.raw(s) - real).abs() <= bound);
                            }
                            None => forbidden += pr,
                        }
                    }
                    let null = v.null();
                    for (a, b) in pmf.iter().zip(null.pmf()) {
                        assert!((a - b).abs() < 1e-12, "motif {i} context {c}: {a} vs {b}");
                    }
                    assert!((forbidden - null.forbidden()).abs() < 1e-12);
                }
            }
        }
    }
}

use crate::support::resident::resident;
use omicsweave_ffm::{Factor, FactorKind, Motif, Support};
use omicsweave_motif::{Program, Quantizer};

#[test]
fn tiny_rounded_pvalues_keep_cutoff_search_monotone() {
    let reader = resident(&Motif {
        width: 1,
        support: Support::All,
        sites: None,
        factors: vec![Factor {
            kind: FactorKind::Potential,
            positions: vec![0],
            table: vec![1.0, 1.0, 1.0, 2.0],
        }],
        gaps: None,
    });
    let motif = reader.motif(0).unwrap();
    let mut scoring = Scoring::exact();
    scoring.quantizer = Quantizer::Precision { states: 128 };
    scoring.rounding = omicsweave_motif::Rounding::Significant(6);
    let background = Background::frequencies(&[0.2, 0.3, 0.5, 1e-310]).unwrap();
    let program = Program::compile(
        &[motif.view()],
        &Alphabet::dna(),
        &background,
        &scoring,
        Strands::Forward,
        Some(1),
    )
    .unwrap();
    let variant = program.motif(0).strands()[0].full();
    let pvalues: Vec<f64> = (0..=variant.int_model().max)
        .map(|s| variant.pvalue(s))
        .collect();
    assert!(pvalues.iter().all(|p| p.is_finite()));
    assert!(pvalues.windows(2).all(|p| p[0] >= p[1]));
    assert_eq!(variant.cutoff(Threshold::PValue(2e-310)), Some(1));
    for threshold in [5e-311, 1e-310, 2e-310, 1.0] {
        let expected = pvalues
            .iter()
            .position(|&p| p <= threshold)
            .map(|s| s as i32);
        assert_eq!(variant.cutoff(Threshold::PValue(threshold)), expected);
    }
}
