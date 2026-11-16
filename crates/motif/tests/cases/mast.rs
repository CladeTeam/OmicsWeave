use crate::support::resident::resident;
use omicsweave_core::{Strand, Strings};
use omicsweave_ffm::{Alphabet, Motif};
use omicsweave_motif::mast::{Mast, MastOptions, StrandMode};
use omicsweave_motif::{
    Background, Hits, Program, Quantizer, ScanOptions, Scoring, Strands, Threshold,
};

#[test]
fn mast_keeps_the_best_windows_own_strand_and_context_statistics() {
    let reader = resident(&Motif::pfm(
        &[vec![20.0, 1.0, 4.0, 2.0], vec![1.0, 3.0, 15.0, 1.0]],
        None,
    ));
    let motif = reader.motif(0).unwrap();
    let bg = Background::from_kmer_frequencies(
        4,
        vec![
            vec![1.0, 2.0, 3.0, 5.0],
            (1..=16).map(|x| (x * x) as f64).collect(),
        ],
    )
    .unwrap();
    let mut scoring = Scoring::exact();
    scoring.quantizer = Quantizer::Precision { states: 128 };
    let program = Program::compile(
        &[motif.view()],
        &Alphabet::dna(),
        &bg,
        &scoring,
        Strands::Both,
        Some(1),
    )
    .unwrap();
    let mut data = Vec::new();
    let mut offsets = vec![0i64];
    for length in [2, 4] {
        for mut word in 0..4usize.pow(length) {
            for _ in 0..length {
                data.push(b"ACGT"[word % 4]);
                word /= 4;
            }
            offsets.push(data.len() as i64);
        }
    }
    let seqs = Strings {
        offsets: &offsets,
        data: &data,
    };
    for mode in [
        StrandMode::Combine,
        StrandMode::Separate,
        StrandMode::Forward,
    ] {
        let mut options = MastOptions::new(mode);
        options.hit_threshold = 1.0;
        let (report, all) = program
            .scan(
                seqs,
                (
                    Mast::new(options),
                    Hits::new(Threshold::Score(f64::NEG_INFINITY)),
                ),
                &ScanOptions {
                    threads: Some(2),
                    tile_len: 1024,
                },
            )
            .unwrap();
        for (threads, tile_len) in [(1, 1024), (2, 8192)] {
            let mut options = MastOptions::new(mode);
            options.hit_threshold = 1.0;
            let other = program
                .scan(
                    seqs,
                    Mast::new(options),
                    &ScanOptions {
                        threads: Some(threads),
                        tile_len,
                    },
                )
                .unwrap();
            assert_eq!(other, report);
        }
        for seq in report.sequences {
            let mut expected = None;
            for hit in all.iter().filter(|h| h.sequence == seq.sequence) {
                if mode == StrandMode::Forward && hit.strand != Strand::Forward
                    || seq.strand.is_some_and(|s| s != hit.strand)
                {
                    continue;
                }
                if expected.is_none_or(|e: &omicsweave_motif::HitEntry| hit.score > e.score) {
                    expected = Some(hit);
                }
            }
            let expected = expected.unwrap();
            assert_eq!(seq.best.len(), 1);
            assert_eq!(seq.best[0].score, expected.score);
            assert_eq!(seq.best[0].pvalue, expected.pvalue);
            for hit in &seq.hits {
                let expected = all
                    .iter()
                    .find(|h| {
                        h.sequence == seq.sequence && h.start == hit.start && h.strand == hit.strand
                    })
                    .unwrap();
                assert_eq!(hit.score, expected.score);
                assert_eq!(hit.pvalue, expected.pvalue);
            }
        }
    }
}

#[test]
fn supplied_null_letters_follow_strand_policy_without_normalizing_weights() {
    use omicsweave_motif::mast::null_letters;
    let alphabet = Alphabet::dna();
    let given = [1.0, 2.0, 3.0, 4.0];
    assert_eq!(
        null_letters(&alphabet, Some(&given), StrandMode::Forward).unwrap(),
        given
    );
    for mode in [StrandMode::Combine, StrandMode::Separate] {
        assert_eq!(
            null_letters(&alphabet, Some(&given), mode).unwrap(),
            [2.5; 4]
        );
    }
    assert!(null_letters(&alphabet, Some(&given[..3]), StrandMode::Forward).is_err());
}
