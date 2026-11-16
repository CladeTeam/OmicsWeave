use crate::support::scan_fixtures::{Rng, dinucleotide, pwm, reader};
use crate::support::scan_oracle::Collect;
use omicsweave_core::Strand;
use omicsweave_ffm::{Alphabet, Motif, MotifRef};
use omicsweave_motif::{Background, Program, ScanOptions, Scoring, Strands};

#[test]
fn reverse_strand_is_the_reverse_complement() {
    let mut rng = Rng(53);
    let motifs: Vec<Motif> = vec![pwm(&mut rng, 8), dinucleotide(&mut rng, 6)];
    let r = reader(&motifs);
    let refs: Vec<MotifRef> = (0..r.len()).map(|i| r.motif(i).unwrap()).collect();
    let views: Vec<_> = refs.iter().map(|m| m.view()).collect();
    let alphabet = Alphabet::dna();
    let bg = Background::frequencies(&[0.3, 0.2, 0.2, 0.3]).unwrap();
    let program = Program::compile(
        &views,
        &alphabet,
        &bg,
        &Scoring::exact(),
        Strands::Both,
        Some(1),
    )
    .unwrap();
    let seq: Vec<u8> = (0..3000).map(|_| b"ACGT"[rng.below(4)]).collect();
    let mut rc = seq.clone();
    omicsweave_core::orient::reverse_complement(&mut rc);
    let opt = ScanOptions::default();
    let a = program
        .scan_sequence(&seq, Collect::default(), &opt)
        .unwrap();
    let b = program
        .scan_sequence(&rc, Collect::default(), &opt)
        .unwrap();
    let flip = |v: &[(u32, Strand, u32, u32, u8, i32)], strand: Strand| {
        let mut out: Vec<(u32, u32, i64)> = v
            .iter()
            .filter(|w| w.1 == strand)
            .map(|w| {
                let wd = program.motif(w.0 as usize).width();
                let raw = program
                    .motif(w.0 as usize)
                    .strand(strand)
                    .unwrap()
                    .full()
                    .int_model()
                    .raw(w.5);
                (w.0, 3000 - w.3 - wd, (raw * 1e6).round() as i64)
            })
            .collect();
        out.sort();
        out
    };
    let forward_on_rc: Vec<(u32, u32, i64)> = {
        let mut v: Vec<_> = b
            .iter()
            .filter(|w| w.1 == Strand::Forward)
            .map(|w| {
                let raw = program
                    .motif(w.0 as usize)
                    .strand(Strand::Forward)
                    .unwrap()
                    .full()
                    .int_model()
                    .raw(w.5);
                (w.0, w.3, (raw * 1e6).round() as i64)
            })
            .collect();
        v.sort();
        v
    };
    let reverse = flip(&a, Strand::Reverse);
    assert_eq!(reverse.len(), forward_on_rc.len());
    for (x, y) in reverse.iter().zip(&forward_on_rc) {
        assert_eq!((x.0, x.1), (y.0, y.1));
        assert!((x.2 - y.2).abs() <= 2 * 10_000, "{x:?} vs {y:?}");
    }
}
