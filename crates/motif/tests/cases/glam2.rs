use crate::support::glam2::*;
use omicsweave_ffm::{Alphabet, Motif};
use omicsweave_motif::profile::{
    AllAlignments, DirichletMixture, GLAM2_HITS, ProfileScoring, TopAlignments, format_alt_g,
    glam2_dna_alphabet, write_glam2_header,
};
use omicsweave_motif::{ScanOptions, Strands};

#[test]
fn reproduces_a_glam2scan_transcript() {
    let alphabet = glam2_dna_alphabet();
    let motif = parse_glam2_motif(HAND_MOTIF, &alphabet);
    let scoring = ProfileScoring::glam2(&alphabet).unwrap();
    let p = program(&[motif], &alphabet, &scoring, Strands::Both).unwrap();
    let fasta = Fasta::parse(HAND_FASTA);
    for threads in [Some(1), Some(3), None] {
        assert_eq!(hits_text(&p, &fasta, 7, threads), HAND_HITS);
    }
    let mut header = Vec::new();
    write_glam2_header(&mut header, "1056", "glam2scan -t -2 -n 7 n m.glam2 s.fa").unwrap();
    assert_eq!(
        String::from_utf8(header).unwrap(),
        "GLAM2SCAN\nVersion 1056\n\nglam2scan -t -2 -n 7 n m.glam2 s.fa\n\n"
    );
}

#[test]
fn hit_streams_never_rise_and_motifs_scan_independently() {
    let alphabet = glam2_dna_alphabet();
    let mut rng = Rng(99);
    let motifs: Vec<Motif> = (0..3)
        .map(|_| parse_glam2_motif(&random_motif(&mut rng, b"acgt", 8, 12), &alphabet))
        .collect();
    let fasta = Fasta::parse(&random_fasta(&mut rng, b"acgtACGTNn", 20, 120));
    let scoring = ProfileScoring::glam2(&alphabet).unwrap();
    let all = program(&motifs, &alphabet, &scoring, Strands::Both).unwrap();
    let hits = all
        .scan(
            fasta.sequences(),
            AllAlignments::new(None),
            &ScanOptions::default(),
        )
        .unwrap();
    assert!(!hits.is_empty());
    for pair in hits.windows(2) {
        let (a, b) = (&pair[0], &pair[1]);
        if (a.sequence, a.motif, a.strand) == (b.sequence, b.motif, b.strand) {
            assert!(b.score <= a.score);
        }
    }
    let together = all
        .scan(
            fasta.sequences(),
            TopAlignments::glam2(GLAM2_HITS, fasta.names()),
            &ScanOptions::default(),
        )
        .unwrap();
    for (i, m) in motifs.iter().enumerate() {
        let alone = program(std::slice::from_ref(m), &alphabet, &scoring, Strands::Both).unwrap();
        let mut single = alone
            .scan(
                fasta.sequences(),
                TopAlignments::glam2(GLAM2_HITS, fasta.names()),
                &ScanOptions::default(),
            )
            .unwrap();
        for h in &mut single[0] {
            h.motif = i as u32;
        }
        assert_eq!(single[0], together[i]);
    }
}

#[test]
fn ungapped_motifs_are_rejected() {
    let alphabet = Alphabet::dna();
    let pfm = Motif::pfm(&[vec![1.0, 2.0, 3.0, 4.0]], Some(10.0));
    let scoring = ProfileScoring::glam2(&alphabet).unwrap();
    let err = program(&[pfm], &alphabet, &scoring, Strands::Forward).unwrap_err();
    assert!(err.to_string().contains("gap"), "{err}");
}

#[test]
fn priors_and_mixture_files() {
    let p = DirichletMixture::glam2_dna().posterior_mean(&[3, 0, 1, 0]);
    for (x, n) in p.iter().zip([3.0, 0.0, 1.0, 0.0]) {
        assert!((x - (n + 0.4) / 5.6).abs() < 1e-15);
    }
    let protein = DirichletMixture::glam2_protein();
    assert_eq!(protein.weights.len(), 20);
    assert_eq!(protein.dimension(), 20);
    let mean = protein.posterior_mean(&[0; 20]);
    assert!((mean.iter().sum::<f64>() - 1.0).abs() < 1e-12);

    let text = "Mixture= 0.7\nAlpha= 1.2 0.3 0.3 0.3 0.3\nMixture= 0.3\nAlpha= 9 5.0 1.0 2.0 1.0\n";
    let m = DirichletMixture::parse_glam2(text, 4).unwrap();
    assert_eq!(m.weights, vec![0.7, 0.3]);
    assert_eq!(m.alphas[1], vec![5.0, 1.0, 2.0, 1.0]);
    assert!(DirichletMixture::parse_glam2(text, 20).is_err());
    assert!(DirichletMixture::parse_glam2("Mixture= 1\n", 4).is_err());
    assert!(DirichletMixture::parse_glam2("Mixture= 1\nAlpha= 1 1 0 1 1\n", 4).is_err());
    assert!(DirichletMixture::parse_glam2("Mixture= 1\nAlpha= 1\n", 4).is_err());
}

#[test]
fn alternate_g_formatting_matches_c() {
    assert_eq!(format_alt_g(9.5712, 3), "9.57");
    assert_eq!(format_alt_g(-3.7, 3), "-3.70");
    assert_eq!(format_alt_g(0.158, 3), "0.158");
    assert_eq!(format_alt_g(123.4, 3), "123.");
    assert_eq!(format_alt_g(1234.0, 3), "1.23e+03");
    assert_eq!(format_alt_g(999.7, 3), "1.00e+03");
    assert_eq!(format_alt_g(1e-5, 3), "1.00e-05");
    assert_eq!(format_alt_g(0.0, 3), "0.00");
    assert_eq!(format_alt_g(5.0, 1), "5.");
}

#[test]
fn invalid_scoring_and_short_name_tables_are_reported() {
    use omicsweave_core::Strings;
    let alphabet = glam2_dna_alphabet();
    let motif = parse_glam2_motif(HAND_MOTIF, &alphabet);
    for value in [0.0, -1.0, f64::INFINITY, f64::NAN] {
        let mut scoring = ProfileScoring::glam2(&alphabet).unwrap();
        scoring.delete_pseudo = value;
        let error = program(
            std::slice::from_ref(&motif),
            &alphabet,
            &scoring,
            Strands::Both,
        )
        .unwrap_err();
        assert!(
            error.to_string().contains("pseudocount"),
            "delete_pseudo={value}: {error}"
        );
    }
    let p = program(
        &[motif],
        &alphabet,
        &ProfileScoring::glam2(&alphabet).unwrap(),
        Strands::Both,
    )
    .unwrap();
    let fasta = Fasta::parse(HAND_FASTA);
    let error = p
        .scan(
            fasta.sequences(),
            TopAlignments::glam2(
                7,
                Strings {
                    offsets: &[0],
                    data: b"",
                },
            ),
            &ScanOptions::default(),
        )
        .unwrap_err();
    assert!(error.to_string().contains("names"), "{error}");
    assert_eq!(hits_text(&p, &fasta, 7, Some(1)), HAND_HITS);
    let none = p
        .scan(
            fasta.sequences(),
            TopAlignments::glam2(0, fasta.names()),
            &ScanOptions::default(),
        )
        .unwrap();
    assert_eq!(none.len(), 1);
    assert!(none[0].is_empty());
}
