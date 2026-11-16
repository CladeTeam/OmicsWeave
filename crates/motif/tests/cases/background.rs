use omicsweave_motif::Background;

#[test]
fn conditional_rows_normalize_independently_and_survive_large_counts() {
    let background =
        Background::from_kmer_frequencies(2, vec![vec![1e308, 1e308], vec![1.0, 3.0, 6.0, 2.0]])
            .unwrap();
    assert_eq!(background.letters(), &[0.5, 0.5]);
    assert_eq!(background.conditional(1), &[0.25, 0.75, 0.75, 0.25]);
}

#[test]
fn malformed_frequencies_are_rejected() {
    for (n, frequencies) in [
        (0, vec![vec![]]),
        (256, vec![vec![1.0; 256]]),
        (4, vec![]),
        (4, vec![vec![1.0; 3]]),
        (4, vec![vec![1.0; 4], vec![1.0; 15]]),
    ] {
        assert!(
            Background::from_kmer_frequencies(n, frequencies.clone()).is_err(),
            "alphabet={n}, frequencies={frequencies:?}"
        );
    }
    for bad in [0.0, -1.0, f64::INFINITY, f64::NAN] {
        assert!(
            Background::frequencies(&[1.0, 1.0, 1.0, bad]).is_err(),
            "frequency={bad}"
        );
    }
}
