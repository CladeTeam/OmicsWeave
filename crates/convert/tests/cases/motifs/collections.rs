use super::*;

#[test]
fn ffm_sources_merge() {
    let dir = tempfile::tempdir().unwrap();
    let a = write(dir.path(), "a.meme", MEME);
    let (_, _) = run(dir.path(), &[a]);
    std::fs::rename(dir.path().join("out.ffm"), dir.path().join("a.ffm")).unwrap();
    let b = write(dir.path(), "b.pcm", ">other\n1 2 3 4\n4 3 2 1\n3 3 3 3\n");
    let (_, _) = run(dir.path(), &[b]);
    std::fs::rename(dir.path().join("out.ffm"), dir.path().join("b.ffm")).unwrap();
    let (out, _) = run(
        dir.path(),
        &[dir.path().join("a.ffm"), dir.path().join("b.ffm")],
    );
    let ids: Vec<&str> = out.motifs.iter().map(|m| m.0.as_str()).collect();
    assert_eq!(ids, ["MA0002.1", "lone", "other"]);
    assert_eq!(out.motifs[0].1, "RUNX1");
}

#[test]
fn detection_keeps_fasta_and_refuses_square_matrices() {
    let dir = tempfile::tempdir().unwrap();
    let fasta = write(dir.path(), "s.fa", ">chr1\nACGTACGT\n>chr2\nNNNN\n");
    assert_eq!(Format::detect(&fasta).unwrap(), Format::Fasta);
    let square = write(
        dir.path(),
        "q.txt",
        ">q\n1 2 3 4\n1 2 3 4\n1 2 3 4\n1 2 3 4\n",
    );
    assert!(
        Format::detect(&square).is_err(),
        "4 x 4 matrices need their format named"
    );
    assert!(detect_sources(&[square]).is_err());
}

#[test]
fn conflicts_are_errors() {
    let dir = tempfile::tempdir().unwrap();
    let a = write(dir.path(), "a.meme", MEME);
    let b = write(dir.path(), "b.meme", MEME);
    let dst = dir.path().join("o.ffm");
    let mut c = Conversion::between(Format::Meme, Format::Ffm, a.clone(), &dst).unwrap();
    c.add_sources([b]).unwrap();
    assert!(
        convert(c, Budget::new(1, 1 << 30), None).is_err(),
        "same IDs twice"
    );

    let pcm = write(dir.path(), "x.pcm", ">m\n1 2 3 4\n1 2 3 4\n1 2 3 4\n");
    assert!(detect_sources(&[a, pcm]).is_err(), "formats differ");
    assert!(detect_sources(&[dir.path().join("none.meme")]).is_err());

    let bad = write(
        dir.path(),
        "bad.meme",
        "MEME version 4\nMOTIF x\nletter-probability matrix: w= 2\n0.5 0.5 0 0\n",
    );
    let c = Conversion::between(Format::Meme, Format::Ffm, bad, &dst).unwrap();
    assert!(
        convert(c, Budget::new(1, 1 << 30), None).is_err(),
        "w= says 2 rows"
    );

    let mut fasta = Conversion::between(
        Format::Fasta,
        Format::TwoBit,
        dir.path().join("a.fa"),
        dir.path().join("a.2bit"),
    )
    .unwrap();
    assert!(fasta.add_sources([dir.path().join("b.fa")]).is_err());
}

#[test]
fn empty_sources_and_malformed_models_return_errors() {
    let dir = tempfile::tempdir().unwrap();
    let dst = dir.path().join("out.ffm");
    assert!(detect_sources(&[]).is_err());
    let src = write(dir.path(), "valid.meme", MEME);
    let mut empty = Conversion::between(Format::Meme, Format::Ffm, &src, &dst).unwrap();
    if let Conversion::MotifToFfm(c) = &mut empty {
        c.sources.clear();
    }
    assert!(convert(empty, Budget::new(1, 1 << 20), None).is_err());

    for (format, text) in [
        (Format::Jaspar, ">bad\n]A[1\n"),
        (
            Format::Meme,
            "MEME version 5\nALPHABET \"custom\"\n\"label only\"\nEND ALPHABET\n",
        ),
        (
            Format::Meme,
            "MEME version 5\nALPHABET= ACGT\nMOTIF bad\nletter-probability matrix: w= 1.5\n1 0 0 0\n",
        ),
        (
            Format::Meme,
            "MEME version 5\nALPHABET= ACGT\nMOTIF bad\nletter-probability matrix: alength= 4.5 w= 1\n1 0 0 0\n",
        ),
        (
            Format::Glam2,
            "GLAM2: Gapped Local Alignment of Motifs\nScore: 1 Columns: 18446744073709551615 Sequences: 1\n",
        ),
        (
            Format::ProBound,
            r#"{"modelSettings":{"letterOrder":"ACGT","bindingModes":[{"size":18446744073709551615}]},"coefficients":{"bindingModes":[{"mononucleotide":[]}]}}"#,
        ),
    ] {
        let src = write(dir.path(), "invalid.txt", text);
        let c = Conversion::between(format, Format::Ffm, src, &dst).unwrap();
        assert!(
            convert(c, Budget::new(1, 1 << 20), None).is_err(),
            "{format}"
        );
    }
}

#[test]
fn headed_detection_keeps_complete_final_records_without_newlines() {
    let dir = tempfile::tempdir().unwrap();
    let square = write(
        dir.path(),
        "square",
        ">m\n1 2 3 4\n1 2 3 4\n1 2 3 4\n1 2 3 4",
    );
    assert!(Format::detect(&square).is_err());
    let raw = write(dir.path(), "raw", ">m\n1 2\n3 4\n5 6\n7 8");
    assert_eq!(Format::detect(&raw).unwrap(), Format::JasparPfm);
    let mixed = write(dir.path(), "mixed", ">m\n1 2 3 4\n>n\n1 2\n3 4\n5 6\n7 8\n");
    assert!(Format::detect(&mixed).is_err());
}

#[test]
fn gzip_expansion_respects_the_text_budget() {
    use std::io::Write as _;
    let dir = tempfile::tempdir().unwrap();
    let src = dir.path().join("large.pcm.gz");
    let file = std::fs::File::create(&src).unwrap();
    let mut encoder = flate2::write::GzEncoder::new(file, flate2::Compression::default());
    encoder.write_all(b">m\n").unwrap();
    encoder.write_all(&b"1 0 0 0\n".repeat(10_000)).unwrap();
    encoder.finish().unwrap();
    let dst = dir.path().join("out.ffm");
    std::fs::write(&dst, b"preserved").unwrap();
    let c = Conversion::between(Format::Pcm, Format::Ffm, src, &dst).unwrap();
    let error = convert(c, Budget::new(1, 4096), None).unwrap_err();
    assert!(error.to_string().contains("decoded motif text"), "{error}");
    assert_eq!(std::fs::read(&dst).unwrap(), b"preserved");
}
