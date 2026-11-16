use super::*;

#[test]
fn meme_keeps_probabilities_ids_and_nsites() {
    let dir = tempfile::tempdir().unwrap();
    let src = write(dir.path(), "m.meme", MEME);
    assert_eq!(Format::detect(&src).unwrap(), Format::Meme);
    let (out, _) = run(dir.path(), &[src]);
    assert_eq!(out.alphabet, b"ACGT");
    let (id, name, w, sites, f) = &out.motifs[0];
    assert_eq!(
        (id.as_str(), name.as_str(), *w, *sites),
        ("MA0002.1", "RUNX1", 2, Some(26.0))
    );
    assert_eq!(f[1], (FactorKind::Joint, vec![1], vec![0.1, 0.2, 0.3, 0.4]));
    assert_eq!(out.motifs[1].0, "lone");
    assert_eq!(out.motifs[1].3, None);
}

#[test]
fn meme_protein_and_rna_alphabets() {
    let dir = tempfile::tempdir().unwrap();
    let row = vec!["0.05"; 20].join(" ");
    let protein = format!(
        "MEME version 5\n\nALPHABET= ACDEFGHIKLMNPQRSTVWY\n\nMOTIF p\nletter-probability matrix: alength= 20 w= 1\n{row}\n"
    );
    let src = write(dir.path(), "p.meme", &protein);
    let (out, _) = run(dir.path(), &[src]);
    assert_eq!(out.alphabet, b"ACDEFGHIKLMNPQRSTVWY");
    let rna = "MEME version 4\nALPHABET= ACGU\nMOTIF r\nletter-probability matrix: w= 1\n0.1 0.2 0.3 0.4\n";
    let src = write(dir.path(), "r.meme", rna);
    let (out, _) = run(dir.path(), &[src]);
    assert_eq!(out.alphabet, b"ACGU");
}

#[test]
fn jaspar_brackets_and_raw_pfm_store_counts() {
    let dir = tempfile::tempdir().unwrap();
    let brackets = ">MA0001.1\tAGL3\nA  [ 0  3 ]\nC  [ 94 75 ]\nG  [ 1  4 ]\nT  [ 2 15 ]\n";
    let src = write(dir.path(), "a.jaspar", brackets);
    assert_eq!(Format::detect(&src).unwrap(), Format::Jaspar);
    let (out, _) = run(dir.path(), &[src]);
    let (id, name, w, sites, f) = &out.motifs[0];
    assert_eq!(
        (id.as_str(), name.as_str(), *w, *sites),
        ("MA0001.1", "AGL3", 2, Some(97.0))
    );
    assert_eq!(f[0].2, vec![0.0, 94.0, 1.0, 2.0]);

    let raw = ">AHR.H14CORE.0.P.B\n176 188 102\n307 363 133\n392 152 377\n128 300 391\n";
    let src = write(dir.path(), "a.txt", raw);
    assert_eq!(Format::detect(&src).unwrap(), Format::JasparPfm);
    let (out, _) = run(dir.path(), &[src]);
    assert_eq!(out.motifs[0].4[2].2, vec![102.0, 133.0, 377.0, 391.0]);
    assert_eq!(out.motifs[0].3, Some(1003.0));
}

#[test]
fn transfac_reorders_columns_and_takes_accession() {
    let dir = tempfile::tempdir().unwrap();
    let text = "AC  MA0002.1\nXX\nID  RUNX1\nXX\nDE  MA0002.1 RUNX1 ; From JASPAR\n\
                PO      T      G      C      A\n01      1      2      3      4      M\n\
                02      4      3      2      1      K\nXX\n//\n";
    let src = write(dir.path(), "x.transfac", text);
    assert_eq!(Format::detect(&src).unwrap(), Format::Transfac);
    let (out, _) = run(dir.path(), &[src]);
    let (id, name, _, sites, f) = &out.motifs[0];
    assert_eq!(
        (id.as_str(), name.as_str(), *sites),
        ("MA0002.1", "RUNX1", Some(10.0))
    );
    assert_eq!(f[0].2, vec![4.0, 3.0, 2.0, 1.0]);
}

#[test]
fn homer_takes_the_name_column() {
    let dir = tempfile::tempdir().unwrap();
    let text = ">TGTGGTTW\tRUNX(Runt)/HPC7-Runx1-ChIP-Seq/Homer\t6.509\t-1.886e+03\n\
                0.1 0.2 0.3 0.4\n0.4 0.3 0.2 0.1\n";
    let src = write(dir.path(), "r.motif", text);
    assert_eq!(Format::detect(&src).unwrap(), Format::Homer);
    let (out, _) = run(dir.path(), &[src]);
    assert_eq!(out.motifs[0].0, "RUNX(Runt)/HPC7-Runx1-ChIP-Seq/Homer");
    assert_eq!(out.motifs[0].2, 2);
}

#[test]
fn cisbp_files_form_one_collection_and_shells_are_passed_over() {
    let dir = tempfile::tempdir().unwrap();
    let pwms = dir.path().join("pwms");
    std::fs::create_dir(&pwms).unwrap();
    write(
        &pwms,
        "M1_3.10.txt",
        "Pos\tA\tC\tG\tT\n1\t0.1\t0.2\t0.3\t0.4\n2\t0.25\t0.25\t0.25\t0.25\n",
    );
    write(&pwms, "M2_3.10.txt", "Pos\tA\tC\tG\tT\n");
    write(&pwms, "M3_3.10.txt", "Pos\tA\tC\tG\tU\n1\t1\t0\t0\t0\n");
    let (m1, m2) = (pwms.join("M1_3.10.txt"), pwms.join("M2_3.10.txt"));
    let (out, warnings) = run(dir.path(), &[m1.clone(), m2, m1]);
    assert_eq!(
        out.motifs.len(),
        1,
        "M1 once despite being given twice; M2 is a shell"
    );
    assert_eq!(out.motifs[0].0, "M1_3.10");
    assert_eq!(warnings.len(), 1);
    // DNA and RNA files do not form one collection.
    let all: Vec<PathBuf> = ["M1_3.10.txt", "M3_3.10.txt"].map(|n| pwms.join(n)).into();
    let from = detect_sources(&all).unwrap();
    let mut c =
        Conversion::between(from, Format::Ffm, all[0].clone(), dir.path().join("x.ffm")).unwrap();
    c.add_sources(all[1..].iter().cloned()).unwrap();
    assert!(convert(c, Budget::new(1, 1 << 30), None).is_err());
}

#[test]
fn meme_custom_alphabet_block_and_guessed_alphabets() {
    let dir = tempfile::tempdir().unwrap();
    let text = "MEME version 5\n\nALPHABET \"methyl\" DNA-LIKE\n\
                A \"Adenine\" CC0000 ~ T \"Thymine\" 008000\nC ~ G # a pair\nm ~ 1\n\
                N = ACGT\nU = T\nEND ALPHABET\n\n\
                MOTIF x\nletter-probability matrix: alength= 6 w= 1\n0.1 0.1 0.1 0.1 0.2 0.4\n";
    let src = write(dir.path(), "c.meme", text);
    let (out, _) = run(dir.path(), std::slice::from_ref(&src));
    // Columns in MEME's order: uppercase, lowercase, digits.
    assert_eq!(out.alphabet, b"ACGTm1");
    let r = FfmReader::open(std::fs::read(dir.path().join("out.ffm")).unwrap()).unwrap();
    let a = r.alphabet();
    assert_eq!(a.complement(), Some(&[3u8, 2, 1, 0, 5, 4][..]));
    assert_eq!(a.code(b'U'), Some(3));
    assert_eq!(a.code(b'N'), None, "ambiguous symbols are invalid input");

    let row = vec!["0.05"; 20].join(" ");
    let text =
        format!("MEME version 4\nMOTIF p\nletter-probability matrix: alength= 20 w= 1\n{row}\n");
    let src = write(dir.path(), "g.meme", &text);
    assert_eq!(run(dir.path(), &[src]).0.alphabet, b"ACDEFGHIKLMNPQRSTVWY");
}
