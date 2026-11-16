use super::*;

#[test]
fn position_rows_weights_and_dinucleotides() {
    let dir = tempfile::tempdir().unwrap();
    let pcm = write(dir.path(), "a.pcm", ">m1 Foo\n10 0 0 0\n0 5 5 0\n1 1 1 7\n");
    assert_eq!(Format::detect(&pcm).unwrap(), Format::Pcm);
    let (out, _) = run(dir.path(), &[pcm]);
    assert_eq!(
        (out.motifs[0].1.as_str(), out.motifs[0].3),
        ("Foo", Some(10.0))
    );

    let pwm = write(dir.path(), "a.pwm", ">w\n-1 0.5 0 2\n0 0 0 0\n1 1 -1 -1\n");
    assert_eq!(Format::detect(&pwm).unwrap(), Format::Pwm);
    let (out, _) = run(dir.path(), &[pwm]);
    let (kind, _, table) = &out.motifs[0].4[0];
    assert_eq!(*kind, FactorKind::Potential);
    assert!((table[0].ln() + 1.0).abs() < 1e-12);

    let row: Vec<String> = (0..16).map(|i| i.to_string()).collect();
    let row = row.join(" ");
    let dpcm = write(dir.path(), "a.dpcm", &format!(">d\n{row}\n{row}\n"));
    assert_eq!(Format::detect(&dpcm).unwrap(), Format::Dpcm);
    let (out, _) = run(dir.path(), &[dpcm]);
    let m = &out.motifs[0];
    assert_eq!(m.2, 3);
    assert_eq!(m.4[1].1, vec![1, 2]);
    assert_eq!(
        m.4[1].2[5], 5.0,
        "AC AG ... TT, first position varying slowest"
    );
}

#[test]
fn ffm_conversion_keeps_support_factors_gaps_and_names() {
    use omicsweave_ffm::{
        Alphabet, Factor, FfmWriter, Gaps, Key, KeyType, Motif, Schema, Support, WriterOptions,
    };

    let dir = tempfile::tempdir().unwrap();
    let src = dir.path().join("models.ffm");
    let mut gapped = Motif::pfm(&[vec![2.0, 1.0, 1.0, 0.0], vec![0.0; 4]], Some(4.0));
    gapped.gaps = Some(Gaps {
        deletions: vec![0.0, 4.0],
        insertions: vec![2.0],
    });
    let conditional = Motif {
        width: 2,
        support: Support::Codes(vec![0, 2]),
        sites: Some(12.0),
        factors: vec![
            Factor {
                kind: FactorKind::Joint,
                positions: vec![0],
                table: vec![2.0, 1.0, 1.0],
            },
            Factor {
                kind: FactorKind::Conditional,
                positions: vec![0, 1],
                table: vec![1.0, 2.0, 3.0, 3.0, 2.0, 1.0, 1.0, 1.0, 1.0],
            },
        ],
        gaps: None,
    };
    let models = [gapped, conditional];
    let schema = Schema::new(vec![Key::new("name", KeyType::String)]).unwrap();
    let mut writer = FfmWriter::new(
        std::fs::File::create(&src).unwrap(),
        Alphabet::dna(),
        schema,
        WriterOptions::default(),
    )
    .unwrap();
    for (i, model) in models.iter().enumerate() {
        writer
            .push(&format!("model{i}"), model, &["display name"])
            .unwrap();
    }
    drop(writer.finish().unwrap());

    run(dir.path(), &[src]);
    let reader = FfmReader::open(std::fs::read(dir.path().join("out.ffm")).unwrap()).unwrap();
    for (i, model) in models.iter().enumerate() {
        assert_eq!(reader.motif(i).unwrap().view().to_motif(), *model);
        assert_eq!(reader.id(i), Some(format!("model{i}").as_str()));
        assert_eq!(reader.annotation(i).unwrap().get(0), Some("display name"));
    }
}

#[test]
fn bamm_keeps_true_orders_at_the_motif_start() {
    let dir = tempfile::tempdir().unwrap();
    let p0 = "2.5e-01 2.5e-01 2.5e-01 2.5e-01 ";
    let tiled = [p0.trim(); 4].join(" ");
    let order1: Vec<String> = (0..16)
        .map(|i| format!("{:.3e}", [0.1, 0.2, 0.3, 0.4][i % 4]))
        .collect();
    let text = format!(
        "{p0}\n{tiled}\n\n1.000e-01 2.000e-01 3.000e-01 4.000e-01 \n{}\n\n",
        order1.join(" ")
    );
    let src = write(dir.path(), "JunD_motif_1.ihbcp", &text);
    assert_eq!(Format::detect(&src).unwrap(), Format::Bamm);
    let (out, _) = run(dir.path(), &[src]);
    let (id, _, w, _, f) = &out.motifs[0];
    assert_eq!((id.as_str(), *w), ("JunD_motif_1", 2));
    assert_eq!(
        (f[0].0, f[0].1.clone(), f[0].2.len()),
        (FactorKind::Conditional, vec![0], 4)
    );
    assert_eq!((f[1].1.clone(), f[1].2.len()), (vec![0, 1], 16));
    assert_eq!(f[1].2[6], 0.3);
}

#[test]
fn probound_modes_become_potentials() {
    let dir = tempfile::tempdir().unwrap();
    let mono: Vec<String> = (0..12).map(|i| format!("{}", -(i as f64) / 10.0)).collect();
    let di: Vec<String> = (0..32).map(|i| format!("{}", i as f64 / 100.0)).collect();
    let text = format!(
        r#"{{"modelSettings":{{"letterOrder":"ACGT","letterComplement":"C-G,A-T","bindingModes":[{{"size":0}},{{"size":3,"dinucleotideDistance":1,"singleStrand":false}}]}},
"coefficients":{{"bindingModes":[{{"mononucleotide":[],"dinucleotide":[]}},{{"mononucleotide":[{}],"dinucleotide":[[{}]],"activity":[[0]]}}]}}}}"#,
        mono.join(","),
        di.join(",")
    );
    let src = write(dir.path(), "fit_15413.json", &text);
    assert_eq!(Format::detect(&src).unwrap(), Format::ProBound);
    let (out, _) = run(dir.path(), &[src]);
    let (id, _, w, _, f) = &out.motifs[0];
    assert_eq!((id.as_str(), *w, f.len()), ("fit_15413", 3, 5));
    assert!((f[1].2[2].ln() + 0.6).abs() < 1e-12, "position 1, G: -0.6");
    assert_eq!(f[4].1, vec![1, 2]);
    assert!(
        (f[4].2[1].ln() - 0.17).abs() < 1e-12,
        "second pair, A then C"
    );
}

const GLAM2: &str = "GLAM2: Gapped Local Alignment of Motifs
Version 5.5.9

glam2 n seqs.fa
Sequences: 4
Greatest sequence length: 30
Residue counts: a=30 c=30 g=30 t=30 N=0

Score: 21.5  Columns: 3  Sequences: 3

                **..*
s1           2 AC..G   4 + 7.10
s2           5 A.TTG   9 + 6.20
s3           1 ACT.G   4 + 5.00

                ACG

 a  c  g  t Del Ins Score
 3  0  0  0   0      4.75
                  0 -0.102
 0  2  0  0   1      2.10
                  3 -1.50
 0  0  3  0   0      4.75

Score: 9.5  Columns: 2  Sequences: 2

                **
s1           2 A.   3 + 1.00
s2           5 A.   6 + 1.00

                A

 a  c  g  t Del Ins Score
 2  0  0  0   0      3.00
                  0 -0.102
 0  0  0  0   2     -inf

";

#[test]
fn glam2_alignments_become_gapped_profiles() {
    let dir = tempfile::tempdir().unwrap();
    let src = write(dir.path(), "glam2.txt", GLAM2);
    assert_eq!(Format::detect(&src).unwrap(), Format::Glam2);
    let (out, warnings) = run(dir.path(), std::slice::from_ref(&src));
    assert_eq!(out.alphabet, b"ACGT");
    assert_eq!(out.motifs.len(), 2);
    assert!(warnings.is_empty());
    let (id, _, w, sites, f) = &out.motifs[0];
    assert_eq!((id.as_str(), *w, *sites), ("glam2.1", 3, Some(3.0)));
    assert_eq!(f[1], (FactorKind::Joint, vec![1], vec![0.0, 2.0, 0.0, 0.0]));
    let r = FfmReader::open(std::fs::read(dir.path().join("out.ffm")).unwrap()).unwrap();
    let m = r.motif(0).unwrap();
    let gaps = m.view().gaps().expect("a gapped profile");
    assert_eq!(gaps.deletions(), [0.0, 1.0, 0.0]);
    assert_eq!(gaps.insertions(), [0.0, 3.0]);
    // A column deleted in every sequence keeps its all-zero counts.
    let (_, _, _, _, f) = &out.motifs[1];
    assert_eq!(f[1].2, vec![0.0; 4]);
    let m = r.motif(1).unwrap();
    assert_eq!(m.view().gaps().unwrap().deletions(), [0.0, 2.0]);

    // A single alignment takes the file name; a broken table is an error.
    let one = GLAM2.split("\nScore: 9.5").next().unwrap();
    let src = write(dir.path(), "zf.txt", one);
    let (out, _) = run(dir.path(), &[src]);
    assert_eq!(out.motifs[0].0, "zf");
    let bad = write(
        dir.path(),
        "bad.txt",
        &one.replace("   1      2.10", "      2.10"),
    );
    let dst = dir.path().join("bad.ffm");
    let c = Conversion::between(Format::Glam2, Format::Ffm, bad, &dst).unwrap();
    assert!(convert(c, Budget::new(1, 1 << 30), None).is_err());
}
