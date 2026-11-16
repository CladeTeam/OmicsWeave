use omicsweave_ffm::{
    Alphabet, FactorKind, FactorView, FfmReader, FfmWriter, Key, KeyType, MAX_ORDER, MAX_TABLE,
    Motif, MotifSetBuilder, Schema, Support, WriterOptions, table_len,
};

use crate::support::{pair, pfm, pwm};

#[test]
fn alphabet_parts_and_encoding_enforce_code_and_output_bounds() {
    for symbols in [b"".as_slice(), b"AA", &[0; 256]] {
        assert!(
            Alphabet::new("bad", symbols).is_err(),
            "symbols {symbols:?}"
        );
    }
    assert!(Alphabet::new("x".repeat(65536), b"A").is_err());
    let alphabet = Alphabet::new("two", b"AB")
        .unwrap()
        .with_complement(b"BA")
        .unwrap();
    assert!(alphabet.clone().with_alias(b'A', b'B').is_err());
    assert!(alphabet.clone().with_alias(b'a', b'C').is_err());
    assert!(alphabet.clone().with_complement(b"A").is_err());
    assert!(alphabet.clone().with_complement(b"BC").is_err());
    let mut lut = *alphabet.lut();
    lut[b'A' as usize] = 1;
    assert!(Alphabet::from_parts("bad".into(), b"AB".to_vec(), None, lut).is_err());
    lut = *alphabet.lut();
    lut[0] = 2;
    assert!(Alphabet::from_parts("bad".into(), b"AB".to_vec(), None, lut).is_err());
    let symbols: Vec<_> = (0..255).collect();
    let largest = Alphabet::new("all", &symbols).unwrap();
    assert_eq!(largest.len(), 255);
    assert!(!largest.is_empty());
    assert_eq!(largest.code(254), Some(254));
    assert_eq!(largest.code(255), None);
    for preset in [Alphabet::dna(), Alphabet::rna(), Alphabet::protein()] {
        for (code, &symbol) in preset.symbols().iter().enumerate() {
            assert_eq!(preset.code(symbol), Some(code as u8));
            assert_eq!(preset.code(symbol.to_ascii_lowercase()), Some(code as u8));
        }
    }
    for size in [1, 3] {
        let mut out = vec![99; size];
        assert!(alphabet.encode(b"AB", &mut out).is_err());
        assert_eq!(out, vec![99; size]);
    }
    alphabet.encode(&[], &mut []).unwrap();
}

#[test]
fn schemas_preserve_order_and_failed_annotations_leave_both_builders_usable() {
    for name in ["", "a\tb", "a\nb", "a\rb", "a\0b"] {
        assert!(
            Schema::new(vec![Key::new(name, KeyType::String)]).is_err(),
            "name {name:?}"
        );
    }
    assert!(
        Schema::new(vec![
            Key::new("x", KeyType::String),
            Key::new("x", KeyType::Float)
        ])
        .is_err()
    );
    assert!(Schema::new(vec![Key::new("x".repeat(65536), KeyType::String)]).is_err());
    let schema = crate::support::schema();
    assert_eq!(schema.index("p"), Some(2));
    assert_eq!(schema.index("missing"), None);
    assert_eq!(schema.keys()[1], Key::new("n", KeyType::Integer));
    let mut builder = MotifSetBuilder::new(Alphabet::dna(), schema.clone());
    let mut writer = FfmWriter::new(
        Vec::new(),
        Alphabet::dna(),
        schema.clone(),
        WriterOptions::default(),
    )
    .unwrap();
    assert_eq!(builder.schema(), &schema);
    assert_eq!(writer.schema(), &schema);
    assert_eq!(builder.alphabet(), writer.alphabet());
    for values in [
        vec!["short"],
        vec!["tab\tvalue", "1", "0.5"],
        vec!["NUL\0value", "1", "0.5"],
        vec!["ok", "9223372036854775808", "0.5"],
        vec!["ok", "1", "not-float"],
    ] {
        assert!(
            builder.push("retry", &pfm(), &values).is_err(),
            "{values:?}"
        );
        assert!(writer.push("retry", &pfm(), &values).is_err(), "{values:?}");
        assert_eq!(builder.len(), 0);
        assert_eq!(writer.len(), 0);
    }
    for id in ["", "a\tb", "a\nb", "a\rb", "a\0b"] {
        assert!(
            builder.push(id, &pfm(), &["ok", "1", "0.5"]).is_err(),
            "ID {id:?}"
        );
        assert!(
            writer.push(id, &pfm(), &["ok", "1", "0.5"]).is_err(),
            "ID {id:?}"
        );
    }
    // Annotation floats use f64 parsing, independently of model table constraints.
    builder
        .push(
            "retry",
            &pfm(),
            &["text\nline", "-9223372036854775808", "NaN"],
        )
        .unwrap();
    writer
        .push(
            "retry",
            &pfm(),
            &["text\nline", "-9223372036854775808", "NaN"],
        )
        .unwrap();
    let set = builder.finish();
    let reader = FfmReader::open(writer.finish().unwrap()).unwrap();
    assert_eq!(reader.len(), 1);
    assert_eq!(set.annotation(0).unwrap().get(2), Some("NaN"));
    assert_eq!(
        reader.annotation(0).unwrap().get(1),
        Some("-9223372036854775808")
    );
    assert_eq!(
        reader.into_set().unwrap().annotation(0).unwrap().as_str(),
        set.annotation(0).unwrap().as_str()
    );
}

#[test]
fn writer_options_reject_invalid_stream_codecs_and_limits() {
    let default = WriterOptions::default();
    let bad = [
        WriterOptions {
            block_values: 0,
            ..default
        },
        WriterOptions {
            block_records: 0,
            ..default
        },
        WriterOptions {
            structure_codec: omicsweave_ffm::Codec::Zstd,
            ..default
        },
        WriterOptions {
            values_codec: omicsweave_ffm::Codec::Zstd,
            ..default
        },
        WriterOptions {
            index_codec: omicsweave_ffm::Codec::Zstd,
            ..default
        },
        WriterOptions {
            annotation_codec: omicsweave_ffm::Codec::Pco,
            ..default
        },
        WriterOptions {
            ids_codec: omicsweave_ffm::Codec::Pco,
            ..default
        },
        WriterOptions {
            levels: omicsweave_ffm::Levels {
                pco: 13,
                ..default.levels
            },
            ..default
        },
        WriterOptions {
            levels: omicsweave_ffm::Levels {
                zstd: i32::MAX,
                ..default.levels
            },
            ..default
        },
    ];
    for opts in bad {
        assert!(
            FfmWriter::new(Vec::new(), Alphabet::dna(), Schema::default(), opts).is_err(),
            "{opts:?}"
        );
    }
}

#[test]
fn model_validation_and_table_limits_reject_invalid_shapes_without_large_tables() {
    assert_eq!(table_len(4, 12).unwrap(), MAX_TABLE);
    assert!(table_len(4, 13).is_err());
    assert_eq!(table_len(1, MAX_ORDER).unwrap(), 1);
    assert!(table_len(1, MAX_ORDER + 1).is_err());
    assert!(table_len(255, 16).is_err());
    for support in [vec![], vec![1, 0], vec![0, 0], vec![4]] {
        let mut model = pfm();
        model.support = Support::Codes(support);
        assert!(model.validate(4).is_err(), "support {:?}", model.support);
    }
    for sites in [0.0, -1.0, f64::NAN, f64::INFINITY] {
        let mut model = pfm();
        model.sites = Some(sites);
        assert!(model.validate(4).is_err(), "sites {sites}");
    }
    let mut model = pfm();
    model.width = 0;
    assert!(model.validate(4).is_err());
    model = pfm();
    model.factors.clear();
    assert!(model.validate(4).is_err());
    for positions in [vec![], vec![0, 0], vec![1, 0], vec![2]] {
        let mut model = pfm();
        model.factors[0].positions = positions;
        assert!(
            model.validate(4).is_err(),
            "positions {:?}",
            model.factors[0].positions
        );
    }
    for value in [-1.0, f64::NAN, f64::INFINITY] {
        let mut model = pfm();
        model.factors[0].table[0] = value;
        assert!(model.validate(4).is_err(), "table value {value}");
    }
    let mut conditional = pair();
    conditional.factors[0].table = vec![1.0, 1.0, 0.0, 0.0];
    assert!(conditional.validate(4).is_err());
    conditional.factors[0].kind = FactorKind::Potential;
    conditional.factors[0].table.fill(0.0);
    conditional.validate(4).unwrap();
}

#[test]
fn normalization_errors_do_not_partially_write_the_output() {
    for table in [
        &[1.0, 1.0, 0.0, 0.0][..],
        &[1.0, -1.0, 1.0, 1.0],
        &[1.0, f64::NAN, 1.0, 1.0],
        &[1.0, f64::INFINITY, 1.0, 1.0],
        &[1.0, 1.0, 1.0],
    ] {
        let factor = FactorView {
            kind: FactorKind::Conditional,
            positions: &[0, 1],
            table,
        };
        let mut out = vec![73.0; table.len()];
        assert!(factor.normalize(2, &mut out).is_err(), "{table:?}");
        assert_eq!(out, vec![73.0; table.len()]);
    }
}

#[test]
fn invalid_models_are_rejected() {
    let mut w = FfmWriter::new(
        Vec::new(),
        Alphabet::dna(),
        Schema::default(),
        WriterOptions::default(),
    )
    .unwrap();
    let mut m = pwm(4, 1);
    m.factors[0].table[1] = f64::NAN;
    assert!(w.push("a", &m, &[]).is_err());
    let mut m = pwm(4, 1);
    m.factors[1].table = vec![0.0; 4];
    assert!(w.push("a", &m, &[]).is_err(), "joint group summing to 0");
    let mut m = pwm(4, 1);
    m.factors[2].positions = vec![4];
    assert!(w.push("a", &m, &[]).is_err(), "position past the width");
    let mut m = pwm(4, 1);
    m.support = Support::Codes(vec![0, 1, 2, 3]);
    assert!(w.push("a", &m, &[]).is_err(), "a full support must be All");
    w.push("a", &pwm(4, 1), &[]).unwrap();
    assert!(w.push("a", &pwm(4, 2), &[]).is_err(), "duplicate ID");
    assert!(w.push("b\tc", &pwm(4, 2), &[]).is_err());
}

#[test]
fn alphabets_map_bytes_and_complements() {
    let dna = Alphabet::dna();
    let mut out = [0u8; 6];
    dna.encode(b"ACgtNx", &mut out).unwrap();
    assert_eq!(out, [0, 1, 2, 3, 0xFF, 0xFF]);
    assert_eq!(dna.complement(), Some(&[3u8, 2, 1, 0][..]));
    assert!(
        Alphabet::new("x", b"AB")
            .unwrap()
            .with_complement(b"AA")
            .is_err()
    );
    let rna_like = Alphabet::new("x", b"ACGT")
        .unwrap()
        .with_alias(b'U', b'T')
        .unwrap();
    assert_eq!(rna_like.code(b'U'), Some(3));
}

#[test]
fn reverse_encoding_complements_before_collapsing_rest() {
    let mut builder = MotifSetBuilder::new(Alphabet::dna(), Schema::default());
    builder.push("subset", &pair(), &[]).unwrap();
    let set = builder.finish();
    let encoding = set.motif(0).unwrap().encoding(set.alphabet()).unwrap();
    assert_eq!(encoding.forward()[b'A' as usize], 0);
    assert_eq!(encoding.forward()[b'T' as usize], 1);
    assert_eq!(encoding.forward()[b'C' as usize], 1);
    let reverse = encoding.reverse().unwrap();
    assert_eq!(reverse[b'T' as usize], 0);
    assert_eq!(reverse[b'C' as usize], 1);
    assert_eq!(reverse[b't' as usize], 0);
    assert_eq!(reverse[b'N' as usize], 255);
    assert_eq!(encoding.forward()[b'N' as usize], 255);
    assert!(
        set.motif(0)
            .unwrap()
            .encoding(&Alphabet::protein())
            .is_err()
    );
    let mut protein = MotifSetBuilder::new(Alphabet::protein(), Schema::default());
    protein
        .push("p", &Motif::pfm(&[vec![1.0; 20]], None), &[])
        .unwrap();
    let protein = protein.finish();
    assert!(
        protein
            .motif(0)
            .unwrap()
            .encoding(protein.alphabet())
            .unwrap()
            .reverse()
            .is_none()
    );
}

#[test]
fn normalization_handles_extreme_counts_and_zero_contexts() {
    let positions = [0];
    let table = [f64::MAX; 4];
    let joint = FactorView {
        kind: FactorKind::Joint,
        positions: &positions,
        table: &table,
    };
    let mut out = [0.0; 4];
    joint.normalize(4, &mut out).unwrap();
    assert_eq!(out, [0.25; 4]);
    assert!(pfm().validate(0).is_err());
    assert!(Motif::pfm(&[vec![]], None).validate(0).is_err());
    assert!(Motif::pfm(&[vec![1.0; 256]], None).validate(256).is_err());
    assert!(table_len(0, 1).is_err());
    assert!(table_len(4, 0).is_err());
    assert!(table_len(1, usize::MAX).is_err());
    let conditional = FactorView {
        kind: FactorKind::Conditional,
        positions: &[0, 1],
        table: &[f64::MAX, f64::MAX, 0.0, 2.0],
    };
    conditional.normalize(2, &mut out).unwrap();
    assert_eq!(out, [0.5, 0.5, 0.0, 1.0]);
    let zero = FactorView {
        table: &[0.0; 4],
        ..conditional
    };
    assert!(zero.normalize(2, &mut out).is_err());
    assert_eq!(out, [0.5, 0.5, 0.0, 1.0]);
    assert!(
        FactorView {
            kind: FactorKind::Potential,
            ..joint
        }
        .normalize(4, &mut out)
        .is_err()
    );
    assert!(joint.normalize(4, &mut [0.0; 3]).is_err());
}
