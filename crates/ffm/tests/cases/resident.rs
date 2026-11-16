use omicsweave_ffm::{
    Alphabet, AnnotationValue, Codec, FfmReader, FfmWriter, Motif, MotifSet, MotifSetBuilder,
    Schema, WriterOptions,
};

use crate::support::{pair, pfm, schema};

#[test]
fn entry_annotation_and_factor_iterators_track_remaining_items() {
    let mut builder = MotifSetBuilder::new(Alphabet::dna(), schema());
    for (id, model) in [("pfm", pfm()), ("pair", pair()), ("tail", pfm())] {
        builder.push(id, &model, &["label", "", "0.25"]).unwrap();
    }
    let set = builder.finish();
    let bytes = set.write_to(Vec::new(), WriterOptions::default()).unwrap();
    let reader = FfmReader::open(bytes).unwrap();
    let mut entries = set.iter();
    assert_eq!(entries.len(), 3);
    assert_eq!(entries.next().unwrap().id(), "pfm");
    assert_eq!(entries.next_back().unwrap().index(), 2);
    assert_eq!(entries.len(), 1);
    assert_eq!(entries.next_back().unwrap().id(), "pair");
    assert!(entries.next().is_none());
    assert!(entries.next_back().is_none());
    let annotation = set.annotation(0).unwrap();
    let mut values = annotation.values();
    assert_eq!(values.len(), 3);
    assert_eq!(values.next_back(), Some("0.25"));
    assert_eq!(values.next(), Some("label"));
    assert_eq!(values.next_back(), Some(""));
    assert_eq!(values.len(), 0);
    assert_eq!(values.next(), None);
    assert_eq!(annotation.get_by_name("missing"), None);
    assert_eq!(annotation.value(usize::MAX), None);
    assert!(set.annotation(usize::MAX).is_none());

    // Both PFM and general factors have streamed and resident layouts.
    for i in 0..set.len() {
        let streamed = reader.motif(i).unwrap();
        let expected = set.motif(i).unwrap().to_motif();
        for model in [streamed.view(), set.motif(i).unwrap()] {
            let mut factors = model.factors();
            for expected in &expected.factors {
                assert_eq!(factors.size_hint(), (factors.len(), Some(factors.len())));
                let mut cloned = factors.clone();
                let got = factors.next().unwrap();
                assert_eq!(got.positions, expected.positions);
                assert_eq!(got.table, expected.table);
                assert_eq!(cloned.next().unwrap().table, got.table);
            }
            assert_eq!(factors.len(), 0);
            assert!(factors.next().is_none());
            let mut skipped = model.factors();
            assert!(skipped.nth(usize::MAX).is_none());
            assert_eq!(skipped.len(), 0);
            assert!(skipped.next().is_none());
        }
    }
}

#[test]
fn borrowed_append_rejects_alphabet_mismatch_without_claiming_the_id() {
    let mut source = MotifSetBuilder::new(Alphabet::dna(), Schema::default());
    source.push("source", &pfm(), &[]).unwrap();
    let source = source.finish();
    let view = source.motif(0).unwrap();
    let mut builder = MotifSetBuilder::new(Alphabet::protein(), Schema::default());
    assert!(builder.push_view("new", view, &[]).is_err());
    assert!(builder.is_empty());
    let protein = Motif::pfm(&[vec![1.0; 20]], None);
    builder.push("new", &protein, &[]).unwrap();
    assert_eq!(builder.finish().motif(0).unwrap().to_motif(), protein);
    let mut writer = FfmWriter::new(
        Vec::new(),
        Alphabet::protein(),
        Schema::default(),
        WriterOptions::default(),
    )
    .unwrap();
    assert!(writer.push_view("new", view, &[]).is_err());
    assert!(writer.is_empty());
    writer.push("new", &protein, &[]).unwrap();
    let reader = FfmReader::open(writer.finish().unwrap()).unwrap();
    assert_eq!(reader.motif(0).unwrap().view().to_motif(), protein);
}

#[test]
fn resident_views_share_tables_and_index_fields() {
    let mut builder = MotifSetBuilder::new(Alphabet::dna(), schema());
    builder.push("α", &pfm(), &["TF\n一", "42", "0.5"]).unwrap();
    builder.push("β", &pair(), &["", "", "-0"]).unwrap();
    let set = builder.finish();
    assert_eq!(set.len(), 2);
    assert_eq!(set.factor_count(), 3);
    assert_eq!(set.value_count(), 12);
    let first = set.get_by_id("α").unwrap();
    assert_eq!(first.index(), 0);
    let model = first.model();
    assert_eq!(model.to_motif(), pfm());
    let matrix = model.pfm().unwrap();
    assert_eq!(matrix.width(), 2);
    assert_eq!(matrix.arity(), 4);
    assert_eq!(matrix.row(1), Some(&[4.0, 3.0, 2.0, 1.0][..]));
    assert_eq!(matrix.row(usize::MAX), None);
    assert_eq!(matrix.rows().count(), 2);
    assert_eq!(
        matrix.values().as_ptr(),
        model.factor(0).unwrap().table.as_ptr()
    );
    assert_eq!(matrix.values().as_ptr(), model.table_values().as_ptr());
    assert_eq!(first.annotation().as_str(), "TF\n一\t42\t0.5");
    assert_eq!(
        first.annotation().value(0),
        Some(AnnotationValue::String("TF\n一"))
    );
    assert_eq!(
        first.annotation().value(1),
        Some(AnnotationValue::Integer(42))
    );
    assert_eq!(
        first.annotation().value(2),
        Some(AnnotationValue::Float(0.5))
    );
    assert_eq!(first.annotation().get(usize::MAX), None);
    assert_eq!(set.annotation(1).unwrap().get(0), Some(""));
    assert_eq!(set.annotation(1).unwrap().value(1), None);
    assert_eq!(set.get(1).unwrap().model().to_motif(), pair());
    assert!(set.motif(1).unwrap().pfm().is_none());
    assert!(set.motif(0).unwrap().factor(usize::MAX).is_none());
    assert!(set.get(usize::MAX).is_none());
    assert!(set.get_by_id("missing").is_none());
    let mut factors = model.factors();
    assert_eq!(factors.nth(1).unwrap().positions, &[1]);
    assert_eq!(factors.len(), 0);
    let cloned = set.clone();
    let ptr = model.table_values().as_ptr();
    drop(set);
    assert_eq!(cloned.motif(0).unwrap().table_values().as_ptr(), ptr);
    assert!(cloned.heap_bytes() > cloned.value_count() * 8);
}

#[test]
fn load_and_write_borrowed_models_across_blocks() {
    for codec in [Codec::Raw, Codec::Pco] {
        let mut builder = MotifSetBuilder::new(Alphabet::dna(), schema());
        builder.push("first", &pfm(), &["A", "-1", "0.25"]).unwrap();
        builder.push("second", &pair(), &["", "", ""]).unwrap();
        let expected = builder.finish();
        let opts = WriterOptions {
            block_records: 1,
            structure_codec: codec,
            values_codec: codec,
            index_codec: codec,
            ..WriterOptions::default()
        };
        let bytes = expected.write_to(Vec::new(), opts).unwrap();
        let reader = FfmReader::open(bytes.clone()).unwrap();
        let held = reader.motif(0).unwrap();
        assert_eq!(
            held.view().pfm().unwrap().values(),
            expected.motif(0).unwrap().table_values()
        );
        assert_eq!(held.view().table_values().len(), 8);
        assert!(held.view().factor(usize::MAX).is_none());
        assert_eq!(
            reader.annotation(0).unwrap().value(1),
            Some(AnnotationValue::Integer(-1))
        );
        let loaded = reader.into_set().unwrap();
        assert_eq!(held.view().to_motif(), pfm());
        let reloaded = MotifSet::load(loaded.write_to(Vec::new(), opts).unwrap()).unwrap();
        assert_eq!(bytes, reloaded.write_to(Vec::new(), opts).unwrap());
        for (a, b) in expected.iter().zip(reloaded.iter()) {
            assert_eq!(a.id(), b.id());
            assert_eq!(a.model().to_motif(), b.model().to_motif());
            assert_eq!(a.annotation().as_str(), b.annotation().as_str());
        }
        let mut builder = MotifSetBuilder::new(loaded.alphabet().clone(), loaded.schema().clone());
        for r in loaded.iter().rev() {
            let annotation: Vec<_> = r.annotation().values().collect();
            builder.push_view(r.id(), r.model(), &annotation).unwrap();
        }
        assert_eq!(builder.finish().id(0), Some("second"));
    }
}

#[test]
fn empty_collections_and_missing_annotation() {
    let empty = MotifSetBuilder::new(Alphabet::dna(), Schema::default()).finish();
    let empty = MotifSet::load(
        empty
            .write_to(Vec::new(), WriterOptions::default())
            .unwrap(),
    )
    .unwrap();
    assert!(empty.is_empty());
    assert_eq!(empty.iter().count(), 0);
    let mut builder = MotifSetBuilder::new(Alphabet::dna(), Schema::default());
    builder.push("x", &pfm(), &[]).unwrap();
    assert!(builder.push("x", &pfm(), &[]).is_err());
    assert!(builder.push("bad", &pfm(), &["extra"]).is_err());
    assert_eq!(builder.len(), 1);
    let set = builder.finish();
    let loaded =
        MotifSet::load(set.write_to(Vec::new(), WriterOptions::default()).unwrap()).unwrap();
    assert_eq!(loaded.annotation(0).unwrap().as_str(), "");
    assert_eq!(loaded.annotation(0).unwrap().values().count(), 0);
    assert!(loaded.annotation(0).unwrap().get(0).is_none());
}

#[test]
fn flat_matrix_builder_validates_shape_and_mass() {
    let mut builder = MotifSetBuilder::new(Alphabet::dna(), Schema::default());
    assert!(builder.push_pfm("x", &[], None, &[]).is_err());
    assert!(builder.push_pfm("x", &[1.0; 3], None, &[]).is_err());
    assert!(builder.push_pfm("x", &[0.0; 4], None, &[]).is_err());
    assert!(builder.push_pfm("x", &[f64::NAN; 4], None, &[]).is_err());
    assert!(builder.is_empty());
    builder.push_pfm("x", &[1.0; 8], Some(4.0), &[]).unwrap();
    builder.push_pfm("wide", &[2.0; 20], None, &[]).unwrap();
    builder.push_pfm("narrow", &[3.0; 4], None, &[]).unwrap();
    let set = builder.finish();
    assert_eq!(set.motif(0).unwrap().pfm().unwrap().width(), 2);
    assert_eq!(&set.table_values()[..8], &[1.0; 8]);
    assert_eq!(set.factor_count(), 8);
    assert_eq!(set.motif(1).unwrap().factor(4).unwrap().positions, &[4]);
    assert_eq!(set.motif(2).unwrap().factor(0).unwrap().positions, &[0]);
    fn send_sync<T: Send + Sync>() {}
    send_sync::<MotifSet>();
}

#[test]
fn gapped_profiles_round_trip() {
    use omicsweave_ffm::{FfmReader, FfmWriter, Gaps, MotifSet, WriterOptions};
    let mut m = Motif::pfm(
        &[
            vec![3.0, 0.0, 1.0, 0.0],
            vec![0.0, 4.0, 0.0, 0.0],
            vec![1.0, 1.0, 1.0, 0.0],
        ],
        Some(4.0),
    );
    m.gaps = Some(Gaps {
        deletions: vec![0.0, 0.0, 1.0],
        insertions: vec![2.0, 0.0],
    });
    let mut w = FfmWriter::new(
        Vec::new(),
        Alphabet::dna(),
        Schema::default(),
        WriterOptions::default(),
    )
    .unwrap();
    w.push("g", &m, &[]).unwrap();
    w.push("plain", &pfm(), &[]).unwrap();
    let bytes = w.finish().unwrap();
    let reader = FfmReader::open(bytes.clone()).unwrap();
    let back = reader.motif(0).unwrap();
    assert_eq!(back.view().to_motif(), m);
    let g = back.view().gaps().unwrap();
    assert_eq!(g.deletions(), &[0.0, 0.0, 1.0]);
    assert_eq!(g.insertions(), &[2.0, 0.0]);
    assert!(reader.motif(1).unwrap().view().gaps().is_none());
    let set = MotifSet::load(bytes).unwrap();
    assert_eq!(set.motif(0).unwrap().to_motif(), m);
    assert!(set.motif(1).unwrap().gaps().is_none());

    // A column every sequence deletes may have no residues.
    let mut empty = m.clone();
    empty.factors[2].table = vec![0.0; 4];
    empty.gaps.as_mut().unwrap().deletions[2] = 4.0;
    assert!(empty.validate(4).is_ok());
    let mut w = FfmWriter::new(
        Vec::new(),
        Alphabet::dna(),
        Schema::default(),
        WriterOptions::default(),
    )
    .unwrap();
    w.push("empty", &empty, &[]).unwrap();
    let back = FfmReader::open(w.finish().unwrap()).unwrap();
    assert_eq!(back.motif(0).unwrap().view().to_motif(), empty);
    empty.gaps.as_mut().unwrap().deletions[2] = 0.0;
    assert!(empty.validate(4).is_err());

    // Gaps need one single-position joint factor per column.
    let mut bad = m.clone();
    bad.gaps.as_mut().unwrap().insertions.push(1.0);
    assert!(bad.validate(4).is_err());
    let mut bad = m.clone();
    bad.factors.swap(0, 1);
    assert!(bad.validate(4).is_err());
}
