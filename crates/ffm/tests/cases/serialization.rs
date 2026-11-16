use std::sync::Arc;
use std::sync::atomic::{AtomicUsize, Ordering};

use omicsweave_core::{ReadAt, cache::BlockCache};
use omicsweave_ffm::{
    Alphabet, Codec, Factor, FactorKind, FfmReader, Key, KeyType, Motif, Schema, Support,
    WriterOptions,
};

use crate::support::{pwm, write};

#[test]
fn independently_encoded_file_has_expected_models_and_annotations() {
    use omicsweave_ffm::{AnnotationValue, FfmWriter, MotifSet};

    let bytes = crate::support::wire::independent_file();
    let reader = FfmReader::open(&bytes).unwrap();
    let expected = Motif {
        width: 2,
        support: Support::All,
        sites: Some(5.0),
        factors: vec![
            Factor {
                kind: FactorKind::Joint,
                positions: vec![0],
                table: vec![1.0, 3.0],
            },
            Factor {
                kind: FactorKind::Conditional,
                positions: vec![0, 1],
                table: vec![2.0, 6.0, 4.0, 4.0],
            },
        ],
        gaps: None,
    };
    assert_eq!(reader.len(), 1);
    assert_eq!(reader.block_count(), 1);
    assert_eq!(reader.minor_version(), 0);
    assert_eq!(reader.id(0), Some("ref"));
    assert_eq!(reader.find("ref"), Some(0));
    assert_eq!(reader.alphabet().symbols(), b"AB");
    assert_eq!(reader.alphabet().complement(), Some(&[1, 0][..]));
    assert_eq!(reader.motif(0).unwrap().view().to_motif(), expected);
    let annotation = reader.annotation(0).unwrap();
    assert_eq!(annotation.values().collect::<Vec<_>>(), ["line\nβ", "-7"]);
    assert_eq!(annotation.value(1), Some(AnnotationValue::Integer(-7)));
    assert_eq!(annotation.get_by_name("unknown"), None);
    let resident = MotifSet::load(&bytes).unwrap();
    assert_eq!(resident.motif(0).unwrap().to_motif(), expected);
    assert_eq!(resident.annotation(0).unwrap().as_str(), "line\nβ\t-7");

    let opts = WriterOptions {
        structure_codec: Codec::Raw,
        values_codec: Codec::Raw,
        index_codec: Codec::Raw,
        annotation_codec: Codec::Raw,
        ids_codec: Codec::Raw,
        ..WriterOptions::default()
    };
    let mut writer = FfmWriter::new(
        Vec::new(),
        reader.alphabet().clone(),
        reader.schema().clone(),
        opts,
    )
    .unwrap();
    writer.push("ref", &expected, &["line\nβ", "-7"]).unwrap();
    assert_eq!(writer.finish().unwrap(), bytes);
}

#[test]
fn oversized_record_stays_in_one_block_and_reader_bounds_are_checked() {
    let rows = vec![
        ("wide".into(), pwm(5, 1), vec![]),
        ("small".into(), pwm(1, 2), vec![]),
    ];
    let bytes = write(
        Alphabet::dna(),
        Schema::default(),
        WriterOptions {
            block_values: 1,
            ..WriterOptions::default()
        },
        &rows,
    );
    let reader = FfmReader::open(bytes).unwrap();
    assert_eq!(reader.block_count(), 2);
    for (i, (_, expected, _)) in rows.iter().enumerate() {
        assert_eq!(reader.motif(i).unwrap().view().to_motif(), *expected);
    }
    for i in [reader.len(), usize::MAX] {
        assert_eq!(reader.id(i), None);
        assert!(reader.motif(i).is_err(), "index {i}");
        assert!(reader.annotation(i).is_err(), "index {i}");
    }
    let mut iter = reader.motifs();
    assert_eq!(iter.size_hint(), (2, Some(2)));
    iter.next().unwrap().unwrap();
    assert_eq!(iter.size_hint(), (1, Some(1)));
    iter.next().unwrap().unwrap();
    assert!(iter.next().is_none());
    assert!(iter.next().is_none());
}

#[test]
fn shared_cache_namespaces_separate_sources_and_detach_preserves_leases() {
    let bytes = |seed| {
        write(
            Alphabet::dna(),
            Schema::default(),
            WriterOptions::default(),
            &[("same-id".into(), pwm(2, seed), vec![])],
        )
    };
    let mut first = FfmReader::open(bytes(1)).unwrap();
    let mut second = FfmReader::open(bytes(2)).unwrap();
    let cache = Arc::new(BlockCache::new(1 << 20, 16));
    first.set_cache(Some(cache.clone()));
    second.set_cache(Some(cache.clone()));
    assert_ne!(first.cache_namespace(), second.cache_namespace());
    assert!(Arc::ptr_eq(first.cache().unwrap(), &cache));
    let held = first.motif(0).unwrap();
    let cloned = held.clone();
    assert_eq!(second.motif(0).unwrap().view().to_motif(), pwm(2, 2));
    cache.clear_source(first.cache_namespace());
    first.set_cache(None);
    assert!(first.cache().is_none());
    drop(first);
    drop(held);
    assert_eq!(cloned.view().to_motif(), pwm(2, 1));
    assert_eq!(second.motif(0).unwrap().view().to_motif(), pwm(2, 2));
    assert_eq!(cache.stats().loads, 2);
}

fn mixed_models() -> Vec<Motif> {
    // Protein motif over D,E,K,R with a rest category: arity 5.
    let protein_pair = Motif {
        width: 3,
        support: Support::Codes(vec![2, 3, 8, 14]),
        sites: None,
        factors: vec![
            Factor {
                kind: FactorKind::Joint,
                positions: vec![0, 2],
                table: (0..25).map(|i| i as f64).collect(),
            },
            Factor {
                kind: FactorKind::Potential,
                positions: vec![1],
                table: vec![0.0, 1.0, 2.0, 0.5, 3.0],
            },
        ],
        gaps: None,
    };
    // Second-order conditional chain over the whole protein alphabet.
    let markov = Motif {
        width: 2,
        support: Support::All,
        sites: Some(7.5),
        factors: vec![Factor {
            kind: FactorKind::Conditional,
            positions: vec![0, 1],
            table: vec![0.05; 400],
        }],
        gaps: None,
    };
    vec![protein_pair, markov]
}

#[test]
fn round_trip_across_blocks_and_codecs() {
    let schema = Schema::new(vec![
        Key::new("name", KeyType::String),
        Key::new("rank", KeyType::Integer),
        Key::new("score", KeyType::Float),
    ])
    .unwrap();
    let motifs: Vec<_> = (0..300)
        .map(|i| {
            (
                format!("MA{i:04}.1"),
                pwm(6 + i % 15, i as u64),
                vec![
                    format!("TF{i}"),
                    if i % 3 == 0 {
                        String::new()
                    } else {
                        i.to_string()
                    },
                    "0.5".into(),
                ],
            )
        })
        .collect();
    for (sc, vc, ac) in [
        (Codec::Pco, Codec::Pco, Codec::Zstd),
        (Codec::Raw, Codec::Raw, Codec::Raw),
    ] {
        let opts = WriterOptions {
            block_values: 500,
            structure_codec: sc,
            values_codec: vc,
            annotation_codec: ac,
            index_codec: sc,
            ids_codec: ac,
            ..WriterOptions::default()
        };
        let bytes = write(Alphabet::dna(), schema.clone(), opts, &motifs);
        let reader = FfmReader::open(bytes).unwrap();
        assert_eq!(reader.len(), 300);
        assert!(reader.block_count() > 10);
        // Out of order access exercises both block reuse and reloads.
        for i in (0..300).rev().step_by(7).chain(0..300) {
            let (id, model, ann) = &motifs[i];
            assert_eq!(reader.id(i), Some(id.as_str()));
            assert_eq!(reader.find(id), Some(i));
            assert_eq!(&reader.motif(i).unwrap().view().to_motif(), model);
            let a = reader.annotation(i).unwrap();
            assert_eq!(
                a.values().collect::<Vec<_>>(),
                ann.iter().map(String::as_str).collect::<Vec<_>>()
            );
            assert_eq!(a.get_by_name("name"), Some(ann[0].as_str()));
        }
        assert_eq!(reader.find("absent"), None);
    }
}

#[test]
fn support_subsets_and_higher_order_factors() {
    let models = mixed_models();
    let motifs: Vec<_> = models
        .iter()
        .enumerate()
        .map(|(i, m)| (format!("m{i}"), m.clone(), vec![]))
        .collect();
    let bytes = write(
        Alphabet::protein(),
        Schema::default(),
        WriterOptions::default(),
        &motifs,
    );
    let reader = FfmReader::open(bytes).unwrap();
    let v0 = reader.motif(0).unwrap();
    let view = v0.view();
    assert_eq!(view.arity(), 5);
    assert_eq!(view.support(), Some(&[2u32, 3, 8, 14][..]));
    let f: Vec<_> = view.factors().collect();
    assert_eq!(f[0].positions, &[0, 2]);
    assert_eq!(f[1].kind, FactorKind::Potential);
    assert_eq!(reader.motif(1).unwrap().view().to_motif(), models[1]);
    assert_eq!(reader.annotation(1).unwrap().values().count(), 0);
}

#[test]
fn shared_cache_serves_the_same_models() {
    let motifs: Vec<_> = (0..50)
        .map(|i| (format!("id{i}"), pwm(10, i), vec![]))
        .collect();
    let opts = WriterOptions {
        block_records: 8,
        ..WriterOptions::default()
    };
    let bytes = write(Alphabet::dna(), Schema::default(), opts, &motifs);
    let mut reader = FfmReader::open(bytes).unwrap();
    let cache = Arc::new(BlockCache::new(1 << 20, 64));
    reader.set_cache(Some(cache.clone()));
    for _ in 0..2 {
        for (i, (_, m, _)) in motifs.iter().enumerate() {
            assert_eq!(&reader.motif(i).unwrap().view().to_motif(), m);
        }
    }
    let stats = cache.stats();
    assert_eq!(stats.loads, 7);
    assert!(stats.hits > 0);
}

#[test]
fn traversal_and_shared_cache_bound_block_reads_and_preserve_held_views() {
    struct Counted {
        bytes: Vec<u8>,
        reads: Arc<AtomicUsize>,
    }
    impl ReadAt for Counted {
        fn len(&self) -> u64 {
            self.bytes.len() as u64
        }
        fn read_exact_at(&self, out: &mut [u8], at: u64) -> omicsweave_core::Result<()> {
            self.reads.fetch_add(1, Ordering::Relaxed);
            self.bytes.read_exact_at(out, at)
        }
    }
    let models: Vec<_> = (0..5)
        .map(|i| (format!("id{i}"), pwm(i + 2, i as u64), vec![]))
        .collect();
    let opts = WriterOptions {
        block_records: 2,
        ..WriterOptions::default()
    };
    let bytes = write(Alphabet::dna(), Schema::default(), opts, &models);
    let reads = Arc::new(AtomicUsize::new(0));
    let mut reader = FfmReader::open(Counted {
        bytes,
        reads: reads.clone(),
    })
    .unwrap();
    reads.store(0, Ordering::Relaxed);
    let held: Vec<_> = reader.motifs().collect::<Result<_, _>>().unwrap();
    assert_eq!(reads.load(Ordering::Relaxed), reader.block_count());
    for (got, (_, expected, _)) in held.iter().zip(&models) {
        let view = got.view();
        assert_eq!(view.to_motif(), *expected);
        assert_eq!(view.pfm().unwrap().values(), view.table_values());
        let mut factors = view.factors();
        let last = factors.nth(view.factor_count() - 1).unwrap();
        assert_eq!(last.positions, &[view.width() - 1]);
        assert!(factors.next().is_none());
        assert!(factors.nth(usize::MAX).is_none());
        assert!(factors.next().is_none());
        let mut factors = view.factors();
        factors.next().unwrap();
        let mut clone = factors.clone();
        assert_eq!(
            factors.nth(view.factor_count() - 2).unwrap().table,
            clone.nth(view.factor_count() - 2).unwrap().table
        );
    }
    reads.store(0, Ordering::Relaxed);
    reader.motif(0).unwrap();
    reader.motif(1).unwrap();
    assert_eq!(reads.load(Ordering::Relaxed), 2);
    let cache = Arc::new(BlockCache::new(1 << 20, 16));
    reader.set_cache(Some(cache.clone()));
    reads.store(0, Ordering::Relaxed);
    let lease = reader.motif(0).unwrap();
    reader.motif(1).unwrap();
    assert_eq!(reads.load(Ordering::Relaxed), 1);
    cache.clear();
    assert_eq!(cache.stats().resident_bytes, 0);
    assert!(cache.stats().live_bytes > 0);
    assert_eq!(lease.view().to_motif(), models[0].1);
    reader.motif(0).unwrap();
    assert_eq!(reads.load(Ordering::Relaxed), 2);
    drop(lease);
    drop(reader);
    cache.clear();
    assert_eq!(cache.stats().live_bytes, 0);
    assert_eq!(held[0].view().to_motif(), models[0].1);
}

#[test]
fn lazy_annotation_retains_its_charged_lease_after_cache_eviction() {
    let text = "注释".repeat(4096);
    let schema = Schema::new(vec![Key::new("description", KeyType::String)]).unwrap();
    let bytes = write(
        Alphabet::dna(),
        schema,
        WriterOptions::default(),
        &[("annotated".into(), pwm(2, 1), vec![text.clone()])],
    );
    let mut reader = FfmReader::open(bytes).unwrap();
    let cache = Arc::new(BlockCache::new(1 << 20, 16));
    reader.set_cache(Some(cache.clone()));
    let model = reader.motif(0).unwrap();
    let charged = cache.stats().live_bytes;
    assert!(charged > text.len());
    let annotation = reader.annotation(0).unwrap();
    assert_eq!(annotation.get(0), Some(text.as_str()));
    assert_eq!(cache.stats().live_bytes, charged);
    cache.clear_source(reader.cache_namespace());
    assert_eq!(cache.stats().resident_bytes, 0);
    assert_eq!(cache.stats().live_bytes, charged);
    drop(model);
    drop(reader);
    assert_eq!(annotation.as_str(), text);
    drop(annotation);
    assert_eq!(cache.stats().live_bytes, 0);
}

#[test]
fn empty_collection() {
    let bytes = write(
        Alphabet::dna(),
        Schema::default(),
        WriterOptions::default(),
        &[],
    );
    let reader = FfmReader::open(bytes).unwrap();
    assert!(reader.is_empty());
    assert_eq!(reader.find("x"), None);
    assert!(reader.motif(0).is_err());
}
