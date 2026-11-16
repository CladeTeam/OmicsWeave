use crate::support::*;
use omicsweave_core::cache::BlockCache;
use omicsweave_core::{AlignmentFields, AlignmentRead, AlignmentsOut, Error, FileSource, Region};
use omicsweave_cram::CramReader;
use omicsweave_fasta::FastaReader;
use std::sync::Arc;

#[test]
fn external_cache_preserves_fields_budgets_clone_and_default_reads() {
    let f = fixture(
        71,
        vec![("chr1".into(), 30_000)],
        SamSpec {
            per_ref: 300,
            gap: 20,
            long_cigars: false,
            unplaced: 0,
        },
    );
    let (cram, crai) = f.cram(
        "cache",
        &["version=3.1", "seqs_per_slice=50", "slices_per_container=3"],
        true,
    );
    let hts = Hts::open(&cram, &crai, Some(&f.fasta));
    let region = Region::new("chr1", 1200, 2600);
    let expected = hts.query(0, 1200, 2600);
    let mut reader = open(&cram, &crai, Some(&f.fasta));
    assert_eq!(
        query(&mut reader, 0, region, AlignmentFields::ALL),
        expected
    );
    reader.reset_stats();
    assert_eq!(
        query(&mut reader, 0, region, AlignmentFields::ALL),
        expected
    );
    assert!(reader.stats().read_calls > 0);
    assert!(reader.stats().slices_decoded > 0);
    let cache = Arc::new(BlockCache::new(64 << 20, 4096));
    reader.set_cache(Some(cache.clone()));
    let namespace = reader.cache_namespace();
    for bits in 0..=AlignmentFields::ALL.bits() {
        let fields = AlignmentFields::from_bits(bits).unwrap();
        let mut keys = Vec::new();
        reader.describe_cache(region, fields, &mut keys).unwrap();
        assert!(keys.iter().all(|k| k.source == namespace && k.kind == 6));
        let before = reader.stats().slices_decoded;
        reader.prefetch_cache(region, fields).unwrap();
        assert!(reader.stats().slices_decoded > before);
        reader.reset_stats();
        assert_eq!(
            query(&mut reader, 0, region, fields),
            only(&expected, fields)
        );
        assert_eq!(reader.stats().read_calls, 0);
        assert_eq!(reader.stats().slices_decoded, 0);
    }
    let mut cloned = reader
        .try_clone_with(FileSource::try_clone, FastaReader::try_clone)
        .unwrap();
    assert_eq!(cloned.cache_namespace(), namespace);
    assert_eq!(
        query(&mut cloned, 0, region, AlignmentFields::ALL),
        expected
    );
    assert_eq!(cloned.stats().read_calls, 0);
    assert_eq!(cloned.stats().slices_decoded, 0);
    let replacement = reader.with_parts(
        FileSource::open(&cram).unwrap(),
        Some(fasta_reader(&f.fasta)),
    );
    assert_ne!(replacement.cache_namespace(), namespace);
    cache.clear_source(namespace);
    cloned.reset_stats();
    assert_eq!(
        query(&mut cloned, 0, region, AlignmentFields::ALL),
        expected
    );
    assert!(cloned.stats().slices_decoded > 0);
    for budget in [0, 1, 1024] {
        let cache = Arc::new(BlockCache::new(budget, 64));
        reader.set_cache(Some(cache.clone()));
        assert_eq!(
            query(&mut reader, 0, region, AlignmentFields::ALL),
            expected
        );
        reader.reset_stats();
        assert_eq!(
            query(&mut reader, 0, region, AlignmentFields::ALL),
            expected
        );
        assert!(reader.stats().slices_decoded > 0);
        assert!(cache.stats().live_bytes <= budget);
    }
}

#[test]
fn invalid_cache_payload_is_rejected_and_staged_output_is_not_ready() {
    let f = fixture(
        72,
        vec![("chr1".into(), 20_000)],
        SamSpec {
            per_ref: 50,
            gap: 20,
            long_cigars: false,
            unplaced: 0,
        },
    );
    let (cram, crai) = f.cram("typed", &["version=3.0"], true);
    let mut reader = open(&cram, &crai, Some(&f.fasta));
    let region = Region::new("chr1", 0, 3000);
    let mut keys = Vec::new();
    reader
        .describe_cache(region, AlignmentFields::ALL, &mut keys)
        .unwrap();
    let cache = Arc::new(BlockCache::new(1 << 20, 64));
    cache.get_or_load(keys[0], || Ok((17u32, 0))).unwrap();
    reader.set_cache(Some(cache.clone()));
    assert!(reader.query(region, AlignmentFields::ALL).is_err());
    assert!(matches!(
        reader.write(AlignmentsOut {
            entries: &mut [],
            names: None,
            seqs: None,
            quals: None,
            cigars: None,
            aux: None,
            modifications: None,
        }),
        Err(Error::NotQueried)
    ));
    cache.clear();
    assert!(reader.query(region, AlignmentFields::ALL).is_ok());
}

#[test]
fn replacement_reference_and_checksum_options_do_not_share_decoded_slices() {
    let f = fixture(
        73,
        vec![("chr1".into(), 20_000)],
        SamSpec {
            per_ref: 50,
            gap: 20,
            long_cigars: false,
            unplaced: 0,
        },
    );
    let (cram, crai) = f.cram("reference", &["version=3.1"], true);
    let region = Region::new("chr1", 0, 3000);
    let cache = Arc::new(BlockCache::new(8 << 20, 256));
    let mut reader = open(&cram, &crai, Some(&f.fasta));
    reader.set_cache(Some(cache.clone()));
    let correct = query(&mut reader, 0, region, AlignmentFields::ALL);
    let mut rng = Rng::new(1729);
    let wrong = f.dir.path().join("wrong.fa");
    write_fasta(&wrong, &f.refs, &reference(&mut rng, &f.refs));
    let mut changed =
        reader.with_parts(FileSource::open(&cram).unwrap(), Some(fasta_reader(&wrong)));
    let error = changed.query(region, AlignmentFields::ALL).unwrap_err();
    assert!(error.to_string().contains("MD5"));
    let mut options = reader.options().clone();
    options.verify_checksums = false;
    let mut permissive = CramReader::open(
        FileSource::open(&cram).unwrap(),
        reader.index().clone(),
        Some(fasta_reader(&wrong)),
        options,
    )
    .unwrap();
    let expected = query(&mut permissive, 0, region, AlignmentFields::ALL);
    assert_ne!(correct, expected);
    permissive.set_cache(Some(cache));
    permissive
        .prefetch_cache(region, AlignmentFields::ALL)
        .unwrap();
    permissive.reset_stats();
    assert_eq!(
        query(&mut permissive, 0, region, AlignmentFields::ALL),
        expected
    );
    assert_eq!(permissive.stats().slices_decoded, 0);
    assert_ne!(permissive.cache_namespace(), reader.cache_namespace());
}
