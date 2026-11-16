use crate::support::{fixture, open, query_all};
use omicsweave_core::{
    AlignmentFields, AlignmentRead, AlignmentsOut, Error, Region, cache::BlockCache,
};
use std::sync::Arc;

#[test]
fn external_blocks_prefetch_clone_and_budget_fallback_preserve_records() {
    let f = fixture(73, 4_000, 5, 6);
    let region = Region::new("chr1", 0, 100_000);
    let mut reader = open(&f.bam, &f.bai);
    let expected = query_all(&mut reader, region);
    assert!(!expected.is_empty());
    reader.reset_stats();
    assert_eq!(query_all(&mut reader, region), expected);
    assert!(
        reader.stats().inflate_calls > 0,
        "default queries do not retain decoded blocks"
    );

    let cache = Arc::new(BlockCache::new(64 << 20, 8192));
    reader.set_cache(Some(Arc::clone(&cache)));
    let mut keys = Vec::new();
    reader.describe_cache(region, &mut keys).unwrap();
    assert!(!keys.is_empty());
    assert!(
        keys.iter()
            .all(|key| key.source == reader.cache_namespace() && key.kind == 5)
    );
    assert!(keys.windows(2).all(|pair| pair[0] != pair[1]));
    reader.prefetch_cache(region).unwrap();
    assert!(cache.stats().loads > 0);
    reader.reset_stats();
    assert_eq!(query_all(&mut reader, region), expected);
    assert_eq!(reader.stats().read_calls, 0);
    assert_eq!(reader.stats().inflate_calls, 0);
    let mut clone = reader.try_clone().unwrap();
    assert!(std::ptr::eq(reader.header(), clone.header()));
    assert!(std::ptr::eq(reader.index(), clone.index()));
    assert_eq!(reader.cache_namespace(), clone.cache_namespace());
    drop(reader);
    assert_eq!(query_all(&mut clone, region), expected);
    assert_eq!(clone.stats().read_calls, 0);
    assert_eq!(clone.stats().inflate_calls, 0);
    cache.clear_source(clone.cache_namespace());
    clone.reset_stats();
    assert_eq!(query_all(&mut clone, region), expected);
    assert!(clone.stats().inflate_calls > 0);
    for (bytes, entries) in [(0, 0), (1, 1)] {
        let limited = Arc::new(BlockCache::new(bytes, entries));
        clone.set_cache(Some(Arc::clone(&limited)));
        clone.reset_stats();
        assert_eq!(query_all(&mut clone, region), expected);
        assert!(clone.stats().inflate_calls > 0);
        assert_eq!(limited.stats().resident_entries, 0);
    }
}

#[test]
fn wrong_cached_payload_fails_safely_and_can_be_evicted() {
    let f = fixture(77, 500, 30, 6);
    let mut reader = open(&f.bam, &f.bai);
    let region = Region::new("chr1", 0, 100_000);
    let expected = query_all(&mut reader, region);
    let cache = Arc::new(BlockCache::new(64 << 20, 8192));
    reader.set_cache(Some(Arc::clone(&cache)));
    let mut keys = Vec::new();
    reader.describe_cache(region, &mut keys).unwrap();
    assert!(!keys.is_empty());
    cache.get_or_load(keys[0], || Ok((7u64, 0))).unwrap();
    assert!(reader.query(region, AlignmentFields::ALL).is_err());
    let mut entries = Vec::new();
    assert!(matches!(
        reader.write(AlignmentsOut {
            entries: &mut entries,
            names: None,
            seqs: None,
            quals: None,
            cigars: None,
            aux: None,
            modifications: None,
        }),
        Err(Error::NotQueried)
    ));
    assert!(cache.evict(keys[0]));
    assert_eq!(query_all(&mut reader, region), expected);
}

#[test]
fn independent_sources_do_not_share_namespaces_and_clones_do_not_share_staged_results() {
    let f = fixture(85, 24, 30, 0);
    let region = Region::new("chr1", 0, 3_000_000);
    let mut original = open(&f.bam, &f.bai);
    let expected = query_all(&mut original, region);
    let cache = Arc::new(BlockCache::new(64 << 20, 8192));
    original.set_cache(Some(Arc::clone(&cache)));
    original.prefetch_cache(region).unwrap();
    let mut cloned = original.try_clone().unwrap();
    assert!(matches!(
        cloned.write(AlignmentsOut {
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
    assert_eq!(
        cloned
            .query(Region::new("empty", 0, 50_000), AlignmentFields::NONE)
            .unwrap()
            .records,
        0
    );
    // Query scratch in another execution context cannot replace the original staged records.
    let mut entries = vec![Default::default(); expected.len()];
    let size = original.query(region, AlignmentFields::NONE).unwrap();
    assert_eq!(size.records, expected.len());
    cloned
        .query(Region::new("empty", 0, 50_000), AlignmentFields::NONE)
        .unwrap();
    original
        .write(AlignmentsOut {
            entries: &mut entries,
            names: None,
            seqs: None,
            quals: None,
            cigars: None,
            aux: None,
            modifications: None,
        })
        .unwrap();
    for (got, want) in entries.iter().zip(&expected) {
        assert_eq!(
            (got.start, got.end, got.flag),
            (want.pos, want.end, want.flag)
        );
    }
    let mut independent = open(&f.bam, &f.bai);
    independent.set_cache(Some(Arc::clone(&cache)));
    assert_ne!(original.cache_namespace(), independent.cache_namespace());
    let (mut shared_keys, mut own_keys) = (Vec::new(), Vec::new());
    original.describe_cache(region, &mut shared_keys).unwrap();
    independent.describe_cache(region, &mut own_keys).unwrap();
    assert!(!own_keys.is_empty());
    assert!(own_keys.iter().all(|key| !shared_keys.contains(key)));
    independent.reset_stats();
    assert_eq!(query_all(&mut independent, region), expected);
    assert!(independent.stats().inflate_calls > 0);
    cache.clear_source(original.cache_namespace());
    independent.reset_stats();
    assert_eq!(query_all(&mut independent, region), expected);
    assert_eq!(independent.stats().inflate_calls, 0);
}
