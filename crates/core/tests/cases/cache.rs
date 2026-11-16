//! Public cache invalidation keeps leases usable while accounting stays live.

use omicsweave_core::cache::{BlockCache, BlockKey, SourceToken};

#[test]
fn source_invalidation_preserves_other_sources_and_held_payloads() {
    let cache = BlockCache::new(16_384, 8);
    let source = SourceToken::new();
    let unrelated = SourceToken::new();
    let key = BlockKey::new(source, 4, [1, 2], 3, 0);
    let other_key = BlockKey {
        source: unrelated,
        ..key
    };
    let held = cache.get_or_load(key, || Ok((vec![1u8, 2], 2))).unwrap();
    let other = cache.get_or_load(other_key, || Ok((vec![9u8], 1))).unwrap();
    let before = cache.stats();

    cache.clear_source(source);
    assert!(cache.get::<Vec<u8>>(key).unwrap().is_none());
    assert_eq!(&**cache.get::<Vec<u8>>(other_key).unwrap().unwrap(), &[9]);
    assert_eq!(&**held, &[1, 2]);
    assert_eq!(cache.stats().resident_entries, 1);
    assert_eq!(cache.stats().live_entries, 2);
    assert_eq!(cache.stats().live_bytes, before.live_bytes);

    let replacement = cache.get_or_load(key, || Ok((vec![3u8], 1))).unwrap();
    assert_eq!(&**replacement, &[3]);
    assert_eq!(&**held, &[1, 2]);
    cache.clear();
    assert!(cache.keys().is_empty());
    assert_eq!(cache.stats().live_entries, 3);
    drop((held, other, replacement));
    assert_eq!(cache.stats().live_entries, 0);
    assert_eq!(cache.stats().live_bytes, 0);
}

#[test]
fn disabling_admission_keeps_results_usable_and_reenabling_allows_caching() {
    let cache = BlockCache::new(0, 0);
    let key = BlockKey::new(SourceToken::new(), 1, [0, 0], 0, 0);
    let direct = cache.get_or_load(key, || Ok((vec![7u8], 1))).unwrap();
    assert!(!direct.is_cached());
    assert_eq!(&**direct, &[7]);
    assert!(cache.get::<Vec<u8>>(key).unwrap().is_none());
    assert_eq!(cache.stats().live_entries, 0);

    cache.resize(4096, 2);
    assert_eq!(cache.limits(), (4096, 2));
    let admitted = cache.get_or_load(key, || Ok((vec![8u8], 1))).unwrap();
    assert!(admitted.is_cached());
    assert_eq!(&**admitted, &[8]);
    assert_eq!(&**direct, &[7]);
    assert_eq!(cache.keys(), vec![key]);
}
