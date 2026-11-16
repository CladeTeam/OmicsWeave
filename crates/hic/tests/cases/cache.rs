use crate::support::io::Counted;
use omicsweave_core::cache::BlockCache;
use omicsweave_core::{
    ContactEntry, ContactQuery, ContactRead, ContactUnit, Error, ReadAt, Region,
};
use omicsweave_hic::{Chrom, HicCacheLimits, HicHeader, HicReader, HicWriter, Record};
use std::cell::RefCell;
use std::io::Cursor;
use std::rc::Rc;
use std::sync::Arc;

fn fixture(version: i32, value: f32) -> Vec<u8> {
    let mut header = HicHeader::new(
        "g",
        vec![
            Chrom {
                name: "a".into(),
                len: 256,
            },
            Chrom {
                name: "b".into(),
                len: 256,
            },
        ],
        vec![1, 2],
    );
    header.version = version;
    header.frag_resolutions = vec![1];
    header.sites = vec![(0..256).collect(), (0..256).collect()];
    let mut writer = HicWriter::create(Cursor::new(Vec::new()), header, 1).unwrap();
    for (unit, resolution) in [
        (ContactUnit::Bp, 1),
        (ContactUnit::Bp, 2),
        (ContactUnit::Frag, 1),
    ] {
        for (first, second) in [(0, 0), (0, 1)] {
            let records: Vec<_> = (0..128)
                .flat_map(|y| {
                    (0..128)
                        .filter(move |&x| first != second || x <= y)
                        .map(move |x| Record {
                            x,
                            y,
                            value: value + ((x * 3 + y * 7) % 17) as f32 / 8.0,
                        })
                })
                .collect();
            writer
                .write_matrix(first, second, unit, resolution, &records)
                .unwrap();
        }
    }
    writer.finish().unwrap().into_inner()
}

fn query<'a>(first: &'a str, second: &'a str, lo: u32, hi: u32) -> ContactQuery<'a> {
    ContactQuery::new(Region::new(first, lo, hi), Region::new(second, lo, hi), 1)
}

fn matrix(reader: &mut HicReader<Counted>, query: ContactQuery<'_>) -> Vec<f32> {
    let (rows, cols) = query.shape();
    let mut output = vec![0.0; rows * cols];
    reader.matrix(query, &mut output).unwrap();
    output
}

fn pixels(reader: &mut HicReader<Counted>, query: ContactQuery<'_>) -> Vec<ContactEntry> {
    let count = reader.query_entries(query).unwrap();
    let mut output = vec![ContactEntry::default(); count];
    reader.write_entries(&mut output).unwrap();
    output
}

#[test]
fn prepared_blocks_share_across_shapes_contexts_and_output_types_without_reads() {
    for version in [8, 9] {
        let source = Counted::new(fixture(version, 0.5));
        let cache = Arc::new(BlockCache::new(64 << 20, 1000));
        let mut uncached = HicReader::open(source.clone()).unwrap();
        let mut cached = HicReader::open(source.clone())
            .unwrap()
            .with_cache(cache.clone());
        for (first, second) in [("a", "a"), ("a", "b"), ("b", "a")] {
            matrix(&mut cached, query(first, second, 0, 128));
            let narrow = query(first, second, 37, 42);
            let expected = matrix(&mut uncached, narrow);
            let expected_pixels = pixels(&mut uncached, narrow);
            let before = source.reads();
            assert_eq!(matrix(&mut cached, narrow), expected);
            assert_eq!(pixels(&mut cached, narrow), expected_pixels);
            assert_eq!(
                source.reads(),
                before,
                "hot queries must not read even compressed bytes"
            );
            assert!(
                cached.cached_records_visited() <= 50,
                "index must skip unrelated records"
            );
            let mut clone = cached.clone_context();
            assert_eq!(clone.cache_namespace(), cached.cache_namespace());
            assert_eq!(matrix(&mut clone, narrow), expected);
            assert_eq!(
                source.reads(),
                before,
                "context clone shares source session, metadata and cache"
            );
        }
        for q in [
            query("a", "b", 12, 20).with_unit(ContactUnit::Frag),
            ContactQuery::new(Region::new("a", 12, 20), Region::new("b", 12, 20), 2),
        ] {
            assert_eq!(matrix(&mut cached, q), matrix(&mut uncached, q));
        }
        let kinds: std::collections::HashSet<_> =
            cache.keys().iter().map(|k| k.stream[1]).collect();
        assert_eq!(
            kinds.len(),
            3,
            "unit and resolution partition cached payloads"
        );
    }
}

#[test]
fn external_plans_eviction_source_isolation_and_budget_bypass_preserve_results() {
    let source = Counted::new(fixture(9, 1.0));
    let cache = Arc::new(BlockCache::new(32 << 20, 1000));
    let mut reader = HicReader::open(source.clone())
        .unwrap()
        .with_cache(cache.clone());
    let q = query("a", "b", 10, 15);
    let mut requests = Vec::new();
    reader.cache_plan(q, &mut requests).unwrap();
    assert!(!requests.is_empty());
    let leases: Vec<_> = requests
        .iter()
        .map(|r| reader.load_cached_block(r).unwrap().unwrap())
        .collect();
    let before = source.reads();
    let expected = matrix(&mut reader, q);
    assert_eq!(
        source.reads(),
        before,
        "explicit prefetch must satisfy the query"
    );
    let pinned = leases[0].len();
    cache.clear_source(reader.cache_namespace());
    assert_eq!(
        leases[0].len(),
        pinned,
        "eviction cannot invalidate a live read lease"
    );
    assert_eq!(matrix(&mut reader, q), expected);
    assert!(source.reads() > before);
    let mut other = HicReader::open(Counted::new(fixture(9, 7.0)))
        .unwrap()
        .with_cache(cache.clone());
    assert_ne!(other.cache_namespace(), reader.cache_namespace());
    assert!(other.load_cached_block(&requests[0]).is_err());
    assert_ne!(matrix(&mut other, q), expected);
    let tiny = Arc::new(BlockCache::new(1, 1));
    reader.set_cache(Some(tiny.clone()));
    let prepared = reader.cache_preparations();
    assert_eq!(matrix(&mut reader, q), expected);
    assert_eq!(matrix(&mut reader, q), expected);
    assert_eq!(tiny.stats().resident_entries, 0);
    assert_eq!(reader.cache_preparations(), prepared);
}

#[test]
fn zero_and_small_cache_budgets_skip_preparation_and_repeated_inflation() {
    for version in [8, 9] {
        let source = Counted::new(fixture(version, 1.0));
        let q = query("a", "b", 10, 15);
        let mut uncached = HicReader::open(source.clone()).unwrap();
        let expected = matrix(&mut uncached, q);
        for (bytes, entries) in [(0, 1000), (1 << 20, 0), (1, 1), (4096, 1)] {
            let cache = Arc::new(BlockCache::new(bytes, entries));
            let mut reader = HicReader::open(source.clone())
                .unwrap()
                .with_cache(cache.clone());
            let mut requests = Vec::new();
            reader.cache_plan(q, &mut requests).unwrap();
            assert!(!requests.is_empty());
            for _ in 0..2 {
                let before = reader.block_inflations();
                assert_eq!(matrix(&mut reader, q), expected);
                assert_eq!(
                    reader.block_inflations() - before,
                    requests.len() as u64,
                    "a rejected preparation reuses the inflated block"
                );
            }
            assert_eq!(
                reader.cache_preparations(),
                0,
                "an unusable budget must never build or sort prepared indices"
            );
            assert_eq!(cache.stats().resident_entries, 0);
            reader
                .set_cache_limits(HicCacheLimits {
                    max_inflated_bytes: 1,
                    max_records: 1,
                })
                .unwrap();
            assert!(
                reader.matrix(q, &mut [0.0; 25]).is_err(),
                "zero/low cache capacity cannot disable explicit resource limits"
            );
        }
    }
}

#[test]
fn existing_small_block_hits_survive_a_low_nonzero_budget() {
    let mut writer = HicWriter::create(
        Cursor::new(Vec::new()),
        HicHeader::new(
            "g",
            vec![Chrom {
                name: "a".into(),
                len: 20,
            }],
            vec![1],
        ),
        1,
    )
    .unwrap();
    writer
        .write_matrix(
            0,
            0,
            ContactUnit::Bp,
            1,
            &[Record {
                x: 1,
                y: 2,
                value: 1.0,
            }],
        )
        .unwrap();
    let source = Counted::new(writer.finish().unwrap().into_inner());
    let cache = Arc::new(BlockCache::new(1 << 20, 10));
    let mut reader = HicReader::open(source.clone())
        .unwrap()
        .with_cache(cache.clone());
    let q = query("a", "a", 0, 4);
    let expected = matrix(&mut reader, q);
    let charge = cache.stats().resident_bytes;
    cache.resize(charge, 1);
    let reads = source.reads();
    let preparations = reader.cache_preparations();
    assert_eq!(matrix(&mut reader, q), expected);
    assert_eq!(source.reads(), reads);
    assert_eq!(reader.cache_preparations(), preparations);
    assert_eq!(cache.stats().resident_entries, 1);
}

#[test]
fn explicit_decode_limits_fail_cleanly_and_can_be_retried() {
    let source = Counted::new(fixture(8, 1.0));
    let cache = Arc::new(BlockCache::new(32 << 20, 1000));
    let mut reader = HicReader::open(source).unwrap().with_cache(cache);
    let q = query("a", "a", 1, 4);
    reader
        .set_cache_limits(HicCacheLimits {
            max_inflated_bytes: 1,
            max_records: 1,
        })
        .unwrap();
    assert!(reader.matrix(q, &mut [0.0; 9]).is_err());
    reader.set_cache_limits(HicCacheLimits::default()).unwrap();
    assert_eq!(matrix(&mut reader, q).len(), 9);
    reader
        .set_cache_limits(HicCacheLimits {
            max_inflated_bytes: 1,
            max_records: 1,
        })
        .unwrap();
    assert!(
        reader.matrix(q, &mut [0.0; 9]).is_err(),
        "limits apply to existing leases too"
    );
}

#[test]
fn cache_preparation_limits_do_not_reject_larger_valid_blocks() {
    let source = Counted::new(fixture(9, 1.0));
    let cache = Arc::new(BlockCache::new(32 << 20, 1000));
    let mut expected_reader = HicReader::open(source.clone()).unwrap();
    let q = query("a", "a", 2, 7);
    let expected = matrix(&mut expected_reader, q);
    for limits in [
        HicCacheLimits {
            max_inflated_bytes: 1,
            max_records: 1,
        },
        HicCacheLimits {
            max_inflated_bytes: 1 << 20,
            max_records: 1,
        },
    ] {
        let mut reader = HicReader::open(source.clone())
            .unwrap()
            .with_cache(cache.clone());
        reader.set_cache_admission_limits(limits).unwrap();
        let mut requests = Vec::new();
        reader.cache_plan(q, &mut requests).unwrap();
        assert!(reader.load_cached_block(&requests[0]).unwrap().is_none());
        let inflations = reader.block_inflations();
        assert_eq!(matrix(&mut reader, q), expected);
        assert_eq!(
            reader.block_inflations() - inflations,
            requests.len() as u64,
            "admission fallback inflates each dependency once"
        );
        assert_eq!(cache.stats().resident_entries, 0);
    }
}

#[test]
fn strict_limits_return_the_first_decode_error_without_retrying() {
    for limits in [
        HicCacheLimits {
            max_inflated_bytes: 1,
            max_records: 1,
        },
        HicCacheLimits {
            max_inflated_bytes: 1 << 20,
            max_records: 1,
        },
    ] {
        let source = Counted::new(fixture(9, 1.0));
        let mut reader = HicReader::open(source)
            .unwrap()
            .with_cache(Arc::new(BlockCache::new(32 << 20, 1000)));
        reader.set_cache_limits(limits).unwrap();
        let before = reader.block_inflations();
        assert!(
            reader
                .matrix(query("a", "a", 2, 7), &mut [0.0; 25])
                .is_err()
        );
        assert_eq!(reader.block_inflations() - before, 1);
    }
}
#[derive(Clone)]
struct EditableSource(Rc<RefCell<Vec<u8>>>);

impl ReadAt for EditableSource {
    fn len(&self) -> u64 {
        self.0.borrow().len() as u64
    }

    fn read_exact_at(&self, buf: &mut [u8], offset: u64) -> omicsweave_core::Result<()> {
        self.0.borrow().as_slice().read_exact_at(buf, offset)
    }
}

#[test]
fn conflicting_matrix_identity_invalidates_entries_and_allows_retry() {
    let header = HicHeader::new(
        "g",
        ["a", "b"]
            .into_iter()
            .map(|name| Chrom {
                name: name.into(),
                len: 16,
            })
            .collect(),
        vec![1],
    );
    let mut writer = HicWriter::create(Cursor::new(Vec::new()), header, 1).unwrap();
    for (c2, x, y, value) in [(0, 1, 1, 5.0), (1, 2, 3, 7.0)] {
        writer
            .write_matrix(0, c2, ContactUnit::Bp, 1, &[Record { x, y, value }])
            .unwrap();
    }
    let original = writer.finish().unwrap().into_inner();
    let own = ContactQuery::new(Region::new("a", 0, 8), Region::new("a", 0, 8), 1);
    let inter = ContactQuery::new(Region::new("a", 0, 8), Region::new("b", 0, 8), 1);
    for bad_pair in [(0i32, 0i32), (1, 0), (0, 77)] {
        let source = EditableSource(Rc::new(RefCell::new(original.clone())));
        let mut reader = HicReader::open(source.clone()).unwrap();
        assert_eq!(reader.query_entries(own).unwrap(), 1);
        let matrix = reader.footer().matrices[&(0, 1)] as usize;
        {
            let mut bytes = source.0.borrow_mut();
            bytes[matrix..matrix + 4].copy_from_slice(&bad_pair.0.to_le_bytes());
            bytes[matrix + 4..matrix + 8].copy_from_slice(&bad_pair.1.to_le_bytes());
        }
        assert!(matches!(
            reader.query_entries(inter),
            Err(Error::Invalid { .. })
        ));
        assert!(matches!(
            reader.write_entries(&mut [ContactEntry::default()]),
            Err(Error::NotQueried)
        ));
        source.0.borrow_mut()[matrix..matrix + 8].copy_from_slice(&original[matrix..matrix + 8]);
        assert_eq!(reader.query_entries(inter).unwrap(), 1);
        let mut entries = [ContactEntry::default()];
        reader.write_entries(&mut entries).unwrap();
        assert_eq!(
            entries[0],
            ContactEntry {
                row: 2,
                col: 3,
                value: 7.0
            }
        );
    }
}

#[test]
fn adjacent_compressed_blocks_share_one_positional_read() {
    let header = HicHeader::new(
        "g",
        vec![Chrom {
            name: "s".into(),
            len: 2_000_000,
        }],
        vec![1000],
    );
    let mut writer = HicWriter::create(Cursor::new(Vec::new()), header, 1).unwrap();
    let records: Vec<_> = (0..2000u32)
        .map(|x| Record {
            x,
            y: x,
            value: 1.0,
        })
        .collect();
    writer
        .write_matrix(0, 0, ContactUnit::Bp, 1000, &records)
        .unwrap();
    let bytes = writer.finish().unwrap().into_inner();
    let source = Counted::new(bytes);
    let mut reader = HicReader::open(source.clone()).unwrap();
    let blocks = reader
        .zoom(0, 0, ContactUnit::Bp, 1000)
        .unwrap()
        .unwrap()
        .blocks
        .len();
    assert!(blocks > 1);
    source.reset();
    let region = Region::new("s", 0, 2_000_000);
    let n = reader
        .query_entries(ContactQuery::new(region, region, 1000))
        .unwrap();
    assert_eq!(n, 2000);
    assert_eq!(source.reads(), 1);
}
