use std::io::Cursor;
use std::sync::Arc;
use std::sync::atomic::{AtomicBool, AtomicUsize, Ordering};

use crate::support::bigbed::{Record, Spec, write};
use omicsweave_bbi::{BbiReader, BigWigWriter, Chrom, WriteOptions};
use omicsweave_core::cache::BlockCache;
use omicsweave_core::{
    Error, FileSource, IntervalEntry, IntervalRead, IntervalsOut, Locus, ReadAt, Region,
    StringColumn, Strings, TrackEntry, TrackRead,
};

struct Counted {
    bytes: Vec<u8>,
    reads: Arc<AtomicUsize>,
    fail: Arc<AtomicBool>,
}

impl ReadAt for Counted {
    fn len(&self) -> u64 {
        self.bytes.len() as u64
    }
    fn read_exact_at(&self, out: &mut [u8], at: u64) -> omicsweave_core::Result<()> {
        self.reads.fetch_add(1, Ordering::Relaxed);
        if self.fail.swap(false, Ordering::Relaxed) {
            return Err(Error::invalid("counted source", "injected read failure"));
        }
        self.bytes.read_exact_at(out, at)
    }
}

fn source(bytes: Vec<u8>) -> (Counted, Arc<AtomicUsize>, Arc<AtomicBool>) {
    let reads = Arc::new(AtomicUsize::new(0));
    let fail = Arc::new(AtomicBool::new(false));
    (
        Counted {
            bytes,
            reads: reads.clone(),
            fail: fail.clone(),
        },
        reads,
        fail,
    )
}

fn wig() -> Vec<u8> {
    let options = WriteOptions {
        items_per_slot: Some(4),
        zoom_levels: 0,
        ..Default::default()
    };
    let mut writer = BigWigWriter::create(
        Cursor::new(Vec::new()),
        vec![Chrom {
            name: "c".into(),
            len: 1024,
        }],
        options,
    )
    .unwrap();
    for i in 0..64 {
        writer.add(0, i * 4, i * 4 + 3, i as f32).unwrap();
    }
    writer.finish().unwrap().0.into_inner()
}

fn assert_wig(reader: &mut BbiReader<impl ReadAt>, start: u32, end: u32) {
    let mut actual = vec![0.; (end - start) as usize];
    reader
        .values(&Locus::from(Region::new("c", start, end)), &mut actual)
        .unwrap();
    let expected: Vec<_> = (start..end)
        .map(|i| {
            if i < 256 && i % 4 < 3 {
                (i / 4) as f32
            } else {
                f32::NAN
            }
        })
        .collect();
    assert!(
        actual
            .iter()
            .zip(expected)
            .all(|(a, b)| a.to_bits() == b.to_bits())
    );
    let n = reader
        .query_entries(Region::new("c", start, end), true)
        .unwrap();
    let mut entries = vec![TrackEntry::default(); n];
    reader.write_entries(&mut entries).unwrap();
    let expected: Vec<_> = (0..64u32)
        .filter_map(|i| {
            let lo = (i * 4).max(start);
            let hi = (i * 4 + 3).min(end);
            (lo < hi).then_some(TrackEntry {
                start: lo,
                end: hi,
                value: i as f32,
            })
        })
        .collect();
    assert_eq!(entries, expected);

    let n = reader
        .query_entries(Region::new("c", start, end), false)
        .unwrap();
    let mut stored = vec![TrackEntry::default(); n];
    reader.write_entries(&mut stored).unwrap();
    let expected: Vec<_> = (0..64u32)
        .filter(|i| i * 4 < end && i * 4 + 3 > start)
        .map(|i| TrackEntry {
            start: i * 4,
            end: i * 4 + 3,
            value: i as f32,
        })
        .collect();
    assert_eq!(stored, expected);
}

#[test]
fn prefetched_wig_hits_do_not_read_and_zero_budget_uses_direct_queries() {
    let (source, reads, _) = source(wig());
    let mut reader = BbiReader::open(source).unwrap();
    let cache = Arc::new(BlockCache::new(1 << 20, 128));
    reader.set_cache(Some(cache.clone())).unwrap();
    reads.store(0, Ordering::Relaxed);
    reader
        .prefetch_cache(&Locus::from(Region::new("c", 0, 256)))
        .unwrap();
    assert_eq!(
        reads.load(Ordering::Relaxed),
        1,
        "adjacent misses share one read"
    );
    let initial = cache.stats();
    let mut keys = Vec::new();
    reader
        .describe_cache(&Locus::from(Region::new("c", 20, 100)), &mut keys)
        .unwrap();
    assert!(!keys.is_empty());
    reads.store(0, Ordering::Relaxed);
    for (start, end) in [(20, 100), (0, 3), (180, 256), (256, 260)] {
        assert_wig(&mut reader, start, end);
    }
    assert_eq!(reads.load(Ordering::Relaxed), 0);
    assert_eq!(cache.stats().loads, initial.loads);
    assert!(cache.evict(keys[0]));
    assert_wig(&mut reader, 20, 100);
    assert!(reads.load(Ordering::Relaxed) > 0);
    cache.resize(256, 128);
    let loads = cache.stats().loads;
    assert_wig(&mut reader, 20, 100);
    assert_eq!(
        cache.stats().loads,
        loads,
        "undersized cache bypasses preparation"
    );
    cache.resize(0, 0);
    assert_wig(&mut reader, 20, 100);
    assert_eq!(cache.stats().resident_entries, 0);
    reader.set_cache(None).unwrap();
    assert_wig(&mut reader, 20, 100);
}

fn bed(compress: bool) -> (Vec<u8>, Vec<Record>) {
    let records: Vec<_> = (0..128u32)
        .map(|i| Record {
            chrom: 0,
            start: i / 2,
            end: if i % 3 == 0 { 300 - i / 2 } else { i / 2 + 1 },
            rest: if i % 7 == 0 {
                Vec::new()
            } else {
                format!("row{i}\t+\t").into_bytes()
            },
        })
        .collect();
    let bytes = write(&Spec {
        chroms: vec![("c".into(), 512)],
        records: records.clone(),
        items_per_block: 128,
        split_chroms: true,
        chrom_fanout: 2,
        index_fanout: 2,
        auto_sql: None,
        compress,
        field_count: 6,
        defined_field_count: 3,
    });
    (bytes, records)
}

fn bed_output(
    reader: &BbiReader<impl ReadAt>,
    size: omicsweave_core::IntervalSize,
) -> Vec<(u32, u32, Vec<u8>)> {
    let mut entries = vec![IntervalEntry::default(); size.entries];
    let mut offsets = vec![0; size.entries + 1];
    let mut bytes = vec![0; size.annotation_bytes];
    reader
        .write(IntervalsOut {
            entries: &mut entries,
            annotations: StringColumn {
                offsets: &mut offsets,
                data: &mut bytes,
            },
        })
        .unwrap();
    let text = Strings {
        offsets: &offsets,
        data: &bytes,
    };
    entries
        .iter()
        .map(|r| (r.start, r.end, text.get(r.annotation as usize).to_vec()))
        .collect()
}

#[test]
fn bigbed_prepared_hits_preserve_nested_records_and_leases_survive_clear() {
    for compressed in [false, true] {
        let (bytes, records) = bed(compressed);
        let (source, reads, _) = source(bytes);
        let mut reader = BbiReader::open(source).unwrap();
        let cache = Arc::new(BlockCache::new(1 << 20, 128));
        reader.set_cache(Some(cache.clone())).unwrap();
        reader
            .prefetch_cache(&Locus::from(Region::new("c", 0, 300)))
            .unwrap();
        reads.store(0, Ordering::Relaxed);
        for (lo, hi) in [(33, 34), (20, 70), (90, 130), (0, 300), (10, 10)] {
            let size = reader.query(Region::new("c", lo, hi), false).unwrap();
            let got = bed_output(&reader, size);
            let expected: Vec<_> = records
                .iter()
                .filter(|r| lo < hi && r.end > lo && r.start < hi)
                .map(|r| (r.start, r.end, r.rest.clone()))
                .collect();
            assert_eq!(got, expected);
            let size = reader.query(Region::new("c", lo, hi), true).unwrap();
            let clipped: Vec<_> = expected
                .iter()
                .map(|(s, e, rest)| ((*s).max(lo), (*e).min(hi), rest.clone()))
                .collect();
            assert_eq!(bed_output(&reader, size), clipped);
        }
        assert_eq!(reads.load(Ordering::Relaxed), 0);
        let size = reader.query(Region::new("c", 20, 70), false).unwrap();
        let expected = bed_output(&reader, size);
        cache.clear_source(reader.cache_namespace());
        assert_eq!(cache.stats().resident_entries, 0);
        assert!(cache.stats().live_bytes > 0);
        reader.set_cache(None).unwrap();
        assert_eq!(bed_output(&reader, size), expected);
        reader.query(Region::new("c", 10, 11), false).unwrap();
        assert_eq!(cache.stats().live_bytes, 0);
    }
}

#[test]
fn failed_loads_publish_nothing_and_retry_with_the_same_source_session() {
    let (source, reads, fail) = source(wig());
    let mut reader = BbiReader::open(source).unwrap();
    let cache = Arc::new(BlockCache::new(1 << 20, 128));
    reader.set_cache(Some(cache.clone())).unwrap();
    fail.store(true, Ordering::Relaxed);
    assert!(
        reader
            .prefetch_cache(&Locus::from(Region::new("c", 0, 10)))
            .is_err()
    );
    assert!(cache.keys().is_empty());
    assert_wig(&mut reader, 0, 10);
    reads.store(0, Ordering::Relaxed);
    assert_wig(&mut reader, 1, 9);
    assert_eq!(reads.load(Ordering::Relaxed), 0);
    let mut keys = Vec::new();
    reader
        .describe_cache(&Locus::from(Region::new("c", 0, 10)), &mut keys)
        .unwrap();
    cache.clear();
    cache
        .get_or_load(keys[0], || Ok((vec![0u8; 32], 32)))
        .unwrap();
    assert!(
        reader
            .values(&Locus::from(Region::new("c", 0, 10)), &mut [0.; 10])
            .is_err()
    );
    cache.clear();
    assert_wig(&mut reader, 0, 10);
}

#[test]
fn cloned_file_contexts_share_cache_but_independent_opens_are_isolated() {
    let dir = tempfile::tempdir().unwrap();
    let path = dir.path().join("data.bw");
    std::fs::write(&path, wig()).unwrap();
    let mut reader = BbiReader::open(FileSource::open(&path).unwrap()).unwrap();
    let cache = Arc::new(BlockCache::new(1 << 20, 128));
    reader.set_cache(Some(cache.clone())).unwrap();
    reader
        .prefetch_cache(&Locus::from(Region::new("c", 0, 100)))
        .unwrap();
    let before = cache.stats().loads;
    let mut cloned = reader.try_clone().unwrap();
    assert_eq!(cloned.cache_namespace(), reader.cache_namespace());
    drop(reader);
    assert_wig(&mut cloned, 20, 90);
    assert_eq!(cache.stats().loads, before);
    let independent = BbiReader::open(FileSource::open(&path).unwrap()).unwrap();
    assert_ne!(independent.cache_namespace(), cloned.cache_namespace());
}
