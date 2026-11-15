use std::sync::Arc;
use std::sync::atomic::{AtomicBool, AtomicUsize, Ordering};

use omicsweave_bgzf::{
    BLOCK_TEXT_LEN, BinIndex, Chunk, EOF_BLOCK, SeekReader, check_level, deflate_block,
};
use omicsweave_core::cache::{BlockCache, SourceToken};
use omicsweave_core::{Error, ReadAt, Result};

struct Source {
    bytes: Vec<u8>,
    reads: AtomicUsize,
    fail: AtomicBool,
}

impl ReadAt for Source {
    fn len(&self) -> u64 {
        self.bytes.len() as u64
    }
    fn read_exact_at(&self, out: &mut [u8], at: u64) -> Result<()> {
        self.reads.fetch_add(1, Ordering::Relaxed);
        if self.fail.swap(false, Ordering::Relaxed) {
            return Err(Error::invalid("fixture", "transient read failure"));
        }
        self.bytes.read_exact_at(out, at)
    }
}

fn fixture() -> (Source, Vec<u8>, Chunk, Vec<u64>) {
    let mut compressor = libdeflater::Compressor::new(check_level(6).unwrap());
    let mut bytes = Vec::new();
    let mut text = Vec::new();
    let mut addresses = Vec::new();
    for (i, len) in [3000, 5000, 7000].into_iter().enumerate() {
        addresses.push(bytes.len() as u64);
        let payload: Vec<_> = (0..len).map(|j| b"ACGTN"[(j + i) % 5]).collect();
        deflate_block(&mut compressor, &payload, &mut bytes).unwrap();
        text.extend_from_slice(&payload);
        if i == 0 {
            addresses.push(bytes.len() as u64);
            bytes.extend_from_slice(&EOF_BLOCK);
        }
    }
    let end = bytes.len() as u64;
    bytes.extend_from_slice(&EOF_BLOCK);
    (
        Source {
            bytes,
            reads: AtomicUsize::new(0),
            fail: AtomicBool::new(false),
        },
        text,
        Chunk {
            start: 7,
            end: end << 16,
        },
        addresses,
    )
}

fn query(reader: &mut SeekReader, source: &Source, chunk: Chunk, expected: &[u8]) {
    reader.reset_query();
    reader.seek(source, chunk.start, chunk.end).unwrap();
    let range = reader.take(source, expected.len(), "fixture text").unwrap();
    assert_eq!(reader.slice(range), expected);
    assert_eq!(reader.tell(), chunk.end);
}

#[test]
fn every_uncached_query_reloads_but_one_query_batches_members() {
    let (source, text, chunk, _) = fixture();
    let mut reader = SeekReader::new();
    for i in 1..=3 {
        query(&mut reader, &source, chunk, &text[7..]);
        assert_eq!(reader.stats().read_calls, i);
        assert_eq!(reader.stats().inflate_calls, i * 4);
    }
    assert_eq!(source.reads.load(Ordering::Relaxed), 3);
}

#[test]
fn header_description_prefetch_clone_and_clear_share_verified_members() {
    let (source, text, chunk, addresses) = fixture();
    let cache = Arc::new(BlockCache::new(1 << 20, 64));
    let mut reader = SeekReader::new();
    let token = SourceToken::new();
    reader.set_cache(Some(cache.clone()), token);
    let mut keys = Vec::new();
    reader
        .describe_cache(&source, &[chunk, chunk], &mut keys)
        .unwrap();
    assert_eq!(keys.iter().map(|k| k.block).collect::<Vec<_>>(), addresses);
    assert!(keys.iter().all(|k| k.kind == 5 && k.source == token));
    assert_eq!(reader.stats().inflate_calls, 0);
    assert_eq!(reader.stats().read_calls, 1);
    reader.prefetch_cache(&source, &[chunk]).unwrap();
    assert_eq!(reader.stats().inflate_calls, 4);
    assert_eq!(reader.stats().read_calls, 1);
    let mut clone = reader.clone_context();
    assert_eq!(clone.cache_namespace(), token);
    query(&mut clone, &source, chunk, &text[7..]);
    assert_eq!(clone.stats().read_calls, 0);
    assert_eq!(clone.stats().inflate_calls, 0);
    cache.clear_source(token);
    assert_eq!(clone.peek(), b"");
    query(&mut clone, &source, chunk, &text[7..]);
    assert_eq!(clone.stats().read_calls, 1);
    assert_eq!(clone.stats().inflate_calls, 4);
    let mut independent = SeekReader::new();
    independent.set_cache(Some(cache), SourceToken::new());
    query(&mut independent, &source, chunk, &text[7..]);
    assert_eq!(independent.stats().inflate_calls, 4);
}

#[test]
fn zero_budget_uses_direct_decode_and_wrong_payload_is_a_safe_error() {
    let (source, text, chunk, _) = fixture();
    let cache = Arc::new(BlockCache::new(0, 0));
    let mut reader = SeekReader::new();
    reader.set_cache(Some(cache.clone()), SourceToken::new());
    reader.prefetch_cache(&source, &[chunk]).unwrap();
    assert_eq!(reader.stats().inflate_calls, 0);
    assert!(
        reader
            .prefetch_cache(&source, &[Chunk { start: 8, end: 7 }])
            .is_err()
    );
    for _ in 0..2 {
        query(&mut reader, &source, chunk, &text[7..]);
    }
    assert_eq!(reader.stats().inflate_calls, 8);
    assert_eq!(cache.stats().loads, 0);
    cache.resize(1 << 20, 64);
    let mut keys = Vec::new();
    reader.describe_cache(&source, &[chunk], &mut keys).unwrap();
    cache.get_or_load(keys[0], || Ok(("wrong", 0))).unwrap();
    reader.reset_query();
    assert!(reader.seek(&source, chunk.start, chunk.end).is_err());
    cache.clear();
    query(&mut reader, &source, chunk, &text[7..]);
}

#[test]
fn sorted_overlapping_prefetch_visits_each_physical_member_once_without_admission() {
    let (source, _, chunk, addresses) = fixture();
    let cache = Arc::new(BlockCache::new(1, 64));
    let mut reader = SeekReader::new();
    reader.set_cache(Some(cache.clone()), SourceToken::new());
    let mut chunks = vec![chunk; 128];
    chunks.push(Chunk {
        start: addresses[2] << 16 | 3,
        end: chunk.end,
    });
    let mut keys = Vec::new();
    reader.describe_cache(&source, &chunks, &mut keys).unwrap();
    assert_eq!(keys.iter().map(|k| k.block).collect::<Vec<_>>(), addresses);
    assert_eq!(reader.stats().inflate_calls, 0);
    reader.prefetch_cache(&source, &chunks).unwrap();
    assert_eq!(reader.stats().inflate_calls, addresses.len() as u64);
    assert_eq!(cache.stats().resident_entries, 0);

    // A later request in the middle of a compressed member must not be
    // hidden by an earlier request that covers that member.
    let invalid = Chunk {
        start: 1 << 16,
        end: chunk.end,
    };
    assert!(
        reader
            .describe_cache(&source, &[chunk, invalid], &mut keys)
            .is_err()
    );
    assert!(reader.prefetch_cache(&source, &[chunk, invalid]).is_err());
}

#[test]
fn unsorted_cache_requests_preserve_first_seen_order() {
    let (source, _, chunk, addresses) = fixture();
    let tail = Chunk {
        start: addresses[3] << 16,
        end: chunk.end,
    };
    let head = Chunk {
        start: chunk.start,
        end: addresses[1] << 16,
    };
    let cache = Arc::new(BlockCache::new(1, 64));
    let mut reader = SeekReader::new();
    reader.set_cache(Some(cache), SourceToken::new());
    let mut keys = Vec::new();
    reader
        .describe_cache(&source, &[tail, head, tail], &mut keys)
        .unwrap();
    assert_eq!(
        keys.iter().map(|k| k.block).collect::<Vec<_>>(),
        [addresses[3], addresses[0]]
    );
    reader.prefetch_cache(&source, &[tail, head, tail]).unwrap();
    assert_eq!(reader.stats().inflate_calls, 3);
}

#[test]
fn crc_errors_never_publish_and_read_failures_can_retry() {
    let (mut source, text, chunk, _) = fixture();
    source.fail.store(true, Ordering::Relaxed);
    let cache = Arc::new(BlockCache::new(1 << 20, 64));
    let mut reader = SeekReader::new();
    reader.set_cache(Some(cache.clone()), SourceToken::new());
    assert!(reader.seek(&source, chunk.start, chunk.end).is_err());
    assert_eq!(cache.stats().resident_entries, 0);
    reader.seek(&source, chunk.start, chunk.end).unwrap();
    let range = reader.take(&source, text.len() - 7, "text").unwrap();
    assert_eq!(reader.slice(range), &text[7..]);
    query(&mut reader, &source, chunk, &text[7..]);
    let first_len = u16::from_le_bytes([source.bytes[16], source.bytes[17]]) as usize + 1;
    source.bytes[first_len - 8] ^= 1;
    reader.set_cache(Some(cache.clone()), SourceToken::new());
    for _ in 0..2 {
        reader.reset_query();
        assert!(reader.seek(&source, chunk.start, chunk.end).is_err());
    }
    assert_eq!(cache.keys().len(), 4);
}

#[test]
fn empty_members_do_not_change_within_offsets_and_chunks_may_be_empty() {
    let (source, text, chunk, addresses) = fixture();
    let mut reader = SeekReader::new();
    reader
        .seek(&source, addresses[1] << 16 | 3, chunk.end)
        .unwrap();
    let range = reader.take(&source, 12, "text").unwrap();
    assert_eq!(reader.slice(range), &text[3000 + 3..3000 + 15]);
    assert!(reader.seek(&source, 3001, 0).is_err());
    let at = source.len() << 16;
    reader.seek(&source, at, at).unwrap();
    assert_eq!(reader.fill(&source, 1).unwrap(), 0);
    assert!(reader.seek(&source, at | 1, 0).is_err());
    let mut keys = Vec::new();
    reader
        .describe_cache(&source, &[Chunk { start: 7, end: 7 }], &mut keys)
        .unwrap();
    assert!(keys.is_empty());
}

#[test]
fn speculative_fill_compaction_moves_no_more_text_than_it_discards() {
    let mut compressor = libdeflater::Compressor::new(check_level(1).unwrap());
    let payload = vec![b'x'; BLOCK_TEXT_LEN];
    let mut member = Vec::new();
    deflate_block(&mut compressor, &payload, &mut member).unwrap();
    let mut source = member.repeat(512);
    source.extend_from_slice(&EOF_BLOCK);
    let mut reader = SeekReader::new();
    reader.seek(&source, 0, 0).unwrap();
    for _ in 0..256 {
        reader.fill(&source, 12 << 20).unwrap();
        let range = reader.take(&source, 65536, "text").unwrap();
        assert!(reader.slice(range).iter().all(|&b| b == b'x'));
    }
    let stats = reader.stats();
    assert!(stats.discarded_bytes > 0);
    assert!(stats.compacted_bytes <= stats.discarded_bytes);
}

#[test]
fn full_signed_coordinate_csi_geometry_does_not_turn_negative() {
    for (shift, levels) in [(63i32, 0i32), (33, 10)] {
        let mut bytes = b"CSI\x01".to_vec();
        for n in [shift, levels, 0, 1, 1, 0] {
            bytes.extend_from_slice(&n.to_le_bytes());
        }
        bytes.extend_from_slice(&0u64.to_le_bytes());
        bytes.extend_from_slice(&1i32.to_le_bytes());
        bytes.extend_from_slice(&65536u64.to_le_bytes());
        bytes.extend_from_slice(&131072u64.to_le_bytes());
        let index = BinIndex::parse(&bytes).unwrap();
        let expected = vec![Chunk {
            start: 65536,
            end: 131072,
        }];
        assert_eq!(index.chunks(0, 0, 1), expected);
        assert_eq!(index.chunks(0, i64::MAX - 2, i64::MAX), expected);
        assert!(index.chunks(0, 0, 0).is_empty());
    }
}

#[test]
fn an_indexed_tiny_chunk_reads_one_member_window_not_an_unhinted_window() {
    let mut compressor = libdeflater::Compressor::new(check_level(6).unwrap());
    let mut source = Vec::new();
    deflate_block(&mut compressor, b"small", &mut source).unwrap();
    let end = source.len() as u64;
    source.extend_from_slice(&EOF_BLOCK);
    source.resize(1 << 20, 0);
    let mut reader = SeekReader::new();
    reader.seek(&source, 0, end << 16).unwrap();
    let range = reader.take(&source, 5, "text").unwrap();
    assert_eq!(reader.slice(range), b"small");
    assert_eq!(reader.stats().read_calls, 1);
    assert_eq!(reader.stats().read_bytes, 65536);
    let mut unhinted = SeekReader::new();
    unhinted.seek(&source, 0, 0).unwrap();
    assert_eq!(unhinted.stats().read_bytes, 256 << 10);
}

#[test]
fn physical_member_end_must_fit_a_virtual_offset() {
    struct Tail(u64);
    impl ReadAt for Tail {
        fn len(&self) -> u64 {
            self.0 + EOF_BLOCK.len() as u64
        }
        fn read_exact_at(&self, out: &mut [u8], at: u64) -> Result<()> {
            EOF_BLOCK.as_slice().read_exact_at(out, at - self.0)
        }
    }
    let source = Tail((u64::MAX >> 16) - 1);
    let mut reader = SeekReader::new();
    assert!(reader.seek(&source, source.0 << 16, 0).is_err());
    assert_eq!(reader.stats().inflate_calls, 0);
}
