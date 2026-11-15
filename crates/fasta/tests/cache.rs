//! Cache lifetime, actual positional reads and independent sequence oracles.
mod common;

use std::cell::{Cell, RefCell};
use std::sync::Arc;

use omicsweave_core::cache::BlockCache;
use omicsweave_core::{FileSource, Locus, ReadAt, Region, Result, SequenceRead};
use omicsweave_fasta::{EOF_BLOCK, Fai, FaiEntry, FastaReader, Gzi, TwoBitReader, compress_fasta};

struct Source {
    bytes: RefCell<Vec<u8>>,
    reads: Cell<usize>,
    read_bytes: Cell<usize>,
}

impl Source {
    fn new(bytes: Vec<u8>) -> Self {
        Self {
            bytes: RefCell::new(bytes),
            reads: Cell::new(0),
            read_bytes: Cell::new(0),
        }
    }
}

impl ReadAt for Source {
    fn len(&self) -> u64 {
        self.bytes.borrow().len() as u64
    }

    fn read_exact_at(&self, out: &mut [u8], at: u64) -> Result<()> {
        self.reads.set(self.reads.get() + 1);
        self.read_bytes.set(self.read_bytes.get() + out.len());
        self.bytes.borrow().as_slice().read_exact_at(out, at)
    }
}

fn fixture(n: usize) -> (Vec<u8>, Vec<u8>, Fai, Gzi) {
    let bases: Vec<_> = (0..n).map(|i| b"ACgtNNTA"[(i * 7 + i / 29) % 8]).collect();
    let mut text = b">s\n".to_vec();
    for line in bases.chunks(60) {
        text.extend_from_slice(line);
        text.push(b'\n');
    }
    let (bytes, fai, gzi) = compress_fasta(text.as_slice(), Vec::new(), 1).unwrap();
    (bases, bytes, fai, gzi)
}

#[test]
fn bgzf_adjacent_queries_hit_without_source_reads_and_external_clear_reloads() {
    let (bases, bytes, fai, gzi) = fixture(180_000);
    let source = Source::new(bytes);
    let cache = Arc::new(BlockCache::new(1 << 20, 32));
    let mut reader = FastaReader::new(&source, fai, Some(gzi));
    reader.set_cache(Some(cache.clone())).unwrap();
    let mut out = vec![0; 130_000];
    reader
        .sequence(&Locus::from(Region::new("s", 17, 130_017)), &mut out)
        .unwrap();
    assert_eq!(out, bases[17..130_017]);
    let reads = source.reads.get();
    let loads = cache.stats().loads;
    assert!(loads >= 2);
    // Compressed bytes of adjacent misses are fetched as a run.
    assert_eq!(reads, 1);
    reader
        .sequence(&Locus::from(Region::new("s", 19, 130_019)), &mut out)
        .unwrap();
    assert_eq!(out, bases[19..130_019]);
    assert_eq!(source.reads.get(), reads);
    assert_eq!(cache.stats().loads, loads);
    cache.clear_source(reader.cache_namespace());
    reader
        .sequence(&Locus::from(Region::new("s", 19, 130_019)), &mut out)
        .unwrap();
    assert_eq!(out, bases[19..130_019]);
    assert!(source.reads.get() > reads);
    reader.set_cache(None).unwrap();
    let reads = source.reads.get();
    reader
        .sequence(&Locus::from(Region::new("s", 19, 130_019)), &mut out)
        .unwrap();
    assert!(source.reads.get() > reads);
}

#[test]
fn plain_fasta_rejects_data_cache_and_empty_bgzf_query_needs_no_blocks() {
    let (bases, bytes, fai, gzi) = fixture(100);
    let cache = Arc::new(BlockCache::new(4096, 8));
    let mut plain = FastaReader::new(bases, fai.clone(), None);
    assert!(plain.set_cache(Some(cache.clone())).is_err());
    let mut keys = Vec::new();
    assert!(
        plain
            .describe_cache(&Locus::from(Region::new("s", 0, 1)), &mut keys)
            .is_err()
    );
    let mut reader = FastaReader::new(bytes, fai, Some(gzi));
    reader.set_cache(Some(cache.clone())).unwrap();
    reader
        .describe_cache(&Locus::from(Region::new("s", 17, 17)), &mut keys)
        .unwrap();
    assert!(keys.is_empty());
    reader
        .sequence(&Locus::from(Region::new("s", 17, 17)), &mut [])
        .unwrap();
    assert_eq!(cache.stats().loads, 0);
}

#[test]
fn admission_rejection_and_wrong_payload_do_not_bypass_checks() {
    let (bases, bytes, fai, gzi) = fixture(1000);
    let source = Source::new(bytes);
    let cache = Arc::new(BlockCache::new(0, 0));
    let mut reader = FastaReader::new(&source, fai, Some(gzi));
    reader.set_cache(Some(cache.clone())).unwrap();
    let mut out = [0; 128];
    for _ in 0..2 {
        reader
            .sequence(&Locus::from(Region::new("s", 17, 145)), &mut out)
            .unwrap();
        assert_eq!(out.as_slice(), &bases[17..145]);
    }
    assert_eq!(source.reads.get(), 2);
    assert_eq!(cache.stats().resident_entries, 0);
    assert!(cache.stats().admission_rejections >= 2);
    cache.resize(1 << 20, 10);
    let mut keys = Vec::new();
    reader
        .describe_cache(&Locus::from(Region::new("s", 17, 145)), &mut keys)
        .unwrap();
    cache
        .get_or_load(keys[0], || Ok((vec![b'X'; 65536], 65536)))
        .unwrap();
    assert!(
        reader
            .sequence(&Locus::from(Region::new("s", 17, 145)), &mut out)
            .is_err()
    );
    cache.clear();
    reader
        .sequence(&Locus::from(Region::new("s", 17, 145)), &mut out)
        .unwrap();
    assert_eq!(out.as_slice(), &bases[17..145]);
}

#[test]
fn leading_empty_members_are_validated_and_reused_as_the_indexed_data_block() {
    let bases = b"ACGTacgtNNNN";
    let mut text = b">s\n".to_vec();
    text.extend_from_slice(bases);
    text.push(b'\n');
    let (data, _, _) = compress_fasta(text.as_slice(), Vec::new(), 1).unwrap();
    let mut bytes = EOF_BLOCK.repeat(3000);
    bytes.extend_from_slice(&data);
    let source = Source::new(bytes);
    let fai = Fai::from_entries([FaiEntry {
        name: "s".into(),
        len: 12,
        offset: 3,
        line_bases: 12,
        line_width: 13,
    }]);
    let gzi = Gzi::parse(&0u64.to_le_bytes()).unwrap();
    let cache = Arc::new(BlockCache::new(4096, 8));
    let mut reader = FastaReader::new(&source, fai, Some(gzi));
    reader.set_cache(Some(cache.clone())).unwrap();
    // A corrupt empty member must not publish any usable data block.
    source.bytes.borrow_mut()[2339 * EOF_BLOCK.len() + 20] ^= 1;
    assert!(
        reader
            .sequence(&Locus::from(Region::new("s", 0, 12)), &mut [0; 12])
            .is_err()
    );
    assert_eq!(cache.stats().resident_entries, 0);
    source.bytes.borrow_mut()[2339 * EOF_BLOCK.len() + 20] ^= 1;
    let mut out = [0; 12];
    reader
        .sequence(&Locus::from(Region::new("s", 0, 12)), &mut out)
        .unwrap();
    assert_eq!(&out, bases);
    let reads = source.reads.get();
    reader
        .sequence(&Locus::from(Region::new("s", 0, 12)), &mut out)
        .unwrap();
    assert_eq!(&out, bases);
    assert_eq!(source.reads.get(), reads);
}

#[test]
fn twobit_tiles_match_independent_mask_and_n_oracle_across_boundaries() {
    let seq = common::TwoBitSeq {
        name: "s".into(),
        bases: (0..101).map(|i| b"ACTG"[i % 4]).collect(),
        n_blocks: vec![(5, 3), (30, 5), (95, 6)],
        mask_blocks: vec![(2, 10), (28, 15), (60, 39)],
    };
    let bytes = common::write_twobit(std::slice::from_ref(&seq), 1);
    let source = Source::new(bytes);
    let cache = Arc::new(BlockCache::new(4096, 16));
    let mut reader = TwoBitReader::open(&source, true).unwrap();
    reader.set_cache_with_tile(Some(cache.clone()), 32).unwrap();
    let before = source.reads.get();
    let mut out = [0; 101];
    reader
        .sequence(&Locus::from(Region::new("s", 0, 101)), &mut out)
        .unwrap();
    for (p, &value) in out.iter().enumerate() {
        let covers = |&(s, n): &(u32, u32)| s as usize <= p && p < (s + n) as usize;
        let expected = if seq.n_blocks.iter().any(covers) {
            b'N'
        } else if seq.mask_blocks.iter().any(covers) {
            seq.bases[p].to_ascii_lowercase()
        } else {
            seq.bases[p]
        };
        assert_eq!(value, expected, "base {p}");
    }
    assert!(source.reads.get() > before);
    let reads = source.reads.get();
    for start in 0..=101 {
        for end in start..=101 {
            let mut actual = vec![0; (end - start) as usize];
            reader
                .sequence(&Locus::from(Region::new("s", start, end)), &mut actual)
                .unwrap();
            assert_eq!(actual, out[start as usize..end as usize]);
        }
    }
    assert_eq!(source.reads.get(), reads);
    assert_eq!(cache.stats().loads, 4);
    assert!(reader.set_cache_with_tile(Some(cache), 3).is_err());
}

#[test]
fn file_context_clone_shares_namespace_and_blocks_but_not_lifetime() {
    let (bases, bytes, fai, gzi) = fixture(1000);
    let dir = tempfile::tempdir().unwrap();
    let path = dir.path().join("ref.fa.gz");
    std::fs::write(&path, bytes).unwrap();
    let cache = Arc::new(BlockCache::new(1 << 20, 10));
    let mut reader = FastaReader::new(FileSource::open(&path).unwrap(), fai, Some(gzi));
    reader.set_cache(Some(cache.clone())).unwrap();
    let mut out = [0; 128];
    reader
        .sequence(&Locus::from(Region::new("s", 17, 145)), &mut out)
        .unwrap();
    let namespace = reader.cache_namespace();
    let mut clone = reader.try_clone().unwrap();
    assert_eq!(namespace, clone.cache_namespace());
    drop(reader);
    clone
        .sequence(&Locus::from(Region::new("s", 19, 147)), &mut out)
        .unwrap();
    assert_eq!(out.as_slice(), &bases[19..147]);
    assert_eq!(cache.stats().loads, 1);
}

#[test]
fn prefetch_populates_both_formats_without_any_query_output() {
    let (bases, bytes, fai, gzi) = fixture(1000);
    let source = Source::new(bytes);
    let cache = Arc::new(BlockCache::new(1 << 20, 100));
    let mut reader = FastaReader::new(&source, fai, Some(gzi));
    assert!(
        reader
            .prefetch_cache(&Locus::from(Region::new("s", 0, 1000)))
            .is_err()
    );
    reader.set_cache(Some(cache.clone())).unwrap();
    reader
        .prefetch_cache(&Locus::from(Region::new("s", 0, 1000)))
        .unwrap();
    let reads = source.reads.get();
    let mut out = [0; 128];
    reader
        .sequence(&Locus::from(Region::new("s", 23, 151)), &mut out)
        .unwrap();
    assert_eq!(out.as_slice(), &bases[23..151]);
    assert_eq!(source.reads.get(), reads);

    let seq = common::TwoBitSeq {
        name: "s".into(),
        bases: b"ACGT".repeat(1000),
        n_blocks: vec![(10, 3)],
        mask_blocks: vec![(8, 10)],
    };
    let source = Source::new(common::write_twobit(std::slice::from_ref(&seq), 1));
    let mut reader = TwoBitReader::open(&source, true).unwrap();
    reader.set_cache(Some(cache)).unwrap();
    reader
        .prefetch_cache(&Locus::from(Region::new("s", 0, 100)))
        .unwrap();
    let reads = source.reads.get();
    let mut out = [0; 24];
    reader
        .sequence(&Locus::from(Region::new("s", 0, 24)), &mut out)
        .unwrap();
    assert_eq!(&out, b"ACGTACGTacNNNcgtacGTACGT");
    assert_eq!(source.reads.get(), reads);
}

#[test]
fn undersized_twobit_cache_reads_only_requested_packed_range() {
    let seq = common::TwoBitSeq {
        name: "s".into(),
        bases: b"ACGT".repeat(25_000),
        n_blocks: Vec::new(),
        mask_blocks: Vec::new(),
    };
    let source = Source::new(common::write_twobit(std::slice::from_ref(&seq), 1));
    let cache = Arc::new(BlockCache::new(64, 10));
    let mut reader = TwoBitReader::open(&source, false).unwrap();
    reader.set_cache(Some(cache.clone())).unwrap();
    let before = source.read_bytes.get();
    reader
        .prefetch_cache(&Locus::from(Region::new("s", 3, 131)))
        .unwrap();
    assert_eq!(source.read_bytes.get(), before);
    let mut out = [0; 128];
    reader
        .sequence(&Locus::from(Region::new("s", 3, 131)), &mut out)
        .unwrap();
    assert_eq!(out.as_slice(), &seq.bases[3..131]);
    assert_eq!(source.read_bytes.get() - before, 33);
    assert_eq!(cache.stats().loads, 0);
}
