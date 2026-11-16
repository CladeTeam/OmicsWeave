//! Per-operation block reuse preserves direct reads and external cache budgets.

use std::cell::Cell;
use std::io::Cursor;
use std::sync::Arc;

use omicsweave_core::cache::{BlockCache, BlockKey, SourceToken};
use omicsweave_core::{Locus, ReadAt, Result, SequenceRead, Span, Strand};
use omicsweave_fasta::{
    FastaReader, TwoBitPiece, TwoBitReader, TwoBitRecord, TwoBitWriter, compress_fasta,
};

struct Source {
    data: Vec<u8>,
    reads: Cell<usize>,
    bytes: Cell<usize>,
    fail: Cell<Option<usize>>,
}

impl Source {
    fn new(data: Vec<u8>) -> Self {
        Self {
            data,
            reads: Cell::new(0),
            bytes: Cell::new(0),
            fail: Cell::new(None),
        }
    }

    fn reset(&self) {
        self.reads.set(0);
        self.bytes.set(0);
    }
}

impl ReadAt for Source {
    fn len(&self) -> u64 {
        self.data.len() as u64
    }

    fn read_exact_at(&self, out: &mut [u8], at: u64) -> Result<()> {
        let call = self.reads.get() + 1;
        self.reads.set(call);
        self.bytes.set(self.bytes.get() + out.len());
        if self.fail.get() == Some(call) {
            out.fill(0);
            return Err(omicsweave_core::Error::invalid(
                "source",
                "injected I/O failure",
            ));
        }
        self.data.as_slice().read_exact_at(out, at)
    }
}

fn genome(n: usize) -> Vec<u8> {
    (0..n).map(|i| b"ACgtN"[(i * 7 + i / 29) % 5]).collect()
}

fn expected(bases: &[u8], spans: &[Span], strand: Strand) -> Vec<u8> {
    let mut out = Vec::new();
    for s in spans {
        out.extend_from_slice(&bases[s.start as usize..s.end as usize]);
    }
    if strand == Strand::Reverse {
        out.reverse();
        for b in &mut out {
            *b = match *b {
                b'A' => b'T',
                b'C' => b'G',
                b'G' => b'C',
                b'T' => b'A',
                b'a' => b't',
                b'c' => b'g',
                b'g' => b'c',
                b't' => b'a',
                x => x,
            };
        }
    }
    out
}

fn fasta(bases: &[u8]) -> Vec<u8> {
    let mut text = b">s\n".to_vec();
    text.extend_from_slice(bases);
    text.push(b'\n');
    text
}

fn twobit(bases: &[u8]) -> Vec<u8> {
    let mut writer = TwoBitWriter::create(Cursor::new(Vec::new()), vec!["s".into()], 0).unwrap();
    let piece = TwoBitPiece::encode(0, bases).unwrap();
    let mut record = TwoBitRecord::new();
    record.push(&piece).unwrap();
    writer.add("s", &record, &mut &piece.packed[..]).unwrap();
    writer.finish().unwrap().into_inner()
}

#[test]
fn bgzf_spans_share_a_block_even_when_admission_fails() {
    let bases = genome(180_000);
    let (bytes, fai, gzi) = compress_fasta(fasta(&bases).as_slice(), Vec::new(), 1).unwrap();
    let spans: Vec<_> = (0..128).map(|i| Span::new(i * 33, i * 33 + 9)).collect();
    for strand in [Strand::Forward, Strand::Reverse] {
        let locus = Locus::new("s", &spans, strand).unwrap();
        for budget in [0, 1, 1 << 20] {
            let cache = Arc::new(BlockCache::new(budget, 32));
            let source = Source::new(bytes.clone());
            let mut reader = FastaReader::new(&source, fai.clone(), Some(gzi.clone()));
            reader.set_cache(Some(cache.clone())).unwrap();
            let mut out = vec![0xa5; locus.len() as usize + 2];
            reader
                .sequence(&locus, &mut out[1..locus.len() as usize + 1])
                .unwrap();
            assert_eq!(&out[1..out.len() - 1], expected(&bases, &spans, strand));
            assert_eq!((out[0], *out.last().unwrap()), (0xa5, 0xa5));
            assert_eq!(source.reads.get(), 1);
            assert_eq!(cache.stats().loads, 1);
            let before = source.reads.get();
            reader
                .sequence(&locus, &mut out[1..locus.len() as usize + 1])
                .unwrap();
            assert_eq!(source.reads.get() - before, usize::from(budget <= 1));
            cache.clear();
            let before = source.reads.get();
            let loads = cache.stats().loads;
            reader.prefetch_cache(&locus).unwrap();
            assert_eq!(source.reads.get() - before, 1);
            assert_eq!(cache.stats().loads - loads, 1);
            if budget <= 1 {
                assert_eq!(cache.stats().live_bytes, 0);
            }
            let mut keys = Vec::new();
            reader.describe_cache(&locus, &mut keys).unwrap();
            assert_eq!(keys.len(), 1);
        }
    }
}

#[test]
fn failed_spliced_reads_do_not_reuse_partial_input_on_retry() {
    let bases = genome(180_000);
    let (bytes, fai, gzi) = compress_fasta(fasta(&bases).as_slice(), Vec::new(), 1).unwrap();
    let spans = [Span::new(10, 20), Span::new(170_000, 170_010)];
    let locus = Locus::new("s", &spans, Strand::Reverse).unwrap();
    for budget in [0, 1, 1 << 20] {
        let source = Source::new(bytes.clone());
        source.fail.set(Some(2));
        let mut reader = FastaReader::new(&source, fai.clone(), Some(gzi.clone()));
        reader
            .set_cache(Some(Arc::new(BlockCache::new(budget, 32))))
            .unwrap();
        let mut out = [0u8; 20];
        assert!(reader.sequence(&locus, &mut out).is_err());
        source.fail.set(None);
        reader.sequence(&locus, &mut out).unwrap();
        assert_eq!(out.as_slice(), expected(&bases, &spans, Strand::Reverse));
    }
}

#[test]
fn twobit_direct_and_tiny_cache_keep_exact_packed_reads() {
    let bases = genome(180_000);
    let bytes = twobit(&bases);
    for budget in [None, Some(0), Some(1)] {
        let source = Source::new(bytes.clone());
        let mut reader = TwoBitReader::open(&source, true).unwrap();
        if let Some(budget) = budget {
            reader
                .set_cache(Some(Arc::new(BlockCache::new(budget, 32))))
                .unwrap();
        }
        for spans in [
            vec![Span::new(3, 11)],
            vec![Span::new(3, 131)],
            vec![Span::new(3, 11), Span::new(19, 27), Span::new(51, 59)],
            vec![Span::new(3, 11), Span::new(65_000, 65_008)],
        ] {
            let locus = Locus::new("s", &spans, Strand::Forward).unwrap();
            let mut out = vec![0; locus.len() as usize];
            source.reset();
            reader.sequence(&locus, &mut out).unwrap();
            assert_eq!(out, expected(&bases, &spans, Strand::Forward));
            let packed_bytes: usize = spans
                .iter()
                .map(|s| s.end.div_ceil(4) as usize - s.start as usize / 4)
                .sum();
            assert_eq!(source.reads.get(), spans.len());
            assert_eq!(source.bytes.get(), packed_bytes);
        }
    }
}

#[test]
fn twobit_operation_lease_reuses_unadmitted_tile_without_pinning_another() {
    let bases = genome(180_000);
    let source = Source::new(twobit(&bases));
    let mut reader = TwoBitReader::open(&source, true).unwrap();
    let cache = Arc::new(BlockCache::new(128 << 10, 32));
    let key = BlockKey::new(SourceToken::new(), 100, [0, 0], 0, 0);
    let pinned = cache
        .get_or_load(key, || Ok((vec![0u8; 100_000], 100_000)))
        .unwrap();
    cache.evict(key);
    let pinned_bytes = cache.stats().live_bytes;
    reader.set_cache(Some(cache.clone())).unwrap();
    let spans = [
        Span::new(3, 11),
        Span::new(65_000, 65_008),
        Span::new(66_000, 66_008),
    ];
    let locus = Locus::new("s", &spans, Strand::Reverse).unwrap();
    source.reset();
    let mut out = vec![0; locus.len() as usize];
    reader.sequence(&locus, &mut out).unwrap();
    assert_eq!(out, expected(&bases, &spans, Strand::Reverse));
    assert_eq!(source.reads.get(), 2);
    assert_eq!(source.bytes.get(), 2 * (65_536 / 4));
    assert_eq!(cache.stats().live_bytes, pinned_bytes);
    drop(pinned);
    assert_eq!(cache.stats().live_bytes, 0);
}

#[test]
fn sparse_uncached_bgzf_locus_keeps_output_order_across_bounded_runs() {
    let bases = genome(4 << 20);
    let (bytes, fai, gzi) = compress_fasta(fasta(&bases).as_slice(), Vec::new(), 1).unwrap();
    let offset = fai.get("s").unwrap().offset;
    let mut spans = Vec::new();
    for &(_, text) in gzi.starts().iter().skip(1) {
        if text > offset && text + 1 < offset + bases.len() as u64 {
            let p = (text - offset) as u32;
            spans.push(Span::new(p - 1, p + 1));
        }
    }
    spans.push(Span::new(bases.len() as u32, bases.len() as u32));
    let source = Source::new(bytes);
    let mut reader = FastaReader::new(&source, fai, Some(gzi));
    for strand in [Strand::Forward, Strand::Reverse] {
        let locus = Locus::new("s", &spans, strand).unwrap();
        let mut out = vec![0; locus.len() as usize];
        reader.sequence(&locus, &mut out).unwrap();
        assert_eq!(out, expected(&bases, &spans, strand));
    }
}
