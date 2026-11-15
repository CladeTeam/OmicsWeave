//! Multi-span and reverse-strand loci read the concatenated spans, reverse
//! complemented on `-`, from every sequence format and read path.

use std::io::Cursor;
use std::sync::Arc;

use omicsweave_core::cache::BlockCache;
use omicsweave_core::orient::{complement, reverse_complement};
use omicsweave_core::{Error, Locus, SequenceRead, Span, Strand};
use omicsweave_fasta::{
    Fai, FastaReader, TwoBitPiece, TwoBitReader, TwoBitRecord, TwoBitWriter, compress_fasta,
};
use proptest::prelude::*;

fn genome(n: usize, seed: u64) -> Vec<u8> {
    let mut x = seed | 1;
    (0..n)
        .map(|_| {
            x ^= x << 13;
            x ^= x >> 7;
            x ^= x << 17;
            b"ACGTacgtNRYKMSWBDHVn"[(x % 20) as usize]
        })
        .collect()
}

fn fasta_text(bases: &[u8], width: usize) -> Vec<u8> {
    let mut text = b">s\n".to_vec();
    for line in bases.chunks(width) {
        text.extend_from_slice(line);
        text.push(b'\n');
    }
    text
}

fn twobit(bases: &[u8]) -> Vec<u8> {
    let mut w = TwoBitWriter::create(Cursor::new(Vec::new()), vec!["s".into()], 0).unwrap();
    let piece = TwoBitPiece::encode(0, bases).unwrap();
    let mut record = TwoBitRecord::new();
    record.push(&piece).unwrap();
    w.add("s", &record, &mut &piece.packed[..]).unwrap();
    w.finish().unwrap().into_inner()
}

fn expected(stored: &[u8], spans: &[Span], strand: Strand) -> Vec<u8> {
    let mut out: Vec<u8> = spans
        .iter()
        .flat_map(|s| stored[s.start as usize..s.end as usize].iter().copied())
        .collect();
    if strand == Strand::Reverse {
        reverse_complement(&mut out);
    }
    out
}

fn spans_strategy(len: u32) -> impl Strategy<Value = (Vec<Span>, Strand)> {
    (
        prop::collection::vec((0u32..20_000, 0u32..4_000), 0..6),
        prop_oneof![
            Just(Strand::Forward),
            Just(Strand::Reverse),
            Just(Strand::Unstranded)
        ],
    )
        .prop_map(move |(parts, strand)| {
            let mut spans = Vec::new();
            let mut at = 0u32;
            for (gap, width) in parts {
                let start = (at + gap).min(len);
                let end = (start + width).min(len);
                spans.push(Span::new(start, end));
                at = end;
            }
            (spans, strand)
        })
}

const LEN: usize = 300_000;

proptest! {
    #![proptest_config(ProptestConfig::with_cases(48))]

    #[test]
    fn every_reader_concatenates_and_orients(
        seed in any::<u64>(),
        (spans, strand) in spans_strategy(LEN as u32),
        width in prop_oneof![Just(60usize), Just(61), Just(4096)],
    ) {
        let bases = genome(LEN, seed);
        let locus = Locus::new("s", &spans, strand).unwrap();
        let text = fasta_text(&bases, width);
        let (bgzf, fai, gzi) = compress_fasta(text.as_slice(), Vec::new(), 1).unwrap();
        let want = expected(&bases, &spans, strand);
        let mut out = vec![0u8; locus.len() as usize];

        let mut plain = FastaReader::new(&text[..], Fai::build(&text[..]).unwrap(), None);
        plain.sequence(&locus, &mut out).unwrap();
        prop_assert_eq!(&out, &want);

        let mut zipped = FastaReader::new(&bgzf[..], fai.clone(), Some(gzi.clone()));
        out.fill(0);
        zipped.sequence(&locus, &mut out).unwrap();
        prop_assert_eq!(&out, &want);

        let mut cached = FastaReader::new(&bgzf[..], fai, Some(gzi));
        cached.set_cache(Some(Arc::new(BlockCache::new(8 << 20, 64)))).unwrap();
        out.fill(0);
        cached.sequence(&locus, &mut out).unwrap();
        prop_assert_eq!(&out, &want);

        // 2bit keeps ACGT, N runs and soft-mask runs only; compare with its
        // own whole-sequence read, which the 2bit differential tests cover.
        let packed = twobit(&bases);
        let mut two = TwoBitReader::open(&packed[..], true).unwrap();
        let mut stored2 = vec![0u8; LEN];
        two.sequence(&Locus::from(omicsweave_core::Region::new("s", 0, LEN as u32)), &mut stored2)
            .unwrap();
        out.fill(0);
        two.sequence(&locus, &mut out).unwrap();
        prop_assert_eq!(&out, &expected(&stored2, &spans, strand));
    }
}

#[test]
fn complement_keeps_case_and_iupac() {
    let forward = b"ACGTUacgtuRYKMBVDHSWNrykmbvdhswn-*.";
    let mapped: Vec<u8> = forward.iter().map(|&b| complement(b)).collect();
    assert_eq!(&mapped, b"TGCAAtgcaaYRMKVBHDSWNyrmkvbhdswn-*.");
}

#[test]
fn loci_are_checked() {
    assert!(matches!(
        Locus::new(
            "s",
            &[Span::new(10, 20), Span::new(15, 30)],
            Strand::Forward
        ),
        Err(Error::Invalid { .. })
    ));
    assert!(matches!(
        Locus::new("s", &[Span::new(20, 10)], Strand::Forward),
        Err(Error::InvalidRegion { .. })
    ));
    let bases = genome(1000, 3);
    let text = fasta_text(&bases, 60);
    let mut reader = FastaReader::new(&text[..], Fai::build(&text[..]).unwrap(), None);
    let spans = [Span::new(0, 10), Span::new(990, 1001)];
    let locus = Locus::new("s", &spans, Strand::Reverse).unwrap();
    let mut out = vec![0u8; 21];
    assert!(matches!(
        reader.sequence(&locus, &mut out),
        Err(Error::RegionOutOfRange { .. })
    ));
}
