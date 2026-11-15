//! 2bit reads against lib2bit, on random files from the test writer.

mod common;

use common::{OpenTwoBit, Rng, lib2bit, random_twobit, write_twobit};
use omicsweave_core::{Locus, Region, SequenceRead};
use omicsweave_fasta::TwoBitReader;
use proptest::prelude::*;

proptest! {
    #![proptest_config(ProptestConfig { cases: 64, ..ProptestConfig::default() })]

    #[test]
    fn matches_lib2bit(seed in any::<u64>(), version in 0u32..=1) {
        let mut rng = Rng::new(seed);
        let seqs = random_twobit(&mut rng);
        let bytes = write_twobit(&seqs, version);
        let dir = tempfile::tempdir().unwrap();
        let path = dir.path().join("t.2bit");
        std::fs::write(&path, &bytes).unwrap();

        let ours = OpenTwoBit::open(&bytes);
        prop_assert!(ours.rounds.iter().all(|&r| r <= 3), "rounds {:?}", ours.rounds);
        let names: Vec<&str> = ours.directory.entries().iter().map(|(n, _)| n.as_str()).collect();
        let expected: Vec<&str> = seqs.iter().map(|s| s.name.as_str()).collect();
        prop_assert_eq!(names, expected);

        for mask in [false, true] {
            let theirs = lib2bit::TwoBit::open(&path, mask).unwrap();
            let mut reader = TwoBitReader::open(&bytes[..], mask).unwrap();
            for (i, s) in seqs.iter().enumerate() {
                let len = ours.records[i].len();
                prop_assert_eq!(len, theirs.len(&s.name));
                prop_assert_eq!(reader.chrom_len(&s.name).unwrap(), len);
                if len == 0 {
                    continue;
                }
                let mut ranges = vec![(0, len)];
                for _ in 0..30 {
                    let start = rng.coord(0, len - 1);
                    ranges.push((start, rng.coord(start + 1, len)));
                }
                for (start, end) in ranges {
                    let expected = theirs.sequence(&s.name, start, end).unwrap();
                    prop_assert_eq!(
                        &ours.fetch(i, start, end, mask), &expected,
                        "{}:{}-{} mask={}", s.name, start, end, mask
                    );
                    let mut out = vec![0u8; (end - start) as usize];
                    reader.sequence(&Locus::from(Region::new(&s.name, start, end)), &mut out).unwrap();
                    prop_assert_eq!(&out, &expected);
                }
            }
        }
    }
}

#[test]
fn packed_bytes_decode_at_every_alignment_without_touching_guards() {
    use omicsweave_fasta::{Parsed, Record};
    let packed: Vec<_> = (0..=255u8).collect();
    let mut header = Vec::new();
    for n in [1024u32, 0, 0, 0] {
        header.extend(n.to_le_bytes());
    }
    let Parsed::Done(record) = Record::parse(&header, 0).unwrap() else {
        panic!("complete header")
    };
    for start in 0..1024u32 {
        for len in [0usize, 1, 2, 3, 4, 5, 15, 16, 17, 63, 64, 65] {
            if start as usize + len > 1024 {
                continue;
            }
            let mut out = vec![0xAA; len + 2];
            record
                .decode(
                    start,
                    &packed[start as usize / 4..],
                    &mut out[1..len + 1],
                    false,
                )
                .unwrap();
            assert_eq!(out[0], 0xAA);
            assert_eq!(out[len + 1], 0xAA);
            for (i, &base) in out[1..len + 1].iter().enumerate() {
                let pos = start as usize + i;
                assert_eq!(
                    base,
                    b"TCAG"[((packed[pos / 4] >> (6 - 2 * (pos % 4))) & 3) as usize]
                );
            }
        }
    }
}
