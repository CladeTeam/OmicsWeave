//! 2bit reads against lib2bit, on random files from the test writer.

use crate::support::Source;
use crate::support::{OpenTwoBit, Rng, TwoBitSeq, lib2bit, random_twobit, write_twobit};
use omicsweave_core::{Locus, Region, SequenceRead};
use omicsweave_fasta::{Directory, Error, Parsed, Record, TwoBitHeader, TwoBitReader};
use proptest::prelude::*;

proptest! {
    #![proptest_config(ProptestConfig { cases: 64, failure_persistence: Some(Box::new(proptest::test_runner::FileFailurePersistence::Direct(concat!(env!("CARGO_MANIFEST_DIR"), "/tests/cases/twobit.proptest-regressions")))), ..ProptestConfig::default() })]

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

#[test]
fn twobit_header_and_records() {
    let mut h = 0x1A41_2743u32.to_be_bytes().to_vec();
    h.extend_from_slice(&[0u8; 12]);
    assert!(matches!(TwoBitHeader::parse(&h), Err(Error::BigEndian)));

    // 8 bases, one N block that runs past the end.
    let mut r = Vec::new();
    for v in [8u32, 1, 6, 4, 0, 0] {
        r.extend_from_slice(&v.to_le_bytes());
    }
    assert_eq!(Record::parse(&r[..4], 0).unwrap(), Parsed::Need(8));
    assert_eq!(Record::parse(&r[..8], 0).unwrap(), Parsed::Need(20));
    assert_eq!(Record::parse(&r[..20], 0).unwrap(), Parsed::Need(24));
    assert!(matches!(Record::parse(&r, 0), Err(Error::Invalid { .. })));
}

#[test]
fn directory_versions_and_prefix_needs_match_wire_widths() {
    for version in [0, 1] {
        let bytes: &[u8] = if version == 0 {
            &[1, b's', 0x78, 0x56, 0x34, 0x12]
        } else {
            &[1, b's', 0x78, 0x56, 0x34, 0x12, 0xEF, 0xCD, 0xAB, 0x90]
        };
        let header = TwoBitHeader { version, count: 1 };
        for prefix in 0..bytes.len() {
            let need = if prefix == 0 { 1 } else { bytes.len() };
            assert_eq!(
                Directory::parse(&header, &bytes[..prefix]).unwrap(),
                Parsed::Need(need)
            );
        }
        let Parsed::Done(directory) = Directory::parse(&header, bytes).unwrap() else {
            panic!("complete wire directory");
        };
        let expected = if version == 0 {
            0x1234_5678
        } else {
            0x90AB_CDEF_1234_5678
        };
        assert_eq!(directory.entries(), &[("s".into(), expected)]);
        assert_eq!(directory.offset("s"), Some(expected));
    }
    // A version determines how offsets are interpreted even in an empty file.
    for version in [2, 3, 255, u32::MAX] {
        for count in [0, 1, u32::MAX] {
            assert!(Directory::parse(&TwoBitHeader { version, count }, &[]).is_err());
        }
    }
}

#[test]
fn minimal_record_headers_match_lib2bit_across_directory_read_boundary() {
    for version in [0, 1] {
        let seqs: Vec<_> = (0..8)
            .map(|i| TwoBitSeq {
                name: format!("s{i}"),
                bases: b"TCAG".repeat(4096),
                n_blocks: Vec::new(),
                mask_blocks: Vec::new(),
            })
            .collect();
        let bytes = write_twobit(&seqs, version);
        let dir = tempfile::tempdir().unwrap();
        let path = dir.path().join("records.2bit");
        std::fs::write(&path, &bytes).unwrap();
        let reference = lib2bit::TwoBit::open(&path, true).unwrap();
        let source = Source::new(bytes.clone());
        let mut reader = TwoBitReader::open(&source, true).unwrap();
        // The directory bound includes s0; each other record needs one
        // complete minimal header, with no extra packed-base prefetch.
        assert_eq!(source.calls.get(), 2 + seqs.len() - 1);
        assert_eq!(
            source.total.get(),
            16 + 264 * seqs.len() + 16 * (seqs.len() - 1)
        );
        for seq in &seqs {
            for (start, end) in [(0, 17), (123, 1024), (16380, 16384)] {
                let mut out = vec![0; (end - start) as usize];
                reader
                    .sequence(&Locus::from(Region::new(&seq.name, start, end)), &mut out)
                    .unwrap();
                assert_eq!(out, reference.sequence(&seq.name, start, end).unwrap());
            }
        }
    }
}

#[test]
fn staged_record_rejects_count_and_offset_overflow_before_decoding() {
    for fields in [[1u32, 2, 0, 0], [1, 0, 2, 0]] {
        let bytes: Vec<_> = fields.into_iter().flat_map(u32::to_le_bytes).collect();
        assert!(Record::parse(&bytes, 0).is_err(), "{fields:?}");
    }
    let header: Vec<_> = [8u32, 0, 0, 0]
        .into_iter()
        .flat_map(u32::to_le_bytes)
        .collect();
    assert!(Record::parse(&header, u64::MAX - 16).is_err());
    let Parsed::Done(record) = Record::parse(&header, 0).unwrap() else {
        panic!("complete header")
    };
    for mask in [false, true] {
        let mut out = [0xA5; 10];
        assert!(record.decode(0, &[0], &mut out[1..9], mask).is_err());
        assert_eq!(out, [0xA5; 10]);
        record
            .decode(0, &[0x9C, 0x9C], &mut out[1..9], mask)
            .unwrap();
        assert_eq!(&out, b"\xA5ACGTACGT\xA5");
    }
}
