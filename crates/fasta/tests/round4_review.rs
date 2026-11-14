//! Independent byte fixtures and sequence oracles for format boundary checks.

use std::io::Write;

use omicsweave_core::{Region, SequenceRead};
use omicsweave_fasta::{
    BgzfWriter, BlockHeader, Directory, Error, Fai, FaiEntry, FastaReader, Parsed, TwoBitHeader,
};

fn entry() -> FaiEntry {
    FaiEntry {
        name: "s".into(),
        len: 12,
        offset: 3,
        line_bases: 4,
        line_width: 6,
    }
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
    }
    // A version determines how offsets are interpreted even in an empty file.
    for version in [2, 3, 255, u32::MAX] {
        for count in [0, 1, u32::MAX] {
            assert!(Directory::parse(&TwoBitHeader { version, count }, &[]).is_err());
        }
    }
}

#[test]
fn every_gap_byte_is_classified_against_an_explicit_c_locale_oracle() {
    let metadata = entry();
    for byte in 0..=u8::MAX {
        let mut raw = b"ACGT\0\nacgt\0\nNNNN".to_vec();
        raw[4] = byte;
        raw[10] = byte;
        let allowed = [9, 10, 11, 12, 13, 32].contains(&byte);
        // One crossing uses the short path; two crossings use streaming lines.
        for length in [5, 12] {
            let mut guarded = vec![0xAA; length + 2];
            let result = metadata.decode(0, &raw, &mut guarded[1..length + 1]);
            assert_eq!(result.is_ok(), allowed, "gap byte {byte}, length {length}");
            if allowed {
                assert_eq!(&guarded[1..length + 1], &b"ACGTacgtNNNN"[..length]);
            }
            assert_eq!((guarded[0], guarded[length + 1]), (0xAA, 0xAA));
        }
    }
}

#[test]
fn vertical_tab_gaps_decode_at_every_member_split_and_window() {
    let text = b">s\nACGT\x0B\nacgt\x0B\nNNNN\x0B\n";
    let expected = b"ACGTacgtNNNN";
    let fai = Fai::from_entries([entry()]);
    // Coordinates and expected output come from the sequence constant, not a
    // second invocation of the FASTA builder or text-range planner.
    for split in 1..text.len() {
        let mut writer = BgzfWriter::new(Vec::new());
        writer.write_all(&text[..split]).unwrap();
        writer.flush().unwrap();
        writer.write_all(&text[split..]).unwrap();
        let (compressed, gzi) = writer.finish().unwrap();
        let mut reader = FastaReader::new(compressed.as_slice(), fai.clone(), Some(gzi));
        for start in 0..=expected.len() {
            for end in start..=expected.len() {
                let mut guarded = vec![0xAA; end - start + 2];
                reader
                    .sequence(
                        Region::new("s", start as u32, end as u32),
                        &mut guarded[1..end - start + 1],
                    )
                    .unwrap();
                assert_eq!(&guarded[1..end - start + 1], &expected[start..end]);
                assert_eq!((guarded[0], guarded[end - start + 1]), (0xAA, 0xAA));
            }
        }
    }
}

#[test]
fn plain_reads_exclude_unrequested_bases_and_trailing_layout_bytes() {
    use std::cell::RefCell;

    use omicsweave_core::{ReadAt, Result};

    struct Traced<'a> {
        bytes: &'a [u8],
        requests: RefCell<Vec<(u64, usize)>>,
    }
    impl ReadAt for Traced<'_> {
        fn len(&self) -> u64 {
            self.bytes.len() as u64
        }

        fn read_exact_at(&self, out: &mut [u8], at: u64) -> Result<()> {
            self.requests.borrow_mut().push((at, out.len()));
            self.bytes.read_exact_at(out, at)
        }
    }

    let source = Traced {
        bytes: b">s\nACGT\x0B\nacgt\x0B\nNNNN\x0B\n",
        requests: RefCell::new(Vec::new()),
    };
    let mut reader = FastaReader::new(&source, Fai::from_entries([entry()]), None);
    for (start, end, expected_read) in [
        (0, 4, (3, 4)),
        (3, 5, (6, 4)),
        (4, 8, (9, 4)),
        (0, 12, (3, 16)),
    ] {
        source.requests.borrow_mut().clear();
        let mut out = vec![0; (end - start) as usize];
        reader
            .sequence(Region::new("s", start, end), &mut out)
            .unwrap();
        assert_eq!(&out, &b"ACGTacgtNNNN"[start as usize..end as usize]);
        assert_eq!(source.requests.borrow().as_slice(), &[expected_read]);
    }
    source.requests.borrow_mut().clear();
    reader.sequence(Region::new("s", 4, 4), &mut []).unwrap();
    assert!(source.requests.borrow().is_empty());
}

#[test]
fn decoding_and_crc_failures_do_not_poison_later_queries() {
    let good = b">s\nACGT\x0B\nacgt\x0B\nNNNN\x0B\n";
    for corrupt_crc in [false, true] {
        let mut text = *good;
        if !corrupt_crc {
            text[13] = b'X';
        }
        let mut writer = BgzfWriter::new(Vec::new());
        writer.write_all(&text[..9]).unwrap();
        writer.flush().unwrap();
        writer.write_all(&text[9..]).unwrap();
        let (mut compressed, gzi) = writer.finish().unwrap();
        if corrupt_crc {
            let second = gzi.starts()[1].0 as usize;
            let header = BlockHeader::parse(&compressed[second..]).unwrap();
            compressed[second + header.block_len - 8] ^= 1;
        }
        let mut reader = FastaReader::new(
            compressed.as_slice(),
            Fai::from_entries([entry()]),
            Some(gzi),
        );
        let mut guarded = [0xAA; 14];
        for _ in 0..2 {
            reader
                .sequence(Region::new("s", 0, 4), &mut guarded[1..5])
                .unwrap();
            assert_eq!(&guarded[1..5], b"ACGT");
            let error = reader
                .sequence(Region::new("s", 0, 12), &mut guarded[1..13])
                .unwrap_err();
            if corrupt_crc {
                assert!(matches!(error, Error::Inflate(_)));
            } else {
                assert!(matches!(error, Error::Invalid { .. }));
            }
            assert_eq!((guarded[0], guarded[13]), (0xAA, 0xAA));
            if !corrupt_crc {
                reader
                    .sequence(Region::new("s", 8, 12), &mut guarded[1..5])
                    .unwrap();
                assert_eq!(&guarded[1..5], b"NNNN");
            }
        }
    }
}
