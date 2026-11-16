//! Layout bytes, bounded reads and recoverable decode failures.
use crate::support::htslib;
use omicsweave_core::{Locus, Region, SequenceRead};
use omicsweave_fasta::{
    BgzfWriter, BlockHeader, Error, Fai, FaiEntry, FastaReader, compress_fasta,
};
use std::io::Write;

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
fn every_gap_byte_matches_explicit_and_htslib_oracles() {
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
    check_htslib_line_gaps();
}

#[test]
fn plain_reads_exclude_unrequested_bases_and_trailing_layout_bytes() {
    let source = crate::support::Source::new(b">s\nACGT\x0B\nacgt\x0B\nNNNN\x0B\n".to_vec());
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
            .sequence(&Locus::from(Region::new("s", start, end)), &mut out)
            .unwrap();
        assert_eq!(&out, &b"ACGTacgtNNNN"[start as usize..end as usize]);
        assert_eq!(source.requests.borrow().as_slice(), &[expected_read]);
    }
    source.requests.borrow_mut().clear();
    reader
        .sequence(&Locus::from(Region::new("s", 4, 4)), &mut [])
        .unwrap();
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
                .sequence(&Locus::from(Region::new("s", 0, 4)), &mut guarded[1..5])
                .unwrap();
            assert_eq!(&guarded[1..5], b"ACGT");
            let error = reader
                .sequence(&Locus::from(Region::new("s", 0, 12)), &mut guarded[1..13])
                .unwrap_err();
            if corrupt_crc {
                assert!(matches!(error, Error::Inflate(_)));
            } else {
                assert!(matches!(error, Error::Invalid { .. }));
            }
            assert_eq!((guarded[0], guarded[13]), (0xAA, 0xAA));
            if !corrupt_crc {
                reader
                    .sequence(&Locus::from(Region::new("s", 8, 12)), &mut guarded[1..5])
                    .unwrap();
                assert_eq!(&guarded[1..5], b"NNNN");
            }
        }
    }
}

fn check_htslib_line_gaps() {
    for whitespace in [b' ', b'\t', 0x0B, 0x0C, b'\r'] {
        let mut text = b">s\n".to_vec();
        for line in [b"AC".as_slice(), b"GT", b"N"] {
            text.extend_from_slice(line);
            text.extend_from_slice(&[whitespace, b'\n']);
        }
        let dir = tempfile::tempdir().unwrap();
        let path = dir.path().join("spaces.fa");
        std::fs::write(&path, &text).unwrap();
        htslib::quiet();
        assert!(htslib::build_index(&path));
        let index_bytes = std::fs::read(dir.path().join("spaces.fa.fai")).unwrap();
        let fai = Fai::parse(&index_bytes).unwrap();
        let reference = htslib::Faidx::load(&path).unwrap();
        let (compressed, built_fai, gzi) = compress_fasta(text.as_slice(), Vec::new(), 1).unwrap();
        assert_eq!(built_fai, fai);
        let mut plain = FastaReader::new(text.as_slice(), fai.clone(), None);
        let mut bgzf = FastaReader::new(compressed.as_slice(), fai, Some(gzi));
        for (start, end) in [(0, 1), (0, 3), (0, 5), (1, 4), (4, 5), (0, 5)] {
            let expected = reference.fetch("s", start, end).unwrap();
            assert_eq!(&expected, &b"ACGTN"[start as usize..end as usize]);
            for reader in [&mut plain, &mut bgzf] {
                let mut guarded = vec![0xAA; expected.len() + 2];
                reader
                    .sequence(
                        &Locus::from(Region::new("s", start, end)),
                        &mut guarded[1..expected.len() + 1],
                    )
                    .unwrap();
                assert_eq!(&guarded[1..expected.len() + 1], expected.as_slice());
                assert_eq!((guarded[0], guarded[expected.len() + 1]), (0xAA, 0xAA));
            }
        }
    }
}

#[test]
fn line_endings_decode_at_every_member_split_and_window_and_recover() {
    let expected = b"ACGTacgtNNNN";
    for (text, offset) in [
        (b">s\nACGT\x0B\nacgt\x0B\nNNNN\x0B\n".as_slice(), 3),
        (b">s\r\nACGT\r\nacgt\r\nNNNN\r\n".as_slice(), 4),
    ] {
        let fai = Fai::from_entries([FaiEntry { offset, ..entry() }]);
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
                            &Locus::from(Region::new("s", start as u32, end as u32)),
                            &mut guarded[1..end - start + 1],
                        )
                        .unwrap();
                    assert_eq!(
                        &guarded[1..end - start + 1],
                        &expected[start..end],
                        "text={text:?} split={split} range={start}..{end}"
                    );
                    assert_eq!((guarded[0], guarded[end - start + 1]), (0xAA, 0xAA));
                }
            }
            let mut guarded = [0xAA; 14];
            assert!(
                reader
                    .sequence(&Locus::from(Region::new("s", 0, 13)), &mut guarded[1..13])
                    .is_err()
            );
            assert_eq!(guarded, [0xAA; 14]);
            reader
                .sequence(&Locus::from(Region::new("s", 3, 5)), &mut guarded[1..3])
                .unwrap();
            assert_eq!(&guarded[1..3], b"Ta");
            reader
                .sequence(&Locus::from(Region::new("s", 0, 12)), &mut guarded[1..13])
                .unwrap();
            assert_eq!(&guarded, b"\xAAACGTacgtNNNN\xAA");
        }
    }
}
