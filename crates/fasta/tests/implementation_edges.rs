//! Independent format-oracle coverage for metadata and line-layout edges.

mod common;

use common::htslib;
use omicsweave_core::{Locus, Region, SequenceRead};
use omicsweave_fasta::{Directory, Fai, FastaReader, Parsed, TwoBitHeader, compress_fasta};

#[test]
fn directory_validates_public_header_versions_before_interpreting_offsets() {
    let name = b's';
    for version in [0, 1] {
        let offset = if version == 1 {
            0x1234_5678_9ABC_DEF0u64
        } else {
            0x1234_5678u64
        };
        let width = if version == 1 { 8 } else { 4 };
        let mut bytes = vec![1, name];
        bytes.extend_from_slice(&offset.to_le_bytes()[..width]);
        let header = TwoBitHeader { version, count: 1 };
        let Parsed::Done(directory) = Directory::parse(&header, &bytes).unwrap() else {
            panic!("complete directory");
        };
        assert_eq!(directory.offset("s"), Some(offset));
        for bad_version in [2, u32::MAX] {
            let header = TwoBitHeader {
                version: bad_version,
                count: 1,
            };
            assert!(Directory::parse(&header, &bytes).is_err());
            assert!(Directory::parse(&header, &[]).is_err());
        }
    }
}

#[test]
fn all_c_locale_line_gaps_match_htslib_in_plain_and_bgzf_reads() {
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
