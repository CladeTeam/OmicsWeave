//! Original block metadata and an independent per-base coverage oracle.
mod common;

use common::{OpenTwoBit, TwoBitSeq, lib2bit, write_twobit};
use omicsweave_core::{Locus, Region, SequenceRead};
use omicsweave_fasta::{Block, TwoBitReader};

fn expected(seq: &TwoBitSeq, start: u32, end: u32, mask: bool) -> Vec<u8> {
    (start..end)
        .map(|p| {
            let covers = |&(s, n): &(u32, u32)| s <= p && p < s + n;
            if seq.n_blocks.iter().any(covers) {
                b'N'
            } else if mask && seq.mask_blocks.iter().any(covers) {
                seq.bases[p as usize].to_ascii_lowercase()
            } else {
                seq.bases[p as usize]
            }
        })
        .collect()
}

#[test]
fn empty_and_adjacent_blocks_preserve_wire_metadata_and_lib2bit_bases() {
    let seq = TwoBitSeq {
        name: "s".into(),
        bases: b"ACGTACGTACGTACGTACGTACGTACGTACGT".to_vec(),
        n_blocks: vec![
            (0, 0),
            (0, 1),
            (1, 0),
            (1, 2),
            (3, 2),
            (8, 0),
            (10, 1),
            (32, 0),
        ],
        mask_blocks: vec![
            (0, 0),
            (1, 1),
            (2, 0),
            (2, 5),
            (7, 1),
            (12, 0),
            (13, 2),
            (32, 0),
        ],
    };
    let dir = tempfile::tempdir().unwrap();
    for version in [0, 1] {
        let bytes = write_twobit(std::slice::from_ref(&seq), version);
        let opened = OpenTwoBit::open(&bytes);
        let record = &opened.records[0];
        for (raw, wanted) in [
            (record.n_blocks(), &seq.n_blocks),
            (record.mask_blocks(), &seq.mask_blocks),
        ] {
            assert_eq!(
                raw,
                wanted
                    .iter()
                    .map(|&(start, len)| Block { start, len })
                    .collect::<Vec<_>>()
            );
        }
        let path = dir.path().join(format!("runs-{version}.2bit"));
        std::fs::write(&path, &bytes).unwrap();
        for mask in [false, true] {
            let oracle = lib2bit::TwoBit::open(&path, mask).unwrap();
            let mut reader = TwoBitReader::open(bytes.as_slice(), mask).unwrap();
            for start in 0..=seq.bases.len() as u32 {
                for end in start..=seq.bases.len() as u32 {
                    let wanted = expected(&seq, start, end, mask);
                    let mut guarded = vec![0xA5; wanted.len() + 2];
                    reader
                        .sequence(
                            &Locus::from(Region::new("s", start, end)),
                            &mut guarded[1..wanted.len() + 1],
                        )
                        .unwrap();
                    assert_eq!(&guarded[1..wanted.len() + 1], wanted);
                    assert_eq!((guarded[0], guarded[wanted.len() + 1]), (0xA5, 0xA5));
                    if start != end {
                        assert_eq!(oracle.sequence("s", start, end).unwrap(), wanted);
                    }
                    if !wanted.is_empty() && wanted.iter().all(|&c| c == b'N') {
                        assert!(record.packed_range(start, end).unwrap().is_empty());
                    }
                }
            }
            assert!(
                reader
                    .sequence(&Locus::from(Region::new("s", 0, 33)), &mut [0; 33])
                    .is_err()
            );
            let mut out = [0; 6];
            reader
                .sequence(&Locus::from(Region::new("s", 0, 6)), &mut out)
                .unwrap();
            assert_eq!(out.as_slice(), expected(&seq, 0, 6, mask));
        }
    }
}

#[test]
fn short_queries_cross_many_empty_blocks_and_adjacent_n_runs() {
    for all_n in [false, true] {
        let mut n_blocks = vec![(0, 1)];
        n_blocks.extend(std::iter::repeat_n((1, 0), 4096));
        if all_n {
            n_blocks.push((1, 1));
        }
        let mut mask_blocks = vec![(0, 1)];
        mask_blocks.extend(std::iter::repeat_n((1, 0), 4096));
        mask_blocks.push((1, 2));
        let seq = TwoBitSeq {
            name: "s".into(),
            bases: vec![b'A'; 8192],
            n_blocks,
            mask_blocks,
        };
        let bytes = write_twobit(std::slice::from_ref(&seq), 0);
        let mut reader = TwoBitReader::open(bytes.as_slice(), true).unwrap();
        for (start, end) in [(0, 2), (0, 3), (1, 3), (2, 4), (0, 2)] {
            let mut out = vec![0; (end - start) as usize];
            reader
                .sequence(&Locus::from(Region::new("s", start, end)), &mut out)
                .unwrap();
            assert_eq!(out, expected(&seq, start, end, true));
        }
    }
}
