//! BAI and CSI wire layouts, coordinate bounds and failed builders.
use omicsweave_bgzf::{
    BinIndex, BinningIndex, Blocks, Chunk, IndexFormat, Logical, Reader, try_csi_levels,
    try_reg2bin,
};
use omicsweave_core::exec::Inline;
use std::io::Read;

fn u32le(out: &mut Vec<u8>, v: u32) {
    out.extend_from_slice(&v.to_le_bytes());
}

fn u64le(out: &mut Vec<u8>, v: u64) {
    out.extend_from_slice(&v.to_le_bytes());
}

fn layout() -> Blocks {
    Blocks {
        starts: vec![(0, 0), (100, 65_280)],
        end: 200,
        text_len: 65_281,
    }
}

#[test]
fn bai_one_record_matches_independent_byte_oracle() {
    let blocks = layout();
    let (start, end) = (Logical::new(0, 3), Logical::new(1, 5));
    let mut index = BinningIndex::new(IndexFormat::Bai, 1, start);
    index.push(0, 0, 1, end, true).unwrap();
    let actual = index.finish(end, &blocks).unwrap();
    let end = (100u64 << 16) | 5;
    let mut expected = b"BAI\x01".to_vec();
    for v in [1, 2, 4681, 1] {
        u32le(&mut expected, v);
    }
    for v in [3, end] {
        u64le(&mut expected, v);
    }
    for v in [37450, 2] {
        u32le(&mut expected, v);
    }
    for v in [3, end, 1, 0] {
        u64le(&mut expected, v);
    }
    u32le(&mut expected, 1);
    for v in [3, 0] {
        u64le(&mut expected, v);
    }
    assert_eq!(actual, expected);
}

#[test]
fn csi_one_record_matches_independent_decompressed_byte_oracle() {
    let blocks = layout();
    let (start, end) = (Logical::new(0, 3), Logical::new(1, 5));
    let mut index = BinningIndex::new(
        IndexFormat::Csi {
            min_shift: 14,
            levels: 5,
        },
        1,
        start,
    );
    index.push(0, 0, 1, end, true).unwrap();
    let actual = index.finish(end, &blocks).unwrap();
    let mut decoded = Vec::new();
    Reader::new(&actual[..], Box::new(Inline::new()))
        .read_to_end(&mut decoded)
        .unwrap();
    let end = (100u64 << 16) | 5;
    let mut expected = b"CSI\x01".to_vec();
    for v in [14, 5, 0, 1, 2, 4681] {
        u32le(&mut expected, v);
    }
    u64le(&mut expected, 3);
    u32le(&mut expected, 1);
    for v in [3, end] {
        u64le(&mut expected, v);
    }
    u32le(&mut expected, 37450);
    u64le(&mut expected, 0);
    u32le(&mut expected, 2);
    for v in [3, end, 1, 0, 0] {
        u64le(&mut expected, v);
    }
    assert_eq!(decoded, expected);
}

#[test]
fn repeated_wide_intervals_and_gaps_preserve_first_window_offsets() {
    let mut index = BinningIndex::new(IndexFormat::Bai, 1, Logical::new(0, 1));
    for (beg, end, at) in [
        (0, 32_768, 2),
        (0, 16_384, 3),
        (16_384, 49_152, 4),
        (81_920, 98_304, 5),
    ] {
        index.push(0, beg, end, Logical::new(0, at), true).unwrap();
    }
    let actual = index.finish(Logical::new(0, 5), &layout()).unwrap();
    // Linear index is the last counted array before n_no_coor. Holes are
    // backfilled with the next occupied window's first-record offset.
    let mut expected_tail = Vec::new();
    u32le(&mut expected_tail, 6);
    for v in [1, 1, 3, 4, 4, 4, 0] {
        u64le(&mut expected_tail, v);
    }
    assert!(actual.ends_with(&expected_tail));
}

#[test]
fn all_negative_tids_are_coordinate_free() {
    let mut index = BinningIndex::new(IndexFormat::Bai, 1, Logical::new(0, 0));
    index.push(-1, 0, 0, Logical::new(0, 1), false).unwrap();
    index.push(-2, 0, 0, Logical::new(0, 2), false).unwrap();
    assert!(index.finish(Logical::new(0, 2), &layout()).is_ok());
}

#[test]
fn checked_csi_parameters_cover_arithmetic_and_bin_id_boundaries() {
    assert!(try_csi_levels(u64::MAX, 14).is_err());
    assert!(try_csi_levels(0, u32::MAX).is_err());
    assert!(try_reg2bin(0, i64::MIN, 14, 5).is_err());
    assert!(try_reg2bin(0, 1, u32::MAX, u32::MAX).is_err());
    assert!(
        BinningIndex::try_new(
            IndexFormat::Csi {
                min_shift: 0,
                levels: 11
            },
            0,
            Logical(0)
        )
        .is_err()
    );
    assert!(
        BinningIndex::try_new(
            IndexFormat::Csi {
                min_shift: 34,
                levels: 10
            },
            0,
            Logical(0)
        )
        .is_err()
    );
    assert_eq!(try_reg2bin(0, 1, 33, 10).unwrap(), 153_391_689);
    let mut index = BinningIndex::try_new(
        IndexFormat::Csi {
            min_shift: 33,
            levels: 10,
        },
        1,
        Logical(0),
    )
    .unwrap();
    index.push(0, 0, 1, Logical(1), true).unwrap();
    assert!(index.finish(Logical(1), &layout()).is_ok());
    let mut invalid = BinningIndex::new(
        IndexFormat::Csi {
            min_shift: u32::MAX,
            levels: u32::MAX,
        },
        1,
        Logical(0),
    );
    invalid.amend_last(Logical(1));
    assert!(invalid.push(0, 0, 1, Logical(2), true).is_err());
    assert!(invalid.finish(Logical(2), &layout()).is_err());
}

#[test]
fn failed_index_rejects_later_records_and_decreasing_offsets() {
    let mut index = BinningIndex::new(IndexFormat::Bai, 2, Logical(0));
    assert!(index.push(1, 10, 9, Logical(1), true).is_err());
    assert!(index.push(1, 10, 11, Logical(2), true).is_err());
    assert!(index.finish(Logical(2), &layout()).is_err());
    let mut index = BinningIndex::new(IndexFormat::Bai, 1, Logical(2));
    assert!(index.push(0, 0, 1, Logical(1), true).is_err());
}

#[test]
fn checked_regions_and_amended_offsets_reject_invalid_boundaries() {
    assert!(try_reg2bin(i64::MAX - 1, i64::MAX, 0, 10).is_err());
    assert!(try_reg2bin(i64::MIN, i64::MIN, 14, 5).is_err());
    assert!(try_reg2bin(0, -1, 14, 5).is_err());
    assert!(try_reg2bin(10, 9, 14, 5).is_err());
    assert_eq!(try_reg2bin(-1, 0, 14, 5).unwrap(), 4680);
    assert_eq!(try_reg2bin(0, 1 << 29, 14, 5).unwrap(), 0);
    let mut index = BinningIndex::new(IndexFormat::Bai, 1, Logical(2));
    index.amend_last(Logical(1));
    assert!(index.push(0, 0, 1, Logical(3), true).is_err());
    assert!(index.finish(Logical(3), &Blocks::default()).is_err());
    let format = IndexFormat::Csi {
        min_shift: 63,
        levels: 0,
    };
    let mut index = BinningIndex::new(format, 1, Logical(0));
    assert!(
        index
            .push(0, i64::MAX, i64::MAX - 1, Logical(1), true)
            .is_err()
    );
    let mut index = BinningIndex::new(format, 1, Logical(0));
    index.push(0, i64::MAX, i64::MAX, Logical(1), true).unwrap();
    assert!(index.push(0, 0, 1, Logical(2), true).is_err());
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
