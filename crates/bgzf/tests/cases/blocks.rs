//! Block framing, inflation and independent gzip compatibility.
use crate::support::compression::{decode_gzip_block, text};
use crate::support::wire::{crc32, full_zero_block, stored};
use omicsweave_bgzf::{
    BLOCK_TEXT_LEN, EOF_BLOCK, Inflater, check_level, deflate_block, inflate_block_owned,
    read_block,
};

#[test]
fn stored_deflate_oracle_and_crc_rejection() {
    assert_eq!(crc32(b"123456789"), 0xcbf4_3926);
    let raw = stored(b"123456789");
    let mut out = [0; 9];
    assert_eq!(
        Inflater::new().inflate_block(&raw, &mut out).unwrap(),
        (40, 9)
    );
    assert_eq!(&out, b"123456789");
    let mut corrupt = raw.clone();
    let at = corrupt.len() - 8;
    corrupt[at] ^= 1;
    assert!(Inflater::new().inflate_block(&corrupt, &mut out).is_err());
    let mut joined = raw.clone();
    joined.extend_from_slice(&EOF_BLOCK);
    let mut input = &joined[..];
    let mut block = Vec::new();
    assert!(read_block(&mut input, &mut block).unwrap());
    assert_eq!(block, raw);
    assert!(read_block(&mut input, &mut block).unwrap());
    assert_eq!(block, EOF_BLOCK);
    assert!(!read_block(&mut input, &mut block).unwrap());
}

#[test]
fn owned_inflate_validates_stored_empty_small_and_large_payloads() {
    for len in [0, 1, 37, 65_280] {
        let text: Vec<u8> = (0..len).map(|i| (i % 251) as u8).collect();
        let raw = stored(&text);
        assert_eq!(inflate_block_owned(&raw).unwrap(), text);
        let mut oversized = raw.clone();
        let last = oversized.len();
        oversized[last - 4..].copy_from_slice(&65_537u32.to_le_bytes());
        assert!(inflate_block_owned(&oversized).is_err());
        let mut undersized = raw.clone();
        undersized[last - 4..].copy_from_slice(&1u32.to_le_bytes());
        if len != 1 {
            assert!(inflate_block_owned(&undersized).is_err());
        }
        if len != 0 {
            undersized[last - 4..].copy_from_slice(&0u32.to_le_bytes());
            assert!(inflate_block_owned(&undersized).is_err());
        }
        let mut trailing = raw;
        trailing.push(0);
        assert!(inflate_block_owned(&trailing).is_err());
    }
    assert!(inflate_block_owned(&EOF_BLOCK).unwrap().is_empty());
    assert_eq!(
        inflate_block_owned(&full_zero_block()).unwrap(),
        vec![0; 65_536]
    );
}

#[test]
fn compression_levels_sizes_and_append_preserve_gzip_members() {
    for level in 0..=12 {
        let mut compressor = libdeflater::Compressor::new(check_level(level).unwrap());
        for len in [0, 1, 64, 4096, BLOCK_TEXT_LEN] {
            for random in [false, true] {
                let expected = text(len, random);
                let mut raw = b"preserved-prefix".to_vec();
                let base = raw.len();
                deflate_block(&mut compressor, &expected, &mut raw).unwrap();
                assert_eq!(&raw[..base], b"preserved-prefix");
                decode_gzip_block(&raw[base..], &expected);
            }
        }
    }
    let mut compressor = libdeflater::Compressor::new(check_level(6).unwrap());
    let mut raw = b"preserved-prefix".to_vec();
    assert!(deflate_block(&mut compressor, &text(BLOCK_TEXT_LEN + 1, false), &mut raw).is_err());
    assert_eq!(raw, b"preserved-prefix");
}

#[test]
fn inflate_respects_output_capacity_and_can_be_reused_after_rejection() {
    let raw = stored(b"payload");
    let mut inflater = Inflater::new();
    for capacity in 0..=9 {
        let mut out = [0xa5; 13];
        let result = inflater.inflate_block(&raw, &mut out[2..2 + capacity]);
        assert_eq!(&out[..2], &[0xa5; 2], "capacity={capacity}");
        assert!(out[2 + capacity..].iter().all(|&b| b == 0xa5));
        if capacity < 7 {
            assert!(
                matches!(result, Err(omicsweave_core::Error::Truncated {
                what: "bgzf output", need: 7, have, ..
            }) if have == capacity as u64),
                "capacity={capacity}: {result:?}"
            );
            assert!(out.iter().all(|&b| b == 0xa5));
        } else {
            assert_eq!(result.unwrap(), (raw.len(), 7));
            assert_eq!(&out[2..9], b"payload");
            assert!(out[9..].iter().all(|&b| b == 0xa5));
        }
    }
    let mut out = [0; 7];
    assert_eq!(
        inflater.inflate_block(&raw, &mut out).unwrap(),
        (raw.len(), 7)
    );
    assert_eq!(&out, b"payload");
}

#[test]
fn multi_block_inflation_counts_members_and_preserves_output_guards() {
    let mut raw = stored(b"ab");
    raw.extend_from_slice(&EOF_BLOCK);
    raw.extend_from_slice(&stored(b"cde"));
    raw.extend_from_slice(b"ignored trailing bytes");
    let mut inflater = Inflater::new();
    for (blocks, expected) in [(0, b"".as_slice()), (1, b"ab"), (2, b"ab"), (3, b"abcde")] {
        let mut out = [0xa5; 9];
        let count = inflater
            .inflate_blocks(&raw, blocks, &mut out[2..7])
            .unwrap();
        assert_eq!(count, expected.len(), "members={blocks}");
        assert_eq!(&out[2..2 + count], expected, "members={blocks}");
        assert_eq!(&out[..2], &[0xa5; 2]);
        assert!(out[2 + count..].iter().all(|&b| b == 0xa5));
    }
}

proptest::proptest! {
    #![proptest_config(proptest::test_runner::Config {
        cases: 64,
        failure_persistence: Some(Box::new(proptest::test_runner::FileFailurePersistence::Direct(
            concat!(env!("CARGO_MANIFEST_DIR"), "/tests/blocks.proptest-regressions")
        ))),
        ..proptest::test_runner::Config::default()
    })]

    #[test]
    fn structured_members_decode_in_order_with_guarded_output(
        members in proptest::collection::vec(
            proptest::collection::vec(proptest::prelude::any::<u8>(), 0..257), 0..8
        )
    ) {
        let raw: Vec<u8> = members.iter().flat_map(|member| stored(member)).collect();
        let expected = members.concat();
        let mut out = vec![0xa5; expected.len() + 4];
        let end = out.len() - 2;
        let count = Inflater::new().inflate_blocks(&raw, members.len(), &mut out[2..end]).unwrap();
        proptest::prop_assert_eq!(count, expected.len());
        proptest::prop_assert_eq!(&out[2..end], expected.as_slice());
        proptest::prop_assert_eq!(&out[..2], &[0xa5; 2]);
        proptest::prop_assert_eq!(&out[end..], &[0xa5; 2]);
    }
}
