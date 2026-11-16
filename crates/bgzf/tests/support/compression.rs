//! Deterministic payloads and an independent gzip decoder.
use omicsweave_bgzf::{EOF_BLOCK, MAX_BLOCK_LEN};

pub fn text(len: usize, random: bool) -> Vec<u8> {
    let mut state = 0x1234_5678u32;
    (0..len)
        .map(|i| {
            if random {
                state ^= state << 13;
                state ^= state >> 17;
                state ^= state << 5;
                state as u8
            } else {
                b"ACGT\n"[i % 5]
            }
        })
        .collect()
}

pub fn decode_gzip_block(raw: &[u8], expected: &[u8]) {
    assert!((28..=MAX_BLOCK_LEN).contains(&raw.len()));
    assert_eq!(&raw[..16], &EOF_BLOCK[..16]);
    assert_eq!(
        u16::from_le_bytes([raw[16], raw[17]]) as usize + 1,
        raw.len()
    );
    let mut decoded = vec![0; expected.len()];
    let n = libdeflater::Decompressor::new()
        .gzip_decompress(raw, &mut decoded)
        .unwrap();
    assert_eq!(n, expected.len());
    assert_eq!(decoded, expected);
}
