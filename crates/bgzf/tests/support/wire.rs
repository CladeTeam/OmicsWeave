//! BGZF members encoded without the crate's compressor.
use omicsweave_bgzf::{EOF_BLOCK, MAX_BLOCK_LEN};

pub fn crc32(text: &[u8]) -> u32 {
    let mut crc = !0u32;
    for &byte in text {
        crc ^= u32::from(byte);
        for _ in 0..8 {
            crc = (crc >> 1) ^ (0xedb8_8320 & 0u32.wrapping_sub(crc & 1));
        }
    }
    !crc
}

pub fn stored(text: &[u8]) -> Vec<u8> {
    stored_with_extra(text, false)
}

// The extra field exercises headers whose BC subfield is not first.
pub fn stored_with_extra(text: &[u8], extended: bool) -> Vec<u8> {
    assert!(text.len() <= u16::MAX as usize);
    let mut raw = EOF_BLOCK[..12].to_vec();
    if extended {
        raw.extend_from_slice(&[b'X', b'Y', 12, 0]);
        raw.extend_from_slice(&[0xab; 12]);
    }
    let at = raw.len();
    raw.extend_from_slice(&[b'B', b'C', 2, 0, 0, 0]);
    let xlen = (raw.len() - 12) as u16;
    raw[10..12].copy_from_slice(&xlen.to_le_bytes());
    let n = text.len() as u16;
    raw.push(1);
    raw.extend_from_slice(&n.to_le_bytes());
    raw.extend_from_slice(&(!n).to_le_bytes());
    raw.extend_from_slice(text);
    raw.extend_from_slice(&crc32(text).to_le_bytes());
    raw.extend_from_slice(&(text.len() as u32).to_le_bytes());
    assert!(raw.len() <= MAX_BLOCK_LEN);
    let bsize = (raw.len() - 1) as u16;
    raw[at + 4..at + 6].copy_from_slice(&bsize.to_le_bytes());
    raw
}

// Raw DEFLATE and CRC of 65,536 zero bytes from independent zlib encoding.
pub fn full_zero_block() -> Vec<u8> {
    let mut raw = EOF_BLOCK[..18].to_vec();
    raw.extend_from_slice(&[237, 193, 1, 1, 0, 0, 0, 128, 144, 254, 175, 238, 8, 10]);
    raw.extend_from_slice(&[0; 63]);
    raw.push(0x6a);
    raw.extend_from_slice(&0xd797_8eebu32.to_le_bytes());
    raw.extend_from_slice(&65_536u32.to_le_bytes());
    let bsize = (raw.len() - 1) as u16;
    raw[16..18].copy_from_slice(&bsize.to_le_bytes());
    raw
}
