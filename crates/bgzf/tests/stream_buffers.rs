use std::io::{self, Read};

use omicsweave_bgzf::{EOF_BLOCK, Reader, read_block};
use omicsweave_core::exec::{Inline, Pool};

fn stored(text: &[u8], extended: bool) -> Vec<u8> {
    let mut crc = !0u32;
    for &byte in text {
        crc ^= u32::from(byte);
        for _ in 0..8 {
            crc = (crc >> 1) ^ (0xedb8_8320 & 0u32.wrapping_sub(crc & 1));
        }
    }
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
    raw.extend_from_slice(&(!crc).to_le_bytes());
    raw.extend_from_slice(&(text.len() as u32).to_le_bytes());
    let bsize = (raw.len() - 1) as u16;
    raw[at + 4..at + 6].copy_from_slice(&bsize.to_le_bytes());
    raw
}

struct Pieces<'a> {
    input: &'a [u8],
    calls: usize,
}

impl Read for Pieces<'_> {
    fn read(&mut self, out: &mut [u8]) -> io::Result<usize> {
        self.calls += 1;
        if self.calls.is_multiple_of(5) {
            return Err(io::ErrorKind::Interrupted.into());
        }
        let n = out.len().min(3).min(self.input.len());
        out[..n].copy_from_slice(&self.input[..n]);
        self.input = &self.input[n..];
        Ok(n)
    }
}

#[test]
fn block_headers_and_bodies_survive_short_reads_with_reused_storage() {
    let blocks = [
        stored(b"one", false),
        stored(b"extended", true),
        EOF_BLOCK.to_vec(),
    ];
    let stream = blocks.concat();
    for capacity in [0, 17, 4096] {
        let mut source = Pieces {
            input: &stream,
            calls: 0,
        };
        let mut raw = Vec::with_capacity(capacity);
        for expected in &blocks {
            assert!(read_block(&mut source, &mut raw).unwrap());
            assert_eq!(&raw, expected);
        }
        assert!(!read_block(&mut source, &mut raw).unwrap());
        assert!(raw.is_empty());
    }
}

#[test]
fn consumed_text_storage_can_cross_empty_and_different_sized_members() {
    let large: Vec<_> = (0..60_000).map(|n| (n % 251) as u8).collect();
    let texts = [b"a".to_vec(), large, b"last".to_vec()];
    let mut raw = Vec::new();
    for (i, text) in texts.iter().enumerate() {
        raw.extend_from_slice(&stored(text, i % 2 != 0));
        raw.extend_from_slice(&EOF_BLOCK);
    }
    let expected = texts.concat();
    let pool = Pool::new(3);
    for queue in [Box::new(Inline::new()) as _, Box::new(pool.ordered(3)) as _] {
        let mut reader = Reader::new(&raw[..], queue);
        let mut actual = Vec::new();
        let mut piece = [0; 37];
        loop {
            let n = reader.read(&mut piece).unwrap();
            if n == 0 {
                break;
            }
            actual.extend_from_slice(&piece[..n]);
        }
        assert_eq!(actual, expected);
        assert_eq!(reader.compressed_consumed(), raw.len() as u64);
        assert_eq!(reader.read(&mut piece).unwrap(), 0);
    }
}

#[test]
fn a_failed_member_remains_failed_after_consumed_storage_is_recycled() {
    let prefix = stored(b"prefix", false);
    let mut invalid = stored(b"bad", true);
    let crc_at = invalid.len() - 8;
    invalid[crc_at] ^= 1;
    let mut raw = prefix.clone();
    raw.extend_from_slice(&EOF_BLOCK);
    raw.extend_from_slice(&invalid);
    raw.extend_from_slice(&stored(b"unreachable", false));
    let mut reader = Reader::new(&raw[..], Box::new(Inline::new()));
    let mut out = [0; 64];
    assert_eq!(reader.read(&mut out).unwrap(), 6);
    assert_eq!(&out[..6], b"prefix");
    assert!(reader.read(&mut out).is_err());
    assert!(reader.read(&mut out).is_err());
    assert_eq!(
        reader.compressed_consumed(),
        (prefix.len() + EOF_BLOCK.len() + invalid.len()) as u64
    );
}
