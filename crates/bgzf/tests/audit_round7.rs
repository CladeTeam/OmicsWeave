//! Independent stored-DEFLATE fixtures and queue/sink fault injection.

use std::collections::VecDeque;
use std::io::{self, Read, Write};
use std::sync::{Arc, Mutex};

use omicsweave_bgzf::{
    BinningIndex, Blocks, EOF_BLOCK, IndexFormat, Logical, Reader, Writer, inflate_block_owned,
    try_reg2bin,
};
use omicsweave_core::exec::{Inline, Job, Ordered};
use omicsweave_core::{Error, Result};

fn stored(text: &[u8]) -> Vec<u8> {
    assert!(text.len() <= 65_505);
    let mut crc = !0u32;
    for &byte in text {
        crc ^= u32::from(byte);
        for _ in 0..8 {
            crc = (crc >> 1) ^ (0xedb8_8320 & 0u32.wrapping_sub(crc & 1));
        }
    }
    let n = text.len() as u16;
    let mut raw = vec![
        31, 139, 8, 4, 0, 0, 0, 0, 0, 255, 6, 0, 66, 67, 2, 0, 0, 0, 1,
    ];
    raw.extend_from_slice(&n.to_le_bytes());
    raw.extend_from_slice(&(!n).to_le_bytes());
    raw.extend_from_slice(text);
    raw.extend_from_slice(&(!crc).to_le_bytes());
    raw.extend_from_slice(&(text.len() as u32).to_le_bytes());
    let bsize = (raw.len() - 1) as u16;
    raw[16..18].copy_from_slice(&bsize.to_le_bytes());
    raw
}

// Raw DEFLATE of 65,536 zero bytes, fixed from Python's zlib (wbits=-15).
// Its trailer CRC is the independent zlib.crc32 result.
fn full_zero_block() -> Vec<u8> {
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

struct RejectSecond {
    accepted: bool,
    done: VecDeque<Result<Vec<u8>>>,
}

impl Ordered<Vec<u8>> for RejectSecond {
    fn submit(&mut self, job: Job<Vec<u8>>) -> Result<()> {
        if self.accepted {
            return Err(Error::invalid("fixture", "submission rejected"));
        }
        self.accepted = true;
        self.done.push_back(job());
        Ok(())
    }

    fn next(&mut self, _: bool) -> Result<Option<Vec<u8>>> {
        self.done.pop_front().transpose()
    }

    fn pending(&self) -> usize {
        self.done.len()
    }

    fn full(&self) -> bool {
        false
    }
}

#[test]
fn reader_delivers_accepted_prefix_before_submission_failure() {
    let first = stored(b"first");
    let mut raw = first.clone();
    raw.extend_from_slice(&stored(b"rejected"));
    let queue = RejectSecond {
        accepted: false,
        done: VecDeque::new(),
    };
    let mut reader = Reader::new(&raw[..], Box::new(queue));
    let mut out = [0; 20];
    assert_eq!(reader.read(&mut out).unwrap(), 5);
    assert_eq!(&out[..5], b"first");
    assert_eq!(reader.compressed_consumed(), first.len() as u64);
    assert!(reader.read(&mut out).is_err());
    assert!(reader.read(&mut out).is_err());
    assert_eq!(reader.compressed_consumed(), first.len() as u64);
}

#[test]
fn flush_try_large_hint_cuts_current_block_without_overflow() {
    let mut writer = Writer::new(Vec::new(), 6, Box::new(Inline::new())).unwrap();
    writer.write(b"a").unwrap();
    writer.flush_try(usize::MAX).unwrap();
    assert_eq!(writer.tell(), Logical::new(1, 0));
    writer.write(b"b").unwrap();
    let file = writer.finish().unwrap();
    assert_eq!(file.blocks.starts.len(), 2);
    let mut text = Vec::new();
    Reader::new(&file.out[..], Box::new(Inline::new()))
        .read_to_end(&mut text)
        .unwrap();
    assert_eq!(text, b"ab");
}

struct PartialFailure(Arc<Mutex<Vec<u8>>>);

impl Write for PartialFailure {
    fn write(&mut self, input: &[u8]) -> io::Result<usize> {
        let mut bytes = self.0.lock().unwrap();
        if bytes.is_empty() {
            bytes.extend_from_slice(&input[..3]);
            Ok(3)
        } else {
            Err(io::Error::other("injected partial write failure"))
        }
    }

    fn flush(&mut self) -> io::Result<()> {
        Ok(())
    }
}

#[test]
fn partial_sink_failure_prevents_further_output() {
    let bytes = Arc::new(Mutex::new(Vec::new()));
    let mut writer =
        Writer::new(PartialFailure(bytes.clone()), 6, Box::new(Inline::new())).unwrap();
    writer.write(b"a").unwrap();
    assert!(writer.flush_block().is_err());
    assert!(writer.write(b"later").is_err());
    assert!(writer.sync().is_err());
    assert!(writer.finish().is_err());
    assert_eq!(bytes.lock().unwrap().len(), 3);
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
#[ignore = "manual timing workload: run only on the designated CPU pod"]
fn audit_owned_inflate_workload() {
    for len in [0, 1, 4096, 65_536] {
        let text = vec![0; len];
        let raw = if len == 65_536 {
            full_zero_block()
        } else {
            stored(&text)
        };
        assert_eq!(inflate_block_owned(&raw).unwrap(), text);
        for sample in 0..5 {
            let mut rounds = 0u64;
            let begin = std::time::Instant::now();
            while rounds < 100 || begin.elapsed().as_millis() < 100 {
                let decoded = inflate_block_owned(std::hint::black_box(&raw)).unwrap();
                std::hint::black_box(decoded);
                rounds += 1;
            }
            let ns = begin.elapsed().as_nanos();
            println!(
                "name=owned_inflate_{len} sample={sample} iterations={rounds} total_ns={ns} ns_per_iter={:.3}",
                ns as f64 / rounds as f64
            );
        }
    }
}
