//! Independent byte fixtures and fault injection for the BGZF audit.

use std::collections::VecDeque;
use std::io::{self, Read, Write};

use omicsweave_bgzf::{
    BinningIndex, Blocks, EOF_BLOCK, IndexFormat, Inflater, Logical, Reader, Writer, read_block,
    try_csi_levels, try_reg2bin,
};
use omicsweave_core::Result;
use omicsweave_core::exec::{Inline, Job, Ordered};

fn crc32(text: &[u8]) -> u32 {
    let mut crc = !0u32;
    for &b in text {
        crc ^= b as u32;
        for _ in 0..8 {
            crc = (crc >> 1) ^ (0xedb8_8320 & 0u32.wrapping_sub(crc & 1));
        }
    }
    !crc
}

// One final stored DEFLATE block, constructed without the compression API.
fn stored(text: &[u8]) -> Vec<u8> {
    assert!(text.len() <= 65_505);
    let len = text.len() as u16;
    let mut raw = vec![
        31, 139, 8, 4, 0, 0, 0, 0, 0, 255, 6, 0, b'B', b'C', 2, 0, 0, 0, 1,
    ];
    raw.extend_from_slice(&len.to_le_bytes());
    raw.extend_from_slice(&(!len).to_le_bytes());
    raw.extend_from_slice(text);
    raw.extend_from_slice(&crc32(text).to_le_bytes());
    raw.extend_from_slice(&(text.len() as u32).to_le_bytes());
    let bsize = (raw.len() - 1) as u16;
    raw[16..18].copy_from_slice(&bsize.to_le_bytes());
    raw
}

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

struct Batch {
    done: VecDeque<Result<Vec<u8>>>,
    limit: usize,
}

impl Batch {
    fn boxed(limit: usize) -> Box<dyn Ordered<Vec<u8>>> {
        Box::new(Self {
            done: VecDeque::new(),
            limit,
        })
    }
}

impl Ordered<Vec<u8>> for Batch {
    fn submit(&mut self, job: Job<Vec<u8>>) -> Result<()> {
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
        self.pending() >= self.limit
    }
}

#[test]
fn reader_delivers_valid_queued_prefix_before_later_framing_error() {
    let mut raw = stored(b"valid prefix");
    raw.extend_from_slice(&[31, 139, 8]);
    let mut reader = Reader::new(&raw[..], Batch::boxed(3));
    let mut out = [0; 12];
    assert_eq!(reader.read(&mut out).unwrap(), 12);
    assert_eq!(&out, b"valid prefix");
    assert!(reader.read(&mut out).is_err());
}

#[test]
fn reader_remains_failed_after_corrupt_block() {
    let mut raw = stored(b"a");
    let mut bad = stored(b"broken");
    let at = bad.len() - 8;
    bad[at] ^= 1;
    raw.extend_from_slice(&bad);
    raw.extend_from_slice(&stored(b"later"));
    let mut reader = Reader::new(&raw[..], Batch::boxed(2));
    let mut out = [0; 10];
    assert_eq!(reader.read(&mut out).unwrap(), 1);
    assert!(reader.read(&mut out).is_err());
    assert!(reader.read(&mut out).is_err());
}

struct FlushOnce(bool);

impl Write for FlushOnce {
    fn write(&mut self, text: &[u8]) -> io::Result<usize> {
        Ok(text.len())
    }

    fn flush(&mut self) -> io::Result<()> {
        if self.0 {
            self.0 = false;
            Err(io::Error::other("injected flush failure"))
        } else {
            Ok(())
        }
    }
}

#[test]
fn sync_flush_error_is_terminal() {
    let mut writer = Writer::new(FlushOnce(true), 6, Box::new(Inline::new())).unwrap();
    writer.write(b"a").unwrap();
    assert!(writer.sync().is_err());
    assert!(writer.write(b"later").is_err());
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

struct InterruptedOnce<'a>(&'a [u8], bool);

impl Read for InterruptedOnce<'_> {
    fn read(&mut self, out: &mut [u8]) -> io::Result<usize> {
        if self.1 {
            self.1 = false;
            Err(io::Error::from(io::ErrorKind::Interrupted))
        } else {
            self.0.read(out)
        }
    }
}

#[test]
fn read_block_retries_interrupts_and_reports_actual_truncated_bytes() {
    let raw = stored(b"abc");
    let mut block = Vec::new();
    assert!(read_block(&mut InterruptedOnce(&raw, true), &mut block).unwrap());
    assert_eq!(block, raw);
    let short = &raw[..raw.len() - 2];
    assert!(matches!(
        read_block(&mut &short[..], &mut block),
        Err(omicsweave_core::Error::Truncated { need, have, .. })
        if need == raw.len() as u64 && have == short.len() as u64
    ));
}
