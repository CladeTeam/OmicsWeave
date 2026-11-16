//! Sequential delivery, short reads and terminal decode failures.
use crate::support::{
    execution::Ready,
    io::Fragmented,
    wire::{stored, stored_with_extra},
};
use omicsweave_bgzf::{EOF_BLOCK, Reader, read_block};
use omicsweave_core::exec::{Inline, Pool};
use std::io::{self, Read};

#[test]
fn block_headers_and_bodies_survive_short_reads_with_reused_storage() {
    let blocks = [
        stored_with_extra(b"one", false),
        stored_with_extra(b"extended", true),
        EOF_BLOCK.to_vec(),
    ];
    let stream = blocks.concat();
    for capacity in [0, 17, 4096] {
        let mut source = Fragmented {
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
        raw.extend_from_slice(&stored_with_extra(text, i % 2 != 0));
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
    let prefix = stored_with_extra(b"prefix", false);
    let mut invalid = stored_with_extra(b"bad", true);
    let crc_at = invalid.len() - 8;
    invalid[crc_at] ^= 1;
    let mut raw = prefix.clone();
    raw.extend_from_slice(&EOF_BLOCK);
    raw.extend_from_slice(&invalid);
    raw.extend_from_slice(&stored_with_extra(b"unreachable", false));
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

#[test]
fn reader_delivers_valid_queued_prefix_before_later_framing_error() {
    let mut raw = stored(b"valid prefix");
    raw.extend_from_slice(&[31, 139, 8]);
    let mut reader = Reader::new(&raw[..], Ready::boxed(3, None));
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
    let mut reader = Reader::new(&raw[..], Ready::boxed(2, None));
    let mut out = [0; 10];
    assert_eq!(reader.read(&mut out).unwrap(), 1);
    assert!(reader.read(&mut out).is_err());
    assert!(reader.read(&mut out).is_err());
}

#[test]
fn reader_delivers_accepted_prefix_before_submission_failure() {
    let first = stored(b"first");
    let mut raw = first.clone();
    raw.extend_from_slice(&stored(b"rejected"));
    let queue = Ready::boxed(usize::MAX, Some(1));
    let mut reader = Reader::new(&raw[..], queue);
    let mut out = [0; 20];
    assert_eq!(reader.read(&mut out).unwrap(), 5);
    assert_eq!(&out[..5], b"first");
    assert_eq!(reader.compressed_consumed(), first.len() as u64);
    assert!(reader.read(&mut out).is_err());
    assert!(reader.read(&mut out).is_err());
    assert_eq!(reader.compressed_consumed(), first.len() as u64);
}

struct InterruptedOnce<'a>(&'a [u8], bool);
impl Read for InterruptedOnce<'_> {
    fn read(&mut self, out: &mut [u8]) -> io::Result<usize> {
        if self.1 {
            self.1 = false;
            Err(io::ErrorKind::Interrupted.into())
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

#[test]
fn each_incomplete_member_reports_truncation_and_reused_storage_can_read_again() {
    for extended in [false, true] {
        let raw = stored_with_extra(b"payload", extended);
        let mut buffer = Vec::with_capacity(raw.len() * 2);
        for cut in 0..raw.len() {
            let result = read_block(&mut &raw[..cut], &mut buffer);
            if cut == 0 {
                assert!(!result.unwrap());
                assert!(buffer.is_empty());
            } else {
                assert!(
                    matches!(result, Err(omicsweave_core::Error::Truncated { .. })),
                    "extended={extended} cut={cut}: {result:?}"
                );
            }
            assert!(read_block(&mut &raw[..], &mut buffer).unwrap());
            assert_eq!(buffer, raw, "extended={extended} cut={cut}");
        }
    }
}

#[test]
fn empty_reads_do_not_consume_input_or_hide_later_failure() {
    let mut raw = stored(b"prefix");
    raw.extend_from_slice(&[31, 139, 8]);
    let source = Fragmented {
        input: &raw,
        calls: 0,
    };
    let mut reader = Reader::new(source, Ready::boxed(3, None));
    assert_eq!(reader.read(&mut []).unwrap(), 0);
    assert_eq!(reader.compressed_consumed(), 0);
    let mut out = [0; 6];
    reader.read_exact(&mut out).unwrap();
    assert_eq!(&out, b"prefix");
    assert!(reader.read(&mut out).is_err());
    assert_eq!(reader.read(&mut []).unwrap(), 0);
    assert!(reader.read(&mut out).is_err());
}
