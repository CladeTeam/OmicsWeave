//! Catalog and streaming writers share delivery and failure contracts.
use crate::support::compression::{decode_gzip_block, text};
use crate::support::{
    execution::{Deferred, Ready},
    io::FailAfter,
};
use omicsweave_bgzf::{BLOCK_TEXT_LEN, EOF_BLOCK, Logical, Reader, Writer};
use omicsweave_core::exec::{Inline, Ordered, Pool};
use std::collections::VecDeque;
use std::io::{self, Read, Write};
use std::sync::{Arc, Mutex};

fn output(
    records: &[Vec<u8>],
    level: i32,
    streaming: bool,
    exec: Box<dyn Ordered<Vec<u8>>>,
) -> (Vec<u8>, Vec<Logical>) {
    let mut writer = if streaming {
        Writer::new_streaming(Vec::new(), level, exec)
    } else {
        Writer::new(Vec::new(), level, exec)
    }
    .unwrap();
    let mut positions = vec![writer.tell()];
    let mut total = 0;
    for (i, record) in records.iter().enumerate() {
        if i % 5 == 0 {
            writer.flush_try(usize::MAX).unwrap();
        }
        writer.write_record(record).unwrap();
        if i % 7 == 1 {
            writer.flush_block().unwrap();
        }
        if i % 11 == 0 {
            writer.sync().unwrap();
        }
        total += record.len();
        assert_eq!(writer.text_len(), total as u64);
        positions.push(writer.tell());
    }
    let raw = if streaming {
        writer.finish_streaming().unwrap()
    } else {
        let finished = writer.finish().unwrap();
        assert_eq!(finished.blocks.text_len, total as u64);
        finished.out
    };
    assert!(raw.ends_with(&EOF_BLOCK));
    let mut decoded = Vec::new();
    Reader::new(&raw[..], Box::new(Inline::new()))
        .read_to_end(&mut decoded)
        .unwrap();
    assert_eq!(decoded, records.concat());
    (raw, positions)
}

#[test]
fn streaming_matches_catalog_bytes_and_positions_for_ready_delayed_and_pool_queues() {
    let records: Vec<_> = (0..80usize)
        .map(|i| {
            let len = if i % 13 == 0 {
                BLOCK_TEXT_LEN * 2 + 17
            } else {
                i * 127 % 9000
            };
            (0..len).map(|j| ((i * 23 + j * 97) % 251) as u8).collect()
        })
        .collect();
    let pool = Pool::new(3);
    for level in [0, 6, 12] {
        let expected = output(&records, level, false, Box::new(Inline::new()));
        assert_eq!(
            output(&records, level, true, Box::new(Inline::new())),
            expected
        );
        assert_eq!(
            output(&records, level, true, Box::new(Deferred(VecDeque::new()))),
            expected
        );
        assert_eq!(
            output(&records, level, true, Box::new(pool.ordered(3))),
            expected
        );
    }
}

#[test]
fn finish_requires_the_explicit_streaming_contract_and_supports_empty_output() {
    let writer = Writer::new_streaming(Vec::new(), 6, Box::new(Inline::new())).unwrap();
    let error = writer.finish().err().unwrap().to_string();
    assert!(error.contains("finish_streaming"));
    let writer = Writer::new_streaming(Vec::new(), 6, Box::new(Inline::new())).unwrap();
    assert_eq!(writer.finish_streaming().unwrap(), EOF_BLOCK);
    let mut writer = Writer::new(Vec::new(), 6, Box::new(Inline::new())).unwrap();
    writer.write(b"full mode may discard the catalog").unwrap();
    let raw = writer.finish_streaming().unwrap();
    let mut text = Vec::new();
    Reader::new(&raw[..], Box::new(Inline::new()))
        .read_to_end(&mut text)
        .unwrap();
    assert_eq!(text, b"full mode may discard the catalog");
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

struct FlushOnce(bool);
impl Write for FlushOnce {
    fn write(&mut self, input: &[u8]) -> io::Result<usize> {
        Ok(input.len())
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
fn sync_flush_error_is_terminal_in_both_modes() {
    for streaming in [false, true] {
        let queue = Box::new(Inline::new());
        let mut writer = if streaming {
            Writer::new_streaming(FlushOnce(true), 6, queue)
        } else {
            Writer::new(FlushOnce(true), 6, queue)
        }
        .unwrap();
        writer.write(b"a").unwrap();
        assert!(writer.sync().is_err(), "streaming={streaming}");
        assert!(writer.write(b"later").is_err());
        assert!(writer.flush_block().is_err());
        assert!(writer.finish_streaming().is_err());
    }
}

#[test]
fn partial_sink_failure_prevents_further_output_in_both_modes() {
    for streaming in [false, true] {
        let bytes = Arc::new(Mutex::new(Vec::new()));
        let sink = FailAfter {
            bytes: bytes.clone(),
            limit: 3,
        };
        let queue = Box::new(Inline::new());
        let mut writer = if streaming {
            Writer::new_streaming(sink, 6, queue)
        } else {
            Writer::new(sink, 6, queue)
        }
        .unwrap();
        writer.write(b"a").unwrap();
        assert!(writer.flush_block().is_err(), "streaming={streaming}");
        assert!(writer.write(b"later").is_err());
        assert!(writer.sync().is_err());
        if streaming {
            assert!(writer.finish_streaming().is_err());
        } else {
            assert!(writer.finish().is_err());
        }
        assert_eq!(bytes.lock().unwrap().len(), 3, "streaming={streaming}");
    }
}

#[test]
fn rejected_compression_submission_is_terminal_in_both_modes() {
    for streaming in [false, true] {
        let queue = Ready::boxed(usize::MAX, Some(0));
        let mut writer = if streaming {
            Writer::new_streaming(Vec::new(), 6, queue)
        } else {
            Writer::new(Vec::new(), 6, queue)
        }
        .unwrap();
        writer.write(b"a").unwrap();
        assert!(writer.flush_block().is_err(), "streaming={streaming}");
        assert!(writer.write(b"later").is_err());
        assert!(writer.sync().is_err());
        assert!(writer.finish_streaming().is_err());
    }
}

struct FinishFailure {
    eof: bool,
}
impl Write for FinishFailure {
    fn write(&mut self, input: &[u8]) -> io::Result<usize> {
        if self.eof && input == EOF_BLOCK {
            Err(io::Error::other("EOF write failure"))
        } else {
            Ok(input.len())
        }
    }
    fn flush(&mut self) -> io::Result<()> {
        if self.eof {
            Ok(())
        } else {
            Err(io::Error::other("final flush failure"))
        }
    }
}

#[test]
fn finish_propagates_eof_and_final_flush_errors_in_both_modes() {
    for streaming in [false, true] {
        for eof in [false, true] {
            let sink = FinishFailure { eof };
            let queue = Box::new(Inline::new());
            let writer = if streaming {
                Writer::new_streaming(sink, 6, queue)
            } else {
                Writer::new(sink, 6, queue)
            }
            .unwrap();
            let failed = if streaming {
                writer.finish_streaming().is_err()
            } else {
                writer.finish().is_err()
            };
            assert!(failed, "streaming={streaming} eof={eof}");
        }
    }
}

fn write_records(input: &[u8], blocks: usize) -> Vec<u8> {
    let mut writer = Writer::new(Vec::new(), 6, Box::new(Inline::new())).unwrap();
    for _ in 0..blocks {
        writer.write(input).unwrap();
        writer.flush_block().unwrap();
    }
    let done = writer.finish().unwrap();
    assert_eq!(done.blocks.starts.len(), blocks);
    assert_eq!(done.blocks.text_len, (input.len() * blocks) as u64);
    done.out
}

#[test]
fn short_flushed_records_and_full_records_are_independent_gzip_members() {
    for len in [1, 64, 4096, BLOCK_TEXT_LEN] {
        let input = text(len, false);
        let raw = write_records(&input, 3);
        let mut at = 0;
        for _ in 0..3 {
            let len = u16::from_le_bytes([raw[at + 16], raw[at + 17]]) as usize + 1;
            decode_gzip_block(&raw[at..at + len], &input);
            at += len;
        }
        assert_eq!(&raw[at..], &EOF_BLOCK);
    }
}
