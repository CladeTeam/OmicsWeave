use std::collections::VecDeque;
use std::io::{self, Read, Write};
use std::sync::{Arc, Mutex};

use omicsweave_bgzf::{BLOCK_TEXT_LEN, EOF_BLOCK, Logical, Reader, Writer};
use omicsweave_core::exec::{Inline, Job, Ordered, Pool};
use omicsweave_core::{Error, Result};

struct Deferred(VecDeque<Job<Vec<u8>>>);

impl Ordered<Vec<u8>> for Deferred {
    fn submit(&mut self, job: Job<Vec<u8>>) -> Result<()> {
        self.0.push_back(job);
        Ok(())
    }
    fn next(&mut self, wait: bool) -> Result<Option<Vec<u8>>> {
        if wait {
            self.0.pop_front().map(|job| job()).transpose()
        } else {
            Ok(None)
        }
    }
    fn pending(&self) -> usize {
        self.0.len()
    }
    fn full(&self) -> bool {
        self.0.len() >= 3
    }
}

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

struct PartialFailure(Arc<Mutex<Vec<u8>>>);
impl Write for PartialFailure {
    fn write(&mut self, input: &[u8]) -> io::Result<usize> {
        let mut bytes = self.0.lock().unwrap();
        if bytes.is_empty() {
            bytes.extend_from_slice(&input[..3]);
            Ok(3)
        } else {
            Err(io::Error::other("partial sink failure"))
        }
    }
    fn flush(&mut self) -> io::Result<()> {
        Ok(())
    }
}

struct FlushFailure;
impl Write for FlushFailure {
    fn write(&mut self, input: &[u8]) -> io::Result<usize> {
        Ok(input.len())
    }
    fn flush(&mut self) -> io::Result<()> {
        Err(io::Error::other("flush failure"))
    }
}

struct EofFailure;
impl Write for EofFailure {
    fn write(&mut self, input: &[u8]) -> io::Result<usize> {
        if input == EOF_BLOCK {
            Err(io::Error::other("EOF write failure"))
        } else {
            Ok(input.len())
        }
    }
    fn flush(&mut self) -> io::Result<()> {
        Ok(())
    }
}

struct Reject;
impl Ordered<Vec<u8>> for Reject {
    fn submit(&mut self, _: Job<Vec<u8>>) -> Result<()> {
        Err(Error::invalid("queue", "submission rejected"))
    }
    fn next(&mut self, _: bool) -> Result<Option<Vec<u8>>> {
        Ok(None)
    }
    fn pending(&self) -> usize {
        0
    }
    fn full(&self) -> bool {
        false
    }
}

#[test]
fn streaming_partial_writes_and_sync_failures_are_terminal() {
    let bytes = Arc::new(Mutex::new(Vec::new()));
    let mut writer =
        Writer::new_streaming(PartialFailure(bytes.clone()), 6, Box::new(Inline::new())).unwrap();
    writer.write(b"a").unwrap();
    assert!(writer.flush_block().is_err());
    assert!(writer.write(b"later").is_err());
    assert!(writer.sync().is_err());
    assert!(writer.finish_streaming().is_err());
    assert_eq!(bytes.lock().unwrap().len(), 3);
    let mut writer = Writer::new_streaming(FlushFailure, 6, Box::new(Inline::new())).unwrap();
    writer.write(b"a").unwrap();
    assert!(writer.sync().is_err());
    assert!(writer.write(b"later").is_err());
    assert!(writer.flush_block().is_err());
    assert!(writer.finish_streaming().is_err());
    let mut writer = Writer::new_streaming(Vec::new(), 6, Box::new(Reject)).unwrap();
    writer.write(b"a").unwrap();
    assert!(writer.flush_block().is_err());
    assert!(writer.write(b"later").is_err());
    assert!(writer.sync().is_err());
    assert!(writer.finish_streaming().is_err());
}

#[test]
fn streaming_finish_propagates_eof_and_final_flush_errors() {
    let writer = Writer::new_streaming(EofFailure, 6, Box::new(Inline::new())).unwrap();
    assert!(writer.finish_streaming().is_err());
    let writer = Writer::new_streaming(FlushFailure, 6, Box::new(Inline::new())).unwrap();
    assert!(writer.finish_streaming().is_err());
}
