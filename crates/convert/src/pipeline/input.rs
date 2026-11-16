//! Reading sources: telling their compression, inflating them, and cutting
//! text into batches of whole lines.
//!
//! Inputs are read with plain sequential `read`s (no mmap, §1.4 of the
//! design): BGZF is inflated block by block on the pool, plain gzip on a
//! thread of its own (it cannot be split), text as is.

use std::fs::File;
use std::io::{self, BufRead, BufReader, Read};
use std::path::Path;
use std::sync::Arc;
use std::sync::atomic::{AtomicU64, Ordering};
use std::sync::mpsc::{Receiver, sync_channel};

use omicsweave_core::exec::Pool;
use omicsweave_core::{Error, Result};

/// Bytes per read from input files.
pub const READ_CHUNK: usize = 1 << 20;

/// How an input file is compressed.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Compression {
    None,
    Gzip,
    Bgzf,
}

/// Detects the compression of the file at `path` from its first bytes.
pub fn detect(path: &Path) -> Result<Compression> {
    let mut head = [0u8; omicsweave_bgzf::HEADER_LEN];
    let mut f = File::open(path).map_err(|e| with_path(e, path))?;
    let mut n = 0;
    while n < head.len() {
        let got = f.read(&mut head[n..])?;
        if got == 0 {
            break;
        }
        n += got;
    }
    Ok(if omicsweave_bgzf::is_bgzf(&head[..n]) {
        Compression::Bgzf
    } else if n >= 2 && head[0] == 0x1f && head[1] == 0x8b {
        Compression::Gzip
    } else {
        Compression::None
    })
}

pub fn with_path(e: io::Error, path: &Path) -> Error {
    Error::Io {
        kind: e.kind(),
        detail: format!("{}: {e}", path.display()),
    }
}

/// Counts the bytes read from a file, for progress.
struct Counted<R> {
    inner: R,
    count: Arc<AtomicU64>,
}

impl<R: Read> Read for Counted<R> {
    fn read(&mut self, buf: &mut [u8]) -> io::Result<usize> {
        let n = self.inner.read(buf)?;
        self.count.fetch_add(n as u64, Ordering::Relaxed);
        Ok(n)
    }
}

/// The text of an input on a background thread: a stream that cannot be
/// split (plain gzip) still overlaps its decoding with the consumer.
struct Threaded {
    rx: Receiver<io::Result<Vec<u8>>>,
    buf: Vec<u8>,
    at: usize,
    done: bool,
}

impl Threaded {
    fn spawn(mut inner: impl Read + Send + 'static) -> Threaded {
        let (tx, rx) = sync_channel(4);
        std::thread::Builder::new()
            .name("omicsweave-gunzip".into())
            .spawn(move || {
                loop {
                    let mut buf = vec![0u8; READ_CHUNK];
                    let mut n = 0;
                    let r = loop {
                        match inner.read(&mut buf[n..]) {
                            Ok(0) => break Ok(()),
                            Ok(k) => {
                                n += k;
                                if n == buf.len() {
                                    break Ok(());
                                }
                            }
                            Err(e) if e.kind() == io::ErrorKind::Interrupted => {}
                            Err(e) => break Err(e),
                        }
                    };
                    buf.truncate(n);
                    let last = n == 0 || r.is_err();
                    let msg = r.map(|()| buf);
                    if tx.send(msg).is_err() || last {
                        return;
                    }
                }
            })
            .expect("spawn decompression thread");
        Threaded {
            rx,
            buf: Vec::new(),
            at: 0,
            done: false,
        }
    }
}

impl Read for Threaded {
    fn read(&mut self, out: &mut [u8]) -> io::Result<usize> {
        while self.at == self.buf.len() {
            if self.done {
                return Ok(0);
            }
            match self.rx.recv() {
                Ok(Ok(b)) if b.is_empty() => self.done = true,
                Ok(Ok(b)) => {
                    self.buf = b;
                    self.at = 0;
                }
                Ok(Err(e)) => {
                    self.done = true;
                    return Err(e);
                }
                Err(_) => self.done = true,
            }
        }
        let n = (self.buf.len() - self.at).min(out.len());
        out[..n].copy_from_slice(&self.buf[self.at..self.at + n]);
        self.at += n;
        Ok(n)
    }
}

/// An opened input: its decompressed text and the count of file bytes read.
pub struct Input {
    pub text: Box<dyn BufRead + Send>,
    pub file_len: u64,
    pub read: Arc<AtomicU64>,
}

/// Opens `path` as text, inflating BGZF on `pool` (at most `limit` blocks
/// ahead) and plain gzip on a thread of its own.
pub fn open_text(path: &Path, pool: &Pool, limit: usize) -> Result<Input> {
    let compression = detect(path)?;
    let file = File::open(path).map_err(|e| with_path(e, path))?;
    let file_len = file.metadata()?.len();
    let read = Arc::new(AtomicU64::new(0));
    let counted = Counted {
        inner: file,
        count: Arc::clone(&read),
    };
    let text: Box<dyn Read + Send> = match compression {
        Compression::None => Box::new(counted),
        Compression::Bgzf => Box::new(omicsweave_bgzf::Reader::new(
            BufReader::with_capacity(READ_CHUNK, counted),
            Box::new(pool.ordered(limit)),
        )),
        Compression::Gzip => Box::new(Threaded::spawn(flate2::read::MultiGzDecoder::new(
            BufReader::with_capacity(READ_CHUNK, counted),
        ))),
    };
    Ok(Input {
        text: Box::new(BufReader::with_capacity(READ_CHUNK, text)),
        file_len,
        read,
    })
}

/// A batch of whole lines (each ended by `\n`, except perhaps the last of
/// the input) and the 1-based number of its first line.
pub struct LineBatch {
    pub text: Vec<u8>,
    pub first_line: u64,
}

/// Cuts a text stream into batches of whole lines of about `target` bytes.
pub struct LineBatches<R> {
    input: R,
    target: usize,
    carry: Vec<u8>,
    next_line: u64,
    done: bool,
}

impl<R: BufRead> LineBatches<R> {
    pub fn new(input: R, target: usize, first_line: u64) -> LineBatches<R> {
        Self::with_prefix(input, target, first_line, Vec::new())
    }

    /// Starts with bytes already read while inspecting the header. They
    /// become the first batch's owned buffer, without a replay copy.
    pub fn with_prefix(
        input: R,
        target: usize,
        first_line: u64,
        prefix: Vec<u8>,
    ) -> LineBatches<R> {
        LineBatches {
            input,
            target: target.max(1),
            carry: prefix,
            next_line: first_line,
            done: false,
        }
    }

    pub fn next_batch(&mut self) -> Result<Option<LineBatch>> {
        if self.done {
            return Ok(None);
        }
        let mut text = std::mem::take(&mut self.carry);
        text.reserve(self.target.saturating_sub(text.len()));
        while text.len() < self.target {
            let buf = self.input.fill_buf()?;
            if buf.is_empty() {
                self.done = true;
                break;
            }
            // Queue sizing budgets one target-sized batch. A large reader
            // buffer must not silently turn every small batch into 1 MiB.
            let n = buf.len().min(self.target - text.len());
            text.extend_from_slice(&buf[..n]);
            self.input.consume(n);
        }
        if !self.done {
            // Keep the incomplete last line for the next batch; a line
            // longer than the target keeps growing this one.
            let mut searched = 0;
            loop {
                match memchr::memrchr(b'\n', &text[searched..]).map(|at| searched + at) {
                    Some(at) => {
                        self.carry.extend_from_slice(&text[at + 1..]);
                        text.truncate(at + 1);
                        break;
                    }
                    None => {
                        searched = text.len();
                        let buf = self.input.fill_buf()?;
                        if buf.is_empty() {
                            self.done = true;
                            break;
                        }
                        // This batch holds a line longer than its target.
                        // Stop at its first newline instead of also taking
                        // unrelated lines from a large input buffer.
                        let n = memchr::memchr(b'\n', buf).map_or(buf.len(), |at| at + 1);
                        text.extend_from_slice(&buf[..n]);
                        self.input.consume(n);
                    }
                }
            }
        }
        if text.is_empty() {
            return Ok(None);
        }
        let first_line = self.next_line;
        self.next_line += memchr::memchr_iter(b'\n', &text).count() as u64;
        if text.last() != Some(&b'\n') {
            self.next_line += 1;
        }
        Ok(Some(LineBatch { text, first_line }))
    }
}

/// Lines of a batch without their `\n` (a `\r` before it is kept), with
/// their 1-based numbers.
pub fn lines(batch: &LineBatch) -> impl Iterator<Item = (u64, &[u8])> {
    let text = &batch.text[..];
    let mut at = 0;
    let mut no = batch.first_line;
    std::iter::from_fn(move || {
        if at >= text.len() {
            return None;
        }
        let end = memchr::memchr(b'\n', &text[at..]).map_or(text.len(), |i| at + i);
        let line = &text[at..end];
        at = end + 1;
        no += 1;
        Some((no - 1, line))
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn short_line_batches_obey_target_with_large_reader_buffers() {
        let input = b"abcdefghij\n".repeat(20_000);
        for target in [11, 37, 1024] {
            let reader = BufReader::with_capacity(READ_CHUNK, input.as_slice());
            let mut batches = LineBatches::new(reader, target, 1);
            let mut actual = Vec::new();
            let mut numbers = Vec::new();
            while let Some(batch) = batches.next_batch().unwrap() {
                assert!(batch.text.len() <= target);
                numbers.extend(lines(&batch).map(|(number, _)| number));
                actual.extend_from_slice(&batch.text);
            }
            assert_eq!(actual, input);
            assert_eq!(numbers, (1..=20_000).collect::<Vec<_>>());
        }
    }

    #[test]
    fn oversized_line_does_not_absorb_following_short_lines() {
        let mut input = b"first\n".to_vec();
        input.extend(std::iter::repeat_n(b'x', 70_003));
        input.extend_from_slice(b"\nlast\nnext\n");
        let mut batches = LineBatches::new(input.as_slice(), 16, 8);
        let first = batches.next_batch().unwrap().unwrap();
        let long = batches.next_batch().unwrap().unwrap();
        let last = batches.next_batch().unwrap().unwrap();
        assert_eq!(
            (first.first_line, first.text.as_slice()),
            (8, &b"first\n"[..])
        );
        assert_eq!(long.first_line, 9);
        assert_eq!(long.text.len(), 70_004);
        assert_eq!(long.text.last(), Some(&b'\n'));
        assert_eq!(
            (last.first_line, last.text.as_slice()),
            (10, &b"last\nnext\n"[..])
        );
        assert!(batches.next_batch().unwrap().is_none());
        assert_eq!([first.text, long.text, last.text].concat(), input);
    }

    #[test]
    fn long_lines_across_small_reads_preserve_bytes_and_numbers() {
        let input = [
            vec![b'x'; 100_003],
            b"\nshort\n".to_vec(),
            vec![b'y'; 70_001],
        ]
        .concat();
        let reader = BufReader::with_capacity(31, input.as_slice());
        let mut batches = LineBatches::new(reader, 37, 5);
        let mut actual = Vec::new();
        let mut numbers = Vec::new();
        while let Some(batch) = batches.next_batch().unwrap() {
            numbers.extend(lines(&batch).map(|(number, _)| number));
            actual.extend(batch.text);
        }
        assert_eq!(actual, input);
        assert_eq!(numbers, [5, 6, 7]);
    }

    #[test]
    fn batches_hold_whole_lines_and_count_them() {
        let text: Vec<u8> = (0..1000)
            .map(|i| format!("line {i} {}\n", "x".repeat(i % 37)))
            .collect::<String>()
            .into_bytes();
        let mut with_tail = text.clone();
        with_tail.extend_from_slice(b"last");
        for input in [text, with_tail] {
            for target in [1, 10, 333, 1 << 20] {
                let mut b = LineBatches::new(&input[..], target, 1);
                let mut all = Vec::new();
                let mut numbers = Vec::new();
                while let Some(batch) = b.next_batch().unwrap() {
                    assert!(
                        batch.text.last() == Some(&b'\n')
                            || all.len() + batch.text.len() == input.len()
                    );
                    for (no, line) in lines(&batch) {
                        numbers.push(no);
                        let _ = line;
                    }
                    all.extend_from_slice(&batch.text);
                }
                assert_eq!(all, input);
                assert_eq!(numbers, (1..=numbers.len() as u64).collect::<Vec<_>>());
            }
        }
    }
}
