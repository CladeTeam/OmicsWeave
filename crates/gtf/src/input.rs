//! Annotation text from a file or bytes: compression told by content
//! (plain, gzip, BGZF), and the text cut into batches of whole lines.
//!
//! Files are read with plain sequential reads. BGZF is inflated block by
//! block on the worker pool; plain gzip cannot be split and is inflated on a
//! thread of its own so that it overlaps parsing.

use std::fs::File;
use std::io::{self, BufRead, BufReader, Read};
use std::path::Path;
use std::sync::mpsc::{Receiver, sync_channel};

use omicsweave_core::exec::Exec;
use omicsweave_core::{Error, Result};

const READ_CHUNK: usize = 1 << 20;

/// How a source is compressed.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Compression {
    None,
    Gzip,
    Bgzf,
}

/// The compression of a stream from its first bytes.
pub fn detect_compression(head: &[u8]) -> Compression {
    if omicsweave_bgzf::is_bgzf(head) {
        Compression::Bgzf
    } else if head.len() >= 2 && head[0] == 0x1f && head[1] == 0x8b {
        Compression::Gzip
    } else {
        Compression::None
    }
}

fn with_path(e: io::Error, path: &Path) -> Error {
    Error::Io {
        kind: e.kind(),
        detail: format!("{}: {e}", path.display()),
    }
}

/// Opens `path` as text.
pub(crate) fn open_file(path: &Path, exec: Exec<'_>) -> Result<Box<dyn BufRead>> {
    let mut head = [0u8; omicsweave_bgzf::HEADER_LEN];
    let mut f = File::open(path).map_err(|e| with_path(e, path))?;
    let mut n = 0;
    while n < head.len() {
        let got = f.read(&mut head[n..]).map_err(|e| with_path(e, path))?;
        if got == 0 {
            break;
        }
        n += got;
    }
    let file = File::open(path).map_err(|e| with_path(e, path))?;
    Ok(match detect_compression(&head[..n]) {
        Compression::None => Box::new(BufReader::with_capacity(READ_CHUNK, file)),
        Compression::Bgzf => Box::new(BufReader::with_capacity(
            READ_CHUNK,
            omicsweave_bgzf::Reader::new(BufReader::with_capacity(READ_CHUNK, file), exec.queue()),
        )),
        Compression::Gzip => Box::new(BufReader::with_capacity(
            READ_CHUNK,
            Threaded::spawn(flate2::read::MultiGzDecoder::new(BufReader::with_capacity(
                READ_CHUNK, file,
            )))?,
        )),
    })
}

/// Text of in-memory bytes, inflated when they are gzip or BGZF.
pub(crate) fn open_bytes(bytes: &[u8]) -> Box<dyn BufRead + '_> {
    match detect_compression(&bytes[..bytes.len().min(omicsweave_bgzf::HEADER_LEN)]) {
        Compression::None => Box::new(bytes),
        _ => Box::new(BufReader::with_capacity(
            READ_CHUNK,
            flate2::bufread::MultiGzDecoder::new(bytes),
        )),
    }
}

/// A stream decoded on a background thread.
struct Threaded {
    rx: Receiver<io::Result<Vec<u8>>>,
    buf: Vec<u8>,
    at: usize,
    done: bool,
}

impl Threaded {
    fn spawn(mut inner: impl Read + Send + 'static) -> Result<Threaded> {
        let (tx, rx) = sync_channel(4);
        std::thread::Builder::new()
            .name("omicsweave-gtf-gunzip".into())
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
                    if tx.send(r.map(|()| buf)).is_err() || last {
                        return;
                    }
                }
            })?;
        Ok(Threaded {
            rx,
            buf: Vec::new(),
            at: 0,
            done: false,
        })
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

/// Whole lines (each ended by `\n`, except perhaps the last of the input)
/// and the 1-based number of the first.
pub(crate) struct LineBatch {
    pub text: Vec<u8>,
    pub first_line: u64,
}

/// Cuts a text stream into batches of whole lines of about `target` bytes;
/// a line longer than the target makes a batch of its own.
pub(crate) struct LineBatches<R> {
    input: R,
    target: usize,
    carry: Vec<u8>,
    next_line: u64,
    done: bool,
}

impl<R: BufRead> LineBatches<R> {
    pub fn new(input: R, target: usize) -> LineBatches<R> {
        LineBatches {
            input,
            target: target.max(1),
            carry: Vec::new(),
            next_line: 1,
            done: false,
        }
    }

    pub fn next_batch(&mut self) -> Result<Option<LineBatch>> {
        if self.done && self.carry.is_empty() {
            return Ok(None);
        }
        let mut text = std::mem::take(&mut self.carry);
        text.reserve(self.target.saturating_sub(text.len()));
        while !self.done && text.len() < self.target {
            let buf = self.input.fill_buf()?;
            if buf.is_empty() {
                self.done = true;
                break;
            }
            let n = buf.len().min(self.target - text.len());
            text.extend_from_slice(&buf[..n]);
            self.input.consume(n);
        }
        if !self.done {
            let mut searched = 0;
            loop {
                if let Some(at) = memchr::memrchr(b'\n', &text[searched..]).map(|i| searched + i) {
                    self.carry.extend_from_slice(&text[at + 1..]);
                    text.truncate(at + 1);
                    break;
                }
                searched = text.len();
                let buf = self.input.fill_buf()?;
                if buf.is_empty() {
                    self.done = true;
                    break;
                }
                let n = memchr::memchr(b'\n', buf).map_or(buf.len(), |i| i + 1);
                text.extend_from_slice(&buf[..n]);
                self.input.consume(n);
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

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn batches_keep_every_byte_and_line_number() {
        let text: Vec<u8> = (0..2000)
            .map(|i| format!("line {i} {}\n", "x".repeat(i % 53)))
            .collect::<String>()
            .into_bytes();
        let mut tail = text.clone();
        tail.extend_from_slice(b"last");
        for input in [text, tail] {
            for target in [1, 17, 1000, 1 << 20] {
                let mut b = LineBatches::new(BufReader::with_capacity(64, &input[..]), target);
                let mut all = Vec::new();
                let mut line = 1;
                while let Some(batch) = b.next_batch().unwrap() {
                    assert_eq!(batch.first_line, line);
                    line += batch.text.split(|&c| c == b'\n').count() as u64
                        - u64::from(batch.text.last() == Some(&b'\n'));
                    all.extend_from_slice(&batch.text);
                }
                assert_eq!(all, input);
            }
        }
    }
}
