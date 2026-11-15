//! [`Reader`]: the text of a BGZF stream, inflated block by block ahead of
//! the consumer.

use std::collections::VecDeque;
use std::io::{self, Read};

use omicsweave_core::exec::Ordered;
use omicsweave_core::{Error, Result};

use crate::{BlockHeader, HEADER_LEN, inflate_block_owned};

/// Whether `head` (the first bytes of a stream, at least [`HEADER_LEN`]
/// when the stream is that long) starts with a BGZF block rather than
/// plain gzip or text.
pub fn is_bgzf(head: &[u8]) -> bool {
    head.len() >= HEADER_LEN
        && head[0] == 0x1f
        && head[1] == 0x8b
        && head[2] == 8
        && head[3] & 4 != 0
        && BlockHeader::parse(head).is_ok()
}

/// Reads into `buf[from..]` until it is full or the input ends, retrying
/// interrupted reads; returns the bytes now in `buf`.
fn fill(input: &mut impl Read, buf: &mut [u8], from: usize) -> Result<usize> {
    let mut got = from;
    while got < buf.len() {
        match input.read(&mut buf[got..]) {
            Ok(0) => break,
            Ok(n) => got += n,
            Err(e) if e.kind() == io::ErrorKind::Interrupted => {}
            Err(e) => return Err(e.into()),
        }
    }
    Ok(got)
}

/// Reads the next whole block from `input` into `buf` (cleared first);
/// `false` at a clean end of the stream.
pub fn read_block(input: &mut impl Read, buf: &mut Vec<u8>) -> Result<bool> {
    read_block_into(input, buf, false)
}

fn read_block_into(input: &mut impl Read, buf: &mut Vec<u8>, recycled_text: bool) -> Result<bool> {
    buf.clear();
    let mut head = [0; HEADER_LEN];
    let got = fill(input, &mut head, 0)?;
    if got == 0 {
        buf.clear();
        return Ok(false);
    }
    let truncated = |what, need: usize, have: usize| Error::Truncated {
        what,
        need: need as u64,
        have: have as u64,
    };
    if got < HEADER_LEN {
        buf.extend_from_slice(&head);
        return Err(truncated("bgzf block header", HEADER_LEN, got));
    }
    let header = match BlockHeader::parse(&head) {
        Err(Error::Truncated { need, .. }) => {
            let need = need as usize;
            buf.reserve(need);
            buf.extend_from_slice(&head);
            buf.resize(need, 0);
            let got = fill(input, buf, HEADER_LEN)?;
            if got < need {
                return Err(truncated("bgzf block header", need, got));
            }
            BlockHeader::parse(buf)?
        }
        Ok(header) => {
            // The common header describes the whole allocation before any
            // heap storage is needed, avoiding an initial header allocation.
            if recycled_text {
                // Incompressible framing may be slightly larger than the
                // previous text. It must not double a nearly-full block's
                // capacity just to accommodate a few trailer bytes.
                buf.reserve_exact(header.block_len);
            } else {
                buf.reserve(header.block_len);
            }
            buf.extend_from_slice(&head);
            header
        }
        Err(error) => {
            buf.extend_from_slice(&head);
            return Err(error);
        }
    };
    let have = buf.len();
    if recycled_text {
        buf.reserve_exact(header.block_len - have);
    }
    buf.resize(header.block_len, 0);
    let got = fill(input, buf, have)?;
    if got < header.block_len {
        return Err(truncated("bgzf block", header.block_len, got));
    }
    Ok(true)
}

/// The decompressed text of a BGZF stream. Blocks are read from the input
/// in order and inflated on the queue given, up to its limit ahead of what
/// has been consumed; text comes out in stream order.
pub struct Reader<R: Read> {
    input: R,
    exec: Box<dyn Ordered<Vec<u8>>>,
    text: Vec<u8>,
    at: usize,
    input_done: bool,
    /// An input or submission error met while reading ahead; reported once
    /// the accepted blocks before it have been handed out.
    pending: Option<Error>,
    /// Once set, every read fails with it.
    failed: Option<Error>,
    /// Compressed bytes read so far.
    consumed: u64,
    /// Compressed sizes of submitted blocks, for progress.
    sizes: VecDeque<u64>,
}

impl<R: Read> Reader<R> {
    pub fn new(input: R, exec: Box<dyn Ordered<Vec<u8>>>) -> Reader<R> {
        Reader {
            input,
            exec,
            text: Vec::new(),
            at: 0,
            input_done: false,
            pending: None,
            failed: None,
            consumed: 0,
            sizes: VecDeque::new(),
        }
    }

    /// Compressed bytes whose text has been handed out.
    pub fn compressed_consumed(&self) -> u64 {
        self.consumed
    }

    fn fill(&mut self) -> Result<bool> {
        if let Some(e) = &self.failed {
            return Err(e.clone());
        }
        let r = self.fill_next();
        if let Err(e) = &r {
            self.failed = Some(e.clone());
        }
        r
    }

    fn fill_next(&mut self) -> Result<bool> {
        loop {
            while !self.input_done && self.pending.is_none() && !self.exec.full() {
                // All previous text has been copied to the consumer. Its
                // allocation can hold a later compressed block; no decoded
                // bytes survive the read_block clear or enter a data cache.
                let mut raw = std::mem::take(&mut self.text);
                self.at = 0;
                match read_block_into(&mut self.input, &mut raw, true) {
                    Ok(true) => {}
                    Ok(false) => {
                        self.input_done = true;
                        break;
                    }
                    Err(e) => {
                        self.pending = Some(e);
                        break;
                    }
                }
                let size = raw.len() as u64;
                if let Err(e) = self
                    .exec
                    .submit(Box::new(move || inflate_block_owned(&raw)))
                {
                    self.pending = Some(e);
                    break;
                }
                self.sizes.push_back(size);
            }
            let next = self.exec.next(true);
            if !matches!(next, Ok(None)) {
                self.consumed += self.sizes.pop_front().unwrap_or(0);
            }
            match next? {
                Some(text) if text.is_empty() => continue,
                Some(text) => {
                    self.text = text;
                    self.at = 0;
                    return Ok(true);
                }
                None => match self.pending.take() {
                    Some(e) => return Err(e),
                    None => return Ok(false),
                },
            }
        }
    }
}

impl<R: Read> Read for Reader<R> {
    fn read(&mut self, out: &mut [u8]) -> io::Result<usize> {
        if out.is_empty() {
            return Ok(0);
        }
        if self.at == self.text.len() && !self.fill().map_err(io::Error::other)? {
            return Ok(0);
        }
        let n = (self.text.len() - self.at).min(out.len());
        out[..n].copy_from_slice(&self.text[self.at..self.at + n]);
        self.at += n;
        Ok(n)
    }
}
