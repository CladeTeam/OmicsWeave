//! A whole BGZF stream read front to back, with the position `bgzf_tell`
//! reports after every byte consumed: what sequential decoders (BAM and
//! BCF records) and index builders read through.

use std::collections::VecDeque;
use std::io::Read;
use std::ops::Range;

use omicsweave_core::exec::{Inline, Job, Ordered, Pool, PoolQueue};
use omicsweave_core::{Error, Result};

use crate::{Blocks, Logical, inflate_block_owned, read_block};

/// Consumed text kept before it is dropped from the front.
const KEEP: usize = 4 << 20;

#[derive(Debug, Clone, Copy)]
struct Inflated {
    /// The block's number in the file.
    seq: u64,
    text: usize,
    len: usize,
}

/// The text of a whole BGZF stream, front to back: blocks inflate ahead of
/// the reader on an [`Ordered`] queue, and every position consumed comes
/// with the [`Logical`] offset htslib's `bgzf_tell` reports there. With its
/// block directory ([`Scanner::new`]) those positions resolve to virtual
/// offsets, as an index needs them.
pub struct Scanner<R: Read> {
    input: R,
    exec: Box<dyn Ordered<Vec<u8>>>,
    input_done: bool,
    /// An input error met while reading ahead, reported once the blocks
    /// before it are used.
    pending: Option<Error>,
    failed: Option<Error>,
    /// Compressed address of the next block read from the input.
    address: u64,
    /// Blocks taken from the queue so far.
    taken: u64,
    text: Vec<u8>,
    loaded: VecDeque<Inflated>,
    /// The loaded block holding the last consumed byte; advances monotonically.
    cur: usize,
    #[cfg(test)]
    cursor_steps: usize,
    #[cfg(test)]
    peak_loaded: usize,
    pos: usize,
    /// Where every block starts, compressed and in the text.
    blocks: Blocks,
    keep_blocks: bool,
    text_total: u64,
    /// `bgzf_tell` after the bytes consumed, with block numbers.
    here: Logical,
}

impl<R: Read> Scanner<R> {
    /// Reads `input` from its start; blocks inflate on `exec`.
    pub fn new(input: R, exec: Box<dyn Ordered<Vec<u8>>>) -> Scanner<R> {
        Self::with_catalog(input, exec, true)
    }

    /// Reads records without retaining the file-wide block directory.
    ///
    /// Use this for one-pass consumers that need records and input progress,
    /// but do not build an index. The logical positions returned by `tell`
    /// are unchanged; this scanner has no directory to resolve them through.
    pub fn new_streaming(input: R, exec: Box<dyn Ordered<Vec<u8>>>) -> Scanner<R> {
        Self::with_catalog(input, exec, false)
    }

    pub(crate) fn with_catalog(
        input: R,
        exec: Box<dyn Ordered<Vec<u8>>>,
        keep_blocks: bool,
    ) -> Scanner<R> {
        Scanner {
            input,
            exec,
            input_done: false,
            pending: None,
            failed: None,
            address: 0,
            taken: 0,
            text: Vec::new(),
            loaded: VecDeque::new(),
            cur: 0,
            #[cfg(test)]
            cursor_steps: 0,
            #[cfg(test)]
            peak_loaded: 0,
            pos: 0,
            blocks: Blocks::default(),
            keep_blocks,
            text_total: 0,
            here: Logical::new(0, 0),
        }
    }

    /// Takes the next inflated block; `false` at the end of the stream.
    fn load_next(&mut self) -> Result<bool> {
        while !self.input_done && self.pending.is_none() && !self.exec.full() {
            let mut raw = Vec::new();
            match read_block(&mut self.input, &mut raw) {
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
            let len = raw.len() as u64;
            self.exec
                .submit(Box::new(move || inflate_block_owned(&raw)))?;
            if self.keep_blocks {
                self.blocks.starts.push((self.address, 0));
            }
            self.address += len;
        }
        let Some(block) = self.exec.next(true)? else {
            self.blocks.end = self.address;
            self.blocks.text_len = self.text_total;
            return match self.pending.take() {
                Some(e) => Err(e),
                None => Ok(false),
            };
        };
        let seq = self.taken;
        self.taken += 1;
        if self.keep_blocks {
            self.blocks.starts[seq as usize].1 = self.text_total;
        }
        self.text_total += block.len() as u64;
        // Empty members have no byte that advance can land on. Their file
        // positions still belong to the optional directory, but retaining
        // them in the live-text queue would grow with an empty-member run.
        if block.is_empty() {
            return Ok(true);
        }
        if self.pos >= KEEP && !self.loaded.is_empty() {
            // Drop consumed text, keeping the block holding the position.
            let keep_from = self
                .loaded
                .iter()
                .rposition(|b| b.text <= self.pos)
                .unwrap_or(0);
            let drop = self.loaded[keep_from].text;
            self.text.drain(..drop);
            self.loaded.drain(..keep_from);
            self.cur = self.cur.saturating_sub(keep_from);
            for b in &mut self.loaded {
                b.text -= drop;
            }
            self.pos -= drop;
        }
        let at = self.text.len();
        self.text.extend_from_slice(&block);
        self.loaded.push_back(Inflated {
            seq,
            text: at,
            len: block.len(),
        });
        #[cfg(test)]
        {
            self.peak_loaded = self.peak_loaded.max(self.loaded.len());
        }
        Ok(true)
    }

    /// Reads ahead until `n` unread bytes are available or the stream
    /// ends; returns how many are. An error is reported again by every
    /// later call.
    pub fn fill(&mut self, n: usize) -> Result<usize> {
        if let Some(error) = &self.failed {
            return Err(error.clone());
        }
        let result: Result<usize> = (|| {
            while self.text.len() - self.pos < n {
                if !self.load_next()? {
                    break;
                }
            }
            Ok(self.text.len() - self.pos)
        })();
        if let Err(error) = &result {
            self.failed = Some(error.clone());
        }
        result
    }

    fn advance(&mut self, n: usize) {
        if n == 0 {
            return;
        }
        self.pos += n;
        let last = self.pos - 1;
        while self.loaded[self.cur].text + self.loaded[self.cur].len <= last {
            self.cur += 1;
            #[cfg(test)]
            {
                self.cursor_steps += 1;
            }
        }
        let b = self.loaded[self.cur];
        self.here = if self.pos == b.text + b.len {
            Logical::new(b.seq + 1, 0)
        } else {
            Logical::new(b.seq, (self.pos - b.text) as u16)
        };
        if self.cur != 0 {
            // Slice ranges borrow text, not this directory. Retaining only
            // the current and unread frames bounds metadata by the active
            // fill even when each member contributes a single byte.
            self.loaded.drain(..self.cur);
            self.cur = 0;
        }
    }

    /// `bgzf_tell` after what has been read, with the block's number in
    /// place of its address. A retained [`Scanner::blocks`] directory
    /// resolves it; streaming scanners report the same logical position.
    pub fn tell(&self) -> Logical {
        self.here
    }

    /// Compressed bytes taken from the input so far.
    pub fn compressed_read(&self) -> u64 {
        self.address
    }

    /// Where the blocks read so far start; complete once the stream has
    /// ended.
    ///
    /// # Panics
    ///
    /// Panics for a scanner created with [`Scanner::new_streaming`]. Use
    /// [`Scanner::try_blocks`] when the scanning mode is not known.
    pub fn blocks(&self) -> &Blocks {
        self.try_blocks()
            .expect("streaming scanners do not retain a block directory")
    }

    /// The retained block directory, or `None` in streaming mode.
    pub fn try_blocks(&self) -> Option<&Blocks> {
        self.keep_blocks.then_some(&self.blocks)
    }

    /// The unread bytes already inflated (at least what the last
    /// [`Scanner::fill`] made available).
    pub fn peek(&self) -> &[u8] {
        // SAFETY: only advance changes pos, after fill has established
        // pos + n <= text.len(). Compaction subtracts the same prefix
        // from text and pos, so pos remains within this live buffer.
        unsafe { self.text.get_unchecked(self.pos..) }
    }

    /// Consumes the next `n` bytes and returns where they lie, for
    /// [`Scanner::slice`]; [`Error::Truncated`] when the stream ends first.
    pub fn take_range(&mut self, n: usize, what: &'static str) -> Result<Range<usize>> {
        let have = self.fill(n)?;
        if have < n {
            return Err(Error::Truncated {
                what,
                need: n as u64,
                have: have as u64,
            });
        }
        let at = self.pos;
        self.advance(n);
        Ok(at..at + n)
    }

    /// Consumes the next `n` bytes; [`Error::Truncated`] when the stream
    /// ends first.
    pub fn take(&mut self, n: usize, what: &'static str) -> Result<&[u8]> {
        let range = self.take_range(n, what)?;
        // SAFETY: take_range succeeds only after fill makes n bytes
        // available at its start. advance changes only position metadata,
        // so the entire returned range still lies in initialized text.
        Ok(unsafe { self.text.get_unchecked(range) })
    }

    /// Consumes `n` bytes without assembling them into one contiguous
    /// frame. On truncation the available prefix has been consumed, and
    /// the error reports its length; a decoding error remains sticky.
    pub fn skip(&mut self, n: usize, what: &'static str) -> Result<()> {
        if n == 0 {
            self.fill(0)?;
            return Ok(());
        }
        let mut left = n;
        while left > 0 {
            let have = self.fill(1)?;
            if have == 0 {
                return Err(Error::Truncated {
                    what,
                    need: n as u64,
                    have: (n - left) as u64,
                });
            }
            let consumed = left.min(have);
            self.advance(consumed);
            left -= consumed;
        }
        Ok(())
    }

    /// The bytes of a range [`Scanner::take_range`] returned, valid until
    /// the next read.
    pub fn slice(&self, range: Range<usize>) -> &[u8] {
        &self.text[range]
    }
}

/// A queue on a pool of its own, which lives as long as the queue.
struct OwnPool {
    queue: PoolQueue<Vec<u8>>,
    _pool: Pool,
}

impl Ordered<Vec<u8>> for OwnPool {
    fn submit(&mut self, job: Job<Vec<u8>>) -> Result<()> {
        self.queue.submit(job)
    }

    fn next(&mut self, wait: bool) -> Result<Option<Vec<u8>>> {
        self.queue.next(wait)
    }

    fn pending(&self) -> usize {
        self.queue.pending()
    }

    fn full(&self) -> bool {
        self.queue.full()
    }
}

/// An inflating queue of `threads` workers of its own (inline for one).
pub fn inflate_queue(threads: usize) -> Box<dyn Ordered<Vec<u8>>> {
    if threads <= 1 {
        return Box::new(Inline::new());
    }
    let pool = Pool::new(threads);
    Box::new(OwnPool {
        queue: pool.ordered(4 * threads),
        _pool: pool,
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{EOF_BLOCK, Writer, deflate_block_at};

    /// The next length-prefixed frame (an `i32` length, then the bytes).
    fn frame<R: Read>(scan: &mut Scanner<R>) -> Option<Range<usize>> {
        if scan.fill(4).unwrap() == 0 {
            return None;
        }
        let n = i32::from_le_bytes(scan.take(4, "frame").unwrap().try_into().unwrap());
        Some(scan.take_range(n as usize, "frame").unwrap())
    }

    fn framed(payload: &[u8]) -> Vec<u8> {
        let mut out = (payload.len() as i32).to_le_bytes().to_vec();
        out.extend_from_slice(payload);
        out
    }

    #[test]
    fn empty_prefix_and_buffer_pruning_advance_each_block_once() {
        const EMPTY: usize = 4096;
        const RECORDS: usize = 700;
        const SIZE: usize = 8192;
        let record = vec![b'x'; SIZE];
        let mut text = Vec::new();
        for _ in 0..RECORDS {
            text.extend_from_slice(&framed(&record));
        }
        assert!(text.len() > KEEP);
        let mut writer = Writer::new(Vec::new(), 1, Box::new(Inline::new())).unwrap();
        writer.write(&text).unwrap();
        let compressed = writer.finish().unwrap().out;
        let mut raw = EOF_BLOCK.repeat(EMPTY);
        raw.extend_from_slice(&compressed);
        for keep_blocks in [true, false] {
            let mut scan = Scanner::with_catalog(&raw[..], Box::new(Inline::new()), keep_blocks);
            let mut records = 0;
            while let Some(range) = frame(&mut scan) {
                assert_eq!(scan.slice(range), record);
                records += 1;
            }
            assert_eq!(records, RECORDS);
            assert!(scan.text.len() < text.len(), "consumed prefix was pruned");
            assert!(scan.cursor_steps <= scan.taken as usize);
            assert!(scan.loaded.iter().all(|block| block.len != 0));
            assert_eq!(scan.try_blocks().is_some(), keep_blocks);
        }
    }

    #[test]
    fn streaming_positions_match_the_directory_without_retaining_it() {
        const EMPTY_PREFIX: usize = 4096;
        const RECORDS: usize = 8192;
        const SIZE: usize = 1024;
        let record = vec![b'x'; SIZE];
        let one = framed(&record);
        let mut member = Vec::new();
        deflate_block_at(1, &one, &mut member).unwrap();
        let mut raw = EOF_BLOCK.repeat(EMPTY_PREFIX);
        for _ in 0..RECORDS {
            raw.extend_from_slice(&member);
            raw.extend_from_slice(&EOF_BLOCK);
        }
        let mut indexed = Scanner::new(&raw[..], Box::new(Inline::new()));
        let mut streaming = Scanner::new_streaming(&raw[..], Box::new(Inline::new()));
        let mut peak_loaded = 0;
        for _ in 0..RECORDS {
            let a = frame(&mut indexed).unwrap();
            let b = frame(&mut streaming).unwrap();
            assert_eq!(indexed.slice(a), record);
            assert_eq!(streaming.slice(b), record);
            assert_eq!(streaming.tell(), indexed.tell());
            assert_eq!(streaming.compressed_read(), indexed.compressed_read());
            assert_eq!(streaming.blocks.starts.capacity(), 0);
            peak_loaded = peak_loaded.max(streaming.loaded.len());
        }
        assert!(frame(&mut indexed).is_none());
        assert!(frame(&mut streaming).is_none());
        assert_eq!(streaming.compressed_read(), raw.len() as u64);
        assert!(streaming.try_blocks().is_none());
        assert_eq!(indexed.blocks().starts.len(), EMPTY_PREFIX + RECORDS * 2);
        assert!(peak_loaded <= KEEP / one.len() + 2);
        assert!(streaming.text.len() < RECORDS * one.len());
    }

    #[test]
    fn empty_prefix_does_not_fill_the_live_block_queue() {
        const EMPTY: usize = 20_000;
        let mut raw = EOF_BLOCK.repeat(EMPTY);
        deflate_block_at(1, b"head", &mut raw).unwrap();
        let mut scanner = Scanner::new_streaming(&raw[..], Box::new(Inline::new()));
        assert_eq!(scanner.take(4, "head").unwrap(), b"head");
        assert_eq!(scanner.loaded.len(), 1);
        assert_eq!(scanner.blocks.starts.capacity(), 0);
        assert_eq!(scanner.tell(), Logical::new(EMPTY as u64 + 1, 0));
        assert_eq!(scanner.fill(1).unwrap(), 0);
    }

    #[test]
    fn byte_sized_members_keep_only_the_active_frame_metadata() {
        const RECORDS: usize = 1024;
        let data = vec![b'x'; 64];
        let one = framed(&data);
        let members: Vec<Vec<u8>> = (0..=255)
            .map(|byte| {
                let mut encoded = Vec::new();
                deflate_block_at(1, &[byte], &mut encoded).unwrap();
                encoded
            })
            .collect();
        let mut raw = Vec::new();
        for _ in 0..RECORDS {
            for &byte in &one {
                raw.extend_from_slice(&members[byte as usize]);
            }
        }
        let mut indexed = Scanner::new(&raw[..], Box::new(Inline::new()));
        let mut streaming = Scanner::new_streaming(&raw[..], Box::new(Inline::new()));
        for _ in 0..RECORDS {
            let range = frame(&mut indexed).unwrap();
            assert_eq!(indexed.slice(range), data);
            let range = frame(&mut streaming).unwrap();
            assert_eq!(streaming.slice(range), data);
            assert_eq!(streaming.tell(), indexed.tell());
            assert_eq!(streaming.loaded.len(), 1);
        }
        assert!(frame(&mut indexed).is_none());
        assert!(frame(&mut streaming).is_none());
        assert!(streaming.peak_loaded <= one.len() + 1);
        assert!(indexed.peak_loaded <= one.len() + 1);
        assert_eq!(indexed.blocks().starts.len(), RECORDS * one.len());
        assert!(streaming.try_blocks().is_none());
    }

    #[test]
    fn short_takes_after_a_truncated_fill_keep_unread_frame_positions() {
        let data: Vec<u8> = (0..30).collect();
        let mut raw = Vec::new();
        for part in data.chunks(10) {
            deflate_block_at(1, part, &mut raw).unwrap();
        }
        for keep_blocks in [false, true] {
            let mut scan = Scanner::with_catalog(&raw[..], Box::new(Inline::new()), keep_blocks);
            assert!(scan.take(100, "test input").is_err());
            assert_eq!(scan.take(4, "test input").unwrap(), &data[..4]);
            assert_eq!(scan.tell(), Logical::new(0, 4));
            assert_eq!(scan.take(10, "test input").unwrap(), &data[4..14]);
            assert_eq!(scan.tell(), Logical::new(1, 4));
            assert_eq!(scan.take(16, "test input").unwrap(), &data[14..]);
            assert_eq!(scan.tell(), Logical::new(3, 0));
        }
    }

    #[test]
    fn skipping_a_large_frame_keeps_only_a_bounded_text_window() {
        let text = vec![b'x'; 3 * KEEP + 17];
        let mut writer = Writer::new(Vec::new(), 1, Box::new(Inline::new())).unwrap();
        writer.write(&text).unwrap();
        let raw = writer.finish().unwrap().out;
        for keep_blocks in [false, true] {
            let mut skip = Scanner::with_catalog(&raw[..], Box::new(Inline::new()), keep_blocks);
            let mut take = Scanner::with_catalog(&raw[..], Box::new(Inline::new()), keep_blocks);
            skip.skip(text.len() - 1, "frame").unwrap();
            take.take(text.len() - 1, "frame").unwrap();
            assert_eq!(skip.tell(), take.tell());
            assert!(skip.text.len() <= KEEP + crate::MAX_BLOCK_LEN);
            assert_eq!(skip.take(1, "tail").unwrap(), b"x");
            assert_eq!(take.take(1, "tail").unwrap(), b"x");
            assert_eq!(skip.tell(), take.tell());
            assert_eq!(skip.fill(1).unwrap(), 0);
            assert_eq!(take.fill(1).unwrap(), 0);
            if keep_blocks {
                assert_eq!(skip.blocks().starts, take.blocks().starts);
                assert_eq!(skip.blocks().text_len, take.blocks().text_len);
            }
        }
    }

    #[test]
    fn truncated_skips_consume_the_prefix_and_decoding_errors_stay_sticky() {
        let mut raw = Vec::new();
        for bytes in [b"abc", b"def", b"ghi"] {
            deflate_block_at(1, bytes, &mut raw).unwrap();
        }
        for keep_blocks in [false, true] {
            let mut scan = Scanner::with_catalog(&raw[..], Box::new(Inline::new()), keep_blocks);
            scan.skip(4, "frame").unwrap();
            assert_eq!(scan.tell(), Logical::new(1, 1));
            assert_eq!(scan.take(2, "frame").unwrap(), b"ef");
            assert_eq!(
                scan.skip(10, "frame").unwrap_err(),
                Error::Truncated {
                    what: "frame",
                    need: 10,
                    have: 3,
                }
            );
            assert_eq!(scan.tell(), Logical::new(3, 0));
            scan.skip(0, "frame").unwrap();
            assert_eq!(scan.fill(1).unwrap(), 0);
        }
        let crc = raw.len() - crate::TRAILER_LEN;
        raw[crc] ^= 1;
        let mut scan = Scanner::new_streaming(&raw[..], Box::new(Inline::new()));
        let error = scan.skip(9, "frame").unwrap_err();
        assert_eq!(scan.skip(1, "frame").unwrap_err(), error);
        assert_eq!(scan.skip(0, "frame").unwrap_err(), error);
        assert_eq!(scan.take(1, "frame").unwrap_err(), error);
    }
}
