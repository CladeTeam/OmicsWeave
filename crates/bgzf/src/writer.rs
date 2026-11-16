//! [`Writer`]: BGZF output cut as htslib cuts it, compressed through an
//! ordered queue.

use std::collections::VecDeque;
use std::io::Write;

use omicsweave_core::exec::Ordered;
use omicsweave_core::{Error, Result};

use crate::{BLOCK_TEXT_LEN, EOF_BLOCK, check_level, deflate_block_at};

/// A position in the text of a file being written: the block's sequence
/// number in the high 48 bits, the byte within its text in the low 16 —
/// a virtual offset with the block's number in place of its compressed
/// address, which is not known until the block is compressed and written.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct Logical(pub u64);

impl Logical {
    pub fn new(block: u64, within: u16) -> Logical {
        Logical(block << 16 | within as u64)
    }

    pub fn block(self) -> u64 {
        self.0 >> 16
    }

    pub fn within(self) -> u16 {
        self.0 as u16
    }
}

/// Where each written block starts, compressed and in the text.
#[derive(Debug, Clone, Default, PartialEq, Eq)]
pub struct Blocks {
    /// `(compressed offset, text offset)` of each block, in order.
    pub starts: Vec<(u64, u64)>,
    /// Compressed offset of the end-of-file block (the length of the data
    /// blocks).
    pub end: u64,
    /// Text bytes written.
    pub text_len: u64,
}

impl Blocks {
    /// The virtual offset of a [`Logical`] position. A position at the start
    /// of the block after the last one is the end-of-file block's address,
    /// as `bgzf_tell` reports it after the last flush.
    pub fn resolve(&self, at: Logical) -> u64 {
        let block = at.block() as usize;
        let address = match self.starts.get(block) {
            Some(&(c, _)) => c,
            None => self.end,
        };
        address << 16 | at.within() as u64
    }
}

/// The sink and block layout of a finished [`Writer`].
pub struct Finished<W> {
    pub out: W,
    pub blocks: Blocks,
}

/// Writes BGZF. Text is cut into blocks exactly as htslib's `bgzf_write`,
/// `bgzf_flush_try` and `bgzf_flush` cut it, so with the same calls the
/// uncompressed blocks, and every position [`Writer::tell`] reports, match
/// htslib's. Blocks compress on the queue given and are written in order,
/// so the output is the same for any queue.
pub struct Writer<W: Write> {
    out: W,
    level: i32,
    exec: Box<dyn Ordered<Vec<u8>>>,
    text: Vec<u8>,
    /// Blocks cut so far: the sequence number of the block being filled.
    cut: u64,
    /// Text offsets of blocks cut but not yet written.
    queued: VecDeque<u64>,
    text_pos: u64,
    blocks: Blocks,
    catalog: bool,
    failed: bool,
}

impl<W: Write> Writer<W> {
    /// A writer at libdeflate compression `level` (0 to 12) compressing on
    /// `exec`.
    pub fn new(out: W, level: i32, exec: Box<dyn Ordered<Vec<u8>>>) -> Result<Writer<W>> {
        Self::with_catalog(out, level, exec, true)
    }

    /// A writer that keeps no per-block directory. Use this for outputs
    /// that do not need virtual-offset indexes and finish it with
    /// [`Writer::finish_streaming`]. In-flight buffers still follow the
    /// supplied queue's limit; metadata does not grow with the file.
    pub fn new_streaming(out: W, level: i32, exec: Box<dyn Ordered<Vec<u8>>>) -> Result<Writer<W>> {
        Self::with_catalog(out, level, exec, false)
    }

    fn with_catalog(
        out: W,
        level: i32,
        exec: Box<dyn Ordered<Vec<u8>>>,
        catalog: bool,
    ) -> Result<Writer<W>> {
        check_level(level)?;
        Ok(Writer {
            out,
            level,
            exec,
            text: Vec::with_capacity(BLOCK_TEXT_LEN),
            cut: 0,
            queued: VecDeque::new(),
            text_pos: 0,
            blocks: Blocks::default(),
            catalog,
            failed: false,
        })
    }

    /// `bgzf_write`: appends `data`, cutting a block whenever one fills.
    pub fn write(&mut self, mut data: &[u8]) -> Result<()> {
        self.check()?;
        while !data.is_empty() {
            let take = (BLOCK_TEXT_LEN - self.text.len()).min(data.len());
            self.text.extend_from_slice(&data[..take]);
            data = &data[take..];
            if self.text.len() == BLOCK_TEXT_LEN {
                self.cut_block()?;
            }
        }
        Ok(())
    }

    /// `bgzf_flush_try(len)`: starts a new block unless `len` more bytes
    /// fit in the current one.
    pub fn flush_try(&mut self, len: usize) -> Result<()> {
        self.check()?;
        if len > BLOCK_TEXT_LEN - self.text.len() {
            self.cut_block()?;
        }
        Ok(())
    }

    /// [`Writer::flush_try`] then [`Writer::write`], as BAM records are
    /// written.
    pub fn write_record(&mut self, data: &[u8]) -> Result<()> {
        self.flush_try(data.len())?;
        self.write(data)
    }

    /// `bgzf_flush`: ends the current block, if it holds any text.
    pub fn flush_block(&mut self) -> Result<()> {
        self.check()?;
        self.cut_block()
    }

    /// Ends the current block, writes every pending block and flushes the
    /// sink.
    pub fn sync(&mut self) -> Result<()> {
        let r = self
            .finish_blocks()
            .and_then(|()| self.out.flush().map_err(Into::into));
        if r.is_err() {
            self.failed = true;
        }
        r
    }

    /// `bgzf_tell`: the position the next byte will be written at.
    pub fn tell(&self) -> Logical {
        Logical::new(self.cut, self.text.len() as u16)
    }

    /// Text bytes written so far.
    pub fn text_len(&self) -> u64 {
        self.text_pos + self.text.len() as u64
    }

    fn check(&self) -> Result<()> {
        if self.failed {
            return Err(Error::invalid("bgzf writer", "an earlier write failed"));
        }
        Ok(())
    }

    fn cut_block(&mut self) -> Result<()> {
        if self.text.is_empty() {
            return Ok(());
        }
        let r = self.submit();
        if r.is_err() {
            self.failed = true;
        }
        r
    }

    fn submit(&mut self) -> Result<()> {
        while self.exec.full() {
            self.take(true)?;
        }
        let text = std::mem::take(&mut self.text);
        self.queued.push_back(self.text_pos);
        self.text_pos += text.len() as u64;
        self.cut += 1;
        let level = self.level;
        self.exec.submit(Box::new(move || {
            let mut block = Vec::new();
            deflate_block_at(level, &text, &mut block)?;
            Ok(block)
        }))?;
        while self.take(false)? {}
        self.text.reserve_exact(BLOCK_TEXT_LEN);
        Ok(())
    }

    /// Writes the next compressed block if it is ready (or, with `wait`,
    /// once it is); returns whether one was written.
    fn take(&mut self, wait: bool) -> Result<bool> {
        let Some(mut block) = self.exec.next(wait)? else {
            return Ok(false);
        };
        let text = self.queued.pop_front().expect("a queued block per result");
        self.out.write_all(&block)?;
        if self.catalog {
            self.blocks.starts.push((self.blocks.end, text));
        }
        self.blocks.end += block.len() as u64;
        if self.text.is_empty() && block.capacity() > self.text.capacity() {
            // The sink has consumed these compressed bytes. Their allocation
            // can become the next text buffer without retaining their data.
            block.clear();
            self.text = block;
        }
        Ok(true)
    }

    /// Writes the last block and the end-of-file block, and returns the sink
    /// and the layout of the blocks.
    pub fn finish(mut self) -> Result<Finished<W>> {
        if !self.catalog {
            return Err(Error::invalid(
                "bgzf writer",
                "block directory disabled; use finish_streaming",
            ));
        }
        self.complete()?;
        Ok(Finished {
            out: self.out,
            blocks: self.blocks,
        })
    }

    /// Writes the final data block and EOF marker and flushes the sink.
    /// This also accepts a catalog writer, discarding its block directory.
    pub fn finish_streaming(mut self) -> Result<W> {
        self.complete()?;
        Ok(self.out)
    }

    fn complete(&mut self) -> Result<()> {
        let r = self.finish_blocks();
        if r.is_err() {
            self.failed = true;
        }
        r?;
        self.out.write_all(&EOF_BLOCK)?;
        self.out.flush()?;
        self.blocks.text_len = self.text_pos;
        Ok(())
    }

    fn finish_blocks(&mut self) -> Result<()> {
        self.check()?;
        self.cut_block()?;
        while self.take(true)? {}
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::Reader;
    use omicsweave_core::exec::{Inline, Pool};
    use std::io::Read;

    fn write_all(
        exec: Box<dyn Ordered<Vec<u8>>>,
        records: &[Vec<u8>],
    ) -> (Vec<u8>, Vec<Logical>, Blocks) {
        let mut w = Writer::new(Vec::new(), 6, exec).unwrap();
        let mut tells = vec![w.tell()];
        for (i, r) in records.iter().enumerate() {
            if i % 3 == 0 {
                w.write_record(r).unwrap();
            } else {
                w.write(r).unwrap();
            }
            if i % 17 == 5 {
                w.flush_block().unwrap();
            }
            tells.push(w.tell());
        }
        let f = w.finish().unwrap();
        (f.out, tells, f.blocks)
    }

    fn records() -> Vec<Vec<u8>> {
        (0..400u32)
            .map(|i| {
                let n = (i * 7919 % 3000) as usize + if i % 50 == 0 { 70_000 } else { 0 };
                (0..n)
                    .map(|j| (j as u32 * 31 + i) as u8 % 7 + b'a')
                    .collect()
            })
            .collect()
    }

    #[test]
    fn output_does_not_depend_on_the_queue() {
        let recs = records();
        let (a, ta, ba) = write_all(Box::new(Inline::new()), &recs);
        let pool = Pool::new(4);
        let (b, tb, bb) = write_all(Box::new(pool.ordered(3)), &recs);
        assert_eq!(a, b);
        assert_eq!(ta, tb);
        assert_eq!(ba, bb);
        let mut text = Vec::new();
        Reader::new(&a[..], Box::new(Inline::new()))
            .read_to_end(&mut text)
            .unwrap();
        assert_eq!(text, recs.concat());
        assert_eq!(ba.text_len, text.len() as u64);
    }

    #[test]
    fn positions_resolve_to_the_text_they_name() {
        let recs = records();
        let (file, tells, blocks) = write_all(Box::new(Inline::new()), &recs);
        let text = recs.concat();
        let mut at = 0usize;
        for (i, t) in tells.iter().enumerate() {
            let v = blocks.resolve(*t);
            let (addr, within) = ((v >> 16) as usize, (v & 0xffff) as usize);
            // The block at `addr` holds text from its recorded start.
            let start = blocks
                .starts
                .iter()
                .find(|s| s.0 as usize == addr)
                .map(|s| s.1 as usize)
                .unwrap_or(text.len());
            assert_eq!(start + within, at, "tell {i}");
            assert!(addr <= file.len());
            if i < recs.len() {
                at += recs[i].len();
            }
        }
    }

    #[test]
    fn streaming_directory_capacity_does_not_grow_with_flushed_blocks() {
        for streaming in [false, true] {
            let queue = Box::new(Inline::new());
            let mut writer = if streaming {
                Writer::new_streaming(std::io::sink(), 0, queue)
            } else {
                Writer::new(std::io::sink(), 0, queue)
            }
            .unwrap();
            for i in 0..4096 {
                writer.write(b"x").unwrap();
                writer.flush_block().unwrap();
                assert_eq!(writer.tell(), Logical::new(i + 1, 0));
                assert_eq!(writer.text_len(), i + 1);
            }
            if streaming {
                assert_eq!(writer.blocks.starts.len(), 0);
                assert_eq!(writer.blocks.starts.capacity(), 0);
            } else {
                assert_eq!(writer.blocks.starts.len(), 4096);
                assert!(writer.blocks.starts.capacity() >= 4096);
            }
            writer.finish_streaming().unwrap();
        }
    }
}
