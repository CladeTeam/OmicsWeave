//! A whole BAM read front to back: blocks are inflated ahead of the reader
//! on an [`Ordered`] queue, and each record comes with the position htslib's
//! `bgzf_tell` reports after it. Index building uses it.

use std::collections::VecDeque;
use std::fs::File;
use std::io::{BufReader, Read};
use std::ops::Range;
use std::path::{Path, PathBuf};

use omicsweave_bgzf::{
    BinningIndex, Blocks, IndexFormat, IndexKind, Logical, inflate_block_owned, read_block,
};
use omicsweave_core::exec::{Inline, Job, Ordered, Pool, PoolQueue};
use omicsweave_core::{Error, Result, sidecar};

use crate::header::{BamHeader, Take};
use crate::record::{Record, flag};

/// Consumed text kept before it is dropped from the front.
const KEEP: usize = 4 << 20;

#[derive(Debug, Clone, Copy)]
struct Inflated {
    /// The block's number in the file.
    seq: u64,
    text: usize,
    len: usize,
}

/// The records of a BAM stream in file order.
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

    fn with_catalog(input: R, exec: Box<dyn Ordered<Vec<u8>>>, keep_blocks: bool) -> Scanner<R> {
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

    fn fill(&mut self, n: usize) -> Result<usize> {
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
            .expect("streaming BAM scanners do not retain a block directory")
    }

    /// The retained block directory, or `None` in streaming mode.
    pub fn try_blocks(&self) -> Option<&Blocks> {
        self.keep_blocks.then_some(&self.blocks)
    }

    /// Reads the header.
    pub fn header(&mut self) -> Result<BamHeader> {
        BamHeader::read(self)
    }

    /// The next record after its `block_size`; `None` at a clean end.
    pub fn next_record(&mut self) -> Result<Option<Range<usize>>> {
        let have = self.fill(4)?;
        if have == 0 {
            return Ok(None);
        }
        if have < 4 {
            return Err(Error::Truncated {
                what: "BAM record",
                need: 4,
                have: have as u64,
            });
        }
        let size = i32::from_le_bytes(self.text[self.pos..self.pos + 4].try_into().unwrap());
        if size < 32 {
            return Err(Error::invalid(
                "BAM record",
                format!("block_size {size} below the 32 bytes of fixed fields"),
            ));
        }
        let total = 4 + size as usize;
        let have = self.fill(total)?;
        if have < total {
            return Err(Error::Truncated {
                what: "BAM record",
                need: total as u64,
                have: have as u64,
            });
        }
        let start = self.pos + 4;
        self.advance(total);
        Ok(Some(start..start + size as usize))
    }

    /// The text of [`Scanner::next_record`]'s range.
    pub fn slice(&self, range: Range<usize>) -> &[u8] {
        &self.text[range]
    }
}

impl<R: Read> Take for Scanner<R> {
    fn take(&mut self, n: usize, what: &'static str) -> Result<&[u8]> {
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
        Ok(&self.text[at..at + n])
    }
}

/// Which index to build.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum IndexChoice {
    /// BAI when every reference fits (at most 2^29 bases), else CSI.
    Auto,
    Bai,
    /// CSI with bins of `2^min_shift` bases.
    Csi {
        min_shift: u32,
    },
}

impl IndexChoice {
    /// The choice called `kind` (`auto`, `bai` or `csi`); `min_shift` sizes
    /// the bins of a CSI.
    pub fn from_name(kind: &str, min_shift: u32) -> Result<IndexChoice> {
        match kind {
            "auto" => Ok(IndexChoice::Auto),
            "bai" => Ok(IndexChoice::Bai),
            "csi" => Ok(IndexChoice::Csi { min_shift }),
            _ => Err(Error::invalid(
                "bam index",
                format!("unknown index {kind:?}: use auto, bai or csi"),
            )),
        }
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

/// [`build_index`] of the BAM file at `bam`, inflated on `threads`
/// workers, saved beside it as `<bam>.bai` or `<bam>.csi` (replacing
/// one there). Returns the index's kind and path.
pub fn build_index_file(
    bam: &Path,
    choice: IndexChoice,
    threads: usize,
) -> Result<(IndexKind, PathBuf)> {
    let input = BufReader::with_capacity(1 << 20, File::open(bam)?);
    let (kind, bytes) = build_index(input, inflate_queue(threads), choice)?;
    let mut path = bam.as_os_str().to_owned();
    path.push(match kind {
        IndexKind::Bai => ".bai",
        IndexKind::Csi => ".csi",
    });
    let path = PathBuf::from(path);
    sidecar::write_file(&path, &bytes)?;
    Ok((kind, path))
}

/// Builds the index of a coordinate-sorted BAM read from `input`, as
/// `samtools index` builds it (`sam_index_build3`); returns its kind and
/// the bytes of the index file (CSI BGZF-compressed, as htslib writes it).
pub fn build_index<R: Read>(
    input: R,
    exec: Box<dyn Ordered<Vec<u8>>>,
    choice: IndexChoice,
) -> Result<(IndexKind, Vec<u8>)> {
    let mut scan = Scanner::new(input, exec);
    let header = scan.header()?;
    let max_len = header
        .references()
        .iter()
        .map(|r| r.len as u64)
        .max()
        .unwrap_or(0);
    let fits_bai = max_len <= 1 << 29;
    let (format, kind) = match choice {
        IndexChoice::Bai if !fits_bai => {
            return Err(Error::invalid(
                "bam index",
                format!("a reference of {max_len} bases does not fit a BAI index; use CSI"),
            ));
        }
        IndexChoice::Bai => (IndexFormat::Bai, IndexKind::Bai),
        IndexChoice::Auto if fits_bai => (IndexFormat::Bai, IndexKind::Bai),
        IndexChoice::Auto => {
            let (min_shift, levels) = omicsweave_bgzf::try_csi_levels(max_len, 14)?;
            (IndexFormat::Csi { min_shift, levels }, IndexKind::Csi)
        }
        IndexChoice::Csi { min_shift } => {
            let (min_shift, levels) = omicsweave_bgzf::try_csi_levels(max_len, min_shift)?;
            (IndexFormat::Csi { min_shift, levels }, IndexKind::Csi)
        }
    };
    let n_refs = header.references().len();
    let mut index = BinningIndex::try_new(format, n_refs, scan.tell())?;
    while let Some(range) = scan.next_record()? {
        let end_at = scan.tell();
        let rec = Record::parse(scan.slice(range), n_refs)?;
        index.push(
            rec.tid,
            rec.pos as i64,
            rec.end(),
            end_at,
            rec.flag & flag::UNMAPPED == 0,
        )?;
    }
    let end = scan.tell();
    Ok((kind, index.finish(end, scan.blocks())?))
}

#[cfg(test)]
mod tests {
    use super::*;
    use omicsweave_bgzf::{EOF_BLOCK, Writer};

    #[test]
    fn empty_prefix_and_buffer_pruning_advance_each_block_once() {
        const EMPTY: usize = 4096;
        const RECORDS: usize = 700;
        const SIZE: usize = 8192;
        let mut text = b"BAM\x01\0\0\0\0\0\0\0\0".to_vec();
        let mut record = vec![b'x'; SIZE];
        record[..32].fill(0);
        record[..4].copy_from_slice(&(-1i32).to_le_bytes());
        record[4..8].copy_from_slice(&(-1i32).to_le_bytes());
        record[8..12].copy_from_slice(&2u32.to_le_bytes());
        record[12..16].copy_from_slice(&((flag::UNMAPPED as u32) << 16).to_le_bytes());
        record[20..24].copy_from_slice(&(-1i32).to_le_bytes());
        record[24..28].copy_from_slice(&(-1i32).to_le_bytes());
        record[32..37].copy_from_slice(b"r\0AAZ");
        record[SIZE - 1] = 0;
        for _ in 0..RECORDS {
            text.extend_from_slice(&(SIZE as i32).to_le_bytes());
            text.extend_from_slice(&record);
        }
        assert!(text.len() > KEEP);
        let mut writer = Writer::new(Vec::new(), 1, Box::new(Inline::new())).unwrap();
        writer.write(&text).unwrap();
        let compressed = writer.finish().unwrap().out;
        let mut raw = EOF_BLOCK.repeat(EMPTY);
        raw.extend_from_slice(&compressed);
        for keep_blocks in [true, false] {
            let mut scan = Scanner::with_catalog(&raw[..], Box::new(Inline::new()), keep_blocks);
            assert!(scan.header().unwrap().references().is_empty());
            let mut records = 0;
            while let Some(range) = scan.next_record().unwrap() {
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
    fn corrupt_member_is_not_reinterpreted_as_clean_eof_after_error() {
        let mut raw = Vec::new();
        omicsweave_bgzf::deflate_block_at(1, b"BAM\x01\0\0\0\0\0\0\0\0", &mut raw).unwrap();
        omicsweave_bgzf::deflate_block_at(1, &[1u8; 100], &mut raw).unwrap();
        let crc = raw.len() - 8;
        raw[crc] ^= 1;
        raw.extend_from_slice(&EOF_BLOCK);
        for streaming in [false, true] {
            for exec in [
                Box::new(Inline::new()) as Box<dyn Ordered<Vec<u8>>>,
                inflate_queue(2),
            ] {
                let mut scan = Scanner::with_catalog(&raw[..], exec, !streaming);
                scan.header().unwrap();
                let first = scan.next_record().unwrap_err();
                assert_eq!(scan.next_record().unwrap_err(), first);
            }
        }
    }

    #[test]
    fn streaming_records_and_positions_match_without_file_wide_metadata() {
        const EMPTY_PREFIX: usize = 4096;
        const RECORDS: usize = 8192;
        const SIZE: usize = 1024;
        let mut record = vec![b'x'; SIZE];
        record[..32].fill(0);
        record[..4].copy_from_slice(&(-1i32).to_le_bytes());
        record[4..8].copy_from_slice(&(-1i32).to_le_bytes());
        record[8..12].copy_from_slice(&2u32.to_le_bytes());
        record[12..16].copy_from_slice(&((flag::UNMAPPED as u32) << 16).to_le_bytes());
        record[20..24].copy_from_slice(&(-1i32).to_le_bytes());
        record[24..28].copy_from_slice(&(-1i32).to_le_bytes());
        record[32..37].copy_from_slice(b"r\0AAZ");
        record[SIZE - 1] = 0;
        let mut framed = (SIZE as i32).to_le_bytes().to_vec();
        framed.extend_from_slice(&record);
        let mut member = Vec::new();
        omicsweave_bgzf::deflate_block_at(1, &framed, &mut member).unwrap();
        let mut raw = EOF_BLOCK.repeat(EMPTY_PREFIX);
        omicsweave_bgzf::deflate_block_at(1, b"BAM\x01\0\0\0\0\0\0\0\0", &mut raw).unwrap();
        for _ in 0..RECORDS {
            raw.extend_from_slice(&member);
            raw.extend_from_slice(&EOF_BLOCK);
        }
        let mut indexed = Scanner::new(&raw[..], Box::new(Inline::new()));
        let mut streaming = Scanner::new_streaming(&raw[..], Box::new(Inline::new()));
        assert_eq!(
            indexed.header().unwrap().text(),
            streaming.header().unwrap().text()
        );
        let mut peak_loaded = 0;
        for _ in 0..RECORDS {
            let a = indexed.next_record().unwrap().unwrap();
            let b = streaming.next_record().unwrap().unwrap();
            assert_eq!(indexed.slice(a), record);
            assert_eq!(streaming.slice(b), record);
            assert_eq!(streaming.tell(), indexed.tell());
            assert_eq!(streaming.compressed_read(), indexed.compressed_read());
            assert_eq!(streaming.blocks.starts.capacity(), 0);
            peak_loaded = peak_loaded.max(streaming.loaded.len());
        }
        assert!(indexed.next_record().unwrap().is_none());
        assert!(streaming.next_record().unwrap().is_none());
        assert_eq!(streaming.compressed_read(), raw.len() as u64);
        assert!(streaming.try_blocks().is_none());
        assert_eq!(
            indexed.blocks().starts.len(),
            EMPTY_PREFIX + 1 + RECORDS * 2
        );
        assert!(peak_loaded <= KEEP / framed.len() + 2);
        assert!(streaming.text.len() < RECORDS * framed.len());
    }

    #[test]
    fn empty_prefix_does_not_fill_the_live_block_queue() {
        const EMPTY: usize = 20_000;
        let mut raw = EOF_BLOCK.repeat(EMPTY);
        omicsweave_bgzf::deflate_block_at(1, b"BAM\x01\0\0\0\0\0\0\0\0", &mut raw).unwrap();
        let mut scanner = Scanner::new_streaming(&raw[..], Box::new(Inline::new()));
        scanner.header().unwrap();
        assert_eq!(scanner.loaded.len(), 1);
        assert_eq!(scanner.blocks.starts.capacity(), 0);
        assert_eq!(scanner.tell(), Logical::new(EMPTY as u64 + 1, 0));
        assert!(scanner.next_record().unwrap().is_none());
    }

    #[test]
    fn byte_sized_members_keep_only_the_active_record_metadata() {
        const RECORDS: usize = 1024;
        let mut data = vec![b'x'; 64];
        data[..32].fill(0);
        data[..4].copy_from_slice(&(-1i32).to_le_bytes());
        data[4..8].copy_from_slice(&(-1i32).to_le_bytes());
        data[8..12].copy_from_slice(&2u32.to_le_bytes());
        data[12..16].copy_from_slice(&((flag::UNMAPPED as u32) << 16).to_le_bytes());
        data[20..24].copy_from_slice(&(-1i32).to_le_bytes());
        data[24..28].copy_from_slice(&(-1i32).to_le_bytes());
        data[32..37].copy_from_slice(b"r\0AAZ");
        data[63] = 0;
        let mut framed = 64i32.to_le_bytes().to_vec();
        framed.extend_from_slice(&data);
        let members: Vec<Vec<u8>> = (0..=255)
            .map(|byte| {
                let mut encoded = Vec::new();
                omicsweave_bgzf::deflate_block_at(1, &[byte], &mut encoded).unwrap();
                encoded
            })
            .collect();
        let mut raw = Vec::new();
        omicsweave_bgzf::deflate_block_at(1, b"BAM\x01\0\0\0\0\0\0\0\0", &mut raw).unwrap();
        for _ in 0..RECORDS {
            for &byte in &framed {
                raw.extend_from_slice(&members[byte as usize]);
            }
        }
        let mut indexed = Scanner::new(&raw[..], Box::new(Inline::new()));
        let mut streaming = Scanner::new_streaming(&raw[..], Box::new(Inline::new()));
        indexed.header().unwrap();
        streaming.header().unwrap();
        for _ in 0..RECORDS {
            let range = indexed.next_record().unwrap().unwrap();
            assert_eq!(indexed.slice(range), data);
            let range = streaming.next_record().unwrap().unwrap();
            assert_eq!(streaming.slice(range), data);
            assert_eq!(streaming.tell(), indexed.tell());
            assert_eq!(streaming.loaded.len(), 1);
        }
        assert!(indexed.next_record().unwrap().is_none());
        assert!(streaming.next_record().unwrap().is_none());
        assert!(streaming.peak_loaded <= framed.len() + 1);
        assert!(indexed.peak_loaded <= framed.len() + 1);
        assert_eq!(indexed.blocks().starts.len(), RECORDS * framed.len() + 1);
        assert!(streaming.try_blocks().is_none());
    }

    #[test]
    fn short_takes_after_a_truncated_fill_keep_unread_frame_positions() {
        let data: Vec<u8> = (0..30).collect();
        let mut raw = Vec::new();
        for part in data.chunks(10) {
            omicsweave_bgzf::deflate_block_at(1, part, &mut raw).unwrap();
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
}
