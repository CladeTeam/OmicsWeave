//! [`SeekReader`]: the text of a BGZF file read at virtual offsets, as
//! htslib's `bgzf_seek`, `bgzf_read` and `bgzf_tell` see it, over a
//! [`ReadAt`] source. Compressed bytes are fetched in windows sized by the
//! caller's hint, so a chunk of an index query is one positional read.

use std::borrow::Cow;
use std::collections::HashSet;
use std::sync::Arc;

use omicsweave_core::cache::{BlockCache, BlockKey, CacheLease, SourceToken};
use omicsweave_core::{Error, ReadAt, Result};

use crate::{BlockHeader, Chunk, Inflater, MAX_BLOCK_LEN, block_text_len};

/// Minimum read window for an indexed chunk.
const MIN_WINDOW: u64 = MAX_BLOCK_LEN as u64;
/// Read ahead for an unhinted header or sequential scan.
const DEFAULT_WINDOW: u64 = 256 << 10;
/// Most compressed bytes fetched by one read.
const MAX_WINDOW: u64 = 16 << 20;
/// Minimum consumed prefix considered for amortized compaction.
const KEEP: usize = 4 << 20;
const MAX_ADDRESS: u64 = u64::MAX >> 16;

fn window(address: u64, hint: u64) -> u64 {
    if hint == 0 {
        return DEFAULT_WINDOW;
    }
    let last = u64::from(hint & 0xffff != 0) * MAX_BLOCK_LEN as u64;
    (hint >> 16)
        .saturating_sub(address)
        .saturating_add(last)
        .clamp(MIN_WINDOW, MAX_WINDOW)
}

/// Coalesces sorted physical coverage while preserving the first-seen key
/// order. Unsorted requests retain their caller-specified traversal order.
fn cache_chunks(chunks: &[Chunk]) -> Result<Cow<'_, [Chunk]>> {
    if chunks.iter().any(|chunk| chunk.start > chunk.end) {
        return Err(Error::invalid("bgzf chunk", "decreasing virtual offsets"));
    }
    if !chunks.windows(2).all(|pair| pair[0].start <= pair[1].start) {
        return Ok(Cow::Borrowed(chunks));
    }
    if !chunks
        .windows(2)
        .any(|pair| pair[1].start >> 16 <= pair[0].end >> 16)
    {
        return Ok(Cow::Borrowed(chunks));
    }
    let mut merged: Vec<Chunk> = Vec::with_capacity(chunks.len());
    for &chunk in chunks {
        if chunk.start == chunk.end {
            continue;
        }
        if let Some(last) = merged.last_mut()
            && chunk.start >> 16 <= last.end >> 16
        {
            last.end = last.end.max(chunk.end);
        } else {
            merged.push(chunk);
        }
    }
    Ok(Cow::Owned(merged))
}

fn check_chunk_starts(chunks: &[Chunk], cursor: &mut usize, at: u64, next: u64) -> Result<()> {
    while let Some(chunk) = chunks.get(*cursor) {
        if chunk.start == chunk.end {
            *cursor += 1;
            continue;
        }
        let start = chunk.start >> 16;
        if start >= next {
            break;
        }
        if start != at {
            return Err(Error::invalid(
                "bgzf chunk",
                "start is not a member boundary",
            ));
        }
        *cursor += 1;
    }
    Ok(())
}

/// Actual source reads and decompressor calls performed by one context.
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub struct SeekStats {
    pub read_calls: u64,
    pub read_bytes: u64,
    pub inflate_calls: u64,
    pub compacted_bytes: u64,
    pub discarded_bytes: u64,
}

struct Member {
    key: BlockKey,
    next: u64,
    bytes: Vec<u8>,
}

#[derive(Debug, Clone, Copy)]
struct Block {
    address: u64,
    /// Address of the next physical block.
    end_address: u64,
    /// Where its text starts in [`SeekReader::text`].
    text: usize,
    len: usize,
}

/// Inflated text of consecutive blocks with a read position.
pub struct SeekReader {
    inflater: Inflater,
    raw: Vec<u8>,
    raw_at: u64,
    text: Vec<u8>,
    blocks: Vec<Block>,
    pos: usize,
    /// The block holding the last byte consumed, or the block sought.
    cur: usize,
    /// Address of the block after the last loaded one.
    next: u64,
    /// `bgzf_tell` at `pos`.
    here: u64,
    window: u64,
    cache: Option<Arc<BlockCache>>,
    namespace: SourceToken,
    stats: SeekStats,
}

impl Default for SeekReader {
    fn default() -> SeekReader {
        SeekReader::new()
    }
}

impl SeekReader {
    pub fn new() -> SeekReader {
        SeekReader {
            inflater: Inflater::new(),
            raw: Vec::new(),
            raw_at: 0,
            text: Vec::new(),
            blocks: Vec::new(),
            pos: 0,
            cur: 0,
            next: 0,
            here: 0,
            window: DEFAULT_WINDOW,
            cache: None,
            namespace: SourceToken::new(),
            stats: SeekStats::default(),
        }
    }

    /// Bind this context to an immutable source session and external cache.
    /// Clearing buffered data also makes rebinding to another source safe.
    pub fn set_cache(&mut self, cache: Option<Arc<BlockCache>>, namespace: SourceToken) {
        self.reset_query();
        self.cache = cache;
        self.namespace = namespace;
    }

    pub fn cache_namespace(&self) -> SourceToken {
        self.namespace
    }

    /// Starts a separate query without retaining compressed or inflated data.
    /// Capacities and the external cache remain available for reuse.
    pub fn reset_query(&mut self) {
        self.raw.clear();
        self.raw_at = 0;
        self.text.clear();
        self.blocks.clear();
        self.pos = 0;
        self.cur = 0;
        self.next = 0;
        self.here = 0;
        self.window = DEFAULT_WINDOW;
    }

    /// Shares the source namespace/cache with independent decoder and buffers.
    pub fn clone_context(&self) -> Self {
        let mut context = Self::new();
        context.cache = self.cache.clone();
        context.namespace = self.namespace;
        context
    }

    pub fn stats(&self) -> SeekStats {
        self.stats
    }

    pub fn reset_stats(&mut self) {
        self.stats = SeekStats::default();
    }

    fn key(&self, address: u64) -> BlockKey {
        BlockKey::new(self.namespace, 5, [0, 0], address, 1)
    }

    /// The virtual offset of the next byte, as `bgzf_tell` reports it: a
    /// read that ends a block leaves it at the next block's address.
    pub fn tell(&self) -> u64 {
        self.here
    }

    /// Makes `[at, at + len)` of the file available in `raw` and returns
    /// where it starts there.
    fn fetch<S: ReadAt + ?Sized>(&mut self, source: &S, at: u64, len: usize) -> Result<usize> {
        let have = at >= self.raw_at
            && at
                .checked_add(len as u64)
                .zip(self.raw_at.checked_add(self.raw.len() as u64))
                .is_some_and(|(end, loaded)| end <= loaded);
        if !have {
            let file_len = source.len();
            if at.checked_add(len as u64).is_none_or(|e| e > file_len) {
                return Err(Error::Truncated {
                    what: "bgzf block",
                    need: at.saturating_add(len as u64),
                    have: file_len,
                });
            }
            let n = (len as u64).max(self.window).min(file_len - at) as usize;
            self.raw.resize(n, 0);
            self.stats.read_calls = self.stats.read_calls.saturating_add(1);
            self.stats.read_bytes = self.stats.read_bytes.saturating_add(n as u64);
            if let Err(error) = source.read_exact_at(&mut self.raw, at) {
                // A failing ReadAt may have written only a prefix. No part of
                // that window is valid for a later seek or load retry.
                self.raw.clear();
                return Err(error);
            }
            self.raw_at = at;
        }
        Ok((at - self.raw_at) as usize)
    }

    fn member<S: ReadAt + ?Sized>(
        &mut self,
        source: &S,
        at: u64,
    ) -> Result<(usize, BlockHeader, u64)> {
        let file_len = source.len();
        if at >= file_len {
            return Err(Error::Truncated {
                what: "bgzf block",
                need: at.saturating_add(1),
                have: file_len,
            });
        }
        let head = (BlockHeader::MIN_LEN as u64).min(file_len - at) as usize;
        let from = self.fetch(source, at, head)?;
        let header = match BlockHeader::parse(&self.raw[from..from + head]) {
            Err(Error::Truncated { need, .. }) => {
                let from = self.fetch(source, at, need as usize)?;
                BlockHeader::parse(&self.raw[from..from + need as usize])?
            }
            r => r?,
        };
        let from = self.fetch(source, at, header.block_len)?;
        let next = at
            .checked_add(header.block_len as u64)
            .filter(|&end| end <= MAX_ADDRESS && end <= file_len)
            .ok_or_else(|| {
                Error::invalid("bgzf block", "block end exceeds virtual-offset range")
            })?;
        Ok((from, header, next))
    }

    fn cached_member<S: ReadAt + ?Sized>(
        &mut self,
        source: &S,
        at: u64,
        cache: &BlockCache,
    ) -> Result<CacheLease<Member>> {
        let key = self.key(at);
        let value = cache.get_or_load(key, || {
            let (from, header, next) = self.member(source, at)?;
            let raw = &self.raw[from..from + header.block_len];
            let mut bytes = vec![0; block_text_len(raw, header)?];
            self.stats.inflate_calls = self.stats.inflate_calls.saturating_add(1);
            let (used, written) = self.inflater.inflate_block(raw, &mut bytes)?;
            if used != header.block_len || written != bytes.len() {
                return Err(Error::invalid(
                    "bgzf block",
                    "framing changed while reading member",
                ));
            }
            let charge = bytes.capacity();
            Ok((Member { key, next, bytes }, charge))
        })?;
        if value.key != key
            || value.next <= at
            || value.next > source.len()
            || value.next > MAX_ADDRESS
            || value.bytes.len() > MAX_BLOCK_LEN
        {
            return Err(Error::invalid(
                "bgzf cache",
                "member does not match source range",
            ));
        }
        Ok(value)
    }

    fn compact(&mut self) {
        if self.cur == 0 {
            return;
        }
        let drop = self.blocks[self.cur].text;
        let kept = self.text.len() - drop;
        // Copy at most as many text bytes as have been consumed since the
        // preceding compaction, even after a large speculative fill.
        if drop < KEEP || drop < kept {
            return;
        }
        self.text.copy_within(drop.., 0);
        self.text.truncate(kept);
        self.blocks.copy_within(self.cur.., 0);
        self.blocks.truncate(self.blocks.len() - self.cur);
        for b in &mut self.blocks {
            b.text -= drop;
        }
        self.pos -= drop;
        self.cur = 0;
        self.stats.compacted_bytes = self.stats.compacted_bytes.saturating_add(kept as u64);
        self.stats.discarded_bytes = self.stats.discarded_bytes.saturating_add(drop as u64);
    }

    /// Inflates the next physical member, or borrows its explicitly cached
    /// contents. Empty members do not grow the text-position index.
    fn load_next<S: ReadAt + ?Sized>(&mut self, source: &S) -> Result<bool> {
        if self.next >= source.len() {
            return Ok(false);
        }
        self.compact();
        let at = self.next;
        let start = self.text.len();
        let cache = self
            .cache
            .as_ref()
            .filter(|c| {
                let (bytes, entries) = c.limits();
                bytes > 0 && entries > 0
            })
            .cloned();
        let (n, end_address) = match cache {
            Some(cache) => {
                let member = self.cached_member(source, at, &cache)?;
                self.text.extend_from_slice(&member.bytes);
                (member.bytes.len(), member.next)
            }
            None => {
                let (from, header, next) = self.member(source, at)?;
                let raw = &self.raw[from..from + header.block_len];
                let size = block_text_len(raw, header)?;
                let end = start
                    .checked_add(size)
                    .ok_or_else(|| Error::invalid("bgzf text", "decoded length overflows"))?;
                self.text.resize(end, 0);
                self.stats.inflate_calls = self.stats.inflate_calls.saturating_add(1);
                match self.inflater.inflate_block(raw, &mut self.text[start..]) {
                    Ok((used, n)) if used == header.block_len && n == size => (n, next),
                    Ok(_) => {
                        self.text.truncate(start);
                        return Err(Error::invalid(
                            "bgzf block",
                            "framing changed while reading member",
                        ));
                    }
                    Err(error) => {
                        self.text.truncate(start);
                        return Err(error);
                    }
                }
            }
        };
        if n != 0 {
            self.blocks.push(Block {
                address: at,
                end_address,
                text: start,
                len: n,
            });
        }
        self.next = end_address;
        Ok(true)
    }

    /// `bgzf_seek` to virtual offset `voffset`. `hint` is the virtual
    /// offset up to which reading is expected to continue (0: unknown),
    /// so the compressed bytes are fetched in one read.
    pub fn seek<S: ReadAt + ?Sized>(&mut self, source: &S, voffset: u64, hint: u64) -> Result<()> {
        let address = voffset >> 16;
        let within = (voffset & 0xffff) as usize;
        self.window = window(address, hint);
        let loaded = self
            .blocks
            .binary_search_by_key(&address, |b| b.address)
            .ok();
        let i = match loaded {
            Some(i) => i,
            None => {
                self.blocks.clear();
                self.text.clear();
                self.cur = 0;
                self.pos = 0;
                self.next = address;
                // htslib skips empty blocks at the target, keeping the
                // offset within the block.
                loop {
                    if !self.load_next(source)? {
                        if within == 0 {
                            self.here = voffset;
                            return Ok(());
                        }
                        return Err(Error::invalid(
                            "bgzf seek",
                            format!("virtual offset {voffset:#x} is past the end of the file"),
                        ));
                    }
                    if self.blocks.last().is_some_and(|b| b.len > 0) {
                        break;
                    }
                }
                self.blocks.len() - 1
            }
        };
        let b = self.blocks[i];
        if within > b.len {
            return Err(Error::invalid(
                "bgzf seek",
                format!("virtual offset {voffset:#x} is beyond its block"),
            ));
        }
        self.cur = i;
        self.pos = b.text + within;
        self.here = voffset;
        Ok(())
    }

    /// Names every physical member touched by the indexed virtual chunks.
    /// This only parses headers, without inflating or materializing records.
    pub fn describe_cache<S: ReadAt + ?Sized>(
        &mut self,
        source: &S,
        chunks: &[Chunk],
        out: &mut Vec<BlockKey>,
    ) -> Result<()> {
        out.clear();
        let mut seen = HashSet::new();
        let original = chunks;
        let chunks = cache_chunks(chunks)?;
        let merged = matches!(&chunks, Cow::Owned(_));
        let mut cursor = 0;
        for chunk in chunks.iter() {
            if chunk.start == chunk.end {
                continue;
            }
            let mut at = chunk.start >> 16;
            self.window = window(at, chunk.end);
            while (at << 16) < chunk.end {
                let key = self.key(at);
                if seen.insert(key) {
                    out.push(key);
                }
                let cached = match &self.cache {
                    Some(cache) => cache.get::<Member>(key)?,
                    None => None,
                };
                let next = match cached {
                    Some(member) => {
                        if member.key != key || member.next <= at || member.next > source.len() {
                            return Err(Error::invalid("bgzf cache", "invalid member boundary"));
                        }
                        member.next
                    }
                    None => self.member(source, at)?.2,
                };
                if merged {
                    check_chunk_starts(original, &mut cursor, at, next)?;
                }
                at = next;
            }
        }
        Ok(())
    }

    /// Synchronously populates the attached cache without producing records.
    pub fn prefetch_cache<S: ReadAt + ?Sized>(
        &mut self,
        source: &S,
        chunks: &[Chunk],
    ) -> Result<()> {
        let cache = self
            .cache
            .clone()
            .ok_or_else(|| Error::invalid("bgzf cache", "no cache attached"))?;
        let (bytes, entries) = cache.limits();
        let original = chunks;
        let chunks = cache_chunks(chunks)?;
        let merged = matches!(&chunks, Cow::Owned(_));
        let mut cursor = 0;
        for chunk in chunks.iter() {
            if bytes == 0 || entries == 0 || chunk.start == chunk.end {
                continue;
            }
            let mut at = chunk.start >> 16;
            self.window = window(at, chunk.end);
            while (at << 16) < chunk.end {
                let next = self.cached_member(source, at, &cache)?.next;
                if merged {
                    check_chunk_starts(original, &mut cursor, at, next)?;
                }
                at = next;
            }
        }
        Ok(())
    }

    /// Makes at least `n` bytes after the read position available unless
    /// the file ends first; returns how many there are.
    pub fn fill<S: ReadAt + ?Sized>(&mut self, source: &S, n: usize) -> Result<usize> {
        while self.text.len() - self.pos < n {
            if !self.load_next(source)? {
                break;
            }
        }
        Ok(self.text.len() - self.pos)
    }

    /// The available bytes after the read position.
    pub fn peek(&self) -> &[u8] {
        &self.text[self.pos..]
    }

    /// The text of a [`SeekReader::take`] range, valid until the next
    /// [`SeekReader::fill`] or [`SeekReader::seek`].
    pub fn slice(&self, range: std::ops::Range<usize>) -> &[u8] {
        &self.text[range]
    }

    /// Consumes `n` available bytes.
    pub fn advance(&mut self, n: usize) {
        assert!(
            n <= self.text.len() - self.pos,
            "advance exceeds available BGZF text"
        );
        if n == 0 {
            return;
        }
        self.pos += n;
        let last = self.pos - 1;
        while !(self.blocks[self.cur].text <= last
            && last < self.blocks[self.cur].text + self.blocks[self.cur].len)
        {
            self.cur += 1;
        }
        let b = self.blocks[self.cur];
        self.here = if self.pos == b.text + b.len {
            b.end_address << 16
        } else {
            b.address << 16 | (self.pos - b.text) as u64
        };
    }

    /// Consumes the next `n` bytes and returns their [`SeekReader::slice`]
    /// range; [`Error::Truncated`] when the file ends first.
    pub fn take<S: ReadAt + ?Sized>(
        &mut self,
        source: &S,
        n: usize,
        what: &'static str,
    ) -> Result<std::ops::Range<usize>> {
        let have = self.fill(source, n)?;
        if have < n {
            return Err(Error::Truncated {
                what,
                need: n as u64,
                have: have as u64,
            });
        }
        let start = self.pos;
        self.advance(n);
        Ok(start..start + n)
    }
}
