//! UCSC 2bit: bases packed four per byte, with runs of N and soft-masked
//! (lowercase) runs stored as blocks.
//!
//! Metadata comes in three layers whose sizes are only known from the layer
//! before: the 16-byte header, the directory of names and offsets, and one
//! record header per sequence (its N and mask blocks). Parsers return
//! [`Parsed::Need`] with the exact byte count they require when the buffer is
//! too short, so callers can read in rounds and batch them across sequences.
//! Sequence data then needs one read per query.

use std::collections::HashMap;
use std::ops::Range;

use crate::error::{Error, Result, invalid, require, u32_at, u64_at};

pub const TWOBIT_MAGIC: u32 = 0x1A41_2743;
pub const HEADER_LEN: usize = 16;
/// Longest directory entry: name length byte, name, 64-bit offset.
const MAX_ENTRY_LEN: u64 = 1 + 255 + 8;

/// The outcome of parsing a variable-length structure from a prefix.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum Parsed<T> {
    Done(T),
    /// The structure needs at least this many bytes from its start; read
    /// them and parse again.
    Need(usize),
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct TwoBitHeader {
    /// 0: 32-bit directory offsets; 1: 64-bit offsets, for files over 4 GiB.
    pub version: u32,
    pub count: u32,
}

impl TwoBitHeader {
    pub fn parse(buf: &[u8]) -> Result<TwoBitHeader> {
        require(buf, HEADER_LEN, "2bit header")?;
        let magic = u32_at(buf, 0);
        if magic != TWOBIT_MAGIC {
            if magic.swap_bytes() == TWOBIT_MAGIC {
                return Err(Error::BigEndian);
            }
            return Err(Error::BadMagic {
                what: "2bit header",
                found: magic,
            });
        }
        let version = u32_at(buf, 4);
        if version > 1 {
            return Err(invalid("2bit header", format!("unknown version {version}")));
        }
        Ok(TwoBitHeader {
            version,
            count: u32_at(buf, 8),
        })
    }

    fn offset_len(&self) -> usize {
        if self.version == 1 { 8 } else { 4 }
    }

    /// File bytes that certainly hold the directory: it starts after the
    /// header and no entry is longer than 264 bytes. One read of this range
    /// lets [`Directory::parse`] succeed.
    pub fn directory_bound(&self, file_len: u64) -> Range<u64> {
        let start = HEADER_LEN as u64;
        start
            ..(start + self.count as u64 * MAX_ENTRY_LEN)
                .min(file_len)
                .max(start)
    }
}

/// Sequence names and their record offsets, in file order.
#[derive(Debug, Clone, Default, PartialEq, Eq)]
pub struct Directory {
    entries: Vec<(String, u64)>,
    by_name: HashMap<String, usize>,
}

impl Directory {
    /// Parses the directory from `buf`, which starts at file offset
    /// [`HEADER_LEN`].
    pub fn parse(header: &TwoBitHeader, buf: &[u8]) -> Result<Parsed<Directory>> {
        if header.version > 1 {
            return Err(invalid(
                "2bit header",
                format!("unknown version {}", header.version),
            ));
        }
        let width = header.offset_len();
        let capacity = (header.count as usize).min(buf.len() / (1 + width));
        let mut entries = Vec::with_capacity(capacity);
        let mut at = 0usize;
        for _ in 0..header.count {
            if at >= buf.len() {
                return Ok(Parsed::Need(at + 1));
            }
            let name_len = buf[at] as usize;
            let end = at + 1 + name_len + width;
            if end > buf.len() {
                return Ok(Parsed::Need(end));
            }
            let name = std::str::from_utf8(&buf[at + 1..at + 1 + name_len])
                .map_err(|e| invalid("2bit name", e.to_string()))?
                .to_owned();
            let offset_at = at + 1 + name_len;
            let offset = if width == 8 {
                u64_at(buf, offset_at)
            } else {
                u32_at(buf, offset_at) as u64
            };
            entries.push((name, offset));
            at = end;
        }
        let mut by_name = HashMap::with_capacity(entries.len());
        for (i, (name, _)) in entries.iter().enumerate() {
            if by_name.insert(name.clone(), i).is_some() {
                return Err(invalid(
                    "2bit directory",
                    format!("duplicate name {name:?}"),
                ));
            }
        }
        Ok(Parsed::Done(Directory { entries, by_name }))
    }

    /// Record offset of `name`.
    pub fn offset(&self, name: &str) -> Option<u64> {
        self.by_name.get(name).map(|&i| self.entries[i].1)
    }

    pub(crate) fn index(&self, name: &str) -> Option<usize> {
        self.by_name.get(name).copied()
    }

    /// (name, record offset) in file order.
    pub fn entries(&self) -> &[(String, u64)] {
        &self.entries
    }

    pub fn len(&self) -> usize {
        self.entries.len()
    }

    pub fn is_empty(&self) -> bool {
        self.entries.is_empty()
    }
}

/// A run of bases `[start, start + len)`.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Block {
    pub start: u32,
    pub len: u32,
}

impl Block {
    fn end(&self) -> u64 {
        self.start as u64 + self.len as u64
    }
}

/// One sequence's record header: length, N blocks, mask blocks, and where
/// its packed bases start.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Record {
    len: u32,
    n_blocks: BlockList,
    mask_blocks: BlockList,
    packed_offset: u64,
}

/// Retain the wire metadata for inspection, with a compact query view only
/// when empty or adjacent blocks would cause redundant work per query.
#[derive(Debug, Clone, PartialEq, Eq)]
struct BlockList {
    raw: Vec<Block>,
    runs: Option<Box<[Block]>>,
}

impl BlockList {
    fn query(&self) -> &[Block] {
        self.runs.as_deref().unwrap_or(&self.raw)
    }
}

/// Bases of each packed byte, most significant bits first.
const fn byte_table() -> [[u8; 4]; 256] {
    let bases = *b"TCAG";
    let mut table = [[0u8; 4]; 256];
    let mut b = 0;
    while b < 256 {
        table[b] = [
            bases[b >> 6],
            bases[(b >> 4) & 3],
            bases[(b >> 2) & 3],
            bases[b & 3],
        ];
        b += 1;
    }
    table
}
static BYTE_BASES: [[u8; 4]; 256] = byte_table();

fn blocks(buf: &[u8], count: usize, at: usize, len: u32, what: &'static str) -> Result<BlockList> {
    let blocks: Vec<Block> = (0..count)
        .map(|i| Block {
            start: u32_at(buf, at + 4 * i),
            len: u32_at(buf, at + 4 * (count + i)),
        })
        .collect();
    let mut needs_compaction = false;
    for (i, b) in blocks.iter().enumerate() {
        if b.end() > len as u64 {
            return Err(invalid(what, format!("block {i} ends past the sequence")));
        }
        if i > 0 {
            let previous_end = blocks[i - 1].end();
            if (b.start as u64) < previous_end {
                return Err(invalid(
                    what,
                    format!("block {i} is unsorted or overlaps the previous one"),
                ));
            }
            needs_compaction |= b.start as u64 == previous_end;
        }
        needs_compaction |= b.len == 0;
    }
    let compact = if !needs_compaction {
        None
    } else {
        let mut compact: Vec<Block> = Vec::new();
        for b in blocks.iter().filter(|b| b.len > 0) {
            if let Some(last) = compact
                .last_mut()
                .filter(|last| last.end() == b.start as u64)
            {
                // Both endpoints were checked against the u32 sequence length.
                last.len = (b.end() - last.start as u64) as u32;
            } else {
                compact.push(*b);
            }
        }
        Some(compact.into_boxed_slice())
    };
    Ok(BlockList {
        raw: blocks,
        runs: compact,
    })
}

impl Record {
    /// Bytes to read first from a record's offset: the length and the N
    /// block count.
    pub const FIRST_READ: usize = 8;

    /// Parses the record header in `buf`, read from file offset `offset`.
    /// Needs at most three rounds: [`Record::FIRST_READ`] bytes, then up to
    /// the mask block count, then the whole header.
    pub fn parse(buf: &[u8], offset: u64) -> Result<Parsed<Record>> {
        if buf.len() < 8 {
            return Ok(Parsed::Need(8));
        }
        let len = u32_at(buf, 0);
        let n_count = u32_at(buf, 4) as u64;
        if n_count > len as u64 {
            return Err(invalid(
                "2bit record",
                format!("{n_count} N blocks in {len} bases"),
            ));
        }
        let mask_at = usize::try_from(8 + 8 * n_count)
            .map_err(|_| invalid("2bit record", "header length exceeds usize"))?;
        let mask_need = mask_at
            .checked_add(4)
            .ok_or_else(|| invalid("2bit record", "header length exceeds usize"))?;
        if buf.len() < mask_need {
            return Ok(Parsed::Need(mask_need));
        }
        let mask_count = u32_at(buf, mask_at) as u64;
        if mask_count > len as u64 {
            return Err(invalid(
                "2bit record",
                format!("{mask_count} mask blocks in {len} bases"),
            ));
        }
        let header_len = usize::try_from(mask_at as u64 + 8 + 8 * mask_count)
            .map_err(|_| invalid("2bit record", "header length exceeds usize"))?;
        if buf.len() < header_len {
            return Ok(Parsed::Need(header_len));
        }
        let packed_offset = offset
            .checked_add(header_len as u64)
            .filter(|&at| at.checked_add((len as u64).div_ceil(4)).is_some())
            .ok_or_else(|| invalid("2bit record", "file range overflows u64"))?;
        Ok(Parsed::Done(Record {
            len,
            n_blocks: blocks(buf, n_count as usize, 8, len, "2bit N blocks")?,
            mask_blocks: blocks(
                buf,
                mask_count as usize,
                mask_at + 4,
                len,
                "2bit mask blocks",
            )?,
            packed_offset,
        }))
    }

    pub fn len(&self) -> u32 {
        self.len
    }

    pub fn is_empty(&self) -> bool {
        self.len == 0
    }

    pub fn n_blocks(&self) -> &[Block] {
        &self.n_blocks.raw
    }

    pub fn mask_blocks(&self) -> &[Block] {
        &self.mask_blocks.raw
    }

    pub(crate) fn check_file_len(&self, file_len: u64) -> Result<()> {
        let end = self.packed_offset + (self.len as u64).div_ceil(4);
        if end > file_len {
            return Err(Error::Truncated {
                what: "2bit packed bases",
                need: end,
                have: file_len,
            });
        }
        Ok(())
    }

    fn all_n(&self, start: u32, end: u32) -> bool {
        let blocks = self.n_blocks.query();
        let first = blocks.partition_point(|b| b.end() <= start as u64);
        // Query blocks are nonempty, disjoint and separated by a gap, so
        // complete coverage must come from exactly one run.
        blocks
            .get(first)
            .is_some_and(|b| b.start <= start && b.end() >= end as u64)
    }

    /// File range of the packed bytes holding bases `[start, end)`.
    /// Empty and entirely N-covered ranges need no packed bytes.
    pub fn packed_range(&self, start: u32, end: u32) -> Result<Range<u64>> {
        if start > end {
            return Err(Error::InvalidRegion { start, end });
        }
        if end > self.len {
            return Err(Error::RegionOutOfRange {
                start,
                end,
                len: self.len,
            });
        }
        if start == end || self.all_n(start, end) {
            return Ok(self.packed_offset..self.packed_offset);
        }
        Ok(self.packed_offset + start as u64 / 4..self.packed_offset + (end as u64).div_ceil(4))
    }

    /// Writes bases `[start, start + out.len())` into `out` from `packed`,
    /// the bytes of [`Record::packed_range`]. Bases are uppercase; N blocks
    /// read `N`; with `mask`, soft-masked bases are lowercase (as
    /// `twoBitToFa`; without it, as `twoBitToFa -noMask`).
    pub fn decode(&self, start: u32, packed: &[u8], out: &mut [u8], mask: bool) -> Result<()> {
        let end = u32::try_from(start as u64 + out.len() as u64).map_err(|_| {
            Error::RegionOutOfRange {
                start,
                end: u32::MAX,
                len: self.len,
            }
        })?;
        let range = self.packed_range(start, end)?;
        let need = (range.end - range.start) as usize;
        if packed.len() < need {
            return Err(Error::Truncated {
                what: "2bit packed bases",
                need: need as u64,
                have: packed.len() as u64,
            });
        }
        if out.is_empty() {
            return Ok(());
        }
        if range.is_empty() {
            out.fill(b'N');
            return Ok(());
        }

        let skip = (start % 4) as usize;
        let mut dst = 0;
        let mut byte = 0;
        if skip != 0 {
            let first = &BYTE_BASES[packed[0] as usize][skip..];
            let take = first.len().min(out.len());
            out[..take].copy_from_slice(&first[..take]);
            dst = take;
            byte = 1;
        }
        let whole = (out.len() - dst) / 4;
        let mut i = 0;
        while i + 4 <= whole {
            for lane in 0..4 {
                // SAFETY: packed_range and the length check guarantee every
                // input byte. whole counts complete groups of four output
                // bases after dst, and i+lane < whole. The lookup index is a
                // u8, so it is < 256. Reads/writes use alignment-one arrays;
                // exclusive out cannot alias packed or the static lookup.
                unsafe {
                    let b = *packed.get_unchecked(byte + i + lane);
                    out.as_mut_ptr()
                        .add(dst + 4 * (i + lane))
                        .cast::<[u8; 4]>()
                        .write(*BYTE_BASES.get_unchecked(b as usize));
                }
            }
            i += 4;
        }
        while i < whole {
            // SAFETY: The same validated input/output and lookup bounds hold
            // for each remaining complete group, without crossing the tail.
            unsafe {
                let b = *packed.get_unchecked(byte + i);
                out.as_mut_ptr()
                    .add(dst + 4 * i)
                    .cast::<[u8; 4]>()
                    .write(*BYTE_BASES.get_unchecked(b as usize));
            }
            i += 1;
        }
        dst += 4 * whole;
        byte += whole;
        if dst < out.len() {
            let rest = out.len() - dst;
            out[dst..].copy_from_slice(&BYTE_BASES[packed[byte] as usize][..rest]);
        }

        if mask {
            for (lo, hi) in overlaps(self.mask_blocks.query(), start as u64, end as u64) {
                out[lo..hi].make_ascii_lowercase();
            }
        }
        for (lo, hi) in overlaps(self.n_blocks.query(), start as u64, end as u64) {
            out[lo..hi].fill(b'N');
        }
        Ok(())
    }
}

/// Parts of sorted, disjoint `blocks` within `[start, end)`, as output
/// positions relative to `start`.
fn overlaps(blocks: &[Block], start: u64, end: u64) -> impl Iterator<Item = (usize, usize)> + '_ {
    let first = blocks.partition_point(|b| b.end() <= start);
    blocks[first..]
        .iter()
        .take_while(move |b| (b.start as u64) < end)
        .map(move |b| {
            ((b.start as u64).max(start) - start) as usize..(b.end().min(end) - start) as usize
        })
        .map(|r| (r.start, r.end))
}
