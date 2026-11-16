//! BGZF, the blocked gzip of bgzip, BAM and BCF: independent gzip members of
//! at most 64 KiB, each recording its compressed size, so any block
//! inflates on its own and blocks compress in parallel.
//!
//! - [`BlockHeader`], [`Inflater`], [`deflate_block`]: one block at a time.
//! - [`Writer`]: cuts text into blocks exactly where htslib's `bgzf_write`
//!   cuts them, compresses them through an [`omicsweave_core::exec`] queue
//!   and writes them in order; positions it reports are
//!   [`Logical`] offsets that [`Blocks::resolve`] turns into virtual
//!   offsets once the file is written.
//!   [`Writer::new_streaming`] omits that directory for non-indexed output.
//! - [`Reader`]: inflates consecutive blocks of a stream ahead of the
//!   reader through an [`omicsweave_core::exec`] queue.
//! - [`Scanner`]: a whole stream front to back with the [`Logical`]
//!   position of every byte consumed (and the block directory resolving
//!   them), for sequential record decoders and index builders.
//! - [`BinningIndex`]: the BAI / CSI index htslib builds from the same
//!   records.
//! - [`BinIndex`], [`SeekReader`], [`ChunkWalk`]: a BAI / CSI index loaded,
//!   and the records of a region read through it at virtual offsets, as
//!   htslib's region iterators read BAM and BCF.

mod bins;
mod index;
mod reader;
mod scan;
mod seek;
mod writer;

use std::cell::RefCell;

use omicsweave_core::{Error, Result};

pub use bins::{BinIndex, Chunk, ChunkWalk, IndexKind, locate_index};
pub use index::{BinningIndex, IndexFormat, try_csi_levels, try_reg2bin};
pub use reader::{Reader, is_bgzf, read_block};
pub use scan::{Scanner, inflate_queue};
pub use seek::{SeekReader, SeekStats};
pub use writer::{Blocks, Finished, Logical, Writer};

/// Upper bound of a block's compressed and decompressed size.
pub const MAX_BLOCK_LEN: usize = 0x10000;
/// Text bytes per full block, as htslib's `BGZF_BLOCK_SIZE`; leaves room
/// for incompressible data to fit [`MAX_BLOCK_LEN`].
pub const BLOCK_TEXT_LEN: usize = 0xff00;
/// The empty block that ends every BGZF file.
pub const EOF_BLOCK: [u8; 28] = [
    0x1f, 0x8b, 0x08, 0x04, 0x00, 0x00, 0x00, 0x00, 0x00, 0xff, 0x06, 0x00, 0x42, 0x43, 0x02, 0x00,
    0x1b, 0x00, 0x03, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
];
/// Header of a block whose only extra subfield is `BC`.
pub const HEADER_LEN: usize = 18;
/// CRC32 and ISIZE.
pub const TRAILER_LEN: usize = 8;
const FEXTRA: u8 = 4;
/// Compression level of htslib's default (`Z_DEFAULT_COMPRESSION`).
pub const DEFAULT_LEVEL: i32 = 6;

#[inline]
fn u16_at(buf: &[u8], at: usize) -> u16 {
    u16::from_le_bytes([buf[at], buf[at + 1]])
}

#[inline]
fn u32_at(buf: &[u8], at: usize) -> u32 {
    u32::from_le_bytes(buf[at..at + 4].try_into().unwrap())
}

fn require(buf: &[u8], need: usize, what: &'static str) -> Result<()> {
    if buf.len() < need {
        return Err(Error::Truncated {
            what,
            need: need as u64,
            have: buf.len() as u64,
        });
    }
    Ok(())
}

fn block_text_len(raw: &[u8], header: BlockHeader) -> Result<usize> {
    require(raw, header.block_len, "bgzf block")?;
    let size = u32_at(raw, header.block_len - 4) as usize;
    if size > MAX_BLOCK_LEN {
        return Err(Error::invalid(
            "bgzf block",
            "uncompressed size exceeds 64 KiB",
        ));
    }
    Ok(size)
}

/// Layout of one block from its gzip header.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct BlockHeader {
    /// Bytes of the whole block.
    pub block_len: usize,
    /// Bytes before the deflate payload.
    pub header_len: usize,
}

impl BlockHeader {
    /// Bytes to read before [`BlockHeader::parse`] can succeed on a block
    /// carrying only the `BC` subfield (the common case; others report a
    /// larger need through [`Error::Truncated`]).
    pub const MIN_LEN: usize = HEADER_LEN;

    pub fn parse(buf: &[u8]) -> Result<BlockHeader> {
        require(buf, 12, "bgzf header")?;
        if buf[0] != 0x1f || buf[1] != 0x8b || buf[2] != 8 || buf[3] & FEXTRA == 0 {
            return Err(Error::invalid("bgzf header", "not a BGZF block"));
        }
        if buf[3] != FEXTRA {
            return Err(Error::Unsupported("BGZF with optional gzip header fields"));
        }
        let xlen = u16_at(buf, 10) as usize;
        require(buf, 12 + xlen, "bgzf header")?;
        let mut at = 12;
        let end = 12 + xlen;
        let mut block_len = None;
        while at < end {
            if end - at < 4 {
                return Err(Error::invalid("bgzf header", "incomplete extra subfield"));
            }
            let len = u16_at(buf, at + 2) as usize;
            if len > end - at - 4 {
                return Err(Error::invalid("bgzf header", "extra subfield exceeds XLEN"));
            }
            if buf[at] == b'B' && buf[at + 1] == b'C' {
                if len != 2 || block_len.is_some() {
                    return Err(Error::invalid(
                        "bgzf header",
                        "invalid or duplicate BC subfield",
                    ));
                }
                block_len = Some(u16_at(buf, at + 4) as usize + 1);
            }
            at += 4 + len;
        }
        let block_len =
            block_len.ok_or_else(|| Error::invalid("bgzf header", "missing BC subfield"))?;
        if block_len < end + TRAILER_LEN {
            return Err(Error::invalid(
                "bgzf header",
                "block size below its header and trailer",
            ));
        }
        Ok(BlockHeader {
            block_len,
            header_len: end,
        })
    }
}

/// Reusable inflate state; create one per thread and keep it.
pub struct Inflater {
    inner: libdeflater::Decompressor,
}

impl Default for Inflater {
    fn default() -> Self {
        Self::new()
    }
}

impl Inflater {
    pub fn new() -> Inflater {
        Inflater {
            inner: libdeflater::Decompressor::new(),
        }
    }

    /// Inflates the block at the start of `raw` into `out` (at most
    /// [`MAX_BLOCK_LEN`] bytes are needed) and returns (block bytes
    /// consumed, text bytes written). The CRC and size in the trailer are
    /// verified.
    pub fn inflate_block(&mut self, raw: &[u8], out: &mut [u8]) -> Result<(usize, usize)> {
        let header = BlockHeader::parse(raw)?;
        self.inflate_parsed(raw, header, out)
    }

    fn inflate_parsed(
        &mut self,
        raw: &[u8],
        header: BlockHeader,
        out: &mut [u8],
    ) -> Result<(usize, usize)> {
        let size = block_text_len(raw, header)?;
        let block = &raw[..header.block_len];
        let crc = u32_at(block, header.block_len - 8);
        if size > out.len() {
            return Err(Error::Truncated {
                what: "bgzf output",
                need: size as u64,
                have: out.len() as u64,
            });
        }
        let payload = &block[header.header_len..header.block_len - TRAILER_LEN];
        let n = self
            .inner
            .deflate_decompress(payload, &mut out[..size])
            .map_err(|_| Error::Inflate("corrupt deflate data"))?;
        if n != size {
            return Err(Error::Inflate("decompressed size differs from the trailer"));
        }
        if libdeflater::crc32(&out[..size]) != crc {
            return Err(Error::Inflate("CRC mismatch"));
        }
        Ok((header.block_len, size))
    }

    /// Inflates `blocks` consecutive blocks from `raw` into `out` and
    /// returns the text bytes written. Bytes after the last block are
    /// ignored.
    pub fn inflate_blocks(&mut self, raw: &[u8], blocks: usize, out: &mut [u8]) -> Result<usize> {
        let (mut src, mut dst) = (0, 0);
        for _ in 0..blocks {
            let (used, written) = self.inflate_block(&raw[src..], &mut out[dst..])?;
            src += used;
            dst += written;
        }
        Ok(dst)
    }
}

/// Checks a libdeflate compression level (0 to 12).
pub fn check_level(level: i32) -> Result<libdeflater::CompressionLvl> {
    libdeflater::CompressionLvl::new(level)
        .map_err(|_| Error::invalid("bgzf", format!("compression level {level} outside 0..=12")))
}

/// Compresses `text` (at most [`BLOCK_TEXT_LEN`] bytes) into one complete
/// block appended to `out`.
pub fn deflate_block(
    compressor: &mut libdeflater::Compressor,
    text: &[u8],
    out: &mut Vec<u8>,
) -> Result<()> {
    if text.len() > BLOCK_TEXT_LEN {
        return Err(Error::invalid("bgzf", "block text exceeds BGZF_BLOCK_SIZE"));
    }
    let base = out.len();
    let room = compressor
        .deflate_compress_bound(text.len())
        .min(MAX_BLOCK_LEN - HEADER_LEN - TRAILER_LEN);
    out.resize(base + HEADER_LEN + room + TRAILER_LEN, 0);
    let body = &mut out[base + HEADER_LEN..base + HEADER_LEN + room];
    let n = compressor
        .deflate_compress(text, body)
        .map_err(|_| Error::invalid("bgzf", "block does not fit after compression"))?;
    let len = HEADER_LEN + n + TRAILER_LEN;
    let block = &mut out[base..base + len];
    block[..HEADER_LEN].copy_from_slice(&EOF_BLOCK[..HEADER_LEN]);
    block[16..18].copy_from_slice(&((len - 1) as u16).to_le_bytes());
    block[len - 8..len - 4].copy_from_slice(&libdeflater::crc32(text).to_le_bytes());
    block[len - 4..len].copy_from_slice(&(text.len() as u32).to_le_bytes());
    out.truncate(base + len);
    Ok(())
}

thread_local! {
    /// Jobs run on whichever worker picks them up and cannot carry a
    /// compressor with them; each thread keeps one per level instead of
    /// allocating libdeflate's state for every 64 KiB block.
    static COMPRESSORS: RefCell<Vec<(i32, libdeflater::Compressor)>> = const { RefCell::new(Vec::new()) };
    static INFLATER: RefCell<Option<Inflater>> = const { RefCell::new(None) };
}

/// [`deflate_block`] with this thread's compressor for `level` (already
/// validated by [`check_level`]).
pub fn deflate_block_at(level: i32, text: &[u8], out: &mut Vec<u8>) -> Result<()> {
    COMPRESSORS.with(|cell| {
        let mut all = cell.borrow_mut();
        let at = match all.iter().position(|(l, _)| *l == level) {
            Some(at) => at,
            None => {
                all.push((level, libdeflater::Compressor::new(check_level(level)?)));
                all.len() - 1
            }
        };
        deflate_block(&mut all[at].1, text, out)
    })
}

/// [`Inflater::inflate_block`] of one whole block with this thread's
/// inflater, returning its text.
pub fn inflate_block_owned(raw: &[u8]) -> Result<Vec<u8>> {
    INFLATER.with(|cell| {
        let mut slot = cell.borrow_mut();
        let inflater = slot.get_or_insert_with(Inflater::new);
        let header = BlockHeader::parse(raw)?;
        let mut text = vec![0u8; block_text_len(raw, header)?];
        let (used, n) = inflater.inflate_parsed(raw, header, &mut text)?;
        if used != raw.len() {
            return Err(Error::invalid("bgzf block", "bytes after the block"));
        }
        text.truncate(n);
        Ok(text)
    })
}
