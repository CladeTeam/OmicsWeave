//! BGZF, the blocked gzip of bgzip: independent gzip members of at most
//! 64 KiB, each recording its compressed size, so any block inflates on its
//! own. Blocks themselves are `omicsweave-bgzf`'s; this module adds what
//! FASTA needs on top: detecting the storage of a source, the header rule
//! of index building, and a writer that returns the `.gzi`.

use std::io::{self, Write};

use omicsweave_core::ReadAt;
use omicsweave_core::exec::{Inline, Ordered};

pub use omicsweave_bgzf::{BLOCK_TEXT_LEN, BlockHeader, EOF_BLOCK, Inflater, MAX_BLOCK_LEN};

use crate::error::{Error, Result, invalid, require, u16_at};
use crate::gzi::Gzi;

/// Header of a block whose only extra subfield is `BC`.
const HEADER_LEN: usize = 18;
const FEXTRA: u8 = 4;

/// Header of a block as htslib's reader accepts it: `XLEN` 6 holding only
/// the `BC` subfield. Index building applies this rule so that it succeeds
/// and fails on the same files as `samtools faidx` and `bgzip -r`; htslib
/// takes any other gzip member, including a BGZF block with more
/// subfields, for plain gzip, which it cannot index.
pub(crate) fn indexable_header(buf: &[u8]) -> Result<BlockHeader> {
    require(buf, HEADER_LEN, "bgzf header")?;
    if buf[0] != 0x1f || buf[1] != 0x8b || buf[2] != 8 {
        return Err(invalid("bgzf header", "not a gzip member"));
    }
    if buf[3] & FEXTRA == 0 || u16_at(buf, 10) != 6 || &buf[12..16] != b"BC\x02\x00" {
        return Err(Error::Unsupported(
            "gzip members that are not BGZF blocks cannot be indexed; recompress with bgzip",
        ));
    }
    BlockHeader::parse(buf)
}

/// How a FASTA source is stored.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) enum Storage {
    Plain,
    Bgzf,
}

/// Detects the storage from the first bytes as htslib does: a source of at
/// least one block header that starts with the gzip magic is compressed,
/// and must then be BGZF ([`Error::Unsupported`] for plain gzip); anything
/// else is plain text.
pub(crate) fn detect<S: ReadAt + ?Sized>(source: &S) -> Result<Storage> {
    if source.len() < HEADER_LEN as u64 {
        return Ok(Storage::Plain);
    }
    let mut magic = [0u8; HEADER_LEN];
    source.read_exact_at(&mut magic, 0)?;
    if magic[0] != 0x1f || magic[1] != 0x8b {
        return Ok(Storage::Plain);
    }
    if magic[3] & FEXTRA == 0 || &magic[12..16] != b"BC\x02\x00" {
        return Err(Error::Unsupported(
            "gzip that is not BGZF cannot be indexed; recompress with bgzip",
        ));
    }
    Ok(Storage::Bgzf)
}

/// Writes BGZF to `inner` and records the `.gzi` of what it writes. Blocks
/// hold [`BLOCK_TEXT_LEN`] text bytes, as bgzip writes them, so block
/// boundaries match bgzip's.
pub struct BgzfWriter<W: Write> {
    inner: omicsweave_bgzf::Writer<W>,
}

impl<W: Write> BgzfWriter<W> {
    /// A writer at libdeflate's default compression level.
    pub fn new(inner: W) -> BgzfWriter<W> {
        Self::with_level(inner, omicsweave_bgzf::DEFAULT_LEVEL).expect("default level is valid")
    }

    /// A writer at compression `level` (0 to 12, libdeflate's scale).
    pub fn with_level(inner: W, level: i32) -> Result<BgzfWriter<W>> {
        Self::with_exec(inner, level, Box::new(Inline::new()))
    }

    /// A writer compressing blocks on `exec`; the output does not depend on
    /// it.
    pub fn with_exec(
        inner: W,
        level: i32,
        exec: Box<dyn Ordered<Vec<u8>>>,
    ) -> Result<BgzfWriter<W>> {
        Ok(BgzfWriter {
            inner: omicsweave_bgzf::Writer::new(inner, level, exec)?,
        })
    }

    /// Writes the last block and the end-of-file marker, and returns the
    /// sink with the index of the written blocks.
    pub fn finish(self) -> io::Result<(W, Gzi)> {
        let done = self.inner.finish().map_err(io::Error::other)?;
        let mut starts = done.blocks.starts;
        if starts.is_empty() {
            starts.push((0, 0));
        }
        Ok((done.out, Gzi::from_starts(starts)))
    }
}

impl<W: Write> Write for BgzfWriter<W> {
    fn write(&mut self, data: &[u8]) -> io::Result<usize> {
        self.inner.write(data).map_err(io::Error::other)?;
        Ok(data.len())
    }

    /// Ends the current block early (like htslib's `bgzf_flush`), writes
    /// it and flushes the sink.
    fn flush(&mut self) -> io::Result<()> {
        self.inner.sync().map_err(io::Error::other)
    }
}
