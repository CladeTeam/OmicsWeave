//! The bgzip index (`.gzi`): where each BGZF block starts in the compressed
//! file and in the uncompressed text. It turns a text range (from
//! [`crate::FaiEntry::text_range`]) into the blocks to read.

use std::ops::Range;
use std::path::{Path, PathBuf};

use omicsweave_core::{FileSource, ReadAt, sidecar};

use crate::bgzf::{MAX_BLOCK_LEN, Storage, detect, indexable_header};
use crate::error::{Error, Result, invalid, require, u32_at, u64_at};

/// Appended to a bgzip file's path to name its `.gzi`.
pub const GZI_SUFFIX: &str = ".gzi";

/// The `.gzi` path of a bgzip file, as bgzip and samtools name it
/// (`a.fa.gz` -> `a.fa.gz.gzi`).
pub fn gzi_path(bgzf: impl AsRef<Path>) -> PathBuf {
    let mut name = bgzf.as_ref().as_os_str().to_owned();
    name.push(GZI_SUFFIX);
    PathBuf::from(name)
}

/// Bytes of a block header as index building reads it.
const HEADER_LEN: usize = 18;
/// Bytes of the `ISIZE` field that ends every block.
const ISIZE_LEN: usize = 4;

/// One BGZF block to read and inflate.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct BlockRef {
    /// Position of the block in the index; stable, usable as a cache key.
    pub index: usize,
    /// Compressed bytes to read. For the last block the index knows of, the
    /// end is an upper bound and may run into the next block; indexes built
    /// by htslib can also omit empty BGZF members between indexed blocks.
    /// The block's own header gives its exact length (at most 64 KiB).
    pub compressed: Range<u64>,
    /// Text offset of the block's first decompressed byte.
    pub text_start: u64,
    /// Decompressed length, when the index knows where the next block starts.
    pub text_len: Option<u64>,
}

/// Consecutive blocks covering a text range, for a single read.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Span {
    pub compressed: Range<u64>,
    /// Text offset of the first block; the requested range starts
    /// `text.start - text_start` bytes into the decompressed span.
    pub text_start: u64,
    /// Number of indexed entries. An index built by htslib may omit empty
    /// BGZF members, so this need not equal the physical member count
    /// expected by [`crate::Inflater::inflate_blocks`]. Use
    /// [`crate::FastaReader`] or follow physical headers until the requested
    /// text range is covered when the span contains skipped empty members.
    pub blocks: usize,
    /// Bytes the decompressed blocks can take.
    pub max_text_len: usize,
}

/// A parsed `.gzi`. Index metadata is small and kept by the value.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Gzi {
    /// (compressed, text) start of every block, beginning with (0, 0), which
    /// the file leaves implicit.
    starts: Vec<(u64, u64)>,
}

impl Gzi {
    /// Parses `.gzi` bytes: a little-endian u64 count followed by that many
    /// (compressed, uncompressed) u64 pairs.
    pub fn parse(buf: &[u8]) -> Result<Gzi> {
        require(buf, 8, "gzi")?;
        let n = u64_at(buf, 0);
        if n > (buf.len() as u64 - 8) / 16 {
            return Err(Error::Truncated {
                what: "gzi",
                need: n.saturating_mul(16).saturating_add(8),
                have: buf.len() as u64,
            });
        }
        let mut starts = Vec::with_capacity(n as usize + 1);
        starts.push((0, 0));
        for i in 0..n as usize {
            let at = 8 + 16 * i;
            let point = (u64_at(buf, at), u64_at(buf, at + 8));
            let last = *starts.last().unwrap();
            if point.0 <= last.0 || point.1 < last.1 {
                return Err(invalid("gzi", format!("entry {i} is out of order")));
            }
            if point.1 - last.1 > MAX_BLOCK_LEN as u64 {
                return Err(invalid(
                    "gzi",
                    format!("entry {i} exceeds the 64 KiB BGZF block limit"),
                ));
            }
            starts.push(point);
        }
        Ok(Gzi { starts })
    }

    /// Index from block starts, the first being (0, 0).
    pub(crate) fn from_starts(starts: Vec<(u64, u64)>) -> Gzi {
        debug_assert_eq!(starts.first(), Some(&(0, 0)));
        Gzi { starts }
    }

    /// Index from the (compressed, text) start of every non-empty block, as
    /// htslib records them while reading: empty blocks get no entry, and
    /// the first entry becomes the implicit (0, 0) even when empty blocks
    /// precede it.
    pub(crate) fn from_nonempty_blocks(mut starts: Vec<(u64, u64)>) -> Gzi {
        match starts.first_mut() {
            Some(first) => *first = (0, 0),
            None => starts.push((0, 0)),
        }
        Gzi { starts }
    }

    /// Builds the index of a BGZF source from its block headers and `ISIZE`
    /// trailers, without inflating; the result is the `.gzi` that
    /// `bgzip -r` and `samtools faidx` write for it. Reads one small window
    /// per block. Truncated or malformed blocks are errors, a missing
    /// end-of-file block is not (htslib only warns), and plain gzip is
    /// [`Error::Unsupported`].
    pub fn build<S: ReadAt + ?Sized>(source: &S) -> Result<Gzi> {
        if detect(source)? != Storage::Bgzf {
            return Err(invalid("bgzf", "not a BGZF file"));
        }
        let len = source.len();
        // The trailer of one block and the header of the next, read together.
        let mut window = [0u8; ISIZE_LEN + HEADER_LEN];
        source.read_exact_at(&mut window[ISIZE_LEN..], 0)?;
        let (mut pos, mut text) = (0u64, 0u64);
        let mut starts = Vec::new();
        while pos < len {
            let block_len = indexable_header(&window[ISIZE_LEN..])?.block_len as u64;
            let end = pos + block_len;
            if end > len {
                return Err(Error::Truncated {
                    what: "bgzf block",
                    need: block_len,
                    have: len - pos,
                });
            }
            let next = (len - end).min(HEADER_LEN as u64) as usize;
            if next > 0 && next < HEADER_LEN {
                return Err(Error::Truncated {
                    what: "bgzf header",
                    need: HEADER_LEN as u64,
                    have: next as u64,
                });
            }
            source.read_exact_at(&mut window[..ISIZE_LEN + next], end - ISIZE_LEN as u64)?;
            let size = u32_at(&window, 0) as u64;
            if size > MAX_BLOCK_LEN as u64 {
                return Err(invalid("bgzf block", "uncompressed size exceeds 64 KiB"));
            }
            if size > 0 {
                starts.push((pos, text));
                text += size;
            }
            pos = end;
        }
        Ok(Gzi::from_nonempty_blocks(starts))
    }

    /// Writes [`gzi_path`]`(bgzf)` atomically.
    pub fn save(&self, bgzf: impl AsRef<Path>) -> Result<()> {
        sidecar::write_file(&gzi_path(bgzf), &self.to_bytes())
    }

    /// Reads and parses [`gzi_path`]`(bgzf)`.
    pub fn load(bgzf: impl AsRef<Path>) -> Result<Gzi> {
        Gzi::parse(&std::fs::read(gzi_path(bgzf))?)
    }

    /// Builds the index of the bgzip file at `bgzf` and saves it beside it.
    pub fn build_file(bgzf: impl AsRef<Path>) -> Result<Gzi> {
        let bgzf = bgzf.as_ref();
        let gzi = Gzi::build(&FileSource::open(bgzf)?)?;
        gzi.save(bgzf)?;
        Ok(gzi)
    }

    /// The `.gzi` bytes, as bgzip writes them.
    pub fn to_bytes(&self) -> Vec<u8> {
        let mut out = Vec::with_capacity(8 + 16 * self.starts.len());
        out.extend_from_slice(&(self.starts.len() as u64 - 1).to_le_bytes());
        for &(c, u) in &self.starts[1..] {
            out.extend_from_slice(&c.to_le_bytes());
            out.extend_from_slice(&u.to_le_bytes());
        }
        out
    }

    /// (compressed, text) start of every indexed block.
    pub fn starts(&self) -> &[(u64, u64)] {
        &self.starts
    }

    /// The block holding text byte `offset`.
    fn block_of(&self, offset: u64) -> usize {
        self.starts.partition_point(|&(_, u)| u <= offset) - 1
    }

    pub(crate) fn block(&self, index: usize, file_len: u64) -> BlockRef {
        let (c, u) = self.starts[index];
        let next = self.starts.get(index + 1);
        BlockRef {
            index,
            compressed: c..next.map_or(c.saturating_add(MAX_BLOCK_LEN as u64).min(file_len), |n| {
                n.0
            }),
            text_start: u,
            text_len: next.map(|n| n.1 - u),
        }
    }

    /// Blocks holding the text bytes `text`, in file order. `file_len` is the
    /// compressed file's length, which bounds the last block.
    pub(crate) fn indices(&self, text: Range<u64>) -> Range<usize> {
        if text.start >= text.end {
            0..0
        } else {
            self.block_of(text.start)..self.block_of(text.end - 1) + 1
        }
    }

    pub fn blocks(
        &self,
        text: Range<u64>,
        file_len: u64,
    ) -> impl ExactSizeIterator<Item = BlockRef> + DoubleEndedIterator + std::iter::FusedIterator + '_
    {
        let range = self.indices(text);
        range.map(move |i| self.block(i, file_len))
    }

    /// The blocks of [`Gzi::blocks`] as one contiguous read; `None` for an
    /// empty range. The returned [`Span::blocks`] counts indexed entries;
    /// empty physical members omitted by the index are not counted.
    pub fn span(&self, text: Range<u64>, file_len: u64) -> Option<Span> {
        let mut blocks = self.blocks(text, file_len);
        let first = blocks.next()?;
        let last = blocks.last().unwrap_or_else(|| first.clone());
        let text_end = last
            .text_len
            .map_or(last.text_start.saturating_add(MAX_BLOCK_LEN as u64), |n| {
                last.text_start.saturating_add(n)
            });
        Some(Span {
            compressed: first.compressed.start..last.compressed.end,
            text_start: first.text_start,
            blocks: last.index - first.index + 1,
            max_text_len: usize::try_from(text_end - first.text_start).ok()?,
        })
    }
}
