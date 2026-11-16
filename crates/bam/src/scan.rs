//! A whole BAM read front to back: records framed over a BGZF
//! [`omicsweave_bgzf::Scanner`], each with the position htslib's
//! `bgzf_tell` reports after it. Index building uses it.

use std::fs::File;
use std::io::{BufReader, Read};
use std::ops::Range;
use std::path::{Path, PathBuf};

use omicsweave_bgzf::{BinningIndex, Blocks, IndexFormat, IndexKind, Logical, inflate_queue};
use omicsweave_core::exec::Ordered;
use omicsweave_core::{Error, Result, sidecar};

use crate::header::{BamHeader, Take};
use crate::record::{Record, flag};

/// The records of a BAM stream in file order.
pub struct Scanner<R: Read> {
    stream: omicsweave_bgzf::Scanner<R>,
}

impl<R: Read> Scanner<R> {
    /// Reads `input` from its start; blocks inflate on `exec`.
    pub fn new(input: R, exec: Box<dyn Ordered<Vec<u8>>>) -> Scanner<R> {
        Scanner {
            stream: omicsweave_bgzf::Scanner::new(input, exec),
        }
    }

    /// Reads records without retaining the file-wide block directory.
    ///
    /// Use this for one-pass consumers that need records and input progress,
    /// but do not build an index. The logical positions returned by `tell`
    /// are unchanged; this scanner has no directory to resolve them through.
    pub fn new_streaming(input: R, exec: Box<dyn Ordered<Vec<u8>>>) -> Scanner<R> {
        Scanner {
            stream: omicsweave_bgzf::Scanner::new_streaming(input, exec),
        }
    }

    /// `bgzf_tell` after what has been read, with the block's number in
    /// place of its address.
    pub fn tell(&self) -> Logical {
        self.stream.tell()
    }

    /// Compressed bytes taken from the input so far.
    pub fn compressed_read(&self) -> u64 {
        self.stream.compressed_read()
    }

    /// Where the blocks read so far start; complete once the stream has
    /// ended. Panics for a streaming scanner.
    pub fn blocks(&self) -> &Blocks {
        self.stream.blocks()
    }

    /// The retained block directory, or `None` in streaming mode.
    pub fn try_blocks(&self) -> Option<&Blocks> {
        self.stream.try_blocks()
    }

    /// Reads the header.
    pub fn header(&mut self) -> Result<BamHeader> {
        BamHeader::read(self)
    }

    /// The next record after its `block_size`; `None` at a clean end.
    pub fn next_record(&mut self) -> Result<Option<Range<usize>>> {
        let have = self.stream.fill(4)?;
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
        // SAFETY: fill returned at least four unread initialized bytes.
        // peek borrows that same buffer and the unaligned load reads
        // exactly the length prefix without retaining its pointer.
        let size =
            unsafe { i32::from_le(self.stream.peek().as_ptr().cast::<i32>().read_unaligned()) };
        if size < 32 {
            return Err(Error::invalid(
                "BAM record",
                format!("block_size {size} below the 32 bytes of fixed fields"),
            ));
        }
        let range = self.stream.take_range(4 + size as usize, "BAM record")?;
        Ok(Some(range.start + 4..range.end))
    }

    /// The text of [`Scanner::next_record`]'s range.
    pub fn slice(&self, range: Range<usize>) -> &[u8] {
        self.stream.slice(range)
    }
}

impl<R: Read> Take for Scanner<R> {
    fn take(&mut self, n: usize, what: &'static str) -> Result<&[u8]> {
        self.stream.take(n, what)
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
    use omicsweave_bgzf::EOF_BLOCK;
    use omicsweave_core::exec::Inline;

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
                let mut scan = if streaming {
                    Scanner::new_streaming(&raw[..], exec)
                } else {
                    Scanner::new(&raw[..], exec)
                };
                scan.header().unwrap();
                let first = scan.next_record().unwrap_err();
                assert_eq!(scan.next_record().unwrap_err(), first);
            }
        }
    }
}
