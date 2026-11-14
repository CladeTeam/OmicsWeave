//! Building the indexes of an existing FASTA file: the `.fai` of plain or
//! bgzip FASTA and, for bgzip, the `.gzi`, in one bounded streaming pass.

use std::path::Path;

use omicsweave_core::{FileSource, ReadAt};

use crate::bgzf::{Inflater, MAX_BLOCK_LEN, Storage, detect, indexable_header};
use crate::error::{Error, Result};
use crate::fai::Fai;
use crate::faibuild::FaiBuilder;
use crate::gzi::Gzi;

/// Bytes read from the source at a time.
const READ_LEN: usize = 1 << 20;

/// Builds the `.fai` and, for BGZF sources, the `.gzi`, in one pass over the
/// source. The storage is detected from the first bytes; plain gzip is
/// [`Error::Unsupported`], as htslib cannot index it either. BGZF blocks are
/// inflated once each and checked (CRC, size); memory stays bounded by one
/// read buffer and one block. The results are byte-identical to the files
/// `samtools faidx` writes, and the `.gzi` equals [`Gzi::build`].
pub fn build_index<S: ReadAt + ?Sized>(source: &S) -> Result<(Fai, Option<Gzi>)> {
    let mut fai = FaiBuilder::new();
    let mut buf = Buffered::new(source);
    match detect(source)? {
        Storage::Plain => {
            let mut pos = 0;
            while pos < source.len() {
                let n = (source.len() - pos).min(READ_LEN as u64) as usize;
                fai.push(buf.get(pos, n, "FASTA text")?)?;
                pos += n as u64;
            }
            Ok((fai.finish()?, None))
        }
        Storage::Bgzf => {
            let mut inflater = Inflater::new();
            let mut block = vec![0u8; MAX_BLOCK_LEN];
            let (mut pos, mut text) = (0u64, 0u64);
            let mut starts = Vec::new();
            while pos < source.len() {
                let header = indexable_header(buf.get(pos, BLOCK_HEADER_LEN, "bgzf header")?)?;
                let raw = buf.get(pos, header.block_len, "bgzf block")?;
                let (_, n) = inflater.inflate_block(raw, &mut block)?;
                if n > 0 {
                    starts.push((pos, text));
                    text += n as u64;
                    fai.push(&block[..n])?;
                }
                pos += header.block_len as u64;
            }
            Ok((fai.finish()?, Some(Gzi::from_nonempty_blocks(starts))))
        }
    }
}

/// [`build_index`] on the file at `fasta`, saving the `.fai` (and the
/// `.gzi`) beside it, as `samtools faidx` names them.
pub fn build_index_files(fasta: impl AsRef<Path>) -> Result<(Fai, Option<Gzi>)> {
    let fasta = fasta.as_ref();
    let (fai, gzi) = build_index(&FileSource::open(fasta)?)?;
    if let Some(gzi) = &gzi {
        gzi.save(fasta)?;
    }
    fai.save(fasta)?;
    Ok((fai, gzi))
}

const BLOCK_HEADER_LEN: usize = 18;

/// Forward reads of a source through one reusable window.
struct Buffered<'a, S: ?Sized> {
    source: &'a S,
    start: u64,
    buf: Vec<u8>,
}

impl<'a, S: ReadAt + ?Sized> Buffered<'a, S> {
    fn new(source: &'a S) -> Buffered<'a, S> {
        Buffered {
            source,
            start: 0,
            buf: Vec::new(),
        }
    }

    /// The `n` bytes at `pos` (`n` at most [`READ_LEN`]); refills the
    /// window from `pos` when they are not all in it. Bytes past the end
    /// of the source are [`Error::Truncated`] `what`.
    fn get(&mut self, pos: u64, n: usize, what: &'static str) -> Result<&[u8]> {
        debug_assert!(n <= READ_LEN);
        let len = self.source.len();
        if pos.checked_add(n as u64).is_none_or(|end| end > len) {
            return Err(Error::Truncated {
                what,
                need: n as u64,
                have: len.saturating_sub(pos),
            });
        }
        let end = self.start + self.buf.len() as u64;
        if pos < self.start || pos + n as u64 > end {
            let fill = (len - pos).min(READ_LEN as u64) as usize;
            self.buf.resize(fill, 0);
            self.source.read_exact_at(&mut self.buf, pos)?;
            self.start = pos;
        }
        let at = (pos - self.start) as usize;
        Ok(&self.buf[at..at + n])
    }
}
