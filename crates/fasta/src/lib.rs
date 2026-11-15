//! Random access to sequences in FASTA (`.fa` + `.fai`), bgzip FASTA
//! (`.fa.gz` + `.fai` + `.gzi`) and UCSC 2bit (`.2bit`) files.
//!
//! The query steps perform no reads and never hold sequence data. Index metadata
//! ([`Fai`], [`Gzi`], the 2bit [`Directory`] and [`Record`]s) is small,
//! parsed once and kept by its value. Every query first computes the file
//! bytes it needs, which the caller reads (or a loader prefetches and
//! caches), then decodes from the caller's bytes into the caller's output:
//!
//! | File | Describe | Decode |
//! |---|---|---|
//! | `.fa` | [`FaiEntry::text_range`] | [`FaiEntry::decode`] |
//! | `.fa.gz` | [`FaiEntry::text_range`], then [`Gzi::blocks`] or [`Gzi::span`] | [`Inflater::inflate_block`]s, then [`FaiEntry::decode`] |
//! | `.2bit` | [`Record::packed_range`] | [`Record::decode`] |
//!
//! [`Gzi::blocks`] yields one [`BlockRef`] per BGZF block with a stable
//! index, so a loader can cache inflated blocks and share them between
//! overlapping queries. 2bit metadata is parsed in rounds: [`Parsed::Need`]
//! reports the bytes still required.
//!
//! [`FastaReader`] (`.fa`, `.fa.gz`) and [`TwoBitReader`] (`.2bit`) compose
//! these steps over a [`omicsweave_core::ReadAt`] source and serve the
//! standard [`omicsweave_core::SequenceRead`] access: a
//! [`omicsweave_core::Region`] in, ASCII bases out.
//!
//! Repeated BGZF reads can explicitly attach an externally owned
//! [`omicsweave_core::cache::BlockCache`] with [`FastaReader::set_cache`].
//! [`TwoBitReader::set_cache_with_tile`] enables decoded fixed tiles. Both
//! expose `describe_cache` and `cache_namespace` for external management;
//! ordinary constructors keep caching disabled. Cached sources and their
//! indexes must stay unchanged for the lifetime of a shared source session.
//! Plain FASTA has no decoded-block cache. File-backed readers can
//! `try_clone` their execution context, sharing immutable metadata and the
//! external cache while keeping work buffers independent.
//!
//! Coordinates are 0-based, half-open. Ranges outside a sequence are errors;
//! padding is the caller's decision.
//!
//! # Indexes
//!
//! `.fai` and `.gzi` are the standard samtools/htslib files, byte-identical
//! to what `samtools faidx` and `bgzip` write, and named beside their FASTA
//! by [`fai_path`] (`a.fa.gz` -> `a.fa.gz.fai`) and [`gzi_path`]
//! (`a.fa.gz.gzi`). [`Fai`] and [`Gzi`] share one shape:
//!
//! | | `.fai` | `.gzi` |
//! |---|---|---|
//! | From a source | [`Fai::build`] | [`Gzi::build`] (block headers only) |
//! | Bytes | [`Fai::parse`], [`Fai::to_bytes`] | [`Gzi::parse`], [`Gzi::to_bytes`] |
//! | Beside a file | [`Fai::save`], [`Fai::load`], [`Fai::build_file`] | [`Gzi::save`], [`Gzi::load`], [`Gzi::build_file`] |
//!
//! [`build_index`] builds both in one streaming pass, detecting plain or
//! BGZF storage from the first bytes (plain gzip cannot be indexed, as in
//! htslib); [`build_index_files`] saves them beside the file. Saving is
//! atomic (temporary file, then rename).
//!
//! Indexes can also come from writing: [`FastaWriter`] writes plain or bgzip
//! FASTA and computes the `.fai` (and `.gzi`) from what it writes, with no
//! re-parsing; [`FaiBuilder`] indexes FASTA text fed in chunks;
//! [`BgzfWriter`] writes a bgzip file and its `.gzi`; [`compress_fasta`]
//! compresses existing FASTA text and indexes it in one pass.
//!
//! 2bit is written by [`TwoBitWriter`]: sequences are encoded in pieces
//! ([`TwoBitPiece`], independent, so they can run in parallel) whose runs
//! [`TwoBitRecord`] joins, then written record by record.

mod bgzf;
mod cache;
mod error;
mod fai;
mod faibuild;
mod gzi;
mod index;
mod reader;
mod twobit;
mod twobitwrite;
mod writer;

use std::io::{Read, Write};

pub use bgzf::{BLOCK_TEXT_LEN, BgzfWriter, BlockHeader, EOF_BLOCK, Inflater, MAX_BLOCK_LEN};
pub use error::{Error, Result};
pub use fai::{FAI_SUFFIX, Fai, FaiEntry, fai_path};
pub use faibuild::FaiBuilder;
pub use gzi::{BlockRef, GZI_SUFFIX, Gzi, Span, gzi_path};
pub use index::{build_index, build_index_files};
pub use reader::{FastaReader, TwoBitReader};
/// Default size, in bases, of an explicitly enabled 2bit decoded tile.
pub const DEFAULT_CACHE_TILE_BASES: u32 = 1 << 16;
pub use twobit::{
    Block, Directory, HEADER_LEN as TWOBIT_HEADER_LEN, Parsed, Record, TWOBIT_MAGIC, TwoBitHeader,
};
pub use twobitwrite::{TwoBitPiece, TwoBitRecord, TwoBitWriter, encode_twobit, twobit_base_size};
pub use writer::FastaWriter;

/// Compresses FASTA text from `input` into bgzip `output`, indexing it on the
/// way. Returns the sink, the `.fai` and the `.gzi`.
pub fn compress_fasta<R: Read, W: Write>(
    mut input: R,
    output: W,
    level: i32,
) -> Result<(W, Fai, Gzi)> {
    let mut fai = FaiBuilder::new();
    let mut writer = BgzfWriter::with_level(output, level)?;
    let mut chunk = vec![0u8; 1 << 20];
    loop {
        let n = input.read(&mut chunk)?;
        if n == 0 {
            break;
        }
        fai.push(&chunk[..n])?;
        writer.write_all(&chunk[..n])?;
    }
    let (output, gzi) = writer.finish()?;
    Ok((output, fai.finish()?, gzi))
}
