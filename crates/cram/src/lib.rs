//! CRAM 3.0 and 3.1, read as htslib reads them.
//!
//! - [`CramReader`]: a CRAM file with its `.crai` ([`CramIndex`]) and,
//!   for slices compressed against one, a reference read through any
//!   [`omicsweave_core::SequenceRead`]; queried by region through the
//!   standard [`omicsweave_core::AlignmentRead`] access. A region returns
//!   the records `sam_itr_queryi` returns, in file order, with the fields
//!   htslib decodes for them (MD and NM generated where it generates
//!   them). Only the blocks of the data series a query needs are
//!   decompressed: the fixed fields and CIGARs need no reference and no
//!   bases, names, qualities or tags.
//! - [`build_index`]: the `.crai` `samtools index` builds.
//! - [`Scanner`]: every container in file order; [`Container::decode_bam`]
//!   decodes one slice into the BAM records htslib reads from it, so
//!   conversions decode slices in parallel.
//!
//! Containers, slices and data series are decoded here; block entropy
//! codecs come from htscodecs (`htscodecs-sys`), gzip from libdeflate,
//! bzip2 and lzma from their `-sys` crates. CRAM 2.x and earlier, and
//! 4.0, are [`omicsweave_core::Error::Unsupported`].

mod codec;
mod comp;
mod container;
mod header;
mod index;
mod md5;
mod reader;
mod record;
mod reference;
mod scan;
mod slice;
mod varint;

pub use container::CRAM_MAGIC;
pub use header::{CramHeader, Reference};
pub use index::{CraiEntry, CramIndex, locate_index};
pub use reader::{CramOptions, CramReader, CramStats};
pub use reference::NoReference;
pub use scan::{
    BamEntry, BamRecords, Container, Scanner, SharedReference, build_index, build_index_file,
    index_entries,
};
pub use slice::References;
