//! BAM, read as htslib reads it.
//!
//! - [`BamReader`]: a BAM file with its BAI or CSI index ([`BinIndex`]),
//!   queried by region through the standard
//!   [`omicsweave_core::AlignmentRead`] access. The
//!   records of a region are exactly those `sam_itr_queryi` returns, in
//!   file order. Indexed chunks determine the physical BGZF blocks read;
//!   a caller-owned block cache can share them across execution contexts.
//! - [`build_index`]: the BAI or CSI `samtools index` builds for a sorted
//!   BAM.
//! - [`Scanner`]: sequential records; [`Scanner::new_streaming`] omits the
//!   file-wide block directory for bounded one-pass consumers, while
//!   [`Scanner::new`] retains it for index building.
//! - [`Record`]: one record read in place, checked as `bam_read1` checks
//!   it; a long CIGAR held in a `CG` tag replaces its placeholder.
//!
//! Records are returned as stored: nothing is filtered, counted or
//! realigned.

mod header;
mod reader;
mod record;
mod scan;
mod stream;

pub use header::{BAM_MAGIC, BamHeader, Reference};
pub use omicsweave_bgzf::{BinIndex, Chunk, IndexKind, locate_index};
pub use reader::BamReader;
pub use record::{NT16, Record, consumes_query, consumes_reference, flag, op, unpack_seq};
pub use scan::{IndexChoice, Scanner, build_index, build_index_file, inflate_queue};
