//! The standard every OmicsWeave format reader follows.
//!
//! - **Access keys**: [`Region`] — a chromosome name and a 0-based, half-open
//!   `u32` base range — keys sparse queries (features, reads, records);
//!   [`Locus`] — ascending spans of one chromosome on a [`Strand`] — keys
//!   dense queries, whose output runs 5' to 3' of the strand over the
//!   concatenated spans (a [`Region`] converts to a one-span locus).
//! - **Output layouts**, written into caller-allocated memory:
//!
//!   | Family | Reader | Output |
//!   |---|---|---|
//!   | Sequence (`.fa`, `.fa.gz`, `.2bit`) | [`SequenceRead`] | `[u8]` ASCII bases |
//!   | Track (bigWig) | [`TrackRead`] | `[f32]` values, `[f32]` bin statistics ([`BinSpec`]), or `[`[`TrackEntry`]`]` clipped intervals |
//!   | Interval (bigBed) | [`IntervalRead`] | `[`[`IntervalEntry`]`]` + annotation [`StringColumn`] |
//!   | Contact (`.hic`) | [`ContactRead`] | `[f32]` row-major tile, or `[`[`ContactEntry`]`]` nonzero cells |
//!   | Alignment (BAM) | [`AlignmentRead`] | `[`[`AlignmentEntry`]`]` + the [`AlignmentFields`] asked for, as string columns |
//!
//!   Entries are `#[repr(C)]` rows matching numpy structured dtypes; string
//!   columns use Arrow's `LargeString` layout. Variable-size outputs are
//!   decoded once by a query ([`TrackRead::query_entries`],
//!   [`IntervalRead::query`], [`AlignmentRead::query`]) that reports their
//!   size, then copied into the caller's buffers.
//! - **Errors**: one [`Error`] for every reader.
//! - **Bytes**: readers read through [`ReadAt`], which the caller provides.
//! - **Ordered execution**: [`exec`], how writers hand pure work (block
//!   compression, record encoding) to worker threads and take results back
//!   in order.
//! - **Sidecars**: the binary container of derived per-file facts
//!   ([`sidecar`]).
//!
//! Format crates keep their own lower-level steps (computing byte ranges,
//! decoding blocks) for loaders that schedule I/O themselves.

mod alignment;
mod bins;
pub mod cache;
mod contact;
mod error;
pub mod exec;
mod layout;
mod locus;
pub mod orient;
mod read;
mod region;
pub mod sidecar;
mod source;

pub use alignment::{
    AlignmentEntry, AlignmentFields, AlignmentRead, AlignmentSize, AlignmentsOut, CigarColumn,
};
pub use bins::{BinAccumulator, BinLayout, BinSource, BinSpec, Stat, ZoomPiece};
pub use contact::{ContactEntry, ContactQuery, ContactRead, ContactUnit};
pub use error::{Error, Result, check_output};
pub use layout::{
    IntervalEntry, IntervalSize, IntervalWriter, IntervalsOut, StringColumn, Strings, TrackEntry,
};
pub use locus::{Locus, Span, Strand};
pub use read::{IntervalRead, SequenceRead, TrackRead};
pub use region::Region;
pub use source::{FileSource, ReadAt};
