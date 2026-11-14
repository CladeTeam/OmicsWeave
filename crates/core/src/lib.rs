//! The standard every OmicsWeave format reader follows.
//!
//! - **Access key**: [`Region`] — a chromosome name and a 0-based, half-open
//!   `u32` base range.
//! - **Output layouts**, written into caller-allocated memory:
//!
//!   | Family | Reader | Output |
//!   |---|---|---|
//!   | Sequence (`.fa`, `.fa.gz`, `.2bit`) | [`SequenceRead`] | `[u8]` ASCII bases |
//!   | Track (bigWig) | [`TrackRead`] | `[f32]` values, or `[`[`TrackEntry`]`]` clipped intervals |
//!   | Interval (bigBed, BED) | [`IntervalRead`] | `[`[`IntervalEntry`]`]` + annotation [`StringColumn`] |
//!   | Contact (`.hic`) | [`ContactRead`] | `[f32]` row-major tile, or `[`[`ContactEntry`]`]` nonzero cells |
//!
//!   Entries are `#[repr(C)]` rows matching numpy structured dtypes; string
//!   columns use Arrow's `LargeString` layout. Variable-size outputs are
//!   decoded once by a query ([`TrackRead::query_entries`],
//!   [`IntervalRead::query`]) that reports their size, then copied into the
//!   caller's buffers.
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

mod contact;
mod error;
pub mod exec;
mod layout;
mod read;
mod region;
pub mod sidecar;
mod source;

pub use contact::{ContactEntry, ContactQuery, ContactRead, ContactUnit};
pub use error::{Error, Result, check_output};
pub use layout::{
    IntervalEntry, IntervalSize, IntervalWriter, IntervalsOut, StringColumn, Strings, TrackEntry,
};
pub use read::{IntervalRead, SequenceRead, TrackRead};
pub use region::Region;
pub use source::{FileSource, ReadAt};
