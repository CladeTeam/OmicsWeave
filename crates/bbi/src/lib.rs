//! bigWig and bigBed (BBI) parsing on caller-provided buffers.
//!
//! The low-level parsers perform no I/O and own no numeric data. Each step either
//! describes what the caller must read and allocate, or operates on bytes the
//! caller already holds:
//!
//! 1. Read [`HEADER_LEN`] bytes at offset 0 and [`Header::parse`] them.
//! 2. Read [`Header::zoom_headers`] and [`ZoomLevel::parse_all`].
//! 3. Build a [`Layout`] from the header, zoom levels and file length; it
//!    gives the byte ranges of the remaining sections.
//! 4. Read [`Layout::chrom_tree`] and [`Chroms::parse`] it.
//! 5. Read [`Layout::index`] into a buffer kept for as long as queries run,
//!    and wrap it in an [`Index`].
//! 6. [`Index::query`] a range to get the [`BlockRef`]s to read.
//! 7. Read each block; if [`Header::compressed`], inflate it with an
//!    [`Inflater`] into a buffer of [`Header::block_buffer_len`] bytes.
//! 8. Parse the bytes as a [`WigBlock`], [`BedBlock`] or [`ZoomBlock`] and
//!    read items or [`WigBlock::fill`] a caller array.
//!
//! Steps 1-5 run once per file; 6-8 run per query.
//! [`BbiReader`] composes
//! them over a [`omicsweave_core::ReadAt`] source and serves the standard
//! [`omicsweave_core::TrackRead`] (bigWig) and
//! [`omicsweave_core::IntervalRead`] (bigBed) access. Coordinates are 0-based,
//! half-open, as stored in the file. Only little-endian files are supported.
//!
//! **Writing and sidecars.** [`BigWigWriter`] and [`BigBedWriter`] write
//! standard version-4 files in one pass and return the [`BbiSidecar`] of the
//! file: every block's span, file range, decompressed size and item count,
//! which the format itself does not record. [`BbiSidecar::build`] makes the
//! same sidecar from any existing file; it is stored next to the file as
//! `<file>.sidecar` ([`omicsweave_core::sidecar`]).
//!
//! # Reusing index and block bytes
//!
//! Parsing stays cheap for one-shot reads. When index bytes or decompressed
//! blocks are retained for many windows, prepare their borrowed views once:
//! [`Index::prepare`], [`WigBlock::prepare`], [`BedBlock::prepare`] and
//! [`ZoomBlock::prepare`]. Preparation certifies actual ordering rather than
//! assuming that starts or ends are monotone. Unsorted and nested intervals
//! retain stored-order behavior; long enclosing intervals are never skipped.
//!
//! Index node metadata and bigBed record offsets go into caller-provided
//! vectors, cleared and filled in the same single pass that validates the
//! bytes; reusing them avoids allocating after warm-up. Their mutable borrows
//! last through preparation; the resulting view holds immutable borrows, so
//! bytes and metadata cannot diverge while in use. No prepared view owns
//! index bytes, records or numeric output.
//! [`BbiReader`] prepares and retains its file-level node metadata when opened;
//! its ordinary read methods reuse that preparation automatically.
//!
//! ```
//! use omicsweave_bbi::{Index, NodeInfo, Result};
//!
//! fn repeated_queries(bytes: &[u8], base: u64) -> Result<()> {
//!     let raw = Index::new(bytes, base)?;
//!     let mut nodes: Vec<NodeInfo> = Vec::new();
//!     let index = raw.prepare(&mut nodes)?;
//!     for start in [100, 200, 300] {
//!         for block in index.query(0, start, start + 64) {
//!             let range = block?.range();
//!             // The caller schedules I/O for this range.
//!             std::hint::black_box(range);
//!         }
//!     }
//!     Ok(())
//! }
//! ```
//!
//! Keep the prepared view beside the cached bytes across queries; preparing
//! anew for every window still incurs a full scan. Preparation never takes a
//! separate counting pass. The original unprepared APIs remain available for
//! one-shot operations.
//!
//! # Work bounds
//!
//! Let N be block records, L output bases, K matching records, C conservative
//! candidates, W total covered lengths including overwrites, B node fanout,
//! T tree nodes, V visited nodes, and Z zoom levels.
//!
//! | Path | Work and storage |
//! | --- | --- |
//! | Header / zoom headers | O(1) / O(Z), checked metadata reads |
//! | Layout | O(Z log Z) to order section boundaries, O(Z) metadata |
//! | Chromosome tree | O(encoded items + names + T log T), O(T + depth) validation state |
//! | Unprepared R-tree query | Linear in items of visited nodes; no ordering assumed |
//! | Index preparation | O(items + T log T), T caller-owned [`NodeInfo`] slots |
//! | Prepared R-tree query | O((V-1) log T + sum(log B + C)), direct root lookup, no per-query allocation |
//! | Prepared Wig / Zoom / bigBed range | O(log N + C); C = K for monotone starts and ends |
//! | Arbitrary nested bounds | Max-span pruning where possible; C can still be N |
//! | fixedStep range | O(1 + K) using the arithmetic progression |
//! | fixedStep fill | O(1 + changed values + written bases), each base written once |
//! | Prepared monotone Wig fill | O(log N + C + L), overwritten suffixes skipped |
//! | General [`WigBlock::fill`] | O(N + W); accepts arbitrary stored order |
//! | [`WigBlock::fill_with_scratch`] | O((N + L) alpha(L)), 3(L+1) caller-owned usize slots |
//! | Prepared bigBed item / iterator nth | O(1), independent of preceding text length |
//! | Inflate / full record traversal | Linear in bytes / records, unavoidable output work |
//!
//! For general heavily overlapping input, `fill_with_scratch` traverses items
//! backwards and uses a disjoint-set forest to skip already-written bases.
//! It preserves last-item-wins values and the logical write count. Its O(L)
//! scratch initialization can cost more than `fill` for sparse windows.
//!
//! Run `cargo bench -p omicsweave-bbi --bench paths` for warmed in-memory
//! measurements. It reports preparation, one-shot parsing and repeated-query
//! paths separately at multiple sizes, including overlap-heavy inputs.
//! These measurements exclude file I/O and decompression.

mod block;
mod bytes;
mod chroms;
mod error;
mod header;
mod index;
mod inflate;
mod layout;
mod reader;
mod search;
mod sidecar;
mod writer;

pub use block::{
    BedBlock, BedEntries, BedEntry, PreparedBedBlock, Summary, WIG_HEADER_LEN, WigBlock, WigKind,
    ZOOM_RECORD_LEN, ZoomBlock,
};
pub use chroms::{CHROM_TREE_MAGIC, Chrom, ChromId, Chroms};
pub use error::{Error, Result};
pub use header::{
    BIGBED_MAGIC, BIGWIG_MAGIC, HEADER_LEN, Header, Kind, TOTAL_SUMMARY_LEN, TotalSummary,
    ZOOM_HEADER_LEN, ZoomLevel, parse_auto_sql,
};
pub use index::{BlockRef, INDEX_MAGIC, Index, IndexHeader, Leaf, MAX_DEPTH, NodeInfo, Query};
pub use inflate::Inflater;
pub use layout::{Layout, Level};
pub use reader::BbiReader;
pub use sidecar::{BBI_SIDECAR_VERSION, BbiSidecar, BlockInfo, LevelBlocks};
pub use writer::{
    BBI_WRITE_VERSION, BedSchema, BigBedWriter, BigWigWriter, WriteOptions, ZOOM_SAMPLE,
};
