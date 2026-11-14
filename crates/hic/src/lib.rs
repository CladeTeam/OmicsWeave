//! Juicer `.hic` contact maps, versions 8 and 9 (v9 is the current
//! specification, <https://github.com/aidenlab/hic-format>): random access
//! to tiles of raw contact counts, and writing. Other versions are
//! [`omicsweave_core::Error::Unsupported`]. Normalization and expected
//! vectors are neither read nor computed: they are analysis, not file
//! content this crate needs to provide.
//!
//! A `.hic` file stores each chromosome pair's contact map at several bin
//! sizes, in base pairs or restriction fragments, cut into zlib-compressed
//! blocks: square grid tiles, or (v9, within a chromosome) tiles laid along
//! the diagonal by distance from it. Only one half of each symmetric
//! intra-chromosomal map is stored.
//!
//! **Reading.** [`HicReader`] serves the standard
//! [`omicsweave_core::ContactRead`] access: a tile of any two regions in, a
//! dense row-major matrix or its nonzero cells out, in base pairs or
//! fragments, as straw reads raw (`NONE`) counts. Tiles come back whole in
//! the orientation asked for, whichever half the file stores. For loaders
//! that schedule I/O themselves, the lower-level pieces are public:
//! [`ZoomData::blocks_for`] gives the blocks a tile needs (stable block
//! numbers, usable as cache keys), [`Inflater`] and [`decode_block`] turn a
//! block's bytes into records.
//!
//! **Writing.** [`HicWriter`] writes binned maps the caller provides,
//! cut into blocks of a configurable size ([`BlockLayout`]; by default
//! [`DEFAULT_INTRA_BLOCK_BINS`]-bin blocks within a chromosome, the Juicer
//! grid between chromosomes) and encoding each block in its smallest
//! representation; [`HicBuilder`] bins
//! contacts at every resolution first, as the Juicer tools' `pre` does.
//! Files carry empty expected and normalization lists. A map can also be
//! streamed block by block ([`HicWriter::matrix_layout`] places each cell
//! in its block, [`HicWriter::add_block`] takes one block's records), so it
//! never has to be in memory whole; blocks are encoded and compressed on an
//! [`omicsweave_core::exec::Exec`] and written in order, the same file for
//! every executor.
//!
//! **Sidecars.** A `.hic` file does not record how large a block is once
//! inflated or how many records it holds. [`HicSidecar`] records both for
//! every block: [`HicWriter::finish_with_sidecar`] and
//! [`HicBuilder::finish_with_sidecar`] return it for the file they write,
//! [`HicSidecar::build`] measures any existing file, and it is stored next
//! to the file as `<file>.sidecar` ([`omicsweave_core::sidecar`]).

mod block;
mod build;
mod cursor;
mod matrix;
mod meta;
mod reader;
mod sidecar;
mod writer;

pub use block::{Deflater, Inflater, Record, decode_block, encode_block};
pub use build::{BuildOptions, HicBuilder};
pub use matrix::{BlockEntry, Grid, ZoomData};
pub use meta::{
    Chrom, DEPTH_BASE_ATTRIBUTE, Depth, Footer, Header, MAX_VERSION, MIN_VERSION, SiteList,
};
pub use reader::HicReader;
pub use sidecar::{BlockSize, HIC_SIDECAR_VERSION, HicSidecar, MatrixBlocks};
pub use writer::{
    BlockBins, BlockLayout, DEFAULT_INTRA_BLOCK_BINS, HicHeader, HicWriter, MatrixLayout,
    MatrixLayouts,
};
