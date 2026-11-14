//! cooler contact maps (`.cool`, and each resolution of a `.mcool`), stored
//! in HDF5: their bin layout and raw pixels, streamed in parallel for
//! conversion. Balancing weights are neither read nor computed: they are
//! analysis, not file content this crate needs to provide.
//!
//! A cooler map (layout versions 2 and 3) is an HDF5 group holding
//! `chroms/{name,length}`, a bin table `bins/{chrom,start,end}`, the
//! sparse matrix as three pixel columns `pixels/{bin1_id,bin2_id,count}`
//! sorted by `(bin1_id, bin2_id)`, and `indexes/chrom_offset`. A `.mcool`
//! holds one such group per bin size under `/resolutions/<res>`
//! ([`mcool_resolutions`]). Only fixed-width bins are supported; variable
//! bins are [`omicsweave_core::Error::Unsupported`].
//!
//! **Reading.** [`Cool::open`] reads the metadata and bin table through
//! libhdf5 (statically linked HDF5 1.14) and checks the bins are the
//! fixed-width tiling of each chromosome. Pixels are never read whole:
//! [`Cool::pixel_spans`] cuts them into row ranges on chunk boundaries and
//! [`Cool::read_pixels`] decodes one range into a [`PixelBatch`].
//!
//! **Parallelism.** libhdf5 is not reentrant and its Rust bindings hold one
//! process-wide lock around every call, so reading pixels through it would
//! run one thread at a time. libhdf5 is therefore used only for metadata
//! and to locate each pixel chunk (file offset, stored size, filter mask);
//! the chunks themselves are read with `pread` and their filters (deflate,
//! shuffle, Fletcher-32) undone in the calling thread. `&Cool` is shared
//! across threads, each keeping its own [`PixelScratch`]. Columns this path
//! cannot decode (other filters, compact or virtual layout, chunks never
//! written, unusual number types) fall back to libhdf5's reader
//! ([`ReadPath::Library`], serialized); [`Cool::read_paths`] tells which
//! path each column takes.

mod column;
mod cool;
mod filter;
mod h5;

pub use column::{PixelScratch, ReadPath};
pub use cool::{Cool, CoolChrom, PixelBatch, StorageMode, is_mcool, mcool_resolutions};
