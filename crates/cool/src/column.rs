//! One pixel column (`bin1_id`, `bin2_id` or `count`): where its bytes are
//! in the file and how to turn a row range of it into values.
//!
//! Chunked columns whose filters are all ones [`crate::filter`] undoes, and
//! contiguous columns, are read with `pread` and decoded here, outside
//! libhdf5 and its global lock. Anything else (compact or virtual layout,
//! other filters, chunks never written, unusual element types) goes through
//! libhdf5's own reader, which is serialized.

use std::ops::Range;

use hdf5_metno as hdf5;
use libdeflater::Decompressor;
use omicsweave_core::{Error, FileSource, ReadAt, Result};

use crate::filter::{self, Stage};
use crate::h5::{self, Elem, ElemKind};

/// How a pixel column is read.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ReadPath {
    /// Chunks read with `pread` and decoded in the calling thread.
    Chunked,
    /// One contiguous extent read with `pread`.
    Contiguous,
    /// libhdf5's reader, serialized by its global lock.
    Library,
}

/// A chunk's place in the file.
#[derive(Debug, Clone, Copy)]
struct ChunkLoc {
    addr: u64,
    size: u64,
    filter_mask: u32,
}

enum Store {
    Chunked {
        chunk_len: u64,
        chunk_bytes: usize,
        pipeline: Vec<Stage>,
        filter_bits: u32,
        chunks: Vec<ChunkLoc>,
    },
    Contiguous {
        offset: u64,
    },
    Library,
}

pub(crate) struct Column {
    name: &'static str,
    ds: hdf5::Dataset,
    len: u64,
    elem: Option<Elem>,
    store: Store,
}

/// Per-thread working state for decoding chunks.
pub struct PixelScratch {
    inflater: Decompressor,
    raw: Vec<u8>,
    spare: Vec<u8>,
}

impl PixelScratch {
    pub fn new() -> PixelScratch {
        PixelScratch {
            inflater: Decompressor::new(),
            raw: Vec::new(),
            spare: Vec::new(),
        }
    }
}

impl Default for PixelScratch {
    fn default() -> PixelScratch {
        PixelScratch::new()
    }
}

/// Bytes per piece when reading a contiguous extent, which bounds the
/// scratch buffer however many rows are asked for.
const CONTIGUOUS_PIECE: u64 = 1 << 22;

/// Values a column decodes into.
pub(crate) trait Value: Copy + hdf5::H5Type {
    /// Appends the elements in `bytes` (whole elements of `elem`).
    fn extend(elem: Elem, bytes: &[u8], out: &mut Vec<Self>, what: &'static str) -> Result<()>;
    /// Converts only the selected elements from HDF5 byte-shuffle planes.
    fn extend_shuffled(
        elem: Elem,
        bytes: &[u8],
        rows: Range<usize>,
        out: &mut Vec<Self>,
        what: &'static str,
    ) -> Result<()>;
    /// Converts a value read by libhdf5 as `Self::Wide`.
    type Wide: Copy + hdf5::H5Type;
    fn from_wide(w: Self::Wide, what: &'static str) -> Result<Self>;
}

macro_rules! element_type {
    ($elem:expr, $run:ident) => {
        match ($elem.kind, $elem.size) {
            (ElemKind::Signed, 1) => $run!(i8, 1),
            (ElemKind::Signed, 2) => $run!(i16, 2),
            (ElemKind::Signed, 4) => $run!(i32, 4),
            (ElemKind::Signed, 8) => $run!(i64, 8),
            (ElemKind::Unsigned, 1) => $run!(u8, 1),
            (ElemKind::Unsigned, 2) => $run!(u16, 2),
            (ElemKind::Unsigned, 4) => $run!(u32, 4),
            (ElemKind::Unsigned, 8) => $run!(u64, 8),
            (ElemKind::Float, 4) => $run!(f32, 4),
            (ElemKind::Float, 8) => $run!(f64, 8),
            _ => unreachable!("element types are checked at open"),
        }
    };
}

macro_rules! each_elem {
    ($elem:expr, $bytes:expr, |$v:ident| $body:expr) => {{
        let e = $elem;
        macro_rules! run {
            ($t:ty, $n:literal) => {
                for b in $bytes.chunks_exact($n) {
                    let b: [u8; $n] = b.try_into().unwrap();
                    let $v = if e.big_endian {
                        <$t>::from_be_bytes(b)
                    } else {
                        <$t>::from_le_bytes(b)
                    };
                    $body;
                }
            };
        }
        element_type!(e, run)
    }};
}

macro_rules! each_shuffled_elem {
    ($elem:expr, $bytes:expr, $rows:expr, |$v:ident| $body:expr) => {{
        let e = $elem;
        let bytes = $bytes;
        let rows = $rows;
        let n = bytes.len().checked_div(e.size).unwrap_or(0);
        if rows.start > rows.end || rows.end > n {
            return Err(Error::invalid(
                "hdf5 shuffle",
                "element range outside byte planes",
            ));
        }
        macro_rules! run {
            ($t:ty, $size:literal) => {
                for i in rows {
                    let b: [u8; $size] = std::array::from_fn(|j| bytes[j * n + i]);
                    let $v = if e.big_endian {
                        <$t>::from_be_bytes(b)
                    } else {
                        <$t>::from_le_bytes(b)
                    };
                    $body;
                }
            };
        }
        element_type!(e, run)
    }};
}

fn negative(what: &'static str) -> Error {
    Error::invalid(what, "negative bin id")
}

impl Value for u64 {
    type Wide = i64;

    #[allow(clippy::unnecessary_cast)]
    fn extend(elem: Elem, bytes: &[u8], out: &mut Vec<u64>, what: &'static str) -> Result<()> {
        if elem.kind == ElemKind::Float {
            return Err(Error::invalid(what, "bin ids are not integers"));
        }
        let mut bad = false;
        each_elem!(elem, bytes, |v| {
            bad |= (v as i128) < 0;
            out.push(v as u64)
        });
        if bad {
            return Err(negative(what));
        }
        Ok(())
    }

    fn from_wide(w: i64, what: &'static str) -> Result<u64> {
        u64::try_from(w).map_err(|_| negative(what))
    }

    #[allow(clippy::unnecessary_cast)]
    fn extend_shuffled(
        elem: Elem,
        bytes: &[u8],
        rows: Range<usize>,
        out: &mut Vec<u64>,
        what: &'static str,
    ) -> Result<()> {
        if elem.kind == ElemKind::Float {
            return Err(Error::invalid(what, "bin ids are not integers"));
        }
        let mut bad = false;
        each_shuffled_elem!(elem, bytes, rows, |v| {
            bad |= (v as i128) < 0;
            out.push(v as u64)
        });
        if bad {
            return Err(negative(what));
        }
        Ok(())
    }
}

impl Value for f64 {
    type Wide = f64;

    fn extend(elem: Elem, bytes: &[u8], out: &mut Vec<f64>, _: &'static str) -> Result<()> {
        each_elem!(elem, bytes, |v| out.push(v as f64));
        Ok(())
    }

    fn from_wide(w: f64, _: &'static str) -> Result<f64> {
        Ok(w)
    }

    fn extend_shuffled(
        elem: Elem,
        bytes: &[u8],
        rows: Range<usize>,
        out: &mut Vec<f64>,
        _: &'static str,
    ) -> Result<()> {
        each_shuffled_elem!(elem, bytes, rows, |v| out.push(v as f64));
        Ok(())
    }
}

impl Column {
    /// Opens a pixel column and works out how to read it.
    pub(crate) fn open(
        group: &hdf5::Group,
        path: &'static str,
        file_len: u64,
        integer: bool,
    ) -> Result<Column> {
        let ds = h5::dataset(group, path)?;
        let len = h5::len_1d(&ds, path)?;
        if integer {
            h5::check_integer(&ds, path)?;
        }
        let elem = h5::elem(&ds, path)?;
        let store = match elem {
            Some(e) if len > 0 => store(&ds, path, len, e, file_len)?,
            _ => Store::Library,
        };
        Ok(Column {
            name: path,
            ds,
            len,
            elem,
            store,
        })
    }

    pub(crate) fn len(&self) -> u64 {
        self.len
    }

    pub(crate) fn path(&self) -> ReadPath {
        match self.store {
            Store::Chunked { .. } => ReadPath::Chunked,
            Store::Contiguous { .. } => ReadPath::Contiguous,
            Store::Library => ReadPath::Library,
        }
    }

    /// Rows per chunk when the column is read chunk by chunk.
    pub(crate) fn chunk_len(&self) -> Option<u64> {
        match self.store {
            Store::Chunked { chunk_len, .. } => Some(chunk_len),
            _ => None,
        }
    }

    /// Appends rows `rows` (within the column) to `out`. Each chunk the
    /// range touches is read and decoded once.
    pub(crate) fn read<T: Value>(
        &self,
        source: &FileSource,
        rows: Range<u64>,
        scratch: &mut PixelScratch,
        out: &mut Vec<T>,
    ) -> Result<()> {
        if rows.is_empty() {
            return Ok(());
        }
        let what = self.name;
        match &self.store {
            Store::Library => self.read_library(rows, out),
            Store::Contiguous { offset } => {
                let elem = self.elem.expect("fast paths have a decodable type");
                let size = elem.size as u64;
                let per_piece = (CONTIGUOUS_PIECE / size).max(1);
                let mut at = rows.start;
                while at < rows.end {
                    let n = per_piece.min(rows.end - at);
                    scratch.raw.resize((n * size) as usize, 0);
                    source.read_exact_at(&mut scratch.raw, offset + at * size)?;
                    T::extend(elem, &scratch.raw, out, what)?;
                    at += n;
                }
                Ok(())
            }
            Store::Chunked {
                chunk_len,
                chunk_bytes,
                pipeline,
                filter_bits,
                chunks,
            } => {
                let elem = self.elem.expect("fast paths have a decodable type");
                let size = elem.size as u64;
                let first = rows.start / chunk_len;
                let last = (rows.end - 1) / chunk_len;
                for k in first..=last {
                    let loc = chunks[k as usize];
                    let base = k * chunk_len;
                    let lo = rows.start.max(base) - base;
                    let hi = (rows.end - base).min(*chunk_len);
                    if loc.filter_mask & *filter_bits == *filter_bits {
                        // Unfiltered chunks are raw elements: read just the
                        // rows needed.
                        if hi * size > loc.size {
                            return Err(Error::Truncated {
                                what: "hdf5 chunk",
                                need: hi * size,
                                have: loc.size,
                            });
                        }
                        scratch.raw.resize(((hi - lo) * size) as usize, 0);
                        source.read_exact_at(&mut scratch.raw, loc.addr + lo * size)?;
                        T::extend(elem, &scratch.raw, out, what)?;
                        continue;
                    }
                    scratch.raw.resize(loc.size as usize, 0);
                    source.read_exact_at(&mut scratch.raw, loc.addr)?;
                    // The first write-stage shuffle is the last decode
                    // stage. Read its byte planes directly into output,
                    // avoiding a full-chunk reorder and intermediate buffer.
                    // Earlier checksums forbid this fusion because they
                    // must see the complete unshuffled byte stream.
                    let shuffled =
                        pipeline.first() == Some(&Stage::Shuffle) && loc.filter_mask & 1 == 0;
                    let (stages, mask) = if shuffled {
                        (&pipeline[1..], loc.filter_mask >> 1)
                    } else {
                        (pipeline.as_slice(), loc.filter_mask)
                    };
                    let in_spare = filter::decode(
                        stages,
                        mask,
                        *chunk_bytes,
                        elem.size,
                        &mut scratch.raw,
                        &mut scratch.spare,
                        &mut scratch.inflater,
                    )?;
                    let decoded = if in_spare {
                        &scratch.spare
                    } else {
                        &scratch.raw
                    };
                    if shuffled {
                        T::extend_shuffled(elem, decoded, lo as usize..hi as usize, out, what)?;
                    } else {
                        let bytes = &decoded[(lo * size) as usize..(hi * size) as usize];
                        T::extend(elem, bytes, out, what)?;
                    }
                }
                Ok(())
            }
        }
    }

    /// Reads rows through libhdf5, converting to a wide type first.
    pub(crate) fn read_library<T: Value>(&self, rows: Range<u64>, out: &mut Vec<T>) -> Result<()> {
        if rows.is_empty() {
            return Ok(());
        }
        let values = self
            .ds
            .read_slice_1d::<T::Wide, _>(rows.start as usize..rows.end as usize)
            .map_err(h5::err(self.name))?;
        out.reserve(values.len());
        for &w in values.iter() {
            out.push(T::from_wide(w, self.name)?);
        }
        Ok(())
    }
}

/// Works out the fast path for a column, or [`Store::Library`] when there
/// is none.
fn store(
    ds: &hdf5::Dataset,
    what: &'static str,
    len: u64,
    elem: Elem,
    file_len: u64,
) -> Result<Store> {
    use hdf5::dataset::Layout;
    let size = elem.size as u64;
    match ds.layout() {
        Layout::Contiguous => {
            let Some(offset) = ds.offset() else {
                return Ok(Store::Library);
            };
            let bytes = len
                .checked_mul(size)
                .ok_or_else(|| Error::invalid(what, "column byte size overflow"))?;
            check_extent(what, offset, bytes, file_len)?;
            Ok(Store::Contiguous { offset })
        }
        Layout::Chunked => {
            let dcpl = ds.dcpl().map_err(h5::err(what))?;
            if dcpl
                .get_chunk_opts()
                .map_err(h5::err(what))?
                .is_some_and(|o| !o.is_empty())
            {
                // Partial edge chunks may be stored unfiltered.
                return Ok(Store::Library);
            }
            let Some(chunk_len) = ds.chunk().and_then(|c| c.first().copied()) else {
                return Ok(Store::Library);
            };
            let chunk_len = chunk_len as u64;
            if chunk_len == 0 {
                return Err(Error::invalid(what, "zero chunk size"));
            }
            let chunk_bytes = chunk_len
                .checked_mul(size)
                .and_then(|n| usize::try_from(n).ok())
                .filter(|&n| n <= isize::MAX as usize)
                .ok_or_else(|| Error::invalid(what, "chunk byte size overflow"))?;
            let mut pipeline = Vec::new();
            let mut has_checksum = false;
            for f in dcpl.get_filters().map_err(h5::err(what))? {
                pipeline.push(match f {
                    hdf5::filters::Filter::Deflate(_) if has_checksum => {
                        Stage::DeflateWithChecksums
                    }
                    hdf5::filters::Filter::Deflate(_) => Stage::Deflate,
                    hdf5::filters::Filter::Shuffle => Stage::Shuffle,
                    hdf5::filters::Filter::Fletcher32 => {
                        has_checksum = true;
                        Stage::Fletcher32
                    }
                    _ => return Ok(Store::Library),
                });
            }
            let filter_bits = match pipeline.len() {
                0..32 => (1u32 << pipeline.len()) - 1,
                32 => u32::MAX,
                _ => return Ok(Store::Library),
            };
            let n = usize::try_from(len.div_ceil(chunk_len))
                .map_err(|_| Error::invalid(what, "chunk index size overflow"))?;
            // An absent chunk requires library fill-value handling for
            // this column. Detect that before allocating a directory for
            // a potentially enormous, mostly unwritten logical extent.
            if ds.num_chunks().is_some_and(|allocated| allocated < n) {
                return Ok(Store::Library);
            }
            // HDF5 reserves HADDR_UNDEF for an absent address; it cannot
            // name any allocated chunk. Keep holes in the final index so
            // constructing it needs only one allocation.
            const MISSING: ChunkLoc = ChunkLoc {
                addr: hdf5_metno_sys::h5::HADDR_UNDEF,
                size: 0,
                filter_mask: 0,
            };
            let mut chunks = Vec::new();
            chunks
                .try_reserve_exact(n)
                .map_err(|_| Error::invalid(what, "chunk index allocation failed"))?;
            chunks.resize(n, MISSING);
            let mut misplaced = false;
            let mut duplicate = false;
            ds.chunks_visit(|info| {
                let k = info.offset[0] / chunk_len;
                if info.offset[0] % chunk_len != 0 {
                    misplaced = true;
                } else if let Some(slot) = usize::try_from(k).ok().and_then(|k| chunks.get_mut(k)) {
                    duplicate |= slot.addr != MISSING.addr;
                    *slot = ChunkLoc {
                        addr: info.addr,
                        size: info.size,
                        filter_mask: info.filter_mask,
                    };
                }
                0
            })
            .map_err(h5::err(what))?;
            if misplaced {
                return Err(Error::invalid(what, "chunk offset off the chunk grid"));
            }
            if duplicate {
                return Err(Error::invalid(what, "duplicate chunk offset"));
            }
            for loc in &chunks {
                // A chunk never written holds the fill value, which only
                // libhdf5 supplies.
                if loc.addr == MISSING.addr {
                    return Ok(Store::Library);
                }
                if usize::try_from(loc.size).is_err() || loc.size > isize::MAX as u64 {
                    return Err(Error::invalid(what, "stored chunk size overflow"));
                }
                check_extent(what, loc.addr, loc.size, file_len)?;
            }
            Ok(Store::Chunked {
                chunk_len,
                chunk_bytes,
                pipeline,
                filter_bits,
                chunks,
            })
        }
        _ => Ok(Store::Library),
    }
}

fn check_extent(what: &'static str, offset: u64, len: u64, file_len: u64) -> Result<()> {
    if offset.checked_add(len).is_none_or(|e| e > file_len) {
        return Err(Error::OutOfRange {
            what,
            offset,
            lo: 0,
            hi: file_len,
        });
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    fn shuffle(bytes: &[u8], size: usize) -> Vec<u8> {
        let n = bytes.len() / size;
        let mut shuffled = Vec::with_capacity(bytes.len());
        for j in 0..size {
            for i in 0..n {
                shuffled.push(bytes[i * size + j]);
            }
        }
        shuffled.extend_from_slice(&bytes[n * size..]);
        shuffled
    }

    #[test]
    fn shuffled_integer_planes_select_exact_rows_and_both_byte_orders() {
        let values = [0u64, 1, 7, 41, 120];
        for kind in [ElemKind::Signed, ElemKind::Unsigned] {
            for size in [1, 2, 4, 8] {
                for big_endian in [false, true] {
                    let elem = Elem {
                        kind,
                        size,
                        big_endian,
                    };
                    let mut encoded = Vec::new();
                    for v in values {
                        if big_endian {
                            encoded.extend_from_slice(&v.to_be_bytes()[8 - size..]);
                        } else {
                            encoded.extend_from_slice(&v.to_le_bytes()[..size]);
                        }
                    }
                    // Byte-shuffle leaves an incomplete final element intact.
                    if size > 1 {
                        encoded.push(0xa5);
                    }
                    let bytes = shuffle(&encoded, size);
                    for rows in [0..0, 0..1, 1..4, 4..5, 0..5] {
                        let mut ids = Vec::new();
                        u64::extend_shuffled(elem, &bytes, rows.clone(), &mut ids, "ids").unwrap();
                        assert_eq!(ids, values[rows.clone()]);
                        let mut counts = Vec::new();
                        f64::extend_shuffled(elem, &bytes, rows.clone(), &mut counts, "count")
                            .unwrap();
                        assert_eq!(
                            counts,
                            values[rows].iter().map(|&v| v as f64).collect::<Vec<_>>()
                        );
                    }
                    assert!(
                        u64::extend_shuffled(elem, &bytes, 0..6, &mut Vec::new(), "ids").is_err()
                    );
                }
            }
        }
    }

    #[test]
    fn shuffled_float_planes_preserve_values_and_reject_float_ids() {
        for size in [4, 8] {
            for big_endian in [false, true] {
                let elem = Elem {
                    kind: ElemKind::Float,
                    size,
                    big_endian,
                };
                let values = [0.0f64, -7.25, 1.5, f64::INFINITY, f64::NAN];
                let mut encoded = Vec::new();
                for v in values {
                    if size == 4 {
                        encoded.extend_from_slice(&if big_endian {
                            (v as f32).to_be_bytes()
                        } else {
                            (v as f32).to_le_bytes()
                        });
                    } else {
                        encoded.extend_from_slice(&if big_endian {
                            v.to_be_bytes()
                        } else {
                            v.to_le_bytes()
                        });
                    }
                }
                let bytes = shuffle(&encoded, size);
                let mut out = Vec::new();
                f64::extend_shuffled(elem, &bytes, 1..5, &mut out, "count").unwrap();
                for (actual, expected) in out.iter().zip(&values[1..]) {
                    assert_eq!(actual.to_bits(), expected.to_bits());
                }
                assert!(u64::extend_shuffled(elem, &bytes, 0..1, &mut Vec::new(), "ids").is_err());
            }
        }
    }

    #[test]
    fn shuffled_negative_id_is_rejected_even_when_other_rows_are_valid() {
        let elem = Elem {
            kind: ElemKind::Signed,
            size: 4,
            big_endian: false,
        };
        let bytes = shuffle(
            &[1i32, -1, 2]
                .into_iter()
                .flat_map(i32::to_le_bytes)
                .collect::<Vec<_>>(),
            4,
        );
        assert!(u64::extend_shuffled(elem, &bytes, 0..3, &mut Vec::new(), "ids").is_err());
        let mut out = Vec::new();
        u64::extend_shuffled(elem, &bytes, 2..3, &mut out, "ids").unwrap();
        assert_eq!(out, [2]);
    }
}
