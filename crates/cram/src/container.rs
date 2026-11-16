//! The CRAM file definition, container headers and blocks, read as
//! htslib's `cram_read_file_def`, `cram_read_container` and
//! `cram_read_block` read them, and block decompression
//! (`cram_uncompress_block`).

use std::ffi::{c_char, c_int, c_uint};
use std::ptr::NonNull;

use omicsweave_core::{Error, Result};

use crate::varint::{int32, itf8, ltf8};

// The bzip2 decoder below is libbz2's; referencing the crate links it.
use bzip2_sys as _;

/// The magic of a CRAM file.
pub const CRAM_MAGIC: [u8; 4] = *b"CRAM";

/// Bytes of the file definition.
pub(crate) const FILE_DEF_LEN: usize = 26;

/// Block content types.
pub(crate) mod content {
    pub const COMPRESSION_HEADER: u8 = 1;
    pub const MAPPED_SLICE: u8 = 2;
    pub const UNMAPPED_SLICE: u8 = 3;
    pub const EXTERNAL: u8 = 4;
    pub const CORE: u8 = 5;
}

/// Block compression methods.
pub(crate) mod method {
    pub const RAW: u8 = 0;
    pub const GZIP: u8 = 1;
    pub const BZIP2: u8 = 2;
    pub const LZMA: u8 = 3;
    pub const RANS: u8 = 4;
    pub const RANS_NX16: u8 = 5;
    pub const ARITH: u8 = 6;
    pub const FQZ: u8 = 7;
    pub const TOK3: u8 = 8;
}

/// The CRAM version of a file definition. Only 3.0 and 3.1 are read.
pub(crate) fn file_version(def: &[u8]) -> Result<(u8, u8)> {
    if def.len() < FILE_DEF_LEN {
        return Err(Error::Truncated {
            what: "CRAM file definition",
            need: FILE_DEF_LEN as u64,
            have: def.len() as u64,
        });
    }
    if def[..4] != CRAM_MAGIC {
        return Err(Error::BadMagic {
            what: "CRAM file definition",
            found: u32::from_le_bytes(def[..4].try_into().unwrap()),
        });
    }
    match (def[4], def[5]) {
        (3, minor @ (0 | 1)) => Ok((3, minor)),
        (major, _) if major < 3 => Err(Error::Unsupported("CRAM before version 3.0")),
        _ => Err(Error::Unsupported("CRAM versions after 3.1")),
    }
}

/// A container header (CRAM 3).
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct ContainerHeader {
    /// Bytes of the container after its header.
    pub length: i32,
    /// Reference id; -1 unmapped, -2 several.
    pub ref_seq_id: i32,
    /// 1-based first position of the container.
    pub ref_seq_start: i64,
    pub ref_seq_span: i64,
    pub num_records: i32,
    pub record_counter: i64,
    pub num_bases: i64,
    pub num_blocks: i32,
    /// Slice offsets from the end of the header.
    pub landmarks: Vec<i32>,
    /// Bytes of the header itself.
    pub header_len: usize,
}

impl ContainerHeader {
    /// Parses the header at the start of `buf`. `Ok(None)` when `buf` is
    /// shorter than the 4-byte length (the end of the file, as htslib
    /// reads it); a header cut short after that is [`Error::Truncated`],
    /// which tells a caller that read only part of the file to read more.
    pub(crate) fn parse(buf: &[u8]) -> Result<Option<ContainerHeader>> {
        if buf.len() < 4 {
            return Ok(None);
        }
        let what = "CRAM container header";
        let mut b = buf;
        let length = int32(&mut b, what)?;
        let ref_seq_id = itf8(&mut b, what)?;
        let ref_seq_start = itf8(&mut b, what)? as i64;
        let ref_seq_span = itf8(&mut b, what)? as i64;
        let num_records = itf8(&mut b, what)?;
        let record_counter = ltf8(&mut b, what)?;
        let num_bases = ltf8(&mut b, what)?;
        let num_blocks = itf8(&mut b, what)?;
        let num_landmarks = itf8(&mut b, what)?;
        if !(0..=1_000_000).contains(&num_landmarks) {
            return Err(Error::invalid(what, format!("{num_landmarks} landmarks")));
        }
        let mut landmarks = Vec::with_capacity(num_landmarks as usize);
        for _ in 0..num_landmarks {
            landmarks.push(itf8(&mut b, what)?);
        }
        let covered = buf.len() - b.len();
        let crc = int32(&mut b, what)? as u32;
        if libdeflater::crc32(&buf[..covered]) != crc {
            return Err(Error::Inflate("CRAM container header CRC32 failure"));
        }
        Ok(Some(ContainerHeader {
            length,
            ref_seq_id,
            ref_seq_start,
            ref_seq_span,
            num_records,
            record_counter,
            num_bases,
            num_blocks,
            landmarks,
            header_len: covered + 4,
        }))
    }

    /// The 1-based last position, `start + span - 1`.
    pub fn ref_seq_end(&self) -> i64 {
        self.ref_seq_start + self.ref_seq_span - 1
    }
}

/// A block as stored: its header fields and its (compressed) bytes.
#[derive(Debug, Clone, Copy)]
pub(crate) struct Block<'a> {
    pub method: u8,
    pub content_type: u8,
    pub content_id: i32,
    pub data: &'a [u8],
    pub uncomp_size: usize,
    /// The bytes the CRC covers: the header and the data.
    covered: &'a [u8],
    crc: u32,
}

impl<'a> Block<'a> {
    /// Parses the block at `buf[*at..]`, advancing `at` past it.
    pub(crate) fn parse(buf: &'a [u8], at: &mut usize) -> Result<Block<'a>> {
        let what = "CRAM block";
        let start = *at;
        let mut b = buf.get(start..).unwrap_or_default();
        if b.len() < 2 {
            return Err(Error::Truncated {
                what,
                need: 2,
                have: b.len() as u64,
            });
        }
        let method = b[0];
        let content_type = b[1];
        if method > method::TOK3 {
            return Err(Error::invalid(
                what,
                format!("unknown compression method {method}"),
            ));
        }
        b = &b[2..];
        let content_id = itf8(&mut b, what)?;
        let comp_size = itf8(&mut b, what)?;
        let uncomp_size = itf8(&mut b, what)?;
        if comp_size < 0 || uncomp_size < 0 || (method == method::RAW && comp_size != uncomp_size) {
            return Err(Error::invalid(
                what,
                format!("sizes {comp_size} compressed, {uncomp_size} uncompressed"),
            ));
        }
        let head = buf.len() - start - b.len();
        let data_at = start + head;
        let end = data_at + comp_size as usize;
        if end + 4 > buf.len() {
            return Err(Error::Truncated {
                what,
                need: (end + 4 - start) as u64,
                have: (buf.len() - start) as u64,
            });
        }
        let crc = u32::from_le_bytes(buf[end..end + 4].try_into().unwrap());
        *at = end + 4;
        Ok(Block {
            method,
            content_type,
            content_id,
            data: &buf[data_at..end],
            uncomp_size: uncomp_size as usize,
            covered: &buf[start..end],
            crc,
        })
    }

    /// Fails unless the block's CRC32 matches its bytes.
    pub(crate) fn check_crc(&self) -> Result<()> {
        if libdeflater::crc32(self.covered) != self.crc {
            return Err(Error::Inflate("CRAM block CRC32 failure"));
        }
        Ok(())
    }

    /// The uncompressed bytes: the block's own for raw (and empty) blocks,
    /// otherwise decompressed into `scratch`. `crc` checks the CRC32 first.
    pub(crate) fn uncompress<'s>(&self, scratch: &'s mut Vec<u8>, crc: bool) -> Result<&'s [u8]>
    where
        'a: 's,
    {
        if crc {
            self.check_crc()?;
        }
        if self.uncomp_size == 0 {
            return Ok(&[]);
        }
        if self.method == method::RAW {
            return Ok(self.data);
        }
        decompress(self.method, self.data, self.uncomp_size, scratch)?;
        Ok(&scratch[..])
    }
}

// libbz2's one-shot decoder, which bzip2-sys does not declare.
unsafe extern "C" {
    fn BZ2_bzBuffToBuffDecompress(
        dest: *mut c_char,
        dest_len: *mut c_uint,
        source: *mut c_char,
        source_len: c_uint,
        small: c_int,
        verbosity: c_int,
    ) -> c_int;
}

fn failed(method: u8) -> Error {
    Error::Inflate(match method {
        method::GZIP => "CRAM gzip block does not inflate",
        method::BZIP2 => "CRAM bzip2 block does not decompress",
        method::LZMA => "CRAM lzma block does not decompress",
        method::RANS => "CRAM rANS 4x8 block does not decode",
        method::RANS_NX16 => "CRAM rANS Nx16 block does not decode",
        method::ARITH => "CRAM arithmetic-coded block does not decode",
        method::FQZ => "CRAM fqzcomp block does not decode",
        _ => "CRAM name tokeniser block does not decode",
    })
}

/// htscodecs' big-endian base-128 uint32, checked before C pointer arithmetic.
fn codec_uint(data: &[u8], at: &mut usize, method: u8) -> Result<usize> {
    let mut value = 0u32;
    for _ in 0..5 {
        let b = *data.get(*at).ok_or_else(|| failed(method))?;
        *at += 1;
        value = value
            .checked_mul(128)
            .and_then(|v| v.checked_add((b & 127) as u32))
            .ok_or_else(|| failed(method))?;
        if b & 128 == 0 {
            return Ok(value as usize);
        }
    }
    Err(failed(method))
}

fn entropy_size(data: &[u8], method: u8) -> Result<Option<usize>> {
    let flags = *data.first().ok_or_else(|| failed(method))?;
    // STRIPE has its own length even when NOSZ is also set. Ordinary NOSZ
    // inherits the caller's capacity; the _to decoders enforce that bound.
    if flags & 8 == 0 && flags & 16 != 0 {
        return Ok(None);
    }
    codec_uint(data, &mut 1, method).map(Some)
}

fn entropy_extent(data: &[u8], size: usize, arithmetic: bool, method: u8) -> Result<()> {
    let mut pending = vec![(data, size, 0u8)];
    while let Some((data, size, depth)) = pending.pop() {
        let flags = *data.first().ok_or_else(|| failed(method))?;
        let mut at = 1;
        if flags & 8 != 0 {
            if depth == 32 {
                return Err(Error::Unsupported("CRAM entropy stripe nesting exceeds 32"));
            }
            if codec_uint(data, &mut at, method)? != size {
                return Err(failed(method));
            }
            let count = *data.get(at).ok_or_else(|| failed(method))? as usize;
            at += 1;
            if count == 0 {
                return Err(failed(method));
            }
            let mut lengths = [0usize; 256];
            for length in &mut lengths[..count] {
                *length = codec_uint(data, &mut at, method)?;
                if *length == 0 {
                    return Err(failed(method));
                }
            }
            for (i, &length) in lengths[..count].iter().enumerate() {
                let end = at
                    .checked_add(length)
                    .filter(|&e| e <= data.len())
                    .ok_or_else(|| failed(method))?;
                pending.push((
                    &data[at..end],
                    size / count + usize::from(i < size % count),
                    depth + 1,
                ));
                at = end;
            }
            continue;
        }
        if flags & 16 == 0 && codec_uint(data, &mut at, method)? != size {
            return Err(failed(method));
        }
        let mut packed = size;
        if flags & 128 != 0 {
            let alphabet = *data.get(at).ok_or_else(|| failed(method))? as usize;
            at += if (1..=16).contains(&alphabet) {
                alphabet + 1
            } else {
                1
            };
            if at > data.len() {
                return Err(failed(method));
            }
            packed = codec_uint(data, &mut at, method)?;
            if packed > size {
                return Err(failed(method));
            }
        }
        if !arithmetic && flags & 64 != 0 {
            let metadata = codec_uint(data, &mut at, method)?;
            let literals = codec_uint(data, &mut at, method)?;
            // RLE metadata is a symbol table (at most 256 plus its count)
            // followed by at most one uint32 varint per literal/run.
            if literals > packed
                || metadata / 2 < 2
                || metadata / 2 > literals.saturating_mul(5).saturating_add(257)
            {
                return Err(failed(method));
            }
            let compressed = if metadata & 1 == 0 {
                codec_uint(data, &mut at, method)?
            } else {
                metadata / 2
            };
            if at
                .checked_add(compressed)
                .is_none_or(|end| end > data.len())
            {
                return Err(failed(method));
            }
        }
    }
    Ok(())
}

fn codec_extent(method: u8, data: &[u8], size: usize) -> Result<()> {
    let word = |at| -> Result<usize> {
        let b = data.get(at..at + 4).ok_or_else(|| failed(method))?;
        Ok(u32::from_le_bytes(b.try_into().unwrap()) as usize)
    };
    let stored = match method {
        method::RANS => Some(word(5)?),
        method::FQZ => Some(codec_uint(data, &mut 0, method)?),
        method::TOK3 => Some(word(0)?),
        method::RANS_NX16 | method::ARITH => entropy_size(data, method)?,
        _ => None,
    };
    if stored.is_some_and(|stored| stored != size) {
        return Err(failed(method));
    }
    if matches!(method, method::RANS_NX16 | method::ARITH) {
        entropy_extent(data, size, method == method::ARITH, method)?;
    }
    if method != method::TOK3 {
        return Ok(());
    }
    let reads = word(4)?;
    if reads == 0 || reads > size || data.len() < 9 {
        return Err(failed(method));
    }
    // Name descriptors hold at most one byte or uint32 per record; ALPHA
    // holds strings and one terminator per record. Validate these allocation
    // extents before tok3 allocates its contexts and nested entropy outputs.
    let mut lengths = [usize::MAX; 128 * 16];
    let mut token = None;
    let mut at = 9;
    while at < data.len() {
        let flags = data[at];
        at += 1;
        if flags & 128 != 0 {
            token = Some(token.map_or(0usize, |n| n + 1));
        }
        let token = token.filter(|&n| n < 128).ok_or_else(|| failed(method))?;
        let ty = (flags & 15) as usize;
        let index = token * 16 + ty;
        if flags & 128 != 0 && ty != 0 {
            lengths[token * 16] = reads;
        }
        if flags & 64 != 0 {
            let pair = data.get(at..at + 2).ok_or_else(|| failed(method))?;
            at += 2;
            let from = ((pair[0] as usize) << 4) + pair[1] as usize;
            if from >= index || lengths[from] == usize::MAX {
                return Err(failed(method));
            }
            lengths[index] = lengths[from];
        } else {
            let bytes = codec_uint(data, &mut at, method)?;
            let end = at
                .checked_add(bytes)
                .filter(|&n| n <= data.len())
                .ok_or_else(|| failed(method))?;
            let length = entropy_size(&data[at..end], method)?.ok_or_else(|| failed(method))?;
            let limit = if ty == 1 {
                size.saturating_add(reads)
            } else if matches!(ty, 0 | 2 | 4 | 8 | 9 | 10 | 11 | 12) {
                reads
            } else {
                reads.saturating_mul(4)
            };
            if length > limit {
                return Err(failed(method));
            }
            entropy_extent(&data[at..end], length, data[8] != 0, method)?;
            lengths[index] = length;
            at = end;
        }
    }
    Ok(())
}

/// Copies a `malloc`ed htscodecs result into `out` and frees it.
///
/// # Safety
/// `ptr` is NULL or a uniquely owned malloc allocation. When `len` fits
/// Rust's slice-size bound it holds at least that many initialized bytes;
/// larger extents are rejected and freed without reading the allocation.
unsafe fn take_malloced(ptr: *mut u8, len: usize, out: &mut Vec<u8>) -> bool {
    let Some(ptr) = NonNull::new(ptr) else {
        return false;
    };
    struct Allocation(NonNull<u8>);
    impl Drop for Allocation {
        fn drop(&mut self) {
            // SAFETY: the private owner receives one codec malloc result
            // and frees it exactly once, independently of Rust's allocator.
            unsafe { libc::free(self.0.as_ptr().cast()) };
        }
    }
    let allocation = Allocation(ptr);
    if len > isize::MAX as usize {
        return false;
    }
    out.clear();
    if out.try_reserve_exact(len).is_err() {
        return false;
    }
    // SAFETY: the caller guarantees `len` readable bytes at `ptr`.
    // The extent also fits Rust's slice-size bound, and the C allocation
    // remains owned until copying finishes or unwinds.
    out.extend_from_slice(unsafe { std::slice::from_raw_parts(allocation.0.as_ptr(), len) });
    true
}

/// Decompresses `data` by `method` into `out`, which ends up holding the
/// uncompressed bytes. Methods that record the size must produce exactly
/// `size` bytes, as htslib requires; bzip2, fqzcomp and the name tokeniser
/// keep what they produce.
pub(crate) fn decompress(method: u8, data: &[u8], size: usize, out: &mut Vec<u8>) -> Result<()> {
    if data.is_empty() && method != method::RAW {
        return Err(failed(method));
    }
    codec_extent(method, data, size)?;
    let input = data.as_ptr() as *mut u8;
    let in_len = c_uint::try_from(data.len()).map_err(|_| failed(method))?;
    let exact = |out: &mut Vec<u8>| {
        out.clear();
        out.resize(size, 0);
    };
    match method {
        method::GZIP => {
            // The largest deflate ratio is 1032:1.
            if size / 2048 > data.len() {
                return Err(Error::Inflate(
                    "CRAM gzip block has an impossibly large compression ratio",
                ));
            }
            exact(out);
            let mut d = libdeflater::Decompressor::new();
            let n = if data.starts_with(&[0x1f, 0x8b]) {
                d.gzip_decompress(data, out)
            } else {
                d.zlib_decompress(data, out)
            }
            .map_err(|_| failed(method))?;
            if n != size {
                return Err(failed(method));
            }
        }
        method::BZIP2 => {
            exact(out);
            let mut got = c_uint::try_from(size).map_err(|_| failed(method))?;
            // SAFETY: both buffers are valid for the lengths passed.
            let r = unsafe {
                BZ2_bzBuffToBuffDecompress(
                    out.as_mut_ptr().cast(),
                    &mut got,
                    input.cast(),
                    in_len,
                    0,
                    0,
                )
            };
            if r != 0 {
                return Err(failed(method));
            }
            out.truncate(got as usize);
        }
        method::LZMA => lzma(data, size, out)?,
        method::RANS => {
            let mut got: c_uint = 0;
            // SAFETY: `input` holds `in_len` bytes; the result is malloced.
            let ptr = unsafe { htscodecs_sys::rans_uncompress(input, in_len, &mut got) };
            if got as usize != size {
                // SAFETY: free accepts NULL or the decoder's malloc result.
                unsafe { libc::free(ptr.cast()) };
                return Err(failed(method));
            }
            // SAFETY: on success the buffer holds `got` bytes.
            if !unsafe { take_malloced(ptr, got as usize, out) } || got as usize != size {
                return Err(failed(method));
            }
        }
        method::RANS_NX16 | method::ARITH => {
            exact(out);
            let mut got = c_uint::try_from(size).map_err(|_| failed(method))?;
            // SAFETY: `out` has room for `got` bytes, `input` holds `in_len`.
            let ptr = unsafe {
                if method == method::RANS_NX16 {
                    htscodecs_sys::rans_uncompress_to_4x16(
                        input,
                        in_len,
                        out.as_mut_ptr(),
                        &mut got,
                    )
                } else {
                    htscodecs_sys::arith_uncompress_to(input, in_len, out.as_mut_ptr(), &mut got)
                }
            };
            if ptr.is_null() || got as usize != size {
                return Err(failed(method));
            }
        }
        method::FQZ => {
            let mut got: usize = size;
            // SAFETY: `input` holds `data.len()` bytes; the result is malloced.
            let ptr = unsafe {
                htscodecs_sys::fqz_decompress(
                    input.cast(),
                    data.len(),
                    &mut got,
                    std::ptr::null_mut(),
                    0,
                )
            };
            // SAFETY: on success the buffer holds `got` bytes.
            if !unsafe { take_malloced(ptr.cast(), got, out) } {
                return Err(failed(method));
            }
        }
        method::TOK3 => {
            let mut got: u32 = 0;
            // SAFETY: `input` holds `in_len` bytes; the result is malloced.
            let ptr = unsafe { htscodecs_sys::tok3_decode_names(input, in_len, &mut got) };
            // SAFETY: on success the buffer holds `got` bytes.
            if !unsafe { take_malloced(ptr, got as usize, out) } {
                return Err(failed(method));
            }
        }
        _ => {
            return Err(Error::invalid(
                "CRAM block",
                format!("unknown compression method {method}"),
            ));
        }
    }
    Ok(())
}

/// An xz stream, decoded as htslib's `lzma_mem_inflate` decodes it.
fn lzma(data: &[u8], size: usize, out: &mut Vec<u8>) -> Result<()> {
    use lzma_sys::*;
    out.clear();
    // One guard byte detects a stream larger than its declared extent,
    // while still letting liblzma consume the footer of an exact stream.
    let limit = size.checked_add(1).ok_or_else(|| failed(method::LZMA))?;
    out.resize(limit, 0);
    // SAFETY: a zeroed lzma_stream is LZMA_STREAM_INIT.
    let mut strm: lzma_stream = unsafe { std::mem::zeroed() };
    // SAFETY: plain FFI calls on a stream this function owns and ends; the
    // output pointer always addresses `avail_out` bytes of `out`.
    let (ok, total) = unsafe {
        if lzma_stream_decoder(&mut strm, lzma_easy_decoder_memusage(9), 0) != LZMA_OK {
            return Err(failed(method::LZMA));
        }
        strm.next_in = data.as_ptr();
        strm.avail_in = data.len();
        let ok = loop {
            let pos = strm.total_out as usize;
            if pos > size {
                break false;
            }
            strm.next_out = out.as_mut_ptr().add(pos);
            strm.avail_out = out.len() - pos;
            let finish = strm.avail_in == 0;
            let r = lzma_code(&mut strm, if finish { LZMA_FINISH } else { LZMA_RUN });
            if r == LZMA_STREAM_END {
                break true;
            }
            if r != LZMA_OK || (finish && strm.avail_out != 0) {
                break false;
            }
        };
        let total = strm.total_out as usize;
        lzma_end(&mut strm);
        (ok, total)
    };
    if !ok || total != size {
        return Err(failed(method::LZMA));
    }
    out.truncate(total);
    Ok(())
}

#[cfg(test)]
mod block_safety_tests {
    use super::*;

    #[test]
    fn empty_compressed_payload_is_rejected_without_entering_c() {
        for method in method::GZIP..=method::TOK3 {
            let mut out = vec![9; 8];
            assert!(decompress(method, &[], 8, &mut out).is_err());
        }
    }

    #[test]
    fn codec_size_headers_are_checked_before_foreign_allocations() {
        let huge = 0x7fff_0000u32;
        let mut rans = vec![0; 26];
        rans[1..5].copy_from_slice(&17u32.to_le_bytes());
        rans[5..9].copy_from_slice(&huge.to_le_bytes());
        let mut tok = vec![0; 9];
        tok[..4].copy_from_slice(&huge.to_le_bytes());
        tok[4..8].copy_from_slice(&1u32.to_le_bytes());
        for (method, data) in [(method::RANS, rans), (method::TOK3, tok)] {
            let mut output = vec![7; 4];
            assert!(decompress(method, &data, 4, &mut output).is_err());
            assert_eq!(output, [7; 4]);
        }
        let large_var = [0x87, 0xff, 0xfc, 0x80, 0x00];
        assert!(codec_extent(method::FQZ, &large_var, 4).is_err());
        for method in [method::RANS_NX16, method::ARITH] {
            let mut data = vec![0];
            data.extend(large_var);
            assert!(codec_extent(method, &data, 4).is_err());
            assert!(codec_extent(method, &[0, 128], 4).is_err());
            assert!(codec_extent(method, &[16], 4).is_ok());
        }
        let mut nested = vec![0; 9];
        nested[..4].copy_from_slice(&4u32.to_le_bytes());
        nested[4..8].copy_from_slice(&1u32.to_le_bytes());
        nested.extend([128, 6, 0]);
        nested.extend(large_var);
        assert!(codec_extent(method::TOK3, &nested, 4).is_err());
    }

    #[test]
    fn nested_entropy_metadata_and_recursion_are_bounded() {
        fn uint(out: &mut Vec<u8>, value: usize) {
            if value >= 128 {
                uint(out, value >> 7);
            }
            out.push((value as u8 & 127) | 128);
        }
        fn encode_uint(value: usize) -> Vec<u8> {
            let mut out = Vec::new();
            uint(&mut out, value);
            *out.last_mut().unwrap() &= 127;
            out
        }
        let mut rle = vec![64, 4];
        rle.extend(encode_uint(0x7fff_0000));
        rle.extend([1, 1, 0]);
        assert!(codec_extent(method::RANS_NX16, &rle, 4).is_err());
        let mut data = vec![48, b'A'];
        for level in 0..40 {
            let mut stripe = vec![8, 1, 1];
            stripe.extend(encode_uint(data.len()));
            stripe.extend(data);
            data = stripe;
            if level == 1 {
                let mut out = Vec::new();
                decompress(method::RANS_NX16, &data, 1, &mut out).unwrap();
                assert_eq!(out, b"A");
            }
        }
        assert!(matches!(
            codec_extent(method::RANS_NX16, &data, 1),
            Err(Error::Unsupported(_))
        ));
    }

    #[test]
    fn malloc_result_is_released_on_slice_size_rejection() {
        // SAFETY: malloc returns either NULL or an allocation owned here.
        let ptr: *mut u8 = unsafe { libc::malloc(1) }.cast();
        assert!(!ptr.is_null());
        let mut out = vec![1, 2, 3];
        // SAFETY: the extent is intentionally rejected before a read; the
        // function owns and releases the valid malloc allocation on this path.
        assert!(!unsafe { take_malloced(ptr, usize::MAX, &mut out) });
        assert_eq!(out, [1, 2, 3]);
    }

    #[test]
    fn lzma_output_is_bounded_by_the_declared_extent() {
        let text = vec![b'A'; 64 * 1024];
        // SAFETY: this pure size function receives the actual input extent.
        let bound = unsafe { lzma_sys::lzma_stream_buffer_bound(text.len()) };
        let mut compressed = vec![0; bound];
        let mut count = 0;
        // SAFETY: source and destination are independent initialized buffers
        // with their real extents, and the default allocator is selected.
        let status = unsafe {
            lzma_sys::lzma_easy_buffer_encode(
                1,
                lzma_sys::LZMA_CHECK_CRC64,
                std::ptr::null(),
                text.as_ptr(),
                text.len(),
                compressed.as_mut_ptr(),
                &mut count,
                compressed.len(),
            )
        };
        assert_eq!(status, lzma_sys::LZMA_OK);
        compressed.truncate(count);
        let mut output = Vec::new();
        assert!(lzma(&compressed, 16, &mut output).is_err());
        assert_eq!(output.len(), 17);
        lzma(&compressed, text.len(), &mut output).unwrap();
        assert_eq!(output, text);
    }
}
