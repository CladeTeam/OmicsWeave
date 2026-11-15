//! The CRAM file definition, container headers and blocks, read as
//! htslib's `cram_read_file_def`, `cram_read_container` and
//! `cram_read_block` read them, and block decompression
//! (`cram_uncompress_block`).

use std::ffi::{c_char, c_int, c_uint};

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

/// Copies a `malloc`ed htscodecs result into `out` and frees it.
///
/// # Safety
/// `ptr` is NULL or a `malloc`ed buffer of at least `len` bytes.
unsafe fn take_malloced(ptr: *mut u8, len: usize, out: &mut Vec<u8>) -> bool {
    if ptr.is_null() {
        return false;
    }
    out.clear();
    // SAFETY: the caller guarantees `len` readable bytes at `ptr`.
    out.extend_from_slice(unsafe { std::slice::from_raw_parts(ptr, len) });
    // SAFETY: `ptr` came from htscodecs' malloc.
    unsafe { libc::free(ptr.cast()) };
    true
}

/// Decompresses `data` by `method` into `out`, which ends up holding the
/// uncompressed bytes. Methods that record the size must produce exactly
/// `size` bytes, as htslib requires; bzip2, fqzcomp and the name tokeniser
/// keep what they produce.
pub(crate) fn decompress(method: u8, data: &[u8], size: usize, out: &mut Vec<u8>) -> Result<()> {
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
    out.resize(size.max(1), 0);
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
            if pos == out.len() {
                out.resize(out.len() * 2 + 32768, 0);
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
