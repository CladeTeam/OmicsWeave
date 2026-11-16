//! Streams: a 32-byte descriptor (codec, decoded length, stored length,
//! checksum) followed by the stored bytes.

use pco::ChunkConfig;
use pco::data_types::Number;
use xxhash_rust::xxh3::Xxh3;

use crate::error::{Error, Result, invalid};
use crate::layout::Cursor;
use crate::options::{allocate, reserve};

pub const DESCRIPTOR_LEN: usize = 32;

/// How a stream's bytes are stored.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub enum Codec {
    /// Little-endian elements as they are.
    Raw,
    /// pcodec standalone format; numeric streams only.
    Pco,
    /// zstd frames; byte streams only.
    Zstd,
}

impl Codec {
    fn id(self) -> u8 {
        match self {
            Codec::Raw => 0,
            Codec::Pco => 1,
            Codec::Zstd => 2,
        }
    }

    fn from_id(id: u8) -> Option<Codec> {
        match id {
            0 => Some(Codec::Raw),
            1 => Some(Codec::Pco),
            2 => Some(Codec::Zstd),
            _ => None,
        }
    }
}

/// Compression effort for the codecs that take one.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct Levels {
    pub pco: usize,
    pub zstd: i32,
}

/// A numeric element stored little-endian.
pub(crate) trait Element: Number {
    const WIDTH: usize;
    fn put(self, out: &mut Vec<u8>);
    fn get(bytes: &[u8]) -> Self;
}

impl Element for u32 {
    const WIDTH: usize = 4;
    fn put(self, out: &mut Vec<u8>) {
        out.extend_from_slice(&self.to_le_bytes());
    }
    fn get(bytes: &[u8]) -> Self {
        u32::from_le_bytes(bytes.try_into().expect("4"))
    }
}

impl Element for f64 {
    const WIDTH: usize = 8;
    fn put(self, out: &mut Vec<u8>) {
        out.extend_from_slice(&self.to_le_bytes());
    }
    fn get(bytes: &[u8]) -> Self {
        f64::from_le_bytes(bytes.try_into().expect("8"))
    }
}

fn push_stream(out: &mut Vec<u8>, codec: Codec, decoded_len: usize, data: &[u8]) -> Result<()> {
    if decoded_len > u32::MAX as usize {
        return Err(invalid("FFM stream", "more than 2^32 - 1 elements"));
    }
    let mut desc = [0u8; DESCRIPTOR_LEN];
    desc[0] = codec.id();
    desc[8..16].copy_from_slice(&(decoded_len as u64).to_le_bytes());
    desc[16..24].copy_from_slice(&(data.len() as u64).to_le_bytes());
    let mut h = Xxh3::new();
    h.update(&desc[..24]);
    h.update(data);
    desc[24..].copy_from_slice(&h.digest().to_le_bytes());
    out.extend_from_slice(&desc);
    out.extend_from_slice(data);
    Ok(())
}

/// Appends a numeric stream.
pub(crate) fn push_numbers<T: Element>(
    out: &mut Vec<u8>,
    values: &[T],
    codec: Codec,
    levels: Levels,
) -> Result<()> {
    #[cfg(target_endian = "little")]
    if codec == Codec::Raw {
        // SAFETY: Element is private and implemented only by u32 and f64,
        // both initialized, padding-free values with WIDTH == size_of::<T>().
        // Their native bytes are the raw little-endian stream. The source
        // slice bounds the byte length and cannot alias mutable out.
        let data = unsafe {
            std::slice::from_raw_parts(values.as_ptr().cast::<u8>(), std::mem::size_of_val(values))
        };
        return push_stream(out, codec, values.len(), data);
    }
    let data = match codec {
        Codec::Raw => {
            let mut b = Vec::with_capacity(values.len() * T::WIDTH);
            for &v in values {
                v.put(&mut b);
            }
            b
        }
        Codec::Pco => {
            let config = ChunkConfig::default().with_compression_level(levels.pco);
            pco::standalone::simple_compress(values, &config)
                .map_err(|e| invalid("FFM stream", format!("pcodec: {e}")))?
        }
        Codec::Zstd => return Err(invalid("FFM stream", "zstd is for byte streams")),
    };
    push_stream(out, codec, values.len(), &data)
}

/// Appends a byte stream.
pub(crate) fn push_bytes(
    out: &mut Vec<u8>,
    bytes: &[u8],
    codec: Codec,
    levels: Levels,
) -> Result<()> {
    if codec == Codec::Raw {
        return push_stream(out, codec, bytes.len(), bytes);
    }
    let data = match codec {
        Codec::Raw => unreachable!("raw stream handled above"),
        Codec::Zstd => zstd::bulk::compress(bytes, levels.zstd)?,
        Codec::Pco => return Err(invalid("FFM stream", "pcodec is for numeric streams")),
    };
    push_stream(out, codec, bytes.len(), &data)
}

/// A checked stream inside a buffer.
pub(crate) struct Stream<'a> {
    pub codec: Codec,
    pub decoded_len: usize,
    pub data: &'a [u8],
}

/// Reads the next stream from `r`, verifying its checksum.
pub(crate) fn next_stream<'a>(r: &mut Cursor<'a>, what: &'static str) -> Result<Stream<'a>> {
    let desc = r.bytes(DESCRIPTOR_LEN)?;
    let codec = Codec::from_id(desc[0])
        .ok_or_else(|| invalid(what, format!("unknown codec {}", desc[0])))?;
    if desc[1..8].iter().any(|&b| b != 0) {
        return Err(invalid(what, "reserved bytes are not zero"));
    }
    let decoded_len = u64::from_le_bytes(desc[8..16].try_into().expect("8"));
    let stored_len = u64::from_le_bytes(desc[16..24].try_into().expect("8"));
    let sum = u64::from_le_bytes(desc[24..32].try_into().expect("8"));
    let stored_len = usize::try_from(stored_len).map_err(|_| invalid(what, "stream too long"))?;
    let data = r.bytes(stored_len)?;
    let mut h = Xxh3::new();
    h.update(&desc[..24]);
    h.update(data);
    if h.digest() != sum {
        return Err(invalid(what, "checksum mismatch"));
    }
    // Every element count is addressed with u32 offsets.
    let decoded_len = u32::try_from(decoded_len)
        .map_err(|_| invalid(what, "stream holds more than 2^32 - 1 elements"))?
        as usize;
    Ok(Stream {
        codec,
        decoded_len,
        data,
    })
}

impl Stream<'_> {
    pub fn numbers<T: Element>(&self, what: &'static str, limit: usize) -> Result<Vec<T>> {
        match self.codec {
            Codec::Raw => {
                if self.decoded_len.checked_mul(T::WIDTH) != Some(self.data.len()) {
                    return Err(invalid(what, "raw length does not match element count"));
                }
                let mut out = reserve(self.decoded_len, limit, what)?;
                out.extend(self.data.chunks_exact(T::WIDTH).map(T::get));
                Ok(out)
            }
            Codec::Pco => {
                let mut out = allocate(self.decoded_len, limit, what)?;
                let progress = pco::standalone::simple_decompress_into(self.data, &mut out)
                    .map_err(|e| invalid(what, format!("pcodec: {e}")))?;
                if !progress.finished || progress.n_processed != self.decoded_len {
                    return Err(invalid(what, "pcodec element count does not match"));
                }
                Ok(out)
            }
            Codec::Zstd => Err(invalid(what, "zstd is not a numeric codec")),
        }
    }

    pub fn bytes(&self, what: &'static str, limit: usize) -> Result<Vec<u8>> {
        if self.codec == Codec::Raw && self.data.len() != self.decoded_len {
            return Err(invalid(what, "raw length does not match element count"));
        }
        if self.codec == Codec::Pco {
            return Err(invalid(what, "pcodec is not a byte codec"));
        }
        let mut out = reserve(self.decoded_len, limit, what)?;
        let len = match self.codec {
            Codec::Raw => {
                out.extend_from_slice(self.data);
                out.len()
            }
            Codec::Zstd => zstd::bulk::Decompressor::new()?
                .decompress_to_buffer(self.data, &mut out)
                .map_err(|_| Error::Inflate("FFM zstd stream"))?,
            Codec::Pco => return Err(invalid(what, "pcodec is not a byte codec")),
        };
        if len != self.decoded_len {
            return Err(invalid(what, "decoded length does not match"));
        }
        Ok(out)
    }
}
