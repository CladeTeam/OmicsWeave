//! Streaming decompression of whole gzip files read at once, for formats that
//! are loaded entirely rather than accessed randomly (such as `.bed.gz`).
//!
//! Plain gzip has no block structure, so the stream is inflated front to
//! back and handed on in chunks; the decompressed text is never held whole.
//! Concatenated members and bgzip files (which are concatenated members)
//! decompress the same way.

use std::io::Read;

use omicsweave_core::{Error, Result};

const CHUNK: usize = 1 << 18;

/// The input starts with the gzip magic bytes.
pub fn is_gzip(raw: &[u8]) -> bool {
    raw.starts_with(&[0x1f, 0x8b])
}

/// Passes the text of `raw` to `sink` in chunks: inflated if `raw` is gzip
/// (any number of members), as is otherwise.
pub fn text_chunks(raw: &[u8], mut sink: impl FnMut(&[u8]) -> Result<()>) -> Result<()> {
    if !is_gzip(raw) {
        return if raw.is_empty() { Ok(()) } else { sink(raw) };
    }
    let mut decoder = flate2::bufread::MultiGzDecoder::new(raw);
    let mut chunk = vec![0u8; CHUNK];
    loop {
        let n = decoder
            .read(&mut chunk)
            .map_err(|e| Error::invalid("gzip", e.to_string()))?;
        if n == 0 {
            return Ok(());
        }
        sink(&chunk[..n])?;
    }
}
