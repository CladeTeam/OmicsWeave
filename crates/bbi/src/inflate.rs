//! zlib decompression of data blocks into caller buffers.

use crate::error::{Error, Result, invalid};

const INITIAL_ROOM: usize = 64 * 1024;

fn inflate_error(error: libdeflater::DecompressionError) -> Error {
    match error {
        libdeflater::DecompressionError::BadData => Error::Inflate("corrupt zlib stream"),
        libdeflater::DecompressionError::InsufficientSpace => {
            Error::Inflate("output buffer smaller than the decompressed block")
        }
    }
}

/// Reusable decompressor state. Create one per thread and keep it: the state
/// is not data, but rebuilding it per block costs an allocation.
pub struct Inflater {
    inner: libdeflater::Decompressor,
    last_successful_output: usize,
}

impl Default for Inflater {
    fn default() -> Self {
        Self::new()
    }
}

impl Inflater {
    pub fn new() -> Inflater {
        Inflater {
            inner: libdeflater::Decompressor::new(),
            last_successful_output: 0,
        }
    }

    /// Inflates the zlib stream `raw` into `out` and returns the decompressed
    /// length. `out` must hold [`crate::Header::block_buffer_len`] bytes to
    /// fit any block of the file.
    pub fn inflate(&mut self, raw: &[u8], out: &mut [u8]) -> Result<usize> {
        let n = self
            .inner
            .zlib_decompress(raw, out)
            .map_err(inflate_error)?;
        self.last_successful_output = n;
        Ok(n)
    }

    /// Appends a block at `at`, keeping the file's declaration as a ceiling
    /// without allocating that potentially untrusted ceiling in advance.
    pub(crate) fn inflate_growing(
        &mut self,
        raw: &[u8],
        out: &mut Vec<u8>,
        at: usize,
        limit: usize,
    ) -> Result<usize> {
        let available = out
            .len()
            .checked_sub(at)
            .ok_or_else(|| invalid("BBI inflate buffer", "append position outside buffer"))?;
        // Sequential blocks often need similar room even when the arena's
        // initialized tail is shorter than a full block. Only an observed
        // successful length supplies this hint; the declaration stays a cap.
        let mut room = limit
            .min(INITIAL_ROOM)
            .max(self.last_successful_output.min(limit))
            .max(available.min(limit));
        loop {
            let end = at
                .checked_add(room)
                .ok_or_else(|| invalid("BBI inflate buffer", "length overflows usize"))?;
            if out.len() < end {
                out.try_reserve(end - out.len())
                    .map_err(|e| invalid("BBI inflate buffer", e.to_string()))?;
                out.resize(end, 0);
            }
            // resize initializes every byte before libdeflater receives the
            // mutable slice. Retained prefix bytes are outside that slice.
            match self.inner.zlib_decompress(raw, &mut out[at..end]) {
                Ok(n) => {
                    self.last_successful_output = n;
                    return Ok(n);
                }
                Err(libdeflater::DecompressionError::InsufficientSpace) if room < limit => {
                    room = room.checked_mul(2).unwrap_or(limit).min(limit);
                }
                Err(error) => return Err(inflate_error(error)),
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn compressed(bytes: &[u8]) -> Vec<u8> {
        let mut encoder = libdeflater::Compressor::new(libdeflater::CompressionLvl::default());
        let mut out = vec![0; encoder.zlib_compress_bound(bytes.len())];
        let n = encoder.zlib_compress(bytes, &mut out).unwrap();
        out.truncate(n);
        out
    }

    #[test]
    fn tiny_stream_does_not_allocate_declared_ceiling() {
        let expected = b"small actual block";
        let raw = compressed(expected);
        let mut out = Vec::new();
        let n = Inflater::new()
            .inflate_growing(&raw, &mut out, 0, u32::MAX as usize)
            .unwrap();
        assert_eq!(&out[..n], expected);
        assert_eq!(out.len(), INITIAL_ROOM);
        assert!(out.capacity() <= 2 * INITIAL_ROOM);
    }

    #[test]
    fn large_stream_grows_without_changing_retained_prefix() {
        let expected: Vec<_> = (0..200_000).map(|i| (i % 251) as u8).collect();
        let raw = compressed(&expected);
        let mut out = vec![0xA5; 37];
        let n = Inflater::new()
            .inflate_growing(&raw, &mut out, 37, u32::MAX as usize)
            .unwrap();
        assert_eq!(n, expected.len());
        assert_eq!(&out[..37], &[0xA5; 37]);
        assert_eq!(&out[37..37 + n], expected);
        assert!(out.len() >= 37 + expected.len());
        assert!(out.len() < 37 + 2 * expected.len());
    }

    #[test]
    fn hard_ceiling_and_corruption_errors_allow_reuse() {
        let large = vec![b'a'; 100_000];
        let raw = compressed(&large);
        let mut inflater = Inflater::new();
        let mut out = Vec::new();
        assert!(inflater.inflate_growing(&raw, &mut out, 0, 32_768).is_err());
        assert_eq!(out.len(), 32_768);
        let small = b"recovered output";
        let raw = compressed(small);
        let n = inflater.inflate_growing(&raw, &mut out, 0, 32_768).unwrap();
        assert_eq!(&out[..n], small);
        assert!(
            inflater
                .inflate_growing(b"bad zlib", &mut out, 0, 32_768)
                .is_err()
        );
        let n = inflater.inflate_growing(&raw, &mut out, 0, 32_768).unwrap();
        assert_eq!(&out[..n], small);
        assert_eq!(out.len(), 32_768);
    }

    #[test]
    fn alternating_sizes_keep_appended_bytes_and_recover_after_errors() {
        let mut inflater = Inflater::new();
        let mut out = vec![0xA5; 37];
        let mut expected = out.clone();
        for (block, len) in [80_000, 80_000, 13, 180_000, 11, 80_000]
            .into_iter()
            .enumerate()
        {
            let plain: Vec<_> = (0..len).map(|i| ((i + block * 19) % 251) as u8).collect();
            let raw = compressed(&plain);
            let at = expected.len();
            let n = inflater
                .inflate_growing(&raw, &mut out, at, 500_000)
                .unwrap();
            expected.extend_from_slice(&plain);
            assert_eq!(n, plain.len());
            assert_eq!(&out[..expected.len()], expected);
            assert_eq!(inflater.last_successful_output, plain.len());

            let at = expected.len();
            assert!(
                inflater
                    .inflate_growing(b"bad zlib", &mut out, at, 500_000)
                    .is_err()
            );
            assert_eq!(inflater.last_successful_output, plain.len());
            assert_eq!(&out[..at], expected);
            assert!(
                inflater
                    .inflate_growing(&raw, &mut out, at, plain.len() - 1)
                    .is_err()
            );
            assert_eq!(inflater.last_successful_output, plain.len());
            assert_eq!(&out[..at], expected);
        }
        // The public fixed-buffer entry point shares the same successful
        // length hint and leaves it unchanged when decoding fails.
        let plain = b"fixed buffer output";
        let raw = compressed(plain);
        let mut fixed = [0; 64];
        let n = inflater.inflate(&raw, &mut fixed).unwrap();
        assert_eq!(&fixed[..n], plain);
        assert_eq!(inflater.last_successful_output, plain.len());
        assert!(inflater.inflate(&raw, &mut fixed[..1]).is_err());
        assert_eq!(inflater.last_successful_output, plain.len());
    }
}
