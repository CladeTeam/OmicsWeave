//! The HDF5 filters cooler files use, undone outside libhdf5 so chunks
//! decode on any thread: deflate (zlib), byte shuffle and Fletcher-32.

use libdeflater::{DecompressionError, Decompressor};
use omicsweave_core::{Error, Result};

/// One stage of a dataset's filter pipeline, in pipeline (write) order.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) enum Stage {
    Deflate,
    Shuffle,
    Fletcher32,
    /// Deflate whose input includes an earlier checksum stage.
    DeflateWithChecksums,
}

/// Undoes `pipeline` on the stored bytes of one chunk.
///
/// `raw` holds the chunk as stored and `spare` is a second buffer; stages
/// ping-pong between the two, and the decoded chunk is left in whichever
/// holds it, which the return value names (`false`: `raw`, `true`:
/// `spare`). Stage `i` is skipped when bit `i` of `filter_mask` is set
/// (HDF5 records filters that were skipped for a chunk, e.g. an optional
/// deflate that did not shrink it). `chunk_bytes` is the decoded size;
/// `elem_size` is the shuffle's element size.
pub(crate) fn decode(
    pipeline: &[Stage],
    filter_mask: u32,
    chunk_bytes: usize,
    elem_size: usize,
    raw: &mut Vec<u8>,
    spare: &mut Vec<u8>,
    inflater: &mut Decompressor,
) -> Result<bool> {
    let mut in_spare = false;
    for (i, stage) in pipeline.iter().enumerate().rev() {
        if i < 32 && filter_mask & (1 << i) != 0 {
            continue;
        }
        let (src, dst) = if in_spare {
            (&mut *spare, &mut *raw)
        } else {
            (&mut *raw, &mut *spare)
        };
        match stage {
            Stage::Fletcher32 => {
                fletcher32_strip(src)?;
                continue;
            }
            Stage::Deflate => inflate(src, dst, chunk_bytes, inflater)?,
            Stage::DeflateWithChecksums => inflate(
                src,
                dst,
                checksum_limit(&pipeline[..i], filter_mask, chunk_bytes)?,
                inflater,
            )?,
            Stage::Shuffle => unshuffle(src, dst, elem_size),
        }
        in_spare = !in_spare;
    }
    let out = if in_spare { &*spare } else { &*raw };
    if out.len() < chunk_bytes {
        return Err(Error::Truncated {
            what: "hdf5 chunk",
            need: chunk_bytes as u64,
            have: out.len() as u64,
        });
    }
    Ok(in_spare)
}

/// Checksums before deflate are inside its stream; only those not skipped
/// by this chunk's mask must fit before their reverse stages remove them.
#[cold]
#[inline(never)]
fn checksum_limit(prefix: &[Stage], filter_mask: u32, chunk_bytes: usize) -> Result<usize> {
    let checksums = prefix
        .iter()
        .enumerate()
        .filter(|&(j, stage)| {
            *stage == Stage::Fletcher32 && (j >= 32 || filter_mask & (1 << j) == 0)
        })
        .count();
    checksums
        .checked_mul(4)
        .and_then(|n| chunk_bytes.checked_add(n))
        .filter(|&n| n <= isize::MAX as usize)
        .ok_or_else(|| Error::invalid("hdf5 chunk", "decoded size overflow"))
}

/// Inflates a zlib stream; HDF5's deflate filter writes zlib framing. The
/// decoded stage is bounded by `chunk_bytes`, including any checksums
/// written before deflate and not yet stripped.
fn inflate(
    src: &[u8],
    dst: &mut Vec<u8>,
    chunk_bytes: usize,
    inflater: &mut Decompressor,
) -> Result<()> {
    dst.resize(chunk_bytes, 0);
    match inflater.zlib_decompress(src, dst) {
        Ok(n) => {
            dst.truncate(n);
            Ok(())
        }
        Err(DecompressionError::InsufficientSpace) => {
            Err(Error::Inflate("hdf5 chunk inflates past its chunk size"))
        }
        Err(DecompressionError::BadData) => Err(Error::Inflate("hdf5 chunk: bad zlib data")),
    }
}

/// Reverses HDF5's byte shuffle: the stored chunk holds byte 0 of every
/// element, then byte 1 of every element, and so on; bytes past the last
/// whole element are left in place.
fn unshuffle(src: &[u8], dst: &mut Vec<u8>, size: usize) {
    dst.clear();
    dst.resize(src.len(), 0);
    let n = src.len().checked_div(size).unwrap_or(0);
    if size <= 1 || n <= 1 {
        dst.copy_from_slice(src);
        return;
    }
    for (j, plane) in src[..n * size].chunks_exact(n).enumerate() {
        for (i, &b) in plane.iter().enumerate() {
            dst[i * size + j] = b;
        }
    }
    dst[n * size..].copy_from_slice(&src[n * size..]);
}

/// Checks and removes the trailing Fletcher-32 checksum. Like libhdf5, it
/// also accepts the byte-swapped checksum written by HDF5 before 1.6.3.
fn fletcher32_strip(buf: &mut Vec<u8>) -> Result<()> {
    let Some(n) = buf.len().checked_sub(4) else {
        return Err(Error::Truncated {
            what: "hdf5 fletcher32 chunk",
            need: 4,
            have: buf.len() as u64,
        });
    };
    let stored = u32::from_le_bytes(buf[n..].try_into().unwrap());
    let sum = fletcher32(&buf[..n]);
    let reversed = ((sum & 0x00ff_00ff) << 8) | ((sum >> 8) & 0x00ff_00ff);
    if stored != sum && stored != reversed {
        return Err(Error::Inflate("hdf5 chunk: fletcher32 checksum mismatch"));
    }
    buf.truncate(n);
    Ok(())
}

/// HDF5's Fletcher-32 (`H5_checksum_fletcher32`): big-endian 16-bit words,
/// sums starting from 0, an odd trailing byte taken as a word's high byte.
pub(crate) fn fletcher32(data: &[u8]) -> u32 {
    let (mut sum1, mut sum2) = (0u32, 0u32);
    let mut words = data.chunks_exact(2);
    let tail = words.remainder();
    // 360 words is the most that can be summed before the 32-bit sums
    // could overflow.
    loop {
        let mut summed = false;
        for w in words.by_ref().take(360) {
            sum1 += u32::from(u16::from_be_bytes([w[0], w[1]]));
            sum2 += sum1;
            summed = true;
        }
        if !summed {
            break;
        }
        sum1 = (sum1 & 0xffff) + (sum1 >> 16);
        sum2 = (sum2 & 0xffff) + (sum2 >> 16);
    }
    if let [b] = tail {
        sum1 += u32::from(*b) << 8;
        sum2 += sum1;
        sum1 = (sum1 & 0xffff) + (sum1 >> 16);
        sum2 = (sum2 & 0xffff) + (sum2 >> 16);
    }
    sum1 = (sum1 & 0xffff) + (sum1 >> 16);
    sum2 = (sum2 & 0xffff) + (sum2 >> 16);
    (sum2 << 16) | sum1
}

#[cfg(test)]
mod tests {
    use super::*;
    use libdeflater::{CompressionLvl, Compressor};

    #[test]
    fn nested_checksums_and_skip_masks_restore_exact_bytes() {
        // Both sides of deflate exercise distinct intermediate bounds;
        // the odd payload also exercises Fletcher's trailing high byte.
        let pipeline = [
            Stage::Fletcher32,
            Stage::DeflateWithChecksums,
            Stage::Fletcher32,
        ];
        let expected = [1, 3, 5, 7, 11, 13, 17];
        let mut inflater = Decompressor::new();
        let mut spare = vec![0xff; 100];
        for mask in 0..8 {
            let mut raw = expected.to_vec();
            for (i, stage) in pipeline.iter().enumerate() {
                if mask & (1 << i) != 0 {
                    continue;
                }
                match stage {
                    Stage::Fletcher32 => {
                        raw.extend_from_slice(&fletcher32(&raw).to_le_bytes());
                    }
                    Stage::Deflate | Stage::DeflateWithChecksums => {
                        let mut compressor = Compressor::new(CompressionLvl::default());
                        let mut compressed = vec![0; compressor.zlib_compress_bound(raw.len())];
                        let n = compressor.zlib_compress(&raw, &mut compressed).unwrap();
                        compressed.truncate(n);
                        raw = compressed;
                    }
                    Stage::Shuffle => unreachable!(),
                }
            }
            let in_spare = decode(
                &pipeline,
                mask,
                expected.len(),
                4,
                &mut raw,
                &mut spare,
                &mut inflater,
            )
            .unwrap();
            let actual = if in_spare { &spare } else { &raw };
            assert_eq!(actual, &expected, "mask={mask}");
        }
    }
}
