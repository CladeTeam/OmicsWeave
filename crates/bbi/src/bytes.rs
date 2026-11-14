//! Little-endian field reads from unaligned caller bytes.

use crate::error::{Error, Result, invalid};

/// Fails unless `buf` holds at least `need` bytes.
#[inline]
pub(crate) fn require(buf: &[u8], need: usize, what: &'static str) -> Result<()> {
    if buf.len() < need {
        return Err(Error::Truncated {
            what,
            need: need as u64,
            have: buf.len() as u64,
        });
    }
    Ok(())
}

/// Bounds-checked reads for metadata and unvalidated structures.
#[inline]
fn array<const N: usize>(buf: &[u8], at: usize) -> [u8; N] {
    let mut out = [0u8; N];
    out.copy_from_slice(&buf[at..at + N]);
    out
}

#[inline]
pub(crate) fn u16_at(buf: &[u8], at: usize) -> u16 {
    u16::from_le_bytes(array(buf, at))
}

#[inline]
pub(crate) fn u32_at(buf: &[u8], at: usize) -> u32 {
    u32::from_le_bytes(array(buf, at))
}

#[inline]
pub(crate) fn u64_at(buf: &[u8], at: usize) -> u64 {
    u64::from_le_bytes(array(buf, at))
}

#[inline]
pub(crate) fn f64_at(buf: &[u8], at: usize) -> f64 {
    f64::from_le_bytes(array(buf, at))
}

/// Reads a field from a record whose complete byte range was validated.
///
/// # Safety
/// `at..at + N` must be within `buf`, without integer overflow.
#[inline]
unsafe fn array_unchecked<const N: usize>(buf: &[u8], at: usize) -> [u8; N] {
    debug_assert!(at <= buf.len() && N <= buf.len() - at);
    // SAFETY: The caller guarantees N readable bytes at `at`. Byte arrays
    // have alignment 1 and every bit pattern is valid.
    unsafe { buf.as_ptr().add(at).cast::<[u8; N]>().read_unaligned() }
}

/// # Safety
/// `at..at + 4` must be within `buf`, without integer overflow.
#[inline]
pub(crate) unsafe fn u32_unchecked(buf: &[u8], at: usize) -> u32 {
    // SAFETY: The caller guarantees the complete four-byte field.
    u32::from_le_bytes(unsafe { array_unchecked(buf, at) })
}

/// # Safety
/// `at..at + 8` must be within `buf`, without integer overflow.
#[inline]
pub(crate) unsafe fn u64_unchecked(buf: &[u8], at: usize) -> u64 {
    // SAFETY: The caller guarantees the complete eight-byte field.
    u64::from_le_bytes(unsafe { array_unchecked(buf, at) })
}

/// # Safety
/// `at..at + 4` must be within `buf`, without integer overflow.
#[inline]
pub(crate) unsafe fn f32_unchecked(buf: &[u8], at: usize) -> f32 {
    // SAFETY: The caller guarantees the complete four-byte field. All f32
    // bit patterns, including NaNs, are valid.
    f32::from_bits(unsafe { u32_unchecked(buf, at) })
}

/// Position of a file offset inside a buffer that starts at `base`, checked so
/// that `len` bytes from there are present.
#[inline]
pub(crate) fn locate(
    base: u64,
    buf: &[u8],
    offset: u64,
    len: usize,
    what: &'static str,
) -> Result<usize> {
    let hi = base
        .checked_add(buf.len() as u64)
        .ok_or_else(|| invalid(what, "buffer file range overflows u64"))?;
    if offset < base || offset.checked_add(len as u64).is_none_or(|end| end > hi) {
        return Err(Error::OutOfRange {
            what,
            offset,
            lo: base,
            hi,
        });
    }
    Ok((offset - base) as usize)
}
