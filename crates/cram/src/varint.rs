//! CRAM's variable-length integers, read as htslib's `safe_itf8_get` and
//! `safe_ltf8_get` read them: a value whose bytes run past the end of its
//! buffer is an error.

use omicsweave_core::{Error, Result};

/// Bytes of an ITF8 value by the high nibble of its first byte.
const ITF8_BYTES: [usize; 16] = [1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 3, 3, 4, 5];

fn truncated(what: &'static str, need: usize, have: usize) -> Error {
    Error::Truncated {
        what,
        need: need as u64,
        have: have as u64,
    }
}

/// Reads an ITF8 integer from the front of `buf`, advancing it.
#[inline]
pub(crate) fn itf8(buf: &mut &[u8], what: &'static str) -> Result<i32> {
    let b = *buf;
    let Some(&first) = b.first() else {
        return Err(truncated(what, 1, 0));
    };
    let n = ITF8_BYTES[(first >> 4) as usize];
    if b.len() < n {
        return Err(truncated(what, n, b.len()));
    }
    let v = match n {
        1 => first as u32,
        2 => ((first as u32) << 8 | b[1] as u32) & 0x3fff,
        3 => ((first as u32) << 16 | (b[1] as u32) << 8 | b[2] as u32) & 0x1f_ffff,
        4 => {
            ((first as u32) << 24 | (b[1] as u32) << 16 | (b[2] as u32) << 8 | b[3] as u32)
                & 0x0fff_ffff
        }
        _ => {
            ((first as u32) & 0x0f) << 28
                | (b[1] as u32) << 20
                | (b[2] as u32) << 12
                | (b[3] as u32) << 4
                | (b[4] as u32) & 0x0f
        }
    };
    *buf = &b[n..];
    Ok(v as i32)
}

/// Bytes of an LTF8 value by its first byte.
#[inline]
fn ltf8_len(first: u8) -> usize {
    match first {
        0x00..=0x7f => 1,
        0x80..=0xbf => 2,
        0xc0..=0xdf => 3,
        0xe0..=0xef => 4,
        0xf0..=0xf7 => 5,
        0xf8..=0xfb => 6,
        0xfc..=0xfd => 7,
        0xfe => 8,
        0xff => 9,
    }
}

/// Reads an LTF8 integer from the front of `buf`, advancing it.
pub(crate) fn ltf8(buf: &mut &[u8], what: &'static str) -> Result<i64> {
    let b = *buf;
    let Some(&first) = b.first() else {
        return Err(truncated(what, 1, 0));
    };
    let n = ltf8_len(first);
    if b.len() < n {
        return Err(truncated(what, n, b.len()));
    }
    let be = |bytes: &[u8]| bytes.iter().fold(0u64, |v, &x| v << 8 | x as u64);
    let v = match n {
        1 => first as u64,
        2..=7 => be(&b[..n]) & ((1u64 << (8 * n - n)) - 1),
        _ => be(&b[1..n]),
    };
    *buf = &b[n..];
    Ok(v as i64)
}

/// Reads a little-endian `int32` from the front of `buf`, advancing it.
pub(crate) fn int32(buf: &mut &[u8], what: &'static str) -> Result<i32> {
    let b = *buf;
    if b.len() < 4 {
        return Err(truncated(what, 4, b.len()));
    }
    *buf = &b[4..];
    Ok(i32::from_le_bytes(b[..4].try_into().unwrap()))
}

#[cfg(test)]
mod tests {
    use super::*;

    fn put_itf8(v: i32) -> Vec<u8> {
        let u = v as u32;
        if u < 0x80 {
            vec![u as u8]
        } else if u < 0x4000 {
            vec![(u >> 8 | 0x80) as u8, u as u8]
        } else if u < 0x20_0000 {
            vec![(u >> 16 | 0xc0) as u8, (u >> 8) as u8, u as u8]
        } else if u < 0x1000_0000 {
            vec![
                (u >> 24 | 0xe0) as u8,
                (u >> 16) as u8,
                (u >> 8) as u8,
                u as u8,
            ]
        } else {
            vec![
                (u >> 28 | 0xf0) as u8,
                (u >> 20) as u8,
                (u >> 12) as u8,
                (u >> 4) as u8,
                (u & 0x0f) as u8,
            ]
        }
    }

    #[test]
    fn itf8_round_trips_every_width() {
        for v in [
            0,
            1,
            127,
            128,
            0x3fff,
            0x4000,
            0x1f_ffff,
            0x20_0000,
            0x0fff_ffff,
            -1,
            i32::MIN,
            i32::MAX,
        ] {
            let bytes = put_itf8(v);
            let mut b = &bytes[..];
            assert_eq!(itf8(&mut b, "t").unwrap(), v);
            assert!(b.is_empty());
            let mut short = &bytes[..bytes.len() - 1];
            if !short.is_empty() {
                assert!(itf8(&mut short, "t").is_err());
            }
        }
    }

    #[test]
    fn ltf8_widths() {
        let mut b: &[u8] = &[0x81, 0x00];
        assert_eq!(ltf8(&mut b, "t").unwrap(), 0x100);
        let mut b: &[u8] = &[0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff];
        assert_eq!(ltf8(&mut b, "t").unwrap(), -1);
        let mut b: &[u8] = &[0xfe, 0, 0, 0, 0, 0, 1, 0];
        assert_eq!(ltf8(&mut b, "t").unwrap(), 0x100);
        let mut b: &[u8] = &[0xe1, 2, 3, 4];
        assert_eq!(ltf8(&mut b, "t").unwrap(), 0x0102_0304);
    }
}
