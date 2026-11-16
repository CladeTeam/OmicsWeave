use std::cell::Cell;

use omicsweave_core::sidecar::SourceIdentity;
use omicsweave_core::{ReadAt, Result};

struct CountingSource<'a> {
    bytes: &'a [u8],
    calls: Cell<usize>,
    read_bytes: Cell<usize>,
}

impl ReadAt for CountingSource<'_> {
    fn len(&self) -> u64 {
        self.bytes.len() as u64
    }

    fn read_exact_at(&self, out: &mut [u8], offset: u64) -> Result<()> {
        self.calls.set(self.calls.get() + 1);
        self.read_bytes.set(self.read_bytes.get() + out.len());
        self.bytes.read_exact_at(out, offset)
    }
}

// Independently assembles the specified sample stream and its word hash.
fn sampled_fingerprint(bytes: &[u8]) -> u64 {
    fn word(h: &mut u64, v: u64) {
        *h = (*h ^ v).wrapping_mul(0x9e37_79b9_7f4a_7c15).rotate_left(31);
    }
    let len = bytes.len();
    let mut ranges = vec![(0, len.min(65_536)), (len.saturating_sub(65_536), len)];
    ranges.extend((1..=14).map(|i| {
        let lo = len / 15 * i;
        (lo, (lo + 4096).min(len))
    }));
    let mut h = 0x243f_6a88_85a3_08d3 ^ len as u64;
    for (lo, hi) in ranges {
        word(&mut h, lo as u64);
        let b = &bytes[lo..hi];
        let full = b.len() / 8 * 8;
        for at in (0..full).step_by(8) {
            let mut v = 0;
            for i in 0..8 {
                v |= u64::from(b[at + i]) << (8 * i);
            }
            word(&mut h, v);
        }
        let mut v = ((b.len() - full) as u64) << 56;
        for (i, &byte) in b[full..].iter().enumerate() {
            v ^= u64::from(byte) << (8 * i);
        }
        word(&mut h, v);
    }
    h = (h ^ (h >> 30)).wrapping_mul(0xbf58_476d_1ce4_e5b9);
    h = (h ^ (h >> 27)).wrapping_mul(0x94d0_49bb_1331_11eb);
    h ^ (h >> 31)
}

#[test]
fn fingerprint_preserves_sample_stream_with_one_read_for_small_sources() {
    for len in [
        0, 1, 7, 4095, 4096, 65_535, 65_536, 70_000, 131_072, 131_073, 300_000,
    ] {
        let bytes: Vec<u8> = (0..len).map(|i| (i % 251) as u8).collect();
        let source = CountingSource {
            bytes: &bytes,
            calls: Cell::new(0),
            read_bytes: Cell::new(0),
        };
        let id = SourceIdentity::of(&source).unwrap();
        assert_eq!(id.len, len as u64);
        assert_eq!(id.fingerprint, sampled_fingerprint(&bytes), "len={len}");
        if len <= 131_072 {
            assert_eq!(source.calls.get(), 1, "len={len}");
            assert_eq!(source.read_bytes.get(), len);
        } else {
            assert_eq!(source.calls.get(), 16, "len={len}");
            assert!(source.read_bytes.get() <= 2 * 65_536 + 14 * 4096);
        }
    }
}

#[test]
fn maximal_length_source_keeps_every_sample_in_bounds() {
    struct VirtualSource;
    impl ReadAt for VirtualSource {
        fn len(&self) -> u64 {
            u64::MAX
        }
        fn read_exact_at(&self, out: &mut [u8], offset: u64) -> Result<()> {
            assert!(out.len() <= 65_536);
            assert!(offset.checked_add(out.len() as u64).is_some());
            out.fill(7);
            Ok(())
        }
    }
    assert_eq!(SourceIdentity::of(&VirtualSource).unwrap().len, u64::MAX);
}
