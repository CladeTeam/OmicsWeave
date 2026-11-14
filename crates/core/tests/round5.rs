use std::cell::Cell;
use std::path::PathBuf;
use std::sync::atomic::{AtomicUsize, Ordering};
use std::sync::{Arc, Barrier};

use omicsweave_core::sidecar::{BodyReader, SourceIdentity, write_file};
use omicsweave_core::{Error, ReadAt, Result};

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
fn oversized_body_read_reports_error_without_advancing() {
    let mut reader = BodyReader::new(&[3, 4, 5], "test");
    assert_eq!(reader.u8().unwrap(), 3);
    assert!(matches!(
        reader.bytes(usize::MAX),
        Err(Error::Truncated { have: 3, .. })
    ));
    assert_eq!(reader.remaining(), 2);
    assert_eq!(reader.bytes(2).unwrap(), &[4, 5]);
    reader.finish().unwrap();
}

struct TempDirectory(PathBuf);

impl TempDirectory {
    fn new() -> Self {
        static NEXT: AtomicUsize = AtomicUsize::new(0);
        for _ in 0..128 {
            let path = std::env::temp_dir().join(format!(
                "omicsweave-core-r5-{}-{}",
                std::process::id(),
                NEXT.fetch_add(1, Ordering::Relaxed)
            ));
            match std::fs::create_dir(&path) {
                Ok(()) => return Self(path),
                Err(error) if error.kind() == std::io::ErrorKind::AlreadyExists => continue,
                Err(error) => panic!("{error}"),
            }
        }
        panic!("cannot reserve test directory");
    }
}

impl Drop for TempDirectory {
    fn drop(&mut self) {
        let _ = std::fs::remove_dir_all(&self.0);
    }
}

#[test]
fn concurrent_sidecar_publishers_never_expose_partial_contents() {
    const THREADS: usize = 8;
    const LEN: usize = 32_768;
    let directory = TempDirectory::new();
    let target = directory.0.join("source.sidecar");
    write_file(&target, &vec![0; LEN]).unwrap();
    let barrier = Arc::new(Barrier::new(THREADS + 1));
    let active = Arc::new(AtomicUsize::new(THREADS));
    struct Running<'a>(&'a AtomicUsize);
    impl Drop for Running<'_> {
        fn drop(&mut self) {
            self.0.fetch_sub(1, Ordering::Release);
        }
    }
    std::thread::scope(|scope| {
        for id in 0..THREADS {
            let (target, barrier, active) = (&target, &barrier, &active);
            scope.spawn(move || {
                let _running = Running(active);
                let bytes = vec![(id + 1) as u8; LEN];
                barrier.wait();
                for _ in 0..12 {
                    write_file(target, &bytes).unwrap();
                }
            });
        }
        barrier.wait();
        while active.load(Ordering::Acquire) != 0 {
            let bytes = std::fs::read(&target).unwrap();
            assert_eq!(bytes.len(), LEN);
            assert!(bytes.iter().all(|&byte| byte == bytes[0]));
        }
    });
    assert_eq!(std::fs::read_dir(&directory.0).unwrap().count(), 1);
}

#[test]
fn failed_publish_cleans_only_its_temporary_file() {
    let directory = TempDirectory::new();
    let target = directory.0.join("source.sidecar");
    std::fs::create_dir(&target).unwrap();
    let other = directory.0.join("source.sidecar.tmp");
    std::fs::write(&other, b"another writer").unwrap();
    assert!(write_file(&target, b"contents").is_err());
    assert!(target.is_dir());
    assert_eq!(std::fs::read(&other).unwrap(), b"another writer");
    assert_eq!(std::fs::read_dir(&directory.0).unwrap().count(), 2);
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

#[test]
fn long_target_name_does_not_lengthen_temporary_name() {
    let directory = TempDirectory::new();
    let target = directory.0.join("x".repeat(250));
    write_file(&target, b"complete").unwrap();
    assert_eq!(std::fs::read(&target).unwrap(), b"complete");
    assert_eq!(std::fs::read_dir(&directory.0).unwrap().count(), 1);
}
