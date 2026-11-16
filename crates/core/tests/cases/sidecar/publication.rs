use std::sync::atomic::{AtomicUsize, Ordering};
use std::sync::{Arc, Barrier};

use omicsweave_core::sidecar::{sidecar_path, write_file};

use crate::support::TempDirectory;

#[test]
fn concurrent_sidecar_publishers_never_expose_partial_contents() {
    const THREADS: usize = 8;
    const LEN: usize = 32_768;
    let directory = TempDirectory::new();
    let target = directory.path().join("source.sidecar");
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
    assert_eq!(std::fs::read_dir(directory.path()).unwrap().count(), 1);
}

#[test]
fn failed_publish_cleans_only_its_temporary_file() {
    let directory = TempDirectory::new();
    let target = directory.path().join("source.sidecar");
    std::fs::create_dir(&target).unwrap();
    let other = directory.path().join("source.sidecar.tmp");
    std::fs::write(&other, b"another writer").unwrap();
    assert!(write_file(&target, b"contents").is_err());
    assert!(target.is_dir());
    assert_eq!(std::fs::read(&other).unwrap(), b"another writer");
    assert_eq!(std::fs::read_dir(directory.path()).unwrap().count(), 2);
}

#[test]
fn long_target_name_does_not_lengthen_temporary_name() {
    let directory = TempDirectory::new();
    let target = directory.path().join("x".repeat(250));
    write_file(&target, b"complete").unwrap();
    assert_eq!(std::fs::read(&target).unwrap(), b"complete");
    assert_eq!(std::fs::read_dir(directory.path()).unwrap().count(), 1);
}

#[test]
fn companion_path_preserves_the_full_source_name() {
    assert_eq!(sidecar_path("/a/b.hic").to_str(), Some("/a/b.hic.sidecar"));
}
