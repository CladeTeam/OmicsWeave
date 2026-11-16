//! Temporary files shared by publication and isolated executor tests.

use std::path::PathBuf;
use std::sync::atomic::{AtomicUsize, Ordering};

pub struct TempDirectory(PathBuf);

impl TempDirectory {
    pub fn path(&self) -> &std::path::Path {
        &self.0
    }

    pub fn new() -> Self {
        static NEXT: AtomicUsize = AtomicUsize::new(0);
        for _ in 0..128 {
            let path = std::env::temp_dir().join(format!(
                "omicsweave-core-tests-{}-{}",
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
