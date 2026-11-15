//! Writing destinations: every file under a temporary name, renamed into
//! place once the whole conversion succeeded; and directories for spill files.

use std::fs::File;
use std::io::{self, Write};
use std::path::{Path, PathBuf};
use std::sync::atomic::{AtomicU64, Ordering};

use omicsweave_core::{Error, FileSource, Result};

use super::input::with_path;
use crate::format::Format;

fn partial_name(path: &Path, id: u64) -> PathBuf {
    let mut name = path.as_os_str().to_owned();
    name.push(format!(".partial-{}-{id}", std::process::id()));
    PathBuf::from(name)
}

/// The files a conversion writes: each is created under a unique temporary name and
/// renamed into place by [`Outputs::commit`], in the order created; until
/// then a failure removes them. Each rename is atomic; a set of files is
/// not a filesystem transaction. `F` describes each file (its [`Format`]
/// for conversions).
pub struct Outputs<F = Format> {
    files: Vec<(PathBuf, PathBuf, F)>,
    committed: bool,
}

impl<F> Default for Outputs<F> {
    fn default() -> Outputs<F> {
        Outputs {
            files: Vec::new(),
            committed: false,
        }
    }
}

impl<F: Copy> Outputs<F> {
    pub fn new() -> Outputs<F> {
        Outputs::default()
    }

    /// The temporary file standing for `path` until commit, for writers that
    /// open their destination by name.
    pub fn temporary(&self, path: &Path) -> Option<&Path> {
        self.files
            .iter()
            .find(|(_, p, _)| p == path)
            .map(|(tmp, _, _)| tmp.as_path())
    }

    /// Creates the temporary file for `path`, which is to hold `format`.
    pub fn create(&mut self, path: &Path, format: F) -> Result<File> {
        if path.is_dir() {
            return Err(Error::invalid(
                "output",
                format!("output {} is a directory", path.display()),
            ));
        }
        if let Some(dir) = path.parent().filter(|d| !d.as_os_str().is_empty())
            && !dir.is_dir()
        {
            return Err(Error::Io {
                kind: io::ErrorKind::NotFound,
                detail: format!("output directory {} does not exist", dir.display()),
            });
        }
        if self.files.iter().any(|(_, p, _)| p == path) {
            return Err(Error::invalid("output", "duplicate output path"));
        }
        static NEXT_ID: AtomicU64 = AtomicU64::new(0);
        for _ in 0..1000 {
            let tmp = partial_name(path, NEXT_ID.fetch_add(1, Ordering::Relaxed));
            match File::options().write(true).create_new(true).open(&tmp) {
                Ok(file) => {
                    self.files.push((tmp, path.to_owned(), format));
                    return Ok(file);
                }
                Err(e) if e.kind() == io::ErrorKind::AlreadyExists => continue,
                Err(e) => return Err(with_path(e, &tmp)),
            }
        }
        Err(Error::invalid("output", "cannot create a temporary output"))
    }

    /// Writes `bytes` as the whole content of `path`, which holds `format`.
    pub fn write(&mut self, path: &Path, format: F, bytes: &[u8]) -> Result<()> {
        let mut f = self.create(path, format)?;
        f.write_all(bytes)?;
        Ok(())
    }

    /// Reads a completed temporary output before publishing it, for companions
    /// whose identity must be computed from the new source's bytes.
    pub fn source(&self, path: &Path) -> Result<FileSource> {
        let (tmp, _, _) = self
            .files
            .iter()
            .find(|(_, p, _)| p == path)
            .ok_or_else(|| Error::invalid("output", "source has not been created"))?;
        FileSource::open(tmp)
    }

    /// Syncs and renames every output into place; returns the final paths
    /// with their formats.
    pub fn commit(mut self) -> Result<Vec<(F, PathBuf)>> {
        for (tmp, _, _) in &self.files {
            File::open(tmp)?.sync_all()?;
        }
        for (tmp, path, _) in &self.files {
            std::fs::rename(tmp, path).map_err(|e| with_path(e, path))?;
        }
        self.committed = true;
        Ok(self.files.iter().map(|(_, p, f)| (*f, p.clone())).collect())
    }
}

impl<F> Drop for Outputs<F> {
    fn drop(&mut self) {
        if !self.committed {
            for (tmp, _, _) in &self.files {
                let _ = std::fs::remove_file(tmp);
            }
        }
    }
}

/// A directory for spill files, removed with everything in it when dropped.
pub struct TempDir {
    path: PathBuf,
}

impl TempDir {
    /// A fresh directory inside `parent`.
    pub fn new(parent: &Path) -> Result<TempDir> {
        let stamp = std::time::SystemTime::now()
            .duration_since(std::time::UNIX_EPOCH)
            .map(|d| d.as_nanos())
            .unwrap_or(0);
        for i in 0..1000u32 {
            let path = parent.join(format!(
                ".omicsweave-spill-{}-{stamp}-{i}",
                std::process::id()
            ));
            match std::fs::create_dir(&path) {
                Ok(()) => return Ok(TempDir { path }),
                Err(e) if e.kind() == io::ErrorKind::AlreadyExists => continue,
                Err(e) => return Err(with_path(e, &path)),
            }
        }
        Err(Error::invalid("spill", "cannot create a spill directory"))
    }

    pub fn path(&self) -> &Path {
        &self.path
    }
}

impl Drop for TempDir {
    fn drop(&mut self) {
        let _ = std::fs::remove_dir_all(&self.path);
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn output_instances_own_distinct_files_until_commit() {
        let dir = tempfile::tempdir().unwrap();
        let path = dir.path().join("result");
        std::fs::write(&path, b"original").unwrap();
        let legacy_partial = dir.path().join("result.partial");
        std::fs::write(&legacy_partial, b"unrelated").unwrap();
        let mut first = Outputs::new();
        first.write(&path, Format::Bed, b"first").unwrap();
        let mut second = Outputs::new();
        second.write(&path, Format::Bed, b"second").unwrap();
        assert_ne!(first.files[0].0, second.files[0].0);
        assert!(first.create(&path, Format::Bed).is_err());
        assert_eq!(std::fs::read(&path).unwrap(), b"original");
        drop(first);
        second.commit().unwrap();
        assert_eq!(std::fs::read(&path).unwrap(), b"second");
        assert_eq!(std::fs::read(&legacy_partial).unwrap(), b"unrelated");
        assert_eq!(std::fs::read_dir(dir.path()).unwrap().count(), 2);
    }

    #[test]
    fn failed_commit_removes_unpublished_temporary_files() {
        let dir = tempfile::tempdir().unwrap();
        let path = dir.path().join("directory");
        let mut outputs = Outputs::new();
        outputs
            .write(&path, Format::Bed, b"cannot replace directory")
            .unwrap();
        let tmp = outputs.files[0].0.clone();
        std::fs::create_dir(&path).unwrap();
        assert!(outputs.commit().is_err());
        assert!(!tmp.exists());
        assert!(path.is_dir());
    }
}
