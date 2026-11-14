//! Byte sources. Readers fetch the bytes their queries need through
//! [`ReadAt`], so the caller decides where bytes come from: a file, memory, or
//! a loader's cache.

use std::fs::File;
use std::os::unix::fs::FileExt;
use std::path::Path;

use crate::error::{Error, Result};

/// Positional reads of a fixed-length byte source.
pub trait ReadAt {
    /// Total bytes.
    fn len(&self) -> u64;

    fn is_empty(&self) -> bool {
        self.len() == 0
    }

    /// Fills `buf` with the bytes at `offset`; fails if they run past the end.
    fn read_exact_at(&self, buf: &mut [u8], offset: u64) -> Result<()>;
}

impl ReadAt for [u8] {
    fn len(&self) -> u64 {
        <[u8]>::len(self) as u64
    }

    fn read_exact_at(&self, buf: &mut [u8], offset: u64) -> Result<()> {
        let hi = <[u8]>::len(self) as u64;
        let end = offset.checked_add(buf.len() as u64).filter(|&e| e <= hi);
        let Some(end) = end else {
            return Err(Error::OutOfRange {
                what: "read",
                offset,
                lo: 0,
                hi,
            });
        };
        buf.copy_from_slice(&self[offset as usize..end as usize]);
        Ok(())
    }
}

impl ReadAt for Vec<u8> {
    fn len(&self) -> u64 {
        self.as_slice().len() as u64
    }

    fn read_exact_at(&self, buf: &mut [u8], offset: u64) -> Result<()> {
        self.as_slice().read_exact_at(buf, offset)
    }
}

impl<T: ReadAt + ?Sized> ReadAt for &T {
    fn len(&self) -> u64 {
        (**self).len()
    }

    fn read_exact_at(&self, buf: &mut [u8], offset: u64) -> Result<()> {
        (**self).read_exact_at(buf, offset)
    }
}

/// A file read with `pread`; its length is taken once at open.
pub struct FileSource {
    file: File,
    len: u64,
}

impl FileSource {
    pub fn open(path: impl AsRef<Path>) -> Result<FileSource> {
        let file = File::open(path)?;
        let len = file.metadata()?.len();
        Ok(FileSource { file, len })
    }
}

impl ReadAt for FileSource {
    fn len(&self) -> u64 {
        self.len
    }

    fn read_exact_at(&self, buf: &mut [u8], offset: u64) -> Result<()> {
        if offset
            .checked_add(buf.len() as u64)
            .is_none_or(|e| e > self.len)
        {
            return Err(Error::OutOfRange {
                what: "read",
                offset,
                lo: 0,
                hi: self.len,
            });
        }
        Ok(FileExt::read_exact_at(&self.file, buf, offset)?)
    }
}
