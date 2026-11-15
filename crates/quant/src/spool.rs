//! Ordered spill of fixed-size records: a stream of results that may
//! outgrow memory is kept in a buffer of bounded size, appended to a
//! temporary file when the buffer fills, and read back in order once.

use std::fs::File;
use std::io::{BufReader, BufWriter, Read, Write};
use std::path::{Path, PathBuf};

use omicsweave_convert::stage::TempDir;
use omicsweave_core::Result;

/// A record of a fixed number of bytes.
pub trait Fixed: Sized + Send + 'static {
    const SIZE: usize;
    fn put(&self, out: &mut [u8]);
    fn get(b: &[u8]) -> Self;
}

/// Records in arrival order, spilled when the buffer holds `limit` bytes.
pub struct Spool<T: Fixed> {
    mem: Vec<T>,
    limit: usize,
    parent: PathBuf,
    dir: Option<TempDir>,
    file: Option<(PathBuf, BufWriter<File>)>,
    len: u64,
    spilled: u64,
}

impl<T: Fixed> Spool<T> {
    /// A spool holding up to `limit` bytes before it spills into a
    /// directory created in `parent`.
    pub fn new(limit: u64, parent: &Path) -> Spool<T> {
        Spool {
            mem: Vec::new(),
            limit: (limit / std::mem::size_of::<T>().max(1) as u64).clamp(1024, 1 << 30) as usize,
            parent: parent.to_owned(),
            dir: None,
            file: None,
            len: 0,
            spilled: 0,
        }
    }

    #[inline]
    pub fn push(&mut self, r: T) -> Result<()> {
        self.mem.push(r);
        self.len += 1;
        if self.mem.len() >= self.limit {
            self.spill()?;
        }
        Ok(())
    }

    fn spill(&mut self) -> Result<()> {
        if self.file.is_none() {
            let dir = TempDir::new(&self.parent)?;
            let path = dir.path().join("spool");
            let f = File::create(&path)?;
            self.dir = Some(dir);
            self.file = Some((path, BufWriter::with_capacity(1 << 20, f)));
        }
        let (_, w) = self.file.as_mut().expect("created above");
        let mut buf = vec![0u8; T::SIZE];
        for r in self.mem.drain(..) {
            r.put(&mut buf);
            w.write_all(&buf)?;
        }
        self.spilled += (self.limit * T::SIZE) as u64;
        Ok(())
    }

    pub fn len(&self) -> u64 {
        self.len
    }

    pub fn is_empty(&self) -> bool {
        self.len == 0
    }

    /// Bytes written to the spill file so far.
    pub fn spilled(&self) -> u64 {
        self.spilled
    }

    /// The records, to be read once in order.
    pub fn finish(mut self) -> Result<Spooled<T>> {
        let file = match self.file.take() {
            Some((path, w)) => {
                w.into_inner().map_err(|e| e.into_error())?;
                Some(BufReader::with_capacity(1 << 20, File::open(path)?))
            }
            None => None,
        };
        Ok(Spooled {
            file,
            mem: std::mem::take(&mut self.mem),
            len: self.len,
            spilled: self.spilled,
            _dir: self.dir.take(),
        })
    }
}

/// The records of a finished [`Spool`].
pub struct Spooled<T: Fixed> {
    file: Option<BufReader<File>>,
    mem: Vec<T>,
    len: u64,
    spilled: u64,
    _dir: Option<TempDir>,
}

impl<T: Fixed> Spooled<T> {
    /// Records in memory only.
    pub fn from_vec(v: Vec<T>) -> Spooled<T> {
        Spooled {
            file: None,
            len: v.len() as u64,
            mem: v,
            spilled: 0,
            _dir: None,
        }
    }

    pub fn len(&self) -> u64 {
        self.len
    }

    pub fn is_empty(&self) -> bool {
        self.len == 0
    }

    pub fn spilled(&self) -> u64 {
        self.spilled
    }

    /// Calls `f` with every record in order.
    pub fn drain(mut self, mut f: impl FnMut(T) -> Result<()>) -> Result<()> {
        if let Some(mut r) = self.file.take() {
            let mut buf = vec![0u8; T::SIZE * 4096];
            loop {
                let mut got = 0;
                while got < buf.len() {
                    let n = r.read(&mut buf[got..])?;
                    if n == 0 {
                        break;
                    }
                    got += n;
                }
                if got % T::SIZE != 0 {
                    return Err(omicsweave_core::Error::invalid(
                        "spool",
                        "truncated spill file",
                    ));
                }
                for c in buf[..got].chunks_exact(T::SIZE) {
                    f(T::get(c))?;
                }
                if got < buf.len() {
                    break;
                }
            }
        }
        for r in std::mem::take(&mut self.mem) {
            f(r)?;
        }
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    struct V(u64);
    impl Fixed for V {
        const SIZE: usize = 8;
        fn put(&self, out: &mut [u8]) {
            out.copy_from_slice(&self.0.to_le_bytes());
        }
        fn get(b: &[u8]) -> V {
            V(u64::from_le_bytes(b.try_into().unwrap()))
        }
    }

    #[test]
    fn spilled_records_come_back_in_order() {
        let dir = tempfile::tempdir().unwrap();
        for limit in [1, 1 << 30] {
            let mut s = Spool::<V>::new(limit, dir.path());
            for i in 0..10_000 {
                s.push(V(i)).unwrap();
            }
            assert_eq!(s.spilled() > 0, limit == 1);
            let mut got = Vec::new();
            s.finish()
                .unwrap()
                .drain(|v| {
                    got.push(v.0);
                    Ok(())
                })
                .unwrap();
            assert_eq!(got, (0..10_000).collect::<Vec<_>>());
        }
        assert_eq!(std::fs::read_dir(dir.path()).unwrap().count(), 0);
    }
}
