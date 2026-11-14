//! Buffered sequential reads over a [`ReadAt`] source, for the variable-length
//! metadata of a `.hic` file (header, footer, matrix records), which is read
//! once at open and in pieces whose sizes are only known while parsing.

use omicsweave_core::{Error, ReadAt, Result};

const CHUNK: usize = 1 << 16;

pub(crate) struct Cursor<'s, S: ?Sized> {
    source: &'s S,
    buf: Vec<u8>,
    /// File offset of `buf[0]`.
    start: u64,
    at: usize,
}

impl<'s, S: ReadAt + ?Sized> Cursor<'s, S> {
    pub fn new(source: &'s S, position: u64) -> Cursor<'s, S> {
        Cursor {
            source,
            buf: Vec::new(),
            start: position,
            at: 0,
        }
    }

    pub fn position(&self) -> u64 {
        self.start + self.at as u64
    }

    pub fn seek(&mut self, position: u64) {
        let end = self.start + self.buf.len() as u64;
        if (self.start..=end).contains(&position) {
            self.at = (position - self.start) as usize;
        } else {
            self.buf.clear();
            self.start = position;
            self.at = 0;
        }
    }

    /// Moves `n` bytes forward; the new position must lie within the source.
    pub fn skip(&mut self, n: u64) -> Result<()> {
        let position = self
            .position()
            .checked_add(n)
            .filter(|&p| p <= self.source.len())
            .ok_or_else(|| Error::Truncated {
                what: "hic metadata",
                need: n,
                have: self.source.len().saturating_sub(self.position()),
            })?;
        self.seek(position);
        Ok(())
    }

    /// Makes at least `n` bytes available from the current position; fewer
    /// only at the end of the source.
    fn fill(&mut self, n: usize) -> Result<()> {
        if self.buf.len() - self.at >= n {
            return Ok(());
        }
        let position = self.position();
        self.buf.drain(..self.at);
        self.start = position;
        self.at = 0;
        let have = self.buf.len();
        let left = self.source.len().saturating_sub(position + have as u64);
        let want = (n.max(CHUNK) - have).min(left as usize);
        self.buf.resize(have + want, 0);
        self.source
            .read_exact_at(&mut self.buf[have..], position + have as u64)?;
        Ok(())
    }

    pub fn bytes(&mut self, n: usize) -> Result<&[u8]> {
        self.fill(n)?;
        if self.buf.len() - self.at < n {
            return Err(Error::Truncated {
                what: "hic metadata",
                need: n as u64,
                have: (self.buf.len() - self.at) as u64,
            });
        }
        let out = &self.buf[self.at..self.at + n];
        self.at += n;
        Ok(out)
    }

    fn array<const N: usize>(&mut self) -> Result<[u8; N]> {
        Ok(self.bytes(N)?.try_into().unwrap())
    }

    pub fn i32(&mut self) -> Result<i32> {
        Ok(i32::from_le_bytes(self.array()?))
    }

    pub fn i64(&mut self) -> Result<i64> {
        Ok(i64::from_le_bytes(self.array()?))
    }

    pub fn f32(&mut self) -> Result<f32> {
        Ok(f32::from_le_bytes(self.array()?))
    }

    /// A non-negative 32-bit count.
    pub fn count(&mut self, what: &'static str) -> Result<usize> {
        let n = self.i32()?;
        usize::try_from(n).map_err(|_| Error::invalid(what, format!("negative count {n}")))
    }

    /// A count stored as `int` in v8 and `long` in v9.
    pub fn count_v9(&mut self, version: i32, what: &'static str) -> Result<u64> {
        let n = if version > 8 {
            self.i64()?
        } else {
            self.i32()? as i64
        };
        u64::try_from(n).map_err(|_| Error::invalid(what, format!("negative count {n}")))
    }

    pub fn offset(&mut self, what: &'static str) -> Result<u64> {
        let n = self.i64()?;
        u64::try_from(n).map_err(|_| Error::invalid(what, format!("negative offset {n}")))
    }

    /// A NUL-terminated string.
    pub fn string(&mut self) -> Result<String> {
        let mut scanned = 0;
        loop {
            if let Some(nul) = self.buf[self.at + scanned..].iter().position(|&b| b == 0) {
                let len = scanned + nul;
                let s = std::str::from_utf8(&self.buf[self.at..self.at + len])
                    .map_err(|_| Error::invalid("hic string", "not UTF-8"))?
                    .to_owned();
                self.at += len + 1;
                return Ok(s);
            }
            scanned = self.buf.len() - self.at;
            let before = self.buf.len() - self.at;
            self.fill(before + CHUNK)?;
            if self.buf.len() - self.at == before {
                return Err(Error::invalid("hic string", "missing NUL terminator"));
            }
        }
    }
}
