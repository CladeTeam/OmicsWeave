//! Short reads and sink failures independent of BGZF internals.
use std::io::{self, Read, Write};
use std::sync::{Arc, Mutex};

pub struct Fragmented<'a> {
    pub input: &'a [u8],
    pub calls: usize,
}

impl Read for Fragmented<'_> {
    fn read(&mut self, out: &mut [u8]) -> io::Result<usize> {
        self.calls += 1;
        if self.calls.is_multiple_of(5) {
            return Err(io::ErrorKind::Interrupted.into());
        }
        let n = out.len().min(3).min(self.input.len());
        out[..n].copy_from_slice(&self.input[..n]);
        self.input = &self.input[n..];
        Ok(n)
    }
}

pub struct FailAfter {
    pub bytes: Arc<Mutex<Vec<u8>>>,
    pub limit: usize,
}

impl Write for FailAfter {
    fn write(&mut self, input: &[u8]) -> io::Result<usize> {
        let mut bytes = self.bytes.lock().unwrap();
        let n = input.len().min(self.limit.saturating_sub(bytes.len()));
        if n == 0 {
            return Err(io::Error::other("injected partial write failure"));
        }
        bytes.extend_from_slice(&input[..n]);
        Ok(n)
    }
    fn flush(&mut self) -> io::Result<()> {
        Ok(())
    }
}
