//! Shared positional-read instrumentation for cache and physical I/O assertions.
use std::sync::{Arc, Mutex};

use omicsweave_core::{ReadAt, Result};

#[derive(Clone)]
pub struct Counted<B = Vec<u8>> {
    bytes: Arc<B>,
    calls: Arc<Mutex<Vec<(u64, usize)>>>,
}

impl<B: AsRef<[u8]>> Counted<B> {
    pub fn new(bytes: B) -> Self {
        Self {
            bytes: Arc::new(bytes),
            calls: Arc::new(Mutex::new(Vec::new())),
        }
    }

    pub fn reads(&self) -> usize {
        self.calls.lock().unwrap().len()
    }

    pub fn calls(&self) -> Vec<(u64, usize)> {
        self.calls.lock().unwrap().clone()
    }

    pub fn reset(&self) {
        self.calls.lock().unwrap().clear();
    }
}

impl<B: AsRef<[u8]>> ReadAt for Counted<B> {
    fn len(&self) -> u64 {
        self.bytes.as_ref().as_ref().len() as u64
    }

    fn read_exact_at(&self, out: &mut [u8], offset: u64) -> Result<()> {
        self.calls.lock().unwrap().push((offset, out.len()));
        self.bytes.as_ref().as_ref().read_exact_at(out, offset)
    }
}
