#![allow(dead_code)]
use omicsweave_bbi::{CHROM_TREE_MAGIC, INDEX_MAGIC};
use omicsweave_core::ReadAt;
use std::{cell::Cell, rc::Rc};
pub(crate) struct Failable {
    pub(crate) bytes: Vec<u8>,
    pub(crate) fail: Rc<Cell<bool>>,
}
impl ReadAt for Failable {
    fn len(&self) -> u64 {
        self.bytes.len() as u64
    }
    fn read_exact_at(&self, out: &mut [u8], at: u64) -> omicsweave_core::Result<()> {
        if self.fail.get() {
            return Err(omicsweave_core::Error::invalid(
                "test source",
                "injected read failure",
            ));
        }
        self.bytes.as_slice().read_exact_at(out, at)
    }
}

pub(crate) fn put32(out: &mut Vec<u8>, value: u32) {
    out.extend(value.to_le_bytes());
}
pub(crate) fn put64(out: &mut Vec<u8>, value: u64) {
    out.extend(value.to_le_bytes());
}
pub(crate) fn index_header(out: &mut Vec<u8>, count: u64) {
    put32(out, INDEX_MAGIC);
    put32(out, 1024);
    put64(out, count);
    out.extend([0; 32]);
}
pub(crate) fn chrom_header(count: u64) -> Vec<u8> {
    let mut out = Vec::new();
    for value in [CHROM_TREE_MAGIC, 1024, 1, 8] {
        put32(&mut out, value);
    }
    put64(&mut out, count);
    put64(&mut out, 0);
    out
}
pub(crate) fn bounds(out: &mut Vec<u8>) {
    for value in [0, 0, 0, 1000] {
        put32(out, value);
    }
}

/// Minimal container: supplied blocks retain their order and all overlap s:0-1000.
pub(crate) fn container(magic: u32, blocks: &[Vec<u8>], compressed: bool) -> Vec<u8> {
    let mut out = vec![0; 64];
    let tree = out.len() as u64;
    out.extend(chrom_header(1));
    out.extend([1, 0, 1, 0, b's']);
    put32(&mut out, 0);
    put32(&mut out, 1000);
    let data = out.len() as u64;
    put64(&mut out, blocks.len() as u64);
    let mut references = Vec::new();
    let mut compressor = libdeflater::Compressor::new(libdeflater::CompressionLvl::default());
    for block in blocks {
        let stored = if compressed {
            let mut stored = vec![0; compressor.zlib_compress_bound(block.len())];
            let n = compressor.zlib_compress(block, &mut stored).unwrap();
            stored.truncate(n);
            stored
        } else {
            block.clone()
        };
        references.push((out.len() as u64, stored.len() as u64));
        out.extend(stored);
    }
    let index = out.len() as u64;
    index_header(&mut out, blocks.len() as u64);
    out.extend([1, 0]);
    out.extend((blocks.len() as u16).to_le_bytes());
    for (offset, size) in references {
        bounds(&mut out);
        put64(&mut out, offset);
        put64(&mut out, size);
    }
    put32(&mut out, magic);
    out[..4].copy_from_slice(&magic.to_le_bytes());
    out[4..6].copy_from_slice(&4u16.to_le_bytes());
    out[8..16].copy_from_slice(&tree.to_le_bytes());
    out[16..24].copy_from_slice(&data.to_le_bytes());
    out[24..32].copy_from_slice(&index.to_le_bytes());
    if compressed {
        let max = blocks.iter().map(Vec::len).max().unwrap() as u32;
        out[52..56].copy_from_slice(&max.to_le_bytes());
    }
    out
}
pub(crate) fn bed(records: &[(u32, u32, &[u8])]) -> Vec<u8> {
    let mut out = Vec::new();
    for &(start, end, note) in records {
        for value in [0, start, end] {
            put32(&mut out, value);
        }
        out.extend(note);
        out.push(0);
    }
    out
}
pub(crate) fn wig(records: &[(u32, u32, u32)]) -> Vec<u8> {
    let mut out = Vec::new();
    for value in [0, 0, 1000, 0, 0] {
        put32(&mut out, value);
    }
    out.extend([1, 0]);
    out.extend((records.len() as u16).to_le_bytes());
    for &(start, end, bits) in records {
        for value in [start, end, bits] {
            put32(&mut out, value);
        }
    }
    out
}
