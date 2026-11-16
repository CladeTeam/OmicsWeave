#![allow(dead_code)]
use omicsweave_bbi::{BbiReader, CHROM_TREE_MAGIC, INDEX_MAGIC};
use omicsweave_core::{IntervalEntry, IntervalRead, IntervalsOut, ReadAt, Region, StringColumn};
use std::{cell::RefCell, ops::Range, rc::Rc};
fn u32s(out: &mut Vec<u8>, words: &[u32]) {
    for word in words {
        out.extend(word.to_le_bytes());
    }
}

fn index_header(out: &mut Vec<u8>, count: u64) {
    u32s(out, &[INDEX_MAGIC, 4]);
    out.extend(count.to_le_bytes());
    out.extend([0; 32]);
}

fn bounds(out: &mut Vec<u8>, start: u32, end: u32) {
    u32s(out, &[0, start, 0, end]);
}

#[derive(Clone)]
pub(crate) struct MutableSource(pub(crate) Rc<RefCell<Vec<u8>>>);

impl ReadAt for MutableSource {
    fn len(&self) -> u64 {
        self.0.borrow().len() as u64
    }

    fn read_exact_at(&self, out: &mut [u8], at: u64) -> omicsweave_core::Result<()> {
        self.0.borrow().as_slice().read_exact_at(out, at)
    }
}

pub(crate) fn compressed_container(
    magic: u32,
    blocks: &[Vec<u8>],
    limit: u32,
    truncate: bool,
) -> (Vec<u8>, Vec<Range<usize>>) {
    let mut out = vec![0; 64];
    let tree = out.len() as u64;
    u32s(&mut out, &[CHROM_TREE_MAGIC, 4, 1, 8]);
    out.extend(1u64.to_le_bytes());
    out.extend(0u64.to_le_bytes());
    out.extend([1, 0, 1, 0, b's']);
    u32s(&mut out, &[0, 200]);
    let data = out.len() as u64;
    out.extend((blocks.len() as u64).to_le_bytes());
    let mut ranges = Vec::new();
    let mut compressor = libdeflater::Compressor::new(libdeflater::CompressionLvl::default());
    for block in blocks {
        let mut stored = vec![0; compressor.zlib_compress_bound(block.len())];
        let n = compressor.zlib_compress(block, &mut stored).unwrap();
        stored.truncate(n - usize::from(truncate));
        let at = out.len();
        out.extend(stored);
        ranges.push(at..out.len());
    }
    let index = out.len() as u64;
    index_header(&mut out, blocks.len() as u64);
    out.extend([1, 0]);
    out.extend((blocks.len() as u16).to_le_bytes());
    for (i, range) in ranges.iter().enumerate() {
        bounds(&mut out, i as u32 * 100, i as u32 * 100 + 100);
        out.extend((range.start as u64).to_le_bytes());
        out.extend((range.len() as u64).to_le_bytes());
    }
    out[..4].copy_from_slice(&magic.to_le_bytes());
    out[4..6].copy_from_slice(&4u16.to_le_bytes());
    for (at, offset) in [(8, tree), (16, data), (24, index)] {
        out[at..at + 8].copy_from_slice(&offset.to_le_bytes());
    }
    out[52..56].copy_from_slice(&limit.to_le_bytes());
    (out, ranges)
}

pub(crate) fn bed(start: u32, note: &[u8]) -> Vec<u8> {
    let mut out = Vec::new();
    u32s(&mut out, &[0, start, start + 2]);
    out.extend(note);
    out.push(0);
    out
}

pub(crate) fn wig(start: u32, count: u16) -> Vec<u8> {
    let mut out = Vec::new();
    u32s(&mut out, &[0, start, start + 2, 0, 0]);
    out.extend([1, 0]);
    out.extend(count.to_le_bytes());
    for _ in 0..count {
        u32s(&mut out, &[start, start + 2, 0x8000_0000]);
    }
    out
}

pub(crate) fn bed_answer(
    reader: &mut BbiReader<MutableSource>,
    region: Region<'_>,
) -> Vec<Vec<u8>> {
    let size = reader.query(region, false).unwrap();
    let guard = IntervalEntry {
        start: 901,
        end: 902,
        annotation: 903,
    };
    let mut entries = vec![guard; size.entries + 2];
    let mut offsets = vec![-1; size.entries + 3];
    let mut data = vec![0xFF; size.annotation_bytes + 2];
    reader
        .write(IntervalsOut {
            entries: &mut entries[1..size.entries + 1],
            annotations: StringColumn {
                offsets: &mut offsets[1..size.entries + 2],
                data: &mut data[1..size.annotation_bytes + 1],
            },
        })
        .unwrap();
    assert_eq!([entries[0], entries[size.entries + 1]], [guard; 2]);
    assert_eq!([offsets[0], offsets[size.entries + 2]], [-1; 2]);
    assert_eq!([data[0], data[size.annotation_bytes + 1]], [0xFF; 2]);
    (0..size.entries)
        .map(|i| data[1 + offsets[i + 1] as usize..1 + offsets[i + 2] as usize].to_vec())
        .collect()
}
