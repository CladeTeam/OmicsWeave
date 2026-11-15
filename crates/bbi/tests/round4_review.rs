//! Independent byte fixtures exercise traversal and decoder recovery.
use omicsweave_bbi::{
    BIGBED_MAGIC, BIGWIG_MAGIC, BbiReader, BlockRef, CHROM_TREE_MAGIC, INDEX_MAGIC, Index,
};
use omicsweave_core::{
    IntervalEntry, IntervalRead, IntervalsOut, Locus, ReadAt, Region, StringColumn, TrackEntry,
    TrackRead,
};
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

#[test]
fn three_level_tree_skips_subtrees_with_nonmonotone_physical_offsets() {
    let base = u64::MAX - 4096;
    let mut bytes = Vec::new();
    index_header(&mut bytes, 4);
    bytes.extend([0, 0, 2, 0]);
    for (start, end, child) in [(0, 50, 296), (100, 150, 244)] {
        bounds(&mut bytes, start, end);
        bytes.extend((base + child).to_le_bytes());
    }
    let spans = [(0, 25), (25, 50), (100, 125), (125, 150)];
    for id in [0, 2, 1, 3] {
        bytes.extend([1, 0, 1, 0]);
        bounds(&mut bytes, spans[id].0, spans[id].1);
        bytes.extend((10_000u64 + id as u64).to_le_bytes());
        bytes.extend(1u64.to_le_bytes());
    }
    for children in [
        [(100, 125, 136), (125, 150, 208)],
        [(25, 50, 172), (0, 25, 100)],
    ] {
        bytes.extend([0, 0, 2, 0]);
        for (start, end, child) in children {
            bounds(&mut bytes, start, end);
            bytes.extend((base + child).to_le_bytes());
        }
    }
    bytes.insert(0, 0xA5);
    let raw = Index::new(&bytes[1..], base).unwrap();
    let mut nodes = Vec::new();
    let ready = raw.prepare(&mut nodes).unwrap();
    for (start, end) in [(0, 200), (110, 140), (1, 2), (60, 90), (25, 25), (0, 200)] {
        let expected: Vec<_> = [1usize, 0, 2, 3]
            .into_iter()
            .filter(|&id| start < end && spans[id].0 < end && spans[id].1 > start)
            .map(|id| BlockRef {
                offset: 10_000 + id as u64,
                size: 1,
            })
            .collect();
        for index in [raw, ready] {
            assert_eq!(
                index
                    .query(0, start, end)
                    .collect::<Result<Vec<_>, _>>()
                    .unwrap(),
                expected
            );
        }
    }
    let mut all = ready.query(0, 0, 200);
    let mut right = ready.query(0, 110, 140);
    assert_eq!(all.next().unwrap().unwrap().offset, 10_001);
    assert_eq!(right.next().unwrap().unwrap().offset, 10_002);
    assert_eq!(all.next().unwrap().unwrap().offset, 10_000);
    assert_eq!(right.next().unwrap().unwrap().offset, 10_003);
    assert!(right.next().is_none());
    assert_eq!(
        all.map(|r| r.unwrap().offset).collect::<Vec<_>>(),
        [10_002, 10_003]
    );
}

#[derive(Clone)]
struct MutableSource(Rc<RefCell<Vec<u8>>>);

impl ReadAt for MutableSource {
    fn len(&self) -> u64 {
        self.0.borrow().len() as u64
    }

    fn read_exact_at(&self, out: &mut [u8], at: u64) -> omicsweave_core::Result<()> {
        self.0.borrow().as_slice().read_exact_at(out, at)
    }
}

fn compressed_container(
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

fn bed(start: u32, note: &[u8]) -> Vec<u8> {
    let mut out = Vec::new();
    u32s(&mut out, &[0, start, start + 2]);
    out.extend(note);
    out.push(0);
    out
}

fn wig(start: u32, count: u16) -> Vec<u8> {
    let mut out = Vec::new();
    u32s(&mut out, &[0, start, start + 2, 0, 0]);
    out.extend([1, 0]);
    out.extend(count.to_le_bytes());
    for _ in 0..count {
        u32s(&mut out, &[start, start + 2, 0x8000_0000]);
    }
    out
}

fn bed_answer(reader: &mut BbiReader<MutableSource>, region: Region<'_>) -> Vec<Vec<u8>> {
    let size = reader.query(region).unwrap();
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

#[test]
fn bigbed_large_declared_ceiling_preserves_prefix_and_recovers_after_corruption() {
    let long = vec![b'x'; 180_000];
    let blocks = [bed(1, &long), bed(101, b"end")];
    let (bytes, ranges) = compressed_container(BIGBED_MAGIC, &blocks, u32::MAX, false);
    let source = MutableSource(Rc::new(RefCell::new(bytes)));
    let mut reader = BbiReader::open(source.clone()).unwrap();
    let full = Region::new("s", 0, 200);
    let short = Region::new("s", 100, 200);
    assert_eq!(
        bed_answer(&mut reader, full),
        [long.clone(), b"end".to_vec()]
    );
    assert_eq!(bed_answer(&mut reader, short), [b"end".to_vec()]);
    let original = source.0.borrow()[ranges[1].clone()].to_vec();
    source.0.borrow_mut()[ranges[1].clone()].fill(0);
    assert!(reader.query(full).is_err());
    assert!(
        reader
            .write(IntervalsOut {
                entries: &mut [],
                annotations: StringColumn {
                    offsets: &mut [0],
                    data: &mut []
                },
            })
            .is_err()
    );
    source.0.borrow_mut()[ranges[1].clone()].copy_from_slice(&original);
    assert_eq!(bed_answer(&mut reader, full), [long, b"end".to_vec()]);
}

#[test]
fn bigwig_geometric_growth_declared_ceiling_and_truncated_stream() {
    let blocks = [wig(1, 10_000), wig(101, 1)];
    for limit in [64, u32::MAX] {
        let (bytes, _) = compressed_container(BIGWIG_MAGIC, &blocks, limit, false);
        let mut reader = BbiReader::open(bytes).unwrap();
        let short = Region::new("s", 100, 110);
        assert_eq!(reader.query_entries(short).unwrap(), 1);
        let result = reader.query_entries(Region::new("s", 0, 200));
        if limit == 64 {
            assert!(result.is_err());
            assert!(reader.write_entries(&mut [TrackEntry::default()]).is_err());
        } else {
            assert_eq!(result.unwrap(), 10_001);
            let mut entries = vec![TrackEntry::default(); 10_001];
            reader.write_entries(&mut entries).unwrap();
            assert!(entries.iter().all(|r| r.value.to_bits() == 0x8000_0000));
            let mut values = [9.0; 10];
            reader
                .values(&Locus::from(Region::new("s", 0, 10)), &mut values)
                .unwrap();
            assert_eq!(values[1].to_bits(), 0x8000_0000);
            assert_eq!(values[2].to_bits(), 0x8000_0000);
            assert!(values[0].is_nan() && values[3..].iter().all(|v| v.is_nan()));
        }
        assert_eq!(reader.query_entries(short).unwrap(), 1);
        let mut entry = [TrackEntry::default()];
        reader.write_entries(&mut entry).unwrap();
        assert_eq!(
            (entry[0].start, entry[0].end, entry[0].value.to_bits()),
            (101, 103, 0x8000_0000)
        );
    }
    for magic in [BIGBED_MAGIC, BIGWIG_MAGIC] {
        let block = if magic == BIGBED_MAGIC {
            bed(1, b"note")
        } else {
            wig(1, 1)
        };
        let (bytes, _) = compressed_container(magic, &[block], 256, true);
        let mut reader = BbiReader::open(bytes).unwrap();
        let region = Region::new("s", 0, 100);
        if magic == BIGBED_MAGIC {
            assert!(reader.query(region).is_err());
        } else {
            assert!(reader.query_entries(region).is_err());
        }
    }
}
