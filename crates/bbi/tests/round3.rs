//! Adversarial node layouts and high-level order checked against plain records.
use omicsweave_bbi::{
    BIGBED_MAGIC, BIGWIG_MAGIC, BbiReader, CHROM_TREE_MAGIC, Chroms, INDEX_MAGIC, Index,
};
use omicsweave_core::{
    IntervalEntry, IntervalRead, IntervalsOut, Locus, ReadAt, Region, StringColumn, TrackEntry,
    TrackRead,
};
use std::{cell::Cell, rc::Rc};

struct Failable {
    bytes: Vec<u8>,
    fail: Rc<Cell<bool>>,
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

fn put32(out: &mut Vec<u8>, value: u32) {
    out.extend(value.to_le_bytes());
}
fn put64(out: &mut Vec<u8>, value: u64) {
    out.extend(value.to_le_bytes());
}
fn index_header(out: &mut Vec<u8>, count: u64) {
    put32(out, INDEX_MAGIC);
    put32(out, 1024);
    put64(out, count);
    out.extend([0; 32]);
}
fn chrom_header(count: u64) -> Vec<u8> {
    let mut out = Vec::new();
    for value in [CHROM_TREE_MAGIC, 1024, 1, 8] {
        put32(&mut out, value);
    }
    put64(&mut out, count);
    put64(&mut out, 0);
    out
}
fn bounds(out: &mut Vec<u8>) {
    for value in [0, 0, 0, 1000] {
        put32(out, value);
    }
}

/// Minimal container: supplied blocks retain their order and all overlap s:0-1000.
fn container(magic: u32, blocks: &[Vec<u8>], compressed: bool) -> Vec<u8> {
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
fn bed(records: &[(u32, u32, &[u8])]) -> Vec<u8> {
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
fn wig(records: &[(u32, u32, u32)]) -> Vec<u8> {
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

#[test]
fn high_level_bigbed_orders_all_hits_stably_across_blocks() {
    let records: Vec<(u32, u32, &[u8])> = vec![
        (90, 100, b"miss"),
        (20, 30, b"a"),
        (2, 8, b"b"),
        (20, 25, b"c"),
        (1, 3, b"d"),
    ];
    let blocks = [bed(&records[..3]), bed(&records[3..])];
    let mut expected: Vec<_> = records
        .iter()
        .copied()
        .filter(|r| r.0 < 40 && r.1 > 0)
        .collect();
    expected.sort_by_key(|r| r.0);
    for compressed in [false, true] {
        let mut reader = BbiReader::open(container(BIGBED_MAGIC, &blocks, compressed)).unwrap();
        let size = reader.query(Region::new("s", 0, 40)).unwrap();
        let mut entries = vec![IntervalEntry::default(); size.entries];
        let mut offsets = vec![0; size.entries + 1];
        let mut data = vec![0; size.annotation_bytes];
        reader
            .write(IntervalsOut {
                entries: &mut entries,
                annotations: StringColumn {
                    offsets: &mut offsets,
                    data: &mut data,
                },
            })
            .unwrap();
        assert_eq!(entries.len(), expected.len());
        for (row, &(start, end, note)) in expected.iter().enumerate() {
            assert_eq!(
                entries[row],
                IntervalEntry {
                    start,
                    end,
                    annotation: row as u32
                }
            );
            assert_eq!(
                &data[offsets[row] as usize..offsets[row + 1] as usize],
                note
            );
        }
    }
}

#[test]
fn high_level_bigwig_sorts_entries_but_preserves_overwrite_order() {
    let records = [
        (20, 30, 0x8000_0000),
        (2, 8, 0x7FC0_1234),
        (20, 25, 0x4000_0000),
        (1, 3, 0x4040_0000),
    ];
    let blocks = [wig(&records[..2]), wig(&records[2..])];
    for compressed in [false, true] {
        let mut reader = BbiReader::open(container(BIGWIG_MAGIC, &blocks, compressed)).unwrap();
        let region = Region::new("s", 0, 40);
        let mut values = [0.0f32; 40];
        reader.values(&Locus::from(region), &mut values).unwrap();
        let mut expected = [f32::NAN; 40];
        for &(start, end, bits) in &records {
            expected[start as usize..end as usize].fill(f32::from_bits(bits));
        }
        assert!(
            values
                .iter()
                .zip(expected)
                .all(|(a, b)| a.to_bits() == b.to_bits())
        );
        let count = reader.query_entries(region).unwrap();
        let mut entries = vec![TrackEntry::default(); count];
        reader.write_entries(&mut entries).unwrap();
        let mut expected = records.to_vec();
        expected.sort_by_key(|r| r.0);
        assert_eq!(
            entries
                .iter()
                .map(|r| (r.start, r.end, r.value.to_bits()))
                .collect::<Vec<_>>(),
            expected
        );
    }
}

#[test]
fn malformed_bigbed_tail_invalidates_staged_output() {
    let mut bytes = bed(&[(90, 100, b"miss")]);
    bytes.push(0);
    let mut reader = BbiReader::open(container(BIGBED_MAGIC, &[bytes], false)).unwrap();
    assert!(reader.query(Region::new("s", 0, 1)).is_err());
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
}

#[test]
fn failed_reads_revoke_staging_and_allow_recovery() {
    let fail = Rc::new(Cell::new(false));
    let mut track = BbiReader::open(Failable {
        bytes: container(BIGWIG_MAGIC, &[wig(&[(1, 3, 0x4040_0000)])], true),
        fail: fail.clone(),
    })
    .unwrap();
    let region = Region::new("s", 0, 40);
    assert_eq!(track.query_entries(region).unwrap(), 1);
    fail.set(true);
    assert!(track.query_entries(region).is_err());
    assert!(track.write_entries(&mut [TrackEntry::default()]).is_err());
    fail.set(false);
    assert_eq!(track.query_entries(region).unwrap(), 1);
    let mut entries = [TrackEntry::default()];
    track.write_entries(&mut entries).unwrap();
    assert_eq!(
        entries[0],
        TrackEntry {
            start: 1,
            end: 3,
            value: 3.0
        }
    );

    let mut intervals = BbiReader::open(Failable {
        bytes: container(BIGBED_MAGIC, &[bed(&[(1, 3, b"a")])], true),
        fail: fail.clone(),
    })
    .unwrap();
    assert_eq!(intervals.query(region).unwrap().entries, 1);
    fail.set(true);
    assert!(intervals.query(region).is_err());
    assert!(
        intervals
            .write(IntervalsOut {
                entries: &mut [IntervalEntry::default()],
                annotations: StringColumn {
                    offsets: &mut [0; 2],
                    data: &mut [0]
                },
            })
            .is_err()
    );
    fail.set(false);
    assert_eq!(intervals.query(region).unwrap().entries, 1);
    let mut entries = [IntervalEntry::default()];
    let mut offsets = [0; 2];
    let mut data = [0];
    intervals
        .write(IntervalsOut {
            entries: &mut entries,
            annotations: StringColumn {
                offsets: &mut offsets,
                data: &mut data,
            },
        })
        .unwrap();
    assert_eq!(
        entries[0],
        IntervalEntry {
            start: 1,
            end: 3,
            annotation: 0
        }
    );
    assert_eq!(offsets, [0, 1]);
    assert_eq!(&data, b"a");
}

#[test]
fn raw_index_limits_repeated_large_nodes_by_bytes() {
    let count = 512u16;
    let mut bytes = Vec::new();
    index_header(&mut bytes, count as u64);
    bytes.extend([0, 0, 2, 0]);
    for _ in 0..2 {
        bounds(&mut bytes);
        put64(&mut bytes, 100);
    }
    bytes.extend([1, 0]);
    bytes.extend(count.to_le_bytes());
    for i in 0..count {
        bounds(&mut bytes);
        put64(&mut bytes, i as u64);
        put64(&mut bytes, 1);
    }
    let results: Vec<_> = Index::new(&bytes, 0).unwrap().query(0, 0, 1).collect();
    assert_eq!(results.len(), count as usize + 1);
    assert!(results[..count as usize].iter().all(Result::is_ok));
    assert!(results.last().unwrap().is_err());
}

#[test]
fn chromosome_tree_rejects_reused_empty_nodes_and_header_overlap() {
    let mut bytes = chrom_header(0);
    bytes.extend([0, 0, 2, 0]);
    for _ in 0..2 {
        bytes.push(b's');
        put64(&mut bytes, 54);
    }
    bytes.extend([1, 0, 0, 0]);
    assert!(Chroms::parse(&bytes, 0).is_err());
    // An empty leaf encoded in reserved header bytes is not a separate node.
    bytes[37..45].copy_from_slice(&24u64.to_le_bytes());
    assert!(Chroms::parse(&bytes, 0).is_err());
}
