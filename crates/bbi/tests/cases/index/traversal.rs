use crate::support::byte_fixture::*;
use omicsweave_bbi::{Chroms, Index};
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

mod physical_order {
    use omicsweave_bbi::{BlockRef, INDEX_MAGIC, Index};
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
}
