use omicsweave_bbi::{CHROM_TREE_MAGIC, Chroms, Error, INDEX_MAGIC, Index, MAX_DEPTH};
fn index_header(item_count: u64) -> Vec<u8> {
    let mut b = Vec::new();
    b.extend_from_slice(&INDEX_MAGIC.to_le_bytes());
    b.extend_from_slice(&64u32.to_le_bytes());
    b.extend_from_slice(&item_count.to_le_bytes());
    b.extend_from_slice(&[0u8; 16]);
    b.extend_from_slice(&0u64.to_le_bytes());
    b.extend_from_slice(&1u32.to_le_bytes());
    b.extend_from_slice(&0u32.to_le_bytes());
    b
}

#[test]
fn index_cycle_is_bounded() {
    // A branch node whose only child is itself.
    let base = 1000u64;
    let mut b = index_header(1);
    let node = base + b.len() as u64;
    b.extend_from_slice(&[0, 0]);
    b.extend_from_slice(&1u16.to_le_bytes());
    for v in [0u32, 0, 0, 100] {
        b.extend_from_slice(&v.to_le_bytes());
    }
    b.extend_from_slice(&node.to_le_bytes());
    let index = Index::new(&b, base).unwrap();
    let results: Vec<_> = index.query(0, 0, 10).collect();
    assert_eq!(results, vec![Err(Error::TooDeep(MAX_DEPTH))]);
    assert!(matches!(
        index.prepare(&mut Vec::new()),
        Err(Error::TooDeep(MAX_DEPTH))
    ));
}

#[test]
fn index_child_outside_buffer_is_an_error() {
    let base = 0u64;
    let mut b = index_header(1);
    b.extend_from_slice(&[0, 0]);
    b.extend_from_slice(&1u16.to_le_bytes());
    for v in [0u32, 0, 0, 100] {
        b.extend_from_slice(&v.to_le_bytes());
    }
    b.extend_from_slice(&1_000_000u64.to_le_bytes());
    let index = Index::new(&b, base).unwrap();
    let results: Vec<_> = index.query(0, 0, 10).collect();
    assert!(matches!(
        results.as_slice(),
        [Err(Error::OutOfRange { .. })]
    ));
}

#[test]
fn empty_query_yields_nothing() {
    let mut b = index_header(0);
    b.extend_from_slice(&[1, 0, 0, 0]);
    let index = Index::new(&b, 0).unwrap();
    assert_eq!(index.query(0, 5, 5).count(), 0);
    assert_eq!(index.query(0, 0, 10).count(), 0);
}

#[test]
fn chrom_tree_rejects_duplicate_ids() {
    let mut b = Vec::new();
    b.extend_from_slice(&CHROM_TREE_MAGIC.to_le_bytes());
    for v in [8u32, 4, 8] {
        b.extend_from_slice(&v.to_le_bytes());
    }
    b.extend_from_slice(&2u64.to_le_bytes());
    b.extend_from_slice(&0u64.to_le_bytes());
    b.extend_from_slice(&[1, 0]);
    b.extend_from_slice(&2u16.to_le_bytes());
    for name in [b"chrA", b"chrB"] {
        b.extend_from_slice(name);
        b.extend_from_slice(&0u32.to_le_bytes());
        b.extend_from_slice(&10u32.to_le_bytes());
    }
    assert!(matches!(Chroms::parse(&b, 0), Err(Error::Invalid { .. })));
}

#[test]
fn index_and_chrom_file_ranges_cannot_wrap() {
    let mut index = index_header(0);
    index.extend_from_slice(&[1, 0, 0, 0]);
    for base in [u64::MAX, u64::MAX - 47, u64::MAX - 51] {
        assert!(Index::new(&index, base).is_err());
    }
    assert!(Index::new(&index, u64::MAX - index.len() as u64).is_ok());

    let mut chrom = vec![0u8; 36];
    chrom[..4].copy_from_slice(&CHROM_TREE_MAGIC.to_le_bytes());
    chrom[8..12].copy_from_slice(&1u32.to_le_bytes());
    chrom[12..16].copy_from_slice(&8u32.to_le_bytes());
    chrom[32] = 1;
    assert!(Chroms::parse(&chrom, u64::MAX - 35).is_err());
    assert!(Chroms::parse(&chrom, u64::MAX - 36).is_ok());
}

#[test]
fn index_rejects_truncated_nodes_and_overflowing_block_ranges() {
    let mut buf = index_header(1);
    buf.extend_from_slice(&[1, 0, 1, 0]);
    for v in [0u32, 0, 0, 100] {
        buf.extend(v.to_le_bytes());
    }
    buf.extend(u64::MAX.to_le_bytes());
    buf.extend(1u64.to_le_bytes());
    let mut query = Index::new(&buf, 0).unwrap().query(0, 0, 10);
    assert!(matches!(query.next(), Some(Err(Error::Invalid { .. }))));
    assert!(query.next().is_none());
    for len in 52..buf.len() {
        let mut query = Index::new(&buf[..len], 0).unwrap().query(0, 0, 10);
        assert!(matches!(query.next(), Some(Err(Error::OutOfRange { .. }))));
        assert!(query.next().is_none());
    }
}

#[test]
fn chrom_tree_bounds_repeated_empty_subtrees() {
    let base = 0u64;
    let mut buf = vec![0u8; 32];
    buf[..4].copy_from_slice(&CHROM_TREE_MAGIC.to_le_bytes());
    buf[8..12].copy_from_slice(&1u32.to_le_bytes());
    buf[12..16].copy_from_slice(&8u32.to_le_bytes());
    for _ in 0..20 {
        let child = buf.len() as u64 + 22;
        buf.extend([0, 0, 2, 0]);
        for _ in 0..2 {
            buf.push(b'a');
            buf.extend(child.to_le_bytes());
        }
    }
    buf.extend([1, 0, 0, 0]);
    assert!(matches!(
        Chroms::parse(&buf, base),
        Err(Error::Invalid { .. })
    ));
}

#[test]
fn index_bounds_repeated_empty_subtrees() {
    let mut buf = index_header(0);
    for _ in 0..20 {
        let child = buf.len() as u64 + 52;
        buf.extend([0, 0, 2, 0]);
        for _ in 0..2 {
            for v in [0u32, 0, 0, 100] {
                buf.extend(v.to_le_bytes());
            }
            buf.extend(child.to_le_bytes());
        }
    }
    buf.extend([1, 0, 0, 0]);
    let mut query = Index::new(&buf, 0).unwrap().query(0, 0, 10);
    assert!(matches!(query.next(), Some(Err(Error::Invalid { .. }))));
    assert!(query.next().is_none());
}

#[test]
fn unordered_index_entries_preserve_all_hits() {
    let mut buf = index_header(2);
    buf.extend([1, 0, 2, 0]);
    for (s, e, offset) in [(100u32, 200u32, 16u64), (0, 50, 32)] {
        for v in [0, s, 0, e] {
            buf.extend(v.to_le_bytes());
        }
        buf.extend(offset.to_le_bytes());
        buf.extend(16u64.to_le_bytes());
    }
    let raw = Index::new(&buf, 0).unwrap();
    let mut nodes = Vec::new();
    for index in [raw, raw.prepare(&mut nodes).unwrap()] {
        let hits: Vec<_> = index.query(0, 0, 25).collect::<Result<_, _>>().unwrap();
        assert_eq!(
            hits,
            vec![omicsweave_bbi::BlockRef {
                offset: 32,
                size: 16
            }]
        );
    }
}

#[test]
fn prepared_index_rejects_reused_nodes() {
    let mut buf = index_header(0);
    buf.extend([0, 0, 2, 0]);
    for _ in 0..2 {
        for v in [0u32, 0, 0, 100] {
            buf.extend(v.to_le_bytes());
        }
        buf.extend(100u64.to_le_bytes());
    }
    buf.extend([1, 0, 0, 0]);
    let index = Index::new(&buf, 0).unwrap();
    assert!(matches!(
        index.prepare(&mut Vec::new()),
        Err(Error::Invalid { .. })
    ));
}

#[test]
fn prepared_index_rejects_distinct_overlapping_nodes() {
    let mut buf = index_header(1);
    let leaf = 48 + 4 + 2 * 24;
    buf.extend([0, 0, 2, 0]);
    for child in [leaf, leaf + 4] {
        for n in [0u32, 0, 2, 100] {
            buf.extend(n.to_le_bytes());
        }
        buf.extend((child as u64).to_le_bytes());
    }
    buf.extend([1, 0, 1, 0]);
    // Chromosome 1 also encodes a valid empty leaf header at leaf+4.
    for n in [1u32, 0, 1, 100] {
        buf.extend(n.to_le_bytes());
    }
    buf.extend(16u64.to_le_bytes());
    buf.extend(4u64.to_le_bytes());
    let index = Index::new(&buf, 0).unwrap();
    assert!(matches!(
        index.prepare(&mut Vec::new()),
        Err(Error::Invalid { .. })
    ));
}
