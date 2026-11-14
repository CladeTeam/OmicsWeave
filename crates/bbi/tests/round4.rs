//! Prepared metadata hints must agree with a flat scan for arbitrary layouts.
use omicsweave_bbi::{BlockRef, INDEX_MAGIC, Index};
use proptest::prelude::*;

fn put32(out: &mut Vec<u8>, n: u32) {
    out.extend(n.to_le_bytes());
}
fn put64(out: &mut Vec<u8>, n: u64) {
    out.extend(n.to_le_bytes());
}

proptest! {
    #![proptest_config(ProptestConfig { cases: 256, ..ProptestConfig::default() })]

    #[test]
    fn prepared_nodes_match_flat_scan_in_permuted_physical_order(
        leaves in prop::collection::vec(
            prop::collection::vec((0u32..8, 0u32..1024, 0u32..128), 0..24), 1..16),
        keys in any::<[u8; 16]>(),
        chrom in 0u32..8, start in 0u32..2048, end in 0u32..2048,
        high_base in any::<bool>(),
    ) {
        let base = if high_base { u64::MAX - 65_536 } else { 0 };
        let count = leaves.iter().map(Vec::len).sum::<usize>();
        let mut physical: Vec<_> = (0..leaves.len()).collect();
        physical.sort_by_key(|&i| keys[i]);
        let mut offsets = vec![0; leaves.len()];
        let mut at = base + 48 + 4 + leaves.len() as u64 * 24;
        for &i in &physical {
            offsets[i] = at;
            at += 4 + leaves[i].len() as u64 * 32;
        }
        // Prefix one byte so all encoded integer fields are unaligned.
        let mut bytes = vec![0xAB];
        put32(&mut bytes, INDEX_MAGIC);
        put32(&mut bytes, 32);
        put64(&mut bytes, count as u64);
        for n in [0, 0, 8, 0] { put32(&mut bytes, n); }
        put64(&mut bytes, 0);
        put32(&mut bytes, 1);
        put32(&mut bytes, 0);
        bytes.extend([0, 0]);
        bytes.extend((leaves.len() as u16).to_le_bytes());
        for &offset in &offsets {
            for n in [0, 0, 8, 0] { put32(&mut bytes, n); }
            put64(&mut bytes, offset);
        }
        let mut expected = Vec::new();
        let mut ordinal = 0u64;
        for leaf in &leaves {
            for &(c, s, width) in leaf {
                expected.push((c, s, s + width, BlockRef { offset: ordinal * 16, size: 16 }));
                ordinal += 1;
            }
        }
        let mut first = vec![0u64; leaves.len()];
        for i in 1..leaves.len() { first[i] = first[i - 1] + leaves[i - 1].len() as u64; }
        for &i in &physical {
            bytes.extend([1, 0]);
            bytes.extend((leaves[i].len() as u16).to_le_bytes());
            for (j, &(c, s, width)) in leaves[i].iter().enumerate() {
                for n in [c, s, c, s + width] { put32(&mut bytes, n); }
                put64(&mut bytes, (first[i] + j as u64) * 16);
                put64(&mut bytes, 16);
            }
        }
        let raw = Index::new(&bytes[1..], base).unwrap();
        let mut nodes = Vec::new();
        let prepared = raw.prepare(&mut nodes).unwrap();
        for (c, lo, hi) in [(chrom, start, end), (chrom, 0, 2048), (7, 1536, 1600), (3, 0, 0), (chrom, start, end)] {
            let oracle: Vec<_> = expected.iter()
                .filter(|&&(rc, s, e, _)| lo < hi && rc == c && s < hi && e > lo)
                .map(|e| e.3).collect();
            for index in [raw, prepared] {
                let got = index.query(c, lo, hi).collect::<Result<Vec<_>, _>>().unwrap();
                prop_assert_eq!(&got, &oracle);
            }
        }
        let mut broken = bytes.clone();
        broken[1 + 48 + 4 + 16..1 + 48 + 4 + 24].copy_from_slice(&u64::MAX.to_le_bytes());
        let bad = Index::new(&broken[1..], base).unwrap();
        let mut query = bad.query(chrom, 0, 2048);
        prop_assert!(query.next().unwrap().is_err());
        prop_assert!(query.next().is_none());
        let got = prepared.query(chrom, 0, 2048).collect::<Result<Vec<_>, _>>().unwrap();
        let oracle: Vec<_> = expected.iter().filter(|&&(c, s, e, _)| c == chrom && s < 2048 && e > 0).map(|e| e.3).collect();
        prop_assert_eq!(got, oracle);
    }
}
