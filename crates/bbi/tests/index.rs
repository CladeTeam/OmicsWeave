//! Query traversal agrees with a flat scan over overlapping, sorted bounds.

use omicsweave_bbi::{BlockRef, INDEX_MAGIC, Index, NodeInfo};
use proptest::prelude::*;

proptest! {
    #![proptest_config(ProptestConfig { cases: 256, ..ProptestConfig::default() })]

    #[test]
    fn leaf_query_matches_scan(
        mut spans in prop::collection::vec((0u32..8, any::<u32>(), 0u32..8, any::<u32>()), 0..128),
        chrom in 0u32..8, start in any::<u32>(), end in any::<u32>(),
        base in prop_oneof![Just(0u64),Just(u64::MAX - 8192),any::<u32>().prop_map(u64::from)],
    ) {
        for (sc,s,ec,e) in &mut spans {
            if (*sc,*s) > (*ec,*e) { std::mem::swap(sc,ec); std::mem::swap(s,e); }
        }
        spans.sort_unstable();
        let mut buf = Vec::new();
        buf.extend(INDEX_MAGIC.to_le_bytes());
        buf.extend(128u32.to_le_bytes());
        buf.extend((spans.len() as u64).to_le_bytes());
        buf.extend([0u8;32]);
        buf.extend([1,0]);
        buf.extend((spans.len() as u16).to_le_bytes());
        for (i,&(sc,s,ec,e)) in spans.iter().enumerate() {
            for v in [sc,s,ec,e] { buf.extend(v.to_le_bytes()); }
            buf.extend((i as u64 * 16).to_le_bytes());
            buf.extend(16u64.to_le_bytes());
        }
        let expected: Vec<_> = spans.iter().enumerate()
            .filter(|&(_, &(sc,s,ec,e))| start < end && (sc,s) < (chrom,end) && (ec,e) > (chrom,start))
            .map(|(i,_)| BlockRef { offset: i as u64 * 16, size:16 }).collect();
        let results: Vec<_> = Index::new(&buf,base).unwrap().query(chrom,start,end)
            .collect::<Result<_,_>>().unwrap();
        prop_assert_eq!(&results,&expected);
        let raw = Index::new(&buf,base).unwrap();
        let mut nodes = vec![NodeInfo::default(); 2];
        let ready = raw.prepare(&mut nodes).unwrap();
        let results:Vec<_> = ready.query(chrom,start,end).collect::<Result<_,_>>().unwrap();
        prop_assert_eq!(results,expected);
        // A complete node is checked before any item can be read unchecked.
        if !spans.is_empty() {
            buf.pop();
            if start < end {
                let mut query = Index::new(&buf,base).unwrap().query(chrom,start,end);
                prop_assert!(query.next().unwrap().is_err());
                prop_assert!(query.next().is_none());
            }
        }
    }
}
