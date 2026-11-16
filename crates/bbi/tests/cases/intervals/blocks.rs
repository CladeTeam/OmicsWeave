use omicsweave_bbi::{BedBlock, BedEntry};
use proptest::prelude::*;

proptest! {
    #![proptest_config(ProptestConfig { cases: 256, ..ProptestConfig::default() })]

    #[test]
    fn prepared_bed_preserves_records_and_query_order(
        mut records in prop::collection::vec((0u32..4,0u32..512,0u32..1024,prop::collection::vec(1u8..=255,0..40)),0..80),
        sorted in any::<bool>(), chrom in 0u32..5, start in 0u32..1024, end in 0u32..1024,
    ) {
        if sorted { records.sort_by_key(|r|(r.0,r.1)); }
        let mut buf = vec![0xAA];
        for (c,s,w,text) in &records {
            for v in [*c,*s,s+w] { buf.extend(v.to_le_bytes()); }
            buf.extend(text); buf.push(0);
        }
        let raw = BedBlock::parse(&buf[1..]).unwrap();
        let mut offsets = vec![u32::MAX; 3];
        let view = raw.prepare(&mut offsets).unwrap();
        prop_assert_eq!(view.len(),records.len());
        let expected:Vec<_> = records.iter().map(|(c,s,w,text)|BedEntry{chrom:*c,start:*s,end:s+w,rest:text}).collect();
        prop_assert_eq!(view.entries().collect::<Vec<_>>(),expected.clone());
        prop_assert_eq!(view.entries().rev().collect::<Vec<_>>(),expected.iter().rev().copied().collect::<Vec<_>>());
        for (i,entry) in expected.iter().enumerate() { prop_assert_eq!(view.item(i),*entry); }
        let hits:Vec<_> = expected.iter().copied().filter(|e|e.overlaps(chrom,start,end)).collect();
        prop_assert_eq!(view.overlapping(chrom,start,end).collect::<Vec<_>>(),hits);
        prop_assert!(std::panic::catch_unwind(|| view.item(usize::MAX)).is_err());
    }
}

#[test]
fn bed_fields_preserve_bytes_across_word_and_simd_boundaries() {
    for len in 0..80 {
        for byte in [1u8, 0x7F, 0x80, 0xFF] {
            let mut buf = vec![0xAA];
            for v in [3u32, 10, 20] {
                buf.extend(v.to_le_bytes());
            }
            let text = vec![byte; len];
            buf.extend_from_slice(&text);
            buf.push(0);
            // Adjacent 0x01 bytes exercise borrow propagation after NUL.
            buf.extend([1; 12]);
            buf.push(0);
            let mut entries = BedBlock::parse(&buf[1..]).unwrap().entries();
            let entry = entries.next().unwrap().unwrap();
            assert_eq!((entry.chrom, entry.start, entry.end), (3, 10, 20));
            assert_eq!(entry.rest, text);
            assert!(entries.next().unwrap().unwrap().rest.is_empty());
            assert!(entries.next().is_none());
        }
    }
}

#[test]
fn bed_entries_stop_at_a_malformed_record() {
    let mut b = Vec::new();
    for v in [0u32, 1, 2] {
        b.extend_from_slice(&v.to_le_bytes());
    }
    b.extend_from_slice(b"ok\0");
    for v in [0u32, 3, 4] {
        b.extend_from_slice(&v.to_le_bytes());
    }
    b.extend_from_slice(b"unterminated");
    let entries: Vec<_> = BedBlock::parse(&b).unwrap().entries().collect();
    assert_eq!(entries.len(), 2);
    assert_eq!(entries[0].as_ref().unwrap().rest, b"ok");
    assert!(entries[1].is_err());
}
