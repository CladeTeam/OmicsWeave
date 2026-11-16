use omicsweave_bbi::{WigBlock, ZoomBlock};
use proptest::prelude::*;
fn wig(kind: u8, start: u32, step: u32, span: u32, items: &[(u32, u32, u32)]) -> Vec<u8> {
    let mut buf = Vec::new();
    for v in [0, start, u32::MAX, step, span] {
        buf.extend(v.to_le_bytes());
    }
    buf.extend([kind, 0]);
    buf.extend((items.len() as u16).to_le_bytes());
    for &(s, e, bits) in items {
        if kind != 3 {
            buf.extend(s.to_le_bytes());
        }
        if kind == 1 {
            buf.extend(e.to_le_bytes());
        }
        buf.extend(bits.to_le_bytes());
    }
    buf
}

fn read_u32(buf: &[u8], at: usize) -> u32 {
    u32::from_le_bytes(buf[at..at + 4].try_into().unwrap())
}

proptest! {
    #![proptest_config(ProptestConfig { cases: 256, ..ProptestConfig::default() })]

    #[test]
    fn fill_matches_per_base_reference(
        kind in 1u8..=3,
        origin in prop_oneof![0u32..512, (u32::MAX - 2048)..u32::MAX],
        step in 0u32..16, span in 0u32..32,
        items in prop::collection::vec((any::<u32>(), any::<u32>(), any::<u32>()), 0..64),
        query in any::<u32>(), relative in 0u32..1024,
        len in 0usize..256,
    ) {
        let mut buf = wig(kind, origin, step, span, &items);
        // Exercise unaligned input as well as buffers starting at their allocation.
        buf.insert(0, 0xAA);
        let Ok(block) = WigBlock::parse(&buf[1..]) else { return Ok(()); };
        let expected: Vec<_> = items.iter().enumerate().map(|(i, &(s, e, bits))| {
            let (s, e) = match kind {
                1 => (s, e),
                2 => (s, s.checked_add(span).unwrap()),
                _ => {
                    let s = (origin as u64 + i as u64 * step as u64) as u32;
                    (s, s + span)
                }
            };
            (s, e, bits)
        }).collect();
        let decoded: Vec<_> = block.intervals().map(|(s,e,v)| (s,e,v.to_bits())).collect();
        prop_assert_eq!(&decoded, &expected);
        let folded = block.intervals().fold(Vec::new(),|mut out,(s,e,v)|{ out.push((s,e,v.to_bits())); out });
        prop_assert_eq!(&folded,&expected);
        let mut iter = block.intervals();
        if expected.len() > 1 {
            let (s,e,v) = iter.next().unwrap();
            prop_assert_eq!((s,e,v.to_bits()),expected[0]);
            let (s,e,v) = iter.next_back().unwrap();
            prop_assert_eq!((s,e,v.to_bits()),expected[expected.len()-1]);
            prop_assert_eq!(iter.len(),expected.len()-2);
            let middle = iter.fold(Vec::new(),|mut out,(s,e,v)|{ out.push((s,e,v.to_bits())); out });
            prop_assert_eq!(middle,&expected[1..expected.len()-1]);
        }
        prop_assert!(block.intervals().nth(usize::MAX).is_none());
        prop_assert!(block.intervals().nth_back(usize::MAX).is_none());
        for (i, &(s,e,bits)) in expected.iter().enumerate() {
            let (a,b,v) = block.item(i);
            prop_assert_eq!((a,b,v.to_bits()), (s,e,bits));
        }
        for start in [query, origin.saturating_add(relative), origin.saturating_sub(relative)] {
            let mut out = vec![f32::from_bits(0x7FC0_1234); len];
            let mut reference = out.clone();
            let mut writes = 0;
            for &(s,e,bits) in &expected {
                for (j, v) in reference.iter_mut().enumerate() {
                    let pos = start as u64 + j as u64;
                    if s as u64 <= pos && pos < e as u64 {
                        *v = f32::from_bits(bits);
                        writes += 1;
                    }
                }
            }
            prop_assert_eq!(block.fill(0, start, &mut out), writes);
            prop_assert_eq!(out.iter().map(|v| v.to_bits()).collect::<Vec<_>>(),
                reference.iter().map(|v| v.to_bits()).collect::<Vec<_>>());
            let mut prepared_out = vec![f32::from_bits(0x7FC0_1234);len];
            prop_assert_eq!(block.prepare().fill(0,start,&mut prepared_out),writes);
            prop_assert_eq!(prepared_out.iter().map(|v| v.to_bits()).collect::<Vec<_>>(),
                reference.iter().map(|v| v.to_bits()).collect::<Vec<_>>());
            let mut scratch = vec![usize::MAX;3*(len+1)];
            for view in [block,block.prepare()] {
                let mut out = vec![f32::from_bits(0x7FC0_1234);len];
                prop_assert_eq!(view.fill_with_scratch(0,start,&mut out,&mut scratch).unwrap(),writes);
                prop_assert_eq!(out.iter().map(|v| v.to_bits()).collect::<Vec<_>>(),
                    reference.iter().map(|v| v.to_bits()).collect::<Vec<_>>());
                let end = start.saturating_add(len as u32);
                let expected_hits:Vec<_> = expected.iter().copied().filter(|&(s,e,_)| start < end && e > start && s < end).collect();
                let hits:Vec<_> = view.overlapping(0,start,end).map(|(s,e,v)|(s,e,v.to_bits())).collect();
                prop_assert_eq!(hits,expected_hits);
            }
        }
    }
}
proptest! {
    #![proptest_config(ProptestConfig { cases: 256, ..ProptestConfig::default() })]

    #[test]
    fn zoom_decodes_all_fields(buf in prop::collection::vec(any::<u8>(), 0..1024)) {
        let mut padded = vec![0xAA];
        padded.extend_from_slice(&buf[..buf.len() / 32 * 32]);
        let block = ZoomBlock::parse(&padded[1..]).unwrap();
        for (i, r) in block.records().enumerate() {
            let bytes = &padded[1 + i * 32..1 + (i + 1) * 32];
            prop_assert_eq!([r.chrom, r.start, r.end, r.valid, r.min.to_bits(), r.max.to_bits(),
                r.sum.to_bits(), r.sum_squares.to_bits()],
                std::array::from_fn::<_,8,_>(|j| read_u32(bytes,j*4)));
            let direct = block.record(i);
            prop_assert_eq!(direct.sum_squares.to_bits(), r.sum_squares.to_bits());
        }
    }
}
proptest! {
    #![proptest_config(ProptestConfig { cases: 256, ..ProptestConfig::default() })]

    #[test]
    fn prepared_wig_sorted_overlaps_match_stored_order(
        kind in 1u8..=2, mut items in prop::collection::vec((0u32..512,0u32..1024,any::<u32>()),0..80),
        span in 0u32..1024, start in 0u32..1024, len in 0usize..256,
    ) {
        items.sort_by_key(|item|item.0);
        let items:Vec<_> = items.into_iter().map(|(s,w,v)|(s,s+w,v)).collect();
        let buf = wig(kind,0,0,span,&items);
        let raw = WigBlock::parse(&buf).unwrap();
        let view = raw.prepare();
        let mut expected = vec![0.0;len];
        let mut out = expected.clone();
        let count = raw.fill(0,start,&mut expected);
        prop_assert_eq!(view.fill(0,start,&mut out),count);
        prop_assert_eq!(out.iter().map(|v|v.to_bits()).collect::<Vec<_>>(),expected.iter().map(|v|v.to_bits()).collect::<Vec<_>>());
        let mut scratch = vec![0;3*(len+1)];
        out.fill(0.0);
        prop_assert_eq!(view.fill_with_scratch(0,start,&mut out,&mut scratch).unwrap(),count);
        prop_assert_eq!(out.iter().map(|v|v.to_bits()).collect::<Vec<_>>(),expected.iter().map(|v|v.to_bits()).collect::<Vec<_>>());
        let end = start + len as u32;
        let hits:Vec<_> = raw.overlapping(0,start,end).map(|(s,e,v)|(s,e,v.to_bits())).collect();
        prop_assert_eq!(view.overlapping(0,start,end).map(|(s,e,v)|(s,e,v.to_bits())).collect::<Vec<_>>(),hits);
    }
}
proptest! {
    #![proptest_config(ProptestConfig { cases: 256, ..ProptestConfig::default() })]

    #[test]
    fn prepared_zoom_matches_linear_filter(
        mut records in prop::collection::vec((0u32..4,0u32..512,0u32..1024,any::<u32>()),0..80),
        sorted in any::<bool>(),chrom in 0u32..5,start in 0u32..1024,end in 0u32..1024,
    ) {
        if sorted { records.sort_by_key(|r|(r.0,r.1)); }
        let mut buf = vec![0xAA];
        for &(c,s,w,bits) in &records {
            for v in [c,s,s+w,1,bits,bits,bits,bits] { buf.extend(v.to_le_bytes()); }
        }
        let raw = ZoomBlock::parse(&buf[1..]).unwrap();
        let hits:Vec<_> = records.iter().filter(|&&(c,s,w,_)| c == chrom && start < end && s+w > start && s < end)
            .map(|&(c,s,w,v)|(c,s,s+w,v)).collect();
        for view in [raw,raw.prepare()] {
            prop_assert_eq!(view.overlapping(chrom,start,end).map(|r|(r.chrom,r.start,r.end,r.sum.to_bits())).collect::<Vec<_>>(),hits.clone());
        }
    }
}

#[test]
fn fixed_overlap_closed_form_counts_match_all_original_intervals() {
    let n = u16::MAX as usize;
    let items: Vec<_> = (0..n).map(|i| (0, 0, (i as f32).to_bits())).collect();
    for (origin, step, span) in [
        (100, 0, 4096),
        (100, 1, 4096),
        (u32::MAX - 200_000, 2, 4096),
    ] {
        let buf = wig(3, origin, step, span, &items);
        let block = WigBlock::parse(&buf).unwrap();
        for start in [
            origin - 1,
            origin,
            origin + 2048,
            origin + span + step * n as u32 - 10,
            u32::MAX,
        ] {
            let mut expected = [f32::NAN; 64];
            let mut count = 0;
            for (i, _) in items.iter().enumerate() {
                let s = origin as u64 + i as u64 * step as u64;
                let e = s + span as u64;
                let lo = s.max(start as u64);
                let hi = e.min(start as u64 + 64);
                if lo < hi {
                    count += (hi - lo) as usize;
                    expected[(lo - start as u64) as usize..(hi - start as u64) as usize]
                        .fill(i as f32);
                }
            }
            let mut out = [f32::NAN; 64];
            assert_eq!(block.fill(0, start, &mut out), count);
            assert_eq!(out.map(f32::to_bits), expected.map(f32::to_bits));
        }
    }
}

#[test]
fn public_indices_cannot_wrap_or_read_trailing_bytes() {
    let buf = wig(1, 0, 0, 0, &[(10, 20, 1)]);
    let mut padded = buf.clone();
    padded.extend_from_slice(&buf[24..]);
    let block = WigBlock::parse(&padded).unwrap();
    let zoom = ZoomBlock::parse(&[0; 32]).unwrap();
    for i in [1, usize::MAX, usize::MAX / 4 + 1, usize::MAX / 32 + 1] {
        assert!(std::panic::catch_unwind(|| block.item(i)).is_err());
        assert!(std::panic::catch_unwind(|| zoom.record(i)).is_err());
    }
}

#[test]
fn step_coordinates_must_fit_u32() {
    let items = [(u32::MAX, 0, 0); 2];
    assert!(WigBlock::parse(&wig(2, 0, 0, 1, &items)).is_err());
    assert!(WigBlock::parse(&wig(3, u32::MAX, 1, 1, &items)).is_err());
}

#[test]
fn maximum_item_count_and_last_record_are_readable() {
    let items = vec![(7, 8, 1.5f32.to_bits()); u16::MAX as usize];
    for kind in 1..=3 {
        let buf = wig(kind, 7, 0, 1, &items);
        let block = WigBlock::parse(&buf).unwrap();
        assert_eq!(block.item(block.len() - 1), (7, 8, 1.5));
        let mut out = [0.0];
        assert_eq!(block.fill(0, 7, &mut out), items.len());
        assert_eq!(out, [1.5]);
    }
}

#[test]
fn dense_fixed_step_preserves_float_bits_and_output_edges() {
    let items = [0, 0x8000_0000, 0x7FC0_1234, 0x7F80_0000, 0xFF80_0000].map(|bits| (0, 0, bits));
    for origin in [0, 100, u32::MAX - items.len() as u32] {
        let mut buf = vec![0xAA];
        buf.extend(wig(3, origin, 1, 1, &items));
        let block = WigBlock::parse(&buf[1..]).unwrap();
        for start in [origin.saturating_sub(2), origin, origin + 2, u32::MAX] {
            for len in 0..10 {
                let mut out = vec![f32::from_bits(0xDEAD_BEEF); len];
                let mut expected = out.clone();
                let mut count = 0;
                for (i, v) in expected.iter_mut().enumerate() {
                    let pos = start as u64 + i as u64;
                    if pos >= origin as u64 && pos < origin as u64 + items.len() as u64 {
                        *v = f32::from_bits(items[(pos - origin as u64) as usize].2);
                        count += 1;
                    }
                }
                assert_eq!(block.fill(0, start, &mut out), count);
                assert_eq!(
                    out.iter().map(|v| v.to_bits()).collect::<Vec<_>>(),
                    expected.iter().map(|v| v.to_bits()).collect::<Vec<_>>()
                );
            }
        }
    }
}

mod validity {
    use super::wig;
    use omicsweave_bbi::{BedBlock, BedEntry, Error, WigBlock, ZoomBlock};
    use proptest::prelude::*;
    proptest! {
        #![proptest_config(ProptestConfig { cases: 256, ..ProptestConfig::default() })]

        #[test]
        fn arbitrary_blocks_are_bounded(buf in prop::collection::vec(any::<u8>(), 0..1024)) {
            if let Ok(block) = WigBlock::parse(&buf) {
                let mut out = [0.0; 32];
                block.fill(block.chrom(), 0, &mut out);
                prop_assert_eq!(block.intervals().count(), block.len());
            }
            let mut entries = BedBlock::parse(&buf).unwrap().entries();
            let mut count = 0;
            while let Some(entry) = entries.next() {
                count += 1;
                prop_assert!(count <= buf.len() / 13 + 1);
                if entry.is_err() { prop_assert!(entries.next().is_none()); break; }
            }
        }
    }

    #[test]
    fn prepared_storage_is_validated_before_use() {
        // One empty record: the reused vector is cleared, not appended to.
        let mut offsets = vec![7, 7];
        let raw = BedBlock::parse(&[0; 13]).unwrap();
        assert_eq!(raw.prepare(&mut offsets).unwrap().len(), 1);
        // A record cut short is rejected before any unchecked access.
        let raw = BedBlock::parse(&[0; 12]).unwrap();
        assert!(raw.prepare(&mut offsets).is_err());
        let buf = wig(1, 0, 0, 0, &[(0, 8, 1.0f32.to_bits())]);
        let block = WigBlock::parse(&buf).unwrap();
        let mut out = [0.0; 8];
        assert!(
            block
                .fill_with_scratch(0, 0, &mut out, &mut [0; 26])
                .is_err()
        );
        assert_eq!(out, [0.0; 8]);
        assert_eq!(block.fill_with_scratch(1, 0, &mut out, &mut []).unwrap(), 0);
    }

    #[test]
    fn empty_and_reversed_queries_do_not_overlap() {
        let buf = wig(1, 0, 0, 0, &[(10, 20, 1)]);
        let block = WigBlock::parse(&buf).unwrap();
        let entry = BedEntry {
            chrom: 0,
            start: 10,
            end: 20,
            rest: b"",
        };
        for (start, end) in [(15, 15), (16, 14)] {
            assert_eq!(block.overlapping(0, start, end).count(), 0);
            assert!(!entry.overlaps(0, start, end));
        }
        // bigBed zero-width features follow libBigWig's interior-point semantics.
        let point = BedEntry {
            start: 15,
            end: 15,
            ..entry
        };
        assert!(point.overlaps(0, 10, 20));
        assert!(!point.overlaps(0, 15, 15));
    }

    #[test]
    fn wig_block_checks_type_and_item_bytes() {
        let mut b = vec![0u8; 24];
        b[20] = 9;
        assert!(matches!(WigBlock::parse(&b), Err(Error::Invalid { .. })));

        b[20] = 1;
        b[22..24].copy_from_slice(&3u16.to_le_bytes());
        b.extend_from_slice(&[0u8; 12 * 2]);
        assert!(matches!(WigBlock::parse(&b), Err(Error::Truncated { .. })));
    }

    #[test]
    fn fill_ignores_other_chromosomes_and_clips() {
        let mut b = Vec::new();
        for v in [3u32, 100, 110, 0, 0] {
            b.extend_from_slice(&v.to_le_bytes());
        }
        b.extend_from_slice(&[1, 0]);
        b.extend_from_slice(&1u16.to_le_bytes());
        for v in [100u32, 110] {
            b.extend_from_slice(&v.to_le_bytes());
        }
        b.extend_from_slice(&2.5f32.to_le_bytes());
        let block = WigBlock::parse(&b).unwrap();

        let mut out = [f32::NAN; 8];
        assert_eq!(block.fill(4, 104, &mut out), 0);
        assert!(out.iter().all(|v| v.is_nan()));
        assert_eq!(block.fill(3, 106, &mut out), 4);
        assert_eq!(&out[..4], &[2.5; 4]);
        assert!(out[4..].iter().all(|v| v.is_nan()));
    }

    #[test]
    fn zoom_block_requires_whole_records() {
        assert!(ZoomBlock::parse(&[0u8; 33]).is_err());
        assert_eq!(ZoomBlock::parse(&[0u8; 64]).unwrap().len(), 2);
    }
}
