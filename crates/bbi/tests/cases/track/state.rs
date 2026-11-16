use crate::support::byte_fixture::{container, wig};
use omicsweave_bbi::{BIGWIG_MAGIC, BbiReader};
use omicsweave_core::{Locus, Region, TrackEntry, TrackRead};
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
        let count = reader.query_entries(region, true).unwrap();
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

mod recovery {
    use crate::support::compressed_fixture::{bed, compressed_container, wig};
    use omicsweave_bbi::{BIGBED_MAGIC, BIGWIG_MAGIC, BbiReader};
    use omicsweave_core::{IntervalRead, Locus, Region, TrackEntry, TrackRead};
    #[test]
    fn bigwig_geometric_growth_declared_ceiling_and_truncated_stream() {
        let blocks = [wig(1, 10_000), wig(101, 1)];
        for limit in [64, u32::MAX] {
            let (bytes, _) = compressed_container(BIGWIG_MAGIC, &blocks, limit, false);
            let mut reader = BbiReader::open(bytes).unwrap();
            let short = Region::new("s", 100, 110);
            assert_eq!(reader.query_entries(short, true).unwrap(), 1);
            let result = reader.query_entries(Region::new("s", 0, 200), true);
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
            assert_eq!(reader.query_entries(short, true).unwrap(), 1);
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
                assert!(reader.query(region, false).is_err());
            } else {
                assert!(reader.query_entries(region, true).is_err());
            }
        }
    }
}
