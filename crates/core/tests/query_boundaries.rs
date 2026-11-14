use omicsweave_core::{ContactQuery, Error, Region};

#[test]
fn region_overlap_matches_covered_bases_and_interior_points() {
    for query_start in 0..=8 {
        for query_end in 0..=8 {
            let query = Region::new("s", query_start, query_end);
            assert_eq!(query.len(), query_end.saturating_sub(query_start));
            assert_eq!(query.check().is_ok(), query_start <= query_end);
            for start in 0..=8 {
                for end in start..=8 {
                    let expected = if start == end {
                        query_start < start && start < query_end
                    } else {
                        (start..end).any(|base| (query_start..query_end).contains(&base))
                    };
                    assert_eq!(query.overlaps(start, end), expected);
                }
            }
        }
    }
}

#[test]
fn bins_match_bins_of_individual_bases() {
    for resolution in 1..=9 {
        for start in 0..=20 {
            for end in 0..=20 {
                let region = Region::new("s", start, end);
                let mut expected: Vec<_> = (start..end).map(|base| base / resolution).collect();
                expected.dedup();
                assert_eq!(
                    ContactQuery::bins(&region, resolution).collect::<Vec<_>>(),
                    expected
                );
            }
        }
    }
}

#[test]
// The oracle uses widened arithmetic independently of production div_ceil.
#[allow(clippy::manual_div_ceil)]
fn maximum_coordinates_do_not_overflow_rounding() {
    for start in [0, 1, u32::MAX - 2, u32::MAX] {
        let region = Region::new("s", start, u32::MAX);
        assert!(region.check_within(u32::MAX).is_ok());
        for resolution in [1, 2, 3, 1 << 31, u32::MAX] {
            let actual = ContactQuery::bins(&region, resolution);
            let lo = u64::from(start) / u64::from(resolution);
            let hi = if start == u32::MAX {
                lo
            } else {
                (u64::from(u32::MAX) + u64::from(resolution) - 1) / u64::from(resolution)
            };
            assert_eq!(u64::from(actual.start), lo);
            assert_eq!(u64::from(actual.end), hi);
            let query = ContactQuery::new(region, region, resolution);
            let bins = usize::try_from(hi - lo).unwrap();
            assert_eq!(query.shape(), (bins, bins));
        }
    }
    assert!(matches!(
        Region::new("s", u32::MAX, 0).check_within(0),
        Err(Error::InvalidRegion { .. })
    ));
    assert!(matches!(
        Region::new("s", u32::MAX, u32::MAX).check_within(u32::MAX - 1),
        Err(Error::RegionOutOfRange { .. })
    ));
}

#[test]
fn zero_resolution_is_rejected_even_for_empty_regions() {
    for region in [Region::new("s", 1, 1), Region::new("s", 1, 2)] {
        let query = ContactQuery::new(region, region, 0);
        assert!(query.check().is_err());
        assert!(std::panic::catch_unwind(|| query.shape()).is_err());
    }
}
