//! Coordinate validation and spliced locus traversal.

use omicsweave_core::{Error, Locus, Region, Span, Strand};

#[test]
fn loci_validate_and_describe_their_spans() {
    let spans = [
        Span::new(10, 20),
        Span::new(20, 25),
        Span::new(40, 40),
        Span::new(50, 60),
    ];
    let locus = Locus::new("c", &spans, Strand::Reverse).unwrap();
    assert_eq!(locus.len(), 25);
    assert!(locus.is_reverse());
    assert_eq!(locus.hull(), Region::new("c", 10, 60));
    let pieces: Vec<_> = locus.pieces().collect();
    assert_eq!(
        pieces,
        [
            (0, Region::new("c", 10, 20)),
            (10, Region::new("c", 20, 25)),
            (15, Region::new("c", 50, 60))
        ]
    );
    assert!(matches!(
        locus.check_within(59),
        Err(Error::RegionOutOfRange { .. })
    ));
    assert!(locus.check_within(60).is_ok());
    assert!(matches!(
        Locus::new(
            "c",
            &[Span::new(10, 20), Span::new(19, 30)],
            Strand::Forward
        ),
        Err(Error::Invalid { .. })
    ));
    let empty = Locus::new("c", &[], Strand::Forward).unwrap();
    assert!(empty.is_empty() && empty.hull().is_empty());
    assert!(matches!(
        Locus::from(Region::new("c", 5, 3)).check(),
        Err(Error::InvalidRegion { .. })
    ));
    assert!(matches!(
        Locus::from(Region::new("c", 5, 3)).check_within(10),
        Err(Error::InvalidRegion { .. })
    ));
}

#[test]
fn within_checks_converted_keys_and_extreme_spans() {
    for region in [Region::new("c", 7, 3), Region::new("c", u32::MAX, 2)] {
        assert!(matches!(
            Locus::from(region).check_within(0),
            Err(Error::InvalidRegion { .. })
        ));
    }
    let empty = Locus::new("c", &[], Strand::Reverse).unwrap();
    assert!(empty.check_within(0).is_ok());
    let at_end = Locus::from(Region::new("c", u32::MAX, u32::MAX));
    assert!(at_end.check_within(u32::MAX).is_ok());
    assert!(matches!(
        at_end.check_within(u32::MAX - 1),
        Err(Error::RegionOutOfRange { .. })
    ));
    let spans = [
        Span::new(0, 5),
        Span::new(10, 20),
        Span::new(u32::MAX, u32::MAX),
    ];
    let spliced = Locus::new("c", &spans, Strand::Forward).unwrap();
    assert!(spliced.check().is_ok());
    assert!(spliced.check_within(u32::MAX).is_ok());
    assert!(matches!(
        spliced.check_within(20),
        Err(Error::RegionOutOfRange { .. })
    ));
    assert!(matches!(
        spliced.check_within(2),
        Err(Error::RegionOutOfRange {
            start: u32::MAX,
            end: u32::MAX,
            len: 2
        })
    ));
    let whole_span = [Span::new(0, u32::MAX)];
    let whole = Locus::new("c", &whole_span, Strand::Forward).unwrap();
    assert_eq!(whole.len(), u32::MAX);
    assert!(whole.check_within(u32::MAX).is_ok());
}

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
