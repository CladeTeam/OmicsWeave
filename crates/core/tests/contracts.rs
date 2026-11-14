use omicsweave_core::{IntervalSize, IntervalsOut, Region, StringColumn};

#[test]
fn empty_and_reversed_queries_never_overlap() {
    for region in [Region::new("s", 5, 5), Region::new("s", 6, 5)] {
        assert!(!region.overlaps(0, 10));
    }
    let region = Region::new("s", 5, 10);
    assert!(region.overlaps(7, 7));
    assert!(!region.overlaps(5, 5));
    assert!(!region.overlaps(10, 10));
}

#[test]
fn interval_buffers_require_exact_sizes() {
    let mut entries = [];
    let mut offsets = [99];
    let mut data = [];
    let out = IntervalsOut {
        entries: &mut entries,
        annotations: StringColumn {
            offsets: &mut offsets,
            data: &mut data,
        },
    };
    assert!(out.check(IntervalSize::default()).is_ok());
    assert!(
        out.check(IntervalSize {
            entries: usize::MAX,
            annotation_bytes: 0
        })
        .is_err()
    );
    assert!(
        out.check(IntervalSize {
            entries: 0,
            annotation_bytes: 1
        })
        .is_err()
    );
}
