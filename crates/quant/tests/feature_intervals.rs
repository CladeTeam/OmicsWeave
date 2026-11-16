use omicsweave_gtf::{IntervalSet, Strand};
use omicsweave_quant::FeatureSpace;

fn intervals() -> IntervalSet {
    IntervalSet::from_columns(
        &["chr1", "chr2"],
        vec![0, 1, 0],
        vec![10, 20, 40],
        vec![15, 30, 40],
        vec![Strand::Forward, Strand::Reverse, Strand::Unstranded],
        None,
        Some(&["named", "", "point"]),
        None,
    )
    .unwrap()
}

#[test]
fn interval_columns_preserve_rows_names_and_zero_widths() {
    let space = FeatureSpace::of_intervals(intervals(), "Regions").unwrap();
    assert_eq!(space.column, [0, 1, 2]);
    assert_eq!(space.features.kind, "Regions");
    assert_eq!(
        space.features.id,
        ["chr1:10-15", "chr2:20-30", "chr1:40-40"]
    );
    assert_eq!(space.features.name, ["named", "chr2:20-30", "point"]);
    assert_eq!(space.features.length, [5, 10, 0]);
    assert_eq!(space.features.strand, ['+', '-', '.']);
    assert_eq!(
        space.match_references(&["chr2", "chr1"], None),
        [Some(1), Some(0)]
    );
}

#[test]
fn interval_columns_are_checked_before_building_features() {
    let mut malformed = Vec::new();
    let mut set = intervals();
    set.end.pop();
    malformed.push(("end column length", set));
    let mut set = intervals();
    set.seqid[0] = 2;
    malformed.push(("unknown sequence id", set));
    let mut set = intervals();
    set.start[0] = set.end[0] + 1;
    malformed.push(("inverted interval", set));
    let mut set = intervals();
    set.name[0] = set.pool().len() as u32;
    malformed.push(("name outside pool", set));
    for (case, set) in malformed {
        assert!(
            FeatureSpace::of_intervals(set, "Regions").is_err(),
            "{case}"
        );
    }
}

#[test]
fn feature_queries_use_the_owned_coordinate_columns() {
    let mut set = intervals();
    let mut hits = Vec::new();
    set.overlaps(0, 10, 15, &mut hits);
    assert_eq!(hits, [0]);
    set.start[0] = 100;
    set.end[0] = 105;
    let space = FeatureSpace::of_intervals(set, "Regions").unwrap();
    hits.clear();
    space.set.overlaps(0, 10, 15, &mut hits);
    assert!(hits.is_empty());
    space.set.overlaps(0, 100, 105, &mut hits);
    assert_eq!(hits, [0]);
}
