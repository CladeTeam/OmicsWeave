//! External interval columns, builder validation and index invalidation.
use omicsweave_gtf::{IntervalSet, IntervalSetBuilder, NONE, Strand};
fn sorted(mut ids: Vec<u32>) -> Vec<u32> {
    ids.sort_unstable();
    ids
}

#[test]
fn interval_columns_preserve_chromosome_order_and_reject_duplicates() {
    let names: Vec<String> = (0..1024).map(|i| format!("contig_{i}")).collect();
    let mut chroms: Vec<&str> = names.iter().map(String::as_str).collect();
    let make = |chroms: &[&str]| {
        IntervalSet::from_columns(chroms, vec![], vec![], vec![], vec![], None, None, None)
    };
    let set = make(&chroms).unwrap();
    assert_eq!(set.seqid_count(), chroms.len());
    for (i, name) in chroms.iter().enumerate() {
        assert_eq!(set.seqid_name(i as u32), *name);
    }
    for at in [0, 511, 1023] {
        chroms.push(chroms[at]);
        assert!(make(&chroms).unwrap_err().to_string().contains("twice"));
        chroms.pop();
    }
    assert!(make(&["", ""]).is_err());
    assert!(make(&[]).unwrap().is_empty());
}

#[test]
fn interval_sets_from_any_source() {
    let mut b = IntervalSetBuilder::new().chroms(["chr1", "chr2"]).unwrap();
    b.push("chr2", 100, 200, Strand::Unstranded, NONE, Some("peak1"))
        .unwrap();
    b.push("chr1", 5, 50, Strand::Forward, 7, None).unwrap();
    b.push("chr3", 0, 10, Strand::Reverse, NONE, Some("peak3"))
        .unwrap();
    assert!(b.push("chr1", 10, 5, Strand::Forward, NONE, None).is_err());
    let set = b.build();
    assert_eq!(set.seqid_count(), 3);
    assert_eq!(set.seqid_of("chr3"), Some(2));
    assert_eq!(set.name(0), Some("peak1"));
    let mut hits = Vec::new();
    set.overlaps(1, 150, 151, &mut hits);
    assert_eq!(hits, [0]);
    let mut c = set.cursor(0);
    hits.clear();
    c.query(40, 60, &mut hits);
    assert_eq!(hits, [1]);
    let sub = set.filter_owners(|o| o == 7);
    assert_eq!(sub.len(), 1);
    assert_eq!(sub.seqid_name(sub.seqid[0]), "chr1");
}

#[test]
fn edited_interval_columns_can_release_and_rebuild_the_index() {
    let mut set = IntervalSet::from_columns(
        &["chr1", "chr2"],
        vec![0, 0, 1],
        vec![1, 5, 10],
        vec![3, 15, 20],
        vec![Strand::Forward; 3],
        None,
        None,
        None,
    )
    .unwrap();
    assert_eq!(set.index().count_overlaps(0, 2, 8), 2);
    set.start[0] = 30;
    set.end[0] = 35;
    set.seqid[1] = 1;
    set.clear_index();
    for seqid in 0..=2 {
        let mut cursor = set.cursor(seqid);
        for start in 0..40 {
            let end = start + 5;
            let want: Vec<u32> = (0..set.len() as u32)
                .filter(|&r| {
                    let i = r as usize;
                    set.seqid[i] == seqid && set.start[i] < end && set.end[i] > start
                })
                .collect();
            let mut tree = Vec::new();
            set.overlaps(seqid, start, end, &mut tree);
            assert_eq!(sorted(tree), want);
            let mut sweep = Vec::new();
            cursor.query(start, end, &mut sweep);
            assert_eq!(sorted(sweep), want);
        }
    }
}

#[test]
fn interval_columns_reject_mismatched_lengths_invalid_seqids_and_reversed_spans() {
    for (label, seqid, start, end, strand, owner, names) in [
        (
            "seqid length",
            vec![],
            vec![1],
            vec![2],
            vec![Strand::Forward],
            None,
            None,
        ),
        (
            "end length",
            vec![0],
            vec![1],
            vec![],
            vec![Strand::Forward],
            None,
            None,
        ),
        (
            "strand length",
            vec![0],
            vec![1],
            vec![2],
            vec![],
            None,
            None,
        ),
        (
            "owner length",
            vec![0],
            vec![1],
            vec![2],
            vec![Strand::Forward],
            Some(vec![]),
            None,
        ),
        (
            "name length",
            vec![0],
            vec![1],
            vec![2],
            vec![Strand::Forward],
            None,
            Some(vec![]),
        ),
        (
            "unknown seqid",
            vec![1],
            vec![1],
            vec![2],
            vec![Strand::Forward],
            None,
            None,
        ),
        (
            "reversed span",
            vec![0],
            vec![2],
            vec![1],
            vec![Strand::Forward],
            None,
            None,
        ),
    ] {
        let result = IntervalSet::from_columns(
            &["chr1"],
            seqid,
            start,
            end,
            strand,
            owner,
            names.as_deref(),
            None,
        );
        assert!(result.is_err(), "accepted invalid columns: {label}");
    }
    let mut builder = IntervalSetBuilder::new().chroms(["chr1"]).unwrap();
    assert!(
        builder
            .push("invalid_chr", 2, 1, Strand::Forward, NONE, None)
            .is_err()
    );
    builder
        .push("chr1", 0, 1, Strand::Forward, 7, Some("one_base"))
        .unwrap();
    let set = builder.build();
    assert_eq!(set.len(), 1);
    assert_eq!(set.seqid_count(), 1);
    assert_eq!(set.name(0), Some("one_base"));
    assert_eq!(set.owner, [7]);
}

#[test]
fn cloned_interval_sets_rebuild_their_own_index_after_column_edits() {
    let original = IntervalSet::from_columns(
        &["chr1"],
        vec![0],
        vec![1],
        vec![5],
        vec![Strand::Forward],
        None,
        None,
        None,
    )
    .unwrap();
    assert_eq!(original.index().count_overlaps(0, 1, 2), 1);
    let mut changed = original.clone();
    changed.start[0] = 10;
    changed.end[0] = 15;
    changed.clear_index();
    assert_eq!(changed.index().count_overlaps(0, 1, 2), 0);
    assert_eq!(changed.index().count_overlaps(0, 10, 11), 1);
    assert_eq!(original.index().count_overlaps(0, 1, 2), 1);
    assert_eq!(original.index().count_overlaps(0, 10, 11), 0);
    let mut hits = Vec::new();
    changed.cursor(0).query(10, 11, &mut hits);
    assert_eq!(hits, [0]);
}
