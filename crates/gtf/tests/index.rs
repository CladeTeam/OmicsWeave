//! The interval index and sweep cursor against naive scans.

use omicsweave_gtf::{IntervalIndex, IntervalSet, IntervalSetBuilder, NONE, Strand};
use proptest::prelude::*;

type Item = (u32, u32, u32);

fn items() -> impl Strategy<Value = Vec<Item>> {
    prop::collection::vec(
        (
            0u32..3,
            0u32..2000,
            prop_oneof![0u32..5, 0u32..60, 0u32..1500],
        ),
        0..400,
    )
    .prop_map(|v| v.into_iter().map(|(s, st, l)| (s, st, st + l)).collect())
}

fn build(items: &[Item]) -> IntervalIndex {
    let seq: Vec<u32> = items.iter().map(|x| x.0).collect();
    let st: Vec<u32> = items.iter().map(|x| x.1).collect();
    let en: Vec<u32> = items.iter().map(|x| x.2).collect();
    IntervalIndex::build(3, &seq, &st, &en)
}

fn naive(items: &[Item], s: u32, pred: impl Fn(u32, u32) -> bool) -> Vec<u32> {
    (0..items.len() as u32)
        .filter(|&i| {
            let (seq, st, en) = items[i as usize];
            seq == s && pred(st, en)
        })
        .collect()
}

fn overlap(qs: u32, qe: u32) -> impl Fn(u32, u32) -> bool {
    move |st, en| qs < qe && st < qe && qs < en
}

fn sorted(mut v: Vec<u32>) -> Vec<u32> {
    v.sort_unstable();
    v
}

proptest! {
    #![proptest_config(ProptestConfig::with_cases(200))]

    #[test]
    fn queries_match_scans(items in items(), queries in prop::collection::vec((0u32..4, 0u32..2200, 0u32..300), 1..60)) {
        let ix = build(&items);
        for (s, qs, l) in queries {
            let qe = qs + l;
            let mut got = Vec::new();
            ix.overlaps(s, qs, qe, &mut got);
            // Results come in start order.
            let starts: Vec<u32> = got.iter().map(|&i| items[i as usize].1).collect();
            prop_assert!(starts.windows(2).all(|w| w[0] <= w[1]));
            prop_assert_eq!(sorted(got), naive(&items, s, overlap(qs, qe)));
            prop_assert_eq!(ix.count_overlaps(s, qs, qe), naive(&items, s, overlap(qs, qe)).len());

            let mut got = Vec::new();
            ix.contained_in(s, qs, qe, &mut got);
            prop_assert_eq!(sorted(got), naive(&items, s, |st, en| overlap(qs, qe)(st, en) && st >= qs && en <= qe));

            let qe1 = qe.max(qs + 1);
            let mut got = Vec::new();
            ix.containing(s, qs, qe, &mut got);
            prop_assert_eq!(sorted(got), naive(&items, s, |st, en| st <= qs && en >= qe1 && st < en));

            // Nearest: overlap first, else the smallest gap on either side.
            let mut got = Vec::new();
            let d = ix.nearest(s, qs, qe, &mut got);
            let ov = naive(&items, s, overlap(qs, qe1));
            let any: Vec<u32> = naive(&items, s, |_, _| true);
            if any.is_empty() {
                prop_assert_eq!(d, None);
            } else if !ov.is_empty() {
                prop_assert_eq!(d, Some(0));
                prop_assert_eq!(sorted(got), ov);
            } else {
                let gap = |i: u32| {
                    let (_, st, en) = items[i as usize];
                    if st >= qe1 { st - qe1 } else { qs - en }
                };
                let best = any.iter().map(|&i| gap(i)).min().unwrap();
                prop_assert_eq!(d, Some(best));
                let want: Vec<u32> = any.iter().copied().filter(|&i| gap(i) == best).collect();
                prop_assert_eq!(sorted(got), want);
            }
        }
    }

    #[test]
    fn sweep_cursor_matches_scans(items in items(), mut steps in prop::collection::vec((0u32..2200, 0u32..400, prop::collection::vec((0u32..800, 0u32..200), 1..4)), 1..80)) {
        let ix = build(&items);
        steps.sort_by_key(|s| s.0);
        for s in 0..3 {
            // Monotone single queries.
            let mut c = ix.cursor(s);
            let mut got = Vec::new();
            for &(qs, l, _) in &steps {
                got.clear();
                c.query(qs, qs + l, &mut got);
                prop_assert_eq!(sorted(got.clone()), naive(&items, s, overlap(qs, qs + l)));
            }
            // A spliced read: one advance, then blocks at or after it.
            c.seek(s);
            for (qs, _, blocks) in &steps {
                c.advance(*qs);
                let mut at = *qs;
                for &(gap, len) in blocks {
                    let bs = at + gap;
                    got.clear();
                    c.overlaps(bs, bs + len, &mut got);
                    prop_assert_eq!(sorted(got.clone()), naive(&items, s, overlap(bs, bs + len)));
                    at = bs + len;
                }
            }
        }
    }

    #[test]
    fn cursor_reseeds_on_backward_queries(items in items(), queries in prop::collection::vec((0u32..2200, 0u32..300), 1..60)) {
        let ix = build(&items);
        let mut c = ix.cursor(1);
        let mut got = Vec::new();
        for (qs, l) in queries {
            got.clear();
            c.query(qs, qs + l, &mut got);
            prop_assert_eq!(sorted(got.clone()), naive(&items, 1, overlap(qs, qs + l)));
        }
    }

    #[test]
    fn batch_matches_single_queries(items in items(), queries in prop::collection::vec((0u32..3, 0u32..2200, 0u32..300), 0..40)) {
        let ix = build(&items);
        let qs: Vec<(u32, u32, u32)> = queries.iter().map(|&(s, a, l)| (s, a, a + l)).collect();
        let (mut off, mut ids) = (Vec::new(), Vec::new());
        ix.overlaps_batch(qs.iter().copied(), &mut off, &mut ids);
        prop_assert_eq!(off.len(), qs.len() + 1);
        for (q, &(s, a, b)) in qs.iter().enumerate() {
            let mut one = Vec::new();
            ix.overlaps(s, a, b, &mut one);
            prop_assert_eq!(&ids[off[q] as usize..off[q + 1] as usize], one.as_slice());
        }
    }
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
fn nearest_nested_intervals_ties_and_coordinate_edges() {
    let mut items = vec![(0, 0, 1_000_000)];
    items.extend((1..20_000).map(|i| (0, i * 2, i * 2 + 1)));
    items.extend([(0, 900_000, 1_000_000), (0, 1_000_000, 1_000_000)]);
    let right = items.len() as u32;
    items.push((0, 1_000_020, 1_000_030));
    items.extend([(1, 0, 0), (1, u32::MAX - 10, u32::MAX)]);
    let ix = build(&items);
    let mut got = vec![u32::MAX];
    assert_eq!(ix.nearest(0, 1_000_010, 1_000_010, &mut got), Some(9));
    assert_eq!(got, [u32::MAX, right]);
    got.clear();
    assert_eq!(ix.nearest(0, 1_000_009, 1_000_011, &mut got), Some(9));
    assert_eq!(got, [0, right - 2, right - 1, right]);
    for (seqid, start, end) in [
        (0, 1_000_001, 1_000_001),
        (0, 1_000_002, 1_000_001),
        (1, 0, 0),
        (1, u32::MAX - 1, u32::MAX),
        (1, u32::MAX, u32::MAX),
        (2, 0, 1),
        (u32::MAX, 0, 1),
    ] {
        let qe = end.max(start.saturating_add(1));
        let mut want = naive(&items, seqid, overlap(start, qe));
        let distance = if want.is_empty() {
            let candidates = naive(&items, seqid, |_, _| true);
            let gap = |i: u32| {
                let (_, s, e) = items[i as usize];
                if s >= qe { s - qe } else { start - e }
            };
            let best = candidates.iter().map(|&i| gap(i)).min();
            want = candidates
                .into_iter()
                .filter(|&i| Some(gap(i)) == best)
                .collect();
            best
        } else {
            Some(0)
        };
        got.clear();
        assert_eq!(ix.nearest(seqid, start, end, &mut got), distance);
        assert_eq!(sorted(got.clone()), want);
    }
}

#[test]
fn index_rejects_invalid_columns_before_querying() {
    assert!(std::panic::catch_unwind(|| IntervalIndex::build(1, &[0], &[], &[1])).is_err());
    assert!(std::panic::catch_unwind(|| IntervalIndex::build(1, &[0], &[0], &[])).is_err());
    assert!(std::panic::catch_unwind(|| IntervalIndex::build_with(1, [(1, 0, 1, 0)])).is_err());
    assert!(std::panic::catch_unwind(|| IntervalIndex::build_with(1, [(0, 2, 1, 0)])).is_err());
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
