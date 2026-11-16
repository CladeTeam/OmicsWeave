//! The interval index and sweep cursor against naive scans.

use crate::support::{LEN, random};
use omicsweave_gtf::{IntervalIndex, Level};
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
fn cursor_seek_backward_floor_and_equal_boundaries_preserve_appended_results() {
    let items = [
        (0, 10, 20),
        (0, 10, 30),
        (0, 30, 40),
        (1, 10, 20),
        (0, 0, 100),
    ];
    let ix = build(&items);
    let mut cursor = ix.cursor(0);
    for (seqid, start, end) in [
        (0, 50, 60),
        (0, 10, 11),
        (0, 20, 30),
        (0, 30, 31),
        (0, 10, 10),
        (1, 10, 20),
        (3, 0, 10),
        (0, 0, 10),
    ] {
        if cursor.seqid() != seqid {
            cursor.seek(seqid);
            assert_eq!(cursor.floor(), 0);
            assert_eq!(cursor.active_len(), 0);
        }
        let mut got = vec![u32::MAX];
        cursor.query(start, end, &mut got);
        assert_eq!(cursor.floor(), start);
        assert_eq!(got.remove(0), u32::MAX);
        assert_eq!(
            sorted(got),
            naive(&items, seqid, overlap(start, end)),
            "seqid={seqid}, query={start}-{end}"
        );
    }
    cursor.advance(50);
    let mut got = Vec::new();
    cursor.overlaps(10, 11, &mut got);
    assert_eq!(cursor.floor(), 10);
    assert_eq!(sorted(got), [0, 1, 4]);
}

#[test]
fn exon_index_and_cursor_cover_feature_rows() {
    let a = random(11);
    let ix = a.interval_index(Level::Exon);
    assert_eq!(ix.len(), a.exons().len());
    let f = a.features();
    for s in 0..a.seqid_count() as u32 {
        let mut c = a.cursor(Level::Exon, s);
        let mut got = Vec::new();
        for p in (0..LEN).step_by(37) {
            got.clear();
            c.query(p, p + 50, &mut got);
            let mut want: Vec<u32> = a
                .exons()
                .iter()
                .copied()
                .filter(|&e| {
                    let i = e as usize;
                    a.feature_seqid(e) == s && f.start[i] < p + 50 && p < f.end[i]
                })
                .collect();
            want.sort_unstable();
            got.sort_unstable();
            assert_eq!(got, want);
        }
    }
}
