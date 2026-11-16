use omicsweave_core::ContactUnit;
use omicsweave_hic::{BlockEntry, Depth, ZoomData};
use proptest::prelude::*;
use std::collections::BTreeSet;

fn zoom(bb: u32, columns: u32, intra: bool) -> ZoomData {
    ZoomData {
        chr1: 0,
        chr2: u32::from(!intra),
        unit: ContactUnit::Bp,
        bin_size: 1,
        sum_counts: 0.0,
        block_bin_count: bb,
        block_column_count: columns,
        blocks: (0..2048)
            .filter(|i| i % 3 != 0)
            .map(|number| BlockEntry {
                number,
                position: 0,
                size: 1,
            })
            .collect(),
    }
}

proptest! {
    #![proptest_config(ProptestConfig { failure_persistence: Some(Box::new(proptest::test_runner::FileFailurePersistence::Direct(concat!(env!("CARGO_MANIFEST_DIR"), "/tests/cases/index/selection.proptest-regressions")))), ..ProptestConfig::default() })]
    #[test]
    fn block_selection_matches_candidate_enumeration(
        bb in 1u32..16, columns in 1u32..32, x in 0u32..128, y in 0u32..128,
        w in 0u32..64, h in 0u32..64, intra in any::<bool>(), v9 in any::<bool>(), kind in 0u32..3,
    ) {
        let depth = match kind { 0 => Depth::Log(2), 1 => Depth::Log(3), _ => Depth::Linear(4) };
        let version = if v9 {9} else {8};
        let z = zoom(bb,columns,intra);
        let mut numbers = BTreeSet::new();
        if w>0 && h>0 {
            let (x1,x2,y1,y2)=(x as i64,(x+w-1) as i64,y as i64,(y+h-1) as i64);
            let bb=bb as i64;
            let rects = if v9 && intra {
                let (a,b)=(depth.of(x1-y2,bb),depth.of(x2-y1,bb));
                vec![(if x1<=y2 && y1<=x2 {0} else {a.min(b)}, a.max(b)+1,(x1+y1)/2/bb,(x2+y2)/2/bb+1)]
            } else {
                let (cl,ch,rl,rh)=(x1/bb,(x2+1)/bb,y1/bb,(y2+1)/bb);
                let mut rects=vec![(rl,rh,cl,ch)];
                if intra { rects.push((cl,ch,rl,rh)); }
                rects
            };
            for (rl,rh,cl,ch) in rects { for r in rl..=rh { for c in cl..=ch { numbers.insert((r*columns as i64+c) as u32); } } }
        }
        let want:Vec<_>=z.blocks.iter().copied().filter(|b| numbers.contains(&b.number)).collect();
        prop_assert_eq!(z.blocks_for(version,depth,x..x+w,y..y+h), want);
    }

}
#[test]
fn huge_queries_and_invalid_grids_are_bounded() {
    let z = zoom(1, 1, false);
    assert_eq!(
        z.blocks_for(8, Depth::Log(2), 0..u32::MAX, 0..u32::MAX),
        z.blocks
    );
    for (bb, cols) in [(0, 1), (1, 0)] {
        assert!(
            zoom(bb, cols, true)
                .blocks_for(9, Depth::Log(2), 0..1, 0..1)
                .is_empty()
        );
    }
}

fn stored_number(
    version: i32,
    depth: Depth,
    intra: bool,
    bb: u32,
    columns: u32,
    x: u32,
    y: u32,
) -> u32 {
    let bb = u64::from(bb);
    let x = u64::from(x);
    let y = u64::from(y);
    let row = if version == 9 && intra {
        let distance = x.abs_diff(y) as f64 / std::f64::consts::SQRT_2 / bb as f64;
        match depth {
            Depth::Log(base) => (1.0 + distance).log(f64::from(base)).floor() as u64,
            Depth::Linear(step) => (distance / f64::from(step)).floor() as u64,
        }
    } else {
        y / bb
    };
    let col = if version == 9 && intra {
        (x + y) / (2 * bb)
    } else {
        x / bb
    };
    u32::try_from(row * u64::from(columns) + col).unwrap()
}

#[test]
fn candidate_blocks_cover_actual_contacts_and_their_mirrors() {
    for version in [8, 9] {
        for intra in [false, true] {
            for depth in [Depth::Log(2), Depth::Log(7), Depth::Linear(3)] {
                for (base, bb, columns) in [
                    (0, 1, 2),
                    (0, 3, 20),
                    (0, 7, 3),
                    (u32::MAX - 32, 1_000_000_000, 3),
                ] {
                    let records: Vec<_> = (base..base + 24)
                        .flat_map(|x| {
                            (base..base + 24)
                                .filter(move |&y| !intra || x <= y)
                                .map(move |y| (x, y))
                        })
                        .collect();
                    let numbers: BTreeSet<_> = records
                        .iter()
                        .map(|&(x, y)| stored_number(version, depth, intra, bb, columns, x, y))
                        .collect();
                    let zoom = ZoomData {
                        chr1: 0,
                        chr2: u32::from(!intra),
                        unit: ContactUnit::Bp,
                        bin_size: 1,
                        sum_counts: 0.0,
                        block_bin_count: bb,
                        block_column_count: columns,
                        blocks: numbers
                            .into_iter()
                            .map(|number| BlockEntry {
                                number,
                                position: 0,
                                size: 0,
                            })
                            .collect(),
                    };
                    for (x, y, width, height) in [
                        (0, 0, 1, 1),
                        (3, 9, 5, 7),
                        (9, 3, 7, 5),
                        (7, 8, 9, 2),
                        (1, 1, 21, 21),
                    ] {
                        let xr = base + x..base + x + width;
                        let yr = base + y..base + y + height;
                        let blocks = zoom.blocks_for(version, depth, xr.clone(), yr.clone());
                        assert!(blocks.windows(2).all(|w| w[0].number < w[1].number));
                        let selected: BTreeSet<_> = blocks.iter().map(|b| b.number).collect();
                        for &(rx, ry) in &records {
                            if (xr.contains(&rx) && yr.contains(&ry))
                                || (intra && xr.contains(&ry) && yr.contains(&rx))
                            {
                                assert!(
                                    selected.contains(&stored_number(
                                        version, depth, intra, bb, columns, rx, ry
                                    )),
                                    "v{version}, {depth:?}, ({rx},{ry}), {xr:?} x {yr:?}"
                                );
                            }
                        }
                    }
                }
            }
        }
    }
}
#[test]
fn wide_rectangles_match_the_union_of_numeric_row_intervals() {
    let blocks: Vec<_> = (0..200_000)
        .filter(|number| number % 3 != 0)
        .map(|number| BlockEntry {
            number,
            position: 0,
            size: 1,
        })
        .collect();
    for intra in [false, true] {
        let zoom = ZoomData {
            chr1: 0,
            chr2: u32::from(!intra),
            unit: ContactUnit::Bp,
            bin_size: 1,
            sum_counts: 0.0,
            block_bin_count: 1,
            block_column_count: 16,
            blocks: blocks.clone(),
        };
        for (x, y) in [(3..40, 100..200), (3..40, 100..u32::MAX)] {
            let x_end = u64::from(x.end);
            let y_end = u64::from(y.end);
            let expected: Vec<_> = blocks
                .iter()
                .copied()
                .filter(|block| {
                    let number = u64::from(block.number);
                    (u64::from(y.start) * 16 + u64::from(x.start)..=y_end * 16 + x_end)
                        .contains(&number)
                        || (intra
                            && (u64::from(x.start) * 16 + u64::from(y.start)..=x_end * 16 + y_end)
                                .contains(&number))
                })
                .collect();
            assert_eq!(zoom.blocks_for(8, Depth::Log(2), x, y), expected);
        }
    }
    let zoom = ZoomData {
        chr1: 0,
        chr2: 1,
        unit: ContactUnit::Bp,
        bin_size: 1,
        sum_counts: 0.0,
        block_bin_count: 1,
        block_column_count: u32::MAX - 3,
        blocks,
    };
    // Even the first numeric interval starts above i64::MAX; there is no
    // block number in the format's u32 range and no reason to visit rows.
    assert!(
        zoom.blocks_for(8, Depth::Log(2), 0..u32::MAX, u32::MAX - 2..u32::MAX)
            .is_empty()
    );
}
