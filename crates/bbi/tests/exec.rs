//! Writers compressing on a worker pool: for any thread count and bound on
//! blocks in flight, the file and sidecar are byte for byte those written
//! with inline compression.

use std::io::Cursor;
use std::sync::LazyLock;

use omicsweave_bbi::{BbiSidecar, BedSchema, BigBedWriter, BigWigWriter, Chrom, WriteOptions};
use omicsweave_core::exec::{Exec, Pool};
use proptest::prelude::*;

static POOL4: LazyLock<Pool> = LazyLock::new(|| Pool::new(4));
static POOL8: LazyLock<Pool> = LazyLock::new(|| Pool::new(8));

/// Inline, then pools with a single block in flight, a few, and many.
fn executors() -> [Exec<'static>; 5] {
    [
        Exec::Inline,
        Exec::Pool {
            pool: &POOL4,
            limit: 1,
        },
        Exec::Pool {
            pool: &POOL4,
            limit: 3,
        },
        Exec::Pool {
            pool: &POOL8,
            limit: 2,
        },
        Exec::Pool {
            pool: &POOL8,
            limit: 10_000,
        },
    ]
}

fn chroms(lens: &[u32]) -> Vec<Chrom> {
    lens.iter()
        .enumerate()
        .map(|(i, &len)| Chrom {
            name: format!("chr{i}"),
            len,
        })
        .collect()
}

type Written = (Vec<u8>, BbiSidecar);

/// bedGraph items `(chrom, start, end, value)`, sorted and disjoint.
fn bigwig(lens: &[u32], items: &[(u32, u32, u32, f32)], o: &WriteOptions, e: Exec<'_>) -> Written {
    let out = Cursor::new(Vec::new());
    let mut w = BigWigWriter::create_with(out, chroms(lens), o.clone(), e).unwrap();
    for &(c, s, end, v) in items {
        w.add(c, s, end, v).unwrap();
    }
    let (out, sidecar) = w.finish().unwrap();
    (out.into_inner(), sidecar)
}

/// Features `(chrom, start, end, rest)`, sorted by start, overlapping.
fn bigbed(
    lens: &[u32],
    items: &[(u32, u32, u32, Vec<u8>)],
    o: &WriteOptions,
    e: Exec<'_>,
) -> Written {
    let out = Cursor::new(Vec::new());
    let schema = BedSchema::bed(4);
    let mut w = BigBedWriter::create_with(out, chroms(lens), schema, o.clone(), e).unwrap();
    for (c, s, end, rest) in items {
        w.add(*c, *s, *end, rest).unwrap();
    }
    let (out, sidecar) = w.finish().unwrap();
    (out.into_inner(), sidecar)
}

/// Walks the chromosomes, placing items at random gaps and spans; bigBed
/// items may overlap the ones before them.
fn layout_items(
    lens: &[u32],
    steps: &[(u8, u16, u16)],
    overlap: bool,
) -> Vec<(u32, u32, u32, usize)> {
    let mut out = Vec::new();
    let (mut chrom, mut pos) = (0usize, 0u32);
    for (i, &(skip, gap, span)) in steps.iter().enumerate() {
        if skip == 0 {
            chrom += 1;
            pos = 0;
        }
        let Some(&len) = lens.get(chrom) else { break };
        let start = pos.saturating_add(gap as u32);
        let end = start.saturating_add(span as u32 + 1).min(len);
        if start >= end {
            continue;
        }
        out.push((chrom as u32, start, end, i));
        pos = if overlap { start } else { end };
    }
    out
}

fn options() -> impl Strategy<Value = WriteOptions> {
    (
        prop::option::of(1u32..64),
        2u32..9,
        prop::bool::weighted(0.8),
        1i32..=12,
        0usize..10,
        prop::option::weighted(0.2, prop::collection::vec(1u32..5_000, 1..4)),
    )
        .prop_map(
            |(items_per_slot, block_size, compress, level, zoom_levels, zoom_reductions)| {
                WriteOptions {
                    items_per_slot,
                    block_size,
                    compress,
                    level,
                    zoom_levels,
                    zoom_reductions,
                }
            },
        )
}

proptest! {
    #![proptest_config(ProptestConfig::with_cases(64))]

    #[test]
    fn bigwig_is_the_same_on_every_executor(
        lens in prop::collection::vec(100u32..200_000, 1..4),
        steps in prop::collection::vec((0u8..40, 0u16..300, 0u16..200), 0..4_000),
        values in prop::collection::vec(-1e6f32..1e6, 1..64),
        options in options(),
    ) {
        let items: Vec<_> = layout_items(&lens, &steps, false)
            .into_iter()
            .map(|(c, s, e, i)| (c, s, e, values[i % values.len()]))
            .collect();
        let want = bigwig(&lens, &items, &options, Exec::Inline);
        for e in executors() {
            let got = bigwig(&lens, &items, &options, e);
            prop_assert!(got.0 == want.0, "file differs");
            prop_assert_eq!(&got.1, &want.1);
        }
    }

    #[test]
    fn bigbed_is_the_same_on_every_executor(
        lens in prop::collection::vec(100u32..200_000, 1..4),
        steps in prop::collection::vec((0u8..40, 0u16..300, 0u16..2_000), 0..3_000),
        names in prop::collection::vec("[a-z]{0,12}", 1..16),
        options in options(),
    ) {
        let items: Vec<_> = layout_items(&lens, &steps, true)
            .into_iter()
            .map(|(c, s, e, i)| (c, s, e, names[i % names.len()].as_bytes().to_vec()))
            .collect();
        let want = bigbed(&lens, &items, &options, Exec::Inline);
        for e in executors() {
            let got = bigbed(&lens, &items, &options, e);
            prop_assert!(got.0 == want.0, "file differs");
            prop_assert_eq!(&got.1, &want.1);
        }
    }
}

/// Default options over many blocks and zoom levels: every pool keeps far
/// more than one block in flight.
#[test]
fn large_bigwig_is_the_same_on_every_executor() {
    let lens = [5_000_000, 3_000_000];
    let mut items = Vec::new();
    for (c, &len) in lens.iter().enumerate() {
        let mut pos = 0u32;
        let mut i = 0u32;
        while pos + 50 < len {
            let span = 1 + (i.wrapping_mul(2_654_435_761) >> 27);
            items.push((c as u32, pos, pos + span, (i % 97) as f32 * 0.25));
            pos += span + (i % 3);
            i += 1;
        }
    }
    let options = WriteOptions::default();
    let want = bigwig(&lens, &items, &options, Exec::Inline);
    assert!(want.1.levels.len() > 2);
    for e in executors() {
        let got = bigwig(&lens, &items, &options, e);
        assert!(got.0 == want.0, "file differs");
        assert_eq!(got.1, want.1);
    }
}
