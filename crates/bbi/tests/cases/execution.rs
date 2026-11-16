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
        for e in executors().into_iter().skip(1) {
            let got = bigwig(&lens, &items, &options, e);
            prop_assert!(got.0 == want.0, "file differs: options={:?}, lens={:?}, items={:?}", options, lens, items);
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
        for e in executors().into_iter().skip(1) {
            let got = bigbed(&lens, &items, &options, e);
            prop_assert!(got.0 == want.0, "file differs: options={:?}, lens={:?}, items={:?}", options, lens, items);
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
    for e in executors().into_iter().skip(1) {
        let got = bigwig(&lens, &items, &options, e);
        assert!(got.0 == want.0, "file differs: options={options:?}");
        assert_eq!(got.1, want.1);
    }
}

mod variable_blocks {
    use super::*;
    fn chroms() -> Vec<Chrom> {
        vec![Chrom {
            name: "c".into(),
            len: 100_000,
        }]
    }

    fn options(compress: bool, zooms: bool, slots: u32) -> WriteOptions {
        WriteOptions {
            compress,
            items_per_slot: Some(slots),
            zoom_reductions: Some(if zooms { vec![8, 128] } else { vec![] }),
            ..WriteOptions::default()
        }
    }

    fn bed(exec: Exec<'_>, compress: bool, zooms: bool, slots: u32) -> (Vec<u8>, BbiSidecar) {
        let mut writer = BigBedWriter::create_with(
            Cursor::new(Vec::new()),
            chroms(),
            BedSchema::bed(4),
            options(compress, zooms, slots),
            exec,
        )
        .unwrap();
        let notes: Vec<_> = [1, 8192, 2, 257]
            .into_iter()
            .map(|len| vec![b'x'; len])
            .collect();
        for i in 0..257 {
            writer
                .add(0, i * 4, i * 4 + 17, &notes[i as usize % notes.len()])
                .unwrap();
        }
        let (out, sidecar) = writer.finish().unwrap();
        (out.into_inner(), sidecar)
    }

    #[test]
    fn variable_base_and_zoom_blocks_match_ordered_pool_bytes() {
        for compress in [false, true] {
            for zooms in [false, true] {
                for slots in [1, 3, 64] {
                    let expected = bed(Exec::Inline, compress, zooms, slots);
                    assert_eq!(BbiSidecar::build(&expected.0).unwrap(), expected.1);
                    for exec in executors().into_iter().skip(1) {
                        assert_eq!(bed(exec, compress, zooms, slots), expected);
                    }
                }
            }
        }
    }
}

use std::process::{Command, Stdio};
use std::time::{Duration, Instant};
fn closed_pool_chroms() -> Vec<Chrom> {
    vec![
        Chrom {
            name: "a".into(),
            len: 100_000,
        },
        Chrom {
            name: "b".into(),
            len: 100_000,
        },
    ]
}
fn closed_pool_options(compress: bool) -> WriteOptions {
    WriteOptions {
        compress,
        items_per_slot: Some(3),
        block_size: 4,
        zoom_reductions: Some(vec![16, 128]),
        ..WriteOptions::default()
    }
}
#[test]
fn writer_rejects_submissions_after_pool_shutdown_without_hanging() {
    let executable = std::env::current_exe().unwrap();
    let marker = tempfile::NamedTempFile::new().unwrap();
    let mut child = Command::new(executable)
        .env("OMICSWEAVE_BBI_CLOSED_POOL_MARKER", marker.path())
        .args([
            "--ignored",
            "--exact",
            "execution::closed_pool_child",
            "--nocapture",
        ])
        .stdin(Stdio::null())
        .spawn()
        .unwrap();
    let deadline = Instant::now() + Duration::from_secs(5);
    loop {
        if let Some(status) = child.try_wait().unwrap() {
            assert!(status.success(), "closed-pool child failed: {status}");
            assert_eq!(
                std::fs::read(marker.path()).unwrap(),
                b"closed-pool assertions completed"
            );
            return;
        }
        if Instant::now() >= deadline {
            let _ = child.kill();
            let _ = child.wait();
            panic!("writer or queue drop hung after its pool was dropped");
        }
        std::thread::sleep(Duration::from_millis(10));
    }
}

#[test]
#[ignore = "invoked by bounded parent process"]
fn closed_pool_child() {
    let pool = Pool::new(1);
    let mut writer_options = closed_pool_options(false);
    writer_options.items_per_slot = Some(1);
    writer_options.zoom_reductions = Some(vec![]);
    let mut w = BigWigWriter::create_with(
        Cursor::new(Vec::new()),
        closed_pool_chroms(),
        writer_options,
        Exec::Pool {
            pool: &pool,
            limit: 1,
        },
    )
    .unwrap();
    w.add(0, 0, 1, 1.0).unwrap();
    drop(pool);
    assert!(w.add(0, 1, 2, 2.0).is_err());
    assert!(w.add(0, 2, 3, 3.0).is_err());
    assert!(
        w.finish().is_err(),
        "failed submission must permanently fail the writer"
    );

    let pool = Pool::new(1);
    let mut writer_options = closed_pool_options(false);
    writer_options.items_per_slot = Some(1);
    writer_options.zoom_reductions = Some(vec![]);
    let mut w = BigBedWriter::create_with(
        Cursor::new(Vec::new()),
        closed_pool_chroms(),
        BedSchema::bed(4),
        writer_options,
        Exec::Pool {
            pool: &pool,
            limit: 1,
        },
    )
    .unwrap();
    w.add(0, 0, 1, b"a").unwrap();
    drop(pool);
    assert!(w.add(0, 1, 2, b"b").is_err());
    assert!(w.add(0, 2, 3, b"c").is_err());
    assert!(w.finish().is_err());
    if let Some(path) = std::env::var_os("OMICSWEAVE_BBI_CLOSED_POOL_MARKER") {
        std::fs::write(path, b"closed-pool assertions completed").unwrap();
    }
}
