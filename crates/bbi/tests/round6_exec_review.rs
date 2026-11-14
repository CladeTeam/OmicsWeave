//! Ordered writer jobs are compared against the inline byte stream. Pool
//! shutdown is exercised in a subprocess so an executor bug cannot hang tests.
use std::io::Cursor;
use std::process::{Command, Stdio};
use std::time::{Duration, Instant};

use omicsweave_bbi::{BbiSidecar, BedSchema, BigBedWriter, BigWigWriter, Chrom, WriteOptions};
use omicsweave_core::exec::{Exec, Pool};

fn chroms() -> Vec<Chrom> {
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
fn options(compress: bool) -> WriteOptions {
    WriteOptions {
        compress,
        items_per_slot: Some(3),
        block_size: 4,
        zoom_reductions: Some(vec![16, 128]),
        ..WriteOptions::default()
    }
}
fn wig(exec: Exec<'_>, compress: bool) -> (Vec<u8>, BbiSidecar) {
    let mut w =
        BigWigWriter::create_with(Cursor::new(Vec::new()), chroms(), options(compress), exec)
            .unwrap();
    for chrom in 0..2 {
        for i in 0..128 {
            w.add(chrom, i * 8, i * 8 + 4, (i % 31) as f32 - 15.0)
                .unwrap();
        }
    }
    let (out, sidecar) = w.finish().unwrap();
    (out.into_inner(), sidecar)
}
fn bed(exec: Exec<'_>, compress: bool) -> (Vec<u8>, BbiSidecar) {
    let mut w = BigBedWriter::create_with(
        Cursor::new(Vec::new()),
        chroms(),
        BedSchema::bed(4),
        options(compress),
        exec,
    )
    .unwrap();
    for chrom in 0..2 {
        for i in 0..128 {
            // Alternating block sizes give workers different amounts of work.
            let note = vec![b'a' + (i % 20) as u8; if i % 7 == 0 { 17_000 } else { 3 }];
            w.add(chrom, i * 3, i * 3 + 50, &note).unwrap();
        }
    }
    let (out, sidecar) = w.finish().unwrap();
    (out.into_inner(), sidecar)
}

#[test]
fn inline_and_bounded_pool_writers_have_identical_bytes_and_sidecars() {
    let pool = Pool::new(3);
    for compress in [false, true] {
        let expected_wig = wig(Exec::Inline, compress);
        let expected_bed = bed(Exec::Inline, compress);
        for limit in [1, 2, 8] {
            let exec = Exec::Pool { pool: &pool, limit };
            let actual_wig = wig(exec, compress);
            let actual_bed = bed(exec, compress);
            assert_eq!(actual_wig, expected_wig);
            assert_eq!(actual_bed, expected_bed);
            assert_eq!(BbiSidecar::build(&actual_wig.0).unwrap(), actual_wig.1);
            assert_eq!(BbiSidecar::build(&actual_bed.0).unwrap(), actual_bed.1);
        }
    }
}

#[test]
fn writer_rejects_submissions_after_pool_shutdown_without_hanging() {
    let executable = std::env::current_exe().unwrap();
    let mut child = Command::new(executable)
        .args(["--ignored", "--exact", "closed_pool_child", "--nocapture"])
        .stdin(Stdio::null())
        .spawn()
        .unwrap();
    let deadline = Instant::now() + Duration::from_secs(5);
    loop {
        if let Some(status) = child.try_wait().unwrap() {
            assert!(status.success(), "closed-pool child failed: {status}");
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
    let mut writer_options = options(false);
    writer_options.items_per_slot = Some(1);
    writer_options.zoom_reductions = Some(vec![]);
    let mut w = BigWigWriter::create_with(
        Cursor::new(Vec::new()),
        chroms(),
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
    let mut writer_options = options(false);
    writer_options.items_per_slot = Some(1);
    writer_options.zoom_reductions = Some(vec![]);
    let mut w = BigBedWriter::create_with(
        Cursor::new(Vec::new()),
        chroms(),
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
}
