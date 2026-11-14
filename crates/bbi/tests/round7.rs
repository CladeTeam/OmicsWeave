//! Writer output and failure boundaries, with variable block sizes and zooms.
use std::cell::Cell;
use std::hint::black_box;
use std::io::{self, Cursor, Seek, SeekFrom, Write};
use std::rc::Rc;
use std::time::{Duration, Instant};

use omicsweave_bbi::{BbiSidecar, BedSchema, BigBedWriter, BigWigWriter, Chrom, WriteOptions};
use omicsweave_core::exec::{Exec, Pool};
use omicsweave_core::sidecar::checksum;

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
    let pool = Pool::new(2);
    for compress in [false, true] {
        for zooms in [false, true] {
            for slots in [1, 3, 64] {
                let expected = bed(Exec::Inline, compress, zooms, slots);
                assert_eq!(BbiSidecar::build(&expected.0).unwrap(), expected.1);
                for limit in [1, 8] {
                    assert_eq!(
                        bed(Exec::Pool { pool: &pool, limit }, compress, zooms, slots),
                        expected
                    );
                }
            }
        }
    }
}

struct GatedWrite {
    out: Cursor<Vec<u8>>,
    fail: Rc<Cell<bool>>,
    writes: Rc<Cell<usize>>,
}

impl Write for GatedWrite {
    fn write(&mut self, bytes: &[u8]) -> io::Result<usize> {
        self.writes.set(self.writes.get() + 1);
        if self.fail.get() {
            return Err(io::Error::other("injected data-block write failure"));
        }
        self.out.write(bytes)
    }
    fn flush(&mut self) -> io::Result<()> {
        Ok(())
    }
}

impl Seek for GatedWrite {
    fn seek(&mut self, from: SeekFrom) -> io::Result<u64> {
        self.out.seek(from)
    }
}

#[test]
fn inline_data_write_failure_is_permanent_for_both_formats() {
    for compress in [false, true] {
        let fail = Rc::new(Cell::new(false));
        let writes = Rc::new(Cell::new(0));
        let out = || GatedWrite {
            out: Cursor::new(Vec::new()),
            fail: Rc::clone(&fail),
            writes: Rc::clone(&writes),
        };
        let mut wig = BigWigWriter::create(out(), chroms(), options(compress, false, 1)).unwrap();
        wig.add(0, 0, 1, 1.0).unwrap();
        fail.set(true);
        assert!(wig.add(0, 2, 3, 2.0).is_err());
        fail.set(false);
        let count = writes.get();
        assert!(wig.add(0, 4, 5, 3.0).is_err());
        assert!(wig.finish().is_err());
        assert_eq!(writes.get(), count);

        let mut bed = BigBedWriter::create(
            out(),
            chroms(),
            BedSchema::bed(4),
            options(compress, true, 1),
        )
        .unwrap();
        bed.add(0, 0, 1, b"a").unwrap();
        fail.set(true);
        assert!(bed.add(0, 2, 3, b"b").is_err());
        fail.set(false);
        let count = writes.get();
        assert!(bed.add(0, 4, 5, b"c").is_err());
        assert!(bed.finish().is_err());
        assert_eq!(writes.get(), count);
    }
}

#[test]
#[ignore = "root runs identical fixture source on baseline and candidate serially"]
fn paired_variable_block_benchmark() {
    let pool = Pool::new(2);
    for (mode, exec) in [
        ("inline", Exec::Inline),
        (
            "pool",
            Exec::Pool {
                pool: &pool,
                limit: 4,
            },
        ),
    ] {
        for compress in [false, true] {
            for zooms in [false, true] {
                for slots in [1, 64] {
                    let name =
                        format!("writer.bed.{mode}.compressed{compress}.zooms{zooms}.slots{slots}");
                    let expected = bed(exec, compress, zooms, slots);
                    let digest = checksum(&expected.0);
                    for _ in 0..4 {
                        black_box(bed(exec, compress, zooms, slots));
                    }
                    let mut samples = [0.0f64; 5];
                    for sample in &mut samples {
                        let start = Instant::now();
                        let mut n = 0;
                        while start.elapsed() < Duration::from_millis(100) {
                            black_box(bed(exec, compress, zooms, slots));
                            n += 1;
                        }
                        *sample = start.elapsed().as_nanos() as f64 / n as f64;
                    }
                    samples.sort_by(f64::total_cmp);
                    println!("{name}\tns={:.2}\tdigest={digest:016x}", samples[2]);
                }
            }
        }
    }
}
