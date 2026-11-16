use std::cell::Cell;
use std::io::{self, Cursor, Seek, SeekFrom, Write};
use std::rc::Rc;

use omicsweave_bbi::{BedSchema, BigBedWriter, BigWigWriter, Chrom, WriteOptions};

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
