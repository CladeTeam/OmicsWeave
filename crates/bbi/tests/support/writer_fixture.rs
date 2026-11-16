#![allow(dead_code)]
use omicsweave_bbi::{BbiSidecar, BedSchema, BigBedWriter, BigWigWriter, Chrom, WriteOptions};
use std::io::{self, Cursor, Seek, SeekFrom, Write};
pub(crate) fn chroms() -> Vec<Chrom> {
    vec![Chrom {
        name: "c".into(),
        len: 200_100,
    }]
}
pub(crate) fn options(explicit: bool, zooms: usize) -> WriteOptions {
    WriteOptions {
        compress: false,
        block_size: 4,
        zoom_levels: zooms,
        zoom_reductions: explicit.then(Vec::new),
        ..WriteOptions::default()
    }
}
pub(crate) fn wig(n: u32, explicit: bool, zooms: usize) -> (Vec<u8>, BbiSidecar) {
    let mut writer =
        BigWigWriter::create(Cursor::new(Vec::new()), chroms(), options(explicit, zooms)).unwrap();
    for i in 0..n {
        writer
            .add(0, i * 8, i * 8 + 4, (i % 97) as f32 - 48.0)
            .unwrap();
    }
    let (out, sidecar) = writer.finish().unwrap();
    (out.into_inner(), sidecar)
}
pub(crate) fn bed(n: u32, explicit: bool) -> (Vec<u8>, BbiSidecar) {
    let mut writer = BigBedWriter::create(
        Cursor::new(Vec::new()),
        chroms(),
        BedSchema::bed(4),
        options(explicit, 0),
    )
    .unwrap();
    for i in 0..n {
        writer.add(0, i * 3, i * 3 + 50, b"x").unwrap();
    }
    let (out, sidecar) = writer.finish().unwrap();
    (out.into_inner(), sidecar)
}

#[derive(Default)]
pub(crate) struct Discard(u64);
impl Write for Discard {
    fn write(&mut self, bytes: &[u8]) -> io::Result<usize> {
        self.0 += bytes.len() as u64;
        Ok(bytes.len())
    }
    fn flush(&mut self) -> io::Result<()> {
        Ok(())
    }
}
impl Seek for Discard {
    fn seek(&mut self, from: SeekFrom) -> io::Result<u64> {
        match from {
            SeekFrom::Start(at) => self.0 = at,
            _ => return Err(io::Error::other("only absolute seeks are used")),
        }
        Ok(self.0)
    }
}
