//! Copy this identical fixture into the baseline and candidate source trees.
use omicsweave_core::{ReadAt, Region, Result, SequenceRead};
use omicsweave_fasta::{TwoBitReader, TwoBitWriter, encode_twobit};
use std::{
    cell::Cell,
    hint::black_box,
    io::Cursor,
    time::{Duration, Instant},
};

struct Counted<'a> {
    bytes: &'a [u8],
    calls: Cell<usize>,
    read: Cell<usize>,
}
impl ReadAt for Counted<'_> {
    fn len(&self) -> u64 {
        self.bytes.len() as u64
    }
    fn read_exact_at(&self, out: &mut [u8], at: u64) -> Result<()> {
        self.calls.set(self.calls.get() + 1);
        self.read.set(self.read.get() + out.len());
        self.bytes.read_exact_at(out, at)
    }
}
fn measure(name: &str, mut f: impl FnMut()) {
    for _ in 0..20 {
        f();
    }
    let mut samples = [0f64; 5];
    for sample in &mut samples {
        let start = Instant::now();
        let mut count = 0;
        while start.elapsed() < Duration::from_millis(100) {
            for _ in 0..16 {
                f();
            }
            count += 16;
        }
        *sample = start.elapsed().as_nanos() as f64 / count as f64;
    }
    samples.sort_by(f64::total_cmp);
    println!("{name}\t{:.2}", samples[2]);
}
fn main() {
    for (kind, bases) in [
        ("ordinary", b"TCAG".repeat(4096)),
        ("blocks", b"nacGTN".repeat(2731)),
    ] {
        for count in [8, 64] {
            let names: Vec<_> = (0..count).map(|i| format!("s{i}")).collect();
            let mut writer =
                TwoBitWriter::create(Cursor::new(Vec::new()), names.clone(), 1).unwrap();
            let (record, packed) = encode_twobit(&bases).unwrap();
            for name in &names {
                writer.add(name, &record, &mut packed.as_slice()).unwrap();
            }
            let bytes = writer.finish().unwrap().into_inner();
            let source = Counted {
                bytes: &bytes,
                calls: Cell::new(0),
                read: Cell::new(0),
            };
            let mut reader = TwoBitReader::open(&source, true).unwrap();
            println!(
                "io\tmetadata.{kind}.{count}\t{}\t{}",
                source.calls.get(),
                source.read.get()
            );
            let mut checksum = 0u64;
            for name in &names {
                let mut out = vec![0; bases.len()];
                reader
                    .sequence(Region::new(name, 0, bases.len() as u32), &mut out)
                    .unwrap();
                let expected: Vec<_> = bases
                    .iter()
                    .map(|&b| if b == b'n' { b'N' } else { b })
                    .collect();
                assert_eq!(out, expected);
                checksum = out
                    .iter()
                    .fold(checksum, |h, &b| h.wrapping_mul(131).wrapping_add(b as u64));
            }
            println!("checksum\tmetadata.{kind}.{count}\t{checksum:016x}");
            measure(&format!("metadata.{kind}.{count}.open"), || {
                black_box(TwoBitReader::open(black_box(&source), true).unwrap());
            });
        }
    }
}
