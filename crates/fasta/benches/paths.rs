//! Warm-buffer sequence paths and positional I/O counts, excluding fixture setup.
#[path = "../tests/common/mod.rs"]
mod common;

use omicsweave_core::{ReadAt, Region, Result, SequenceRead};
use omicsweave_fasta::{FaiBuilder, FastaReader, TwoBitReader, compress_fasta};
use std::{
    cell::Cell,
    hint::black_box,
    rc::Rc,
    time::{Duration, Instant},
};

#[derive(Clone)]
struct Counted<'a> {
    bytes: &'a [u8],
    calls: Rc<Cell<usize>>,
    read: Rc<Cell<usize>>,
}
impl<'a> Counted<'a> {
    fn new(bytes: &'a [u8]) -> Self {
        Self {
            bytes,
            calls: Rc::new(Cell::new(0)),
            read: Rc::new(Cell::new(0)),
        }
    }
    fn reset(&self) {
        self.calls.set(0);
        self.read.set(0);
    }
    fn report(&self, name: &str) {
        println!("io\t{name}\t{}\t{}", self.calls.get(), self.read.get());
    }
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
    if std::env::var("FASTA_BENCH_FILTER").is_ok_and(|s| !name.contains(&s)) {
        return;
    }
    for _ in 0..20 {
        f();
    }
    let mut samples = [0.0f64; 5];
    for sample in &mut samples {
        let start = Instant::now();
        let mut n = 0;
        while start.elapsed() < Duration::from_millis(100) {
            for _ in 0..16 {
                f();
            }
            n += 16;
        }
        *sample = start.elapsed().as_nanos() as f64 / n as f64;
    }
    samples.sort_by(f64::total_cmp);
    println!("{name}\t{:.2}", samples[2]);
}
fn main() {
    for kind in ["long", "many", "duplicate"] {
        let mut text = Vec::new();
        let records = if kind == "long" { 1 } else { 8192 };
        let len = if kind == "long" { 2_000_000 } else { 64 };
        for i in 0..records {
            let name = if kind == "duplicate" { i % 16 } else { i };
            text.extend(format!(">seq{name} description\n").as_bytes());
            for line in vec![b'A'; len].chunks(60) {
                text.extend(line);
                text.push(b'\n');
            }
        }
        for chunk_len in [37, 1 << 20] {
            measure(&format!("fasta.builder.{kind}.chunk{chunk_len}"), || {
                let mut builder = FaiBuilder::new();
                for chunk in black_box(text.as_slice()).chunks(chunk_len) {
                    builder.push(chunk).unwrap();
                }
                let index = builder.finish().unwrap();
                assert_eq!(index.len(), if kind == "duplicate" { 16 } else { records });
                black_box(index);
            });
        }
    }
    for width in [60, 2_100_000] {
        let mut text = b">s\n".to_vec();
        let bases: Vec<_> = (0..2_000_000).map(|i| b"ACGT"[i % 4]).collect();
        for line in bases.chunks(width) {
            text.extend(line);
            text.push(b'\n');
        }
        let mut builder = FaiBuilder::new();
        builder.push(&text).unwrap();
        let fai = builder.finish().unwrap();
        let (gz, _, gzi) = compress_fasta(&text[..], Vec::new(), 1).unwrap();
        for (kind, bytes, index) in [
            ("plain", text.as_slice(), None),
            ("bgzf", gz.as_slice(), Some(gzi)),
        ] {
            let counted = Counted::new(bytes);
            let mut reader = FastaReader::new(counted.clone(), fai.clone(), index);
            for len in [64, 131072, 1_000_000] {
                let name = format!("fasta.{kind}.width{width}.len{len}");
                let region = Region::new("s", 123, 123 + len);
                let mut out = vec![0; len as usize];
                counted.reset();
                reader.sequence(region, &mut out).unwrap();
                counted.report(&name);
                assert_eq!(&out, &bases[123..123 + len as usize]);
                measure(&name, || {
                    reader
                        .sequence(black_box(region), black_box(&mut out))
                        .unwrap();
                    black_box(&out);
                });
            }
        }
    }
    for count in [64, 1024] {
        let seqs: Vec<_> = (0..count)
            .map(|i| common::TwoBitSeq {
                name: format!("chr{i}"),
                bases: vec![b'A'; 1024],
                n_blocks: Vec::new(),
                mask_blocks: Vec::new(),
            })
            .collect();
        let bytes = common::write_twobit(&seqs, 1);
        let source = Counted::new(&bytes);
        source.reset();
        let reader = TwoBitReader::open(source.clone(), false).unwrap();
        source.report(&format!("twobit.open.{count}"));
        measure(&format!("twobit.open.{count}"), || {
            black_box(TwoBitReader::open(source.clone(), false).unwrap());
        });
        measure(&format!("twobit.names.{count}"), || {
            black_box(reader.names().map(str::len).sum::<usize>());
        });
    }
    for all_n in [false, true] {
        let seq = common::TwoBitSeq {
            name: "s".into(),
            bases: vec![b'A'; 1_000_000],
            n_blocks: if all_n {
                vec![(0, 1_000_000)]
            } else {
                Vec::new()
            },
            mask_blocks: Vec::new(),
        };
        let bytes = common::write_twobit(&[seq], 0);
        let counted = Counted::new(&bytes);
        let mut reader = TwoBitReader::open(counted.clone(), true).unwrap();
        let mut out = vec![0; 131072];
        let region = Region::new("s", 123, 131195);
        counted.reset();
        reader.sequence(region, &mut out).unwrap();
        counted.report(&format!("twobit.sequence.allN{all_n}"));
        measure(&format!("twobit.sequence.allN{all_n}"), || {
            reader
                .sequence(black_box(region), black_box(&mut out))
                .unwrap();
            black_box(&out);
        });
    }
}
