//! Leading empty BGZF members: warm-memory queries and physical I/O counts.
use std::cell::Cell;
use std::hint::black_box;
use std::io::Write;
use std::time::{Duration, Instant};

use omicsweave_core::{ReadAt, Region, Result, SequenceRead};
use omicsweave_fasta::{BgzfWriter, EOF_BLOCK, Fai, FaiEntry, FastaReader, Gzi};

struct Source {
    bytes: Vec<u8>,
    calls: Cell<usize>,
    total: Cell<usize>,
}

impl ReadAt for Source {
    fn len(&self) -> u64 {
        self.bytes.len() as u64
    }
    fn read_exact_at(&self, out: &mut [u8], offset: u64) -> Result<()> {
        self.calls.set(self.calls.get() + 1);
        self.total.set(self.total.get() + out.len());
        self.bytes.as_slice().read_exact_at(out, offset)
    }
}

fn main() {
    for count in [0, 1, 64, 1024, 4096] {
        let mut bytes = Vec::new();
        for _ in 0..count {
            bytes.extend_from_slice(&EOF_BLOCK);
        }
        let mut writer = BgzfWriter::new(bytes);
        writer.write_all(b">s\nACGTacgtNNNN\n").unwrap();
        let (bytes, _) = writer.finish().unwrap();
        let physical = bytes.len();
        let source = Source {
            bytes,
            calls: Cell::new(0),
            total: Cell::new(0),
        };
        let fai = Fai::from_entries([FaiEntry {
            name: "s".into(),
            len: 12,
            offset: 3,
            line_bases: 12,
            line_width: 13,
        }]);
        let gzi = Gzi::parse(&0u64.to_le_bytes()).unwrap();
        let mut reader = FastaReader::new(&source, fai, Some(gzi));
        let region = Region::new("s", 0, 12);
        let mut out = [0; 12];
        reader.sequence(region, &mut out).unwrap();
        assert_eq!(&out, b"ACGTacgtNNNN");
        println!(
            "io\tempty-members.{count}\tphysical={physical}\tcalls={}\tbytes={}\tdigest={:016x}",
            source.calls.get(),
            source.total.get(),
            out.iter()
                .fold(0u64, |h, &b| h.wrapping_mul(131).wrapping_add(b as u64))
        );
        let mut samples = [0.0; 5];
        for sample in &mut samples {
            let start = Instant::now();
            let mut n = 0;
            while start.elapsed() < Duration::from_millis(75) {
                reader
                    .sequence(black_box(region), black_box(&mut out))
                    .unwrap();
                black_box(&out);
                n += 1;
            }
            *sample = start.elapsed().as_nanos() as f64 / n as f64;
        }
        samples.sort_by(f64::total_cmp);
        assert_eq!(&out, b"ACGTacgtNNNN");
        println!("empty-members.{count}\t{:.2}", samples[2]);
    }
}
