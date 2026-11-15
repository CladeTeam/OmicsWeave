//! Physical-member windows, independent metadata and true htslib results.
mod common;

use std::cell::{Cell, RefCell};
use std::io::Write;

use omicsweave_core::{Locus, ReadAt, Region, Result, SequenceRead};
use omicsweave_fasta::{BgzfWriter, EOF_BLOCK, Fai, FaiEntry, FastaReader, Gzi, MAX_BLOCK_LEN};

struct Source {
    bytes: RefCell<Vec<u8>>,
    calls: Cell<usize>,
    total: Cell<usize>,
    largest: Cell<usize>,
}

impl Source {
    fn new(bytes: Vec<u8>) -> Self {
        Self {
            bytes: RefCell::new(bytes),
            calls: Cell::new(0),
            total: Cell::new(0),
            largest: Cell::new(0),
        }
    }
}

impl ReadAt for Source {
    fn len(&self) -> u64 {
        self.bytes.borrow().len() as u64
    }

    fn read_exact_at(&self, out: &mut [u8], offset: u64) -> Result<()> {
        self.calls.set(self.calls.get() + 1);
        self.total.set(self.total.get() + out.len());
        self.largest.set(self.largest.get().max(out.len()));
        self.bytes.borrow().as_slice().read_exact_at(out, offset)
    }
}

fn fixture(empty_members: usize, bases: &[u8]) -> (Vec<u8>, Fai, Gzi) {
    let mut bytes = Vec::new();
    for _ in 0..empty_members {
        bytes.extend_from_slice(&EOF_BLOCK);
    }
    let mut writer = BgzfWriter::new(bytes);
    writer.write_all(b">s\n").unwrap();
    writer.write_all(bases).unwrap();
    writer.write_all(b"\n").unwrap();
    let (bytes, _) = writer.finish().unwrap();
    let fai = Fai::from_entries([FaiEntry {
        name: "s".into(),
        len: bases.len() as u32,
        offset: 3,
        line_bases: bases.len() as u32,
        line_width: bases.len() as u32 + 1,
    }]);
    // htslib leaves the first non-empty member implicit, even after empties.
    let gzi = Gzi::parse(&0u64.to_le_bytes()).unwrap();
    (bytes, fai, gzi)
}

#[test]
fn leading_empty_members_read_bounded_windows_without_repeated_suffixes() {
    let bases = b"ACGTacgtNNNN";
    for count in [1, 2, 17, 2339, 2340, 2341, 2342, 3000, 5000] {
        let (bytes, fai, gzi) = fixture(count, bases);
        let physical = bytes.len();
        let source = Source::new(bytes);
        let mut reader = FastaReader::new(&source, fai, Some(gzi));
        let mut guarded = [0xA5; 14];
        reader
            .sequence(&Locus::from(Region::new("s", 0, 12)), &mut guarded[1..13])
            .unwrap();
        assert_eq!(&guarded, b"\xA5ACGTacgtNNNN\xA5", "count={count}");
        assert!(source.largest.get() <= MAX_BLOCK_LEN);
        assert!(
            source.calls.get() <= 3 + physical / (MAX_BLOCK_LEN - EOF_BLOCK.len()),
            "count={count}"
        );
        assert!(
            source.total.get() <= 2 * physical + MAX_BLOCK_LEN,
            "count={count}, read={}",
            source.total.get()
        );
    }
}

#[test]
fn partial_header_and_data_windows_match_htslib() {
    let bases: Vec<u8> = (0u32..4096)
        .map(|i| {
            let v = i.wrapping_mul(0x9E37_79B9).rotate_left(13) ^ i.wrapping_mul(31);
            b"ACGT"[(v as usize >> 16) & 3]
        })
        .collect();
    common::htslib::quiet();
    let dir = tempfile::tempdir().unwrap();
    for count in [2338, 2339, 2340, 2341, 2342] {
        let (bytes, fai, gzi) = fixture(count, &bases);
        let path = dir.path().join(format!("leading-{count}.fa.gz"));
        std::fs::write(&path, &bytes).unwrap();
        assert!(common::htslib::build_index(&path));
        assert_eq!(
            std::fs::read(omicsweave_fasta::fai_path(&path)).unwrap(),
            fai.to_bytes()
        );
        assert_eq!(
            std::fs::read(omicsweave_fasta::gzi_path(&path)).unwrap(),
            gzi.to_bytes()
        );
        let oracle = common::htslib::Faidx::load(&path).unwrap();
        let source = Source::new(bytes);
        let mut reader = FastaReader::new(&source, fai, Some(gzi));
        for (start, end) in [(0, 4096), (1, 2), (27, 100), (4000, 4096)] {
            let mut out = vec![0; (end - start) as usize];
            reader
                .sequence(&Locus::from(Region::new("s", start, end)), &mut out)
                .unwrap();
            assert_eq!(out, bases[start as usize..end as usize]);
            assert_eq!(Some(out), oracle.fetch("s", start, end));
        }
    }
}

#[test]
fn failed_empty_crc_and_truncated_data_allow_repaired_source_queries() {
    let bases = b"ACGTacgtNNNN";
    let count = 2342;
    let (bytes, fai, gzi) = fixture(count, bases);
    let good = bytes.clone();
    let source = Source::new(bytes);
    let mut reader = FastaReader::new(&source, fai, Some(gzi));
    let mut out = [0xA5; 14];
    let crc = 2340 * EOF_BLOCK.len() + EOF_BLOCK.len() - 8;
    source.bytes.borrow_mut()[crc] ^= 1;
    assert!(
        reader
            .sequence(&Locus::from(Region::new("s", 0, 12)), &mut out[1..13])
            .is_err()
    );
    assert_eq!(out, [0xA5; 14]);
    source.bytes.borrow_mut()[crc] ^= 1;
    reader
        .sequence(&Locus::from(Region::new("s", 0, 12)), &mut out[1..13])
        .unwrap();
    assert_eq!(&out, b"\xA5ACGTacgtNNNN\xA5");

    source
        .bytes
        .borrow_mut()
        .truncate(count * EOF_BLOCK.len() + 20);
    assert!(
        reader
            .sequence(&Locus::from(Region::new("s", 0, 12)), &mut out[1..13])
            .is_err()
    );
    *source.bytes.borrow_mut() = good;
    reader
        .sequence(&Locus::from(Region::new("s", 0, 12)), &mut out[1..13])
        .unwrap();
    assert_eq!(&out, b"\xA5ACGTacgtNNNN\xA5");
}
