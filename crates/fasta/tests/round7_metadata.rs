//! Metadata round trips and transactional rejection of malformed public pieces.
mod common;

use common::{TwoBitSeq, lib2bit, write_twobit};
use omicsweave_core::{Locus, ReadAt, Region, Result, SequenceRead};
use omicsweave_fasta::{TwoBitPiece, TwoBitReader, TwoBitRecord, TwoBitWriter, encode_twobit};
use std::{cell::Cell, io::Cursor};

struct Counted<'a> {
    bytes: &'a [u8],
    calls: Cell<usize>,
    read: Cell<usize>,
}
impl ReadAt for Counted<'_> {
    fn len(&self) -> u64 {
        self.bytes.len() as u64
    }
    fn read_exact_at(&self, out: &mut [u8], offset: u64) -> Result<()> {
        self.calls.set(self.calls.get() + 1);
        self.read.set(self.read.get() + out.len());
        self.bytes.read_exact_at(out, offset)
    }
}

#[test]
fn minimal_record_headers_match_lib2bit_across_directory_read_boundary() {
    for version in [0, 1] {
        let seqs: Vec<_> = (0..8)
            .map(|i| TwoBitSeq {
                name: format!("s{i}"),
                bases: b"TCAG".repeat(4096),
                n_blocks: Vec::new(),
                mask_blocks: Vec::new(),
            })
            .collect();
        let bytes = write_twobit(&seqs, version);
        let dir = tempfile::tempdir().unwrap();
        let path = dir.path().join("records.2bit");
        std::fs::write(&path, &bytes).unwrap();
        let reference = lib2bit::TwoBit::open(&path, true).unwrap();
        let source = Counted {
            bytes: &bytes,
            calls: Cell::new(0),
            read: Cell::new(0),
        };
        let mut reader = TwoBitReader::open(&source, true).unwrap();
        // The directory bound includes s0; each other record needs one
        // complete minimal header, with no extra packed-base prefetch.
        assert_eq!(source.calls.get(), 2 + seqs.len() - 1);
        assert_eq!(
            source.read.get(),
            16 + 264 * seqs.len() + 16 * (seqs.len() - 1)
        );
        for seq in &seqs {
            for (start, end) in [(0, 17), (123, 1024), (16380, 16384)] {
                let mut out = vec![0; (end - start) as usize];
                reader
                    .sequence(&Locus::from(Region::new(&seq.name, start, end)), &mut out)
                    .unwrap();
                assert_eq!(out, reference.sequence(&seq.name, start, end).unwrap());
            }
        }
    }
}

#[test]
fn invalid_public_pieces_leave_the_record_unchanged() {
    let prefix = TwoBitPiece::encode(0, b"acNN").unwrap();
    let valid = TwoBitPiece::encode(4, b"nnAC").unwrap();
    let mut record = TwoBitRecord::new();
    record.push(&prefix).unwrap();
    let original = record.clone();
    let mut invalid = Vec::new();
    let mut piece = valid.clone();
    piece.packed.clear();
    invalid.push(piece);
    for runs in [
        vec![(3, 5)],
        vec![(4, 9)],
        vec![(4, 4)],
        vec![(6, 8), (4, 6)],
        vec![(4, 7), (6, 8)],
    ] {
        for mask in [false, true] {
            let mut piece = valid.clone();
            if mask {
                piece.mask_runs = runs.clone();
            } else {
                piece.n_runs = runs.clone();
            }
            invalid.push(piece);
        }
    }
    for piece in invalid {
        assert!(record.push(&piece).is_err());
        assert_eq!(record, original);
    }
    record.push(&valid).unwrap();
    let mut writer = TwoBitWriter::create(Cursor::new(Vec::new()), vec!["s".into()], 0).unwrap();
    let packed = [prefix.packed[0], valid.packed[0]];
    writer.add("s", &record, &mut packed.as_slice()).unwrap();
    let bytes = writer.finish().unwrap().into_inner();
    let mut reader = TwoBitReader::open(bytes, true).unwrap();
    let mut out = [0; 8];
    reader
        .sequence(&Locus::from(Region::new("s", 0, 8)), &mut out)
        .unwrap();
    assert_eq!(&out, b"acNNNNAC");
    assert!(TwoBitPiece::encode(u64::MAX - 3, b"ACGT").is_err());
}

#[test]
fn short_packed_input_prevents_reusing_a_partially_written_sink() {
    let (record, packed) = encode_twobit(b"ACGTACGT").unwrap();
    let mut writer = TwoBitWriter::create(Cursor::new(Vec::new()), vec!["s".into()], 0).unwrap();
    assert!(writer.add("s", &record, &mut &packed[..1]).is_err());
    assert!(writer.add("s", &record, &mut packed.as_slice()).is_err());
    assert!(writer.finish().is_err());
}

#[test]
fn rejected_record_name_does_not_poison_the_writer() {
    let (record, packed) = encode_twobit(b"ACGT").unwrap();
    let mut writer = TwoBitWriter::create(Cursor::new(Vec::new()), vec!["s".into()], 0).unwrap();
    assert!(
        writer
            .add("wrong", &record, &mut packed.as_slice())
            .is_err()
    );
    writer.add("s", &record, &mut packed.as_slice()).unwrap();
    assert!(writer.finish().is_ok());
}
