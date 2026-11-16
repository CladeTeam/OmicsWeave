//! 2bit written by `TwoBitWriter`, piece by piece: byte-identical to the
//! test writer's `faToTwoBit` layout and read back by lib2bit.

use crate::support::{Rng, TwoBitSeq, lib2bit, write_twobit};
use omicsweave_core::{Locus, Region, SequenceRead};
use omicsweave_fasta::{
    Error, Fai, FastaWriter, TwoBitPiece, TwoBitReader, TwoBitRecord, TwoBitWriter, encode_twobit,
};
use proptest::prelude::*;
use std::{
    cell::Cell,
    io::{self, Cursor, Seek, SeekFrom, Write},
    rc::Rc,
};

fn random_seq(rng: &mut Rng) -> Vec<u8> {
    let len = if rng.chance(5) {
        0
    } else {
        rng.range(1, 6_000) as usize
    };
    let mut seq = Vec::with_capacity(len);
    while seq.len() < len {
        let run = (rng.range(1, 300) as usize).min(len - seq.len());
        let alphabet: &[u8] = match rng.range(0, 6) {
            0 => b"NNNN",
            1 => b"nnnn",
            2 => b"acgt",
            3 => b"ACGTRYKMn-",
            _ => b"ACGT",
        };
        seq.extend((0..run).map(|_| *rng.pick(alphabet)));
    }
    seq
}

fn runs(seq: &[u8], hit: impl Fn(u8) -> bool) -> Vec<(u32, u32)> {
    let mut out = Vec::new();
    let mut i = 0;
    while i < seq.len() {
        if hit(seq[i]) {
            let s = i;
            while i < seq.len() && hit(seq[i]) {
                i += 1;
            }
            out.push((s as u32, (i - s) as u32));
        } else {
            i += 1;
        }
    }
    out
}

fn is_acgt(c: u8) -> bool {
    b"ACGTacgt".contains(&c)
}

proptest! {
    #![proptest_config(ProptestConfig { cases: 48, failure_persistence: Some(Box::new(proptest::test_runner::FileFailurePersistence::Direct(concat!(env!("CARGO_MANIFEST_DIR"), "/tests/cases/writer.proptest-regressions")))), ..ProptestConfig::default() })]

    #[test]
    fn matches_the_reference_layout_and_lib2bit(seed in any::<u64>(), version in 0u32..=1) {
        let mut rng = Rng::new(seed);
        let seqs: Vec<(String, Vec<u8>)> = (0..rng.range(1, 6))
            .map(|i| (format!("s{i}"), random_seq(&mut rng)))
            .collect();
        let piece = 4 * rng.range(1, 400) as usize;

        let names = seqs.iter().map(|(n, _)| n.clone()).collect();
        let mut w = TwoBitWriter::create(Cursor::new(Vec::new()), names, version).unwrap();
        for (name, seq) in &seqs {
            let mut rec = TwoBitRecord::new();
            let mut packed = Vec::new();
            for (i, chunk) in seq.chunks(piece).enumerate() {
                let p = TwoBitPiece::encode((i * piece) as u64, chunk).unwrap();
                rec.push(&p).unwrap();
                packed.extend_from_slice(&p.packed);
            }
            w.add(name, &rec, &mut &packed[..]).unwrap();
        }
        let ours = w.finish().unwrap().into_inner();

        let reference: Vec<TwoBitSeq> = seqs
            .iter()
            .map(|(name, seq)| TwoBitSeq {
                name: name.clone(),
                bases: seq.iter().map(|&c| if is_acgt(c) { c.to_ascii_uppercase() } else { b'T' }).collect(),
                n_blocks: runs(seq, |c| !is_acgt(c)),
                mask_blocks: runs(seq, |c| c.is_ascii_lowercase()),
            })
            .collect();
        prop_assert_eq!(&ours, &write_twobit(&reference, version));

        let dir = tempfile::tempdir().unwrap();
        let path = dir.path().join("t.2bit");
        std::fs::write(&path, &ours).unwrap();
        for mask in [false, true] {
            let theirs = lib2bit::TwoBit::open(&path, mask).unwrap();
            for (name, seq) in &seqs {
                prop_assert_eq!(theirs.len(name), seq.len() as u32);
                if seq.is_empty() {
                    continue;
                }
                let got = theirs.sequence(name, 0, seq.len() as u32).unwrap();
                for (i, (&g, &c)) in got.iter().zip(seq).enumerate() {
                    if is_acgt(c) {
                        let want = if mask { c } else { c.to_ascii_uppercase() };
                        prop_assert_eq!(g, want, "{}:{}", name, i);
                    } else {
                        prop_assert_eq!(g.to_ascii_uppercase(), b'N', "{}:{}", name, i);
                    }
                }
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

#[test]
fn failed_fasta_sink_does_not_commit_index_or_accept_more_records() {
    let (sink, remaining, calls) = sink();
    let mut writer = FastaWriter::new(sink, 0);
    remaining.set(32);
    assert!(writer.add("s", &vec![b'A'; 70_000]).is_err());
    assert!(writer.fai().is_empty());
    let failed_calls = calls.get();
    remaining.set(usize::MAX);
    assert!(writer.add("t", b"ACGT").is_err());
    assert!(writer.finish().is_err());
    assert_eq!(calls.get(), failed_calls);
}

#[test]
fn failed_twobit_sink_cannot_be_finished_or_reused() {
    let (sink, remaining, calls) = sink();
    let mut writer = TwoBitWriter::create(sink, vec!["s".into()], 0).unwrap();
    let (record, packed) = encode_twobit(&vec![b'A'; 1024]).unwrap();
    remaining.set(32);
    assert!(writer.add("s", &record, &mut packed.as_slice()).is_err());
    let failed_calls = calls.get();
    remaining.set(usize::MAX);
    assert!(writer.add("s", &record, &mut packed.as_slice()).is_err());
    assert!(writer.finish().is_err());
    assert_eq!(calls.get(), failed_calls);
}

struct FailAfter {
    bytes: Cursor<Vec<u8>>,
    remaining: Rc<Cell<usize>>,
    calls: Rc<Cell<usize>>,
}
impl Write for FailAfter {
    fn write(&mut self, bytes: &[u8]) -> io::Result<usize> {
        self.calls.set(self.calls.get() + 1);
        let n = bytes.len().min(self.remaining.get());
        if n == 0 {
            return Err(io::Error::other("injected sink failure"));
        }
        self.remaining.set(self.remaining.get() - n);
        self.bytes.write(&bytes[..n])
    }
    fn flush(&mut self) -> io::Result<()> {
        Ok(())
    }
}
impl Seek for FailAfter {
    fn seek(&mut self, pos: SeekFrom) -> io::Result<u64> {
        self.bytes.seek(pos)
    }
}
fn sink() -> (FailAfter, Rc<Cell<usize>>, Rc<Cell<usize>>) {
    let remaining = Rc::new(Cell::new(usize::MAX));
    let calls = Rc::new(Cell::new(0));
    (
        FailAfter {
            bytes: Cursor::new(Vec::new()),
            remaining: remaining.clone(),
            calls: calls.clone(),
        },
        remaining,
        calls,
    )
}

#[test]
fn writer_rejects_records_that_would_not_index_as_written() {
    let mut writer = FastaWriter::new(Vec::new(), 4);
    writer.add("ok", b"ACGTA").unwrap();
    let rejected: [(&str, &[u8]); 13] = [
        ("", b"ACGT"),
        ("a b", b"ACGT"),
        ("a\tb", b"ACGT"),
        ("a\rb", b"ACGT"),
        ("a\0b", b"ACGT"),
        ("ok", b"ACGT"),
        ("e", b""),
        ("g", b"AC GT"),
        ("g", b"AC\nGT"),
        ("g", b"AC>GT"),
        ("g", b">ACGT"),
        ("g", b"AC\x7fGT"),
        ("g", "ACéGT".as_bytes()),
    ];
    for (name, seq) in rejected {
        assert!(
            matches!(writer.add(name, seq), Err(Error::Invalid { .. })),
            "{name:?} {seq:?}"
        );
    }
    writer.add("染色体", b"N").unwrap();
    let (bytes, fai, gzi) = writer.finish().unwrap();
    assert_eq!(bytes, ">ok\nACGT\nA\n>染色体\nN\n".as_bytes());
    assert_eq!(fai, Fai::build(&bytes[..]).unwrap());
    assert_eq!(gzi, None);
    assert!(FastaWriter::new(Vec::new(), 60).finish().is_err());
    assert!(FastaWriter::bgzf(Vec::new(), 60, 13).is_err());
}

#[test]
fn named_records_advance_only_after_success_and_consume_exact_packed_bytes() {
    for version in [0, 1] {
        for names in [
            vec!["".into()],
            vec!["s".into(), "s".into()],
            vec!["x".repeat(256)],
        ] {
            assert!(TwoBitWriter::create(Cursor::new(Vec::new()), names, version).is_err());
        }
        let (record, packed) = encode_twobit(b"AcGTN").unwrap();
        let mut writer = TwoBitWriter::create(
            Cursor::new(Vec::new()),
            vec!["s".into(), "empty".into()],
            version,
        )
        .unwrap();
        assert_eq!(writer.next_name(), Some("s"));
        assert!(
            writer
                .add("empty", &record, &mut packed.as_slice())
                .is_err()
        );
        assert_eq!(writer.next_name(), Some("s"));
        let mut extra = packed.clone();
        extra.extend_from_slice(b"sentinel");
        let mut input = Cursor::new(extra);
        writer.add("s", &record, &mut input).unwrap();
        assert_eq!(input.position(), packed.len() as u64);
        assert_eq!(writer.next_name(), Some("empty"));
        let (empty, packed_empty) = encode_twobit(b"").unwrap();
        writer
            .add("empty", &empty, &mut packed_empty.as_slice())
            .unwrap();
        assert_eq!(writer.next_name(), None);
        assert!(
            writer
                .add("extra", &empty, &mut packed_empty.as_slice())
                .is_err()
        );
        let bytes = writer.finish().unwrap().into_inner();
        let mut reader = TwoBitReader::open(bytes, true).unwrap();
        assert_eq!(reader.names().collect::<Vec<_>>(), ["s", "empty"]);
        assert_eq!(reader.chrom_len("empty").unwrap(), 0);
        let mut out = [0; 5];
        reader
            .sequence(&Locus::from(Region::new("s", 0, 5)), &mut out)
            .unwrap();
        assert_eq!(&out, b"AcGTN");
        reader
            .sequence(&Locus::from(Region::new("empty", 0, 0)), &mut [])
            .unwrap();
        let incomplete =
            TwoBitWriter::create(Cursor::new(Vec::new()), vec!["s".into()], version).unwrap();
        assert!(incomplete.finish().is_err());
    }
}
