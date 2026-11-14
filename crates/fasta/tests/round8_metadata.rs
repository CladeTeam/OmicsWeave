//! Metadata parsing and sink failures at the public writer boundary.
use omicsweave_fasta::{Fai, FastaWriter, TwoBitWriter, encode_twobit};
use std::{
    cell::Cell,
    io::{self, Cursor, Seek, SeekFrom, Write},
    rc::Rc,
};

#[test]
fn parsed_duplicate_names_keep_first_rows_and_still_validate_every_row() {
    let text = b"s\t12\t3\t6\t7\ns\t99\t15\t9\t10\nt\t8\t42\t4\t5\n";
    let fai = Fai::parse(text).unwrap();
    assert_eq!(fai.len(), 2);
    assert_eq!(fai.get("s").unwrap().len, 12);
    assert_eq!(fai.to_bytes(), b"s\t12\t3\t6\t7\nt\t8\t42\t4\t5\n");
    for tail in [
        b"s\tnan\t15\t9\t10\n".as_slice(),
        b"s\t99\t15\t10\t9\n",
        b"s\t99\t15\t9\t10\textra\n",
    ] {
        let mut invalid = text.to_vec();
        invalid.extend_from_slice(tail);
        assert!(Fai::parse(&invalid).is_err());
    }
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
