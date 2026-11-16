use std::cell::Cell;
use std::io::{self, Cursor, Seek, SeekFrom, Write};
use std::rc::Rc;

use omicsweave_core::ContactUnit;
use omicsweave_hic::{Chrom, HicHeader, HicWriter, Record};

/// The fault occurs after a short successful write, and then clears, so an
/// output that resumes writing would otherwise appear to recover.
struct FaultWriter {
    bytes: Cursor<Vec<u8>>,
    remaining: Rc<Cell<Option<usize>>>,
    writes: Rc<Cell<usize>>,
}

impl Write for FaultWriter {
    fn write(&mut self, bytes: &[u8]) -> io::Result<usize> {
        self.writes.set(self.writes.get() + 1);
        match self.remaining.get() {
            Some(0) => {
                self.remaining.set(None);
                Err(io::Error::other("injected output failure"))
            }
            Some(left) => {
                let n = self.bytes.write(&bytes[..bytes.len().min(left)])?;
                self.remaining.set(Some(left - n));
                Ok(n)
            }
            None => self.bytes.write(bytes),
        }
    }

    fn flush(&mut self) -> io::Result<()> {
        Ok(())
    }
}

impl Seek for FaultWriter {
    fn seek(&mut self, position: SeekFrom) -> io::Result<u64> {
        self.bytes.seek(position)
    }
}

#[test]
fn partial_output_failure_prevents_successful_finalization() {
    for second in [false, true] {
        let remaining = Rc::new(Cell::new(None));
        let writes = Rc::new(Cell::new(0));
        let output = FaultWriter {
            bytes: Cursor::new(Vec::new()),
            remaining: Rc::clone(&remaining),
            writes: Rc::clone(&writes),
        };
        let header = HicHeader::new(
            "g",
            vec![Chrom {
                name: "s".into(),
                len: 100,
            }],
            vec![1, 2],
        );
        let mut writer = HicWriter::create(output, header, 1).unwrap();
        let record = [Record {
            x: 0,
            y: 0,
            value: 1.0,
        }];
        // The failure hits the first map, or one written after a good one.
        if second {
            writer
                .write_matrix(0, 0, ContactUnit::Bp, 2, &record)
                .unwrap();
        }
        remaining.set(Some(1));
        let result = writer.write_matrix(0, 0, ContactUnit::Bp, 1, &record);
        assert!(result.is_err());
        assert_eq!(remaining.get(), None, "the injected fault has cleared");
        let attempts = writes.get();
        assert!(
            writer
                .write_matrix(0, 0, ContactUnit::Bp, 1, &record)
                .is_err()
        );
        assert!(writer.finish().is_err());
        assert_eq!(
            writes.get(),
            attempts,
            "failed output must receive no more writes"
        );
    }
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum OutputOperation {
    Write,
    Seek,
    Flush,
}

struct FinalizeWriter {
    bytes: Cursor<Vec<u8>>,
    operations: Rc<std::cell::RefCell<Vec<OutputOperation>>>,
    fail_at: Rc<Cell<Option<usize>>>,
}

impl FinalizeWriter {
    fn operation(&self, operation: OutputOperation) -> io::Result<()> {
        let mut operations = self.operations.borrow_mut();
        let index = operations.len();
        operations.push(operation);
        if self.fail_at.get() == Some(index) {
            return Err(io::Error::other("injected finalization failure"));
        }
        Ok(())
    }
}

impl Write for FinalizeWriter {
    fn write(&mut self, bytes: &[u8]) -> io::Result<usize> {
        self.operation(OutputOperation::Write)?;
        self.bytes.write(bytes)
    }

    fn flush(&mut self) -> io::Result<()> {
        self.operation(OutputOperation::Flush)
    }
}

impl Seek for FinalizeWriter {
    fn seek(&mut self, position: SeekFrom) -> io::Result<u64> {
        self.operation(OutputOperation::Seek)?;
        self.bytes.seek(position)
    }
}

fn finalization_case(
    version: i32,
    with_sidecar: bool,
    fail_step: Option<usize>,
) -> Vec<OutputOperation> {
    let operations = Rc::new(std::cell::RefCell::new(Vec::new()));
    let fail_at = Rc::new(Cell::new(None));
    let out = FinalizeWriter {
        bytes: Cursor::new(Vec::new()),
        operations: operations.clone(),
        fail_at: fail_at.clone(),
    };
    let mut header = HicHeader::new(
        "g",
        vec![Chrom {
            name: "a".into(),
            len: 100,
        }],
        vec![1],
    );
    header.version = version;
    let mut writer = HicWriter::create(out, header, 1).unwrap();
    writer
        .write_matrix(
            0,
            0,
            ContactUnit::Bp,
            1,
            &[Record {
                x: 2,
                y: 3,
                value: 7.0,
            }],
        )
        .unwrap();
    let start = operations.borrow().len();
    fail_at.set(fail_step.map(|step| start + step));
    let result = if with_sidecar {
        writer.finish_with_sidecar().map(|(out, _)| out)
    } else {
        writer.finish()
    };
    if let Some(step) = fail_step {
        let error = result
            .err()
            .expect("finalization must report the injected output failure");
        assert!(
            error.to_string().contains("injected finalization failure"),
            "version={version}, sidecar={with_sidecar}, step={step}: {error}"
        );
        assert_eq!(
            operations.borrow().len(),
            start + step + 1,
            "no operation may follow an output error"
        );
    } else {
        let bytes = result.unwrap().bytes.into_inner();
        assert!(omicsweave_hic::HicReader::open(bytes.as_slice()).is_ok());
    }
    operations.borrow()[start..].to_vec()
}

#[test]
fn finalization_propagates_every_footer_write_backpatch_seek_and_flush_failure() {
    for version in [8, 9] {
        for with_sidecar in [false, true] {
            let operations = finalization_case(version, with_sidecar, None);
            assert!(operations.contains(&OutputOperation::Write));
            assert!(
                operations
                    .iter()
                    .filter(|&&op| op == OutputOperation::Seek)
                    .count()
                    >= 2
            );
            assert_eq!(operations.last(), Some(&OutputOperation::Flush));
            for step in 0..operations.len() {
                let attempted = finalization_case(version, with_sidecar, Some(step));
                assert_eq!(
                    attempted,
                    operations[..=step],
                    "version={version}, sidecar={with_sidecar}, step={step}"
                );
            }
        }
    }
}
