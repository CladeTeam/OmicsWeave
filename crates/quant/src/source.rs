//! Input BAMs: header checks, and the serial reader that inflates BGZF
//! blocks on the worker pool and cuts the records into batches of
//! consecutive records, checking their coordinate order as it goes.

use std::fs::File;
use std::io::{BufReader, Read};
use std::path::{Path, PathBuf};

use omicsweave_bam::{BamHeader, Scanner};
use omicsweave_core::exec::{Inline, Ordered};
use omicsweave_core::{Error, Result};

use crate::prelude::coord_key;

/// The fix for an unsorted input, named in errors.
pub const SORT_HINT: &str = "sort it first: ow convert in.bam sorted.bam --sort";

/// An input BAM, its header read.
#[derive(Debug, Clone)]
pub struct Input {
    pub path: PathBuf,
    /// References `(name, length)` by tid.
    pub refs: Vec<(String, u32)>,
    /// `ID`s of the `@RG` lines, in header order.
    pub read_groups: Vec<String>,
    /// Compressed size.
    pub len: u64,
}

/// The value of `key` (`SO`, `ID`, ...) on a header line's fields.
fn field<'a>(line: &'a str, key: &str) -> Option<&'a str> {
    line.split('\t')
        .skip(1)
        .find_map(|f| f.strip_prefix(key)?.strip_prefix(':'))
}

impl Input {
    /// Reads the header of `path` and checks that it declares coordinate
    /// order (`@HD SO:coordinate`).
    pub fn open(path: &Path) -> Result<Input> {
        let file = File::open(path).map_err(|e| Error::Io {
            kind: e.kind(),
            detail: format!("{}: {e}", path.display()),
        })?;
        let len = file.metadata()?.len();
        let mut scan = Scanner::new(
            BufReader::with_capacity(1 << 16, file),
            Box::new(Inline::new()),
        );
        let header = scan.header().map_err(|e| with_input(e, path))?;
        Input::check(path, &header, len)
    }

    fn check(path: &Path, header: &BamHeader, len: u64) -> Result<Input> {
        let text = String::from_utf8_lossy(header.text());
        let so = text
            .lines()
            .find(|l| l.starts_with("@HD"))
            .and_then(|l| field(l, "SO"));
        if so != Some("coordinate") {
            return Err(Error::invalid(
                "input",
                format!(
                    "{} is not sorted by coordinate ({}); {SORT_HINT}",
                    path.display(),
                    match so {
                        Some(s) => format!("@HD SO:{s}"),
                        None => "no @HD SO".into(),
                    }
                ),
            ));
        }
        let read_groups = text
            .lines()
            .filter(|l| l.starts_with("@RG"))
            .filter_map(|l| field(l, "ID").map(str::to_owned))
            .collect();
        Ok(Input {
            path: path.to_owned(),
            refs: header
                .references()
                .iter()
                .map(|r| (r.name.clone(), r.len))
                .collect(),
            read_groups,
            len,
        })
    }
}

fn with_input(e: Error, path: &Path) -> Error {
    match e {
        Error::Invalid { what, detail } => Error::Invalid {
            what,
            detail: format!("{}: {detail}", path.display()),
        },
        other => other,
    }
}

/// A run of consecutive records of one input: each with its 4-byte
/// `block_size`.
pub struct Batch {
    pub input: u32,
    pub data: Vec<u8>,
    pub records: u32,
    /// Number of the first record in its input.
    pub first: u64,
    /// [`coord_key`] of the last record.
    pub frontier: u64,
}

/// Reads the placed records of an input in batches; stops at the first
/// unplaced record (they sort last).
pub struct BatchReader<R: Read> {
    scan: Scanner<R>,
    input: u32,
    path: PathBuf,
    refs: Vec<String>,
    target: usize,
    next: u64,
    last: u64,
    done: bool,
}

impl BatchReader<BufReader<File>> {
    /// Opens `input`, inflating on `exec`, for batches of about `target`
    /// bytes.
    pub fn open(
        input: &Input,
        index: u32,
        exec: Box<dyn Ordered<Vec<u8>>>,
        target: usize,
    ) -> Result<BatchReader<BufReader<File>>> {
        let file = File::open(&input.path)?;
        let mut scan = Scanner::new(BufReader::with_capacity(1 << 20, file), exec);
        scan.header().map_err(|e| with_input(e, &input.path))?;
        Ok(BatchReader {
            scan,
            input: index,
            path: input.path.clone(),
            refs: input.refs.iter().map(|r| r.0.clone()).collect(),
            target: target.max(1 << 12),
            next: 0,
            last: 0,
            done: false,
        })
    }
}

impl<R: Read> BatchReader<R> {
    /// Compressed bytes read so far.
    pub fn bytes_read(&self) -> u64 {
        self.scan.compressed_read()
    }

    fn position(&self, key: u64) -> String {
        let tid = (key >> 32) as usize;
        let pos = key & 0xffff_ffff;
        match self.refs.get(tid) {
            Some(name) => format!("{name}:{pos}"),
            None => format!("reference {tid}:{pos}"),
        }
    }

    /// The next batch; `None` once every placed record has been read.
    pub fn next_batch(&mut self) -> Result<Option<Batch>> {
        if self.done {
            return Ok(None);
        }
        let mut data = Vec::with_capacity(self.target + (1 << 16));
        let first = self.next;
        let mut records = 0u32;
        while data.len() < self.target {
            let Some(range) = self
                .scan
                .next_record()
                .map_err(|e| with_input(e, &self.path))?
            else {
                self.done = true;
                break;
            };
            let rec = self.scan.slice(range);
            let tid = i32::from_le_bytes(rec[0..4].try_into().unwrap());
            if tid < 0 {
                self.done = true;
                break;
            }
            let pos = i32::from_le_bytes(rec[4..8].try_into().unwrap());
            let key = coord_key(tid, pos);
            if key < self.last {
                return Err(Error::invalid(
                    "input",
                    format!(
                        "{}: record {} at {} comes after {}: the BAM is not sorted by \
                         coordinate; {SORT_HINT}",
                        self.path.display(),
                        self.next + 1,
                        self.position(key),
                        self.position(self.last)
                    ),
                ));
            }
            self.last = key;
            data.extend_from_slice(&(rec.len() as u32).to_le_bytes());
            data.extend_from_slice(rec);
            records += 1;
            self.next += 1;
        }
        if records == 0 {
            return Ok(None);
        }
        Ok(Some(Batch {
            input: self.input,
            data,
            records,
            first,
            frontier: self.last,
        }))
    }
}
