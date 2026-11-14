//! Building a `.fai` from FASTA text fed in chunks of any size.
//!
//! The rules follow htslib's `fai_build_core` (FASTA part), so the index is
//! byte-identical to `samtools faidx`:
//!
//! - A name runs from the first non-space byte after `>` to the next
//!   whitespace; the rest of the header line is ignored.
//! - All lines of a sequence must have the same length, except the last,
//!   which may be shorter; a shorter line must be followed by a header or a
//!   blank line.
//! - Sequences without any sequence line are dropped, except that a file whose
//!   last record has none is rejected.
//! - A repeated name keeps its first sequence.
//! - Names must be UTF-8, matching the string keys used by readers.
//! - Bases are counted as printable non-space bytes (`isgraph`); line widths
//!   count every byte including the terminator, so `\r\n` files work.

use crate::error::{Error, Result};
use crate::fai::{Fai, FaiEntry};

#[inline]
fn is_space(c: u8) -> bool {
    matches!(c, b' ' | b'\t' | b'\n' | 0x0B | 0x0C | b'\r')
}

#[inline]
fn is_graph(c: u8) -> bool {
    (0x21..=0x7E).contains(&c)
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum State {
    /// Between records: expecting `>` or a blank line.
    Between,
    /// After a `\r` between records: expecting `\n`.
    BetweenCr,
    /// Right after `>`; the previous record is stored on the next byte.
    HeaderStart,
    Name,
    /// The rest of a header line after its name.
    HeaderRest,
    /// At the start of a line within a sequence.
    LineStart,
    /// Within a sequence line.
    Line,
}

/// Incremental `.fai` builder: [`push`](FaiBuilder::push) the FASTA text in
/// order, then [`finish`](FaiBuilder::finish).
pub struct FaiBuilder {
    state: State,
    /// Bytes consumed so far.
    pos: u64,
    line: u64,
    name: Vec<u8>,
    has_line: bool,
    seq_offset: u64,
    seq_len: u64,
    line_width: u64,
    line_bases: u64,
    width: u64,
    bases: u64,
    index: Fai,
}

impl Default for FaiBuilder {
    fn default() -> Self {
        Self::new()
    }
}

impl FaiBuilder {
    pub fn new() -> FaiBuilder {
        FaiBuilder {
            state: State::Between,
            pos: 0,
            line: 1,
            name: Vec::new(),
            has_line: false,
            seq_offset: 0,
            seq_len: 0,
            line_width: 0,
            line_bases: 0,
            width: 0,
            bases: 0,
            index: Fai::default(),
        }
    }

    fn error(&self, detail: impl Into<String>) -> Error {
        Error::Format {
            line: self.line,
            detail: detail.into(),
        }
    }

    fn store(&mut self) -> Result<()> {
        let name = std::str::from_utf8(&self.name)
            .map_err(|_| self.error("sequence name is not UTF-8"))?;
        let small = |v: u64| {
            u32::try_from(v).map_err(|_| self.error(format!("sequence '{name}' exceeds 2^32 - 1")))
        };
        let (len, line_bases, line_width) = (
            small(self.seq_len)?,
            small(self.line_bases)?,
            small(self.line_width)?,
        );
        if self.index.get(name).is_none() {
            self.index.insert_missing(FaiEntry {
                name: name.to_owned(),
                len,
                offset: self.seq_offset,
                line_bases,
                line_width,
            });
        }
        self.has_line = false;
        Ok(())
    }

    /// The header line ended; the sequence starts at byte `offset`.
    fn begin_sequence(&mut self, offset: u64) {
        self.state = State::LineStart;
        self.seq_offset = offset;
        self.seq_len = 0;
        self.line_width = 0;
        self.line_bases = 0;
        self.line += 1;
    }

    fn end_line(&mut self) -> Result<()> {
        self.width += 1;
        self.seq_len += self.bases;
        if self.line_width == 0 {
            self.line_width = self.width;
            self.line_bases = self.bases;
            self.state = State::LineStart;
        } else if self.line_width > self.width {
            self.state = State::Between;
        } else if self.line_width < self.width {
            let name = String::from_utf8_lossy(&self.name).into_owned();
            return Err(self.error(format!("different line length in sequence '{name}'")));
        } else {
            self.state = State::LineStart;
        }
        self.line += 1;
        Ok(())
    }

    /// Consumes the next chunk of FASTA text.
    pub fn push(&mut self, chunk: &[u8]) -> Result<()> {
        let mut i = 0;
        while i < chunk.len() {
            if self.state == State::Line {
                let rest = &chunk[i..];
                let stop = memchr::memchr(b'\n', rest);
                let body = &rest[..stop.unwrap_or(rest.len())];
                self.width += body.len() as u64;
                self.bases += body.iter().filter(|&&c| is_graph(c)).count() as u64;
                i += body.len();
                self.pos += body.len() as u64;
                if stop.is_some() {
                    i += 1;
                    self.pos += 1;
                    self.end_line()?;
                }
                continue;
            }
            let c = chunk[i];
            i += 1;
            self.pos += 1;
            self.byte(c)?;
        }
        Ok(())
    }

    fn byte(&mut self, c: u8) -> Result<()> {
        match self.state {
            State::Between => match c {
                b'>' => self.state = State::HeaderStart,
                b'@' => return Err(self.error("FASTQ is not supported")),
                b'\r' => self.state = State::BetweenCr,
                b'\n' => self.line += 1,
                _ => return Err(self.error(format!("unexpected byte {c:#04x}"))),
            },
            State::BetweenCr => {
                if c != b'\n' {
                    return Err(self.error("carriage return not followed by new line"));
                }
                self.line += 1;
                self.state = State::Between;
            }
            State::HeaderStart => {
                if self.has_line {
                    self.store()?;
                }
                self.name.clear();
                self.state = State::Name;
                self.name_byte(c);
            }
            State::Name => self.name_byte(c),
            State::HeaderRest => {
                if c == b'\n' {
                    self.begin_sequence(self.pos);
                }
            }
            State::LineStart => match c {
                b'\n' => {
                    self.state = State::Between;
                    self.line += 1;
                }
                b'>' => self.state = State::HeaderStart,
                _ => {
                    self.has_line = true;
                    self.width = 0;
                    self.bases = 0;
                    self.state = State::Line;
                    self.width += 1;
                    self.bases += is_graph(c) as u64;
                }
            },
            State::Line => unreachable!("lines are consumed in bulk"),
        }
        Ok(())
    }

    fn name_byte(&mut self, c: u8) {
        if !is_space(c) {
            self.name.push(c);
        } else if c == b'\n' {
            self.begin_sequence(self.pos);
        } else if !self.name.is_empty() {
            self.state = State::HeaderRest;
        }
    }

    /// Ends the text and returns the index.
    pub fn finish(mut self) -> Result<Fai> {
        match self.state {
            State::Line => self.end_line()?,
            State::Name | State::HeaderRest => {
                let name = String::from_utf8_lossy(&self.name).into_owned();
                return Err(self.error(format!("the last entry '{name}' has no sequence")));
            }
            State::BetweenCr => return Err(self.error("carriage return not followed by new line")),
            State::Between | State::HeaderStart | State::LineStart => {}
        }
        if !self.has_line {
            return Err(self.error("file truncated"));
        }
        self.store()?;
        Ok(self.index)
    }
}
