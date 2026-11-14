//! Writing FASTA with its indexes: the `.fai` entries follow from the layout
//! the writer chooses, and the `.gzi` from the blocks it compresses, so
//! nothing is parsed back.

use std::io::Write;

use crate::bgzf::BgzfWriter;
use crate::error::{Result, invalid};
use crate::fai::{Fai, FaiEntry};
use crate::gzi::Gzi;

/// Text gathered before it is handed to the sink; larger pieces go through
/// directly.
const BUFFER_LEN: usize = 1 << 16;

enum Sink<W: Write> {
    Plain(W),
    Bgzf(BgzfWriter<W>),
}

impl<W: Write> Sink<W> {
    fn write_all(&mut self, bytes: &[u8]) -> Result<()> {
        match self {
            Sink::Plain(out) => out.write_all(bytes)?,
            Sink::Bgzf(out) => out.write_all(bytes)?,
        }
        Ok(())
    }
}

/// Writes FASTA records, plain or bgzip, and returns their `.fai` (and
/// `.gzi`) when finished. Each record is `>name`, a new line, then the
/// sequence in lines of `line_width` bases, each ended by `\n`.
///
/// Records are checked before anything of them is written, so the file
/// indexes, by [`crate::build_index`] or `samtools faidx`, to exactly the
/// entries returned; a rejected record leaves the output unchanged.
pub struct FastaWriter<W: Write> {
    sink: Sink<W>,
    line_width: usize,
    buf: Vec<u8>,
    /// Text bytes produced so far, buffered ones included.
    pos: u64,
    fai: Fai,
    broken: bool,
}

impl<W: Write> FastaWriter<W> {
    /// Plain FASTA; `line_width` bases per line (0 = each sequence on one
    /// line).
    pub fn new(out: W, line_width: usize) -> FastaWriter<W> {
        FastaWriter::with_sink(Sink::Plain(out), line_width)
    }

    /// BGZF-compressed FASTA at libdeflate `level` (0..=12), in blocks as
    /// bgzip writes them ([`BgzfWriter`]).
    pub fn bgzf(out: W, line_width: usize, level: i32) -> Result<FastaWriter<W>> {
        Ok(FastaWriter::with_sink(
            Sink::Bgzf(BgzfWriter::with_level(out, level)?),
            line_width,
        ))
    }

    fn with_sink(sink: Sink<W>, line_width: usize) -> FastaWriter<W> {
        FastaWriter {
            sink,
            line_width,
            buf: Vec::with_capacity(BUFFER_LEN),
            pos: 0,
            fai: Fai::default(),
            broken: false,
        }
    }

    /// Index of the records written so far.
    pub fn fai(&self) -> &Fai {
        &self.fai
    }

    /// Appends a record `>name\n` and the wrapped sequence. Rejects empty
    /// names, names with ASCII whitespace or control bytes, names already
    /// written, empty sequences, sequences longer than a `.fai` can
    /// describe, and sequence bytes other than printable non-space ASCII
    /// or that are `>`.
    pub fn add(&mut self, name: &str, seq: &[u8]) -> Result<()> {
        self.check_intact()?;
        let entry = self.entry(name, seq)?;
        if let Err(e) = self.write_record(name, seq, entry.line_bases as usize) {
            self.broken = true;
            return Err(e);
        }
        self.pos = entry.offset + text_len(seq.len() as u64, entry.line_bases as u64);
        self.fai.insert_missing(entry);
        Ok(())
    }

    fn write_record(&mut self, name: &str, seq: &[u8], line_bases: usize) -> Result<()> {
        self.put(b">")?;
        self.put(name.as_bytes())?;
        self.put(b"\n")?;
        for line in seq.chunks(line_bases) {
            self.put(line)?;
            self.put(b"\n")?;
        }
        Ok(())
    }

    /// After a failed write the sink holds part of a record, which no
    /// index describes; the writer then refuses to go on.
    fn check_intact(&self) -> Result<()> {
        if self.broken {
            return Err(invalid("FASTA writer", "an earlier write failed"));
        }
        Ok(())
    }

    /// The index entry of a record appended now, once the record is valid.
    fn entry(&self, name: &str, seq: &[u8]) -> Result<FaiEntry> {
        if name.is_empty() {
            return Err(invalid("FASTA record", "empty name"));
        }
        if let Some(c) = name.bytes().find(|&c| c <= b' ' || c == 0x7F) {
            return Err(invalid(
                "FASTA record",
                format!("name {name:?} contains byte {c:#04x}"),
            ));
        }
        if self.fai.get(name).is_some() {
            return Err(invalid("FASTA record", format!("duplicate name {name:?}")));
        }
        if seq.is_empty() {
            return Err(invalid("FASTA record", format!("{name}: empty sequence")));
        }
        if let Some(at) = seq
            .iter()
            .position(|&c| !(0x21..=0x7E).contains(&c) || c == b'>')
        {
            return Err(invalid(
                "FASTA record",
                format!(
                    "{name}: byte {:#04x} at {at} is not a printable non-space base",
                    seq[at]
                ),
            ));
        }
        let too_long = || invalid("FASTA record", format!("{name}: too long for a .fai entry"));
        let len = u32::try_from(seq.len()).map_err(|_| too_long())?;
        let line_bases = match self.line_width {
            0 => len,
            w => len.min(u32::try_from(w).unwrap_or(u32::MAX)),
        };
        let line_width = line_bases.checked_add(1).ok_or_else(too_long)?;
        Ok(FaiEntry {
            name: name.to_owned(),
            len,
            offset: self.pos + name.len() as u64 + 2,
            line_bases,
            line_width,
        })
    }

    fn put(&mut self, bytes: &[u8]) -> Result<()> {
        if self.buf.len() + bytes.len() > BUFFER_LEN {
            self.drain()?;
            if bytes.len() >= BUFFER_LEN {
                return self.sink.write_all(bytes);
            }
        }
        self.buf.extend_from_slice(bytes);
        Ok(())
    }

    fn drain(&mut self) -> Result<()> {
        if !self.buf.is_empty() {
            self.sink.write_all(&self.buf)?;
            self.buf.clear();
        }
        Ok(())
    }

    /// Ends the file and returns the sink, the `.fai` and, for bgzip, the
    /// `.gzi`. A file without records is an error: it has no index
    /// (`samtools faidx` rejects it), and nothing is written for it.
    pub fn finish(mut self) -> Result<(W, Fai, Option<Gzi>)> {
        self.check_intact()?;
        if self.fai.is_empty() {
            return Err(invalid(
                "FASTA writer",
                "no records written; an empty FASTA has no index",
            ));
        }
        self.drain()?;
        match self.sink {
            Sink::Plain(mut out) => {
                out.flush()?;
                Ok((out, self.fai, None))
            }
            Sink::Bgzf(out) => {
                let (out, gzi) = out.finish()?;
                Ok((out, self.fai, Some(gzi)))
            }
        }
    }
}

/// Text bytes of `len` bases in lines of `line_bases`, terminators included.
fn text_len(len: u64, line_bases: u64) -> u64 {
    len + len.div_ceil(line_bases)
}
