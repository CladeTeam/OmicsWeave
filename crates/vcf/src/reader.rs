//! [`VariantReader`]: VCF or BCF, plain, gzip or BGZF, read front to back.

use std::io::{BufRead, BufReader, Cursor, Read};

use omicsweave_core::{Error, Result};

use crate::bcf::{Fixed, parse_record};
use crate::header::Header;
use crate::record::Record;
use crate::text::{Scratch, parse_line};

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Format {
    Vcf,
    Bcf,
}

/// Reads a whole variant file in order. Compression (gzip of any number of
/// members, which includes BGZF) and the format are detected from the
/// first bytes.
pub struct VariantReader<'r> {
    input: Box<dyn BufRead + 'r>,
    header: Header,
    format: Format,
    line: Vec<u8>,
    line_no: usize,
    shared: Vec<u8>,
    indiv: Vec<u8>,
    scratch: Scratch,
    failed: bool,
}

const CHUNK: usize = 1 << 18;

fn io(e: std::io::Error) -> Error {
    Error::from(e)
}

fn prefix(input: &mut dyn Read, bytes: &mut [u8]) -> Result<usize> {
    let mut at = 0;
    while at < bytes.len() {
        match input.read(&mut bytes[at..]) {
            Ok(0) => break,
            Ok(n) => at += n,
            Err(e) if e.kind() == std::io::ErrorKind::Interrupted => continue,
            Err(e) => return Err(io(e)),
        }
    }
    Ok(at)
}

impl<'r> VariantReader<'r> {
    /// Opens `source` and reads the header. Headers and textual dictionary
    /// identifiers must be UTF-8; alleles and annotation values remain bytes.
    pub fn new(mut source: impl Read + 'r) -> Result<VariantReader<'r>> {
        let mut first = [0u8; 2];
        let n = prefix(&mut source, &mut first)?;
        let compressed = first[..n] == [0x1f, 0x8b];
        let replay = Cursor::new(first[..n].to_vec()).chain(source);
        let mut raw: Box<dyn Read + 'r> = if compressed {
            Box::new(flate2::bufread::MultiGzDecoder::new(
                BufReader::with_capacity(CHUNK, replay),
            ))
        } else {
            Box::new(replay)
        };
        let mut magic = [0u8; 5];
        let n = prefix(raw.as_mut(), &mut magic)?;
        let raw = Box::new(BufReader::with_capacity(
            CHUNK,
            Cursor::new(magic[..n].to_vec()).chain(raw),
        ));
        let mut reader = VariantReader {
            input: raw,
            header: Header::default(),
            format: Format::Vcf,
            line: Vec::new(),
            line_no: 0,
            shared: Vec::new(),
            indiv: Vec::new(),
            scratch: Scratch::default(),
            failed: false,
        };
        if magic[..n].starts_with(b"BCF") {
            if n != 5 || magic[..4] != *b"BCF\x02" || !matches!(magic[4], 1 | 2) {
                return Err(Error::invalid(
                    "bcf header",
                    "unsupported or truncated magic",
                ));
            }
            reader.format = Format::Bcf;
            reader.read_bcf_header()?;
        } else {
            reader.read_vcf_header()?;
        }
        Ok(reader)
    }

    pub fn header(&self) -> &Header {
        &self.header
    }

    pub fn format(&self) -> Format {
        self.format
    }

    fn read_line(&mut self) -> Result<bool> {
        self.line.clear();
        let n = self.input.read_until(b'\n', &mut self.line).map_err(io)?;
        if n == 0 {
            return Ok(false);
        }
        self.line_no += 1;
        while matches!(self.line.last(), Some(b'\n' | b'\r')) {
            self.line.pop();
        }
        Ok(true)
    }

    fn read_vcf_header(&mut self) -> Result<()> {
        let mut text = String::new();
        loop {
            if !self.read_line()? {
                return Err(Error::invalid("vcf header", "no #CHROM line"));
            }
            if !self.line.starts_with(b"#") {
                return Err(Error::Format {
                    line: self.line_no as u64,
                    detail: "header line without '#'".into(),
                });
            }
            let line = std::str::from_utf8(&self.line)
                .map_err(|_| Error::invalid("vcf header", "header is not UTF-8"))?;
            text.push_str(line);
            text.push('\n');
            if !self.line.starts_with(b"##") {
                break;
            }
        }
        self.header = Header::parse(&text)?;
        Ok(())
    }

    fn read_exact(&mut self, buf: &mut [u8]) -> Result<bool> {
        let mut at = 0;
        while at < buf.len() {
            let n = match self.input.read(&mut buf[at..]) {
                Err(e) if e.kind() == std::io::ErrorKind::Interrupted => continue,
                r => r.map_err(io)?,
            };
            if n == 0 {
                if at == 0 {
                    return Ok(false);
                }
                return Err(Error::Truncated {
                    what: "bcf",
                    need: buf.len() as u64,
                    have: at as u64,
                });
            }
            at += n;
        }
        Ok(true)
    }

    /// Grow only after observing input bytes, rather than trusting a stored
    /// length to justify an allocation before truncation can be detected.
    fn read_block(&mut self, buf: &mut Vec<u8>, len: usize) -> Result<()> {
        buf.clear();
        while buf.len() < len {
            let bytes = match self.input.fill_buf() {
                Err(e) if e.kind() == std::io::ErrorKind::Interrupted => continue,
                r => r.map_err(io)?,
            };
            if bytes.is_empty() {
                return Err(Error::Truncated {
                    what: "bcf block",
                    need: len as u64,
                    have: buf.len() as u64,
                });
            }
            let n = bytes.len().min(len - buf.len());
            buf.try_reserve(n)
                .map_err(|_| Error::invalid("bcf block", "allocation failed"))?;
            buf.extend_from_slice(&bytes[..n]);
            self.input.consume(n);
        }
        Ok(())
    }

    fn read_bcf_header(&mut self) -> Result<()> {
        let mut head = [0u8; 9];
        if !self.read_exact(&mut head)? {
            return Err(Error::invalid("bcf header", "empty"));
        }
        let len = u32::from_le_bytes(head[5..9].try_into().unwrap()) as usize;
        let mut text = Vec::new();
        self.read_block(&mut text, len)?;
        if text.last() != Some(&0) {
            return Err(Error::invalid("bcf header", "missing terminating NUL"));
        }
        let text = std::str::from_utf8(&text)
            .map_err(|_| Error::invalid("bcf header", "header is not UTF-8"))?;
        self.header = Header::parse(text)?;
        Ok(())
    }

    /// Reads the next record into `rec`; `false` at the end. A malformed
    /// VCF line is an error after which reading continues with the next
    /// line; a malformed BCF record or an I/O error ends the stream. `rec`
    /// is cleared on an error or at the end, retaining its buffer capacities.
    pub fn read(&mut self, rec: &mut Record) -> Result<bool> {
        rec.clear();
        if self.failed {
            return Ok(false);
        }
        let result = self.read_next(rec);
        if let Err(e) = &result {
            rec.clear();
            self.failed = self.format == Format::Bcf || !matches!(e, Error::Format { .. });
        }
        result
    }

    fn read_next(&mut self, rec: &mut Record) -> Result<bool> {
        match self.format {
            Format::Vcf => loop {
                if !self.read_line()? {
                    return Ok(false);
                }
                if self.line.is_empty() {
                    continue;
                }
                return parse_line(&mut self.header, &self.line, rec, &mut self.scratch)
                    .map(|()| true)
                    .map_err(|e| Error::Format {
                        line: self.line_no as u64,
                        detail: e.to_string(),
                    });
            },
            Format::Bcf => {
                let mut head = [0u8; 32];
                if !self.read_exact(&mut head)? {
                    return Ok(false);
                }
                let l_shared = u32::from_le_bytes(head[0..4].try_into().unwrap()) as usize;
                let l_indiv = u32::from_le_bytes(head[4..8].try_into().unwrap()) as usize;
                if l_shared < 24 {
                    return Err(Error::invalid("bcf record", "shared block too short"));
                }
                let fixed = Fixed::parse(head[8..32].try_into().unwrap());
                let mut shared = std::mem::take(&mut self.shared);
                let mut indiv = std::mem::take(&mut self.indiv);
                let result = (|| {
                    self.read_block(&mut shared, l_shared - 24)?;
                    self.read_block(&mut indiv, l_indiv)?;
                    parse_record(&self.header, &fixed, &shared, &indiv, rec)
                })();
                self.shared = shared;
                self.indiv = indiv;
                result.map(|()| true)
            }
        }
    }
}
