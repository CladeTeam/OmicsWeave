//! Byte layout of the fixed header, checksummed sections and the block
//! directory.

use xxhash_rust::xxh3::xxh3_64;

use crate::error::{Error, Result, invalid};

pub const MAGIC: [u8; 8] = *b"\x89FFM\r\n\x1a\n";
pub const MAJOR: u16 = 1;
pub const MINOR: u16 = 0;
pub const HEADER_LEN: usize = 128;
pub const DIR_ENTRY_LEN: usize = 16;
pub const CHECKSUM_LEN: usize = 8;

/// A little-endian reader over a byte slice that reports truncation against
/// the structure it names.
pub(crate) struct Cursor<'a> {
    bytes: &'a [u8],
    pos: usize,
    what: &'static str,
}

impl<'a> Cursor<'a> {
    pub fn new(bytes: &'a [u8], what: &'static str) -> Self {
        Cursor {
            bytes,
            pos: 0,
            what,
        }
    }

    pub fn bytes(&mut self, n: usize) -> Result<&'a [u8]> {
        let end = self.pos.checked_add(n).filter(|&e| e <= self.bytes.len());
        let Some(end) = end else {
            return Err(Error::Truncated {
                what: self.what,
                need: (self.pos as u64).saturating_add(n as u64),
                have: self.bytes.len() as u64,
            });
        };
        let s = &self.bytes[self.pos..end];
        self.pos = end;
        Ok(s)
    }

    pub fn u8(&mut self) -> Result<u8> {
        Ok(self.bytes(1)?[0])
    }

    pub fn u16(&mut self) -> Result<u16> {
        Ok(u16::from_le_bytes(self.bytes(2)?.try_into().expect("2")))
    }

    pub fn u32(&mut self) -> Result<u32> {
        Ok(u32::from_le_bytes(self.bytes(4)?.try_into().expect("4")))
    }

    pub fn u64(&mut self) -> Result<u64> {
        Ok(u64::from_le_bytes(self.bytes(8)?.try_into().expect("8")))
    }

    pub fn zeros(&mut self, n: usize) -> Result<()> {
        if self.bytes(n)?.iter().any(|&b| b != 0) {
            return Err(invalid(self.what, "reserved bytes are not zero"));
        }
        Ok(())
    }

    /// Fails unless every byte was consumed.
    pub fn finish(self) -> Result<()> {
        if self.pos != self.bytes.len() {
            return Err(invalid(
                self.what,
                format!("{} trailing bytes", self.bytes.len() - self.pos),
            ));
        }
        Ok(())
    }
}

/// Appends the checksum of everything in `section`.
pub(crate) fn seal(section: &mut Vec<u8>) {
    let sum = xxh3_64(section);
    section.extend_from_slice(&sum.to_le_bytes());
}

/// The payload of a sealed section, after checking its trailing checksum.
pub(crate) fn unseal<'a>(section: &'a [u8], what: &'static str) -> Result<&'a [u8]> {
    if section.len() < CHECKSUM_LEN {
        return Err(Error::Truncated {
            what,
            need: CHECKSUM_LEN as u64,
            have: section.len() as u64,
        });
    }
    let (payload, sum) = section.split_at(section.len() - CHECKSUM_LEN);
    if xxh3_64(payload) != u64::from_le_bytes(sum.try_into().expect("8")) {
        return Err(invalid(what, "checksum mismatch"));
    }
    Ok(payload)
}

/// A section's place in the file.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub struct Span {
    pub offset: u64,
    pub len: u64,
}

impl Span {
    fn end(self) -> Option<u64> {
        self.offset.checked_add(self.len)
    }
}

#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub struct Header {
    pub minor: u16,
    pub motif_count: u64,
    pub block_count: u64,
    pub alphabet: Span,
    pub schema: Span,
    pub ids: Span,
    pub directory: Span,
    pub index: Span,
}

impl Header {
    pub fn encode(&self) -> [u8; HEADER_LEN] {
        let mut b = Vec::with_capacity(HEADER_LEN);
        b.extend_from_slice(&MAGIC);
        b.extend_from_slice(&MAJOR.to_le_bytes());
        b.extend_from_slice(&self.minor.to_le_bytes());
        b.extend_from_slice(&0u32.to_le_bytes());
        b.extend_from_slice(&self.motif_count.to_le_bytes());
        b.extend_from_slice(&self.block_count.to_le_bytes());
        for span in [
            self.alphabet,
            self.schema,
            self.ids,
            self.directory,
            self.index,
        ] {
            b.extend_from_slice(&span.offset.to_le_bytes());
            b.extend_from_slice(&span.len.to_le_bytes());
        }
        b.resize(HEADER_LEN - CHECKSUM_LEN, 0);
        seal(&mut b);
        b.try_into().expect("header length")
    }

    /// Parses and checks the header of a file of `file_len` bytes.
    pub fn parse(bytes: &[u8], file_len: u64) -> Result<Header> {
        const WHAT: &str = "FFM header";
        if bytes.len() < HEADER_LEN || bytes[..8] != MAGIC {
            return Err(invalid(WHAT, "not an FFM file (bad magic)"));
        }
        let payload = unseal(&bytes[..HEADER_LEN], WHAT)?;
        let mut r = Cursor::new(payload, WHAT);
        r.bytes(8)?;
        if r.u16()? != MAJOR {
            return Err(Error::Unsupported("FFM major version other than 1"));
        }
        let minor = r.u16()?;
        if r.u32()? != 0 {
            return Err(invalid(WHAT, "unknown flags"));
        }
        let motif_count = r.u64()?;
        let block_count = r.u64()?;
        let mut spans = [Span::default(); 5];
        for span in &mut spans {
            span.offset = r.u64()?;
            span.len = r.u64()?;
        }
        r.zeros(8)?;
        r.finish()?;
        let [alphabet, schema, ids, directory, index] = spans;
        if motif_count > crate::MAX_MOTIFS as u64 {
            return Err(invalid(
                WHAT,
                "motif count exceeds the four-column index limit",
            ));
        }
        // Metadata sections follow the header back to back; blocks fill the
        // rest of the file.
        let ends = [alphabet, schema, ids, directory, index].map(Span::end);
        let ordered = alphabet.offset == HEADER_LEN as u64
            && ends[0] == Some(schema.offset)
            && ends[1] == Some(ids.offset)
            && ends[2] == Some(directory.offset)
            && ends[3] == Some(index.offset)
            && ends[4].is_some_and(|e| e <= file_len);
        if !ordered {
            return Err(invalid(WHAT, "sections are not laid out in order"));
        }
        let dir_len = block_count
            .checked_mul(DIR_ENTRY_LEN as u64)
            .and_then(|n| n.checked_add(CHECKSUM_LEN as u64));
        if dir_len != Some(directory.len) {
            return Err(invalid(WHAT, "directory length does not match block count"));
        }
        if (block_count == 0) != (motif_count == 0) || block_count > motif_count {
            return Err(invalid(WHAT, "block count does not fit the motif count"));
        }
        Ok(Header {
            minor,
            motif_count,
            block_count,
            alphabet,
            schema,
            ids,
            directory,
            index,
        })
    }

    /// The metadata sections, which follow the header back to back.
    pub fn metadata(&self) -> Span {
        Span {
            offset: HEADER_LEN as u64,
            len: self.index.offset + self.index.len - HEADER_LEN as u64,
        }
    }

    /// The bytes holding blocks: from the end of the metadata to `file_len`.
    pub fn blocks(&self, file_len: u64) -> Span {
        let offset = self.index.offset + self.index.len;
        Span {
            offset,
            len: file_len - offset,
        }
    }
}

/// One block of the directory.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct DirEntry {
    pub offset: u64,
    pub length: u32,
    pub record_count: u32,
}

impl DirEntry {
    pub fn write_to(&self, out: &mut Vec<u8>) {
        out.extend_from_slice(&self.offset.to_le_bytes());
        out.extend_from_slice(&self.length.to_le_bytes());
        out.extend_from_slice(&self.record_count.to_le_bytes());
    }
}

/// Parses the directory and checks that its blocks follow each other through
/// the block region and hold every motif.
pub(crate) fn parse_directory(
    section: &[u8],
    header: &Header,
    file_len: u64,
) -> Result<Vec<DirEntry>> {
    const WHAT: &str = "FFM directory";
    let payload = unseal(section, WHAT)?;
    let mut r = Cursor::new(payload, WHAT);
    let region = header.blocks(file_len);
    let mut entries = Vec::with_capacity(header.block_count as usize);
    let mut next_byte = region.offset;
    let mut motifs = 0u64;
    for i in 0..header.block_count {
        let e = DirEntry {
            offset: r.u64()?,
            length: r.u32()?,
            record_count: r.u32()?,
        };
        if e.offset != next_byte || e.record_count == 0 || e.length < 8 + 3 * 32 {
            return Err(invalid(
                WHAT,
                format!("block {i} does not follow the previous block"),
            ));
        }
        next_byte = e.offset.saturating_add(e.length as u64);
        motifs += e.record_count as u64;
        entries.push(e);
    }
    r.finish()?;
    if next_byte != region.offset + region.len || motifs != header.motif_count {
        return Err(invalid(
            WHAT,
            "blocks do not fill the block region or hold every motif",
        ));
    }
    Ok(entries)
}
