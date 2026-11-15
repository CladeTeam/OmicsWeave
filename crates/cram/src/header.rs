//! The file definition and SAM header of a CRAM file (`cram_dopen`,
//! `cram_read_SAM_hdr`), and container reads through [`ReadAt`].

use std::io::Cursor;

use omicsweave_core::{Error, ReadAt, Result};
use omicsweave_sam::SamHeader;

use crate::container::{Block, ContainerHeader, FILE_DEF_LEN, file_version};
use crate::varint::int32;

/// A reference sequence of the SAM header (`@SQ`).
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Reference {
    pub name: String,
    pub len: u32,
}

/// The header of a CRAM file: its version and the SAM header it holds.
#[derive(Debug, Clone)]
pub struct CramHeader {
    major: u8,
    minor: u8,
    sam: SamHeader,
    refs: Vec<Reference>,
    pub(crate) ref_lens: Vec<i64>,
    /// `@RG` IDs in header order.
    pub(crate) read_groups: Vec<Vec<u8>>,
    /// The index of a last read group named `UNKNOWN`, which htslib treats
    /// as none; -1 otherwise.
    pub(crate) unknown_rg: i32,
    first_container: u64,
}

/// A container read from a source: its header and the bytes after it.
pub(crate) struct ReadContainer {
    pub header: ContainerHeader,
    /// File offset of the container.
    pub offset: u64,
}

/// Reads the container header at `offset`; `None` at the end of the
/// file.
pub(crate) fn read_container_header<S: ReadAt + ?Sized>(
    source: &S,
    offset: u64,
    buf: &mut Vec<u8>,
) -> Result<Option<ReadContainer>> {
    let len = source.len();
    if offset >= len {
        return Ok(None);
    }
    let mut want = 1024u64;
    loop {
        let n = want.min(len - offset) as usize;
        buf.clear();
        buf.resize(n, 0);
        source.read_exact_at(buf, offset)?;
        match ContainerHeader::parse(buf) {
            Ok(Some(header)) => return Ok(Some(ReadContainer { header, offset })),
            Ok(None) => return Ok(None),
            Err(Error::Truncated { .. }) if (n as u64) < len - offset => want *= 8,
            Err(e) => return Err(e),
        }
    }
}

/// Reads the `length` bytes of a container after its header into `buf`.
pub(crate) fn read_container_body<S: ReadAt + ?Sized>(
    source: &S,
    c: &ReadContainer,
    buf: &mut Vec<u8>,
) -> Result<()> {
    let at = c.offset + c.header.header_len as u64;
    let length = usize::try_from(c.header.length.max(0)).unwrap_or(0);
    if at + length as u64 > source.len() {
        return Err(Error::Truncated {
            what: "CRAM container",
            need: length as u64,
            have: source.len().saturating_sub(at),
        });
    }
    buf.clear();
    buf.resize(length, 0);
    source.read_exact_at(buf, at)
}

impl CramHeader {
    /// Reads the file definition and the header container of `source`.
    /// CRAM versions other than 3.0 and 3.1 are [`Error::Unsupported`].
    pub fn read<S: ReadAt + ?Sized>(source: &S) -> Result<CramHeader> {
        let mut def = [0u8; FILE_DEF_LEN];
        if source.len() < FILE_DEF_LEN as u64 {
            return Err(Error::Truncated {
                what: "CRAM file definition",
                need: FILE_DEF_LEN as u64,
                have: source.len(),
            });
        }
        source.read_exact_at(&mut def, 0)?;
        let (major, minor) = file_version(&def)?;
        let mut buf = Vec::new();
        let c = read_container_header(source, FILE_DEF_LEN as u64, &mut buf)?
            .ok_or_else(|| Error::invalid("CRAM header", "no header container"))?;
        if c.header.num_blocks < 1 {
            return Err(Error::invalid(
                "CRAM header",
                "header container has no blocks",
            ));
        }
        let first_container =
            FILE_DEF_LEN as u64 + c.header.header_len as u64 + c.header.length.max(0) as u64;
        read_container_body(source, &c, &mut buf)?;
        let mut at = 0;
        let block = Block::parse(&buf, &mut at)?;
        let mut scratch = Vec::new();
        let data = block.uncompress(&mut scratch, true)?;
        let mut cp = data;
        let len = int32(&mut cp, "CRAM header")?;
        if len < 0 || len as usize > cp.len() {
            return Err(Error::invalid("CRAM header", format!("text length {len}")));
        }
        let text = &cp[..len as usize];
        Self::from_text(major, minor, text, first_container)
    }

    fn from_text(major: u8, minor: u8, text: &[u8], first_container: u64) -> Result<CramHeader> {
        let end = memchr::memchr(0, text).unwrap_or(text.len());
        let sam = SamHeader::read(&mut Cursor::new(&text[..end]))?;
        let mut refs = Vec::with_capacity(sam.references().len());
        for r in sam.references() {
            let name = String::from_utf8(r.name.clone())
                .map_err(|_| Error::invalid("CRAM header", "a reference name is not UTF-8"))?;
            refs.push(Reference {
                name,
                len: r.len.min(u32::MAX as u64) as u32,
            });
        }
        let ref_lens = sam.references().iter().map(|r| r.len as i64).collect();
        let mut read_groups = Vec::new();
        for line in sam.text().split(|&c| c == b'\n') {
            if let Some(rest) = line.strip_prefix(b"@RG\t")
                && let Some(id) = rest
                    .split(|&c| c == b'\t')
                    .find_map(|t| t.strip_prefix(b"ID:"))
            {
                read_groups.push(id.to_vec());
            }
        }
        let unknown_rg = match read_groups.last() {
            Some(id) if id.as_slice() == b"UNKNOWN" => read_groups.len() as i32 - 1,
            _ => -1,
        };
        Ok(CramHeader {
            major,
            minor,
            sam,
            refs,
            ref_lens,
            read_groups,
            unknown_rg,
            first_container,
        })
    }

    /// `(major, minor)`: 3.0 or 3.1.
    pub fn version(&self) -> (u8, u8) {
        (self.major, self.minor)
    }

    /// The SAM header text as htslib holds it after parsing (what it
    /// writes into a BAM converted from this file).
    pub fn text(&self) -> &[u8] {
        self.sam.text()
    }

    /// The parsed SAM header.
    pub fn sam(&self) -> &SamHeader {
        &self.sam
    }

    /// References in reference-id order.
    pub fn references(&self) -> &[Reference] {
        &self.refs
    }

    /// The reference id of `name` (or of an `AN` alternative name).
    pub fn tid(&self, name: &str) -> Option<usize> {
        self.sam.tid(name.as_bytes()).map(|t| t as usize)
    }

    /// The length of `chrom`; [`Error::UnknownChrom`] if absent.
    pub fn chrom_len(&self, chrom: &str) -> Result<u32> {
        self.tid(chrom)
            .map(|t| self.refs[t].len)
            .ok_or_else(|| Error::UnknownChrom(chrom.to_owned()))
    }

    /// File offset of the first data container.
    pub fn first_container(&self) -> u64 {
        self.first_container
    }
}
