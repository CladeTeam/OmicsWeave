//! The BAM header as htslib's `bam_hdr_read` reads it: the SAM header text
//! and the binary reference list, which is what records' reference ids
//! index.

use std::collections::HashMap;

use omicsweave_core::{Error, Result};

/// The BAM magic, `BAM\1`.
pub const BAM_MAGIC: [u8; 4] = *b"BAM\x01";

/// Decompressed bytes taken in order, as a header is read.
pub(crate) trait Take {
    /// The next `n` bytes; [`Error::Truncated`] when the text ends first.
    fn take(&mut self, n: usize, what: &'static str) -> Result<&[u8]>;
}

/// A reference sequence of the binary header.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Reference {
    /// The name up to its first NUL.
    pub name: String,
    pub len: u32,
}

/// A BAM header. Reference names must be UTF-8; when a name repeats,
/// lookups by name find its first reference.
#[derive(Debug, Clone, Default)]
pub struct BamHeader {
    text: Vec<u8>,
    refs: Vec<Reference>,
    names: HashMap<String, usize>,
}

fn i32_of(b: &[u8]) -> i32 {
    i32::from_le_bytes(b[..4].try_into().unwrap())
}

fn u32_of(b: &[u8]) -> u32 {
    u32::from_le_bytes(b[..4].try_into().unwrap())
}

impl BamHeader {
    pub(crate) fn read(input: &mut impl Take) -> Result<BamHeader> {
        let magic = input.take(4, "BAM header")?;
        if magic != BAM_MAGIC {
            return Err(Error::BadMagic {
                what: "BAM header",
                found: u32_of(magic),
            });
        }
        let l_text = u32_of(input.take(4, "BAM header")?) as usize;
        let text = input.take(l_text, "BAM header text")?.to_vec();
        let n_ref = i32_of(input.take(4, "BAM header")?);
        if n_ref < 0 {
            return Err(Error::invalid(
                "BAM header",
                format!("negative reference count {n_ref}"),
            ));
        }
        let mut h = BamHeader {
            text,
            refs: Vec::with_capacity((n_ref as usize).min(1 << 16)),
            names: HashMap::new(),
        };
        for tid in 0..n_ref as usize {
            let l_name = i32_of(input.take(4, "BAM reference")?);
            if l_name <= 0 {
                return Err(Error::invalid(
                    "BAM reference",
                    format!("reference {tid} has name length {l_name}"),
                ));
            }
            let raw = input.take(l_name as usize, "BAM reference name")?;
            let raw = &raw[..raw.iter().position(|&c| c == 0).unwrap_or(raw.len())];
            let name = String::from_utf8(raw.to_vec()).map_err(|_| {
                Error::invalid(
                    "BAM reference",
                    format!("name of reference {tid} is not UTF-8"),
                )
            })?;
            let len = u32_of(input.take(4, "BAM reference")?);
            h.names.entry(name.clone()).or_insert(tid);
            h.refs.push(Reference { name, len });
        }
        Ok(h)
    }

    /// The SAM header text as stored (`l_text` bytes, NUL padding
    /// included).
    pub fn text(&self) -> &[u8] {
        &self.text
    }

    /// References in reference-id order.
    pub fn references(&self) -> &[Reference] {
        &self.refs
    }

    /// The reference id of `name`.
    pub fn tid(&self, name: &str) -> Option<usize> {
        self.names.get(name).copied()
    }

    /// The length of `chrom`; [`Error::UnknownChrom`] if absent.
    pub fn chrom_len(&self, chrom: &str) -> Result<u32> {
        self.tid(chrom)
            .map(|t| self.refs[t].len)
            .ok_or_else(|| Error::UnknownChrom(chrom.to_owned()))
    }
}
