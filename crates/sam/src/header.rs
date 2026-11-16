//! The SAM header as htslib reads it (`sam_hdr_build_from_sam_file` and
//! `sam_hrecs_update_hashes`, `header.c`) and writes it to BAM
//! (`bam_hdr_write`, `sam.c`).

use std::collections::HashMap;
use std::hash::{BuildHasherDefault, Hasher};
use std::io::BufRead;
use std::sync::Arc;

use omicsweave_core::{Error, Result};

use crate::cnum::strtoll;

/// A reference sequence of the header: an `@SQ` line's `SN` and `LN`.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Reference {
    /// The last `SN` value of the line.
    pub name: Vec<u8>,
    /// `LN` as `strtoll` reads it; always positive. BAM stores it as
    /// `u32`, saturated at `u32::MAX` as htslib does.
    pub len: u64,
}

/// FxHash: reference names are short and trusted, and the lookup runs once
/// or twice per record.
#[derive(Default, Clone, Copy)]
pub(crate) struct FxHasher(u64);

impl FxHasher {
    #[inline]
    fn add(&mut self, w: u64) {
        self.0 = (self.0.rotate_left(5) ^ w).wrapping_mul(0x51_7c_c1_b7_27_22_0a_95);
    }
}

impl Hasher for FxHasher {
    #[inline]
    fn write(&mut self, bytes: &[u8]) {
        let mut chunks = bytes.chunks_exact(8);
        for c in &mut chunks {
            self.add(u64::from_le_bytes(c.try_into().unwrap()));
        }
        let rest = chunks.remainder();
        if !rest.is_empty() {
            let mut w = [0u8; 8];
            w[..rest.len()].copy_from_slice(rest);
            self.add(u64::from_le_bytes(w));
        }
    }

    #[inline]
    fn write_usize(&mut self, n: usize) {
        self.add(n as u64);
    }

    #[inline]
    fn finish(&self) -> u64 {
        self.0
    }
}

/// Reference name (`SN`, or an `AN` alternative name) to tid.
pub(crate) type NameMap = HashMap<Box<[u8]>, i32, BuildHasherDefault<FxHasher>>;

/// A SAM header: every leading line that starts with `@`.
///
/// [`SamHeader::read`] applies htslib's checks, so a header htslib rejects
/// is an error here too: a line without a two-letter type, a type other
/// than `HD`/`SQ`/`RG`/`PG`/`CO`, a field that is not `XX:value`, an
/// `@SQ` without `SN` or with a missing, non-positive or conflicting `LN`,
/// a duplicate `SN`, an `@RG` or `@PG` without `ID`. Header lines holding a
/// NUL byte are rejected, where htslib would silently cut them short.
#[derive(Debug, Clone, Default)]
pub struct SamHeader {
    raw: Vec<u8>,
    text: Vec<u8>,
    lines: u64,
    refs: Vec<Reference>,
    names: Arc<NameMap>,
}

impl SamHeader {
    /// Reads the header from decompressed SAM text, consuming every line
    /// that starts with `@` and nothing after them: the first data line is
    /// left in `input` for the caller.
    pub fn read<R: BufRead + ?Sized>(input: &mut R) -> Result<SamHeader> {
        let mut h = SamHeader::default();
        // Build under one exclusive borrow; readers and encoders share the
        // immutable dictionary after the header is complete.
        let names = Arc::get_mut(&mut h.names).expect("new header dictionary is exclusive");
        let mut first_hd = None;
        let mut line = Vec::new();
        loop {
            let buf = input.fill_buf()?;
            if buf.first() != Some(&b'@') {
                break;
            }
            line.clear();
            input.read_until(b'\n', &mut line)?;
            h.raw.extend_from_slice(&line);
            h.lines += 1;
            // hts_getline: a "\n" ends the line, and a "\r" before it goes too.
            if line.last() == Some(&b'\n') {
                line.pop();
                if line.last() == Some(&b'\r') {
                    line.pop();
                }
            }
            let lno = h.lines;
            let fail = |detail: &str| Error::Format {
                line: lno,
                detail: detail.to_string(),
            };
            if line.len() < 3 || !line[1].is_ascii_alphabetic() || !line[2].is_ascii_alphabetic() {
                return Err(fail("header line does not have a two character key"));
            }
            if memchr::memchr(0, &line).is_some() {
                return Err(fail("NUL byte in header line"));
            }
            // htslib accepts a bare "@CO" and gives it the missing tab.
            if line == b"@CO" {
                line.push(b'\t');
            }
            let ty = [line[1], line[2]];
            let tab = line.get(3) == Some(&b'\t');
            let valid = match &ty {
                b"HD" | b"SQ" | b"RG" | b"PG" => tab,
                _ => ty[0] == b'C' && ty[1] == b'O',
            };
            if !valid {
                return Err(fail("header line type must be HD, SQ, RG, PG or CO"));
            }
            if &ty == b"CO" {
                if !tab {
                    return Err(fail("missing tab after @CO"));
                }
            } else {
                let tags = tags(&line[3..]).ok_or_else(|| fail("malformed key:value pair"))?;
                match &ty {
                    b"SQ" => Self::add_reference(&mut h.refs, names, &tags).map_err(fail)?,
                    b"RG" | b"PG" if !tags.iter().any(|t| t.starts_with(b"ID")) => {
                        return Err(fail("header line has no ID tag"));
                    }
                    _ => {}
                }
            }
            let start = h.text.len();
            h.text.extend_from_slice(&line);
            h.text.push(b'\n');
            if &ty == b"HD" && first_hd.is_none() {
                first_hd = Some(start..h.text.len());
            }
        }
        // htslib rebuilds the text from its parsed lines, with the first
        // @HD moved to the top.
        // Keep one reusable input line and one contiguous text buffer;
        // moving a late @HD needs no separately retained line buffers.
        if let Some(hd) = first_hd {
            h.text[..hd.end].rotate_right(hd.len());
        }
        if h.text.len() > u32::MAX as usize {
            return Err(Error::invalid("SAM header", "header text too long for BAM"));
        }
        Ok(h)
    }

    /// `@SQ` bookkeeping as `sam_hrecs_update_hashes` does it, including
    /// `AN` alternative names, which resolve RNAME and RNEXT as well.
    fn add_reference(
        refs: &mut Vec<Reference>,
        names: &mut NameMap,
        tags: &[&[u8]],
    ) -> std::result::Result<(), &'static str> {
        let mut name: Option<&[u8]> = None;
        let mut alt: Option<&[u8]> = None;
        let mut len: i64 = -1;
        let mut conflicting = false;
        for t in tags {
            match &t[..2] {
                b"SN" => name = Some(&t[3..]),
                b"LN" => {
                    let before = len;
                    len = strtoll(&t[3..]);
                    conflicting |= before != -1 && before != len;
                }
                b"AN" => alt = Some(&t[3..]),
                _ => {}
            }
        }
        let name = name.ok_or("@SQ line with no SN tag")?;
        if len <= 0 {
            return Err("@SQ line with invalid or no LN tag");
        }
        if conflicting {
            return Err("@SQ line with multiple LN tags of different values");
        }
        if let Some(&tid) = names.get(name) {
            // An SN equal to an earlier SN is a duplicate; one equal to an
            // earlier AN takes the name over.
            if refs[tid as usize].name == name {
                return Err("duplicate @SQ SN");
            }
        }
        let tid = refs.len() as i32;
        refs.push(Reference {
            name: name.to_vec(),
            len: len as u64,
        });
        names.insert(name.into(), tid);
        for a in alt.into_iter().flat_map(|a| a.split(|&c| c == b',')) {
            if !a.is_empty() {
                names.entry(a.into()).or_insert(tid);
            }
        }
        Ok(())
    }

    /// The bytes [`SamHeader::read`] consumed, as they were.
    pub fn raw(&self) -> &[u8] {
        &self.raw
    }

    /// The header text as htslib holds and writes it: one `\n` per line
    /// (a `\r` before it dropped), the first `@HD` line first.
    pub fn text(&self) -> &[u8] {
        &self.text
    }

    /// Number of header lines; the first data line is line `lines() + 1`.
    pub fn lines(&self) -> u64 {
        self.lines
    }

    /// The `@SQ` references, in tid order.
    pub fn references(&self) -> &[Reference] {
        &self.refs
    }

    /// The tid of a reference name or `AN` alternative name.
    pub fn tid(&self, name: &[u8]) -> Option<i32> {
        self.names.get(name).copied()
    }

    pub(crate) fn names(&self) -> &Arc<NameMap> {
        &self.names
    }
}

/// Splits the fields after a header line's type: each must be a tab and
/// then at least three bytes with `:` as the third (`parse_noncomment_line`).
fn tags(mut rest: &[u8]) -> Option<Vec<&[u8]>> {
    let mut out = Vec::new();
    loop {
        rest = rest.strip_prefix(b"\t")?;
        let end = memchr::memchr(b'\t', rest).unwrap_or(rest.len());
        let tag = &rest[..end];
        if tag.len() < 3 || tag[2] != b':' {
            return None;
        }
        out.push(tag);
        rest = &rest[end..];
        if rest.is_empty() {
            return Some(out);
        }
    }
}

/// The BAM header htslib writes for this SAM header, before BGZF
/// compression: `BAM\1`, `l_text`, the text, `n_ref`, then each
/// reference's `l_name`, NUL-terminated name and `l_ref`.
pub fn bam_header_bytes(header: &SamHeader) -> Vec<u8> {
    let names: usize = header.refs.iter().map(|r| r.name.len() + 9).sum();
    let mut out = Vec::with_capacity(12 + header.text.len() + names);
    out.extend_from_slice(b"BAM\x01");
    out.extend_from_slice(&(header.text.len() as u32).to_le_bytes());
    out.extend_from_slice(&header.text);
    out.extend_from_slice(&(header.refs.len() as i32).to_le_bytes());
    for r in &header.refs {
        out.extend_from_slice(&(r.name.len() as u32 + 1).to_le_bytes());
        out.extend_from_slice(&r.name);
        out.push(0);
        out.extend_from_slice(&(r.len.min(u64::from(u32::MAX)) as u32).to_le_bytes());
    }
    out
}
