//! VCF text to BCF: the header and each record encoded into the
//! uncompressed bytes htslib writes (`bcf_hdr_write`, `bcf_write`) after
//! reading the same text (`vcf_hdr_read`, `vcf_parse`). BGZF framing,
//! compression and indexing are left to the caller; [`BcfSpan`] is what an
//! index needs of each record.
//!
//! The header is fixed before any record is encoded, so unlike htslib a
//! record that names a contig, FILTER, INFO or FORMAT key the header does
//! not declare is an error: htslib adds a dummy declaration to its header
//! and refuses to write that record, but writes later records that use the
//! new name against a header it has already written, which leaves a BCF
//! whose ids its own header does not define.

use std::io::BufRead;

use omicsweave_core::{Error, Result};

use crate::fields::{
    BT_CHAR, BT_INT8, BT_INT16, BT_INT32, BT_NULL, INT32_MISSING, INT32_VECTOR_END,
};
use crate::header::Header;
use crate::record::Record;
use crate::text::{Declared, Scratch, parse_line};

/// Reads a VCF text header from `input`: lines up to and including the
/// first line that starts with `#` but not `##` (the `#CHROM` line), as
/// htslib's `vcf_hdr_read` collects them, parsed by
/// [`Header::from_vcf_text`]. Data lines after it are left in `input`.
///
/// The text must start with `##fileformat=VCF`, as htslib needs to detect
/// VCF. A line ends at `\n`, and a `\r` before it is dropped; empty lines
/// are skipped. Returns the header and the number of lines consumed, so
/// the first data line is line `consumed + 1`.
pub fn read_vcf_header<R: BufRead + ?Sized>(input: &mut R) -> Result<(Header, u64)> {
    let mut text = Vec::new();
    let mut line = Vec::new();
    let mut n = 0u64;
    loop {
        line.clear();
        if input.read_until(b'\n', &mut line)? == 0 {
            return Err(Error::invalid("vcf header", "no #CHROM line"));
        }
        n += 1;
        if n == 1 && !line.starts_with(b"##fileformat=VCF") {
            return Err(Error::invalid(
                "vcf header",
                "the text does not start with ##fileformat=VCF, which htslib needs to detect VCF",
            ));
        }
        if line.last() == Some(&b'\n') {
            line.pop();
        }
        if line.last() == Some(&b'\r') {
            line.pop();
        }
        if line.is_empty() {
            continue;
        }
        if line[0] != b'#' {
            return Err(Error::Format {
                line: n,
                detail: "header line without '#'".into(),
            });
        }
        text.extend_from_slice(&line);
        text.push(b'\n');
        if line.get(1) != Some(&b'#') {
            break;
        }
    }
    Ok((Header::from_vcf_text(&text)?, n))
}

/// The bytes `bcf_hdr_write` writes, before BGZF compression: the magic
/// `BCF\2\2`, the text length (little-endian `u32`, counting the NUL) and
/// the header text as htslib formats it (with `IDX`), NUL-terminated.
///
/// The header must come from [`Header::from_vcf_text`] or
/// [`read_vcf_header`], which keep the records htslib would print.
pub fn bcf_header_bytes(header: &Header) -> Result<Vec<u8>> {
    let text = header.bcf_text.as_deref().ok_or(Error::Unsupported(
        "BCF header bytes of a header not parsed from VCF text",
    ))?;
    let len = u32::try_from(text.len() + 1)
        .map_err(|_| Error::invalid("bcf header", "header text too long"))?;
    let mut out = Vec::with_capacity(text.len() + 10);
    out.extend_from_slice(b"BCF\x02\x02");
    out.extend_from_slice(&len.to_le_bytes());
    out.extend_from_slice(text);
    out.push(0);
    Ok(out)
}

/// What an index records of one BCF record, as `bcf_write` hands it to
/// `hts_idx_push`: the contig id and the 0-based span `[beg, end)`, with
/// `end = beg + rlen` (`rlen` from REF, INFO/END, INFO/SVLEN and FORMAT/LEN
/// as `vcf_parse` computes it). `beg` is -1 for a record at POS 0;
/// `hts_idx_push` itself clamps a negative `beg` and an empty span.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct BcfSpan {
    pub rid: u32,
    pub beg: i64,
    pub end: i64,
}

/// Per-worker state of [`BcfEncoder::encode_line`]: buffers reused across
/// lines so that encoding allocates only while they grow.
#[derive(Debug, Default)]
pub struct BcfScratch {
    rec: Record,
    sc: Scratch,
}

/// Encodes VCF data lines into BCF records against a fixed header. Shared
/// by reference across threads, each with its own [`BcfScratch`].
#[derive(Debug, Clone)]
pub struct BcfEncoder {
    header: Header,
}

/// `bcf_enc_size`: a typed-value descriptor.
fn put_size(out: &mut Vec<u8>, size: usize, ty: u8) {
    if size < 15 {
        out.push((size as u8) << 4 | ty);
    } else if size < 128 {
        out.extend_from_slice(&[15 << 4 | ty, 1 << 4 | BT_INT8, size as u8]);
    } else if size < 32768 {
        out.extend_from_slice(&[15 << 4 | ty, 1 << 4 | BT_INT16]);
        out.extend_from_slice(&(size as i16).to_le_bytes());
    } else {
        out.extend_from_slice(&[15 << 4 | ty, 1 << 4 | BT_INT32]);
        out.extend_from_slice(&(size as i32).to_le_bytes());
    }
}

/// `bcf_enc_int1`: one integer in the narrowest type.
fn put_int1(out: &mut Vec<u8>, x: i32) {
    match x {
        INT32_VECTOR_END => out.extend_from_slice(&[1 << 4 | BT_INT8, i8::MIN as u8 + 1]),
        INT32_MISSING => out.extend_from_slice(&[1 << 4 | BT_INT8, i8::MIN as u8]),
        -120..=127 => out.extend_from_slice(&[1 << 4 | BT_INT8, x as u8]),
        -32760..=32767 => {
            out.push(1 << 4 | BT_INT16);
            out.extend_from_slice(&(x as i16).to_le_bytes());
        }
        _ => {
            out.push(1 << 4 | BT_INT32);
            out.extend_from_slice(&x.to_le_bytes());
        }
    }
}

/// `bcf_enc_vint` of a FILTER list (`wsize` is the whole list).
fn put_vint(out: &mut Vec<u8>, a: &[i32]) {
    match a {
        [] => return put_size(out, 0, BT_NULL),
        [x] => return put_int1(out, *x),
        _ => {}
    }
    let max = a.iter().copied().max().unwrap();
    let min = a
        .iter()
        .copied()
        .filter(|&v| v > INT32_VECTOR_END)
        .min()
        .unwrap_or(i32::MAX);
    let narrow = |v: i32, missing: i32| match v {
        INT32_VECTOR_END => missing + 1,
        INT32_MISSING => missing,
        v => v,
    };
    if max <= 127 && min >= -120 {
        put_size(out, a.len(), BT_INT8);
        out.extend(a.iter().map(|&v| narrow(v, i8::MIN.into()) as u8));
    } else if max <= 32767 && min >= -32760 {
        put_size(out, a.len(), BT_INT16);
        for &v in a {
            out.extend_from_slice(&(narrow(v, i16::MIN.into()) as i16).to_le_bytes());
        }
    } else {
        put_size(out, a.len(), BT_INT32);
        for &v in a {
            out.extend_from_slice(&v.to_le_bytes());
        }
    }
}

fn put_chars(out: &mut Vec<u8>, s: &[u8]) {
    put_size(out, s.len(), BT_CHAR);
    out.extend_from_slice(s);
}

impl BcfEncoder {
    /// An encoder for records of `header`, which should be the one whose
    /// [`bcf_header_bytes`] precede them.
    pub fn new(header: &Header) -> BcfEncoder {
        BcfEncoder {
            header: header.clone(),
        }
    }

    pub fn header(&self) -> &Header {
        &self.header
    }

    /// Encodes one VCF data line into a complete BCF record (`l_shared`,
    /// `l_indiv`, the fixed fields, the shared and the per-sample blocks)
    /// appended to `out`, and returns its index span.
    ///
    /// `line` excludes the `\n`; a trailing `\r` is dropped. `line_no` is
    /// the 1-based line number reported in errors. Values are encoded as
    /// `vcf_parse` encodes them. Errors are [`Error::Format`] and leave
    /// `out` unchanged: lines htslib fails to parse, and records it parses
    /// but `bcf_write` refuses: a name the header does not declare (see
    /// the module documentation), a repeated FORMAT key, a header with
    /// samples and no FORMAT column, POS beyond 2^31 - 1. A line holding a
    /// NUL byte is also an error.
    pub fn encode_line(
        &self,
        scratch: &mut BcfScratch,
        line: &[u8],
        line_no: u64,
        out: &mut Vec<u8>,
    ) -> Result<BcfSpan> {
        let format_error = |detail: String| Error::Format {
            line: line_no,
            detail,
        };
        let line = line.strip_suffix(b"\r").unwrap_or(line);
        if line.contains(&0) {
            return Err(format_error("NUL byte in a data line".into()));
        }
        if line.len() > i32::MAX as usize {
            return Err(format_error("data line too long".into()));
        }
        let BcfScratch { rec, sc } = scratch;
        parse_line(&mut Declared(&self.header), line, rec, sc)
            .map_err(|e| format_error(e.to_string()))?;
        if rec.pos >= i64::from(i32::MAX) {
            return Err(format_error(
                "position beyond 2^31 - 1, which BCF cannot hold".into(),
            ));
        }

        let start = out.len();
        out.extend_from_slice(&[0; 32]);
        let d = &rec.data;
        put_chars(out, if d.id == b"." { b"" } else { &d.id });
        for allele in rec.alleles() {
            put_chars(out, allele);
        }
        put_vint(out, &d.filters);
        for s in &d.info {
            put_int1(out, s.key as i32);
            put_size(out, s.n as usize, s.ty);
            out.extend_from_slice(&d.info_bytes[s.start..s.end]);
        }
        let shared_end = out.len();
        for s in &d.fmt {
            put_int1(out, s.key as i32);
            put_size(out, s.n as usize, s.ty);
            out.extend_from_slice(&d.fmt_bytes[s.start..s.end]);
        }
        let l_shared = u32::try_from(shared_end - start - 8);
        let l_indiv = u32::try_from(out.len() - shared_end);
        let (Ok(l_shared), Ok(l_indiv)) = (l_shared, l_indiv) else {
            out.truncate(start);
            return Err(format_error("record too large for BCF".into()));
        };

        let mut fixed = [0u8; 32];
        fixed[0..4].copy_from_slice(&l_shared.to_le_bytes());
        fixed[4..8].copy_from_slice(&l_indiv.to_le_bytes());
        fixed[8..12].copy_from_slice(&rec.rid.to_le_bytes());
        fixed[12..16].copy_from_slice(&(rec.pos as u32).to_le_bytes());
        fixed[16..20].copy_from_slice(&(rec.rlen as u32).to_le_bytes());
        fixed[20..24].copy_from_slice(&d.qual.to_le_bytes());
        fixed[24..26].copy_from_slice(&(d.info.len() as u16).to_le_bytes());
        fixed[26..28].copy_from_slice(&(rec.n_alleles() as u16).to_le_bytes());
        let n_fmt_sample = (d.fmt.len() as u32) << 24 | (d.n_sample & 0xff_ffff);
        fixed[28..32].copy_from_slice(&n_fmt_sample.to_le_bytes());
        out[start..start + 32].copy_from_slice(&fixed);

        Ok(BcfSpan {
            rid: rec.rid,
            beg: rec.pos,
            end: rec.pos + rec.rlen,
        })
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn sizes_use_htslib_descriptors() {
        let enc = |n| {
            let mut v = Vec::new();
            put_size(&mut v, n, BT_CHAR);
            v
        };
        assert_eq!(enc(14), [0xe7]);
        assert_eq!(enc(15), [0xf7, 0x11, 15]);
        assert_eq!(enc(128), [0xf7, 0x12, 128, 0]);
        assert_eq!(enc(40000), [0xf7, 0x13, 0x40, 0x9c, 0, 0]);
    }
}
