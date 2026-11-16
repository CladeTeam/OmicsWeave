//! One SAM line to one BAM record, as htslib's `sam_parse1` reads it and
//! `bam_write1` writes it (`sam.c`).

use std::borrow::Cow;
use std::sync::Arc;

use memchr::memchr;
use omicsweave_core::{Error, Result};

use crate::cnum::{hts_str2int, hts_str2uint, strtod, strtoul_zero};
use crate::header::{NameMap, SamHeader};

/// `BAM_FUNMAP`
const FUNMAP: u16 = 4;
const CSOFT_CLIP: u32 = 4;
const CREF_SKIP: u32 = 3;
/// `HTS_POS_MAX`
const POS_MAX: i64 = ((i32::MAX as i64) << 32) | i32::MAX as i64;
/// `BAM_CIGAR_TYPE`: bit 0 consumes the query, bit 1 the reference.
const CIGAR_TYPE: u32 = 0x3C1A7;
/// Most CIGAR operations a BAM record holds in place; longer CIGARs move
/// to a `CG:B:I` tag behind a two-operation placeholder.
const MAX_INLINE_CIGAR: usize = 0xffff;

/// `bam_cigar_table`
const CIGAR_OP: [i8; 256] = {
    let mut t = [-1i8; 256];
    let ops = b"MIDNSHP=XB";
    let mut i = 0;
    while i < ops.len() {
        t[ops[i] as usize] = i as i8;
        i += 1;
    }
    t
};

/// `seq_nt16_table`
const NT16: [u8; 256] = {
    let mut t = [15u8; 256];
    let codes = b"=ACMGRSVTWYHKDBN";
    let mut i = 0;
    while i < codes.len() {
        t[codes[i] as usize] = i as u8;
        t[codes[i].to_ascii_lowercase() as usize] = i as u8;
        i += 1;
    }
    t
};

/// What a BAM index records for a record: htslib's
/// `hts_idx_push(idx, tid, pos, bam_endpos(b), offset, !(flag & BAM_FUNMAP))`.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct BamSpan {
    /// Reference id; -1 when unplaced.
    pub tid: i32,
    /// 0-based leftmost position; -1 when unplaced.
    pub beg: i64,
    /// `bam_endpos`: `beg` plus the CIGAR's reference length, or plus 1
    /// for an unmapped record or a reference length of 0.
    pub end: i64,
    /// FLAG lacks `0x4`.
    pub mapped: bool,
}

/// Per-worker state [`BamEncoder::encode_line`] reuses between lines.
#[derive(Debug, Default)]
pub struct BamScratch {
    cigar: Vec<u32>,
}

impl BamScratch {
    pub fn new() -> BamScratch {
        BamScratch::default()
    }
}

/// Encodes SAM lines as BAM records against one header. Shared by all
/// workers; each worker brings its own [`BamScratch`].
///
/// The record bytes equal what htslib writes for the line after reading it
/// with `sam_read1` and writing it with `sam_write1` to BAM, and the line
/// is rejected wherever htslib rejects it, in reading or in writing.
/// Deliberate departures from htslib:
///
/// - An RNAME or RNEXT not among the header's `@SQ` names (or their `AN`
///   alternatives) is an error. htslib warns and writes the record as
///   unplaced (the mate as unplaced), which silently drops the name.
/// - A line holding a NUL byte is an error; htslib's C strings end there.
/// - One trailing `\r` is always dropped. htslib drops it only before a
///   `\n`, so a last line that ends in `\r` without a `\n` keeps it there.
///
/// htslib's warnings (empty QNAME, a mapped record without coordinate or
/// CIGAR made unmapped, a CG tag moved into the CIGAR) are followed
/// without being reported. Bytes `>= 0x80` where htslib compares C `char`s
/// (aux tag ids, aux value starts, the skip over junk after a value) are
/// treated as negative, as on x86-64.
#[derive(Debug, Clone)]
pub struct BamEncoder {
    names: Arc<NameMap>,
    n_targets: usize,
}

type Fail = Cow<'static, str>;

#[inline]
fn at(s: &[u8], i: usize) -> u8 {
    s.get(i).copied().unwrap_or(0)
}

/// A C `char` compared as signed.
#[inline]
fn sc(c: u8) -> i8 {
    c as i8
}

fn tab(line: &[u8], p: usize) -> std::result::Result<usize, Fail> {
    memchr(b'\t', &line[p..])
        .map(|k| p + k)
        .ok_or(Cow::Borrowed("too few fields"))
}

#[inline]
fn put_u32(out: &mut [u8], at: usize, v: u32) {
    out[at..at + 4].copy_from_slice(&v.to_le_bytes());
}

fn cigar_lens(cigar: &[u32]) -> (i64, i64) {
    let (mut query, mut reference) = (0, 0);
    for &c in cigar {
        let ty = CIGAR_TYPE >> ((c & 0xf) << 1);
        let len = i64::from(c >> 4);
        if ty & 1 != 0 {
            query += len;
        }
        if ty & 2 != 0 {
            reference += len;
        }
    }
    (query, reference)
}

/// `hts_reg2bin(beg, end, 14, 5)`, truncated to BAM's 16 bits.
fn reg2bin(beg: i64, end: i64) -> u16 {
    let end = end - 1;
    let mut s = 14;
    let mut t: i64 = 4681;
    for l in (1..=5).rev() {
        if beg >> s == end >> s {
            return (t + (beg >> s)) as u16;
        }
        s += 3;
        t -= 1 << (3 * (l - 1));
    }
    0
}

/// `parse_sam_flag`: (flag, bytes consumed).
fn parse_flag(s: &[u8], overflow: &mut bool) -> (u16, usize) {
    match at(s, 0) {
        b'1'..=b'9' => {
            let (v, n, of) = hts_str2uint(s, 16);
            *overflow |= of;
            (v as u16, n)
        }
        b'0' if at(s, 1) == b'\t' => (0, 1),
        b'0' => {
            let (v, n) = strtoul_zero(s);
            if v > 65535 {
                *overflow = true;
                (65535, n)
            } else {
                (v as u16, n)
            }
        }
        _ => (0, 0),
    }
}

impl BamEncoder {
    /// Shares the header's immutable reference-name dictionary. The encoder
    /// owns its share, so the header can be dropped independently.
    pub fn new(header: &SamHeader) -> BamEncoder {
        BamEncoder {
            names: header.names().clone(),
            n_targets: header.references().len(),
        }
    }

    /// Appends the BAM record of one SAM data line, its 4-byte
    /// `block_size` included, to `out`, and returns its index span.
    ///
    /// `line` is one line without its `\n`; one trailing `\r` is dropped.
    /// `line_no` (1-based) only labels errors. On error `out` is left as it
    /// was.
    pub fn encode_line(
        &self,
        scratch: &mut BamScratch,
        line: &[u8],
        line_no: u64,
        out: &mut Vec<u8>,
    ) -> Result<BamSpan> {
        let start = out.len();
        let line = line.strip_suffix(b"\r").unwrap_or(line);
        self.encode(scratch, line, out, start).map_err(|detail| {
            out.truncate(start);
            Error::Format {
                line: line_no,
                detail: detail.into_owned(),
            }
        })
    }

    fn lookup(&self, name: &[u8]) -> std::result::Result<i32, Fail> {
        self.names.get(name).copied().ok_or_else(|| {
            Cow::Owned(format!(
                "reference {:?} is not in the header",
                String::from_utf8_lossy(name)
            ))
        })
    }

    fn encode(
        &self,
        scratch: &mut BamScratch,
        line: &[u8],
        out: &mut Vec<u8>,
        start: usize,
    ) -> std::result::Result<BamSpan, Fail> {
        if memchr(0, line).is_some() {
            return Err("NUL byte in line".into());
        }
        let len = line.len();
        out.reserve(len + 64);
        // block_size and the 32-byte fixed part are filled in at the end.
        out.resize(start + 36, 0);

        // QNAME
        let e = tab(line, 0)?;
        if e > 254 {
            return Err("query name too long".into());
        }
        out.extend_from_slice(&line[..e]);
        out.push(0);
        let l_read_name = e as u32 + 1;
        let mut p = e + 1;

        // FLAG
        let mut overflow = false;
        let (mut flag, n) = parse_flag(&line[p..], &mut overflow);
        p += n;
        if at(line, p) != b'\t' {
            return Err("malformed FLAG".into());
        }
        p += 1;

        // RNAME
        let e = tab(line, p)?;
        let mut tid = if &line[p..e] == b"*" {
            -1
        } else {
            if self.n_targets == 0 {
                return Err("no SQ lines present in the header".into());
            }
            self.lookup(&line[p..e])?
        };
        p = e + 1;

        // POS
        let (v, n, of) = hts_str2uint(&line[p..], 62);
        overflow |= of;
        let pos = v as i64 - 1;
        p += n;
        if at(line, p) != b'\t' {
            return Err("malformed POS".into());
        }
        p += 1;
        if pos < 0 && tid >= 0 {
            tid = -1;
        }
        if tid < 0 {
            flag |= FUNMAP;
        }

        // MAPQ
        let (v, n, of) = hts_str2uint(&line[p..], 8);
        overflow |= of;
        let mapq = v as u32;
        p += n;
        if at(line, p) != b'\t' {
            return Err("malformed MAPQ".into());
        }
        p += 1;

        // CIGAR
        let cigar = &mut scratch.cigar;
        cigar.clear();
        let (mut qlen, mut rlen) = (0, 0);
        if at(line, p) != b'*' {
            let field = memchr(b'\t', &line[p..]).map_or(len, |k| p + k);
            let n_ops = line[p..field]
                .iter()
                .filter(|c| !c.is_ascii_digit())
                .count();
            if n_ops == 0 {
                return Err("no CIGAR operations".into());
            }
            // Bound speculative allocation before operation syntax has been
            // validated; longer CIGARs can grow through the normal Vec path.
            cigar.reserve(n_ops.min(MAX_INLINE_CIGAR + 1));
            for _ in 0..n_ops {
                let (v, n, of) = hts_str2uint(&line[p..], 28);
                if n == 0 {
                    return Err("invalid CIGAR length".into());
                }
                if of {
                    return Err("CIGAR length too long".into());
                }
                p += n;
                let op = CIGAR_OP[at(line, p) as usize];
                p += 1;
                if op < 0 {
                    return Err("unrecognized CIGAR operator".into());
                }
                cigar.push((v as u32) << 4 | op as u32);
                // Both lengths are needed immediately after parsing. Keep
                // them with the validated operation rather than traversing
                // the scratch CIGAR again.
                let ty = CIGAR_TYPE >> ((op as u32) << 1);
                if ty & 1 != 0 {
                    qlen += v as i64;
                }
                if ty & 2 != 0 {
                    rlen += v as i64;
                }
            }
            if at(line, p) != b'\t' {
                return Err("malformed CIGAR".into());
            }
            p += 1;
        } else {
            flag |= FUNMAP;
            p = tab(line, p)? + 1;
        }
        let reflen = if flag & FUNMAP != 0 { 1 } else { rlen.max(1) };
        if POS_MAX - reflen <= pos {
            return Err("read ends beyond highest supported position".into());
        }

        // RNEXT
        let e = tab(line, p)?;
        let mut mtid = match &line[p..e] {
            b"=" => tid,
            b"*" => -1,
            name => self.lookup(name)?,
        };
        p = e + 1;

        // PNEXT
        let (v, n, of) = hts_str2uint(&line[p..], 62);
        overflow |= of;
        let mpos = v as i64 - 1;
        p += n;
        if at(line, p) != b'\t' {
            return Err("malformed PNEXT".into());
        }
        p += 1;
        if mpos < 0 && mtid >= 0 {
            mtid = -1;
        }

        // TLEN
        let (tlen, n, of) = hts_str2int(&line[p..], 63);
        overflow |= of;
        p += n;
        if at(line, p) != b'\t' {
            return Err("malformed TLEN".into());
        }
        p += 1;
        if overflow {
            return Err("number outside allowed range".into());
        }

        // The CIGAR's place; filled once the final CIGAR is known.
        let cigar_at = out.len();
        let n_cigar = cigar.len();
        let region = if n_cigar > MAX_INLINE_CIGAR {
            8
        } else {
            4 * n_cigar
        };
        out.resize(cigar_at + region, 0);

        // SEQ
        let e = tab(line, p)?;
        let seq = &line[p..e];
        p = e + 1;
        let l_qseq = if seq == b"*" { 0 } else { seq.len() };
        if l_qseq > i32::MAX as usize {
            return Err("read sequence is too long".into());
        }
        if l_qseq > 0 {
            let pairs = seq.chunks_exact(2);
            let last = pairs.remainder();
            out.extend(pairs.map(|b| NT16[b[0] as usize] << 4 | NT16[b[1] as usize]));
            if let [b] = last {
                out.push(NT16[*b as usize] << 4);
            }
        }

        // QUAL
        if at(line, p) == b'*' && matches!(at(line, p + 1), b'\t' | 0) {
            out.resize(out.len() + l_qseq, 0xff);
            p += 2;
        } else {
            if len - p < l_qseq || !matches!(at(line, p + l_qseq), b'\t' | 0) {
                return Err("SEQ and QUAL are of different length".into());
            }
            let o = out.len();
            out.resize(o + l_qseq, 0);
            let mut bad = 0u8;
            for (d, &s) in out[o..].iter_mut().zip(&line[p..p + l_qseq]) {
                let v = s.wrapping_sub(33);
                bad |= v;
                *d = v;
            }
            if bad & 0x80 != 0 {
                return Err("invalid QUAL character".into());
            }
            p += l_qseq + 1;
        }

        // Optional fields
        // Only the soft-clip placeholder can use CG; specialize the aux
        // parser so ordinary records need no additional tag-id checks.
        let cg_at = if n_cigar > 0
            && cigar[0] == (CSOFT_CLIP | (l_qseq as u32) << 4)
            && tid >= 0
            && pos >= 0
        {
            aux::<true>(line, p, out)?
        } else {
            aux::<false>(line, p, out)?
        };

        // bam_tag2cigar: a CIGAR of one soft clip spanning the read and a
        // CG:B:I tag holding at least as many operations stand for the CG
        // tag's CIGAR.
        let mut cg_moved = false;
        if let Some(tag_at) = cg_at {
            let tag = &out[tag_at..];
            let typed = tag[2] == b'B' && matches!(tag[3], b'I' | b'i');
            let count = if typed {
                u32::from_le_bytes(tag[4..8].try_into().unwrap()) as usize
            } else {
                0
            };
            if typed && count >= n_cigar && count < 1 << 29 {
                let end = tag_at + 8 + 4 * count;
                cigar.clear();
                cigar.extend(
                    out[tag_at + 8..end]
                        .chunks_exact(4)
                        .map(|c| u32::from_le_bytes(c.try_into().unwrap())),
                );
                let new_region = if count > MAX_INLINE_CIGAR {
                    8
                } else {
                    4 * count
                };
                // The replacement region is smaller than the removed CG
                // tag, so its middle destination ends before the tail's
                // source. Move each retained byte once, then truncate;
                // the final CIGAR write fills the replacement region.
                let middle = cigar_at + region..tag_at;
                let middle_len = middle.len();
                let middle_at = cigar_at + new_region;
                out.copy_within(middle, middle_at);
                let tail_at = middle_at + middle_len;
                let new_len = tail_at + (out.len() - end);
                out.copy_within(end.., tail_at);
                out.truncate(new_len);
                (qlen, rlen) = cigar_lens(cigar);
                cg_moved = true;
            }
        }
        if l_qseq > 0 && !cigar.is_empty() && qlen != l_qseq as i64 {
            return Err("CIGAR and query sequence are of different length".into());
        }

        // bam_write1
        let n_cigar = cigar.len();
        let reflen = if !cg_moved {
            reflen
        } else if flag & FUNMAP != 0 {
            1
        } else {
            rlen.max(1)
        };
        if pos > i64::from(i32::MAX)
            || mpos > i64::from(i32::MAX)
            || tlen < i64::from(i32::MIN)
            || tlen > i64::from(i32::MAX)
        {
            return Err("positional data is too large for BAM format".into());
        }
        if n_cigar > MAX_INLINE_CIGAR {
            if rlen >= 1 << 28 {
                return Err("CIGAR reference length too long for BAM".into());
            }
            put_u32(out, cigar_at, (l_qseq as u32) << 4 | CSOFT_CLIP);
            put_u32(out, cigar_at + 4, (rlen as u32) << 4 | CREF_SKIP);
            out.extend_from_slice(b"CGBI");
            out.extend_from_slice(&(n_cigar as u32).to_le_bytes());
            out.extend(cigar.iter().flat_map(|c| c.to_le_bytes()));
        } else {
            for (d, c) in out[cigar_at..cigar_at + 4 * n_cigar]
                .chunks_exact_mut(4)
                .zip(cigar.iter())
            {
                d.copy_from_slice(&c.to_le_bytes());
            }
        }

        let n_cigar_field = if n_cigar > MAX_INLINE_CIGAR {
            2
        } else {
            n_cigar as u32
        };
        let bin = u32::from(reg2bin(pos, pos + reflen));
        let block = (out.len() - start - 4) as u32;
        let fixed = [
            block,
            tid as u32,
            pos as u32,
            bin << 16 | mapq << 8 | l_read_name,
            u32::from(flag) << 16 | n_cigar_field,
            l_qseq as u32,
            mtid as u32,
            mpos as u32,
            tlen as u32,
        ];
        for (i, v) in fixed.into_iter().enumerate() {
            put_u32(out, start + 4 * i, v);
        }
        Ok(BamSpan {
            tid,
            beg: pos,
            end: pos + reflen,
            mapped: flag & FUNMAP == 0,
        })
    }
}

/// `aux_parse` (strict): the optional fields from `q` to the line's end.
/// When requested, remember the first CG tag's absolute output offset;
/// a differently typed first CG still shadows all subsequent CG tags.
fn aux<const FIND_CG: bool>(
    line: &[u8],
    mut q: usize,
    out: &mut Vec<u8>,
) -> std::result::Result<Option<usize>, Fail> {
    let p = line.len();
    let mut overflow = false;
    let mut cg_at = None;
    while q < p {
        if p - q < 5 {
            return Err("incomplete aux field".into());
        }
        if sc(line[q]) < sc(b'!') || sc(line[q + 1]) < sc(b'!') {
            return Err("invalid aux tag id".into());
        }
        if FIND_CG && cg_at.is_none() && &line[q..q + 2] == b"CG" {
            cg_at = Some(out.len());
        }
        out.extend_from_slice(&line[q..q + 2]);
        let ty = line[q + 3];
        q += 5;
        if ty != b'Z' && ty != b'H' && sc(at(line, q)) <= sc(b'\t') {
            return Err("incomplete aux field".into());
        }
        match ty {
            b'A' | b'a' | b'c' | b'C' => {
                out.push(b'A');
                out.push(line[q]);
                q += 1;
            }
            b'i' | b'I' => {
                if line[q] == b'-' {
                    let (x, n, of) = hts_str2int(&line[q..], 32);
                    overflow |= of;
                    q += n;
                    let x = x as i32;
                    if x >= i32::from(i8::MIN) {
                        out.push(b'c');
                        out.push(x as u8);
                    } else if x >= i32::from(i16::MIN) {
                        out.push(b's');
                        out.extend_from_slice(&(x as i16).to_le_bytes());
                    } else {
                        out.push(b'i');
                        out.extend_from_slice(&x.to_le_bytes());
                    }
                } else {
                    let (x, n, of) = hts_str2uint(&line[q..], 32);
                    overflow |= of;
                    q += n;
                    let x = x as u32;
                    if x <= u32::from(u8::MAX) {
                        out.push(b'C');
                        out.push(x as u8);
                    } else if x <= u32::from(u16::MAX) {
                        out.push(b'S');
                        out.extend_from_slice(&(x as u16).to_le_bytes());
                    } else {
                        out.push(b'I');
                        out.extend_from_slice(&x.to_le_bytes());
                    }
                }
            }
            b'f' => {
                let (v, n, _) = strtod(&line[q..]);
                q += n;
                out.push(b'f');
                out.extend_from_slice(&(v as f32).to_le_bytes());
            }
            b'd' => {
                let (v, n, _) = strtod(&line[q..]);
                q += n;
                out.push(b'd');
                out.extend_from_slice(&v.to_le_bytes());
            }
            b'Z' | b'H' => {
                let end = memchr(b'\t', &line[q..]).map_or(p, |k| q + k);
                if ty == b'H' && (end - q) & 1 != 0 {
                    return Err("hex field does not have an even number of digits".into());
                }
                out.push(ty);
                out.extend_from_slice(&line[q..end]);
                out.push(0);
                q = end;
            }
            b'B' => {
                let sub = line[q];
                q += 1;
                if !matches!(at(line, q), 0 | b',' | b'\t') {
                    return Err("B aux field type not followed by ','".into());
                }
                q = b_array(line, sub, q, out, &mut 0)?;
            }
            _ => {
                return Err(Cow::Owned(format!(
                    "unrecognized aux type {:?}",
                    char::from(ty)
                )));
            }
        }
        // Skip any junk to the next tab.
        while sc(at(line, q)) > sc(b'\t') {
            q += 1;
        }
        q += 1;
    }
    if overflow {
        return Err("numeric value out of allowed range".into());
    }
    Ok(cg_at)
}

/// `skip_to_comma_`
#[inline]
fn skip_to_comma(line: &[u8], q: &mut usize) {
    while sc(at(line, *q)) > sc(b'\t') && at(line, *q) != b',' {
        *q += 1;
    }
}

/// `sam_parse_B_vals_r`: a B array of subtype `sub` whose values start at
/// `input` (the first comma). Values out of the subtype's range make
/// htslib retry once with the smallest subtype holding them all.
fn b_array(
    line: &[u8],
    sub: u8,
    input: usize,
    out: &mut Vec<u8>,
    tries: &mut u32,
) -> std::result::Result<usize, Fail> {
    *tries += 1;
    if *tries > 2 {
        return Err("malformed data in B array".into());
    }
    if !matches!(sub, b'A' | b'c' | b'C' | b's' | b'S' | b'i' | b'I' | b'f') {
        return Err(Cow::Owned(format!(
            "unrecognized type B:{:?}",
            char::from(sub)
        )));
    }
    let orig = out.len();
    out.push(b'B');
    out.push(sub);
    let count_at = out.len();
    out.extend_from_slice(&[0; 4]);
    let mut q = input;
    let mut used: u32 = 0;
    let mut overflow = false;
    // Signed subtypes parse every value; unsigned ones flag a minus sign
    // as an overflow without parsing it.
    macro_rules! values {
        (|$s:ident| $parse:expr, $unsigned:expr) => {
            while at(line, q) == b',' {
                used = used.wrapping_add(1);
                if $unsigned && at(line, q + 1) == b'-' {
                    overflow = true;
                    q += 1;
                    skip_to_comma(line, &mut q);
                } else {
                    let $s = &line[q + 1..];
                    let n: usize = $parse;
                    q += 1 + n;
                }
            }
        };
    }
    match sub {
        b'c' => values!(
            |s| {
                let (v, n, of) = hts_str2int(s, 8);
                overflow |= of;
                out.push(v as u8);
                n
            },
            false
        ),
        b'C' => values!(
            |s| {
                let (v, n, of) = hts_str2uint(s, 8);
                overflow |= of;
                out.push(v as u8);
                n
            },
            true
        ),
        b's' => values!(
            |s| {
                let (v, n, of) = hts_str2int(s, 16);
                overflow |= of;
                out.extend_from_slice(&(v as i16).to_le_bytes());
                n
            },
            false
        ),
        b'S' => values!(
            |s| {
                let (v, n, of) = hts_str2uint(s, 16);
                overflow |= of;
                out.extend_from_slice(&(v as u16).to_le_bytes());
                n
            },
            true
        ),
        b'i' => values!(
            |s| {
                let (v, n, of) = hts_str2int(s, 32);
                overflow |= of;
                out.extend_from_slice(&(v as i32).to_le_bytes());
                n
            },
            false
        ),
        b'I' => values!(
            |s| {
                let (v, n, of) = hts_str2uint(s, 32);
                overflow |= of;
                out.extend_from_slice(&(v as u32).to_le_bytes());
                n
            },
            true
        ),
        b'f' => values!(
            |s| {
                let (v, n, _) = strtod(s);
                out.extend_from_slice(&(v as f32).to_le_bytes());
                n
            },
            false
        ),
        // B:A is a valid, necessarily empty, array to htslib.
        _ => {}
    }
    if !matches!(at(line, q), b'\t' | 0) {
        return Err(Cow::Owned(format!("malformed B:{}", char::from(sub))));
    }
    put_u32(out, count_at, used);
    if !overflow {
        return Ok(q);
    }

    out.truncate(orig);
    let end = q;
    let mut q = input;
    let (mut min, mut max) = (0i64, 0i64);
    let mut overflow = false;
    while q < end {
        let (v, n, of) = hts_str2int(&line[q + 1..], 64);
        overflow |= of;
        q += 1 + n;
        max = max.max(v);
        min = min.min(v);
        skip_to_comma(line, &mut q);
    }
    if !overflow {
        let fit = |lo: i64, hi: i64| min >= lo && max <= hi;
        let retry = if min < 0 {
            if fit(i8::MIN.into(), i8::MAX.into()) {
                Some(b'c')
            } else if fit(i16::MIN.into(), i16::MAX.into()) {
                Some(b's')
            } else if fit(i32::MIN.into(), i32::MAX.into()) {
                Some(b'i')
            } else {
                None
            }
        } else if max < i64::from(u8::MAX) {
            Some(b'C')
        } else if max <= i64::from(u16::MAX) {
            Some(b'S')
        } else if max <= i64::from(u32::MAX) {
            Some(b'I')
        } else {
            None
        };
        if let Some(t) = retry {
            return b_array(line, t, input, out, tries);
        }
    }
    Err("numeric value in B array out of allowed range".into())
}
