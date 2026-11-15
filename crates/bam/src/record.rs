//! One BAM record read in place, checked as htslib's `bam_read1` and
//! `sam_read1` check it.

use std::ops::Range;

use omicsweave_core::{AlignmentEntry, Error, Result};

/// FLAG bits (`BAM_F*`).
pub mod flag {
    pub const PAIRED: u16 = 0x1;
    pub const PROPER_PAIR: u16 = 0x2;
    pub const UNMAPPED: u16 = 0x4;
    pub const MATE_UNMAPPED: u16 = 0x8;
    pub const REVERSE: u16 = 0x10;
    pub const MATE_REVERSE: u16 = 0x20;
    pub const READ1: u16 = 0x40;
    pub const READ2: u16 = 0x80;
    pub const SECONDARY: u16 = 0x100;
    pub const QC_FAIL: u16 = 0x200;
    pub const DUPLICATE: u16 = 0x400;
    pub const SUPPLEMENTARY: u16 = 0x800;
}

/// CIGAR operation codes (`BAM_C*`), the low 4 bits of an operation.
pub mod op {
    pub const MATCH: u8 = 0;
    pub const INS: u8 = 1;
    pub const DEL: u8 = 2;
    pub const REF_SKIP: u8 = 3;
    pub const SOFT_CLIP: u8 = 4;
    pub const HARD_CLIP: u8 = 5;
    pub const PAD: u8 = 6;
    pub const EQUAL: u8 = 7;
    pub const DIFF: u8 = 8;
}

/// `BAM_CIGAR_TYPE`: bit 0 consumes the query, bit 1 the reference.
const CIGAR_TYPE: u32 = 0x3C1A7;

/// Whether operation `op` consumes reference bases.
#[inline]
pub fn consumes_reference(op: u8) -> bool {
    (CIGAR_TYPE >> ((op as u32 & 0xf) << 1)) & 2 != 0
}

/// Whether operation `op` consumes query bases.
#[inline]
pub fn consumes_query(op: u8) -> bool {
    (CIGAR_TYPE >> ((op as u32 & 0xf) << 1)) & 1 != 0
}

/// `seq_nt16_str`: the base of each 4-bit code.
pub const NT16: &[u8; 16] = b"=ACMGRSVTWYHKDBN";

const fn base_pairs() -> [[u8; 2]; 256] {
    let mut table = [[0; 2]; 256];
    let mut i = 0;
    while i < table.len() {
        table[i] = [NT16[i >> 4], NT16[i & 15]];
        i += 1;
    }
    table
}

const BASE_PAIRS: [[u8; 2]; 256] = base_pairs();

#[inline]
fn u32_at(b: &[u8], at: usize) -> u32 {
    u32::from_le_bytes(b[at..at + 4].try_into().unwrap())
}

#[inline]
fn i32_at(b: &[u8], at: usize) -> i32 {
    i32::from_le_bytes(b[at..at + 4].try_into().unwrap())
}

/// `aux_type2size`.
fn aux_size(t: u8) -> u64 {
    match t {
        b'A' | b'c' | b'C' => 1,
        b's' | b'S' => 2,
        b'i' | b'I' | b'f' => 4,
        b'd' => 8,
        b'Z' | b'H' | b'B' => t as u64,
        _ => 0,
    }
}

/// `skip_aux`: the end of the value whose type byte is at `s`, `None` when
/// it is malformed.
fn skip_aux(aux: &[u8], s: usize) -> Option<usize> {
    let end = aux.len();
    if s >= end {
        return Some(end);
    }
    let size = aux_size(aux[s]);
    let s = s + 1;
    match size {
        0 => None,
        90 | 72 => Some(memchr(0, &aux[s..]).map_or(end, |p| s + p + 1)),
        66 => {
            if end - s < 5 {
                return None;
            }
            let size = aux_size(aux[s]);
            let n = u32_at(aux, s + 1) as u64;
            let s = s + 5;
            if size == 0 || ((end - s) as u64) < size * n {
                return None;
            }
            Some(s + (size * n) as usize)
        }
        _ => {
            if ((end - s) as u64) < size {
                return None;
            }
            Some(s + size as usize)
        }
    }
}

fn memchr(c: u8, b: &[u8]) -> Option<usize> {
    b.iter().position(|&x| x == c)
}

/// `bam_aux_get`: the position of `tag`'s type byte in `aux`. A malformed
/// field before it, or in it, is an error.
fn aux_get(aux: &[u8], tag: [u8; 2]) -> std::result::Result<Option<usize>, ()> {
    if aux.len() <= 2 {
        return Ok(None);
    }
    let mut s = 2;
    loop {
        if aux[s - 2] == tag[0] && aux[s - 1] == tag[1] {
            let e = skip_aux(aux, s).ok_or(())?;
            if (aux[s] == b'Z' || aux[s] == b'H') && aux[e - 1] != 0 {
                return Err(());
            }
            return Ok(Some(s));
        }
        let next = skip_aux(aux, s).ok_or(())?;
        if aux.len() - next <= 2 {
            return Ok(None);
        }
        s = next + 2;
    }
}

/// One record: the bytes after its `block_size`, with its fields located.
#[derive(Debug, Clone)]
pub struct Record<'a> {
    data: &'a [u8],
    pub tid: i32,
    pub pos: i32,
    pub mapq: u8,
    pub flag: u16,
    pub mate_tid: i32,
    pub mate_pos: i32,
    pub tlen: i32,
    name: Range<usize>,
    /// The operations in use: in place, or the `CG` tag's array.
    cigar: Range<usize>,
    seq: Range<usize>,
    l_seq: usize,
    qual: Range<usize>,
    aux: Range<usize>,
    /// The `CG` tag the CIGAR came from, which is not part of the record's
    /// optional fields (htslib moves it into the CIGAR).
    cg: Option<Range<usize>>,
    end: i64,
}

fn invalid(detail: impl Into<String>) -> Error {
    Error::invalid("BAM record", detail)
}

/// Reference and query lengths of `cigar` bytes.
fn cigar_lens(cigar: &[u8]) -> (i64, i64) {
    let (mut r, mut q) = (0i64, 0i64);
    for c in cigar.chunks_exact(4) {
        let v = u32::from_le_bytes(c.try_into().unwrap());
        let (op, len) = ((v & 0xf) as u8, (v >> 4) as i64);
        if consumes_reference(op) {
            r += len;
        }
        if consumes_query(op) {
            q += len;
        }
    }
    (r, q)
}

impl<'a> Record<'a> {
    /// Locates the fields of `data` (a record after its `block_size`) of a
    /// file with `n_refs` references, rejecting what htslib rejects.
    pub fn parse(data: &'a [u8], n_refs: usize) -> Result<Record<'a>> {
        if data.len() < 32 {
            return Err(invalid(format!(
                "block_size {} below the 32 bytes of fixed fields",
                data.len()
            )));
        }
        let tid = i32_at(data, 0);
        let pos = i32_at(data, 4);
        let x2 = u32_at(data, 8);
        let l_qname = (x2 & 0xff) as usize;
        let mapq = (x2 >> 8) as u8;
        let x3 = u32_at(data, 12);
        let flag = (x3 >> 16) as u16;
        let n_cigar = (x3 & 0xffff) as usize;
        let l_seq = i32_at(data, 16);
        let mate_tid = i32_at(data, 20);
        let mate_pos = i32_at(data, 24);
        let tlen = i32_at(data, 28);
        if l_seq < 0 || l_qname < 1 {
            return Err(invalid(format!(
                "read name length {l_qname}, sequence length {l_seq}"
            )));
        }
        let l_seq = l_seq as usize;
        let body = data.len() - 32;
        let need =
            (n_cigar as u64) * 4 + l_qname as u64 + (l_seq as u64).div_ceil(2) + l_seq as u64;
        if need > body as u64 {
            return Err(invalid(format!(
                "fields need {need} bytes, record has {body}"
            )));
        }
        let n = n_refs as i64;
        let out_of_range = |t: i32| (t as i64) < -1 || t as i64 >= n;
        if out_of_range(tid) || out_of_range(mate_tid) {
            return Err(invalid(format!(
                "reference {tid} or mate reference {mate_tid} outside the header's {n_refs}"
            )));
        }
        let name_at = 32;
        let raw_name = &data[name_at..name_at + l_qname];
        let name_len = memchr(0, raw_name).unwrap_or(l_qname);
        let cigar_at = name_at + l_qname;
        let seq_at = cigar_at + 4 * n_cigar;
        let qual_at = seq_at + l_seq.div_ceil(2);
        let aux_at = qual_at + l_seq;
        let mut r = Record {
            data,
            tid,
            pos,
            mapq,
            flag,
            mate_tid,
            mate_pos,
            tlen,
            name: name_at..name_at + name_len,
            cigar: cigar_at..seq_at,
            seq: seq_at..qual_at,
            l_seq,
            qual: qual_at..aux_at,
            aux: aux_at..data.len(),
            cg: None,
            end: 0,
        };
        let from_tag = r.cigar_from_tag()?;
        let (rlen, qlen) = cigar_lens(&data[r.cigar.clone()]);
        if n_cigar > 0 || from_tag {
            let check = if from_tag {
                l_seq > 0
            } else {
                l_seq > 0 && flag & flag::UNMAPPED == 0
            };
            if check && qlen != l_seq as i64 {
                return Err(invalid(format!(
                    "CIGAR and query sequence lengths differ ({qlen} and {l_seq})"
                )));
            }
        }
        let rlen = if flag & flag::UNMAPPED != 0 { 0 } else { rlen };
        r.end = pos as i64 + if rlen == 0 { 1 } else { rlen };
        Ok(r)
    }

    /// `bam_tag2cigar`: a CIGAR of `<l_seq>S<n>N` stands for the one in the
    /// `CG:B:I` tag when the record is placed.
    fn cigar_from_tag(&mut self) -> Result<bool> {
        let n_cigar = self.cigar.len() / 4;
        if n_cigar == 0 {
            return Ok(false);
        }
        let first = u32_at(self.data, self.cigar.start);
        if first != (op::SOFT_CLIP as u32 | (self.l_seq as u32) << 4) {
            return Ok(false);
        }
        if self.tid < 0 || self.pos < 0 {
            return Ok(false);
        }
        let aux = &self.data[self.aux.clone()];
        let Some(at) = aux_get(aux, *b"CG").map_err(|_| invalid("corrupted aux data"))? else {
            return Ok(false);
        };
        if aux[at] != b'B' || !(aux[at + 1] == b'I' || aux[at + 1] == b'i') {
            return Ok(false);
        }
        let len = u32_at(aux, at + 2) as usize;
        if len < n_cigar || len >= 1 << 29 {
            return Ok(false);
        }
        let values = self.aux.start + at + 6;
        self.cigar = values..values + 4 * len;
        self.cg = Some(self.aux.start + at - 2..values + 4 * len);
        Ok(true)
    }

    /// The bytes the record was read from.
    pub fn data(&self) -> &'a [u8] {
        self.data
    }

    /// The read name, without its NUL.
    pub fn name(&self) -> &'a [u8] {
        &self.data[self.name.clone()]
    }

    /// Number of CIGAR operations (those of a `CG` tag when it holds them).
    pub fn n_cigar(&self) -> usize {
        self.cigar.len() / 4
    }

    /// The CIGAR operations as `(op, len)`.
    pub fn cigar(&self) -> impl ExactSizeIterator<Item = (u8, u32)> + 'a {
        self.data[self.cigar.clone()].chunks_exact(4).map(|c| {
            let v = u32::from_le_bytes(c.try_into().unwrap());
            ((v & 0xf) as u8, v >> 4)
        })
    }

    /// The raw CIGAR bytes (little-endian `len << 4 | op` words).
    pub fn cigar_bytes(&self) -> &'a [u8] {
        &self.data[self.cigar.clone()]
    }

    /// Query length.
    pub fn l_seq(&self) -> usize {
        self.l_seq
    }

    /// Bases packed two per byte, high nibble first.
    pub fn packed_seq(&self) -> &'a [u8] {
        &self.data[self.seq.clone()]
    }

    /// Base qualities as stored.
    pub fn qual(&self) -> &'a [u8] {
        &self.data[self.qual.clone()]
    }

    /// The optional fields as stored, in up to two pieces: a `CG` tag that
    /// held the CIGAR is left out.
    pub fn aux(&self) -> (&'a [u8], &'a [u8]) {
        match &self.cg {
            None => (&self.data[self.aux.clone()], &[]),
            Some(cg) => (
                &self.data[self.aux.start..cg.start],
                &self.data[cg.end..self.aux.end],
            ),
        }
    }

    /// Bytes of the optional fields [`Record::aux`] returns.
    pub fn aux_len(&self) -> usize {
        let (a, b) = self.aux();
        a.len() + b.len()
    }

    /// `bam_endpos`.
    pub fn end(&self) -> i64 {
        self.end
    }

    pub fn is_reverse(&self) -> bool {
        self.flag & flag::REVERSE != 0
    }

    pub fn entry(&self) -> AlignmentEntry {
        AlignmentEntry {
            start: self.pos as i64,
            end: self.end,
            mate_start: self.mate_pos,
            mate_tid: self.mate_tid,
            tlen: self.tlen,
            flag: self.flag,
            mapq: self.mapq,
            _pad: 0,
        }
    }
}

/// Expands packed bases into ASCII, one per query position.
pub fn unpack_seq(packed: &[u8], l_seq: usize, out: &mut [u8]) {
    debug_assert_eq!(out.len(), l_seq);
    let (pairs, odd) = (l_seq / 2, l_seq % 2);
    let out = &mut out[..l_seq];
    for (dst, &b) in out[..2 * pairs].chunks_exact_mut(2).zip(&packed[..pairs]) {
        dst.copy_from_slice(&BASE_PAIRS[b as usize]);
    }
    if odd == 1 {
        out[l_seq - 1] = NT16[(packed[pairs] >> 4) as usize];
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn cigar_types() {
        assert!(consumes_reference(op::MATCH) && consumes_query(op::MATCH));
        assert!(consumes_reference(op::DEL) && !consumes_query(op::DEL));
        assert!(consumes_reference(op::REF_SKIP) && !consumes_query(op::REF_SKIP));
        assert!(!consumes_reference(op::INS) && consumes_query(op::INS));
        assert!(!consumes_reference(op::SOFT_CLIP) && consumes_query(op::SOFT_CLIP));
        assert!(!consumes_reference(op::HARD_CLIP) && !consumes_query(op::HARD_CLIP));
        assert!(consumes_reference(op::EQUAL) && consumes_reference(op::DIFF));
    }

    #[test]
    fn unpack() {
        let mut out = [0u8; 3];
        unpack_seq(&[0x12, 0x40], 3, &mut out);
        assert_eq!(&out, b"ACG");
    }

    #[test]
    fn every_packed_code_and_odd_tail_match_independent_nibbles() {
        let packed: Vec<u8> = (0..=255).collect();
        for len in 0..=512 {
            let mut out = vec![0; len];
            unpack_seq(&packed, len, &mut out);
            let expected: Vec<_> = (0..len)
                .map(|i| {
                    let shift = if i % 2 == 0 { 4 } else { 0 };
                    b"=ACMGRSVTWYHKDBN"[((packed[i / 2] >> shift) & 15) as usize]
                })
                .collect();
            assert_eq!(out, expected, "length {len}");
        }
    }

    #[test]
    fn aux_walk() {
        let aux = b"NMi\x05\x00\x00\x00XZZabc\x00CGBI\x01\x00\x00\x00\x10\x00\x00\x00";
        assert_eq!(aux_get(aux, *b"NM"), Ok(Some(2)));
        assert_eq!(aux_get(aux, *b"XZ"), Ok(Some(9)));
        assert_eq!(aux_get(aux, *b"CG"), Ok(Some(16)));
        assert_eq!(aux_get(aux, *b"XX"), Ok(None));
        assert_eq!(aux_get(b"NMq\x01", *b"XX"), Err(()));
    }
}
