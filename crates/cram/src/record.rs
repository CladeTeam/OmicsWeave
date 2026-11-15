//! A decoded record as the BAM record htslib builds from it
//! (`cram_to_bam`, `bam_set1`, `bam_tag2cigar`), and that record in BAM's
//! binary layout (`bam_write1`).

use omicsweave_core::{AlignmentEntry, Error, Result};

use crate::slice::{Decoded, MAX_NAME};

const UNMAPPED: u16 = 0x4;

/// `seq_nt16_table`: ASCII to 4-bit base codes.
const fn nt16_table() -> [u8; 256] {
    let mut t = [15u8; 256];
    let codes = b"=ACMGRSVTWYHKDBN";
    let mut i = 0;
    while i < 16 {
        t[codes[i] as usize] = i as u8;
        t[codes[i].to_ascii_lowercase() as usize] = i as u8;
        i += 1;
    }
    t
}

pub(crate) const NT16_CODE: [u8; 256] = nt16_table();
pub(crate) const NT16: &[u8; 16] = b"=ACMGRSVTWYHKDBN";

/// Reference and query lengths of a CIGAR.
pub(crate) fn cigar_lens(cigar: &[u32]) -> (i64, i64) {
    let (mut r, mut q) = (0i64, 0i64);
    for &c in cigar {
        let (op, len) = (c & 0xf, (c >> 4) as i64);
        // BAM_CIGAR_TYPE: bit 0 query, bit 1 reference.
        let t = (0x3C1A7u32 >> (op << 1)) & 3;
        if t & 1 != 0 {
            q += len;
        }
        if t & 2 != 0 {
            r += len;
        }
    }
    (r, q)
}

/// `hts_reg2bin(beg, end, 14, 5)`.
pub(crate) fn reg2bin(beg: i64, end: i64) -> u16 {
    let end = end - 1;
    let mut s = 14;
    let mut t = ((1i64 << 15) - 1) / 7;
    let mut l = 5;
    while l > 0 {
        if beg >> s == end >> s {
            return (t + (beg >> s)) as u16;
        }
        l -= 1;
        s += 3;
        t -= 1 << (3 * l);
    }
    0
}

/// `aux_type2size`.
fn aux_size(t: u8) -> usize {
    match t {
        b'A' | b'c' | b'C' => 1,
        b's' | b'S' => 2,
        b'i' | b'I' | b'f' => 4,
        b'd' => 8,
        b'Z' | b'H' | b'B' => t as usize,
        _ => 0,
    }
}

/// `skip_aux`: the end of the value whose type is at `s`.
fn skip_aux(aux: &[u8], s: usize) -> Option<usize> {
    let end = aux.len();
    if s >= end {
        return Some(end);
    }
    let size = aux_size(aux[s]);
    let s = s + 1;
    match size {
        0 => None,
        90 | 72 => Some(memchr::memchr(0, &aux[s..]).map_or(end, |p| s + p + 1)),
        66 => {
            if end - s < 5 {
                return None;
            }
            let size = aux_size(aux[s]);
            let n = u32::from_le_bytes(aux[s + 1..s + 5].try_into().unwrap()) as usize;
            let s = s + 5;
            if size == 0 || (end - s) < size.checked_mul(n)? {
                return None;
            }
            Some(s + size * n)
        }
        _ => (end - s >= size).then_some(s + size),
    }
}

/// `bam_aux_get`: the position of `tag`'s type byte.
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

/// A record's BAM fields, located in its slice's buffers or built for it.
#[derive(Debug, Default)]
pub(crate) struct BamFields {
    pub name: Vec<u8>,
    pub cigar: Vec<u32>,
    pub aux: Vec<u8>,
    pub entry: AlignmentEntry,
    pub tid: i32,
    pub bin: u16,
    /// Query length (0 when the record stores no bases).
    pub l_seq: usize,
}

fn invalid(detail: impl Into<String>) -> Error {
    Error::invalid("CRAM record", detail)
}

/// What naming and read groups take from the file.
pub(crate) struct Naming<'a> {
    /// Prefix of generated read names (htslib: the file's base name).
    pub prefix: &'a [u8],
    pub read_groups: &'a [Vec<u8>],
}

impl Decoded {
    /// `cram_to_bam` for record `i`, into `f`. Names and aux fields are
    /// built when `name` / `aux` (they need those series decoded).
    pub(crate) fn bam_fields(
        &self,
        i: usize,
        naming: &Naming<'_>,
        name: bool,
        aux: bool,
        f: &mut BamFields,
    ) -> Result<()> {
        let cr = &self.recs[i];
        f.name.clear();
        if name {
            if cr.name_len > 0 {
                let at = cr.name as usize;
                f.name
                    .extend_from_slice(&self.names[at..at + cr.name_len as usize]);
            } else {
                let mate = (cr.mate_line >= 0 && (cr.mate_line as usize) < self.recs.len())
                    .then(|| &self.recs[cr.mate_line as usize])
                    .filter(|m| m.name_len > 0);
                match mate {
                    Some(m) => {
                        let at = m.name as usize;
                        f.name
                            .extend_from_slice(&self.names[at..at + m.name_len as usize]);
                    }
                    None => {
                        let keep = naming.prefix.len().min((MAX_NAME - 21) as usize);
                        f.name.extend_from_slice(&naming.prefix[..keep]);
                        f.name.push(b':');
                        let n = if cr.mate_line >= 0 && (cr.mate_line as usize) < i {
                            cr.mate_line as i64
                        } else {
                            i as i64
                        };
                        f.name.extend_from_slice(
                            (self.record_counter + n + 1).to_string().as_bytes(),
                        );
                    }
                }
            }
        }
        if cr.rg < -1 || cr.rg >= naming.read_groups.len() as i32 {
            return Err(invalid(format!("read group {}", cr.rg)));
        }
        f.cigar.clear();
        let at = cr.cigar as usize;
        f.cigar
            .extend_from_slice(&self.cigars[at..at + cr.ncigar as usize]);
        f.aux.clear();
        if aux {
            let at = cr.aux as usize;
            f.aux
                .extend_from_slice(&self.aux[at..at + cr.aux_size as usize]);
            if cr.rg >= 0 {
                f.aux.extend_from_slice(b"RGZ");
                f.aux.extend_from_slice(&naming.read_groups[cr.rg as usize]);
                f.aux.push(0);
            }
        }
        let flag = cr.flags as u16;
        let pos = cr.apos - 1;
        let l_seq = cr.len.max(0) as usize;
        if f.name.len() > MAX_NAME as usize {
            return Err(invalid("read name too long"));
        }
        let mapped = flag & UNMAPPED == 0;
        let (mut rlen, qlen) = if mapped { cigar_lens(&f.cigar) } else { (0, 0) };
        if mapped && l_seq > 0 && f.cigar.is_empty() {
            return Err(invalid("mapped read without a CIGAR"));
        }
        if mapped && l_seq > 0 && l_seq as i64 != qlen {
            return Err(invalid("CIGAR and query sequence lengths differ"));
        }
        if rlen == 0 {
            rlen = 1;
        }
        f.bin = reg2bin(pos, pos + rlen);
        if cr.has_cg && aux {
            tag2cigar(f, cr.ref_id, pos, l_seq)?;
            let (r, _) = if mapped { cigar_lens(&f.cigar) } else { (0, 0) };
            f.bin = reg2bin(pos, pos + r.max(1));
        }
        let (r, _) = if mapped { cigar_lens(&f.cigar) } else { (0, 0) };
        f.tid = cr.ref_id;
        f.l_seq = l_seq;
        f.entry = AlignmentEntry {
            start: pos,
            end: pos + r.max(1),
            mate_start: (cr.mate_pos - 1) as i32,
            mate_tid: cr.mate_ref_id,
            tlen: cr.tlen as i32,
            flag,
            mapq: cr.mqual as u8,
            _pad: 0,
        };
        Ok(())
    }

    /// Bases of record `i` as BAM holds them (through the 4-bit codes),
    /// appended to `out`.
    pub(crate) fn push_seq(&self, i: usize, l_seq: usize, out: &mut Vec<u8>) {
        let at = self.recs[i].seq as usize;
        out.extend(
            self.seqs[at..at + l_seq]
                .iter()
                .map(|&c| NT16[NT16_CODE[c as usize] as usize]),
        );
    }

    /// Qualities of record `i`, appended to `out`.
    pub(crate) fn push_qual(&self, i: usize, l_seq: usize, out: &mut Vec<u8>) {
        let at = self.recs[i].qual as usize;
        out.extend_from_slice(&self.quals[at..at + l_seq]);
    }

    /// Record `i` in BAM's binary layout, `block_size` first, appended to
    /// `out` (`bam_write1`).
    pub(crate) fn push_bam(&self, i: usize, f: &BamFields, out: &mut Vec<u8>) -> Result<()> {
        let cr = &self.recs[i];
        let e = &f.entry;
        let n_cigar = f.cigar.len();
        let long = n_cigar > 0xffff;
        let l_seq = f.l_seq;
        let start = out.len();
        out.extend_from_slice(&[0; 4]);
        if e.start > i32::MAX as i64 || cr.tlen < i32::MIN as i64 || cr.tlen > i32::MAX as i64 {
            return Err(invalid("positional data is too large for BAM"));
        }
        let put = |out: &mut Vec<u8>, v: u32| out.extend_from_slice(&v.to_le_bytes());
        put(out, f.tid as u32);
        put(out, e.start as u32);
        put(
            out,
            (f.bin as u32) << 16 | (e.mapq as u32) << 8 | (f.name.len() as u32 + 1),
        );
        put(
            out,
            (e.flag as u32) << 16 | if long { 2 } else { n_cigar as u32 },
        );
        put(out, l_seq as u32);
        put(out, e.mate_tid as u32);
        put(out, (cr.mate_pos - 1) as u32);
        put(out, cr.tlen as u32);
        out.extend_from_slice(&f.name);
        out.push(0);
        if long {
            let (rlen, _) = cigar_lens(&f.cigar);
            if rlen >= 1 << 28 {
                return Err(invalid("a CIGAR this long cannot be written to BAM"));
            }
            put(out, (l_seq as u32) << 4 | 4);
            put(out, (rlen as u32) << 4 | 3);
        } else {
            for &c in &f.cigar {
                put(out, c);
            }
        }
        let at = cr.seq as usize;
        let seq = &self.seqs[at..at + l_seq.min(self.seqs.len().saturating_sub(at))];
        let mut pairs = seq.chunks_exact(2);
        for p in &mut pairs {
            out.push(NT16_CODE[p[0] as usize] << 4 | NT16_CODE[p[1] as usize]);
        }
        if let [last] = pairs.remainder() {
            out.push(NT16_CODE[*last as usize] << 4);
        }
        if self.has_qual {
            self.push_qual(i, l_seq, out);
        } else {
            out.resize(out.len() + l_seq, 0xff);
        }
        out.extend_from_slice(&f.aux);
        if long {
            out.extend_from_slice(b"CGBI");
            put(out, n_cigar as u32);
            for &c in &f.cigar {
                put(out, c);
            }
        }
        let size = (out.len() - start - 4) as u32;
        out[start..start + 4].copy_from_slice(&size.to_le_bytes());
        Ok(())
    }
}

/// `bam_tag2cigar`: a CIGAR of `<l_seq>S<n>N` is replaced by a `CG:B:I`
/// tag's operations, and the tag leaves the aux fields.
fn tag2cigar(f: &mut BamFields, tid: i32, pos: i64, l_seq: usize) -> Result<()> {
    let placeholder = 4u32 | (l_seq as u32) << 4;
    if f.cigar.first() != Some(&placeholder) || tid < 0 || pos < 0 {
        return Ok(());
    }
    let Some(at) = aux_get(&f.aux, *b"CG").map_err(|_| invalid("corrupted aux data"))? else {
        return Ok(());
    };
    let aux = &f.aux;
    if aux[at] != b'B' || !(aux[at + 1] == b'I' || aux[at + 1] == b'i') {
        return Ok(());
    }
    let len = u32::from_le_bytes(aux[at + 2..at + 6].try_into().unwrap()) as usize;
    if len < f.cigar.len() || len >= 1 << 29 {
        return Ok(());
    }
    let values = at + 6;
    f.cigar = aux[values..values + 4 * len]
        .chunks_exact(4)
        .map(|c| u32::from_le_bytes(c.try_into().unwrap()))
        .collect();
    f.aux.drain(at - 2..values + 4 * len);
    if l_seq > 0 && !f.cigar.is_empty() && cigar_lens(&f.cigar).1 != l_seq as i64 {
        return Err(invalid("CIGAR and query sequence lengths differ"));
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn bins() {
        assert_eq!(reg2bin(-1, 0), 4680);
        assert_eq!(reg2bin(0, 1), 4681);
        assert_eq!(reg2bin(0, 1 << 14), 4681);
        assert_eq!(reg2bin(0, (1 << 14) + 1), 585);
        assert_eq!(reg2bin(0, 1 << 29), 0);
    }
}
