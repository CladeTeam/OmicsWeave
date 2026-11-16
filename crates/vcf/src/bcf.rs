//! BCF (2.1/2.2) records: typed binary values, decoded as htslib reads them.

use omicsweave_core::{Error, Result};

use crate::fields::{FormatField, Typed};
use crate::header::{Header, VCF44, ValueType};
use crate::record::{
    GT_MISSING_VALUE, GT_VECTOR_END, Record, SpanFields, has_gvcf_allele, has_svlen_allele, rlen,
};

fn bad(detail: impl Into<String>) -> Error {
    Error::invalid("bcf record", detail)
}

const NULL: u8 = 0;
const INT8: u8 = 1;
const INT16: u8 = 2;
const INT32: u8 = 3;
const INT64: u8 = 4;
const FLOAT: u8 = 5;
const CHAR: u8 = 7;

pub(crate) fn width(ty: u8) -> Result<usize> {
    match ty {
        NULL => Ok(0),
        INT8 | CHAR => Ok(1),
        INT16 => Ok(2),
        INT32 | FLOAT => Ok(4),
        INT64 => Ok(8),
        t => Err(bad(format!("unknown value type {t}"))),
    }
}

/// An integer read with its width's missing and end-of-vector values
/// mapped to htslib's 32-bit ones.
enum Int {
    Value(i64),
    Missing,
    VectorEnd,
}

fn int_at(ty: u8, b: &[u8]) -> Int {
    match ty {
        INT8 => match b[0] as i8 {
            i8::MIN => Int::Missing,
            -127 => Int::VectorEnd,
            v => Int::Value(v as i64),
        },
        INT16 => match i16::from_le_bytes([b[0], b[1]]) {
            i16::MIN => Int::Missing,
            -32767 => Int::VectorEnd,
            v => Int::Value(v as i64),
        },
        INT32 => match i32::from_le_bytes(b[..4].try_into().unwrap()) {
            i32::MIN => Int::Missing,
            v if v == i32::MIN + 1 => Int::VectorEnd,
            v => Int::Value(v as i64),
        },
        INT64 => match i64::from_le_bytes(b[..8].try_into().unwrap()) {
            i64::MIN => Int::Missing,
            v if v == i64::MIN + 1 => Int::VectorEnd,
            v => Int::Value(v),
        },
        _ => Int::Missing,
    }
}

/// Reads typed values out of a byte slice.
pub(crate) struct Cursor<'a> {
    pub buf: &'a [u8],
    pub at: usize,
}

impl<'a> Cursor<'a> {
    fn take(&mut self, n: usize) -> Result<&'a [u8]> {
        let end = self
            .at
            .checked_add(n)
            .filter(|&e| e <= self.buf.len())
            .ok_or_else(|| bad("truncated"))?;
        let out = &self.buf[self.at..end];
        self.at = end;
        Ok(out)
    }

    /// A type descriptor: (type, count).
    pub(crate) fn descriptor(&mut self) -> Result<(u8, usize)> {
        let d = self.take(1)?[0];
        let (ty, n) = (d & 0xf, (d >> 4) as usize);
        if n != 15 {
            return Ok((ty, n));
        }
        // The extended length is one scalar, never another vector descriptor.
        let d = self.take(1)?[0];
        let length_ty = d & 0xf;
        if d >> 4 != 1 || !matches!(length_ty, INT8 | INT16 | INT32 | INT64) {
            return Err(bad("expected a scalar vector length"));
        }
        let n = match int_at(length_ty, self.take(width(length_ty)?)?) {
            Int::Value(n) => n,
            _ => return Err(bad("missing vector length")),
        };
        usize::try_from(n)
            .map(|n| (ty, n))
            .map_err(|_| bad("negative vector length"))
    }

    /// A single typed integer.
    pub(crate) fn int(&mut self) -> Result<i64> {
        let (ty, n) = self.descriptor()?;
        if ty == NULL && n != 0 {
            return Err(bad("nonempty null vector"));
        }
        if n != 1 || !matches!(ty, INT8 | INT16 | INT32 | INT64) {
            return Err(bad("expected a typed integer"));
        }
        let w = width(ty)?;
        match int_at(ty, self.take(w)?) {
            Int::Value(v) => Ok(v),
            _ => Err(bad("missing typed integer")),
        }
    }

    /// A typed vector: (type, count, bytes).
    pub fn vector(&mut self) -> Result<(u8, usize, &'a [u8])> {
        let (ty, n) = self.descriptor()?;
        if ty == NULL && n != 0 {
            return Err(bad("nonempty null vector"));
        }
        let bytes = self.take(
            n.checked_mul(width(ty)?)
                .ok_or_else(|| bad("vector size"))?,
        )?;
        Ok((ty, n, bytes))
    }
}

/// The 24 fixed bytes after the two lengths.
pub(crate) struct Fixed {
    pub rid: i32,
    pub pos: u32,
    pub rlen: i32,
    pub qual: u32,
    pub n_info: u16,
    pub n_allele: u16,
    pub n_sample: u32,
    pub n_fmt: u8,
}

impl Fixed {
    pub fn parse(b: &[u8; 24]) -> Fixed {
        let u32_at = |i: usize| u32::from_le_bytes(b[i..i + 4].try_into().unwrap());
        Fixed {
            rid: u32_at(0) as i32,
            pos: u32_at(4),
            rlen: u32_at(8) as i32,
            qual: u32_at(12),
            n_info: u16::from_le_bytes([b[16], b[17]]),
            n_allele: u16::from_le_bytes([b[18], b[19]]),
            n_sample: u32_at(20) & 0xff_ffff,
            n_fmt: b[23],
        }
    }
}

/// Decodes a record's shared and per-sample parts into `rec`.
pub(crate) fn parse_record(
    h: &Header,
    fixed: &Fixed,
    shared: &[u8],
    indiv: &[u8],
    rec: &mut Record,
) -> Result<()> {
    parse_record_fields::<true>(h, fixed, shared, indiv, rec).map(|_| ())
}

/// Matrix queries borrow the validated GT encoding and materialize only the
/// selected rows. Complete records own every field and normalize GT phasing.
pub(crate) fn parse_record_fields<'a, const UNPACK_FIELDS: bool>(
    h: &Header,
    fixed: &Fixed,
    shared: &[u8],
    indiv: &'a [u8],
    rec: &mut Record,
) -> Result<Option<FormatField<'a>>> {
    rec.clear();
    if fixed.rid < 0
        || fixed.rid as usize >= h.contigs.len()
        || h.contigs[fixed.rid as usize].name.is_empty()
    {
        return Err(bad(format!("contig id {} not in the header", fixed.rid)));
    }
    rec.rid = fixed.rid as u32;
    rec.pos = if fixed.pos == u32::MAX {
        -1
    } else {
        fixed.pos as i64
    };
    rec.rlen = fixed.rlen as i64;
    rec.data.qual = fixed.qual;
    rec.data.n_sample = fixed.n_sample;

    let mut c = Cursor { buf: shared, at: 0 };
    let (ty, _, id) = c.vector()?; // ID
    if ty != CHAR {
        return Err(bad("ID is not a string"));
    }
    rec.data.set_id(id);
    if fixed.n_allele == 0 {
        return Err(bad("missing REF allele"));
    }
    let repair_span = fixed.rlen < 0;
    // get_rlen inspects the encoded strings before unpacking, including any
    // NUL padding. Retain those bytes only for this exceptional repair path.
    let mut span_alleles = repair_span.then(Record::default);
    for _ in 0..fixed.n_allele {
        let (ty, _, bytes) = c.vector()?;
        if ty != CHAR {
            return Err(bad("allele is not a string"));
        }
        // Alleles are C strings: padding NULs end them.
        let end = bytes.iter().position(|&b| b == 0).unwrap_or(bytes.len());
        if let Some(raw) = &mut span_alleles {
            raw.push_allele(bytes);
        }
        rec.push_allele(&bytes[..end]);
    }
    let (ty, n, bytes) = c.vector()?; // FILTER
    if n != 0 {
        if !matches!(ty, INT8 | INT16 | INT32 | INT64) {
            return Err(bad("FILTER is not an integer vector"));
        }
        for value in bytes.chunks_exact(width(ty)?) {
            match int_at(ty, value) {
                Int::Value(id) => {
                    let key = checked_key(h, id)?;
                    if UNPACK_FIELDS {
                        rec.data.filters.push(key as i32);
                    }
                }
                _ => return Err(bad("missing FILTER id")),
            }
        }
    }
    let mut spans = SpanFields::default();
    let end_id = h.key_id("END");
    let svlen_id = repair_span
        .then(|| has_svlen_allele(span_alleles.as_ref().unwrap()))
        .filter(|v| *v)
        .and_then(|_| h.key_id("SVLEN"));
    let len_id = repair_span
        .then(|| has_gvcf_allele(span_alleles.as_ref().unwrap()))
        .filter(|v| *v)
        .and_then(|_| h.key_id("LEN"));
    let mut end_seen = false;
    for _ in 0..fixed.n_info {
        let key = checked_key(h, c.int()?)?;
        let (ty, n, bytes) = c.vector()?;
        if UNPACK_FIELDS {
            rec.data.push_info(key, ty, n, bytes);
        }
        if Some(key) == end_id && !end_seen {
            end_seen = true;
            if n == 1 && matches!(ty, INT8 | INT16 | INT32 | INT64) {
                rec.end = span_int(ty, bytes);
                spans.end = rec.end;
            }
        } else if repair_span && Some(key) == svlen_id && spans.svlen.is_none() {
            spans.svlen = Some(span_values(ty, bytes));
        }
    }
    if c.at != shared.len() {
        return Err(bad("trailing bytes in shared block"));
    }

    let n_fmt = if indiv.is_empty() || fixed.n_sample == 0 {
        0
    } else {
        fixed.n_fmt
    };
    let n_sample = h.samples.len();
    if n_fmt > 0 && fixed.n_sample as usize != n_sample {
        return Err(bad(format!(
            "{} samples, the header has {n_sample}",
            fixed.n_sample
        )));
    }
    let gt_name = h.key_id("GT");
    let gt_id = gt_name.filter(|&id| {
        h.key(id)
            .is_some_and(|key| key.format == Some(ValueType::String))
    });
    let mut borrowed_gt = None;
    let mut c = Cursor { buf: indiv, at: 0 };
    for _ in 0..n_fmt {
        let key = checked_key(h, c.int()?)?;
        let (ty, n) = c.descriptor()?;
        if ty == NULL && n != 0 {
            return Err(bad("nonempty null FORMAT vector"));
        }
        let w = width(ty)?;
        let size = n
            .checked_mul(w)
            .and_then(|s| s.checked_mul(n_sample))
            .filter(|&s| s <= i32::MAX as usize)
            .ok_or_else(|| bad("FORMAT size"))?;
        let bytes = c.take(size)?;
        if Some(key) == len_id && spans.len.is_none() {
            spans.len = Some(if matches!(ty, INT8 | INT16 | INT32 | INT64) {
                Ok(bytes
                    .chunks_exact(w)
                    .filter_map(|b| span_int(ty, b))
                    .max()
                    .unwrap_or(0))
            } else {
                Err(())
            });
        }
        if UNPACK_FIELDS {
            let start = rec.data.fmt_bytes.len();
            rec.data.push_format(key, ty, n, bytes);
            if Some(key) == gt_name && h.version < VCF44 && n > 0 {
                update_phasing(&mut rec.data.fmt_bytes[start..], n, w.max(1));
            }
        } else if Some(key) == gt_id && borrowed_gt.is_none() {
            if !matches!(ty, INT8 | INT16 | INT32) {
                return Err(bad("GT is not an integer vector"));
            }
            borrowed_gt = Some(FormatField {
                key,
                n,
                values: Typed { ty, bytes },
            });
        }
    }
    if let Some(gt) = gt_id.and_then(|id| rec.data.fmt.iter().find(|s| s.key == id).copied()) {
        if !matches!(gt.ty, INT8 | INT16 | INT32) {
            return Err(bad("GT is not an integer vector"));
        }
        if UNPACK_FIELDS {
            let (out, data) = rec.gt_and_data(gt.n as usize, gt.ty);
            let bytes = &data.fmt_bytes[gt.start..gt.end];
            decode_gt(gt.ty, gt.n as usize, bytes, n_sample, out);
        }
    }
    if n_fmt > 0 && c.at != indiv.len() {
        return Err(bad("trailing bytes in FORMAT block"));
    }
    if repair_span {
        rec.rlen = rlen(rec.pos, span_alleles.as_ref().unwrap(), &spans)?;
    }
    Ok(borrowed_gt)
}

/// htslib reads genotypes of files before VCF 4.4 with the first allele's
/// phase set as 4.4 defines it: a haploid allele is phased unless missing,
/// otherwise the first is phased if all the others are. It works on the
/// low byte of each stored value, end-of-vector values included, of every
/// FORMAT field keyed GT (`updatephasing`).
fn update_phasing(gt: &mut [u8], ploidy: usize, width: usize) {
    for row in gt.chunks_exact_mut(ploidy * width) {
        let phase = match ploidy {
            1 => u8::from(row[0] != 0),
            _ => (1..ploidy).fold(1, |all, k| all & row[k * width]) & 1,
        };
        row[0] |= phase;
    }
}

/// Genotype values of a stored GT field, widened to htslib's 32-bit
/// encoding; values after a sample's end-of-vector are end-of-vector.
pub(crate) fn decode_gt(ty: u8, n: usize, bytes: &[u8], n_sample: usize, out: &mut Vec<i32>) {
    let w = width(ty).unwrap_or(0);
    out.clear();
    out.reserve(n * n_sample);
    for s in 0..n_sample {
        let row = &bytes[s * n * w..(s + 1) * n * w];
        let mut ended = false;
        for j in 0..n {
            out.push(if ended {
                GT_VECTOR_END
            } else {
                match int_at(ty, &row[j * w..]) {
                    Int::Value(v) => v as i32,
                    Int::Missing => GT_MISSING_VALUE,
                    Int::VectorEnd => {
                        ended = true;
                        GT_VECTOR_END
                    }
                }
            });
        }
    }
}

pub(crate) fn checked_key(h: &Header, key: i64) -> Result<u32> {
    let id = u32::try_from(key).map_err(|_| bad("negative dictionary id"))?;
    h.key(id)
        .ok_or_else(|| bad("dictionary id not in the header"))?;
    Ok(id)
}

pub(crate) fn span_int(ty: u8, bytes: &[u8]) -> Option<i64> {
    match int_at(ty, bytes) {
        Int::Value(v) => Some(v),
        Int::Missing => None,
        Int::VectorEnd => Some(match ty {
            INT8 => i8::MIN as i64 + 1,
            INT16 => i16::MIN as i64 + 1,
            INT32 => i32::MIN as i64 + 1,
            _ => i64::MIN + 1,
        }),
    }
}

pub(crate) fn span_values(ty: u8, bytes: &[u8]) -> std::result::Result<Vec<Option<i64>>, ()> {
    if !matches!(ty, INT8 | INT16 | INT32 | INT64) {
        return Err(());
    }
    Ok(bytes
        .chunks_exact(width(ty).map_err(|_| ())?)
        .map(|b| span_int(ty, b))
        .collect())
}

#[cfg(test)]
mod selected_tests {
    use super::*;
    use crate::{BcfEncoder, BcfScratch};

    #[test]
    fn selected_view_borrows_gt_and_preserves_repaired_end_svlen_and_len_spans() {
        let header = Header::from_vcf_text(
            b"##fileformat=VCFv4.5\n\
            ##contig=<ID=chr1,length=100000>\n\
            ##INFO=<ID=END,Number=1,Type=Integer>\n\
            ##INFO=<ID=SVLEN,Number=A,Type=Integer>\n\
            ##FORMAT=<ID=GT,Number=1,Type=String>\n\
            ##FORMAT=<ID=LEN,Number=1,Type=Integer>\n\
            #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ts0\ts1\n",
        )
        .unwrap();
        let encoder = BcfEncoder::new(&header);
        let mut scratch = BcfScratch::default();
        for (line, span) in [
            (
                b"chr1\t10\tv\tA\tC\t.\tPASS\tEND=25\tGT:LEN\t0|1:3\t1:7".as_slice(),
                16,
            ),
            (
                b"chr1\t10\tv\tA\t<DEL>\t.\tPASS\tEND=11;SVLEN=100\tGT:LEN\t0|1:3\t1:7",
                101,
            ),
            (
                b"chr1\t10\tv\tA\t<NON_REF>\t.\tPASS\tEND=11\tGT:LEN\t0|1:3\t1:70",
                70,
            ),
        ] {
            let mut bytes = Vec::new();
            encoder
                .encode_line(&mut scratch, line, 1, &mut bytes)
                .unwrap();
            let l_shared = u32::from_le_bytes(bytes[..4].try_into().unwrap()) as usize;
            let mut fixed = Fixed::parse(bytes[8..32].try_into().unwrap());
            fixed.rlen = -1;
            let shared = &bytes[32..8 + l_shared];
            let indiv = &bytes[8 + l_shared..];
            let (mut full, mut selected) = (Record::default(), Record::default());
            parse_record(&header, &fixed, shared, indiv, &mut full).unwrap();
            let gt = parse_record_fields::<false>(&header, &fixed, shared, indiv, &mut selected)
                .unwrap()
                .unwrap();
            assert_eq!(full.rlen, span);
            assert_eq!(
                (selected.pos, selected.rlen, selected.end),
                (full.pos, full.rlen, full.end)
            );
            assert_eq!(
                selected.alleles().collect::<Vec<_>>(),
                full.alleles().collect::<Vec<_>>()
            );
            assert_eq!(selected.id(), full.id());
            assert_eq!(gt, full.format_field(header.key_id("GT").unwrap()).unwrap());
            assert!(gt.values.bytes.as_ptr() >= indiv.as_ptr());
            assert!(gt.values.bytes.as_ptr_range().end <= indiv.as_ptr_range().end);
            assert_eq!(selected.data.info_bytes.capacity(), 0);
            assert_eq!(selected.data.fmt_bytes.capacity(), 0);
            assert_eq!(selected.data.info.capacity(), 0);
            assert_eq!(selected.data.fmt.capacity(), 0);
            assert!(selected.genotypes().is_none());
        }
    }
}
