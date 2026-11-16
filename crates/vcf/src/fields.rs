//! The fields beyond the variant itself: ID, QUAL, FILTER and every INFO and
//! FORMAT value. Values stay in BCF's typed encoding (one byte buffer per
//! section and a span per field), so their widths, sentinels and padding
//! are exactly what htslib holds after reading the same record, from VCF
//! text or BCF. Typed getters follow `bcf_get_info_values` and
//! `bcf_get_format_values`.

use crate::header::{Header, LineKind, ValueType};
use crate::record::Record;

/// BCF value type codes (`BCF_BT_*`).
pub const BT_NULL: u8 = 0;
pub const BT_INT8: u8 = 1;
pub const BT_INT16: u8 = 2;
pub const BT_INT32: u8 = 3;
pub const BT_INT64: u8 = 4;
pub const BT_FLOAT: u8 = 5;
pub const BT_CHAR: u8 = 7;

/// Bit patterns of the float sentinels (`bcf_float_missing`,
/// `bcf_float_vector_end`).
pub const FLOAT_MISSING: u32 = 0x7F80_0001;
pub const FLOAT_VECTOR_END: u32 = 0x7F80_0002;
/// The 32-bit integer sentinels the getters widen narrower ones to.
pub const INT32_MISSING: i32 = i32::MIN;
pub const INT32_VECTOR_END: i32 = i32::MIN + 1;

/// Bytes per value of a type; 0 for NULL and unknown types.
pub fn type_width(ty: u8) -> usize {
    match ty {
        BT_INT8 | BT_CHAR => 1,
        BT_INT16 => 2,
        BT_INT32 | BT_FLOAT => 4,
        BT_INT64 => 8,
        _ => 0,
    }
}

/// The stored integer of a width, sign-extended; a float's bits; a char
/// as a signed byte.
pub(crate) fn raw_at(ty: u8, b: &[u8]) -> i64 {
    match ty {
        BT_INT8 | BT_CHAR => b[0] as i8 as i64,
        BT_INT16 => i16::from_le_bytes([b[0], b[1]]) as i64,
        BT_INT32 => i32::from_le_bytes([b[0], b[1], b[2], b[3]]) as i64,
        BT_FLOAT => u32::from_le_bytes([b[0], b[1], b[2], b[3]]) as i64,
        BT_INT64 => i64::from_le_bytes(b[..8].try_into().unwrap()),
        _ => 0,
    }
}

fn missing_of(ty: u8) -> i64 {
    match ty {
        BT_INT8 => i8::MIN as i64,
        BT_INT16 => i16::MIN as i64,
        BT_INT32 => i32::MIN as i64,
        BT_INT64 => i64::MIN,
        BT_FLOAT => FLOAT_MISSING as i64,
        _ => i64::MAX,
    }
}

fn vector_end_of(ty: u8) -> i64 {
    match ty {
        BT_INT8 => i8::MIN as i64 + 1,
        BT_INT16 => i16::MIN as i64 + 1,
        BT_INT32 => i32::MIN as i64 + 1,
        BT_INT64 => i64::MIN + 1,
        BT_FLOAT => FLOAT_VECTOR_END as i64,
        _ => i64::MAX,
    }
}

/// Numeric getters validate the type once and decode complete chunks. The
/// public Typed view retains its checked indexing for caller-built slices.
fn numeric_values(ty: u8, bytes: &[u8]) -> impl Iterator<Item = i64> + '_ {
    let width = match ty {
        BT_INT8 => 1,
        BT_INT16 => 2,
        BT_INT32 | BT_FLOAT => 4,
        _ => unreachable!("numeric getter type validation"),
    };
    bytes.chunks_exact(width).map(move |chunk| {
        // SAFETY: chunks_exact supplies width initialized bytes. Each arm
        // reads exactly that width, unaligned, from the live shared slice;
        // no pointer escapes and the bytes remain immutable during the read.
        unsafe {
            match ty {
                BT_INT8 => *chunk.get_unchecked(0) as i8 as i64,
                BT_INT16 => i16::from_le(chunk.as_ptr().cast::<i16>().read_unaligned()) as i64,
                BT_INT32 => i32::from_le(chunk.as_ptr().cast::<i32>().read_unaligned()) as i64,
                BT_FLOAT => u32::from_le(chunk.as_ptr().cast::<u32>().read_unaligned()) as i64,
                _ => unreachable!("numeric getter type validation"),
            }
        }
    })
}

/// One decoded value.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum Value {
    Int(i64),
    Float(f32),
    Char(u8),
    Missing,
    VectorEnd,
}

/// Values of one BCF type, as stored.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Typed<'a> {
    /// `BT_*` code.
    pub ty: u8,
    pub bytes: &'a [u8],
}

impl<'a> Typed<'a> {
    pub fn len(&self) -> usize {
        match type_width(self.ty) {
            0 => 0,
            w => self.bytes.len() / w,
        }
    }

    pub fn is_empty(&self) -> bool {
        self.len() == 0
    }

    /// See [`raw_at`]: the stored integer, float bits or signed char.
    /// Panics if the index arithmetic overflows or the value is out of bounds.
    pub fn raw(&self, i: usize) -> i64 {
        let w = type_width(self.ty);
        let offset = i.checked_mul(w).expect("typed value offset overflow");
        raw_at(self.ty, &self.bytes[offset..])
    }

    pub fn value(&self, i: usize) -> Value {
        let raw = self.raw(i);
        match self.ty {
            BT_CHAR => Value::Char(raw as u8),
            _ if raw == missing_of(self.ty) => Value::Missing,
            _ if raw == vector_end_of(self.ty) => Value::VectorEnd,
            BT_FLOAT => Value::Float(f32::from_bits(raw as u32)),
            _ => Value::Int(raw),
        }
    }
}

/// An INFO field of a record.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct InfoField<'a> {
    /// Dictionary id.
    pub key: u32,
    pub values: Typed<'a>,
}

impl InfoField<'_> {
    /// htslib's `info->v1.i`: the single value of a one-value field,
    /// otherwise 0. Like the C union, a float's bits read as an integer.
    pub fn v1_int(&self) -> i64 {
        if self.values.len() != 1 {
            return 0;
        }
        self.values.raw(0)
    }

    /// htslib's `info->v1.f`: the float stored in the low bytes of `v1`.
    pub fn v1_float(&self) -> f32 {
        f32::from_bits(self.v1_int() as u32)
    }
}

/// A FORMAT field of a record: `n` values per sample.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct FormatField<'a> {
    /// Dictionary id.
    pub key: u32,
    pub n: usize,
    /// All samples' values, sample-major.
    pub values: Typed<'a>,
}

impl<'a> FormatField<'a> {
    /// Bytes per sample (`fmt->size`).
    /// Panics if the stored count cannot be represented as a byte length.
    pub fn size(&self) -> usize {
        self.n
            .checked_mul(type_width(self.values.ty))
            .expect("FORMAT size overflow")
    }

    pub fn sample(&self, i: usize) -> Typed<'a> {
        let size = self.size();
        let start = i.checked_mul(size).expect("FORMAT sample offset overflow");
        let end = start.checked_add(size).expect("FORMAT sample end overflow");
        Typed {
            ty: self.values.ty,
            bytes: &self.values.bytes[start..end],
        }
    }
}

/// Why a typed getter returns no values, with htslib's return codes.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum GetError {
    /// -1: the tag is not declared for this line kind.
    Undeclared,
    /// -2: declared with another type, or stored with an unexpected one.
    WrongType,
    /// -3: the record does not have the tag.
    Absent,
    /// A stored type htslib's FORMAT getter does not handle (it exits), or
    /// a header/record sample-count mismatch with no initialized C result.
    Unsupported,
}

impl GetError {
    pub fn code(self) -> i32 {
        match self {
            GetError::Undeclared => -1,
            GetError::WrongType => -2,
            GetError::Absent => -3,
            GetError::Unsupported => -100,
        }
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub(crate) struct Span {
    pub key: u32,
    pub ty: u8,
    /// Values (INFO) or values per sample (FORMAT).
    pub n: u32,
    pub start: usize,
    pub end: usize,
}

/// Storage of the fields; part of [`Record`].
#[derive(Debug, Clone, Default, PartialEq, Eq)]
pub(crate) struct RecordData {
    /// As htslib's `d.id` after unpacking: `.` when empty.
    pub id: Vec<u8>,
    pub qual: u32,
    pub filters: Vec<i32>,
    pub info: Vec<Span>,
    pub info_bytes: Vec<u8>,
    pub fmt: Vec<Span>,
    pub fmt_bytes: Vec<u8>,
    pub n_sample: u32,
}

impl RecordData {
    pub fn clear(&mut self) {
        self.id.clear();
        self.qual = 0;
        self.filters.clear();
        self.info.clear();
        self.info_bytes.clear();
        self.fmt.clear();
        self.fmt_bytes.clear();
        self.n_sample = 0;
    }

    /// Sets the ID from its stored bytes: a C string, `.` when empty.
    pub fn set_id(&mut self, bytes: &[u8]) {
        self.id.clear();
        if bytes.is_empty() {
            self.id.push(b'.');
        } else {
            let end = bytes.iter().position(|&b| b == 0).unwrap_or(bytes.len());
            self.id.extend_from_slice(&bytes[..end]);
        }
    }

    pub fn push_info(&mut self, key: u32, ty: u8, n: usize, bytes: &[u8]) {
        let start = self.info_bytes.len();
        self.info_bytes.extend_from_slice(bytes);
        self.info.push(Span {
            key,
            ty,
            n: n as u32,
            start,
            end: self.info_bytes.len(),
        });
    }

    pub fn push_format(&mut self, key: u32, ty: u8, n: usize, bytes: &[u8]) {
        let start = self.fmt_bytes.len();
        self.fmt_bytes.extend_from_slice(bytes);
        self.fmt.push(Span {
            key,
            ty,
            n: n as u32,
            start,
            end: self.fmt_bytes.len(),
        });
    }
}

fn info_view<'a>(d: &'a RecordData, s: &Span) -> InfoField<'a> {
    InfoField {
        key: s.key,
        values: Typed {
            ty: s.ty,
            bytes: &d.info_bytes[s.start..s.end],
        },
    }
}

fn format_view<'a>(d: &'a RecordData, s: &Span) -> FormatField<'a> {
    FormatField {
        key: s.key,
        n: s.n as usize,
        values: Typed {
            ty: s.ty,
            bytes: &d.fmt_bytes[s.start..s.end],
        },
    }
}

impl Record {
    /// The ID column as htslib unpacks it (`.` when empty).
    pub fn id(&self) -> &[u8] {
        &self.data.id
    }

    /// QUAL; `None` when missing.
    pub fn qual(&self) -> Option<f32> {
        (self.data.qual != FLOAT_MISSING).then(|| f32::from_bits(self.data.qual))
    }

    /// The bits of QUAL as stored, the missing sentinel included.
    pub fn qual_bits(&self) -> u32 {
        self.data.qual
    }

    /// FILTER as dictionary ids; empty for `.`.
    pub fn filters(&self) -> &[i32] {
        &self.data.filters
    }

    /// htslib's `n_sample`: the header's sample count when the record has
    /// a FORMAT column (or BCF sample data), otherwise 0.
    pub fn n_sample(&self) -> usize {
        self.data.n_sample as usize
    }

    pub fn info_fields(&self) -> impl ExactSizeIterator<Item = InfoField<'_>> + '_ {
        self.data.info.iter().map(|s| info_view(&self.data, s))
    }

    /// The first INFO field with this dictionary id.
    pub fn info_field(&self, key: u32) -> Option<InfoField<'_>> {
        self.data
            .info
            .iter()
            .find(|s| s.key == key)
            .map(|s| info_view(&self.data, s))
    }

    pub fn format_fields(&self) -> impl ExactSizeIterator<Item = FormatField<'_>> + '_ {
        self.data.fmt.iter().map(|s| format_view(&self.data, s))
    }

    /// The first FORMAT field with this dictionary id.
    pub fn format_field(&self, key: u32) -> Option<FormatField<'_>> {
        self.data
            .fmt
            .iter()
            .find(|s| s.key == key)
            .map(|s| format_view(&self.data, s))
    }

    fn info_lookup(
        &self,
        h: &Header,
        tag: &str,
        ty: ValueType,
    ) -> Result<Option<InfoField<'_>>, GetError> {
        let id = h
            .key_id(tag)
            .filter(|&id| h.declared_as(id, LineKind::Info))
            .ok_or(GetError::Undeclared)?;
        if h.key(id).and_then(|k| k.info) != Some(ty) {
            return Err(GetError::WrongType);
        }
        Ok(self.info_field(id))
    }

    /// `bcf_get_info_flag`: whether the record has the flag.
    pub fn info_flag(&self, h: &Header, tag: &str) -> Result<bool, GetError> {
        Ok(self.info_lookup(h, tag, ValueType::Flag)?.is_some())
    }

    /// `bcf_get_info_string`: the raw bytes (all values, comma-separated).
    pub fn info_string(&self, h: &Header, tag: &str, out: &mut Vec<u8>) -> Result<usize, GetError> {
        let f = self
            .info_lookup(h, tag, ValueType::String)?
            .ok_or(GetError::Absent)?;
        // htslib copies `len` bytes whatever the stored width.
        let n = f.values.len();
        out.clear();
        out.extend_from_slice(&f.values.bytes[..n]);
        Ok(n)
    }

    /// The values of an INFO field as htslib's numeric getters write them:
    /// 32-bit patterns up to the first vector end, integers with the 32-bit
    /// sentinels and floats as stored, whichever type was asked for.
    fn info_numeric(
        &self,
        h: &Header,
        tag: &str,
        ty: ValueType,
        mut emit: impl FnMut(u32),
    ) -> Result<usize, GetError> {
        let f = self.info_lookup(h, tag, ty)?.ok_or(GetError::Absent)?;
        if !matches!(f.values.ty, BT_INT8 | BT_INT16 | BT_INT32 | BT_FLOAT) {
            return Err(GetError::WrongType);
        }
        let missing = missing_of(f.values.ty);
        let end = vector_end_of(f.values.ty);
        let mut n = 0;
        for raw in numeric_values(f.values.ty, f.values.bytes) {
            match raw {
                raw if raw == end => break,
                raw if raw == missing && f.values.ty != BT_FLOAT => emit(INT32_MISSING as u32),
                raw => emit(raw as u32),
            }
            n += 1;
        }
        Ok(n)
    }

    /// `bcf_get_info_int32`: values up to the first vector end, missing
    /// ones as [`INT32_MISSING`].
    pub fn info_int32(&self, h: &Header, tag: &str, out: &mut Vec<i32>) -> Result<usize, GetError> {
        out.clear();
        self.info_numeric(h, tag, ValueType::Integer, |b| out.push(b as i32))
    }

    /// `bcf_get_info_float`: values up to the first vector end, missing
    /// ones as the missing bit pattern.
    pub fn info_float(&self, h: &Header, tag: &str, out: &mut Vec<f32>) -> Result<usize, GetError> {
        out.clear();
        self.info_numeric(h, tag, ValueType::Float, |b| out.push(f32::from_bits(b)))
    }

    fn format_lookup(
        &self,
        h: &Header,
        tag: &str,
        ty: ValueType,
    ) -> Result<FormatField<'_>, GetError> {
        let id = h
            .key_id(tag)
            .filter(|&id| h.declared_as(id, LineKind::Format))
            .ok_or(GetError::Undeclared)?;
        let declared = h.key(id).and_then(|k| k.format);
        // GT is declared a String but stored as integers.
        let expected = if tag == "GT" { ValueType::String } else { ty };
        if declared != Some(expected) {
            return Err(GetError::WrongType);
        }
        let field = self.format_field(id).ok_or(GetError::Absent)?;
        if h.samples.len() != self.n_sample() {
            return Err(GetError::Unsupported);
        }
        Ok(field)
    }

    /// `bcf_get_format_char` (`BCF_HT_STR`): the raw bytes of all samples,
    /// `n` per sample, whatever the stored type.
    pub fn format_char(&self, h: &Header, tag: &str, out: &mut Vec<u8>) -> Result<usize, GetError> {
        let f = self.format_lookup(h, tag, ValueType::String)?;
        let n = f.n * h.samples.len();
        out.clear();
        out.extend_from_slice(&f.values.bytes[..n.min(f.values.bytes.len())]);
        out.resize(n, 0);
        Ok(n)
    }

    /// The values of a FORMAT field as htslib's numeric getters write them:
    /// `n` 32-bit patterns per sample, integers with the 32-bit sentinels
    /// and floats as stored, whichever type was asked for.
    fn format_numeric(
        &self,
        h: &Header,
        tag: &str,
        ty: ValueType,
        mut emit: impl FnMut(u32),
    ) -> Result<usize, GetError> {
        let f = self.format_lookup(h, tag, ty)?;
        if !matches!(f.values.ty, BT_INT8 | BT_INT16 | BT_INT32 | BT_FLOAT) {
            return Err(GetError::Unsupported);
        }
        let float = f.values.ty == BT_FLOAT;
        let (missing, end) = if float {
            (FLOAT_MISSING, FLOAT_VECTOR_END)
        } else {
            (INT32_MISSING as u32, INT32_VECTOR_END as u32)
        };
        let nsmpl = self.n_sample();
        let stored_missing = missing_of(f.values.ty);
        let stored_end = vector_end_of(f.values.ty);
        for s in 0..nsmpl {
            let row = f.sample(s);
            let mut j = 0;
            for raw in numeric_values(row.ty, row.bytes) {
                j += 1;
                match raw {
                    raw if raw == stored_end => {
                        emit(end);
                        break;
                    }
                    raw if raw == stored_missing => emit(missing),
                    raw => emit(raw as u32),
                }
            }
            for _ in j..f.n {
                emit(end);
            }
        }
        Ok(nsmpl * f.n)
    }

    /// `bcf_get_format_int32`: `n` values per sample.
    pub fn format_int32(
        &self,
        h: &Header,
        tag: &str,
        out: &mut Vec<i32>,
    ) -> Result<usize, GetError> {
        out.clear();
        self.format_numeric(h, tag, ValueType::Integer, |b| out.push(b as i32))
    }

    /// `bcf_get_format_float`: `n` values per sample.
    pub fn format_float(
        &self,
        h: &Header,
        tag: &str,
        out: &mut Vec<f32>,
    ) -> Result<usize, GetError> {
        out.clear();
        self.format_numeric(h, tag, ValueType::Float, |b| out.push(f32::from_bits(b)))
    }

    /// `bcf_get_genotypes`.
    pub fn genotype_values(&self, h: &Header, out: &mut Vec<i32>) -> Result<usize, GetError> {
        self.format_int32(h, "GT", out)
    }
}

#[cfg(test)]
mod numeric_tests {
    use super::*;

    #[test]
    fn complete_unaligned_chunks_match_checked_public_decoding() {
        for ty in [BT_INT8, BT_INT16, BT_INT32, BT_FLOAT] {
            let mut storage = vec![99];
            for bits in [
                0,
                1,
                127,
                128,
                255,
                32767,
                32768,
                0x7fff_ffff,
                0x8000_0000,
                0x8000_0001,
                FLOAT_MISSING,
                FLOAT_VECTOR_END,
                0xffff_ffff,
            ] {
                storage.extend_from_slice(&bits.to_le_bytes()[..type_width(ty)]);
            }
            let bytes = &storage[1..];
            let view = Typed { ty, bytes };
            assert_eq!(
                numeric_values(ty, bytes).collect::<Vec<_>>(),
                (0..view.len()).map(|i| view.raw(i)).collect::<Vec<_>>()
            );
        }
    }
}
