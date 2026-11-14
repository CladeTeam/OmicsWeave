//! One variant record, decoded from VCF text or BCF into the same form, and
//! the derived values htslib computes from it (reference length, variant
//! types).

/// Genotype values use htslib's encoding: `(allele + 1) << 1 | phased`, with
/// 0/1 for a missing allele.
pub const GT_MISSING_VALUE: i32 = i32::MIN;
/// Pads a sample's genotype when it has fewer alleles than the widest.
pub const GT_VECTOR_END: i32 = i32::MIN + 1;

/// Variant type bits (`bcf_get_variant_types`).
pub const VCF_REF: i32 = 0;
pub const VCF_SNP: i32 = 1;
pub const VCF_MNP: i32 = 2;
pub const VCF_INDEL: i32 = 4;
pub const VCF_OTHER: i32 = 8;
pub const VCF_BND: i32 = 16;
pub const VCF_OVERLAP: i32 = 32;
const VCF_INS: i32 = 64;
const VCF_DEL: i32 = 128;
const ORIG_VAR_TYPES: i32 = VCF_SNP | VCF_MNP | VCF_INDEL | VCF_OTHER | VCF_BND | VCF_OVERLAP;

/// A decoded record. Readers reuse one across records.
#[derive(Debug, Clone, Default, PartialEq, Eq)]
pub struct Record {
    /// Contig id in the header.
    pub rid: u32,
    /// 0-based position of the first REF base (-1 for a telomere at 0).
    pub pos: i64,
    /// Reference span, from REF, INFO/END, INFO/SVLEN and FORMAT/LEN as
    /// htslib computes it.
    pub rlen: i64,
    /// INFO/END (1-based, inclusive) when given as one integer.
    pub end: Option<i64>,
    allele_bytes: Vec<u8>,
    allele_ends: Vec<usize>,
    /// Genotype values, `gt_width` per sample; empty without genotypes.
    gt: Vec<i32>,
    gt_width: usize,
    has_gt: bool,
    /// ID, QUAL, FILTER and all INFO and FORMAT values (see `fields`).
    pub(crate) data: crate::fields::RecordData,
}

/// Genotypes of a record: `width` values per sample.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Genotypes<'a> {
    pub values: &'a [i32],
    pub width: usize,
}

impl Genotypes<'_> {
    /// Values of sample `i`; a zero-width layout yields an empty slice.
    ///
    /// Panics if the index arithmetic overflows or the sample is out of bounds.
    pub fn sample(&self, i: usize) -> &[i32] {
        let start = i.checked_mul(self.width).expect("sample offset overflow");
        let end = start.checked_add(self.width).expect("sample end overflow");
        &self.values[start..end]
    }
}

/// Whether a genotype value is a missing allele (`.`).
pub fn gt_is_missing(v: i32) -> bool {
    (v >> 1) == 0
}

/// The allele index of a genotype value (-1 when missing).
pub fn gt_allele(v: i32) -> i32 {
    (v >> 1) - 1
}

pub fn gt_is_phased(v: i32) -> bool {
    v & 1 != 0
}

impl Record {
    pub(crate) fn clear(&mut self) {
        self.rid = 0;
        self.pos = 0;
        self.rlen = 0;
        self.end = None;
        self.allele_bytes.clear();
        self.allele_ends.clear();
        self.gt.clear();
        self.gt_width = 0;
        self.has_gt = false;
        self.data.clear();
    }

    pub(crate) fn push_allele(&mut self, allele: &[u8]) {
        self.allele_bytes.extend_from_slice(allele);
        self.allele_ends.push(self.allele_bytes.len());
    }

    /// Marks the record as having genotypes of `width` and lends their
    /// buffer next to the stored fields, to decode one into the other.
    pub(crate) fn gt_and_data(
        &mut self,
        width: usize,
    ) -> (&mut Vec<i32>, &crate::fields::RecordData) {
        self.has_gt = true;
        self.gt_width = width;
        (&mut self.gt, &self.data)
    }

    pub fn n_alleles(&self) -> usize {
        self.allele_ends.len()
    }

    pub fn allele(&self, i: usize) -> &[u8] {
        let start = if i == 0 { 0 } else { self.allele_ends[i - 1] };
        &self.allele_bytes[start..self.allele_ends[i]]
    }

    pub fn alleles(&self) -> impl ExactSizeIterator<Item = &[u8]> + '_ {
        (0..self.n_alleles()).map(|i| self.allele(i))
    }

    /// The record's genotypes, if the header has a GT field and the record
    /// gives one.
    pub fn genotypes(&self) -> Option<Genotypes<'_>> {
        self.has_gt.then_some(Genotypes {
            values: &self.gt,
            width: self.gt_width,
        })
    }

    /// The union of the alternate alleles' types, restricted to the bits
    /// `bcf_get_variant_types` reports.
    pub fn variant_types(&self) -> i32 {
        variant_types(self.alleles())
    }
}

/// The union of the types of `alleles[1..]` against `alleles[0]`, as
/// `bcf_get_variant_types`.
pub fn variant_types<'a>(mut alleles: impl Iterator<Item = &'a [u8]>) -> i32 {
    let Some(reference) = alleles.next() else {
        return VCF_REF;
    };
    alleles.fold(0, |acc, alt| acc | variant_type(reference, alt)) & ORIG_VAR_TYPES
}

fn upper(c: u8) -> u8 {
    c.to_ascii_uppercase()
}

/// The type of one alternate allele (`bcf_set_variant_type`).
pub fn variant_type(reference: &[u8], alt: &[u8]) -> i32 {
    variant(reference, alt).0
}

/// The type of one alternate allele and its length change (`var->n` of
/// `bcf_set_variant_type`): bases inserted (positive) or deleted
/// (negative), 1 for a SNP, the longer side's length for an MNP. `None`
/// where htslib leaves it unset (symbolic and breakend alleles). Strings
/// are read as C strings: positions past the end read as NUL.
pub fn variant(reference: &[u8], alt: &[u8]) -> (i32, Option<i64>) {
    let r = |i: usize| reference.get(i).copied().unwrap_or(0);
    let a = |i: usize| alt.get(i).copied().unwrap_or(0);
    if a(1) == 0 {
        match a(0) {
            b'*' => return (VCF_OVERLAP, Some(0)),
            b'.' => return (VCF_REF, Some(0)),
            _ => {}
        }
        if r(1) == 0 {
            if r(0) == a(0) || a(0) == b'X' {
                return (VCF_REF, Some(0));
            }
            return (VCF_SNP, Some(1));
        }
    }
    if a(0) == b'<' {
        if (a(1) == b'X' || a(1) == b'*') && a(2) == b'>' {
            return (VCF_REF, Some(0));
        }
        if &alt[1..] == b"NON_REF>" {
            return (VCF_REF, Some(0));
        }
        return (VCF_OTHER, None);
    }
    if matches!(a(0), b']' | b'[' | b'.') {
        return (VCF_BND, None);
    }
    let (rl, al) = (reference.len(), alt.len());
    let mut i = 0;
    while i < rl && i < al && upper(r(i)) == upper(a(i)) {
        i += 1;
    }
    let breakend = |c: u8| matches!(c, b']' | b'[' | b'.');
    let change = Some(al as i64 - rl as i64);
    if i < al && i == rl {
        if breakend(a(al - 1)) {
            return (VCF_BND, None);
        }
        return (VCF_INDEL | VCF_INS, change);
    }
    if i < rl && i == al {
        return (VCF_INDEL | VCF_DEL, change);
    }
    if i == rl && i == al {
        return (VCF_REF, Some(0));
    }
    // Both continue past the common prefix; trim the common suffix.
    let (mut re, mut ae) = (rl - 1, al - 1);
    if breakend(a(ae)) {
        return (VCF_BND, None);
    }
    while re > i && ae > i && upper(r(re)) == upper(a(ae)) {
        re -= 1;
        ae -= 1;
    }
    let (dr, da) = ((re - i) as i64, (ae - i) as i64);
    if ae == i {
        if re == i {
            return (VCF_SNP, Some(1));
        }
        let ty = if upper(r(re)) == upper(a(ae)) {
            VCF_INDEL | VCF_DEL
        } else {
            VCF_OTHER
        };
        return (ty, Some(-dr));
    }
    if re == i {
        let ty = if upper(r(re)) == upper(a(ae)) {
            VCF_INDEL | VCF_INS
        } else {
            VCF_OTHER
        };
        return (ty, Some(da));
    }
    let ty = if dr == da { VCF_MNP } else { VCF_OTHER };
    (ty, Some(if dr > da { -(dr + 1) } else { da + 1 }))
}

/// A symbolic allele whose SVLEN spans the reference: `<CNV>`, `<DEL>`,
/// `<DUP>`, `<INV>` and their `:`-subtypes.
pub(crate) fn svlen_on_ref(alt: &[u8]) -> bool {
    alt.len() >= 5
        && alt[0] == b'<'
        && (alt[4] == b'>' || alt[4] == b':')
        && matches!(&alt[..4], b"<CNV" | b"<DEL" | b"<DUP" | b"<INV")
        && alt[alt.len() - 1] == b'>'
}

/// The fields htslib derives the reference span from. Integer values are
/// `None` where absent, missing or not of an integer type.
#[derive(Debug, Default)]
pub(crate) struct SpanFields {
    /// INFO/END, when given as one integer.
    pub end: Option<i64>,
    /// INFO/SVLEN (one per ALT), when the record has an SVLEN-spanning
    /// allele and the key is in the header. `Err` for a non-integer type.
    pub svlen: Option<Result<Vec<Option<i64>>, ()>>,
    /// Maximum FORMAT/LEN over all samples, when the record has a gVCF allele and
    /// the key is in the header. `Err` for a non-integer type.
    pub len: Option<Result<i64, ()>>,
}

/// The reference span (`get_rlen`).
pub(crate) fn rlen(pos: i64, record: &Record, f: &SpanFields) -> omicsweave_core::Result<i64> {
    let len_ref = if record.n_alleles() > 0 {
        record.allele(0).len() as i64
    } else {
        0
    };
    let pos = i128::from(pos);
    let end = i128::from(f.end.unwrap_or(0));
    // htslib marks SVLEN alleles in a bit array, then tests value i against
    // byte i/8 but bit (i+1)%8; the same test is kept.
    let mut len = 0i128;
    if let Some(Ok(values)) = &f.svlen {
        for (i, v) in values.iter().enumerate() {
            if i + 1 >= record.n_alleles() {
                break;
            }
            let marked_allele = (i & !7) | ((i + 1) & 7);
            if marked_allele == 0 || !svlen_on_ref(record.allele(marked_allele)) {
                continue;
            }
            len = len.max(i128::from(v.unwrap_or(0)).abs());
        }
    }
    if len == 0 && end != 0 {
        len = if end > pos { end - pos - 1 } else { 0 };
    }
    let end_svlen = pos + len + 1;

    let mut len = 0i128;
    if let Some(Ok(value)) = f.len {
        len = i128::from(value).max(0);
    }
    if len == 0 && end != 0 {
        len = if end > pos { end - pos } else { 0 };
    }
    let end_fmtlen = pos + len;

    let hpos = end.max(end_svlen).max(end_fmtlen);
    i64::try_from((hpos - pos).max(i128::from(len_ref)))
        .map_err(|_| omicsweave_core::Error::invalid("variant span", "reference span too large"))
}

/// Whether the record has a gVCF allele (`<*>`, `<NON_REF>`), for which
/// FORMAT/LEN extends the span.
pub(crate) fn has_gvcf_allele(record: &Record) -> bool {
    (1..record.n_alleles()).any(|i| matches!(record.allele(i), b"<*>" | b"<NON_REF>"))
}

/// Whether the record has an SVLEN-spanning symbolic allele.
pub(crate) fn has_svlen_allele(record: &Record) -> bool {
    (1..record.n_alleles()).any(|i| svlen_on_ref(record.allele(i)))
}
