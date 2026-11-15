//! VCF text records, parsed as htslib's `vcf_parse` reads them: every
//! column is encoded into the typed form htslib builds (integer widths,
//! sentinels, string padding), so a VCF record and the BCF htslib writes
//! from it decode alike.

use omicsweave_core::{Error, Result};

use crate::bcf::{decode_gt, span_int, span_values};
use crate::cnum::{hts_str2dbl, hts_str2int, hts_str2uint, strtod};
use crate::fields::{
    BT_CHAR, BT_FLOAT, BT_INT8, BT_INT16, BT_INT32, BT_NULL, FLOAT_MISSING, FLOAT_VECTOR_END,
    INT32_MISSING, INT32_VECTOR_END,
};
use crate::header::{Header, VCF44, ValueType};
use crate::record::{Record, SpanFields, has_gvcf_allele, has_svlen_allele, rlen};

fn bad(detail: impl Into<String>) -> Error {
    Error::invalid("vcf record", detail)
}

/// `BCF_MIN_BT_INT32`: the smallest value outside the reserved ones.
const MIN_INT32: i64 = i32::MIN as i64 + 8;
/// htslib's limit on FORMAT keys per record (`n_fmt` is one byte).
const MAX_N_FMT: usize = 255;

const HT_INT: u32 = 1;
const HT_REAL: u32 = 2;
const HT_STR: u32 = 3;

/// Working buffers of the FORMAT parser, kept by the reader so parsing
/// does not allocate once warm.
#[derive(Debug, Default)]
pub(crate) struct Scratch {
    /// The sample columns, with sample-separating tabs turned into NULs.
    line: Vec<u8>,
    /// Per-key value arrays, as htslib's `h->mem`.
    mem: Vec<u8>,
    aux: Vec<FmtAux>,
    ints: Vec<i32>,
    enc: Vec<u8>,
    key_order: Vec<(u32, usize)>,
}

/// htslib's `fmt_aux_t`.
#[derive(Debug, Clone, Copy, Default)]
struct FmtAux {
    key: u32,
    is_gt: bool,
    /// `BCF_HT_*` of the declaration.
    htype: u32,
    max_m: u32,
    max_l: u32,
    max_g: u32,
    /// Bytes per sample; `None` for a dropped duplicate.
    size: Option<usize>,
    offset: usize,
}

/// `bcf_enc_int1`: the narrowest type for one value.
fn enc_int1(x: i32, out: &mut Vec<u8>) -> u8 {
    out.clear();
    if x == INT32_VECTOR_END {
        out.push(i8::MIN as u8 + 1);
        BT_INT8
    } else if x == INT32_MISSING {
        out.push(i8::MIN as u8);
        BT_INT8
    } else if (-120..=127).contains(&x) {
        out.push(x as i8 as u8);
        BT_INT8
    } else if (-32760..=32767).contains(&x) {
        out.extend_from_slice(&(x as i16).to_le_bytes());
        BT_INT16
    } else {
        out.extend_from_slice(&x.to_le_bytes());
        BT_INT32
    }
}

/// `bcf_enc_vint`: the narrowest type holding all values (sentinels
/// mapped to the type's own).
fn enc_vint(a: &[i32], out: &mut Vec<u8>) -> u8 {
    out.clear();
    match a.len() {
        0 => return BT_NULL,
        1 => return enc_int1(a[0], out),
        _ => {}
    }
    let max = a.iter().copied().max().unwrap();
    let min = a
        .iter()
        .copied()
        .filter(|&v| v > INT32_VECTOR_END)
        .min()
        .unwrap_or(i32::MAX);
    if max <= 127 && min >= -120 {
        for &v in a {
            out.push(match v {
                INT32_VECTOR_END => i8::MIN as u8 + 1,
                INT32_MISSING => i8::MIN as u8,
                v => v as i8 as u8,
            });
        }
        BT_INT8
    } else if max <= 32767 && min >= -32760 {
        for &v in a {
            let x = match v {
                INT32_VECTOR_END => i16::MIN + 1,
                INT32_MISSING => i16::MIN,
                v => v as i16,
            };
            out.extend_from_slice(&x.to_le_bytes());
        }
        BT_INT16
    } else {
        for &v in a {
            out.extend_from_slice(&v.to_le_bytes());
        }
        BT_INT32
    }
}

/// POS as `vcf_parse` reads it: `hts_str2uint` of 62 bits that must use
/// the whole column.
fn parse_pos(field: &[u8]) -> Result<i64> {
    let (v, used, overflow) = hts_str2uint(field, 62);
    if overflow {
        return Err(bad("position too large"));
    }
    if used != field.len() {
        return Err(bad(format!(
            "could not parse the position '{}'",
            String::from_utf8_lossy(field)
        )));
    }
    Ok(v as i64 - 1)
}

/// How a record's names resolve against the header: a [`Header`] adds
/// the ones it lacks as htslib does while reading; [`Declared`] rejects
/// them.
pub(crate) trait Dictionary {
    /// Parsing for BCF output, which htslib refuses to write otherwise:
    /// a FORMAT key may not repeat, a header with samples needs the
    /// FORMAT column, and genotypes are not decoded into the record.
    const STRICT: bool;
    fn header(&self) -> &Header;
    fn contig(&mut self, name: &str) -> Result<u32>;
    /// Any dictionary key names a FILTER, as htslib looks them up.
    fn filter(&mut self, name: &str) -> Result<u32>;
    fn info(&mut self, name: &str) -> Result<(u32, ValueType)>;
    fn format(&mut self, name: &str) -> Result<(u32, ValueType)>;
}

impl Dictionary for Header {
    const STRICT: bool = false;

    fn header(&self) -> &Header {
        self
    }

    fn contig(&mut self, name: &str) -> Result<u32> {
        Ok(self.add_contig(name))
    }

    fn filter(&mut self, name: &str) -> Result<u32> {
        Ok(self.add_undeclared_filter(name))
    }

    fn info(&mut self, name: &str) -> Result<(u32, ValueType)> {
        if let Some(id) = self.key_id(name)
            && let Some(ty) = self.key(id).and_then(|key| key.info)
        {
            return Ok((id, ty));
        }
        Ok((self.add_undeclared(name, false), ValueType::String))
    }

    fn format(&mut self, name: &str) -> Result<(u32, ValueType)> {
        if let Some(id) = self.key_id(name)
            && let Some(ty) = self.key(id).and_then(|key| key.format)
        {
            return Ok((id, ty));
        }
        Ok((self.add_undeclared(name, true), ValueType::String))
    }
}

/// A header that is complete: names it does not declare are errors.
pub(crate) struct Declared<'h>(pub &'h Header);

impl Dictionary for Declared<'_> {
    const STRICT: bool = true;

    fn header(&self) -> &Header {
        self.0
    }

    fn contig(&mut self, name: &str) -> Result<u32> {
        self.0
            .contig_id(name)
            .ok_or_else(|| bad(format!("contig '{name}' is not declared in the header")))
    }

    fn filter(&mut self, name: &str) -> Result<u32> {
        self.0
            .key_id(name)
            .ok_or_else(|| bad(format!("FILTER '{name}' is not declared in the header")))
    }

    fn info(&mut self, name: &str) -> Result<(u32, ValueType)> {
        let id = self.0.key_id(name);
        id.zip(id.and_then(|id| self.0.key(id)?.info))
            .ok_or_else(|| bad(format!("INFO '{name}' is not declared in the header")))
    }

    fn format(&mut self, name: &str) -> Result<(u32, ValueType)> {
        let id = self.0.key_id(name);
        id.zip(id.and_then(|id| self.0.key(id)?.format))
            .ok_or_else(|| bad(format!("FORMAT '{name}' is not declared in the header")))
    }
}

/// FILTER (`vcf_parse_filter`): names become dictionary ids; a name not
/// in the dictionary is added as a FILTER.
fn parse_filter<D: Dictionary>(d: &mut D, field: &[u8], rec: &mut Record) -> Result<()> {
    let field = field.strip_suffix(b";").unwrap_or(field);
    for name in field.split(|&c| c == b';') {
        let name = std::str::from_utf8(name).map_err(|_| bad("FILTER is not UTF-8"))?;
        let id = d.filter(name)?;
        rec.data.filters.push(id as i32);
    }
    Ok(())
}

/// INFO (`vcf_parse_info`): each key's values encoded by its declared
/// type; undeclared keys are added as `Number=1,Type=String`.
fn parse_info<D: Dictionary>(
    d: &mut D,
    info: &[u8],
    rec: &mut Record,
    sc: &mut Scratch,
) -> Result<()> {
    let info = info.strip_suffix(b";").unwrap_or(info);
    for item in info.split(|&c| c == b';') {
        let (key, value) = match item.iter().position(|&c| c == b'=') {
            Some(i) => (&item[..i], Some(&item[i + 1..])),
            None => (item, None),
        };
        if key.is_empty() {
            continue;
        }
        if rec.data.info.len() == u16::MAX as usize {
            return Err(bad("too many INFO entries"));
        }
        let name = std::str::from_utf8(key).map_err(|_| bad("INFO tag is not UTF-8"))?;
        let (id, ty) = d.info(name)?;
        let Some(value) = value else {
            rec.data.push_info(id, BT_NULL, 0, &[]);
            continue;
        };
        match ty {
            // A Flag with a value is kept as a string.
            ValueType::Flag | ValueType::String => {
                rec.data.push_info(id, BT_CHAR, value.len(), value);
            }
            ValueType::Integer => {
                sc.ints.clear();
                for v in value.split(|&c| c == b',') {
                    let (x, n, overflow) = hts_str2int(v, 64);
                    sc.ints.push(
                        if n == 0 || overflow || !(MIN_INT32..=i32::MAX as i64).contains(&x) {
                            INT32_MISSING
                        } else {
                            x as i32
                        },
                    );
                }
                let ty = if sc.ints.len() == 1 {
                    enc_int1(sc.ints[0], &mut sc.enc)
                } else {
                    enc_vint(&sc.ints, &mut sc.enc)
                };
                rec.data.push_info(id, ty, sc.ints.len(), &sc.enc);
            }
            ValueType::Float => {
                sc.enc.clear();
                let mut n = 0;
                for v in value.split(|&c| c == b',') {
                    let (x, used, failed) = hts_str2dbl(v);
                    let bits = if used == 0 || failed {
                        FLOAT_MISSING
                    } else {
                        (x as f32).to_bits()
                    };
                    sc.enc.extend_from_slice(&bits.to_le_bytes());
                    n += 1;
                }
                rec.data.push_info(id, BT_FLOAT, n, &sc.enc);
            }
        }
    }
    Ok(())
}

/// Parses one data line into `rec`, registering undeclared contigs and keys
/// in the header.
pub(crate) fn parse_line<D: Dictionary>(
    d: &mut D,
    line: &[u8],
    rec: &mut Record,
    sc: &mut Scratch,
) -> Result<()> {
    rec.clear();
    let mut cols = line.splitn(10, |&c| c == b'\t');
    let mut next = |what: &str| {
        cols.next()
            .ok_or_else(|| bad(format!("missing {what} column")))
    };
    let chrom = next("CHROM")?;
    let pos_field = next("POS")?;
    let id = next("ID")?;
    let reference = next("REF")?;
    let alt = next("ALT")?;
    let qual = next("QUAL")?;
    let filter = next("FILTER")?;
    let info = next("INFO")?;
    let format = cols.next();
    let samples = cols.next();

    let chrom = std::str::from_utf8(chrom).map_err(|_| bad("CHROM is not UTF-8"))?;
    rec.rid = d.contig(chrom)?;
    rec.pos = parse_pos(pos_field)?;
    rec.data.set_id(if id == b"." { b"" } else { id });
    rec.push_allele(reference);
    if alt != b"." {
        for a in alt.split(|&c| c == b',') {
            rec.push_allele(a);
        }
    }
    if rec.n_alleles() > u16::MAX as usize {
        return Err(bad("too many ALT alleles"));
    }
    rec.data.qual = if qual == b"." {
        FLOAT_MISSING
    } else {
        (strtod(qual).0 as f32).to_bits()
    };
    if filter != b"." {
        parse_filter(d, filter, rec)?;
    }
    if info != b"." {
        parse_info(d, info, rec, sc)?;
    }
    if !d.header().samples.is_empty() {
        match format {
            Some(format) => parse_format(d, format, samples, rec, sc)?,
            None if D::STRICT => return Err(bad("no FORMAT column; the header has samples")),
            None => {}
        }
    }
    let h = d.header();

    let mut spans = SpanFields::default();
    if let Some(end) = h.key_id("END").and_then(|id| rec.info_field(id)) {
        let v = end.values;
        if v.len() == 1
            && matches!(
                v.ty,
                BT_INT8 | BT_INT16 | BT_INT32 | crate::fields::BT_INT64
            )
        {
            spans.end = span_int(v.ty, v.bytes);
        }
    }
    if has_svlen_allele(rec)
        && let Some(svlen) = h.key_id("SVLEN").and_then(|id| rec.info_field(id))
    {
        spans.svlen = Some(span_values(svlen.values.ty, svlen.values.bytes));
    }
    if has_gvcf_allele(rec)
        && let Some(len) = h.key_id("LEN").and_then(|id| rec.format_field(id))
    {
        spans.len = Some(span_max(len.values.ty, len.values.bytes));
    }
    rec.rlen = rlen(rec.pos, rec, &spans)?;
    rec.end = spans.end;
    Ok(())
}

/// The largest FORMAT/LEN value over all samples, as the span needs it.
fn span_max(ty: u8, bytes: &[u8]) -> std::result::Result<i64, ()> {
    let w = crate::fields::type_width(ty);
    if !matches!(ty, BT_INT8 | BT_INT16 | BT_INT32 | crate::fields::BT_INT64) || w == 0 {
        return Err(());
    }
    Ok(bytes
        .chunks_exact(w)
        .filter_map(|b| span_int(ty, b))
        .max()
        .unwrap_or(0))
}

/// The meta-characters `vcf_parse_format_max3` stops at.
fn meta(c: u8) -> bool {
    matches!(c, 0 | b'\t' | b',' | b'/' | b':' | b'|')
}

/// FORMAT and the sample columns (`vcf_parse_format`): sizes are taken
/// over all samples first, then each sample's values are written at that
/// size and each key encoded for all samples.
fn parse_format<D: Dictionary>(
    d: &mut D,
    format: &[u8],
    samples: Option<&[u8]>,
    rec: &mut Record,
    sc: &mut Scratch,
) -> Result<()> {
    let Some(rest) = samples else {
        return Err(bad("FORMAT column with no sample columns"));
    };
    let nsamples = d.header().samples.len();
    if format == b"." {
        rec.data.n_sample = nsamples as u32;
        return Ok(());
    }

    // Keys from the dictionary.
    sc.aux.clear();
    for name in format.split(|&c| c == b':') {
        if sc.aux.len() >= MAX_N_FMT {
            return Err(bad("too many FORMAT keys"));
        }
        if name.is_empty() || name == b"." {
            return Err(bad("invalid FORMAT tag name '.'"));
        }
        let name = std::str::from_utf8(name).map_err(|_| bad("FORMAT tag is not UTF-8"))?;
        let (key, ty) = d.format(name)?;
        sc.aux.push(FmtAux {
            key,
            is_gt: name == "GT",
            htype: ty.code(),
            ..FmtAux::default()
        });
    }

    // Largest sizes over the samples; sample-separating tabs become NULs.
    sc.line.clear();
    sc.line.extend_from_slice(rest);
    let line = &mut sc.line;
    let aux = &mut sc.aux;
    let end = line.len();
    let c = |line: &[u8], i: usize| if i < end { line[i] } else { 0 };
    let mut n_sample = 0;
    let (mut m, mut g) = (1u32, 1u32);
    let mut r = 0;
    while r < end {
        let mut j = 0;
        let mut r_start = r;
        loop {
            while !meta(c(line, r)) {
                r += 1;
            }
            match c(line, r) {
                b',' => m += 1,
                b'|' | b'/' => {
                    if aux[j].is_gt {
                        g += 1
                    }
                }
                ch => {
                    if ch == 0 && r < end {
                        return Err(bad("unexpected character in FORMAT fields"));
                    }
                    if ch == b'\t' {
                        line[r] = 0;
                    }
                    // Fields after the first count their leading ':'.
                    let l = (r - r_start) as u32;
                    r_start = r;
                    let f = &mut aux[j];
                    f.max_m = f.max_m.max(m);
                    f.max_l = f.max_l.max(l);
                    if f.is_gt {
                        f.max_g = f.max_g.max(g);
                    }
                    m = 1;
                    g = 1;
                    if ch != b':' {
                        break;
                    }
                    j += 1;
                    if j >= aux.len() {
                        return Err(bad("more sample fields than FORMAT keys"));
                    }
                }
            }
            if r >= end {
                break;
            }
            r += 1;
        }
        n_sample += 1;
        if n_sample == nsamples {
            break;
        }
        r += 1;
    }

    // Arrays of the largest size per sample.
    sc.mem.clear();
    for f in aux.iter_mut() {
        if f.max_m == 0 {
            f.max_m = 1;
        }
        let size = match f.htype {
            HT_STR if f.is_gt => (f.max_g as usize) << 2,
            HT_STR => f.max_l as usize,
            HT_INT | HT_REAL => (f.max_m as usize) << 2,
            t => return Err(bad(format!("FORMAT type {t} is not supported"))),
        };
        f.size = Some(size);
        f.offset = sc.mem.len();
        let end = n_sample
            .checked_mul(size)
            .and_then(|bytes| sc.mem.len().checked_add(bytes))
            .filter(|&end| end <= i32::MAX as usize)
            .ok_or_else(|| bad("FORMAT fields too large"))?;
        sc.mem.resize(end, 0);
    }
    if aux.len() > 32 {
        if aux.iter().all(|f| f.key < 512) {
            // All word indices are below eight; each shift is below 64.
            // Marking in input order keeps the first occurrence unchanged.
            let mut seen = [0u64; 8];
            for f in aux.iter_mut() {
                let word = (f.key >> 6) as usize;
                let bit = 1u64 << (f.key & 63);
                if seen[word] & bit != 0 {
                    f.size = None;
                } else {
                    seen[word] |= bit;
                }
            }
        } else {
            // Sorting a side index handles arbitrary dictionary ids without
            // changing field order or retaining a dictionary-sized bitset.
            sc.key_order.clear();
            sc.key_order
                .extend(aux.iter().enumerate().map(|(i, f)| (f.key, i)));
            sc.key_order.sort_unstable();
            for pair in sc.key_order.windows(2) {
                if pair[0].0 == pair[1].0 {
                    aux[pair[1].1].size = None;
                }
            }
        }
    } else {
        for i in 1..aux.len() {
            if (0..i).any(|j| aux[j].size.is_some() && aux[j].key == aux[i].key) {
                aux[i].size = None;
            }
        }
    }

    if D::STRICT && aux.iter().any(|f| f.size.is_none()) {
        return Err(bad("duplicate FORMAT key"));
    }
    let h = d.header();
    fill(h.version, nsamples, line, end, aux, &mut sc.mem)?;

    // Encode each key for all samples.
    if n_sample > 0 {
        for f in aux.iter() {
            let Some(size) = f.size else { continue };
            let data = &sc.mem[f.offset..f.offset + size * n_sample];
            if f.htype == HT_STR && !f.is_gt {
                rec.data.push_format(f.key, BT_CHAR, size, data);
            } else if f.htype == HT_INT || f.is_gt {
                sc.ints.clear();
                sc.ints.extend(
                    data.chunks_exact(4)
                        .map(|b| i32::from_le_bytes(b.try_into().unwrap())),
                );
                let ty = enc_vint(&sc.ints, &mut sc.enc);
                rec.data.push_format(f.key, ty, size >> 2, &sc.enc);
            } else {
                rec.data.push_format(f.key, BT_FLOAT, size >> 2, data);
            }
        }
    }
    if n_sample != nsamples {
        return Err(bad(format!(
            "{n_sample} sample columns, the header has {nsamples}"
        )));
    }
    rec.data.n_sample = n_sample as u32;

    if !D::STRICT && h.has_gt_field() {
        let gt_id = h.key_id("GT");
        if let Some(gt) = gt_id.and_then(|id| rec.data.fmt.iter().find(|s| s.key == id).copied()) {
            let (out, data) = rec.gt_and_data(gt.n as usize);
            let bytes = &data.fmt_bytes[gt.start..gt.end];
            if gt.ty == BT_NULL {
                out.clear();
            } else {
                decode_gt(gt.ty, gt.n as usize, bytes, n_sample, out);
            }
        }
    }
    Ok(())
}

/// Writes each sample's values at the per-key size
/// (`vcf_parse_format_fill5`).
fn fill(
    version: i32,
    nsamples: usize,
    line: &[u8],
    end: usize,
    aux: &[FmtAux],
    mem: &mut [u8],
) -> Result<()> {
    let c = |i: usize| if i < end { line[i] } else { 0 };
    let put = |mem: &mut [u8], base: usize, size: usize, l: usize, v: u32| -> Result<()> {
        if (l + 1) * 4 > size {
            return Err(bad("FORMAT value outside its field"));
        }
        mem[base + 4 * l..base + 4 * l + 4].copy_from_slice(&v.to_le_bytes());
        Ok(())
    };
    let mut t = 0;
    let mut m = 0;
    while t < end {
        if m == nsamples {
            break;
        }
        let mut j = 0;
        while t < end {
            let f = aux
                .get(j)
                .ok_or_else(|| bad("more sample fields than FORMAT keys"))?;
            j += 1;
            match f.size {
                None => {
                    while c(t) != b':' && c(t) != 0 {
                        t += 1;
                    }
                }
                Some(size) => {
                    let base = f.offset + size * m;
                    match f.htype {
                        HT_STR if f.is_gt => {
                            let mut is_phased = 0u32;
                            let mut prefixed = false;
                            if version >= VCF44 && matches!(c(t), b'|' | b'/') {
                                is_phased = u32::from(c(t) == b'|');
                                t += 1;
                                prefixed = true;
                            }
                            let (mut l, mut ploidy) = (0, 0);
                            let (mut unreadable, mut overflow) = (false, false);
                            let (mut any_unphased, mut unknown_first) = (false, false);
                            let mut max = 0u32;
                            let mut first = 0u32;
                            loop {
                                ploidy += 1;
                                let v = if c(t) == b'.' {
                                    t += 1;
                                    if l == 0 {
                                        unknown_first = true;
                                    }
                                    is_phased
                                } else {
                                    let val = if c(t).is_ascii_digit() && !c(t + 1).is_ascii_digit()
                                    {
                                        t += 1;
                                        u32::from(c(t - 1) - b'0')
                                    } else {
                                        let (n, used, of) =
                                            hts_str2uint(&line[t.min(end)..end], 506);
                                        overflow |= of;
                                        unreadable |= used == 0;
                                        t += used;
                                        n as u32
                                    };
                                    max = max.max(val);
                                    (val.wrapping_add(1) << 1) | is_phased
                                };
                                if l == 0 {
                                    first = v;
                                }
                                put(mem, base, size, l, v)?;
                                l += 1;
                                any_unphased |= ploidy != 1 && is_phased == 0;
                                is_phased = u32::from(c(t) == b'|');
                                if !matches!(c(t), b'|' | b'/') {
                                    break;
                                }
                                t += 1;
                            }
                            if !prefixed {
                                // Without a leading indicator, the first
                                // allele is phased when haploid and known,
                                // or when no other allele is unphased.
                                let phased = if ploidy == 1 {
                                    !unknown_first
                                } else {
                                    !any_unphased
                                };
                                if phased {
                                    put(mem, base, size, 0, first | 1)?;
                                }
                            }
                            if overflow || max > (i32::MAX >> 1) as u32 - 1 {
                                return Err(bad("GT value too large"));
                            }
                            if unreadable {
                                return Err(bad("GT value not a number or '.'"));
                            }
                            for k in l..size / 4 {
                                put(mem, base, size, k, INT32_VECTOR_END as u32)?;
                            }
                        }
                        HT_STR => {
                            let mut l = 0;
                            while c(t) != b':' && c(t) != 0 {
                                if l >= size {
                                    return Err(bad("FORMAT value outside its field"));
                                }
                                mem[base + l] = c(t);
                                l += 1;
                                t += 1;
                            }
                            mem[base + l..base + size].fill(0);
                        }
                        HT_INT => {
                            let mut l = 0;
                            loop {
                                let v = if c(t) == b'.' {
                                    t += 1;
                                    INT32_MISSING
                                } else {
                                    let (x, used, of) = hts_str2int(&line[t.min(end)..end], 64);
                                    t += used;
                                    if used == 0
                                        || of
                                        || !(MIN_INT32..=i32::MAX as i64).contains(&x)
                                    {
                                        INT32_MISSING
                                    } else {
                                        x as i32
                                    }
                                };
                                put(mem, base, size, l, v as u32)?;
                                l += 1;
                                if c(t) != b',' {
                                    break;
                                }
                                t += 1;
                            }
                            for k in l..size / 4 {
                                put(mem, base, size, k, INT32_VECTOR_END as u32)?;
                            }
                        }
                        HT_REAL => {
                            let mut l = 0;
                            loop {
                                let v = if c(t) == b'.' && !c(t + 1).is_ascii_digit() {
                                    t += 1;
                                    FLOAT_MISSING
                                } else {
                                    // A failed conversion stores 0.
                                    let (x, used, _) = hts_str2dbl(&line[t.min(end)..end]);
                                    t += used;
                                    (x as f32).to_bits()
                                };
                                put(mem, base, size, l, v)?;
                                l += 1;
                                if c(t) != b',' {
                                    break;
                                }
                                t += 1;
                            }
                            for k in l..size / 4 {
                                put(mem, base, size, k, FLOAT_VECTOR_END)?;
                            }
                        }
                        t => return Err(bad(format!("FORMAT type {t} is not supported"))),
                    }
                }
            }
            match c(t) {
                0 => break,
                b':' => t += 1,
                _ => return Err(bad("invalid character in FORMAT field")),
            }
        }
        // Keys a sample omits: missing, then end-of-vector.
        for f in &aux[j..] {
            let Some(size) = f.size else { continue };
            let base = f.offset + size * m;
            match f.htype {
                HT_STR if f.is_gt => {
                    if size > 0 {
                        put(mem, base, size, 0, INT32_MISSING as u32)?;
                    }
                    for k in 1..size / 4 {
                        put(mem, base, size, k, INT32_VECTOR_END as u32)?;
                    }
                }
                HT_STR => {
                    if size > 0 {
                        mem[base] = b'.';
                        mem[base + 1..base + size].fill(0);
                    }
                }
                HT_INT => {
                    put(mem, base, size, 0, INT32_MISSING as u32)?;
                    for k in 1..size / 4 {
                        put(mem, base, size, k, INT32_VECTOR_END as u32)?;
                    }
                }
                _ => {
                    put(mem, base, size, 0, FLOAT_MISSING)?;
                    for k in 1..size / 4 {
                        put(mem, base, size, k, FLOAT_VECTOR_END)?;
                    }
                }
            }
        }
        m += 1;
        t += 1;
    }
    Ok(())
}
