//! `filter_test` and what it calls: the setters that read a record into
//! value tokens, the functions, the comparators and the operators.

use omicsweave_core::{Error, Result};

use super::stats::{calc_binom_two_sided, kt_fisher_two};
use super::*;
use crate::fields::{
    BT_CHAR, BT_FLOAT, BT_INT8, BT_INT16, BT_INT32, BT_NULL, FLOAT_MISSING, FLOAT_VECTOR_END,
    FormatField, INT32_MISSING, INT32_VECTOR_END,
};
use crate::record::Record;

fn run_error(detail: impl Into<String>) -> Error {
    Error::invalid("filter evaluation", detail)
}

fn one_based_pos(rec: &Record) -> Result<i64> {
    rec.pos
        .checked_add(1)
        .ok_or_else(|| unsupported("filter position overflow"))
}

/// C's `(int)` of a double as x86-64 computes it: truncation, and
/// `INT_MIN` for NaN and values outside the range.
fn c_int(x: f64) -> i32 {
    if x.is_nan() || x >= 2147483648.0 || x <= -2147483649.0 {
        i32::MIN
    } else {
        x as i32
    }
}

/// `1ULL << k` as x86-64 shifts (the count taken modulo 64).
fn bit(k: i64) -> u64 {
    1u64.wrapping_shl(k as u32)
}

/// Mutable access to token `r` and shared access to the others.
fn split(v: &mut [Token], r: usize) -> (&[Token], &mut Token) {
    let (left, right) = v.split_at_mut(r);
    (left, &mut right[0])
}

/// Mutable access to two distinct tokens, `a < r`.
fn pair(v: &mut [Token], a: usize, r: usize) -> (&mut Token, &mut Token) {
    let (left, right) = v.split_at_mut(r);
    (&mut left[a], &mut right[0])
}

fn gt_missing(v: i64) -> bool {
    (v >> 1) == 0
}

fn gt_allele(v: i64) -> i64 {
    (v >> 1) - 1
}

/// The vector-end value of an integer type, as stored.
fn int_vector_end(ty: u8) -> i64 {
    match ty {
        BT_INT8 => i8::MIN as i64 + 1,
        BT_INT16 => i16::MIN as i64 + 1,
        _ => i32::MIN as i64 + 1,
    }
}

/// `bcf_get_fmt(hdr, line, "GT")`.
fn gt_field<'a>(h: &Header, rec: &'a Record) -> Option<FormatField<'a>> {
    let id = h.key_id("GT")?;
    if !h.declared_as(id, LineKind::Format) {
        return None;
    }
    rec.format_field(id)
}

/// `bcf_format_gt(fmt, i, str)` (no header: no VCF 4.4 phase prefix).
fn format_gt(fmt: &FormatField, i: usize, s: &mut KStr) {
    match fmt.values.ty {
        BT_INT8 | BT_INT16 | BT_INT32 => {
            let row = fmt.sample(i);
            let end = int_vector_end(fmt.values.ty);
            let mut k = 0;
            while k < fmt.n {
                let v = row.raw(k);
                if v == end {
                    break;
                }
                if k > 0 {
                    s.putc(if v & 1 != 0 { b'|' } else { b'/' });
                }
                if (v >> 1) == 0 {
                    s.putc(b'.');
                } else {
                    for c in ((v >> 1) - 1).to_string().bytes() {
                        s.putc(c);
                    }
                }
                k += 1;
            }
            if k == 0 {
                s.putc(b'.');
            }
        }
        BT_NULL => s.putc(b'.'),
        _ => {}
    }
}

/// filter.c's `bcf_get_info_value`: one INFO value as it lands in the
/// caller's variable.
enum InfoValue {
    Absent,
    Missing,
    Int(i64),
    Float(f64),
    /// Nothing written to the variable.
    Untouched,
}

fn info_value(rec: &Record, key: i32, ivec: i32) -> Result<InfoValue> {
    let Some(f) = (key >= 0).then(|| rec.info_field(key as u32)).flatten() else {
        return Ok(InfoValue::Absent);
    };
    let v = f.values;
    let len = v.len();
    if len == 1 {
        return Ok(match v.ty {
            BT_FLOAT => InfoValue::Float(f64::from(f.v1_float())),
            BT_INT8 | BT_INT16 | BT_INT32 => InfoValue::Int(f.v1_int()),
            _ => InfoValue::Untouched,
        });
    }
    let ivec = ivec.max(0) as usize;
    if ivec >= len {
        return Ok(InfoValue::Missing);
    }
    match v.ty {
        BT_INT8 | BT_INT16 | BT_INT32 | BT_FLOAT => {}
        _ => return Err(unsupported("filter.c exits on an INFO value of this type")),
    }
    let (mis, end) = if v.ty == BT_FLOAT {
        (FLOAT_MISSING as i64, FLOAT_VECTOR_END as i64)
    } else {
        let m = match v.ty {
            BT_INT8 => i8::MIN as i64,
            BT_INT16 => i16::MIN as i64,
            _ => i32::MIN as i64,
        };
        (m, m + 1)
    };
    for j in 0..ivec {
        if v.raw(j) == end {
            return Ok(InfoValue::Missing);
        }
    }
    let raw = v.raw(ivec);
    if raw == mis {
        return Ok(InfoValue::Missing);
    }
    Ok(if v.ty == BT_FLOAT {
        InfoValue::Float(f64::from(f32::from_bits(raw as u32)))
    } else {
        InfoValue::Int(raw)
    })
}

fn tag_str(tok: &Token) -> &str {
    tok.tag
        .as_deref()
        .and_then(|t| std::str::from_utf8(t).ok())
        .unwrap_or("")
}

/// Keeps the stale contents of filter.c's shared integer buffer: a getter
/// grows it to `need` and writes `values`.
fn fill_ints(buf: &mut Vec<Option<i32>>, need: usize, values: &[i32]) {
    if buf.len() < need {
        buf.resize(need, None);
    }
    if buf.len() < values.len() {
        buf.resize(values.len(), None);
    }
    for (k, &v) in values.iter().enumerate() {
        buf[k] = Some(v);
    }
}

fn fill_floats(buf: &mut Vec<Option<u32>>, need: usize, values: &[f32]) {
    if buf.len() < need {
        buf.resize(need, None);
    }
    if buf.len() < values.len() {
        buf.resize(values.len(), None);
    }
    for (k, &v) in values.iter().enumerate() {
        buf[k] = Some(v.to_bits());
    }
}

fn read_int(buf: &[Option<i32>], i: isize) -> i32 {
    match usize::try_from(i)
        .ok()
        .and_then(|i| buf.get(i))
        .copied()
        .flatten()
    {
        Some(v) => v,
        None => {
            undefined("read of the integer buffer past what was written");
            0
        }
    }
}

fn read_float(buf: &[Option<u32>], i: isize) -> u32 {
    match usize::try_from(i)
        .ok()
        .and_then(|i| buf.get(i))
        .copied()
        .flatten()
    {
        Some(v) => v,
        None => {
            undefined("read of the float buffer past what was written");
            0
        }
    }
}

/// `_has_missing_string`: whether a comma-separated C string has a `.`.
fn has_missing_string(s: &KStr, mut beg: usize) -> bool {
    while s.get(beg) != 0 {
        let mut end = beg;
        while s.get(end) != 0 && s.get(end) != b',' {
            end += 1;
        }
        if end - beg == 1 && s.get(beg) == b'.' {
            return true;
        }
        if s.get(end) == 0 {
            break;
        }
        beg = end + 1;
    }
    false
}

/// `_match_vector_strings`: whether any field of `a` matches any of `b`.
// Preserve the two independent string spans and the C comparison flags.
#[allow(clippy::too_many_arguments)]
fn match_vector_strings(
    a: &KStr,
    mut abeg: usize,
    alen: usize,
    b: &KStr,
    bstr: usize,
    blen: usize,
    logic: i32,
    ml: &[i32; 3],
) -> bool {
    let aend = abeg + alen;
    let bend = bstr + blen;
    while abeg < aend && a.get(abeg) != 0 {
        let mut amid = abeg;
        while amid < aend && a.get(amid) != 0 && a.get(amid) != b',' {
            amid += 1;
        }
        let miss = usize::from(amid - abeg == 1 && a.get(abeg) == b'.');
        let mut bbeg = bstr;
        while bbeg < bend && b.get(bbeg) != 0 {
            let mut bmid = bbeg;
            while bmid < bend && b.get(bmid) != 0 && b.get(bmid) != b',' {
                bmid += 1;
            }
            let nmiss = miss + usize::from(bmid - bbeg == 1 && b.get(bbeg) == b'.');
            if nmiss > 0 {
                if ml[nmiss] != 0 {
                    return true;
                }
            } else {
                let mut m = amid - abeg == bmid - bbeg
                    && (0..amid - abeg).all(|k| {
                        let (x, y) = (a.get(abeg + k), b.get(bbeg + k));
                        x == y
                    });
                if logic == TOK_NE {
                    m = !m;
                }
                if m {
                    return true;
                }
            }
            if b.get(bmid) == 0 {
                break;
            }
            bbeg = bmid + 1;
        }
        if a.get(amid) == 0 {
            break;
        }
        abeg = amid + 1;
    }
    false
}

/// `_regex_vector_strings`: whether any field matches (or, for `!~`, does
/// not match).
fn regex_vector_strings(
    re: &Regex,
    s: &KStr,
    mut str: usize,
    len: usize,
    logic: i32,
    ml: &[i32; 3],
) -> bool {
    let end = str + len;
    while str < end && s.get(str) != 0 {
        let mut mid = str;
        while mid < end && s.get(mid) != 0 && s.get(mid) != b',' {
            mid += 1;
        }
        let miss = usize::from(mid - str == 1 && s.get(str) == b'.');
        let mut m = miss > 0 && ml[miss] != 0;
        if logic == TOK_NLIKE {
            m = !m;
        }
        if m {
            return true;
        }
        let field: Vec<u8> = (str..mid).map(|k| s.get(k)).collect();
        let mut m = re.is_match(&field);
        if logic == TOK_NLIKE {
            m = !m;
        }
        if m {
            return true;
        }
        if s.get(mid) == 0 {
            break;
        }
        str = mid + 1;
    }
    false
}

/// `tok_init_samples`.
fn init_samples(a: &Token, b: &Token, r: &mut Token) {
    if (a.nsamples != 0 || b.nsamples != 0) && r.nsamples == 0 {
        r.nsamples = if a.nsamples != 0 {
            a.nsamples
        } else {
            b.nsamples
        };
        let mut u = vec![0u8; r.nsamples as usize];
        for t in [a, b] {
            for i in 0..t.nsamples as usize {
                let v = u8::from(t.usmpl(i));
                match u.get_mut(i) {
                    Some(x) => *x |= v,
                    None => undefined("sample mask written past its end"),
                }
            }
        }
        r.usmpl = Some(u);
    }
    for i in 0..r.nsamples.max(0) as usize {
        r.set_pass_sample(i, 0);
    }
}

/// `tok_init_values`.
fn init_values(a: &Token, b: &Token, r: &mut Token) {
    let t = if (a.nsamples != 0 || b.nsamples != 0) && (a.nsamples == 0 || b.nsamples == 0) {
        if a.nsamples != 0 { a } else { b }
    } else if a.nvalues > b.nvalues {
        a
    } else {
        b
    };
    r.nvalues = t.nvalues;
    r.nval1 = t.nval1;
    r.values.expand(r.nvalues.max(0) as usize);
}

/// `VECTOR_ARITHMETICS`.
fn arith(a: &Token, b: &Token, r: &mut Token, op: i32) -> Result<()> {
    let apply = |x: f64, y: f64| -> Result<f64> {
        Ok(match op {
            TOK_ADD => x + y,
            TOK_SUB => x - y,
            TOK_MULT => x * y,
            TOK_DIV => x / y,
            _ => {
                let (x, y) = (c_int(x), c_int(y));
                if y == 0 || (x == i32::MIN && y == -1) {
                    return Err(unsupported("integer division by zero in `%`"));
                }
                f64::from(x % y)
            }
        })
    };
    let mut has = false;
    if a.nvalues != 0 && b.nvalues != 0 {
        init_values(a, b, r);
        init_samples(a, b, r);
        if a.nsamples == 0 && b.nsamples == 0 {
            if a.nvalues != b.nvalues && a.nvalues != 1 && b.nvalues != 1 {
                return Err(run_error(format!(
                    "Cannot run numeric operator in -i/-e filtering on vectors of different lengths: {} vs {}",
                    a.nvalues, b.nvalues
                )));
            }
            let (mut ia, mut ib) = (0, 0);
            for ir in 0..r.nvalues as isize {
                if a.nvalues > 1 {
                    ia = ir;
                }
                if b.nvalues > 1 {
                    ib = ir;
                }
                let (x, y) = (a.values.get(ia), b.values.get(ib));
                if is_mv(x) || is_mv(y) {
                    r.values.set(ir, missing());
                    continue;
                }
                has = true;
                r.values.set(ir, apply(x, y)?);
            }
        } else if a.nsamples != 0 && b.nsamples != 0 {
            if a.nsamples != b.nsamples {
                return Err(unsupported("assertion failed in filter.c"));
            }
            if a.nval1 != b.nval1 && a.nval1 != 1 && b.nval1 != 1 {
                return Err(run_error(format!(
                    "Cannot run numeric operator in -i/-e filtering on vectors of different lengths: {} vs {}",
                    a.nval1, b.nval1
                )));
            }
            for i in 0..r.nsamples as isize {
                let (rv, av, bv) = (
                    i * r.nval1 as isize,
                    i * a.nval1 as isize,
                    i * b.nval1 as isize,
                );
                let (mut ia, mut ib) = (0, 0);
                for ir in 0..r.nval1 as isize {
                    if a.nval1 > 1 {
                        ia = ir;
                    }
                    if b.nval1 > 1 {
                        ib = ir;
                    }
                    let (x, y) = (a.values.get(av + ia), b.values.get(bv + ib));
                    if is_mv(x) || is_mv(y) {
                        r.values.set(rv + ir, missing());
                        continue;
                    }
                    has = true;
                    r.values.set(rv + ir, apply(x, y)?);
                }
            }
        } else if a.nsamples != 0 {
            if b.nvalues != 1 {
                return Err(unsupported("assertion failed in filter.c"));
            }
            let y = b.values.get(0);
            if !is_mv(y) {
                for i in 0..a.nvalues as isize {
                    let x = a.values.get(i);
                    if is_mv(x) {
                        r.values.set(i, missing());
                        continue;
                    }
                    has = true;
                    r.values.set(i, apply(x, y)?);
                }
            }
        } else {
            if a.nvalues != 1 {
                return Err(unsupported("assertion failed in filter.c"));
            }
            let x = a.values.get(0);
            if !is_mv(x) {
                for i in 0..b.nvalues as isize {
                    let y = b.values.get(i);
                    if is_mv(y) {
                        r.values.set(i, missing());
                        continue;
                    }
                    has = true;
                    r.values.set(i, apply(x, y)?);
                }
            }
        }
    }
    if !has {
        r.nvalues = 0;
    }
    Ok(())
}

fn compare(op: i32, x: f64, y: f64) -> bool {
    // Values up to 2^24 compare as floats (bcftools #871).
    let (x, y) = if x > 16777216.0 || y > 16777216.0 {
        (x, y)
    } else {
        (f64::from(x as f32), f64::from(y as f32))
    };
    match op {
        TOK_EQ => x == y,
        TOK_NE => x != y,
        TOK_LE => x <= y,
        TOK_LT => x < y,
        TOK_BT => x > y,
        _ => x >= y,
    }
}

/// A relation is existential over all pairs. Its ordered, nonmissing part
/// depends only on the two extrema, including filter.c's float rounding.
fn ordered_vectors(a: &Token, b: &Token, op: i32, ml: &[i32; 3]) -> Option<bool> {
    if !matches!(op, TOK_LT | TOK_LE | TOK_BT | TOK_BE) {
        return None;
    }
    let (x, y) = (*a.values.0.first()?, *b.values.0.first()?);
    if x.to_bits() == UNINIT || y.to_bits() == UNINIT {
        return None;
    }
    let nmiss = usize::from(is_mv(x)) + usize::from(is_mv(y));
    if if nmiss == 0 {
        compare(op, x, y)
    } else {
        ml[nmiss] != 0
    } {
        return Some(true);
    }
    struct Range {
        missing: bool,
        present: bool,
        min: f64,
        max: f64,
    }
    fn range(t: &Token) -> Option<Range> {
        let n = usize::try_from(t.nvalues).ok()?;
        let values = t.values.0.get(..n)?;
        let mut out = Range {
            missing: false,
            present: false,
            min: f64::INFINITY,
            max: f64::NEG_INFINITY,
        };
        for &v in values {
            // Leave partial shadow arrays on the original short-circuit
            // path: an unwritten tail need not be read by filter.c.
            if v.to_bits() == UNINIT {
                return None;
            }
            if is_mv(v) {
                out.missing = true;
            } else {
                out.present = true;
                if !v.is_nan() {
                    out.min = out.min.min(v);
                    out.max = out.max.max(v);
                }
            }
        }
        Some(out)
    }
    let (a, b) = (range(a)?, range(b)?);
    if (ml[2] != 0 && a.missing && b.missing)
        || (ml[1] != 0 && ((a.missing && b.present) || (a.present && b.missing)))
    {
        return Some(true);
    }
    // NaNs have no ordered pairs. The inverted bounds also distinguish
    // all-NaN vectors from actual infinities.
    if a.min > a.max || b.min > b.max {
        return Some(false);
    }
    let (x, y) = match op {
        TOK_LT | TOK_LE => (a.min, b.max),
        _ => (a.max, b.min),
    };
    Some(compare(op, x, y))
}

/// Equality has two disjoint numeric domains: values above 2^24 retain
/// double precision; all others compare after float rounding. Rounding
/// cannot cross that boundary, and signed zeros compare equal.
fn equality_key(value: f64) -> Option<u64> {
    if value.is_nan() {
        return None;
    }
    let value = if value > 16777216.0 {
        value
    } else {
        f64::from(value as f32)
    };
    Some(if value == 0.0 { 0 } else { value.to_bits() })
}

/// Reduce existential equality/inequality without enumerating all pairs.
/// Inspect only initialized shadow slots; incomplete arrays must retain the
/// original access order, including successful matches before unwritten tails.
fn equality_vectors(a: &Token, b: &Token, op: i32, ml: &[i32; 3]) -> Option<bool> {
    let av = a.values.0.get(..usize::try_from(a.nvalues).ok()?)?;
    let bv = b.values.0.get(..usize::try_from(b.nvalues).ok()?)?;
    let &x = av.first()?;
    if x.to_bits() == UNINIT {
        return None;
    }
    // Preserve cheap matches in the first row before building summaries or
    // a hash table. Every visited pair follows the original access order.
    for &y in bv {
        if y.to_bits() == UNINIT {
            return None;
        }
        let nmiss = usize::from(is_mv(x)) + usize::from(is_mv(y));
        if if nmiss == 0 {
            compare(op, x, y)
        } else {
            ml[nmiss] != 0
        } {
            return Some(true);
        }
    }
    #[derive(Default)]
    struct Summary {
        missing: bool,
        present: bool,
        nan: bool,
        key: Option<u64>,
        varied: bool,
    }
    fn summary(values: &[f64]) -> Option<Summary> {
        let mut out = Summary::default();
        for &value in values {
            if value.to_bits() == UNINIT {
                return None;
            }
            if is_mv(value) {
                out.missing = true;
            } else {
                out.present = true;
                match equality_key(value) {
                    None => out.nan = true,
                    Some(key) => match out.key {
                        None => out.key = Some(key),
                        Some(first) => out.varied |= first != key,
                    },
                }
            }
        }
        Some(out)
    }
    let (a, b) = (summary(av)?, summary(bv)?);
    if (ml[2] != 0 && a.missing && b.missing)
        || (ml[1] != 0 && ((a.missing && b.present) || (a.present && b.missing)))
    {
        return Some(true);
    }
    if !a.present || !b.present {
        return Some(false);
    }
    if op == TOK_NE {
        return Some(a.nan || b.nan || a.varied || b.varied || a.key != b.key);
    }
    let (Some(ak), Some(bk)) = (a.key, b.key) else {
        return Some(false);
    };
    if !a.varied {
        return Some(bv.iter().any(|&v| !is_mv(v) && equality_key(v) == Some(ak)));
    }
    if !b.varied {
        return Some(av.iter().any(|&v| !is_mv(v) && equality_key(v) == Some(bk)));
    }
    let (small, large) = if av.len() <= bv.len() {
        (av, bv)
    } else {
        (bv, av)
    };
    let keys: HashSet<u64> = small
        .iter()
        .filter(|&&v| !is_mv(v))
        .filter_map(|&v| equality_key(v))
        .collect();
    Some(
        large
            .iter()
            .any(|&v| !is_mv(v) && equality_key(v).is_some_and(|key| keys.contains(&key))),
    )
}

/// `CMP_VECTORS`.
fn cmp_vectors(a: &Token, b: &Token, r: &mut Token, op: i32, ml: &[i32; 3]) -> Result<()> {
    init_samples(a, b, r);
    let mv = |t: &Token, i: isize| usize::from(is_mv(t.values.get(i)));
    if a.nsamples == 0 && b.nsamples == 0 {
        if a.nvalues >= 64 && b.nvalues >= 64 && matches!(op, TOK_EQ | TOK_NE) {
            if let Some(pass) = equality_vectors(a, b, op, ml) {
                r.pass_site = i32::from(pass);
                return Ok(());
            }
        }
        if a.nvalues >= 16 && b.nvalues >= 16 {
            if let Some(pass) = ordered_vectors(a, b, op, ml) {
                r.pass_site = i32::from(pass);
                return Ok(());
            }
        }
        if a.nvalues == 0 && b.nvalues == 0 {
            r.pass_site = ml[2];
        } else if a.nvalues == 0 || b.nvalues == 0 {
            let t = if a.nvalues != 0 { a } else { b };
            for j in 0..t.nvalues as isize {
                if mv(t, j) != 0 {
                    if ml[2] != 0 {
                        r.pass_site = 1;
                        break;
                    }
                } else if ml[1] != 0 {
                    r.pass_site = 1;
                    break;
                }
            }
        } else {
            'outer: for i in 0..a.nvalues as isize {
                let amiss = mv(a, i);
                for j in 0..b.nvalues as isize {
                    let nmiss = amiss + mv(b, j);
                    if nmiss != 0 {
                        if ml[nmiss] != 0 {
                            r.pass_site = 1;
                            break 'outer;
                        }
                    } else if compare(op, a.values.get(i), b.values.get(j)) {
                        r.pass_site = 1;
                        break 'outer;
                    }
                }
            }
        }
        return Ok(());
    }
    if a.nvalues == 0 && b.nvalues == 0 {
        if ml[2] != 0 {
            for i in 0..r.nsamples as usize {
                if r.usmpl(i) {
                    r.set_pass_sample(i, ml[2] as u8);
                    r.pass_site = 1;
                }
            }
        }
    } else if a.nvalues == 0 || b.nvalues == 0 {
        let t = if a.nvalues != 0 { a } else { b };
        if t.nsamples == 0 {
            let mut miss = 0;
            for j in 0..t.nvalues as isize {
                miss |= mv(t, j);
            }
            let k = miss + 1;
            if ml[k] != 0 {
                for i in 0..r.nsamples as usize {
                    if r.usmpl(i) {
                        r.set_pass_sample(i, ml[k] as u8);
                        r.pass_site = 1;
                    }
                }
            }
        } else {
            for i in 0..t.nsamples as usize {
                if !r.usmpl(i) {
                    continue;
                }
                let p = (i * t.nval1 as usize) as isize;
                let mut miss = 0;
                for j in 0..t.nval1 as isize {
                    miss |= mv(t, p + j);
                }
                let k = miss + 1;
                if ml[k] != 0 {
                    r.set_pass_sample(i, ml[k] as u8);
                    r.pass_site = 1;
                }
            }
        }
    } else if a.nsamples != 0 && b.nsamples != 0 {
        if a.nval1 != b.nval1 {
            return Err(run_error(format!(
                "Incompatible number of per-sample values in comparison: {} vs {}",
                a.nval1, b.nval1
            )));
        }
        if a.nsamples != b.nsamples {
            return Err(run_error(format!(
                "Incompatible number samples in comparison: {} vs {}",
                a.nsamples, b.nsamples
            )));
        }
        for i in 0..a.nsamples as usize {
            if !a.usmpl(i) || !b.usmpl(i) {
                match r.usmpl.as_mut().and_then(|u| u.get_mut(i)) {
                    Some(x) => *x = 0,
                    None => undefined("sample mask written past its end"),
                }
                continue;
            }
            let ap = (i * a.nval1 as usize) as isize;
            let bp = (i * b.nval1 as usize) as isize;
            for j in 0..a.nval1 as isize {
                let (x, y) = (a.values.get(ap + j), b.values.get(bp + j));
                if a.idx == -3 && is_vector_end(x) {
                    break;
                }
                if b.idx == -3 && is_vector_end(y) {
                    break;
                }
                let mut nmiss = usize::from(is_mv(x));
                if nmiss != 0 && ml[0] == 0 {
                    continue;
                }
                nmiss += usize::from(is_mv(y));
                if nmiss != 0 {
                    if ml[nmiss] != 0 {
                        r.set_pass_sample(i, 1);
                        r.pass_site = 1;
                        break;
                    }
                } else if compare(op, x, y) {
                    r.set_pass_sample(i, 1);
                    r.pass_site = 1;
                    break;
                }
            }
        }
    } else {
        // One side per sample, the other a site-level vector.
        let (s, v, s_first) = if a.nsamples != 0 {
            (a, b, true)
        } else {
            (b, a, false)
        };
        for i in 0..s.nsamples as usize {
            if !r.usmpl(i) {
                continue;
            }
            let p = (i * s.nval1 as usize) as isize;
            'sample: for j in 0..s.nval1 as isize {
                let x = s.values.get(p + j);
                if a.idx == -3 && is_vector_end(x) {
                    break;
                }
                let miss = usize::from(is_mv(x));
                if miss != 0 && ml[0] == 0 {
                    continue;
                }
                for k in 0..v.nvalues as isize {
                    let y = v.values.get(k);
                    let nmiss = miss + usize::from(is_mv(y));
                    let hit = if nmiss != 0 {
                        ml[nmiss] != 0
                    } else if s_first {
                        compare(op, x, y)
                    } else {
                        compare(op, y, x)
                    };
                    if hit {
                        r.set_pass_sample(i, 1);
                        r.pass_site = 1;
                        break 'sample;
                    }
                }
            }
        }
    }
    Ok(())
}

/// `cmp_vector_strings`.
fn cmp_vector_strings(a: &Token, b: &Token, r: &mut Token) -> Result<()> {
    init_samples(a, b, r);
    let logic = r.tok_type;
    let regex = a.regex.as_deref().or(b.regex.as_deref());
    if a.nvalues as usize != a.str_value.l || b.nvalues as usize != b.str_value.l {
        return Err(unsupported("assertion failed in filter.c"));
    }
    if a.nsamples != 0 && b.nsamples != 0 {
        return Err(unsupported("assertion failed in filter.c"));
    }
    let ok = match regex {
        None => logic == TOK_EQ || logic == TOK_NE,
        Some(_) => logic == TOK_LIKE || logic == TOK_NLIKE,
    };
    if !ok {
        return Err(unsupported("assertion failed in filter.c"));
    }
    let mut ml = [0, 0, 0];
    if logic == TOK_EQ || logic == TOK_LIKE {
        ml[0] = 1;
        ml[2] = 1;
    } else if logic == TOK_NE || logic == TOK_NLIKE {
        ml[0] = 1;
        ml[1] = 1;
    }
    if a.nsamples == 0 && b.nsamples == 0 {
        if a.nvalues == 0 && b.nvalues == 0 {
            r.pass_site = ml[2];
            return Ok(());
        }
        if a.nvalues == 0 || b.nvalues == 0 {
            let t = if a.nvalues != 0 { a } else { b };
            let miss = usize::from(has_missing_string(&t.str_value, 0));
            if ml[miss + 1] != 0 {
                r.pass_site = 1;
            }
            return Ok(());
        }
        r.pass_site = i32::from(match regex {
            None => match_vector_strings(
                &a.str_value,
                0,
                a.str_value.l,
                &b.str_value,
                0,
                b.str_value.l,
                logic,
                &ml,
            ),
            Some(re) => {
                let t = if a.regex.is_some() { b } else { a };
                regex_vector_strings(re, &t.str_value, 0, t.str_value.l, logic, &ml)
            }
        });
        return Ok(());
    }
    if a.nvalues == 0 && b.nvalues == 0 {
        if ml[2] != 0 {
            for i in 0..r.nsamples as usize {
                if r.usmpl(i) {
                    r.set_pass_sample(i, ml[2] as u8);
                    r.pass_site = 1;
                }
            }
        }
        return Ok(());
    }
    if a.nvalues == 0 || b.nvalues == 0 {
        let t = if a.nvalues != 0 { a } else { b };
        if t.nsamples == 0 {
            let miss = usize::from(has_missing_string(&t.str_value, 0));
            if ml[miss + 1] == 0 {
                return Ok(());
            }
            for i in 0..r.nsamples as usize {
                if r.usmpl(i) {
                    r.set_pass_sample(i, 1);
                    r.pass_site = 1;
                }
            }
        } else {
            for i in 0..t.nsamples as usize {
                if !r.usmpl(i) {
                    continue;
                }
                let miss = usize::from(has_missing_string(&t.str_value, i * t.nval1 as usize));
                if ml[miss + 1] != 0 {
                    r.set_pass_sample(i, 1);
                    r.pass_site = 1;
                }
            }
        }
        return Ok(());
    }
    let (x, y) = if a.nsamples != 0 { (a, b) } else { (b, a) };
    let same = match (regex, y.regex.as_deref()) {
        (None, None) => true,
        (Some(p), Some(q)) => std::ptr::eq(p, q),
        _ => false,
    };
    if !same {
        return Err(unsupported("assertion failed in filter.c"));
    }
    for i in 0..x.nsamples as usize {
        if !r.usmpl(i) {
            continue;
        }
        let at = i * x.nval1 as usize;
        let m = match regex {
            Some(re) => regex_vector_strings(re, &x.str_value, at, x.nval1 as usize, logic, &ml),
            None => match_vector_strings(
                &x.str_value,
                at,
                x.nval1 as usize,
                &y.str_value,
                0,
                y.str_value.l,
                logic,
                &ml,
            ),
        };
        if m {
            r.set_pass_sample(i, 1);
            r.pass_site = 1;
        }
    }
    Ok(())
}

impl Filter {
    /// Tests a record as `filter_test(filter, rec, NULL)` does: whether the
    /// expression holds at the site. `header` is the reader's current
    /// header (records may have added contigs and keys to it). Rust-evaluated
    /// regex matches share 1,000,000 work units across search starts, with at
    /// most 128 matcher frames. Exhaustion returns Unsupported. Larger or
    /// branching/repeating back-reference patterns use system REG_NOSUB with
    /// checked input/program sizes, without a hard native execution-time bound.
    /// Literal matching uses a separate bounded-stack, linear-time path.
    pub fn test(&mut self, header: &Header, rec: &Record) -> Result<bool> {
        take_undefined();
        self.cached_gt.valid = false;
        self.cached_site.valid = false;
        let r = self.run(header, rec);
        if let Some(what) = take_undefined() {
            return Err(unsupported(what));
        }
        r
    }

    fn shown(&self) -> String {
        String::from_utf8_lossy(&self.str).into_owned()
    }

    fn run(&mut self, h: &Header, rec: &Record) -> Result<bool> {
        let n = self.filters.len();
        let mut nstack: usize = 0;
        for i in 0..n {
            self.filters[i].pass_site = 0;
            if self.filters[i].tok_type == TOK_VAL {
                if let Some(s) = self.filters[i].setter {
                    self.set(s, i, h, rec)?;
                } else {
                    let t = &mut self.filters[i];
                    if let Some(key) = &t.key {
                        t.str_value.l = 0;
                        t.str_value.puts(cstr(key));
                        t.nvalues = t.str_value.l as i32;
                    } else if t.is_constant {
                        t.values.set(0, t.threshold);
                        t.nvalues = 1;
                    }
                }
                self.stack[nstack] = i;
                nstack += 1;
                continue;
            }
            if let Some(f) = self.filters[i].func {
                let nargs = self.func(f, i, h, rec, nstack)?;
                let at = nstack as isize - nargs as isize;
                if at < 0 {
                    return Err(unsupported(
                        "filter function arguments missing from the stack",
                    ));
                }
                self.stack[at as usize] = i;
                if nargs - 1 > 0 {
                    nstack -= (nargs - 1) as usize;
                }
                continue;
            }
            if nstack < 2 {
                return Err(run_error(format!(
                    "Error occurred while processing the filter \"{}\" (1:{nstack})",
                    self.shown()
                )));
            }
            let (a, b) = (self.stack[nstack - 2], self.stack[nstack - 1]);
            let is_str = i32::from(self.filters[b].is_str) + i32::from(self.filters[a].is_str);
            let op = self.filters[i].tok_type;
            if matches!(op, TOK_ADD | TOK_SUB | TOK_MULT | TOK_DIV | TOK_MODULO) {
                let (left, r) = split(&mut self.filters, i);
                arith(&left[a], &left[b], r, op)?;
                self.stack[nstack - 2] = i;
                nstack -= 1;
                continue;
            }
            if let Some(c) = self.filters[i].comparator {
                self.compare_with(c, b, a, i, h, rec)?;
            } else if let Some(c) = self.filters[b].comparator {
                self.compare_with(c, b, a, i, h, rec)?;
            } else if let Some(c) = self.filters[a].comparator {
                self.compare_with(c, a, b, i, h, rec)?;
            } else if is_str == 2 {
                let (left, r) = split(&mut self.filters, i);
                cmp_vector_strings(&left[a], &left[b], r)?;
            } else {
                if is_str == 1 {
                    return Err(run_error(
                        "Error: cannot use arithmetic operators to compare strings and numbers",
                    ));
                }
                let mut ml = [0, 0, 0];
                if op == TOK_EQ {
                    ml = [1, 0, 1];
                }
                if op == TOK_NE {
                    ml = [1, 1, 0];
                }
                if !matches!(op, TOK_EQ | TOK_NE | TOK_LE | TOK_LT | TOK_BT | TOK_BE) {
                    return Err(run_error(format!("todo: filter.c .. type={op}")));
                }
                let (left, r) = split(&mut self.filters, i);
                cmp_vectors(&left[a], &left[b], r, op, &ml)?;
            }
            self.stack[nstack - 2] = i;
            nstack -= 1;
        }
        if nstack > 1 {
            return Err(run_error(format!(
                "Error occurred while processing the filter \"{}\" (2:{nstack})",
                self.shown()
            )));
        }
        if nstack == 0 {
            return Err(unsupported("filter result read from an empty stack"));
        }
        Ok(self.filters[self.stack[0]].pass_site != 0)
    }

    fn compare_with(
        &mut self,
        c: Cmp,
        a: usize,
        b: usize,
        r: usize,
        h: &Header,
        rec: &Record,
    ) -> Result<()> {
        match c {
            Cmp::BitAnd => {
                let (left, rt) = split(&mut self.filters, r);
                let v = |t: &Token| {
                    c_int(if t.nvalues != 0 {
                        t.values.get(0)
                    } else {
                        t.threshold
                    })
                };
                let x = v(&left[a]) & v(&left[b]);
                rt.pass_site = if rt.tok_type == TOK_LIKE {
                    i32::from(x != 0)
                } else {
                    i32::from(x == 0)
                };
                Ok(())
            }
            Cmp::Filter => {
                let (left, rt) = split(&mut self.filters, r);
                cmp_filter(&left[b], rt, rec)
            }
            Cmp::StringHash => {
                let (left, rt) = split(&mut self.filters, r);
                let (mut a, mut b) = (&left[a], &left[b]);
                if b.hash.is_some() {
                    std::mem::swap(&mut a, &mut b);
                }
                cmp_string_hash(a, b, rt)
            }
            Cmp::Id => {
                let (left, rt) = split(&mut self.filters, r);
                let (mut a, mut b) = (&left[a], &left[b]);
                if b.hash.is_some() {
                    std::mem::swap(&mut a, &mut b);
                }
                let _ = h;
                cmp_id(a, b, rt, rec)
            }
        }
    }
}

/// `filters_cmp_filter`: `b` holds the FILTER ids of the expression.
fn cmp_filter(b: &Token, r: &mut Token, rec: &Record) -> Result<()> {
    let flt = rec.filters();
    let has = |id: f64| flt.iter().any(|&f| f64::from(f) == id);
    let nb = b.nvalues.max(0) as isize;
    let all_present = || (0..nb).all(|j| has(b.values.get(j)));
    match r.tok_type {
        TOK_NOT_IN => {
            if nb == 0 {
                if !flt.is_empty() {
                    r.pass_site = 1;
                }
                return Ok(());
            }
            if flt.is_empty() || !all_present() {
                r.pass_site = 1;
            }
        }
        TOK_IN => {
            if nb == 0 {
                if flt.is_empty() {
                    r.pass_site = 1;
                }
                return Ok(());
            }
            if !flt.is_empty() && all_present() {
                r.pass_site = 1;
            }
        }
        TOK_NE => {
            if b.nvalues as usize != flt.len() {
                r.pass_site = 1;
                return Ok(());
            }
            if nb == 0 {
                return Ok(());
            }
            if !all_present() {
                r.pass_site = 1;
            }
        }
        TOK_EQ => {
            if b.nvalues as usize != flt.len() {
                return Ok(());
            }
            if nb == 0 || all_present() {
                r.pass_site = 1;
            }
        }
        _ => {
            return Err(run_error(
                "Only ==, !=, ~, and !~ operators are supported for FILTER",
            ));
        }
    }
    Ok(())
}

/// `filters_cmp_string_hash`: `a` holds the `@file` names.
fn cmp_string_hash(a: &Token, b: &Token, r: &mut Token) -> Result<()> {
    if r.tok_type != TOK_EQ && r.tok_type != TOK_NE {
        return Err(run_error(
            "Only == and != operators are supported for strings read from a file",
        ));
    }
    if !b.is_str {
        return Err(run_error(
            "Error: Lookup of values from a file is currenlty supported only for strings",
        ));
    }
    let hash = a.hash.as_ref().unwrap();
    let ne = r.tok_type == TOK_NE;
    let s = &b.str_value;
    if b.nsamples == 0 {
        if b.idx >= 0 {
            let mut ret = !s.is_null() && hash.contains(&s.cstr(0));
            if ne {
                ret = !ret;
            }
            r.pass_site = i32::from(ret);
            return Ok(());
        }
        let mut ret = false;
        let mut ptr = 0;
        while s.get(ptr) != 0 {
            let mut eptr = ptr + 1;
            while s.get(eptr) != 0 && s.get(eptr) != b',' {
                eptr += 1;
            }
            let field: Vec<u8> = (ptr..eptr).map(|k| s.get(k)).collect();
            ret |= hash.contains(&field);
            if s.get(eptr) == 0 {
                break;
            }
            ptr = eptr + 1;
        }
        if ne {
            ret = !ret;
        }
        r.pass_site = i32::from(ret);
        return Ok(());
    }
    init_samples(a, b, r);
    r.pass_site = 0;
    let nval1 = b.nval1.max(0) as usize;
    for i in 0..b.nsamples as usize {
        if !r.usmpl(i) {
            continue;
        }
        let start = i * nval1;
        // The block is read as a C string ending at its width.
        let byte = |k: usize| if k >= start + nval1 { 0 } else { s.get(k) };
        let mut ret = false;
        if b.idx >= 0 {
            let mut k = start;
            let mut field = Vec::new();
            while byte(k) != 0 {
                field.push(byte(k));
                k += 1;
            }
            ret = hash.contains(&field);
        } else {
            let mut ptr = start;
            while byte(ptr) != 0 {
                let mut eptr = ptr + 1;
                while byte(eptr) != 0 && byte(eptr) != b',' {
                    eptr += 1;
                }
                let field: Vec<u8> = (ptr..eptr).map(byte).collect();
                ret |= hash.contains(&field);
                if byte(eptr) == 0 {
                    break;
                }
                ptr = eptr + 1;
            }
        }
        if ne {
            ret = !ret;
        }
        r.set_pass_sample(i, u8::from(ret));
        r.pass_site |= i32::from(ret);
    }
    Ok(())
}

/// `filters_cmp_id`: ID against a string, a regex or an `@file`.
fn cmp_id(a: &Token, b: &Token, r: &mut Token, rec: &Record) -> Result<()> {
    let ids = rec.id();
    let mut pass;
    let mut start = 0;
    loop {
        let end = ids[start..]
            .iter()
            .position(|&c| c == b';')
            .map(|p| start + p);
        let id = &ids[start..end.unwrap_or(ids.len())];
        if let Some(hash) = &a.hash {
            if r.tok_type != TOK_EQ && r.tok_type != TOK_NE {
                return Err(run_error(
                    "Only == and != operators are supported for strings read from a file",
                ));
            }
            pass = hash.contains(id);
        } else {
            if b.str_value.l == 0 {
                return Err(run_error("Error occurred while evaluating the expression"));
            }
            if r.tok_type == TOK_EQ || r.tok_type == TOK_NE {
                pass = b.str_value.cstr(0) == id;
            } else {
                if r.tok_type != TOK_LIKE && r.tok_type != TOK_NLIKE {
                    return Err(run_error(format!(
                        "Only the following operators are supported for querying ID: ==, !=, ~, !~; the operator type {} is not supported",
                        r.tok_type
                    )));
                }
                let re = a
                    .regex
                    .as_deref()
                    .or(b.regex.as_deref())
                    .ok_or_else(|| run_error("fixme: regex initialization failed"))?;
                pass = re.is_match(id);
            }
        }
        match end {
            Some(e) if !pass => start = e + 1,
            _ => break,
        }
    }
    if r.tok_type == TOK_NE {
        pass = !pass;
    }
    r.pass_site = i32::from(pass);
    Ok(())
}

impl Filter {
    /// `filters_cache_genotypes`: 0, or -1 without usable genotypes.
    fn cache_genotypes(&mut self, h: &Header, rec: &Record) -> i32 {
        let c = &mut self.cached_gt;
        if c.valid {
            return if c.nbuf > 0 { 0 } else { -1 };
        }
        c.valid = true;
        c.nbuf = match rec.genotype_values(h, &mut c.buf) {
            Ok(n) => n as i32,
            Err(e) => e.code(),
        };
        if c.nbuf <= 0 {
            return -1;
        }
        let ns = self.nsamples.max(0) as usize;
        let mask = c.mask.get_or_insert_with(|| vec![0; ns]);
        let n_sample = rec.n_sample();
        let ngt1 = c.nbuf as usize / n_sample;
        for i in 0..n_sample {
            if i >= mask.len() {
                undefined("genotype mask written past its end");
                break;
            }
            mask[i] = 0;
            for j in 0..ngt1 {
                let v = i64::from(c.buf[i * ngt1 + j]);
                if gt_missing(v) {
                    continue;
                }
                if v == i64::from(INT32_VECTOR_END) {
                    break;
                }
                let allele = gt_allele(v);
                if allele > 63 {
                    c.nbuf = 0;
                    return -1;
                }
                mask[i] |= bit(allele);
            }
        }
        0
    }

    fn mask(&self, i: usize) -> u64 {
        match self.cached_gt.mask.as_ref().and_then(|m| m.get(i)) {
            Some(&m) => m,
            None => {
                undefined("genotype mask read past its end");
                0
            }
        }
    }

    /// `filters_cache_AC_stats`: 0, or -1 when the counts are unusable.
    fn cache_ac(&mut self, h: &Header, rec: &Record) -> Result<i32> {
        if self.cached_site.valid {
            return Ok(if self.cached_site.n_samples > 0 {
                0
            } else {
                -1
            });
        }
        self.cached_site.valid = true;
        let n_allele = rec.n_alleles();
        let c = &mut self.cached_site;
        if c.ac.len() < n_allele {
            c.ac.resize(n_allele, 0);
            c.af.resize(n_allele, 0.0);
        }
        if calc_ac_info(h, rec, &mut c.ac)? {
            c.n_samples = rec.n_sample() as i32;
            c.n_missing = -1;
            c.an = c.ac[..n_allele]
                .iter()
                .fold(0i32, |s, &v| s.wrapping_add(v));
            for i in 0..n_allele {
                c.af[i] = if c.an != 0 {
                    f64::from(c.ac[i]) / f64::from(c.an)
                } else {
                    0.0
                };
            }
            return Ok(0);
        }
        let n_sample = rec.n_sample();
        c.n_samples = n_sample as i32;
        c.n_missing = 0;
        c.an = 0;
        c.ac[..n_allele].fill(0);
        let mut gt = Vec::new();
        let ngt = match rec.genotype_values(h, &mut gt) {
            Ok(n) => {
                fill_ints(&mut self.tmpi, n, &gt);
                n as i32
            }
            Err(e) => e.code(),
        };
        if n_sample == 0 {
            return Err(unsupported("integer division by zero in filter.c"));
        }
        let ngt1 = ngt / n_sample as i32;
        let c = &mut self.cached_site;
        for i in 0..n_sample as i32 {
            let mut is_missing = true;
            for j in 0..ngt1 {
                let v = i64::from(read_int(&self.tmpi, (ngt1 * i + j) as isize));
                if gt_missing(v) {
                    continue;
                }
                if v == i64::from(INT32_VECTOR_END) {
                    break;
                }
                let ial = gt_allele(v);
                if ial >= n_allele as i64 {
                    c.n_samples = 0;
                    return Ok(-1);
                }
                if ial < 0 {
                    return Err(unsupported("allele count written before its array"));
                }
                c.ac[ial as usize] += 1;
                c.an += 1;
                is_missing = false;
            }
            if is_missing {
                c.n_missing += 1;
            }
        }
        for i in 0..n_allele {
            c.af[i] = if c.an != 0 {
                f64::from(c.ac[i]) / f64::from(c.an)
            } else {
                0.0
            };
        }
        Ok(0)
    }
}

/// `bcf_calc_ac(hdr, line, ac, BCF_UN_INFO)`: allele counts from INFO/AN
/// and INFO/AC; `false` when they are absent or do not fit.
fn calc_ac_info(h: &Header, rec: &Record, ac: &mut [i32]) -> Result<bool> {
    let n_allele = rec.n_alleles();
    ac[..n_allele].fill(0);
    let (Some(an_id), Some(ac_id)) = (h.key_id("AN"), h.key_id("AC")) else {
        return Ok(false);
    };
    let mut an: i32 = -1;
    let mut acf = None;
    for f in rec.info_fields() {
        if f.key == an_id {
            an = f.v1_int() as i32;
        } else if f.key == ac_id {
            acf = Some(f);
        }
    }
    let Some(acf) = acf else { return Ok(false) };
    if an < 0 {
        return Ok(false);
    }
    let len = acf.values.len();
    if len != n_allele - 1 {
        return Ok(false);
    }
    if !matches!(acf.values.ty, BT_INT8 | BT_INT16 | BT_INT32) {
        return Err(unsupported("htslib exits on INFO/AC of this type"));
    }
    let mut nac: i32 = 0;
    for i in 0..len {
        let v = acf.values.raw(i) as i32;
        ac[i + 1] = v;
        nac = nac.wrapping_add(v);
    }
    if an < nac {
        return Err(unsupported(
            "htslib exits on INFO/AN smaller than the AC sum",
        ));
    }
    ac[0] = an - nac;
    Ok(true)
}

impl Filter {
    fn set(&mut self, s: Setter, i: usize, h: &Header, rec: &Record) -> Result<()> {
        let n_allele = rec.n_alleles();
        match s {
            Setter::Qual => {
                let t = &mut self.filters[i];
                let q = rec.qual_bits();
                t.values.set(
                    0,
                    if q == FLOAT_MISSING {
                        missing()
                    } else {
                        f64::from(f32::from_bits(q))
                    },
                );
                t.nvalues = 1;
            }
            Setter::Type => {
                let t = &mut self.filters[i];
                let v = rec.variant_types();
                t.values.set(0, f64::from(if v == 0 { 1 } else { v << 1 }));
                t.nvalues = 1;
            }
            Setter::Chrom => {
                let t = &mut self.filters[i];
                t.str_value.l = 0;
                let name = h
                    .contigs
                    .get(rec.rid as usize)
                    .map(|c| c.name.as_bytes())
                    .unwrap_or_default();
                t.str_value.puts(cstr(name));
                t.nvalues = t.str_value.l as i32;
                t.is_str = true;
            }
            Setter::Pos => {
                let pos = one_based_pos(rec)?;
                let t = &mut self.filters[i];
                t.values.set(0, pos as f64);
                t.nvalues = 1;
            }
            Setter::Info => self.set_info(i, rec),
            Setter::InfoInt | Setter::InfoFloat => {
                self.set_info_number(i, h, rec, s == Setter::InfoInt)?
            }
            Setter::InfoString => self.set_info_string(i, h, rec),
            Setter::InfoFlag => {
                let t = &mut self.filters[i];
                let present = t.hdr_id >= 0 && rec.info_field(t.hdr_id as u32).is_some();
                t.values.set(0, if present { 1.0 } else { 0.0 });
                t.nvalues = 1;
            }
            Setter::FormatInt | Setter::FormatFloat => {
                self.set_format_number(i, h, rec, s == Setter::FormatInt)?
            }
            Setter::FormatString => self.set_format_string(i, h, rec)?,
            Setter::GenotypeString => self.set_genotype_string(i, h, rec)?,
            Setter::Genotype2 | Setter::Genotype3 | Setter::Genotype4 => {
                let ty = match s {
                    Setter::Genotype2 => 2,
                    Setter::Genotype3 => 3,
                    _ => 4,
                };
                self.set_genotype(i, h, rec, ty)?
            }
            Setter::Ilen => {
                if n_allele == 0 {
                    return Err(unsupported("filter.c reads a missing REF allele"));
                }
                let t = &mut self.filters[i];
                t.nvalues = n_allele as i32 - 1;
                t.values.expand(t.nvalues.max(0) as usize);
                let rlen = cstr(rec.allele(0)).len() as i64;
                for k in 1..n_allele {
                    let a = cstr(rec.allele(k));
                    if a.first() == Some(&b'<') {
                        t.values.set(k as isize - 1, missing());
                        continue;
                    }
                    t.values.set(k as isize - 1, (a.len() as i64 - rlen) as f64);
                }
            }
            Setter::Ref => {
                if n_allele == 0 {
                    return Err(unsupported("filter.c reads a missing REF allele"));
                }
                let t = &mut self.filters[i];
                t.str_value.l = 0;
                t.str_value.puts(cstr(rec.allele(0)));
                t.nvalues = t.str_value.l as i32;
            }
            Setter::Alt => {
                let t = &mut self.filters[i];
                t.str_value.l = 0;
                let na = n_allele as i32;
                if t.idx >= 0 {
                    let allele = t
                        .idx
                        .checked_add(1)
                        .ok_or_else(|| unsupported("filter ALT index overflow"))?;
                    if na > allele {
                        t.str_value.puts(cstr(rec.allele(allele as usize)));
                    } else {
                        t.str_value.putc(b'.');
                    }
                } else if t.idx == -2 {
                    let mut end = if t.nuidxs != 0 { t.nuidxs } else { na - 1 };
                    if end >= na - 1 {
                        end = na - 2;
                    }
                    for k in 0..=end {
                        if k >= t.nidxs() || t.idxs[k as usize] != 0 {
                            if t.str_value.l > 0 {
                                t.str_value.putc(b',');
                            }
                            t.str_value.puts(cstr(rec.allele(k as usize + 1)));
                        }
                    }
                } else if na > 1 {
                    t.str_value.puts(cstr(rec.allele(1)));
                    for k in 2..n_allele {
                        t.str_value.putc(b',');
                        t.str_value.puts(cstr(rec.allele(k)));
                    }
                } else if na == 1 {
                    t.str_value.putc(b'.');
                }
                t.nvalues = t.str_value.l as i32;
            }
            Setter::Nalt => {
                let t = &mut self.filters[i];
                t.nvalues = 1;
                t.values.set(0, (n_allele as i64 - 1) as f64);
            }
            Setter::Nmissing => self.set_nmissing(i, h, rec)?,
            Setter::Ac | Setter::Af => {
                let ok = self.cache_ac(h, rec)? == 0 && self.cached_site.an != 0;
                let c = &self.cached_site;
                let t = &mut self.filters[i];
                if !ok {
                    t.nvalues = 0;
                    return Ok(());
                }
                let get = |k: usize| {
                    if s == Setter::Ac {
                        f64::from(c.ac[k])
                    } else {
                        c.af[k]
                    }
                };
                if t.idx >= 0 {
                    t.nvalues = 1;
                    let k = t.idx as usize + 1;
                    t.values.set(0, if k < n_allele { get(k) } else { 0.0 });
                } else if n_allele == 1 {
                    t.nvalues = 1;
                    t.values.set(0, 0.0);
                } else {
                    t.values.expand(n_allele);
                    for k in 1..n_allele {
                        t.values.set(k as isize - 1, get(k));
                    }
                    t.nvalues = n_allele as i32 - 1;
                }
            }
            Setter::An => {
                let ok = self.cache_ac(h, rec)? == 0;
                let t = &mut self.filters[i];
                t.nvalues = 1;
                t.values.set(
                    0,
                    if ok {
                        f64::from(self.cached_site.an)
                    } else {
                        0.0
                    },
                );
            }
            Setter::Mac => {
                let ok = self.cache_ac(h, rec)? == 0 && self.cached_site.an != 0;
                let c = &self.cached_site;
                let t = &mut self.filters[i];
                if !ok {
                    t.nvalues = 0;
                    return Ok(());
                }
                t.nvalues = 1;
                let mut min = i32::MAX;
                for &v in &c.ac[..n_allele] {
                    if v != 0 && v < min {
                        min = v;
                    }
                }
                if min == 0 {
                    t.nvalues = 0;
                } else if min == c.an {
                    t.values.set(0, 0.0);
                } else {
                    t.values.set(0, f64::from(min));
                }
            }
            Setter::Maf => {
                let ok = self.cache_ac(h, rec)? == 0 && self.cached_site.an != 0;
                let c = &self.cached_site;
                let t = &mut self.filters[i];
                if !ok {
                    t.nvalues = 0;
                    return Ok(());
                }
                t.nvalues = 1;
                let mut min = f64::INFINITY;
                for &v in &c.af[..n_allele] {
                    if v != 0.0 && v < min {
                        min = v;
                    }
                }
                if min == 0.0 {
                    t.nvalues = 0;
                } else if min == 1.0 {
                    t.values.set(0, 0.0);
                } else {
                    t.values.set(0, min);
                }
            }
        }
        Ok(())
    }

    /// `filters_set_info` (Number=1 keys).
    fn set_info(&mut self, i: usize, rec: &Record) {
        let t = &mut self.filters[i];
        let f = (t.hdr_id >= 0)
            .then(|| rec.info_field(t.hdr_id as u32))
            .flatten();
        let Some(f) = f else {
            t.nvalues = 0;
            t.str_value.l = 0;
            return;
        };
        match f.values.ty {
            BT_CHAR => {
                let n = f.values.len();
                if n >= t.str_value.m {
                    t.str_value.grow_exact(n + 1);
                }
                for (k, &c) in f.values.bytes[..n].iter().enumerate() {
                    t.str_value.set(k, c);
                }
                t.str_value.set(n, 0);
                t.str_value.l = n;
                t.nvalues = n as i32;
            }
            BT_FLOAT => {
                let v = f.v1_float();
                if v.to_bits() == FLOAT_MISSING {
                    t.nvalues = 0;
                } else {
                    t.values.set(0, f64::from(v));
                    t.nvalues = 1;
                }
                t.str_value.l = 0;
            }
            ty => {
                t.str_value.l = 0;
                let v1 = f.v1_int();
                let miss = match ty {
                    BT_INT8 => v1 == i8::MIN as i64,
                    BT_INT16 => v1 == i16::MIN as i64,
                    BT_INT32 => v1 == i32::MIN as i64,
                    _ => false,
                };
                if miss {
                    t.nvalues = 0;
                } else {
                    t.values.set(0, v1 as f64);
                    t.nvalues = 1;
                }
            }
        }
    }

    /// `filters_set_info_int` / `filters_set_info_float`.
    fn set_info_number(&mut self, i: usize, h: &Header, rec: &Record, int: bool) -> Result<()> {
        let (tag, idx, hdr_id) = {
            let t = &self.filters[i];
            (tag_str(t).to_owned(), t.idx, t.hdr_id)
        };
        if idx == -2 {
            let stored = (hdr_id >= 0)
                .then(|| rec.info_field(hdr_id as u32))
                .flatten()
                .map_or(0, |f| f.values.len());
            let n = if int {
                let mut v = Vec::new();
                match rec.info_int32(h, &tag, &mut v) {
                    Ok(n) => {
                        fill_ints(&mut self.tmpi, stored, &v);
                        n as i32
                    }
                    Err(e) => e.code(),
                }
            } else {
                let mut v = Vec::new();
                match rec.info_float(h, &tag, &mut v) {
                    Ok(n) => {
                        fill_floats(&mut self.tmpf, stored, &v);
                        n as i32
                    }
                    Err(e) => e.code(),
                }
            };
            let t = &mut self.filters[i];
            if n <= 0 {
                t.nvalues = 0;
                return Ok(());
            }
            t.values.expand(n as usize);
            let nidxs = t.nidxs();
            let mut end = if *t.idxs.last().unwrap_or(&0) < 0 {
                n - 1
            } else {
                nidxs - 1
            };
            if end >= n {
                end = n - 1;
            }
            let mut j = 0;
            for k in 0..=end {
                if k >= nidxs || t.idxs[k as usize] != 0 {
                    let v = if int {
                        f64::from(read_int(&self.tmpi, k as isize))
                    } else {
                        let b = read_float(&self.tmpf, k as isize);
                        if b == FLOAT_MISSING {
                            missing()
                        } else {
                            f64::from(f32::from_bits(b))
                        }
                    };
                    t.values.set(j, v);
                    j += 1;
                }
            }
            t.nvalues = j as i32;
            return Ok(());
        }
        let v = info_value(rec, hdr_id, idx)?;
        let t = &mut self.filters[i];
        let value = match v {
            InfoValue::Absent | InfoValue::Missing => {
                t.nvalues = 0;
                return Ok(());
            }
            InfoValue::Int(x) if int => x as f64,
            InfoValue::Float(x) if !int => x,
            // A value of the other kind lands in the variable bit for bit.
            InfoValue::Int(x) => f64::from_bits(x as u64),
            InfoValue::Float(x) => x.to_bits() as i64 as f64,
            InfoValue::Untouched if int => 0.0,
            InfoValue::Untouched => {
                return Err(unsupported("read of an uninitialised INFO value"));
            }
        };
        t.values.set(0, value);
        t.nvalues = 1;
        Ok(())
    }

    /// `filters_set_info_string`.
    fn set_info_string(&mut self, i: usize, h: &Header, rec: &Record) {
        let Filter { filters, tmps, .. } = self;
        let t = &mut filters[i];
        let tag = tag_str(t).to_owned();
        let mut bytes = Vec::new();
        let n = match rec.info_string(h, &tag, &mut bytes) {
            Ok(n) => {
                if t.str_value.m < n + 1 {
                    t.str_value.grow_exact(n + 1);
                }
                for (k, &c) in bytes.iter().enumerate() {
                    t.str_value.set(k, c);
                }
                t.str_value.set(n, 0);
                n
            }
            Err(_) => {
                t.nvalues = 0;
                t.str_value.l = 0;
                return;
            }
        };
        let s = &mut t.str_value;
        if t.idx >= 0 {
            let idx = t.idx as usize;
            let (mut k, mut ss) = (0, 0);
            while ss < n && k < idx {
                if s.get(ss) == b',' {
                    k += 1;
                }
                ss += 1;
            }
            if ss == n || k != idx {
                t.nvalues = 0;
                s.l = 0;
                return;
            }
            let mut se = ss;
            while se < n && s.get(se) != b',' {
                se += 1;
            }
            if ss == 0 {
                s.set(se, 0);
            } else {
                for k in 0..se - ss {
                    let c = s.get(ss + k);
                    s.set(k, c);
                }
                s.set(se - ss, 0);
            }
            s.l = se - ss;
        } else if t.idx == -2 && t.idxs.first() == Some(&-1) {
            s.l = n;
        } else if t.idx == -2 {
            tmps.l = 0;
            tmps.resize(n);
            let nidxs = t.idxs.len() as i32;
            let mut iend = if *t.idxs.last().unwrap_or(&0) < 0 {
                n as i32 - 1
            } else {
                nidxs - 1
            };
            if iend >= n as i32 {
                iend = n as i32 - 1;
            }
            let (mut beg, mut dst) = (0, 0);
            for k in 0..=iend {
                let mut end = beg;
                while s.get(end) != 0 && s.get(end) != b',' {
                    end += 1;
                }
                if k >= nidxs || t.idxs[k as usize] != 0 {
                    for q in beg..end {
                        let c = s.get(q);
                        tmps.set(dst, c);
                        dst += 1;
                    }
                    tmps.set(dst, b',');
                    dst += 1;
                }
                beg = end + 1;
            }
            tmps.set(dst, 0);
            let l = dst;
            std::mem::swap(&mut tmps.buf, &mut s.buf);
            std::mem::swap(&mut tmps.m, &mut s.m);
            std::mem::swap(&mut tmps.sparse, &mut s.sparse);
            s.l = l;
        }
        t.nvalues = t.str_value.l as i32;
    }

    fn check_nsample(&self, i: usize, h: &Header, rec: &Record) -> Result<()> {
        let t = &self.filters[i];
        if rec.n_sample() as i32 != t.nsamples {
            let chrom = h
                .contigs
                .get(rec.rid as usize)
                .map_or("", |c| c.name.as_str());
            return Err(run_error(format!(
                "Incorrect number of FORMAT fields at {chrom}:{} .. {}, {} vs {}",
                one_based_pos(rec)?,
                tag_str(t),
                rec.n_sample(),
                t.nsamples
            )));
        }
        Ok(())
    }

    /// `filters_set_format_int` / `filters_set_format_float`.
    fn set_format_number(&mut self, i: usize, h: &Header, rec: &Record, int: bool) -> Result<()> {
        self.check_nsample(i, h, rec)?;
        let tag = tag_str(&self.filters[i]).to_owned();
        let nvals = if int {
            let mut v = Vec::new();
            match rec.format_int32(h, &tag, &mut v) {
                Ok(n) => {
                    fill_ints(&mut self.tmpi, n, &v);
                    n as i32
                }
                Err(crate::fields::GetError::Unsupported) => {
                    return Err(unsupported("htslib exits on a FORMAT value of this type"));
                }
                Err(e) => e.code(),
            }
        } else {
            let mut v = Vec::new();
            match rec.format_float(h, &tag, &mut v) {
                Ok(n) => {
                    fill_floats(&mut self.tmpf, n, &v);
                    n as i32
                }
                Err(crate::fields::GetError::Unsupported) => {
                    return Err(unsupported("htslib exits on a FORMAT value of this type"));
                }
                Err(e) => e.code(),
            }
        };
        if nvals < 0 {
            self.filters[i].nvalues = 0;
            return Ok(());
        }
        let nsamples = self.filters[i].nsamples;
        if nsamples == 0 {
            return Err(unsupported("integer division by zero in filter.c"));
        }
        let nsrc1 = nvals / nsamples;
        let idx = self.filters[i].idx;
        // The value at `k` of the shared buffer, as a double.
        let read = |f: &Filter, k: isize, raw: bool| -> f64 {
            if int {
                let v = read_int(&f.tmpi, k);
                if raw {
                    f64::from(v)
                } else if v == INT32_MISSING {
                    missing()
                } else if v == INT32_VECTOR_END {
                    vector_end()
                } else {
                    f64::from(v)
                }
            } else {
                let b = read_float(&f.tmpf, k);
                if b == FLOAT_MISSING {
                    missing()
                } else if b == FLOAT_VECTOR_END {
                    vector_end()
                } else {
                    f64::from(f32::from_bits(b))
                }
            }
        };
        let (nval1, nidxs, idxs, usmpl) = {
            let t = &mut self.filters[i];
            t.nval1 = if t.idx >= 0 {
                1
            } else if t.nuidxs != 0 {
                t.nuidxs
            } else {
                nsrc1
            };
            t.nvalues = t.nval1 * nsamples;
            t.values.expand(t.nvalues.max(0) as usize);
            (
                t.nval1,
                t.nidxs(),
                t.idxs.clone(),
                (0..nsamples as usize)
                    .map(|s| t.usmpl(s))
                    .collect::<Vec<_>>(),
            )
        };
        if idx == -3 && self.cache_genotypes(h, rec) != 0 {
            self.filters[i].nvalues = 0;
            return Ok(());
        }
        let mut out: Vec<(isize, f64)> = Vec::new();
        for s in 0..nsamples as isize {
            if !usmpl[s as usize] {
                continue;
            }
            let src = s * nsrc1 as isize;
            if idx >= 0 {
                let v = if idx >= nsrc1 {
                    missing()
                } else {
                    read(self, src + idx as isize, false)
                };
                out.push((s, v));
            } else if idx == -3 {
                let dst = s * nval1 as isize;
                let m = self.mask(s as usize);
                let mut j = 0;
                for k in 0..nsrc1 as i64 {
                    if m & bit(k) == 0 {
                        continue;
                    }
                    // Integers are copied without mapping their sentinels.
                    out.push((dst + j, read(self, src + k as isize, true)));
                    j += 1;
                }
                if j == 0 {
                    out.push((dst, missing()));
                    j = 1;
                }
                while j < nval1 as isize {
                    out.push((dst + j, vector_end()));
                    j += 1;
                }
            } else {
                let dst = s * nval1 as isize;
                let kend = if *idxs.last().unwrap_or(&0) < 0 {
                    nval1
                } else {
                    nidxs
                };
                let mut j = 0;
                for k in 0..kend {
                    if k < nidxs && idxs[k as usize] == 0 {
                        continue;
                    }
                    out.push((dst + j, read(self, src + k as isize, false)));
                    j += 1;
                }
                if j == 0 {
                    out.push((dst, missing()));
                    j = 1;
                }
                while j < nval1 as isize {
                    out.push((dst + j, vector_end()));
                    j += 1;
                }
            }
        }
        let t = &mut self.filters[i];
        for (k, v) in out {
            t.values.set(k, v);
        }
        Ok(())
    }

    /// `filters_set_format_string`.
    fn set_format_string(&mut self, i: usize, h: &Header, rec: &Record) -> Result<()> {
        self.check_nsample(i, h, rec)?;
        let tag = tag_str(&self.filters[i]).to_owned();
        let mut bytes = Vec::new();
        let got = rec.format_char(h, &tag, &mut bytes);
        {
            let t = &mut self.filters[i];
            let s = &mut t.str_value;
            // GNU C narrows size_t capacity to the getter's signed int.
            let mut ndim = s.m as i32;
            if let Ok(n) = got {
                let n = i32::try_from(n)
                    .map_err(|_| unsupported("filter FORMAT string length overflow"))?;
                if ndim < n {
                    if n == 0 {
                        return Err(unsupported(
                            "filter.c appends after a zero-size realloc failure",
                        ));
                    }
                    s.reallocate(n as usize);
                    ndim = n;
                }
                for (k, &c) in bytes.iter().enumerate() {
                    s.set(k, c);
                }
            }
            if ndim < 0 {
                return Err(unsupported(
                    "filter.c appends using a negative string capacity",
                ));
            }
            s.m = ndim as usize;
            s.l = ndim as usize;
            s.putc(0);
            s.l = 0;
            t.nvalues = 0;
        }
        let Ok(nstr) = got else { return Ok(()) };
        if self.filters[i].idx == -3 && self.cache_genotypes(h, rec) != 0 {
            return Ok(());
        }
        let masks: Vec<u64> = (0..self.filters[i].nsamples.max(0) as usize)
            .map(|s| {
                if self.filters[i].idx == -3 {
                    self.mask(s)
                } else {
                    0
                }
            })
            .collect();
        let t = &mut self.filters[i];
        t.nvalues = nstr as i32;
        t.str_value.l = nstr;
        if t.nsamples == 0 {
            return Err(unsupported("integer division by zero in filter.c"));
        }
        let nval1 = nstr / t.nsamples as usize;
        t.nval1 = nval1 as i32;
        for (smp, &mask) in masks.iter().enumerate() {
            if !t.usmpl(smp) {
                continue;
            }
            let s = &mut t.str_value;
            let src = smp * nval1;
            let mut dst = src;
            let mut ibeg = 0;
            let mut idx: i64 = 0;
            while ibeg < nval1 {
                let mut iend = ibeg;
                while iend < nval1 && s.get(src + iend) != 0 && s.get(src + iend) != b',' {
                    iend += 1;
                }
                let keep = if t.idx >= 0 {
                    i64::from(t.idx) == idx
                } else if t.idx == -3 {
                    mask & bit(idx) != 0
                } else if idx < t.idxs.len() as i64 {
                    t.idxs[idx as usize] != 0
                } else {
                    *t.idxs.last().unwrap_or(&0) < 0
                };
                if keep {
                    let n = iend - ibeg + 1;
                    if ibeg != 0 {
                        for k in 0..n {
                            let c = s.get(src + ibeg + k);
                            s.set(dst + k, c);
                        }
                    }
                    dst += n;
                    if t.idx >= 0 {
                        break;
                    }
                }
                if s.get(src + iend) == 0 {
                    break;
                }
                ibeg = iend + 1;
                idx += 1;
            }
            if dst == src {
                s.set(dst, b'.');
                dst += 2;
            }
            if dst - src < nval1 {
                for k in dst - 1..src + nval1 - 1 {
                    s.set(k, 0);
                }
            }
        }
        Ok(())
    }

    /// `filters_set_genotype_string`: each sample's GT as text, in blocks
    /// of a common width.
    fn set_genotype_string(&mut self, i: usize, h: &Header, rec: &Record) -> Result<()> {
        let t = &mut self.filters[i];
        t.nvalues = 0;
        t.str_value.l = 0;
        let Some(fmt) = gt_field(h, rec) else {
            return Ok(());
        };
        let nsmpl = rec.n_sample();
        let mut blen = 4;
        'again: loop {
            t.str_value.l = 0;
            for s in 0..nsmpl {
                let plen = t.str_value.l;
                format_gt(&fmt, s, &mut t.str_value);
                t.str_value.putc_(0);
                if t.str_value.l - plen > blen {
                    blen *= 2;
                    continue 'again;
                }
                let mut p = t.str_value.l - plen;
                while p < blen {
                    t.str_value.putc_(0);
                    p += 1;
                }
            }
            break;
        }
        if t.nsamples as usize != nsmpl {
            return Err(unsupported("assertion failed in filter.c"));
        }
        t.nvalues = t.str_value.l as i32;
        t.nval1 = blen as i32;
        Ok(())
    }

    /// `_filters_set_genotype`: each sample's GT as a class (`rr`/`ra`/...,
    /// `hom`/`het`/`hap`, `mis`/`ref`/`alt`).
    fn set_genotype(&mut self, i: usize, h: &Header, rec: &Record, ty: usize) -> Result<()> {
        let t = &mut self.filters[i];
        t.nvalues = 0;
        t.str_value.l = 0;
        let Some(fmt) = gt_field(h, rec) else {
            return Ok(());
        };
        let nsmpl = h.samples.len();
        let w = if ty == 2 { 3 } else { 4 };
        if t.str_value.m <= w * nsmpl {
            t.str_value.grow_exact(w * nsmpl + 1);
        }
        if !matches!(fmt.values.ty, BT_INT8 | BT_INT16 | BT_INT32) {
            let chrom = h
                .contigs
                .get(rec.rid as usize)
                .map_or("", |c| c.name.as_str());
            return Err(run_error(format!(
                "The GT type is not lineognised: {} at {chrom}:{}",
                fmt.values.ty,
                one_based_pos(rec)?
            )));
        }
        let end = int_vector_end(fmt.values.ty);
        for s in 0..rec.n_sample() {
            let row = fmt.sample(s);
            let (mut is_het, mut has_ref, mut miss) = (false, false, false);
            let mut jal = 0i64;
            let mut j = 0;
            while j < fmt.n {
                let ial = row.raw(j);
                if ial == end {
                    break;
                }
                if gt_missing(ial) {
                    miss = true;
                    break;
                }
                if gt_allele(ial) == 0 {
                    has_ref = true;
                }
                if j > 0 && gt_allele(ial) != gt_allele(jal) {
                    is_het = true;
                }
                jal = ial;
                j += 1;
            }
            let text: &[u8] = if ty == 4 {
                if j == 0 || miss {
                    b"mis\0"
                } else if !has_ref || is_het {
                    b"alt\0"
                } else {
                    b"ref\0"
                }
            } else if j == 0 || miss {
                b".\0"
            } else if ty == 3 {
                if j == 1 {
                    b"hap\0"
                } else if !is_het {
                    b"hom\0"
                } else {
                    b"het\0"
                }
            } else if j == 1 {
                if has_ref { b"r\0" } else { b"a\0" }
            } else if !is_het {
                if has_ref { b"rr\0" } else { b"aa\0" }
            } else if has_ref {
                b"ra\0"
            } else {
                b"aA\0"
            };
            for (k, &c) in text.iter().enumerate() {
                t.str_value.set(w * s + k, c);
            }
        }
        if t.nsamples as usize != nsmpl {
            return Err(unsupported("assertion failed in filter.c"));
        }
        t.str_value.l = w * nsmpl;
        t.nvalues = t.str_value.l as i32;
        t.str_value.set(t.str_value.l, 0);
        t.nval1 = w as i32;
        Ok(())
    }

    /// `filters_set_nmissing`: N_MISSING and F_MISSING.
    fn set_nmissing(&mut self, i: usize, h: &Header, rec: &Record) -> Result<()> {
        let t = &mut self.filters[i];
        if t.usmpl.is_none() {
            return Err(unsupported("assertion failed in filter.c"));
        }
        if rec.n_sample() == 0 {
            t.nvalues = 1;
            t.values.set(0, 0.0);
            return Ok(());
        }
        let fmt = h.key_id("GT").and_then(|id| rec.format_field(id));
        let Some(fmt) = fmt else {
            t.nvalues = 0;
            return Ok(());
        };
        if !matches!(fmt.values.ty, BT_INT8 | BT_INT16 | BT_INT32) {
            return Err(unsupported("filter.c exits on a GT field of this type"));
        }
        let end = int_vector_end(fmt.values.ty);
        let (mut nsmpl, mut nmissing) = (0, 0);
        for s in 0..rec.n_sample() {
            if !t.usmpl(s) {
                continue;
            }
            nsmpl += 1;
            let row = fmt.sample(s);
            for j in 0..fmt.n {
                let v = row.raw(j);
                if v == end {
                    break;
                }
                if gt_missing(v) {
                    nmissing += 1;
                    break;
                }
            }
        }
        t.nsamples = 0;
        t.nvalues = 1;
        let n_missing = t
            .tag
            .as_deref()
            .is_some_and(|tag| tag.first() == Some(&b'N'));
        t.values.set(
            0,
            if n_missing {
                f64::from(nmissing)
            } else if nsmpl != 0 {
                f64::from(nmissing) / f64::from(nsmpl)
            } else {
                0.0
            },
        );
        Ok(())
    }
}

/// How a per-site reduction combines values.
#[derive(Clone, Copy)]
enum Reduce {
    Max,
    Min,
    Avg,
    Sum,
}

impl Filter {
    fn func(&mut self, f: Func, r: usize, h: &Header, rec: &Record, nstack: usize) -> Result<i32> {
        let top = |s: &[usize]| -> Result<usize> {
            if nstack == 0 {
                return Err(unsupported("filter function without an argument"));
            }
            Ok(s[nstack - 1])
        };
        match f {
            Func::LogicOr | Func::LogicAnd => {
                if nstack < 2 {
                    return Err(run_error(format!(
                        "Error occurred while processing the filter \"{}\"",
                        self.shown()
                    )));
                }
                let (a, b) = (self.stack[nstack - 2], self.stack[nstack - 1]);
                let (left, rt) = split(&mut self.filters, r);
                if f == Func::LogicOr {
                    logic_or(&left[a], &left[b], rt);
                } else {
                    logic_and(&left[a], &left[b], rt);
                }
                Ok(2)
            }
            Func::Npass => {
                if nstack == 0 {
                    return Err(run_error("Error parsing the expression"));
                }
                let a = self.stack[nstack - 1];
                let (left, rt) = split(&mut self.filters, r);
                let t = &left[a];
                if t.nsamples == 0 {
                    return Err(run_error(format!(
                        "The function {} works with FORMAT fields",
                        tag_str(rt)
                    )));
                }
                if t.usmpl.is_none() {
                    return Err(unsupported("assertion failed in filter.c"));
                }
                let (mut nsmpl, mut npass) = (0, 0);
                for i in 0..t.nsamples as usize {
                    if !t.usmpl(i) {
                        continue;
                    }
                    nsmpl += 1;
                    if t.pass_sample(i) != 0 {
                        npass += 1;
                    }
                }
                rt.values.expand(1);
                rt.nsamples = 0;
                rt.nvalues = 1;
                let n = rt
                    .tag
                    .as_deref()
                    .is_some_and(|tag| tag.first() == Some(&b'N'));
                rt.values.set(
                    0,
                    if n {
                        f64::from(npass)
                    } else if nsmpl != 0 {
                        f64::from(npass) / f64::from(nsmpl)
                    } else {
                        0.0
                    },
                );
                Ok(1)
            }
            Func::Max | Func::Min | Func::Avg | Func::Sum => {
                let a = top(&self.stack)?;
                let (left, rt) = split(&mut self.filters, r);
                let red = match f {
                    Func::Max => Reduce::Max,
                    Func::Min => Reduce::Min,
                    Func::Avg => Reduce::Avg,
                    _ => Reduce::Sum,
                };
                reduce(&left[a], rt, red);
                Ok(1)
            }
            Func::SMax | Func::SMin | Func::SAvg | Func::SSum => {
                let a = top(&self.stack)?;
                let red = match f {
                    Func::SMax => Reduce::Max,
                    Func::SMin => Reduce::Min,
                    Func::SAvg => Reduce::Avg,
                    _ => Reduce::Sum,
                };
                let (left, rt) = split(&mut self.filters, r);
                let t = &left[a];
                if t.nsamples == 0 {
                    reduce(t, rt, red);
                    return Ok(1);
                }
                smpl_header(t, rt)?;
                for i in 0..t.nsamples as usize {
                    if !rt.usmpl(i) {
                        continue;
                    }
                    let p = (i * t.nval1.max(0) as usize) as isize;
                    let vals: Vec<f64> = (0..t.nval1 as isize)
                        .map(|j| t.values.get(p + j))
                        .filter(|&v| !is_mv(v))
                        .collect();
                    let v = if vals.is_empty() {
                        missing()
                    } else {
                        match red {
                            Reduce::Max => vals
                                .iter()
                                .fold(f64::NEG_INFINITY, |m, &v| if m < v { v } else { m }),
                            Reduce::Min => vals
                                .iter()
                                .fold(f64::INFINITY, |m, &v| if m > v { v } else { m }),
                            Reduce::Avg => vals.iter().fold(0.0, |s, &v| s + v) / vals.len() as f64,
                            Reduce::Sum => vals.iter().fold(0.0, |s, &v| s + v),
                        }
                    };
                    rt.values.set(i as isize, v);
                }
                Ok(1)
            }
            Func::Median | Func::Stddev => {
                let a = top(&self.stack)?;
                let (t, rt) = pair(&mut self.filters, a, r);
                median_or_stddev(t, rt, f == Func::Median);
                Ok(1)
            }
            Func::SMedian | Func::SStddev => {
                let a = top(&self.stack)?;
                let (t, rt) = pair(&mut self.filters, a, r);
                if t.nsamples == 0 {
                    median_or_stddev(t, rt, f == Func::SMedian);
                    return Ok(1);
                }
                smpl_header(t, rt)?;
                let nval1 = t.nval1.max(0) as isize;
                for i in 0..t.nsamples as usize {
                    if !rt.usmpl(i) {
                        continue;
                    }
                    let p = i as isize * nval1;
                    let mut n = 0isize;
                    for j in 0..nval1 {
                        let v = t.values.get(p + j);
                        if is_mv(v) {
                            continue;
                        }
                        if n < j {
                            t.values.set(p + n, v);
                        }
                        n += 1;
                    }
                    let v = if n == 0 {
                        missing()
                    } else if f == Func::SMedian {
                        if n == 1 {
                            t.values.get(p)
                        } else {
                            let mut vals: Vec<f64> = (0..n).map(|k| t.values.get(p + k)).collect();
                            msort(&mut vals);
                            for (k, &v) in vals.iter().enumerate() {
                                t.values.set(p + k as isize, v);
                            }
                            let n = n as usize;
                            if n % 2 == 1 {
                                vals[n / 2]
                            } else {
                                (vals[n / 2 - 1] + vals[n / 2]) * 0.5
                            }
                        }
                    } else if n == 1 {
                        0.0
                    } else {
                        let vals: Vec<f64> = (0..n).map(|k| t.values.get(p + k)).collect();
                        stddev(&vals)
                    };
                    rt.values.set(i as isize, v);
                }
                Ok(1)
            }
            Func::Abs | Func::Phred => {
                let a = top(&self.stack)?;
                let (left, rt) = split(&mut self.filters, r);
                let t = &left[a];
                if f == Func::Abs {
                    func_abs(t, rt)?;
                } else {
                    func_phred(t, rt)?;
                }
                Ok(1)
            }
            Func::Count | Func::SCount => {
                let a = top(&self.stack)?;
                let (left, rt) = split(&mut self.filters, r);
                let t = &left[a];
                if f == Func::SCount && t.nsamples != 0 {
                    smpl_header(t, rt)?;
                    if t.is_str {
                        return Err(run_error("todo: Type=String for COUNT on FORMAT fields?"));
                    }
                    for i in 0..t.nsamples as usize {
                        if !t.usmpl(i) {
                            continue;
                        }
                        let p = (i * t.nval1.max(0) as usize) as isize;
                        let cnt = (0..t.nval1 as isize)
                            .filter(|&j| !is_mv(t.values.get(p + j)))
                            .count();
                        rt.values.set(i as isize, cnt as f64);
                    }
                    return Ok(1);
                }
                func_count(t, rt)?;
                Ok(1)
            }
            Func::Strlen => {
                let a = top(&self.stack)?;
                let (left, rt) = split(&mut self.filters, r);
                func_strlen(&left[a], rt);
                Ok(1)
            }
            Func::Binom => self.func_binom(r, h, rec, nstack),
            Func::Fisher => self.func_fisher(r, h, rec, nstack),
            Func::Perl => Err(run_error(
                "Perl filtering requires running `configure --enable-perl-filters` at compile time.",
            )),
        }
    }
}

/// The per-sample result shape of the `SMPL_*` functions.
fn smpl_header(t: &Token, rt: &mut Token) -> Result<()> {
    rt.nsamples = t.nsamples;
    rt.nvalues = t.nsamples;
    rt.nval1 = 1;
    rt.values.expand(rt.nvalues.max(0) as usize);
    let Some(u) = &t.usmpl else {
        return Err(unsupported("assertion failed in filter.c"));
    };
    let n = t.nsamples as usize;
    if u.len() < n {
        undefined("sample mask read past its end");
    }
    rt.usmpl = Some(
        u.iter()
            .copied()
            .chain(std::iter::repeat(0))
            .take(n)
            .collect(),
    );
    Ok(())
}

/// `func_max`, `func_min`, `func_avg`, `func_sum` over all used values.
fn reduce(t: &Token, rt: &mut Token, red: Reduce) {
    rt.nvalues = 0;
    if t.nvalues == 0 {
        return;
    }
    let mut vals = Vec::new();
    if t.nsamples != 0 {
        for i in 0..t.nsamples as usize {
            if !t.usmpl(i) {
                continue;
            }
            let p = (i * t.nval1.max(0) as usize) as isize;
            for j in 0..t.nval1 as isize {
                let v = t.values.get(p + j);
                if !is_mv(v) {
                    vals.push(v);
                }
            }
        }
    } else {
        for i in 0..t.nvalues as isize {
            let v = t.values.get(i);
            if !is_mv(v) {
                vals.push(v);
            }
        }
    }
    if vals.is_empty() {
        return;
    }
    let v = match red {
        Reduce::Max => vals
            .iter()
            .fold(f64::NEG_INFINITY, |m, &v| if m < v { v } else { m }),
        Reduce::Min => vals
            .iter()
            .fold(f64::INFINITY, |m, &v| if m > v { v } else { m }),
        Reduce::Avg => vals.iter().fold(0.0, |s, &v| s + v) / vals.len() as f64,
        Reduce::Sum => vals.iter().fold(0.0, |s, &v| s + v),
    };
    rt.values.set(0, v);
    rt.nvalues = 1;
}

fn stddev(vals: &[f64]) -> f64 {
    let n = vals.len() as f64;
    let mut avg = 0.0;
    for &v in vals {
        avg += v;
    }
    avg /= n;
    let mut sdev = 0.0;
    for &v in vals {
        sdev += (v - avg) * (v - avg);
    }
    (sdev / n).sqrt()
}

/// glibc's `qsort` with `compare_doubles`: a merge sort for arrays that
/// fit its temporary buffer, so the order of values that compare equal
/// (NaNs, signed zeros) is reproduced.
fn msort(v: &mut [f64]) {
    fn cmp(a: f64, b: f64) -> i32 {
        if a < b {
            -1
        } else if a > b {
            1
        } else {
            0
        }
    }
    fn rec(b: &mut [f64], tmp: &mut Vec<f64>) {
        let n = b.len();
        if n <= 1 {
            return;
        }
        let n1 = n / 2;
        let (b1, b2) = b.split_at_mut(n1);
        rec(b1, tmp);
        rec(b2, tmp);
        tmp.clear();
        let (mut i, mut j) = (0, 0);
        while i < b1.len() && j < b2.len() {
            if cmp(b1[i], b2[j]) <= 0 {
                tmp.push(b1[i]);
                i += 1;
            } else {
                tmp.push(b2[j]);
                j += 1;
            }
        }
        tmp.extend_from_slice(&b1[i..]);
        let k = tmp.len();
        b[..k].copy_from_slice(&tmp[..k]);
    }
    let mut tmp = Vec::with_capacity(v.len());
    rec(v, &mut tmp);
}

/// `func_median` / `func_stddev`: compacts the argument's used values in
/// place, as filter.c does.
fn median_or_stddev(t: &mut Token, rt: &mut Token, median: bool) {
    rt.nvalues = 0;
    if t.nvalues == 0 {
        return;
    }
    let mut n = 0isize;
    if t.nsamples != 0 {
        let mut k = 0isize;
        for i in 0..t.nsamples as usize {
            if !t.usmpl(i) {
                k += t.nval1 as isize;
                continue;
            }
            for _ in 0..t.nval1 {
                let v = t.values.get(k);
                if !is_mv(v) {
                    if n < k {
                        t.values.set(n, v);
                    }
                    n += 1;
                }
                k += 1;
            }
        }
    } else {
        for i in 0..t.nvalues as isize {
            let v = t.values.get(i);
            if is_mv(v) {
                continue;
            }
            if n < i {
                t.values.set(n, v);
            }
            n += 1;
        }
    }
    if n == 0 {
        return;
    }
    let mut vals: Vec<f64> = (0..n).map(|k| t.values.get(k)).collect();
    let v = if median {
        if n == 1 {
            vals[0]
        } else {
            msort(&mut vals);
            for (k, &v) in vals.iter().enumerate() {
                t.values.set(k as isize, v);
            }
            let n = n as usize;
            if n % 2 == 1 {
                vals[n / 2]
            } else {
                (vals[n / 2 - 1] + vals[n / 2]) * 0.5
            }
        }
    } else if n == 1 {
        0.0
    } else {
        stddev(&vals)
    };
    rt.values.set(0, v);
    rt.nvalues = 1;
}

/// `func_abs`.
fn func_abs(t: &Token, rt: &mut Token) -> Result<()> {
    if t.is_str {
        return Err(run_error("ABS() can be applied only on numeric values"));
    }
    rt.nsamples = t.nsamples;
    rt.nvalues = t.nvalues;
    rt.nval1 = t.nval1;
    rt.values.expand(rt.nvalues.max(0) as usize);
    if let Some(u) = &t.usmpl {
        let n = t.nsamples.max(0) as usize;
        rt.usmpl = Some(
            u.iter()
                .copied()
                .chain(std::iter::repeat(0))
                .take(n)
                .collect(),
        );
    }
    if t.nvalues == 0 {
        return Ok(());
    }
    if t.usmpl.is_some() {
        let mut k = 0isize;
        for i in 0..t.nsamples as usize {
            if !t.usmpl(i) {
                k += t.nval1 as isize;
                continue;
            }
            for _ in 0..t.nval1 {
                let v = t.values.get(k);
                rt.values.set(k, if is_mv(v) { missing() } else { v.abs() });
                k += 1;
            }
        }
    } else {
        for i in 0..t.nvalues as isize {
            let v = t.values.get(i);
            if is_missing(v) {
                rt.values.set(i, missing());
            } else if !is_vector_end(v) {
                rt.values.set(i, v.abs());
            }
        }
    }
    Ok(())
}

/// `func_phred`.
fn func_phred(t: &Token, rt: &mut Token) -> Result<()> {
    if t.is_str {
        return Err(run_error("PHRED() can be applied only on numeric values"));
    }
    rt.nsamples = t.nsamples;
    rt.nval1 = t.nval1;
    for i in 0..rt.nsamples.max(0) as usize {
        let v = t.pass_sample(i);
        rt.set_pass_sample(i, v);
    }
    if rt.usmpl.is_none() {
        let n = t.nsamples.max(0) as usize;
        rt.usmpl = Some(match &t.usmpl {
            Some(u) => u
                .iter()
                .copied()
                .chain(std::iter::repeat(0))
                .take(n)
                .collect(),
            None if n == 0 => Vec::new(),
            None => return Err(unsupported("copy from a NULL sample mask")),
        });
    }
    rt.nvalues = t.nvalues;
    if t.nvalues == 0 {
        return Ok(());
    }
    rt.values.expand(rt.nvalues.max(0) as usize);
    let phred = |v: f64| -4.34294481903 * v.ln();
    if t.usmpl.is_some() {
        let mut k = 0isize;
        for i in 0..t.nsamples as usize {
            if !t.usmpl(i) {
                k += t.nval1 as isize;
                continue;
            }
            for _ in 0..t.nval1 {
                let v = t.values.get(k);
                rt.values
                    .set(k, if is_mv(v) { missing() } else { phred(v) });
                k += 1;
            }
        }
    } else {
        for i in 0..t.nvalues as isize {
            let v = t.values.get(i);
            rt.values
                .set(i, if is_mv(v) { missing() } else { phred(v) });
        }
    }
    Ok(())
}

/// `func_count`.
fn func_count(t: &Token, rt: &mut Token) -> Result<()> {
    let mut cnt = 0i64;
    if t.tag.is_some() && t.nsamples != 0 {
        if t.is_str {
            return Err(run_error("todo: Type=String for COUNT on FORMAT fields?"));
        }
        for i in 0..t.nsamples as usize {
            if !t.usmpl(i) {
                continue;
            }
            let p = (i * t.nval1.max(0) as usize) as isize;
            for j in 0..t.nval1 as isize {
                if !is_mv(t.values.get(p + j)) {
                    cnt += 1;
                }
            }
        }
    } else if t.nsamples != 0 {
        for i in 0..t.nsamples as usize {
            if t.pass_sample(i) != 0 {
                cnt += 1;
            }
        }
    } else if t.is_str {
        if t.str_value.l > 0 {
            cnt = 1;
        }
        for k in 0..t.str_value.l {
            if t.str_value.get(k) == b',' {
                cnt += 1;
            }
        }
    } else {
        cnt = i64::from(t.nvalues);
    }
    rt.nvalues = 1;
    rt.values.set(0, cnt as f64);
    Ok(())
}

/// `func_strlen`.
fn func_strlen(t: &Token, rt: &mut Token) {
    rt.nvalues = 0;
    rt.str_value.l = 0;
    if t.str_value.l == 0 {
        return;
    }
    let s = &t.str_value;
    if t.idx == -2 {
        let mut i = 0;
        let mut ss = 0;
        while s.get(ss) != 0 {
            let mut se = ss;
            while s.get(se) != 0 && s.get(se) != b',' {
                se += 1;
            }
            rt.values.expand(i + 1);
            rt.values.set(i as isize, (se - ss) as f64);
            ss = if s.get(se) != 0 { se + 1 } else { se };
            i += 1;
        }
        rt.nvalues = i as i32;
    } else {
        let v = if s.get(1) == 0 && s.get(0) == b'.' {
            0.0
        } else {
            s.cstrlen(0) as f64
        };
        rt.values.set(0, v);
        rt.nvalues = 1;
    }
}

/// `vector_logic_or`.
fn logic_or(a: &Token, b: &Token, rt: &mut Token) {
    init_samples(a, b, rt);
    if a.pass_site == 0 && b.pass_site == 0 {
        return;
    }
    rt.pass_site = 1;
    if a.nsamples == 0 && b.nsamples == 0 {
        return;
    }
    if rt.tok_type == TOK_OR_VEC {
        if (a.nsamples == 0 && a.pass_site == 0) || (b.nsamples == 0 && b.pass_site == 0) {
            let t = if a.nsamples != 0 { a } else { b };
            for i in 0..rt.nsamples as usize {
                if rt.usmpl(i) {
                    rt.set_pass_sample(i, t.pass_sample(i));
                }
            }
        } else {
            for i in 0..rt.nsamples as usize {
                if rt.usmpl(i) {
                    rt.set_pass_sample(i, 1);
                }
            }
        }
        return;
    }
    if a.nsamples == 0 || b.nsamples == 0 {
        let t = if a.nsamples != 0 { a } else { b };
        for i in 0..rt.nsamples as usize {
            if rt.usmpl(i) {
                rt.set_pass_sample(i, t.pass_sample(i));
            }
        }
        return;
    }
    if a.nsamples != b.nsamples {
        undefined("assertion failed in filter.c");
        return;
    }
    for i in 0..rt.nsamples as usize {
        if rt.usmpl(i) {
            rt.set_pass_sample(i, a.pass_sample(i) | b.pass_sample(i));
        }
    }
}

/// `vector_logic_and`.
fn logic_and(a: &Token, b: &Token, rt: &mut Token) {
    init_samples(a, b, rt);
    if a.pass_site == 0 || b.pass_site == 0 {
        return;
    }
    if a.nsamples == 0 && b.nsamples == 0 {
        rt.pass_site = 1;
        return;
    }
    if a.nsamples == 0 || b.nsamples == 0 {
        let t = if a.nsamples != 0 { a } else { b };
        for i in 0..rt.nsamples as usize {
            if rt.usmpl(i) {
                rt.set_pass_sample(i, t.pass_sample(i));
            }
        }
        rt.pass_site = 1;
        return;
    }
    if a.nsamples != b.nsamples {
        undefined("assertion failed in filter.c");
        return;
    }
    if rt.tok_type == TOK_AND_VEC {
        for i in 0..rt.nsamples as usize {
            if rt.usmpl(i) {
                rt.set_pass_sample(i, a.pass_sample(i) | b.pass_sample(i));
            }
        }
        rt.pass_site = 1;
    } else {
        for i in 0..rt.nsamples as usize {
            if !rt.usmpl(i) {
                continue;
            }
            let v = a.pass_sample(i) & b.pass_sample(i);
            rt.set_pass_sample(i, v);
            if v != 0 {
                rt.pass_site = 1;
            }
        }
    }
}

impl Filter {
    /// The genotypes `binom()` and `fisher()` take allele indices from:
    /// (values, ploidy), or `None` without usable diploid genotypes.
    fn diploid_gts(&mut self, h: &Header, rec: &Record) -> Result<Option<usize>> {
        let mut gt = Vec::new();
        let ngt = match rec.genotype_values(h, &mut gt) {
            Ok(n) => {
                fill_ints(&mut self.tmpi, n, &gt);
                n as i32
            }
            Err(e) => e.code(),
        };
        if rec.n_sample() == 0 {
            return Err(unsupported("integer division by zero in filter.c"));
        }
        let max_ploidy = ngt / rec.n_sample() as i32;
        if ngt <= 0 || max_ploidy < 2 {
            return Ok(None);
        }
        Ok(Some(max_ploidy as usize))
    }

    /// The two allele indices of sample `i`, or `None` when either is
    /// missing.
    fn sample_alleles(
        &self,
        h: &Header,
        rec: &Record,
        i: usize,
        ploidy: usize,
    ) -> Result<Option<(i64, i64)>> {
        let p0 = i64::from(read_int(&self.tmpi, (i * ploidy) as isize));
        let p1 = i64::from(read_int(&self.tmpi, (i * ploidy + 1) as isize));
        if gt_missing(p0) || gt_missing(p1) || p1 == i64::from(INT32_VECTOR_END) {
            return Ok(None);
        }
        let (a1, a2) = (gt_allele(p0), gt_allele(p1));
        let n = rec.n_alleles() as i64;
        if a1 >= n || a2 >= n {
            let chrom = h
                .contigs
                .get(rec.rid as usize)
                .map_or("", |c| c.name.as_str());
            return Err(run_error(format!(
                "Incorrect allele index at {chrom}:{}, sample {}",
                one_based_pos(rec)?,
                h.samples.get(i).map_or("", String::as_str)
            )));
        }
        Ok(Some((a1, a2)))
    }

    /// `func_binom`.
    fn func_binom(&mut self, r: usize, h: &Header, rec: &Record, nstack: usize) -> Result<i32> {
        let nargs = self.filters[r].nargs;
        if nargs != 2 && nargs != 1 {
            return Err(run_error("Error: binom() takes one or two arguments"));
        }
        let istack = nstack as isize - nargs as isize;
        if istack < 0 {
            return Err(unsupported("assertion failed in filter.c"));
        }
        let istack = istack as usize;
        let ta = self.stack[istack];
        let tb = (istack + 1 < nstack).then(|| self.stack[istack + 1]);
        if self.filters[ta].nsamples != 0 {
            let ns = self.filters[ta].nsamples;
            {
                let (left, rt) = split(&mut self.filters, r);
                let t = &left[ta];
                rt.nval1 = 1;
                rt.nvalues = ns;
                rt.nsamples = ns;
                rt.values.expand(ns as usize);
                let Some(u) = &t.usmpl else {
                    return Err(unsupported("assertion failed in filter.c"));
                };
                rt.usmpl = Some(
                    u.iter()
                        .copied()
                        .chain(std::iter::repeat(0))
                        .take(ns as usize)
                        .collect(),
                );
            }
            if istack + 1 == nstack {
                let Some(ploidy) = self.diploid_gts(h, rec)? else {
                    let rt = &mut self.filters[r];
                    for i in 0..ns as usize {
                        if rt.usmpl(i) {
                            rt.values.set(i as isize, missing());
                        }
                    }
                    return Ok(nargs);
                };
                for i in 0..ns as usize {
                    if !self.filters[r].usmpl(i) {
                        continue;
                    }
                    let v = match self.sample_alleles(h, rec, i, ploidy)? {
                        None => missing(),
                        Some((a1, a2)) => {
                            let t = &self.filters[ta];
                            let base = (t.nval1 as i64 * i as i64) as isize;
                            let (x, y) = (
                                t.values.get(base + a1 as isize),
                                t.values.get(base + a2 as isize),
                            );
                            if is_mv(x) || is_mv(y) {
                                missing()
                            } else {
                                let p = calc_binom_two_sided(c_int(x), c_int(y), 0.5);
                                if p < 0.0 { missing() } else { p }
                            }
                        }
                    };
                    self.filters[r].values.set(i as isize, v);
                }
            } else {
                let tb = tb.unwrap();
                let (left, rt) = split(&mut self.filters, r);
                let (t, t2) = (&left[ta], &left[tb]);
                if t.nval1 != 1 || t2.nval1 != 1 {
                    let chrom = h
                        .contigs
                        .get(rec.rid as usize)
                        .map_or("", |c| c.name.as_str());
                    return Err(run_error(format!(
                        "Expected one value per binom() argument, found {} and {} at {chrom}:{}",
                        t.nval1,
                        t2.nval1,
                        one_based_pos(rec)?
                    )));
                }
                for i in 0..ns as usize {
                    if !rt.usmpl(i) {
                        continue;
                    }
                    let (x, y) = (
                        t.values.get(i as isize * t.nval1 as isize),
                        t2.values.get(i as isize * t2.nval1 as isize),
                    );
                    let v = if is_mv(x) || is_mv(y) {
                        missing()
                    } else {
                        let p = calc_binom_two_sided(c_int(x), c_int(y), 0.5);
                        if p < 0.0 { missing() } else { p }
                    };
                    rt.values.set(i as isize, v);
                }
            }
        } else {
            let (left, rt) = split(&mut self.filters, r);
            rt.nvalues = 1;
            rt.values.expand(1);
            let t = &left[ta];
            let pair = if istack + 1 == nstack {
                (t.nvalues == 2).then(|| (t.values.get(0), t.values.get(1)))
            } else {
                let t2 = &left[tb.unwrap()];
                (t.nvalues == 1 && t2.nvalues == 1).then(|| (t.values.get(0), t2.values.get(0)))
            };
            let v = match pair {
                Some((x, y)) if !is_mv(x) && !is_mv(y) => {
                    let p = calc_binom_two_sided(c_int(x), c_int(y), 0.5);
                    if p < 0.0 { missing() } else { p }
                }
                _ => missing(),
            };
            rt.values.set(0, v);
        }
        Ok(nargs)
    }

    /// `func_fisher`.
    fn func_fisher(&mut self, r: usize, h: &Header, rec: &Record, nstack: usize) -> Result<i32> {
        let nargs = self.filters[r].nargs;
        let istack = nstack as isize - nargs as isize;
        if istack < 0 {
            return Err(unsupported(
                "filter function arguments missing from the stack",
            ));
        }
        let istack = istack as usize;
        let ta = self.stack[istack];
        let tb = (istack + 2 == nstack).then(|| self.stack[istack + 1]);
        let fisher = |n: [f64; 4]| -> f64 {
            if n.iter().any(|&x| is_mv(x)) {
                missing()
            } else {
                kt_fisher_two(c_int(n[0]), c_int(n[1]), c_int(n[2]), c_int(n[3])).unwrap_or_else(
                    || {
                        undefined("fisher() on counts too large to sum");
                        missing()
                    },
                )
            }
        };
        if self.filters[ta].nsamples == 0 {
            let (left, rt) = split(&mut self.filters, r);
            rt.nvalues = 1;
            rt.values.expand(1);
            let t = &left[ta];
            let n = if istack + 1 == nstack {
                (t.nvalues == 4).then(|| {
                    [
                        t.values.get(0),
                        t.values.get(1),
                        t.values.get(2),
                        t.values.get(3),
                    ]
                })
            } else if let Some(tb) = tb {
                let t2 = &left[tb];
                (t.nvalues >= 2 && t2.nvalues >= 2).then(|| {
                    [
                        t.values.get(0),
                        t2.values.get(0),
                        t.values.get(1),
                        t2.values.get(1),
                    ]
                })
            } else {
                None
            };
            rt.values.set(0, n.map_or(missing(), fisher));
            return Ok(nargs);
        }
        let ns = self.filters[ta].nsamples;
        {
            let (left, rt) = split(&mut self.filters, r);
            let t = &left[ta];
            rt.nval1 = 1;
            rt.nvalues = ns;
            rt.nsamples = ns;
            rt.values.expand(ns as usize);
            let Some(u) = &t.usmpl else {
                return Err(unsupported("assertion failed in filter.c"));
            };
            rt.usmpl = Some(
                u.iter()
                    .copied()
                    .chain(std::iter::repeat(0))
                    .take(ns as usize)
                    .collect(),
            );
        }
        let (t_nval1, t_vl) = (self.filters[ta].nval1, self.filters[ta].vl_len);
        if istack + 1 == nstack && t_nval1 == 4 {
            let (left, rt) = split(&mut self.filters, r);
            let t = &left[ta];
            for i in 0..ns as usize {
                if !rt.usmpl(i) {
                    continue;
                }
                let b = 4 * i as isize;
                let n = [
                    t.values.get(b),
                    t.values.get(b + 1),
                    t.values.get(b + 2),
                    t.values.get(b + 3),
                ];
                rt.values.set(i as isize, fisher(n));
            }
        } else if let Some(tb) = tb {
            let (t2_vl, nu) = (
                self.filters[tb].vl_len,
                self.filters[ta].nuidxs + self.filters[tb].nuidxs,
            );
            if t_vl == VL_R && t2_vl == VL_R && nu == 4 {
                let (left, rt) = split(&mut self.filters, r);
                let (t, t2) = (&left[ta], &left[tb]);
                for i in 0..ns as usize {
                    if !rt.usmpl(i) {
                        continue;
                    }
                    let (b, b2) = (
                        t.nval1 as isize * i as isize,
                        t2.nval1 as isize * i as isize,
                    );
                    let n = [
                        t.values.get(b),
                        t.values.get(b + 1),
                        t2.values.get(b2),
                        t2.values.get(b2 + 1),
                    ];
                    rt.values.set(i as isize, fisher(n));
                }
            } else if t_vl == VL_R && t2_vl == VL_R {
                let Some(ploidy) = self.diploid_gts(h, rec)? else {
                    let rt = &mut self.filters[r];
                    for i in 0..ns as usize {
                        if rt.usmpl(i) {
                            rt.values.set(i as isize, missing());
                        }
                    }
                    return Ok(nargs);
                };
                for i in 0..ns as usize {
                    if !self.filters[r].usmpl(i) {
                        continue;
                    }
                    let v = match self.sample_alleles(h, rec, i, ploidy)? {
                        None => missing(),
                        Some((a1, a2)) => {
                            let (t, t2) = (&self.filters[ta], &self.filters[tb]);
                            let (b, b2) = (
                                t.nval1 as isize * i as isize,
                                t2.nval1 as isize * i as isize,
                            );
                            fisher([
                                t.values.get(b + a1 as isize),
                                t.values.get(b + a2 as isize),
                                t2.values.get(b2 + a1 as isize),
                                t2.values.get(b2 + a2 as isize),
                            ])
                        }
                    };
                    self.filters[r].values.set(i as isize, v);
                }
            }
        }
        Ok(nargs)
    }
}

#[cfg(test)]
mod ordered_vector_tests {
    use super::*;

    #[test]
    fn extrema_match_all_pairs_and_leave_unwritten_tails_to_short_circuit() {
        let pool = [
            f64::NEG_INFINITY,
            -16777217.0,
            -16777216.0,
            -1.0,
            -0.0,
            0.0,
            1.0,
            16777216.0,
            16777217.0,
            f64::INFINITY,
            f64::NAN,
            missing(),
            vector_end(),
        ];
        for &x in &pool {
            for &y in &pool {
                for mixed in [false, true] {
                    let values = |v| {
                        let mut data = vec![v; 16];
                        if mixed {
                            data.extend(pool);
                        }
                        Token {
                            nvalues: data.len() as i32,
                            values: Vals(data),
                            ..Token::default()
                        }
                    };
                    let (a, b) = (values(x), values(y));
                    for op in [TOK_LT, TOK_LE, TOK_BT, TOK_BE] {
                        for flags in 0..8 {
                            let ml = [flags & 1, (flags >> 1) & 1, (flags >> 2) & 1];
                            let expected = a.values.0.iter().any(|&x| {
                                b.values.0.iter().any(|&y| {
                                    let nmiss = usize::from(is_mv(x)) + usize::from(is_mv(y));
                                    if nmiss != 0 {
                                        ml[nmiss] != 0
                                    } else {
                                        compare(op, x, y)
                                    }
                                })
                            });
                            assert_eq!(ordered_vectors(&a, &b, op, &ml), Some(expected));
                        }
                    }
                }
            }
        }
        let mut a = Token {
            nvalues: 16,
            values: Vals(vec![0.0; 16]),
            ..Token::default()
        };
        let mut b = Token {
            nvalues: 16,
            values: Vals(vec![1.0; 16]),
            ..Token::default()
        };
        b.values.0[0] = -1.0;
        a.values.0[15] = f64::from_bits(UNINIT);
        assert_eq!(ordered_vectors(&a, &b, TOK_LT, &[0; 3]), None);
        take_undefined();
        let mut result = Token::default();
        cmp_vectors(&a, &b, &mut result, TOK_LT, &[0; 3]).unwrap();
        assert_eq!(result.pass_site, 1);
        assert!(take_undefined().is_none());
    }
}

#[cfg(test)]
mod equality_vector_tests {
    use super::*;

    fn token(values: Vec<f64>) -> Token {
        Token {
            nvalues: values.len() as i32,
            values: Vals(values),
            ..Token::default()
        }
    }

    // The site-level, nonempty CMP_VECTORS loop, including Vals::get's
    // initializedness tracking and the original missing-flag lookup.
    fn cartesian(a: &Token, b: &Token, op: i32, ml: &[i32; 3]) -> i32 {
        for i in 0..a.nvalues as isize {
            let amiss = usize::from(is_mv(a.values.get(i)));
            for j in 0..b.nvalues as isize {
                let nmiss = amiss + usize::from(is_mv(b.values.get(j)));
                if if nmiss != 0 {
                    ml[nmiss] != 0
                } else {
                    compare(op, a.values.get(i), b.values.get(j))
                } {
                    return 1;
                }
            }
        }
        0
    }

    #[test]
    fn equality_keys_match_pair_comparison_at_rounding_boundaries() {
        let pool = [
            f64::NEG_INFINITY,
            -f64::MAX,
            -16777217.0,
            -16777216.0,
            -1.00000001,
            -1.0,
            -f64::MIN_POSITIVE,
            -0.0,
            0.0,
            f64::MIN_POSITIVE,
            1.0,
            1.00000001,
            16777215.5,
            16777216.0,
            16777216.25,
            16777217.0,
            f64::INFINITY,
            f64::NAN,
        ];
        for &x in &pool {
            for &y in &pool {
                let same = equality_key(x)
                    .zip(equality_key(y))
                    .is_some_and(|(a, b)| a == b);
                assert_eq!(same, compare(TOK_EQ, x, y), "{x:?}, {y:?}");
            }
        }
    }

    #[test]
    fn equality_reduction_matches_cartesian_oracle_for_all_missing_flags() {
        let pool = [
            f64::NEG_INFINITY,
            -f64::MAX,
            -16777217.0,
            -16777216.0,
            -0.0,
            0.0,
            16777216.0,
            16777216.25,
            f64::INFINITY,
            f64::NAN,
            missing(),
            vector_end(),
        ];
        let mut inputs: Vec<Vec<f64>> = pool.iter().map(|&x| vec![x; 64]).collect();
        inputs.push((0..64).map(|i| f64::from(i * 2)).collect());
        inputs.push((0..64).map(|i| f64::from(i * 2 + 1)).collect());
        inputs.push(pool.repeat(6));
        let mut late = vec![1000.0; 64];
        late[63] = 126.0;
        inputs.push(late);
        for av in &inputs {
            for bv in &inputs {
                let (a, b) = (token(av.clone()), token(bv.clone()));
                for op in [TOK_EQ, TOK_NE] {
                    for flags in 0..8 {
                        let ml = [flags & 1, (flags >> 1) & 1, (flags >> 2) & 1];
                        take_undefined();
                        let expected = cartesian(&a, &b, op, &ml);
                        assert!(take_undefined().is_none());
                        assert_eq!(equality_vectors(&a, &b, op, &ml), Some(expected != 0));
                        let mut result = Token::default();
                        cmp_vectors(&a, &b, &mut result, op, &ml).unwrap();
                        assert_eq!(result.pass_site, expected);
                        assert!(take_undefined().is_none());
                    }
                }
            }
        }
    }

    #[test]
    fn partial_shadows_keep_cartesian_short_circuit_order() {
        for op in [TOK_EQ, TOK_NE] {
            let mut a = token(vec![0.0; 64]);
            let mut b = token(vec![if op == TOK_EQ { 1.0 } else { 0.0 }; 64]);
            a.values.0[1] = 1.0;
            a.values.0[63] = f64::from_bits(UNINIT);
            assert_eq!(equality_vectors(&a, &b, op, &[0; 3]), None);
            take_undefined();
            let mut result = Token::default();
            cmp_vectors(&a, &b, &mut result, op, &[0; 3]).unwrap();
            assert_eq!(result.pass_site, 1);
            assert!(take_undefined().is_none());

            // A hit in the first row must not inspect its unwritten tail.
            b.values.0[0] = if op == TOK_EQ { 0.0 } else { 1.0 };
            b.values.0[63] = f64::from_bits(UNINIT);
            assert_eq!(equality_vectors(&a, &b, op, &[0; 3]), Some(true));

            // An unwritten value before the only match is an error.
            a.values.0[1] = f64::from_bits(UNINIT);
            b.values.0[0] = if op == TOK_EQ { 1.0 } else { 0.0 };
            b.values.0[63] = if op == TOK_EQ { 0.0 } else { 1.0 };
            assert_eq!(equality_vectors(&a, &b, op, &[0; 3]), Some(true));
            b.values.0[1] = f64::from_bits(UNINIT);
            assert_eq!(equality_vectors(&a, &b, op, &[0; 3]), None);
            take_undefined();
            let mut result = Token::default();
            cmp_vectors(&a, &b, &mut result, op, &[0; 3]).unwrap();
            assert!(take_undefined().is_some());
        }
    }

    #[test]
    fn partial_shadows_match_original_accesses_and_all_missing_flags() {
        for op in [TOK_EQ, TOK_NE] {
            for unwritten in [0, 1, 32, 63] {
                for hit in [0, 1, 63] {
                    for flags in 0..8 {
                        let ml = [flags & 1, (flags >> 1) & 1, (flags >> 2) & 1];
                        for mv in [missing(), vector_end(), f64::NAN] {
                            let mut a = token(vec![0.0; 64]);
                            let mut b = token(vec![if op == TOK_EQ { 1.0 } else { 0.0 }; 64]);
                            a.values.0[hit] = 1.0;
                            b.values.0[hit] = mv;
                            for side in [false, true] {
                                let (mut a, mut b) = (a.clone(), b.clone());
                                if side {
                                    a.values.0[unwritten] = f64::from_bits(UNINIT);
                                } else {
                                    b.values.0[unwritten] = f64::from_bits(UNINIT);
                                }
                                take_undefined();
                                let expected = cartesian(&a, &b, op, &ml);
                                let undefined = take_undefined();
                                let mut result = Token::default();
                                cmp_vectors(&a, &b, &mut result, op, &ml).unwrap();
                                assert_eq!(result.pass_site, expected);
                                assert_eq!(take_undefined(), undefined);
                            }
                        }
                    }
                }
            }
        }
    }
}
