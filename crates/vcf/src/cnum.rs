//! Number parsing as the C code this crate matches does it: htslib's
//! `hts_str2int`, `hts_str2uint` and `hts_str2dbl`, and glibc's `strtod`
//! (C locale). Each returns the value and the number of bytes
//! consumed; inputs are read as C strings (the end of the slice and a NUL
//! both end them).

fn at(s: &[u8], i: usize) -> u8 {
    match s.get(i) {
        Some(&c) => c,
        None => 0,
    }
}

/// C locale `isspace`.
pub(crate) fn is_space(c: u8) -> bool {
    matches!(c, b' ' | b'\t' | b'\n' | 0x0b | 0x0c | b'\r')
}

/// `hts_str2int`: (value, consumed, overflowed). On overflow the value is
/// the limit of `bits` and all following digits are consumed.
pub(crate) fn hts_str2int(s: &[u8], bits: u32) -> (i64, usize, bool) {
    let mut limit: u64 = (1u64 << (bits - 1)) - 1;
    let mut fast: u32 = (bits - 1) * 1000 / 3322 + 1;
    let mut i = 0;
    let neg = match at(s, 0) {
        b'-' => {
            limit += 1;
            i = 1;
            true
        }
        b'+' => {
            i = 1;
            false
        }
        _ => false,
    };
    let mut n: u64 = 0;
    loop {
        fast -= 1;
        if fast == 0 || !at(s, i).is_ascii_digit() {
            break;
        }
        n = n.wrapping_mul(10).wrapping_add(u64::from(at(s, i) - b'0'));
        i += 1;
    }
    let mut failed = false;
    if at(s, i) >= b'0' && fast == 0 {
        let (d10, m10) = (limit / 10, limit % 10);
        while at(s, i).is_ascii_digit() {
            let d = u64::from(at(s, i) - b'0');
            if n < d10 || (n == d10 && d <= m10) {
                n = n * 10 + d;
                i += 1;
            } else {
                while at(s, i).is_ascii_digit() {
                    i += 1;
                }
                n = limit;
                failed = true;
                break;
            }
        }
    }
    let v = if neg {
        n.wrapping_neg() as i64
    } else {
        n as i64
    };
    (v, i, failed)
}

/// `hts_str2uint`: (value, consumed, overflowed). Digits read before the
/// overflow check starts wrap, as the C code's do.
pub(crate) fn hts_str2uint(s: &[u8], bits: u32) -> (u64, usize, bool) {
    let limit: u64 = if bits < 64 {
        (1u64 << bits) - 1
    } else {
        u64::MAX
    };
    let mut fast: u32 = bits * 1000 / 3322 + 1;
    let mut i = usize::from(at(s, 0) == b'+');
    let mut n: u64 = 0;
    loop {
        fast -= 1;
        if fast == 0 || !at(s, i).is_ascii_digit() {
            break;
        }
        n = n.wrapping_mul(10).wrapping_add(u64::from(at(s, i) - b'0'));
        i += 1;
    }
    let mut failed = false;
    if at(s, i).is_ascii_digit() && fast == 0 {
        let (d10, m10) = (limit / 10, limit % 10);
        while at(s, i).is_ascii_digit() {
            let d = u64::from(at(s, i) - b'0');
            if n < d10 || (n == d10 && d <= m10) {
                n = n * 10 + d;
                i += 1;
            } else {
                while at(s, i).is_ascii_digit() {
                    i += 1;
                }
                n = limit;
                failed = true;
                break;
            }
        }
    }
    (n, i, failed)
}

/// `hts_str2dbl`: (value, consumed, failed). Short plain decimals are
/// divided out exactly; everything else goes to `strtod`.
pub(crate) fn hts_str2dbl(s: &[u8]) -> (f64, usize, bool) {
    const D: [f64; 22] = [
        1.0, 1.0, 1e1, 1e2, 1e3, 1e4, 1e5, 1e6, 1e7, 1e8, 1e9, 1e10, 1e11, 1e12, 1e13, 1e14, 1e15,
        1e16, 1e17, 1e18, 1e19, 1e20,
    ];
    let fallback = || {
        let (d, n, _) = strtod(s);
        (d, n, n == 0)
    };
    let mut v = 0;
    while is_space(at(s, v)) {
        v += 1;
    }
    let mut neg = false;
    match at(s, v) {
        b'-' => {
            neg = true;
            v += 1;
        }
        b'+' => v += 1,
        _ => {}
    }
    match at(s, v) {
        b'1'..=b'9' => {}
        b'0' if !matches!(at(s, v + 1), b'x' | b'X') => {}
        _ => return fallback(),
    }
    while at(s, v) == b'0' {
        v += 1;
    }
    let start = v;
    let mut max_len = 15;
    let mut n: u64 = 0;
    let mut point: Option<usize> = None;
    loop {
        max_len -= 1;
        if max_len == 0 || !at(s, v).is_ascii_digit() {
            break;
        }
        n = n * 10 + u64::from(at(s, v) - b'0');
        v += 1;
    }
    if max_len != 0 && at(s, v) == b'.' {
        point = Some(v - start);
        v += 1;
        loop {
            max_len -= 1;
            if max_len == 0 || !at(s, v).is_ascii_digit() {
                break;
            }
            n = n * 10 + u64::from(at(s, v) - b'0');
            v += 1;
        }
    }
    let point = point.unwrap_or(v - start);
    if max_len == 0 || matches!(at(s, v), b'e' | b'E') {
        return fallback();
    }
    let d = n as f64 / D[v - start - point];
    (if neg { -d } else { d }, v, false)
}

/// glibc `strtod` in the C locale: (value, consumed, ERANGE). Decimal,
/// hexadecimal, `inf`/`infinity` and `nan`/`nan(chars)` forms; correctly
/// rounded. ERANGE is reported on overflow and when a nonzero input
/// rounds to a subnormal or zero.
pub(crate) fn strtod(s: &[u8]) -> (f64, usize, bool) {
    let mut i = 0;
    while is_space(at(s, i)) {
        i += 1;
    }
    let neg = match at(s, i) {
        b'-' => {
            i += 1;
            true
        }
        b'+' => {
            i += 1;
            false
        }
        _ => false,
    };
    let sign = |x: f64| if neg { -x } else { x };
    let lower = |j: usize| at(s, j).to_ascii_lowercase();
    let word = |j: usize, w: &[u8]| w.iter().enumerate().all(|(k, &c)| lower(j + k) == c);
    if word(i, b"inf") {
        let n = if word(i + 3, b"inity") { i + 8 } else { i + 3 };
        return (sign(f64::INFINITY), n, false);
    }
    if word(i, b"nan") {
        let mut n = i + 3;
        let mut payload = 0;
        if at(s, n) == b'(' {
            let mut k = n + 1;
            while at(s, k).is_ascii_alphanumeric() || at(s, k) == b'_' {
                k += 1;
            }
            if at(s, k) == b')' {
                payload = nan_payload(&s[n + 1..k]);
                n = k + 1;
            }
        }
        return (
            sign(f64::from_bits(0x7ff8_0000_0000_0000 | payload)),
            n,
            false,
        );
    }
    if at(s, i) == b'0' && lower(i + 1) == b'x' {
        let h = i + 2;
        let hex = |j: usize| at(s, j).is_ascii_hexdigit();
        if hex(h) || (at(s, h) == b'.' && hex(h + 1)) {
            let (v, n, erange) = hex_float(s, h);
            return (sign(v), n, erange);
        }
        // Only the `0` is a number.
        return (sign(0.0), i + 1, false);
    }
    // Decimal.
    let mut j = i;
    let mut any = false;
    let mut nonzero = false;
    while at(s, j).is_ascii_digit() {
        nonzero |= at(s, j) != b'0';
        j += 1;
        any = true;
    }
    if at(s, j) == b'.' {
        let mut k = j + 1;
        let mut frac = false;
        while at(s, k).is_ascii_digit() {
            nonzero |= at(s, k) != b'0';
            k += 1;
            frac = true;
        }
        if any || frac {
            j = k;
            any = true;
        }
    }
    if !any {
        return (0.0, 0, false);
    }
    if lower(j) == b'e' {
        let mut k = j + 1;
        if matches!(at(s, k), b'+' | b'-') {
            k += 1;
        }
        if at(s, k).is_ascii_digit() {
            while at(s, k).is_ascii_digit() {
                k += 1;
            }
            j = k;
        }
    }
    // Every byte in this span was checked as decimal syntax above. Borrow
    // the original digits so exponent magnitude and leading zeros are kept.
    let text = std::str::from_utf8(&s[i..j]).expect("decimal syntax is ASCII");
    let v: f64 = text.parse().unwrap_or(0.0);
    let erange = v.is_infinite() || (nonzero && v.abs() < f64::MIN_POSITIVE);
    (sign(v), j, erange)
}

/// glibc interprets an entirely numeric NaN payload as a base-0 unsigned
/// integer. Unrecognised text keeps the default quiet NaN mantissa.
fn nan_payload(s: &[u8]) -> u64 {
    let (base, digits) = if s.starts_with(b"0x") || s.starts_with(b"0X") {
        (16, &s[2..])
    } else if s.starts_with(b"0") {
        (8, s)
    } else {
        (10, s)
    };
    if digits.is_empty() {
        return 0;
    }
    let mut n = 0u64;
    for &c in digits {
        let Some(digit) = (c as char).to_digit(base) else {
            return 0;
        };
        n = n
            .saturating_mul(u64::from(base))
            .saturating_add(u64::from(digit));
    }
    n & 0x000f_ffff_ffff_ffff
}

/// A hexadecimal float from its first digit (after `0x`), correctly
/// rounded to nearest-even.
fn hex_float(s: &[u8], start: usize) -> (f64, usize, bool) {
    let digit = |c: u8| (c as char).to_digit(16);
    let mut j = start;
    // Mantissa: up to 60 significant bits plus a sticky bit.
    let mut m: u64 = 0;
    let mut exp: i64 = 0;
    let mut sticky = false;
    let mut seen_point = false;
    loop {
        let c = at(s, j);
        if c == b'.' && !seen_point {
            seen_point = true;
            j += 1;
            continue;
        }
        let Some(d) = digit(c) else { break };
        if m >> 56 == 0 {
            m = m << 4 | u64::from(d);
            if seen_point {
                exp -= 4;
            }
        } else {
            sticky |= d != 0;
            if !seen_point {
                exp += 4;
            }
        }
        j += 1;
    }
    if at(s, j) | 0x20 == b'p' {
        let mut k = j + 1;
        let neg = match at(s, k) {
            b'-' => {
                k += 1;
                true
            }
            b'+' => {
                k += 1;
                false
            }
            _ => false,
        };
        if at(s, k).is_ascii_digit() {
            let mut e: i64 = 0;
            while at(s, k).is_ascii_digit() {
                e = (e * 10 + i64::from(at(s, k) - b'0')).min(1 << 40);
                k += 1;
            }
            exp += if neg { -e } else { e };
            j = k;
        }
    }
    if m == 0 {
        return (0.0, j, false);
    }
    // Normalise so that bit 63 is the leading one: value = m * 2^exp.
    let lz = m.leading_zeros();
    m <<= lz;
    exp -= i64::from(lz);
    // Unbiased exponent of the leading bit.
    let top = exp + 63;
    if top > 1023 {
        return (f64::INFINITY, j, true);
    }
    // Bits to keep: 53 for normal numbers, fewer for subnormal ones.
    let keep: i64 = if top >= -1022 { 53 } else { 53 - (-1022 - top) };
    if keep < 0 {
        return (0.0, j, true);
    }
    let drop = 64 - keep as u32;
    let (mut q, rest) = if drop >= 64 {
        (0u64, m)
    } else {
        (m >> drop, m & ((1u64 << drop) - 1))
    };
    let half = if drop >= 64 {
        1u64 << 63
    } else {
        1u64 << (drop - 1)
    };
    let inexact = rest != 0 || sticky;
    let round_up = rest > half || (rest == half && (sticky || q & 1 == 1));
    if round_up {
        q += 1;
    }
    // q * 2^(exp + drop)
    let v = q as f64 * 2f64.powi(0) * pow2(exp + i64::from(drop));
    let erange = v.is_infinite() || (inexact && v < f64::MIN_POSITIVE);
    (v, j, erange)
}

/// 2^e as a double, exact over the range a hex float can produce.
fn pow2(e: i64) -> f64 {
    // Split to stay within the normal exponent range at each step.
    let mut v = 1.0f64;
    let mut e = e;
    while e > 1000 {
        v *= 2f64.powi(1000);
        e -= 1000;
    }
    while e < -1000 {
        v *= 2f64.powi(-1000);
        e += 1000;
    }
    v * 2f64.powi(e as i32)
}
