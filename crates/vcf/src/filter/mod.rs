//! bcftools' `-i`/`-e` filter expressions: a rewrite of `filter.c`
//! (bcftools 1.23.1, built without Perl filters) with the same language,
//! the same errors and the same results. [`Filter::new`] compiles an
//! expression as `filter_init` does; [`Filter::test`] evaluates it on a
//! record as `filter_test(filter, rec, NULL)` does, the site-level test
//! `bcftools consensus` uses.
//!
//! The evaluation keeps filter.c's state between records where filter.c
//! keeps it (token buffers, sample masks, string buffers), so quirks that
//! depend on earlier records are reproduced. filter.c caches genotype and
//! allele-count summaries per `bcf1_t` address; [`Filter::test`] behaves
//! as if every call passed a different record object.
//!
//! Where filter.c reports an error (its `error()`), the result is
//! [`Error::Invalid`]; where its behaviour is undefined or it ends the
//! process another way (reads of memory it never wrote, failed
//! assertions, `exit`, division by zero), the result is
//! [`Error::Unsupported`]. Regex compilation and evaluation also return
//! `Unsupported` on the documented recursion or work limits, never a
//! successful false result for an unfinished match.

mod eval;
mod regex;
mod stats;

use std::cell::Cell;
use std::collections::{HashMap, HashSet};
use std::io::Read;

use omicsweave_core::{Error, Result};

use crate::cnum::{is_space, strtod};
use crate::header::{Header, LineKind};

use regex::Regex;

const TOK_VAL: i32 = 0;
const TOK_LFT: i32 = 1;
const TOK_RGT: i32 = 2;
const TOK_LE: i32 = 3;
const TOK_LT: i32 = 4;
const TOK_EQ: i32 = 5;
const TOK_BT: i32 = 6;
const TOK_BE: i32 = 7;
const TOK_NE: i32 = 8;
const TOK_OR: i32 = 9;
const TOK_AND: i32 = 10;
const TOK_ADD: i32 = 11;
const TOK_SUB: i32 = 12;
const TOK_MULT: i32 = 13;
const TOK_DIV: i32 = 14;
const TOK_MAX: i32 = 15;
const TOK_MIN: i32 = 16;
const TOK_AVG: i32 = 17;
const TOK_AND_VEC: i32 = 18;
const TOK_OR_VEC: i32 = 19;
const TOK_LIKE: i32 = 20;
const TOK_NLIKE: i32 = 21;
const TOK_SUM: i32 = 22;
const TOK_ABS: i32 = 23;
const TOK_LEN: i32 = 24;
const TOK_FUNC: i32 = 25;
const TOK_CNT: i32 = 26;
const TOK_PERLSUB: i32 = 27;
const TOK_BINOM: i32 = 28;
const TOK_PHRED: i32 = 29;
const TOK_MEDIAN: i32 = 30;
const TOK_STDEV: i32 = 31;
const TOK_SMAX: i32 = 32;
const TOK_SMIN: i32 = 33;
const TOK_SAVG: i32 = 34;
const TOK_SMEDIAN: i32 = 35;
const TOK_SSTDEV: i32 = 36;
const TOK_SSUM: i32 = 37;
const TOK_IN: i32 = 38;
const TOK_NOT_IN: i32 = 39;
const TOK_MODULO: i32 = 40;
const TOK_EXT: i32 = 41;
const TOK_FISHER: i32 = 42;
const TOK_SCOUNT: i32 = 43;

const OP_PREC: [i32; 44] = [
    0, 1, 1, 5, 5, 5, 5, 5, 5, 2, 3, 6, 6, 7, 7, 8, 8, 8, 3, 2, 5, 5, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8,
    8, 8, 8, 8, 8, 8, 8, 8, 7, 0, 8, 8,
];

const HL_FLT: i32 = 0;
const HL_INFO: i32 = 1;
const HL_FMT: i32 = 2;
const HT_FLAG: i32 = 0;
const HT_INT: i32 = 1;
const HT_REAL: i32 = 2;
const HT_STR: i32 = 3;
const VL_FIXED: i32 = 0;
const VL_R: i32 = 4;
const UN_STR: i32 = 1;
const UN_FLT: i32 = 2;
const UN_INFO: i32 = 4;
const UN_FMT: i32 = 8;

const DBL_MISSING: u64 = 0x7ff0_0000_0000_0001;
const DBL_VECTOR_END: u64 = 0x7ff0_0000_0000_0002;
/// Marks filter buffer slots filter.c has never written.
const UNINIT: u64 = 0x7ff4_0bad_f00d_0001;

fn missing() -> f64 {
    f64::from_bits(DBL_MISSING)
}
fn vector_end() -> f64 {
    f64::from_bits(DBL_VECTOR_END)
}
fn is_missing(x: f64) -> bool {
    x.to_bits() == DBL_MISSING
}
fn is_vector_end(x: f64) -> bool {
    x.to_bits() == DBL_VECTOR_END
}
fn is_mv(x: f64) -> bool {
    is_missing(x) || is_vector_end(x)
}

thread_local! {
    /// Set when the evaluation does what filter.c leaves undefined.
    static UNDEFINED: Cell<Option<&'static str>> = const { Cell::new(None) };
}

fn undefined(what: &'static str) {
    UNDEFINED.with(|u| {
        if u.get().is_none() {
            u.set(Some(what));
        }
    });
}

fn take_undefined() -> Option<&'static str> {
    UNDEFINED.with(|u| u.take())
}

fn unsupported(what: &'static str) -> Error {
    Error::Unsupported(what)
}

fn init_error(detail: impl Into<String>) -> Error {
    Error::invalid("filter expression", detail)
}

/// A C `double` array with its allocated length; slots never written read
/// as undefined.
#[derive(Debug, Clone, Default)]
struct Vals(Vec<f64>);

impl Vals {
    /// `hts_expand`: grows to at least `n` slots, new ones unwritten.
    fn expand(&mut self, n: usize) {
        if n > self.0.len() {
            self.0.resize(n, f64::from_bits(UNINIT));
        }
    }

    /// `hts_expand0`: grows to at least `n` slots, new ones zero.
    fn expand0(&mut self, n: usize) {
        if n > self.0.len() {
            self.0.resize(n, 0.0);
        }
    }

    fn get(&self, i: isize) -> f64 {
        match usize::try_from(i).ok().and_then(|i| self.0.get(i)) {
            Some(v) if v.to_bits() != UNINIT => *v,
            _ => {
                undefined("read of filter values never written");
                0.0
            }
        }
    }

    fn set(&mut self, i: isize, v: f64) {
        match usize::try_from(i).ok().and_then(|i| self.0.get_mut(i)) {
            Some(x) => *x = v,
            None => undefined("write past filter values"),
        }
    }
}

/// A kstring's nominal C capacity and exactly the bytes C has written.
/// The Vec is a contiguous initialized prefix; jump writes stay sparse.
#[derive(Debug, Clone, Default)]
struct KStr {
    buf: Vec<u8>,
    l: usize,
    m: usize,
    // Keep ordinary tokens compact; sparse storage is only allocated on jump writes.
    #[allow(clippy::box_collection)]
    sparse: Option<Box<HashMap<usize, u8>>>,
}

impl KStr {
    fn is_null(&self) -> bool {
        self.m == 0
    }

    /// A real C realloc preserves only offsets below its new size, even
    /// when that size grows relative to a narrowed nominal capacity.
    fn reallocate(&mut self, n: usize) {
        self.buf.truncate(n);
        if let Some(sparse) = &mut self.sparse {
            sparse.retain(|&i, _| i < n);
        }
        self.m = n;
    }

    /// ks_resize changes nominal capacity without materializing its holes.
    fn resize(&mut self, size: usize) {
        if self.m < size {
            let grown = if size > usize::MAX >> 2 {
                Some(size)
            } else {
                size.checked_add(size >> 1)
            };
            let Some(grown) = grown else {
                undefined("filter string capacity overflow");
                return;
            };
            self.reallocate(grown);
        }
    }

    fn grow_exact(&mut self, n: usize) {
        if self.m < n {
            self.reallocate(n);
        }
    }

    fn get(&self, i: usize) -> u8 {
        if let Some(&c) = self.buf.get(i) {
            return c;
        }
        if let Some(&c) = self.sparse.as_ref().and_then(|s| s.get(&i)) {
            return c;
        }
        undefined("read of a filter string byte never written");
        0
    }

    fn set(&mut self, i: usize, c: u8) {
        if i.checked_add(1).is_none() {
            undefined("filter string offset overflow");
            return;
        }
        if let Some(old) = self.buf.get_mut(i) {
            *old = c;
        } else if i == self.buf.len() {
            self.buf.push(c);
            while let Some(c) = self.sparse.as_mut().and_then(|s| s.remove(&self.buf.len())) {
                self.buf.push(c);
            }
        } else {
            self.sparse
                .get_or_insert_with(|| Box::new(HashMap::new()))
                .insert(i, c);
        }
    }

    fn puts(&mut self, s: &[u8]) {
        let Some(end) = self.l.checked_add(s.len()) else {
            undefined("filter string length overflow");
            return;
        };
        let Some(size) = end.checked_add(2) else {
            undefined("filter string length overflow");
            return;
        };
        self.resize(size);
        if self.sparse.is_none() && self.l <= self.buf.len() {
            // The dense prefix proves both the overwrite and append are
            // initialized; bytes beyond end retain C's reusable suffix.
            let overwrite = s.len().min(self.buf.len() - self.l);
            self.buf[self.l..self.l + overwrite].copy_from_slice(&s[..overwrite]);
            self.buf.extend_from_slice(&s[overwrite..]);
            self.l = end;
            self.set(end, 0);
            return;
        }
        for (k, &c) in s.iter().enumerate() {
            self.set(self.l + k, c);
        }
        self.l = end;
        self.set(end, 0);
    }

    fn putc(&mut self, c: u8) {
        let Some(end) = self.l.checked_add(1) else {
            undefined("filter string length overflow");
            return;
        };
        let Some(size) = end.checked_add(1) else {
            undefined("filter string length overflow");
            return;
        };
        self.resize(size);
        self.set(self.l, c);
        self.l = end;
        self.set(end, 0);
    }

    fn putc_(&mut self, c: u8) {
        let Some(end) = self.l.checked_add(1) else {
            undefined("filter string length overflow");
            return;
        };
        self.resize(end);
        self.set(self.l, c);
        self.l = end;
    }

    /// The C string at at: bytes up to the first NUL.
    fn cstr(&self, at: usize) -> Vec<u8> {
        let mut out = Vec::new();
        let mut i = at;
        loop {
            let c = self.get(i);
            if c == 0 {
                return out;
            }
            out.push(c);
            i += 1;
        }
    }

    fn cstrlen(&self, at: usize) -> usize {
        let mut i = at;
        while self.get(i) != 0 {
            i += 1;
        }
        i - at
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum Setter {
    Qual,
    Type,
    Info,
    InfoInt,
    InfoFloat,
    InfoString,
    InfoFlag,
    FormatInt,
    FormatFloat,
    FormatString,
    Genotype2,
    Genotype3,
    Genotype4,
    GenotypeString,
    Ilen,
    Ref,
    Alt,
    Nmissing,
    Nalt,
    Chrom,
    Pos,
    Ac,
    An,
    Mac,
    Af,
    Maf,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum Func {
    Npass,
    Max,
    Min,
    Avg,
    Median,
    Stddev,
    Sum,
    Abs,
    Count,
    Strlen,
    Phred,
    Binom,
    Fisher,
    Perl,
    SCount,
    SMax,
    SMin,
    SAvg,
    SMedian,
    SStddev,
    SSum,
    LogicOr,
    LogicAnd,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum Cmp {
    BitAnd,
    Filter,
    StringHash,
    Id,
}

/// filter.c's `token_t`. A default token is the zeroed one filter.c gets
/// from `hts_expand0`.
#[derive(Debug, Clone, Default)]
struct Token {
    tok_type: i32,
    nargs: i32,
    key: Option<Vec<u8>>,
    tag: Option<Vec<u8>>,
    threshold: f64,
    is_constant: bool,
    hdr_id: i32,
    hl_type: i32,
    ht_type: i32,
    vl_len: i32,
    idx: i32,
    idxs: Vec<i32>,
    nuidxs: i32,
    usmpl: Option<Vec<u8>>,
    nsamples: i32,
    setter: Option<Setter>,
    func: Option<Func>,
    comparator: Option<Cmp>,
    hash: Option<HashSet<Vec<u8>>>,
    regex: Option<Box<Regex>>,
    iext: i32,
    values: Vals,
    str_value: KStr,
    is_str: bool,
    is_missing: bool,
    pass_site: i32,
    pass_samples: Vec<u8>,
    nvalues: i32,
    nval1: i32,
}

impl Token {
    fn nidxs(&self) -> i32 {
        self.idxs.len() as i32
    }

    /// A value token as `filters_init1` starts it.
    fn value() -> Token {
        Token {
            tok_type: TOK_VAL,
            hdr_id: -1,
            hl_type: -1,
            ht_type: -1,
            vl_len: VL_FIXED,
            pass_site: -1,
            ..Token::default()
        }
    }

    fn usmpl(&self, i: usize) -> bool {
        match &self.usmpl {
            Some(u) => match u.get(i) {
                Some(&v) => v != 0,
                None => {
                    undefined("sample mask read past its end");
                    false
                }
            },
            None => {
                undefined("sample mask is NULL");
                false
            }
        }
    }

    fn pass_sample(&self, i: usize) -> u8 {
        match self.pass_samples.get(i) {
            Some(&v) => v,
            None => {
                undefined("sample status read past its end");
                0
            }
        }
    }

    fn set_pass_sample(&mut self, i: usize, v: u8) {
        match self.pass_samples.get_mut(i) {
            Some(x) => *x = v,
            None => undefined("sample status written past its end"),
        }
    }
}

/// filter.c's genotype cache (`cached_GT`).
#[derive(Debug, Clone, Default)]
struct CachedGt {
    valid: bool,
    nbuf: i32,
    buf: Vec<i32>,
    mask: Option<Vec<u64>>,
}

/// filter.c's allele-count cache (`cached_site`).
#[derive(Debug, Clone, Default)]
struct CachedSite {
    valid: bool,
    n_missing: i32,
    n_samples: i32,
    an: i32,
    ac: Vec<i32>,
    af: Vec<f64>,
}

/// A compiled filter expression.
#[derive(Debug, Clone)]
pub struct Filter {
    str: Vec<u8>,
    filters: Vec<Token>,
    stack: Vec<usize>,
    nsamples: i32,
    /// `flt->tmpi`/`tmpf`: shared integer and float buffers, with their
    /// stale contents (`None`: never written).
    tmpi: Vec<Option<i32>>,
    tmpf: Vec<Option<u32>>,
    tmps: KStr,
    cached_gt: CachedGt,
    cached_site: CachedSite,
}

/// Header lookups as filter.c makes them.
struct Hdr<'a>(&'a Header);

fn kind(hl: i32) -> Option<LineKind> {
    match hl {
        HL_FLT => Some(LineKind::Filter),
        HL_INFO => Some(LineKind::Info),
        HL_FMT => Some(LineKind::Format),
        _ => None,
    }
}

impl Hdr<'_> {
    /// `bcf_hdr_id2int(hdr, BCF_DT_ID, name)`.
    fn id(&self, name: &[u8]) -> i32 {
        std::str::from_utf8(name)
            .ok()
            .and_then(|n| self.0.key_id(n))
            .map_or(-1, |i| i as i32)
    }

    /// `bcf_hdr_idinfo_exists`.
    fn exists(&self, hl: i32, id: i32) -> bool {
        id >= 0 && kind(hl).is_some_and(|k| self.0.declared_as(id as u32, k))
    }

    /// The declaration word; reading it for id -1 is undefined in C.
    fn word(&self, hl: i32, id: i32) -> Result<u32> {
        let k = kind(hl).ok_or_else(|| unsupported("header lookup of an invalid line kind"))?;
        if id < 0 {
            return Err(unsupported("header lookup of an undefined tag"));
        }
        self.0
            .declaration_word(id as u32, k)
            .ok_or_else(|| unsupported("header lookup of a removed tag"))
    }

    fn id2type(&self, hl: i32, id: i32) -> Result<i32> {
        Ok((self.word(hl, id)? >> 4 & 0xf) as i32)
    }

    fn id2length(&self, hl: i32, id: i32) -> Result<i32> {
        Ok((self.word(hl, id)? >> 8 & 0xf) as i32)
    }

    fn id2number(&self, hl: i32, id: i32) -> Result<i32> {
        Ok((self.word(hl, id)? >> 12) as i32)
    }

    fn nsamples(&self) -> i32 {
        self.0.samples.len() as i32
    }
}

fn at(s: &[u8], i: usize) -> u8 {
    s.get(i).copied().unwrap_or(0)
}

/// `strncasecmp(s, word, n) == 0` on C strings.
fn ncase_eq(s: &[u8], word: &[u8], n: usize) -> bool {
    for k in 0..n {
        let (a, b) = (at(s, k), at(word, k));
        if !a.eq_ignore_ascii_case(&b) {
            return false;
        }
        if a == 0 {
            return true;
        }
    }
    true
}

/// `strncmp(s, word, n) == 0` on C strings.
fn n_eq(s: &[u8], word: &[u8], n: usize) -> bool {
    for k in 0..n {
        let (a, b) = (at(s, k), at(word, k));
        if a != b {
            return false;
        }
        if a == 0 {
            return true;
        }
    }
    true
}

fn case_eq(a: &[u8], b: &[u8]) -> bool {
    a.eq_ignore_ascii_case(b)
}

/// The C string at the start of `s`.
fn cstr(s: &[u8]) -> &[u8] {
    &s[..s.iter().position(|&c| c == 0).unwrap_or(s.len())]
}

/// glibc `strtol(s, &end, 10)`: (value, consumed). Overflow clamps to the
/// `long` range; no digits consume nothing.
fn strtol(s: &[u8]) -> (i64, usize) {
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
    let digits = i;
    let mut n: i128 = 0;
    while at(s, i).is_ascii_digit() {
        n = (n * 10 + i128::from(at(s, i) - b'0')).min(1 << 70);
        i += 1;
    }
    if i == digits {
        return (0, 0);
    }
    let n = if neg { -n } else { n };
    (n.clamp(i64::MIN as i128, i64::MAX as i128) as i64, i)
}

/// `filters_next_token`: the token at `*pos` (whitespace skipped first)
/// and its length; operators and function names advance `*pos` past
/// themselves as filter.c advances its string. Negative: a function
/// taking a variable number of arguments; -1: unterminated quotes.
fn next_token(s: &[u8], pos: &mut usize) -> (i32, usize) {
    let mut tmp = *pos;
    while at(s, tmp) != 0 && is_space(at(s, tmp)) {
        tmp += 1;
    }
    *pos = tmp;
    let c0 = at(s, *pos);
    if c0.is_ascii_digit() || c0 == b'.' {
        let (_, n, _) = strtod(&s[*pos..]);
        if n > 0 && (at(s, *pos + n) == 0 || !at(s, *pos + n).is_ascii_alphanumeric()) {
            return (TOK_VAL, n);
        }
        tmp = *pos;
    }
    let ci = |p: usize, w: &[u8]| ncase_eq(&s[p.min(s.len())..], w, w.len());
    const FUNCS: &[(&[u8], usize, i32)] = &[
        (b"SMPL_COUNT(", 10, TOK_SCOUNT),
        (b"SMPL_MAX(", 8, TOK_SMAX),
        (b"SMPL_MIN(", 8, TOK_SMIN),
        (b"SMPL_MEAN(", 9, TOK_SAVG),
        (b"SMPL_MEDIAN(", 11, TOK_SMEDIAN),
        (b"SMPL_AVG(", 8, TOK_SAVG),
        (b"SMPL_STDEV(", 10, TOK_SSTDEV),
        (b"SMPL_SUM(", 8, TOK_SSUM),
        (b"sCOUNT(", 6, TOK_SCOUNT),
        (b"sMAX(", 4, TOK_SMAX),
        (b"sMIN(", 4, TOK_SMIN),
        (b"sMEAN(", 5, TOK_SAVG),
        (b"sMEDIAN(", 7, TOK_SMEDIAN),
        (b"sAVG(", 4, TOK_SAVG),
        (b"sSTDEV(", 6, TOK_SSTDEV),
        (b"sSUM(", 4, TOK_SSUM),
        (b"MAX(", 3, TOK_MAX),
        (b"MIN(", 3, TOK_MIN),
        (b"MEAN(", 4, TOK_AVG),
        (b"MEDIAN(", 6, TOK_MEDIAN),
        (b"AVG(", 3, TOK_AVG),
        (b"STDEV(", 5, TOK_STDEV),
        (b"SUM(", 3, TOK_SUM),
        (b"ABS(", 3, TOK_ABS),
        (b"COUNT(", 5, TOK_CNT),
        (b"STRLEN(", 6, TOK_LEN),
        (b"BINOM(", 5, -TOK_BINOM),
        (b"FISHER(", 6, -TOK_FISHER),
        (b"PHRED(", 5, TOK_PHRED),
        (b"%MAX(", 4, TOK_MAX),
        (b"%MIN(", 4, TOK_MIN),
        (b"%AVG(", 4, TOK_AVG),
        (b"%SUM(", 4, TOK_SUM),
    ];
    for &(name, adv, tok) in FUNCS {
        // FISHER is matched on its first six characters, as filter.c does.
        let n = if tok == -TOK_FISHER { 6 } else { name.len() };
        if ncase_eq(&s[tmp.min(s.len())..], name, n) {
            *pos += adv;
            return (tok, 0);
        }
    }
    if ci(tmp, b"INFO/") {
        tmp += 5;
    }
    if ci(tmp, b"FORMAT/") {
        tmp += 7;
    }
    if ci(tmp, b"FMT/") {
        tmp += 4;
    }
    if ci(tmp, b"PERL.") {
        *pos += 5;
        return (-TOK_PERLSUB, 0);
    }
    if ci(tmp, b"N_PASS(") || ci(tmp, b"F_PASS(") {
        *pos += 6;
        return (-TOK_FUNC, 6);
    }
    if ci(tmp, b"%ILEN") {
        return (TOK_VAL, 5);
    }
    for (w, n) in [
        (&b"{}"[..], 2),
        (b"{STR}", 5),
        (b"{INT}", 5),
        (b"{FLOAT}", 7),
    ] {
        if ci(tmp, w) {
            return (TOK_EXT, n);
        }
    }
    if at(s, tmp) == b'@' {
        while at(s, tmp) != 0 && !is_space(at(s, tmp)) && at(s, tmp) != b'=' && at(s, tmp) != b'!' {
            tmp += 1;
        }
        return (TOK_VAL, tmp - *pos);
    }
    let mut brackets = 0;
    while at(s, tmp) != 0 {
        let c = at(s, tmp);
        if brackets == 0
            && (matches!(
                c,
                b'"' | b'\''
                    | b'<'
                    | b'>'
                    | b'='
                    | b'!'
                    | b'&'
                    | b'|'
                    | b'('
                    | b')'
                    | b'+'
                    | b'*'
                    | b'-'
                    | b'/'
                    | b'~'
                    | b'%'
            ) || is_space(c))
        {
            break;
        }
        if c == b']' {
            if brackets != 0 {
                tmp += 1;
            }
            break;
        }
        if c == b'[' {
            brackets += 1;
        }
        tmp += 1;
    }
    if tmp > *pos {
        return (TOK_VAL, tmp - *pos);
    }
    let c = at(s, tmp);
    if c == b'"' || c == b'\'' {
        tmp += 1;
        while at(s, tmp) != 0 && at(s, tmp) != c {
            tmp += 1;
        }
        if at(s, tmp) == 0 {
            return (-1, 0);
        }
        return (TOK_VAL, tmp - *pos + 1);
    }
    let c1 = at(s, tmp + 1);
    let (tok, n) = match c {
        b'!' if c1 == b'=' => (TOK_NE, 2),
        b'!' if c1 == b'~' => (TOK_NLIKE, 2),
        b'<' if c1 == b'=' => (TOK_LE, 2),
        b'<' => (TOK_LT, 1),
        b'>' if c1 == b'=' => (TOK_BE, 2),
        b'>' => (TOK_BT, 1),
        b'=' if c1 == b'=' => (TOK_EQ, 2),
        b'=' => (TOK_EQ, 1),
        b'(' => (TOK_LFT, 1),
        b')' => (TOK_RGT, 1),
        b'&' if c1 == b'&' => (TOK_AND_VEC, 2),
        b'|' if c1 == b'|' => (TOK_OR_VEC, 2),
        b'&' => (TOK_AND, 1),
        b'|' => (TOK_OR, 1),
        b'+' => (TOK_ADD, 1),
        b'-' => (TOK_SUB, 1),
        b'*' => (TOK_MULT, 1),
        b'/' => (TOK_DIV, 1),
        b'~' => (TOK_LIKE, 1),
        b'%' => (TOK_MODULO, 1),
        _ => return (TOK_VAL, tmp - *pos),
    };
    *pos += n;
    (tok, 0)
}

/// `parse_tag_list`: splits on commas outside brackets.
fn parse_tag_list(s: &[u8]) -> Vec<Vec<u8>> {
    let mut out = Vec::new();
    let mut q = 0;
    let mut nest = 0;
    let mut p = 0;
    loop {
        let c = at(s, p);
        if (c == b',' && nest == 0) || c == 0 {
            out.push(s[q..p].to_vec());
            q = p + 1;
        }
        if c == 0 {
            break;
        }
        if c == b'[' || c == b'(' {
            nest += 1;
        }
        if (c == b']' || c == b')') && nest > 0 {
            nest -= 1;
        }
        p += 1;
    }
    out
}

/// `expand_path`: `~`, `~/...` and `$VAR`. A `~user` prefix is kept as
/// written (the password database is not consulted).
fn expand_path(path: &[u8]) -> Vec<u8> {
    use std::os::unix::ffi::OsStrExt;
    if at(path, 0) == b'~' {
        if at(path, 1) == 0 || at(path, 1) == b'/' {
            let mut out = match std::env::var_os("HOME") {
                Some(h) => h.as_bytes().to_vec(),
                None => b".".to_vec(),
            };
            out.extend_from_slice(&path[1..]);
            return out;
        }
        return path.to_vec();
    }
    if at(path, 0) == b'$' {
        if let Some(v) = std::env::var_os(std::ffi::OsStr::from_bytes(&path[1..])) {
            return v.as_bytes().to_vec();
        }
    }
    path.to_vec()
}

/// `hts_readlist(name, 1, &n)`: the non-empty lines of a plain or gzip
/// file, a trailing CR removed, each up to its first NUL.
fn read_list(name: &[u8]) -> Option<Vec<Vec<u8>>> {
    use std::os::unix::ffi::OsStrExt;
    let path = std::path::Path::new(std::ffi::OsStr::from_bytes(name));
    let raw = std::fs::read(path).ok()?;
    let data = if raw.starts_with(&[0x1f, 0x8b]) {
        let mut out = Vec::new();
        flate2::read::MultiGzDecoder::new(&raw[..])
            .read_to_end(&mut out)
            .ok()?;
        out
    } else {
        raw
    };
    let mut out = Vec::new();
    let mut lines = data.split(|&c| c == b'\n').peekable();
    while let Some(line) = lines.next() {
        if lines.peek().is_none() && line.is_empty() {
            break;
        }
        let line = line.strip_suffix(b"\r").unwrap_or(line);
        if line.is_empty() {
            continue;
        }
        out.push(cstr(line).to_vec());
    }
    Some(out)
}

/// `parse_idxs`: (idxs, idx), or `None` for an unparsable list.
fn parse_idxs(s: &[u8]) -> Result<Option<(Vec<i32>, i32)>> {
    if s.is_empty() || s == b"*" {
        return Ok(Some((vec![-1], -2)));
    }
    if s == b"GT" {
        return Ok(Some((vec![-1], -3)));
    }
    let (v, n) = strtol(s);
    let idx = v as i32;
    if idx >= 0 && n == s.len() {
        return Ok(Some((Vec::new(), idx)));
    }
    let mut idxs: Vec<i32> = Vec::new();
    let mut ibeg: i32 = -1;
    let mut beg = 0;
    while beg < s.len() {
        let (v, n) = strtol(&s[beg..]);
        let num = v as i32;
        let end = beg + n;
        match at(s, end) {
            b',' => beg = end + 1,
            0 => beg = end,
            b'-' => {
                beg = end + 1;
                ibeg = num;
                continue;
            }
            _ => return Ok(None),
        }
        if num < 0 {
            return Err(unsupported("negative index in a filter subscript"));
        }
        if num > 1 << 24 {
            return Err(unsupported("filter subscript too large to allocate"));
        }
        if num as usize >= idxs.len() {
            idxs.resize(num as usize + 1, 0);
        }
        if ibeg >= 0 {
            for i in ibeg..=num {
                idxs[i as usize] = 1;
            }
            ibeg = -1;
        }
        idxs[num as usize] = 1;
    }
    if ibeg >= 0 {
        if ibeg > 1 << 24 {
            return Err(unsupported("filter subscript too large to allocate"));
        }
        if ibeg as usize >= idxs.len() {
            idxs.resize(ibeg as usize + 1, 0);
        }
        idxs[ibeg as usize] = -1;
    }
    Ok(Some((idxs, -2)))
}

/// The state `filter_init` builds tokens against.
struct Init<'a> {
    hdr: Hdr<'a>,
    str: &'a [u8],
    max_unpack: i32,
    n_ext: i32,
    ext: Vec<i32>,
}

impl Init<'_> {
    fn init_usmpl(&self, tok: &mut Token) {
        if tok.nsamples != 0 {
            return;
        }
        tok.nsamples = self.hdr.nsamples();
        tok.usmpl = Some(vec![1; tok.nsamples as usize]);
    }

    /// `parse_tag_idx`: the subscript of `tag[...]`.
    fn parse_tag_idx(
        &self,
        is_fmt: bool,
        tag: &[u8],
        tag_idx: &[u8],
        tok: &mut Token,
    ) -> Result<()> {
        let mut tag_idx = tag_idx.to_vec();
        if tag_idx.last() == Some(&b']') {
            tag_idx.pop();
        }
        let ori = String::from_utf8_lossy(&tag_idx).into_owned();
        let nsmpl = self.hdr.nsamples();
        let mut idx1 = 0;
        let mut idxs1: Vec<i32> = Vec::new();
        let mut set_samples = false;
        let colon = tag_idx.iter().rposition(|&c| c == b':');
        if tag_idx.first() == Some(&b'@') {
            if !is_fmt {
                return Err(init_error(format!(
                    "Could not parse \"{ori}\". (Not a FORMAT tag yet a sample list provided.)"
                )));
            }
            let mut fname = expand_path(&tag_idx[1..]);
            let mut list = read_list(&fname);
            if list.is_none() {
                if let Some(c) = colon {
                    let (idxs2, idx2) = parse_idxs(&tag_idx[c + 1..])?
                        .ok_or_else(|| init_error(format!("Could not parse the index: {ori}")))?;
                    tok.idxs = idxs2;
                    tok.idx = idx2;
                    match fname.iter().rposition(|&c| c == b':') {
                        Some(fc) => fname.truncate(fc),
                        None => return Err(unsupported("sample list name without its colon")),
                    }
                    list = read_list(&fname);
                }
            }
            let list = list.ok_or_else(|| {
                init_error(format!(
                    "Could not read: {}",
                    String::from_utf8_lossy(&fname)
                ))
            })?;
            tok.nsamples = nsmpl;
            let mut usmpl = vec![0u8; nsmpl as usize];
            for name in &list {
                let id = std::str::from_utf8(name)
                    .ok()
                    .and_then(|n| self.hdr.0.sample_id(n))
                    .ok_or_else(|| {
                        init_error(format!(
                            "No such sample in the VCF: \"{}\"",
                            String::from_utf8_lossy(name)
                        ))
                    })?;
                usmpl[id] = 1;
            }
            tok.usmpl = Some(usmpl);
            if colon.is_none() {
                tok.idxs = vec![-1];
                tok.idx = -2;
            }
        } else if let Some(c) = colon {
            if !is_fmt {
                return Err(init_error(format!(
                    "Could not parse the index \"{ori}\". (Not a FORMAT tag yet sample index implied.)"
                )));
            }
            let bad = || init_error(format!("Could not parse the index: {ori}"));
            let (i1, x1) = parse_idxs(&tag_idx[..c])?.ok_or_else(bad)?;
            let (i2, x2) = parse_idxs(&tag_idx[c + 1..])?.ok_or_else(bad)?;
            idxs1 = i1;
            idx1 = x1;
            tok.idxs = i2;
            tok.idx = x2;
            set_samples = true;
        } else {
            let (i1, x1) = parse_idxs(&tag_idx)?
                .ok_or_else(|| init_error(format!("Could not parse the index: {ori}")))?;
            if is_fmt {
                if i1.len() == 1 && i1[0] == -1 {
                    tok.idxs = vec![-1];
                    tok.idx = x1;
                } else if self.hdr.id2number(HL_FMT, tok.hdr_id)? != 1 {
                    let t = String::from_utf8_lossy(tag);
                    return Err(init_error(format!(
                        "The FORMAT tag {t} can have multiple subfields, run as {t}[sample:subfield]"
                    )));
                } else {
                    tok.idx = 0;
                }
                idxs1 = i1;
                idx1 = x1;
                set_samples = true;
            } else {
                tok.idxs = i1;
                tok.idx = x1;
            }
        }
        if set_samples {
            tok.nsamples = nsmpl;
            let mut usmpl = vec![0u8; nsmpl as usize];
            if idx1 >= 0 {
                if idx1 >= nsmpl {
                    return Err(init_error(format!("The sample index is too large: {ori}")));
                }
                usmpl[idx1 as usize] = 1;
            } else if idx1 == -2 || idx1 == -3 {
                let mut i = 0;
                while i < idxs1.len() {
                    if idxs1[i] == 0 {
                        i += 1;
                        continue;
                    }
                    if idxs1[i] == -1 {
                        break;
                    }
                    if i as i32 >= nsmpl {
                        return Err(init_error(format!("The sample index is too large: {ori}")));
                    }
                    usmpl[i] = 1;
                    i += 1;
                }
                if idxs1.last() == Some(&-1) {
                    while (i as i32) < nsmpl {
                        usmpl[i] = 1;
                        i += 1;
                    }
                }
            } else {
                return Err(init_error("todo: filter.c sample index"));
            }
            tok.usmpl = Some(usmpl);
        }
        if !tok.idxs.is_empty() && *tok.idxs.last().unwrap() != -1 {
            tok.nuidxs += tok.idxs.iter().filter(|&&v| v != 0).count() as i32;
        }
        Ok(())
    }

    /// `filters_init1_ext`: `{}`, `{str}`, `{int}`, `{float}`.
    fn init1_ext(&mut self, s: &[u8], len: usize) -> Token {
        let mut tok = Token::value();
        self.n_ext += 1;
        tok.iext = self.n_ext;
        if ncase_eq(s, b"{str}", len) {
            tok.ht_type = HT_STR;
            tok.is_str = true;
        } else if ncase_eq(s, b"{int}", len) {
            tok.ht_type = HT_INT;
        } else if ncase_eq(s, b"{float}", len) {
            tok.ht_type = HT_REAL;
        }
        self.ext.push(tok.ht_type);
        tok
    }

    /// `filters_init1`: a value token from `s[..len]` (`s` continues as
    /// the C string does).
    fn init1(&mut self, s: &[u8], len: usize) -> Result<Token> {
        let mut tok = Token::value();
        let c0 = at(s, 0);
        if c0 == b'"' || c0 == b'\'' {
            if at(s, len - 1) != c0 {
                return Err(init_error(format!(
                    "TODO: [{}]",
                    String::from_utf8_lossy(self.str)
                )));
            }
            let key = s[1..len - 1].to_vec();
            tok.values.expand(1);
            tok.values.set(0, (len - 2) as f64);
            tok.is_missing = key == b".";
            tok.nvalues = (len - 2) as i32;
            tok.key = Some(key);
            tok.is_str = true;
            tok.ht_type = HT_STR;
            return Ok(tok);
        }
        if c0 == b'@' {
            let tag = s[..len].to_vec();
            let fname = expand_path(cstr(&tag[1..]));
            let list = read_list(&fname).ok_or_else(|| {
                init_error(format!(
                    "Could not read: {}",
                    String::from_utf8_lossy(&fname)
                ))
            })?;
            let mut hash = HashSet::new();
            for item in list {
                let end = item.iter().position(|&c| is_space(c)).unwrap_or(item.len());
                hash.insert(item[..end].to_vec());
            }
            tok.tag = Some(tag);
            tok.hash = Some(hash);
            return Ok(tok);
        }

        let (mut s, mut len) = (s, len);
        let mut is_fmt: i32 = -1;
        if ncase_eq(s, b"FMT/", 4) {
            s = &s[4..];
            len -= 4;
            is_fmt = 1;
        } else if ncase_eq(s, b"FORMAT/", 7) {
            s = &s[7..];
            len -= 7;
            is_fmt = 1;
        } else {
            if ncase_eq(s, b"INFO/", 5) {
                is_fmt = 0;
                s = &s[5..];
                len -= 5;
            } else if ncase_eq(s, b"QUAL", len) || n_eq(s, b"%QUAL", len) {
                tok.setter = Some(Setter::Qual);
                tok.tag = Some(b"QUAL".to_vec());
                tok.ht_type = HT_REAL;
                return Ok(tok);
            } else if ncase_eq(s, b"TYPE", len) || n_eq(s, b"%TYPE", len) {
                tok.setter = Some(Setter::Type);
                tok.tag = Some(b"TYPE".to_vec());
                tok.ht_type = HT_STR;
                return Ok(tok);
            } else if ncase_eq(s, b"FILTER", len) || n_eq(s, b"%FILTER", len) {
                tok.comparator = Some(Cmp::Filter);
                tok.tag = Some(b"FILTER".to_vec());
                self.max_unpack |= UN_FLT;
                tok.hl_type = HL_FLT;
                tok.ht_type = HT_STR;
                return Ok(tok);
            } else if ncase_eq(s, b"ID", len) || ncase_eq(s, b"%ID", len) {
                tok.comparator = Some(Cmp::Id);
                tok.tag = Some(b"ID".to_vec());
                tok.ht_type = HT_STR;
                return Ok(tok);
            } else if ncase_eq(s, b"CHROM", len) {
                tok.setter = Some(Setter::Chrom);
                tok.tag = Some(b"CHROM".to_vec());
                tok.ht_type = HT_STR;
                return Ok(tok);
            } else if ncase_eq(s, b"POS", len) {
                tok.setter = Some(Setter::Pos);
                tok.tag = Some(b"POS".to_vec());
                tok.ht_type = HT_INT;
                return Ok(tok);
            } else if ncase_eq(s, b"REF", len) {
                tok.setter = Some(Setter::Ref);
                tok.is_str = true;
                tok.tag = Some(b"REF".to_vec());
                tok.ht_type = HT_STR;
                return Ok(tok);
            } else if ncase_eq(s, b"ALT", len) {
                tok.setter = Some(Setter::Alt);
                tok.is_str = true;
                tok.tag = Some(b"ALT".to_vec());
                tok.ht_type = HT_STR;
                tok.idxs = vec![-1];
                tok.idx = -2;
                return Ok(tok);
            } else if ncase_eq(s, b"N_ALT", len) {
                tok.setter = Some(Setter::Nalt);
                tok.tag = Some(b"N_ALT".to_vec());
                tok.ht_type = HT_INT;
                return Ok(tok);
            } else if ncase_eq(s, b"N_SAMPLES", len) {
                tok.tok_type = TOK_VAL;
                tok.threshold = f64::from(self.hdr.nsamples());
                tok.is_constant = true;
                tok.ht_type = HT_INT;
                return Ok(tok);
            } else if ncase_eq(s, b"N_MISSING", len) || ncase_eq(s, b"F_MISSING", len) {
                let n = ncase_eq(s, b"N_MISSING", len);
                self.max_unpack |= UN_FMT;
                tok.setter = Some(Setter::Nmissing);
                tok.tag = Some(if n {
                    b"N_MISSING".to_vec()
                } else {
                    b"F_MISSING".to_vec()
                });
                tok.ht_type = if n { HT_INT } else { HT_REAL };
                self.init_usmpl(&mut tok);
                return Ok(tok);
            }
        }

        let mut name = s[..len].to_vec();
        let mut is_array = 0;
        if name.last() == Some(&b']') {
            if let Some(i) = name.iter().position(|&c| c == b'[') {
                is_array = i + 1;
            }
        }
        let idx_text = if is_array > 0 {
            let t = name[is_array..].to_vec();
            name.truncate(is_array - 1);
            t
        } else {
            Vec::new()
        };
        // The name is a C string from here on.
        let name = cstr(&name).to_vec();
        tok.hdr_id = self.hdr.id(&name);
        if is_fmt == -1 {
            if tok.hdr_id >= 0 {
                let is_info = self.hdr.exists(HL_INFO, tok.hdr_id);
                let fmt = self.hdr.exists(HL_FMT, tok.hdr_id);
                is_fmt = i32::from(fmt);
                if is_info && fmt {
                    let n = String::from_utf8_lossy(&name);
                    return Err(init_error(format!(
                        "Error: ambiguous filtering expression, both INFO/{n} and FORMAT/{n} are defined in the VCF header."
                    )));
                }
            }
            if is_fmt == -1 {
                is_fmt = 0;
            }
        }
        let is_fmt = is_fmt == 1;
        if is_array > 0 {
            self.parse_tag_idx(is_fmt, &name, &idx_text, &mut tok)?;
            if tok.idx == -3 && self.hdr.id2length(HL_FMT, tok.hdr_id)? != VL_R {
                return Err(init_error(
                    "Error: GT subscripts can be used only with Number=R tags",
                ));
            }
        } else if is_fmt {
            self.init_usmpl(&mut tok);
        }
        tok.hl_type = if is_fmt { HL_FMT } else { HL_INFO };
        if tok.hdr_id >= 0 {
            tok.vl_len = self.hdr.id2length(tok.hl_type, tok.hdr_id)?;
        }
        if is_fmt {
            self.max_unpack |= UN_FMT;
        }
        let shown = || String::from_utf8_lossy(&name).into_owned();
        if tok.hdr_id >= 0 {
            if is_fmt && name == b"GT" {
                tok.setter = Some(Setter::GenotypeString);
                tok.is_str = true;
                tok.ht_type = HT_STR;
            } else if is_fmt {
                if !self.hdr.exists(HL_FMT, tok.hdr_id) {
                    return Err(init_error(format!("No such FORMAT field: {}", shown())));
                }
                if self.hdr.id2number(HL_FMT, tok.hdr_id)? != 1 && is_array == 0 {
                    tok.idxs = vec![-1];
                    tok.idx = -2;
                }
                match self.hdr.id2type(HL_FMT, tok.hdr_id)? {
                    HT_INT => {
                        tok.setter = Some(Setter::FormatInt);
                        tok.ht_type = HT_INT;
                    }
                    HT_REAL => {
                        tok.setter = Some(Setter::FormatFloat);
                        tok.ht_type = HT_REAL;
                    }
                    HT_STR => {
                        tok.setter = Some(Setter::FormatString);
                        tok.ht_type = HT_STR;
                        tok.is_str = true;
                    }
                    _ => return Err(init_error("[filter.c] FIXME")),
                }
            } else if !self.hdr.exists(HL_INFO, tok.hdr_id) {
                return Err(init_error(format!("No such INFO field: {}", shown())));
            } else {
                let ty = self.hdr.id2type(HL_INFO, tok.hdr_id)?;
                if ty == HT_FLAG {
                    tok.setter = Some(Setter::InfoFlag);
                    tok.ht_type = HT_INT;
                } else {
                    tok.ht_type = ty;
                    if ty == HT_STR {
                        tok.is_str = true;
                    }
                    if self.hdr.id2number(HL_INFO, tok.hdr_id)? == 1 {
                        tok.setter = Some(Setter::Info);
                    } else {
                        tok.setter = Some(match ty {
                            HT_INT => Setter::InfoInt,
                            HT_REAL => Setter::InfoFloat,
                            _ => {
                                tok.is_str = true;
                                Setter::InfoString
                            }
                        });
                        if is_array == 0 {
                            tok.idx = -2;
                            tok.idxs = vec![-1];
                        }
                    }
                }
                self.max_unpack |= UN_INFO;
            }
            tok.tag = Some(name);
            return Ok(tok);
        }
        if case_eq(&name, b"ALT") {
            tok.setter = Some(Setter::Alt);
            tok.is_str = true;
            tok.ht_type = HT_STR;
            tok.tag = Some(name);
            return Ok(tok);
        }
        if case_eq(&name, b"AN") {
            self.max_unpack |= UN_FMT;
            tok.setter = Some(Setter::An);
            tok.tag = Some(b"AN".to_vec());
            tok.ht_type = HT_INT;
            return Ok(tok);
        }
        if case_eq(&name, b"AC") {
            self.max_unpack |= UN_FMT;
            tok.setter = Some(Setter::Ac);
            tok.tag = Some(b"AC".to_vec());
            tok.idx = -2;
            tok.ht_type = HT_INT;
            return Ok(tok);
        }
        if case_eq(&name, b"MAC") || case_eq(&name, b"AF") || case_eq(&name, b"MAF") {
            self.max_unpack |= self.max_ac_an_unpack();
            let (setter, tag, ht) = if case_eq(&name, b"MAC") {
                (Setter::Mac, &b"MAC"[..], HT_INT)
            } else if case_eq(&name, b"AF") {
                tok.idx = -2;
                (Setter::Af, &b"AF"[..], HT_REAL)
            } else {
                (Setter::Maf, &b"MAF"[..], HT_REAL)
            };
            tok.setter = Some(setter);
            tok.tag = Some(tag.to_vec());
            tok.ht_type = ht;
            return Ok(tok);
        }
        if case_eq(&name, b"ILEN") || case_eq(&name, b"%ILEN") {
            self.max_unpack |= UN_STR;
            tok.setter = Some(Setter::Ilen);
            tok.tag = Some(b"ILEN".to_vec());
            tok.ht_type = HT_INT;
            return Ok(tok);
        }

        // A number: an integer if strtol takes it all, else a double.
        let (v, n) = strtol(&name);
        tok.threshold = v as f64;
        if n != name.len() {
            let (d, n, erange) = strtod(&name);
            tok.threshold = d;
            if erange || n != len {
                return Err(init_error(format!(
                    "[filter.c] Error: the tag \"{}\" is not defined in the VCF header",
                    shown()
                )));
            }
            tok.ht_type = HT_REAL;
        } else {
            tok.ht_type = HT_INT;
        }
        tok.is_constant = true;
        Ok(tok)
    }

    fn max_ac_an_unpack(&self) -> i32 {
        for name in [&b"AC"[..], b"AN"] {
            let id = self.hdr.id(name);
            if id < 0 || !self.hdr.exists(HL_INFO, id) {
                return UN_FMT;
            }
        }
        UN_INFO
    }
}

/// A token of the compiled list at `i`, as filter.c reads `out[i]`: past
/// the end but within the allocation it is a zeroed token; elsewhere the
/// read is undefined.
fn peek(out: &[Token], i: isize) -> Result<&Token> {
    static ZERO: std::sync::OnceLock<Token> = std::sync::OnceLock::new();
    let cap = out.len().next_power_of_two() as isize;
    if i >= 0 && (i as usize) < out.len() {
        Ok(&out[i as usize])
    } else if i >= 0 && i < cap {
        Ok(ZERO.get_or_init(Token::default))
    } else {
        Err(unsupported("filter expression read outside its token list"))
    }
}

impl Filter {
    /// Compiles `expr` against `header` as `filter_init` does; an error
    /// where filter.c reports one. Regex nesting or syntax trees deeper
    /// than 64 nodes return `Error::Unsupported` to bound recursive parsing,
    /// cloning and destruction.
    pub fn new(header: &Header, expr: &str) -> Result<Filter> {
        take_undefined();
        let r = Filter::init(header, expr.as_bytes());
        if let Some(what) = take_undefined() {
            return Err(unsupported(what));
        }
        r
    }

    fn init(header: &Header, expr: &[u8]) -> Result<Filter> {
        let s = cstr(expr);
        let shown = || String::from_utf8_lossy(s).into_owned();
        let mut init = Init {
            hdr: Hdr(header),
            str: s,
            max_unpack: UN_STR,
            n_ext: 0,
            ext: Vec::new(),
        };
        // perl_init: a "perl:" prefix needs Perl filters, which are not built.
        let mut tmp = 0;
        while at(s, tmp) != 0 && is_space(at(s, tmp)) {
            tmp += 1;
        }
        if at(s, tmp) != 0 && ncase_eq(&s[tmp..], b"perl:", 5) {
            return Err(init_error(
                "Perl filtering requires running `configure --enable-perl-filters` at compile time.",
            ));
        }
        let mut tmp = 0;
        let mut out: Vec<Token> = Vec::new();
        let mut ops: Vec<Token> = Vec::new();
        let mut last_op = -1;
        while at(s, tmp) != 0 {
            let (ret, len) = next_token(s, &mut tmp);
            if ret == -1 {
                return Err(init_error(format!("Missing quotes in: {}", shown())));
            }
            if ret == TOK_LFT {
                ops.push(Token {
                    tok_type: ret,
                    ..Token::default()
                });
            } else if ret == TOK_RGT {
                while ops.last().is_some_and(|t| t.tok_type != TOK_LFT) {
                    out.push(ops.pop().unwrap());
                }
                if ops.pop().is_none() {
                    return Err(init_error(format!("Could not parse: {}", shown())));
                }
            } else if ret == TOK_EXT {
                let tok = init.init1_ext(&s[tmp..], len);
                out.push(tok);
                tmp += len;
            } else if ret != TOK_VAL {
                let mut ret = ret;
                if ret == TOK_SUB && last_op != TOK_VAL && last_op != TOK_RGT {
                    // Unary minus: -x is -1*x.
                    out.push(Token {
                        tok_type: TOK_VAL,
                        hdr_id: -1,
                        pass_site: -1,
                        threshold: -1.0,
                        is_constant: true,
                        ..Token::default()
                    });
                    ret = TOK_MULT;
                } else if ret == -TOK_FUNC {
                    let start = tmp - len;
                    let mut tok = Token {
                        tok_type: TOK_FUNC,
                        hdr_id: -1,
                        pass_site: -1,
                        threshold: -1.0,
                        ..Token::default()
                    };
                    let name = &s[start..];
                    if ncase_eq(name, b"N_PASS", 6) || ncase_eq(name, b"F_PASS", 6) {
                        init.max_unpack |= UN_FMT;
                        tok.func = Some(Func::Npass);
                        tok.tag = Some(if ncase_eq(name, b"N_PASS", 6) {
                            b"N_PASS".to_vec()
                        } else {
                            b"F_PASS".to_vec()
                        });
                    } else {
                        return Err(init_error(format!(
                            "The function \"{}\" is not supported",
                            String::from_utf8_lossy(name)
                        )));
                    }
                    ops.push(tok);
                    continue;
                } else if ret < 0 {
                    let ret = -ret;
                    tmp += len;
                    let mut beg = tmp;
                    let mut nargs = 0;
                    if ret == TOK_PERLSUB {
                        while at(s, beg) != 0
                            && ((at(s, beg).is_ascii_alphanumeric()
                                && !at(s, beg).is_ascii_punctuation())
                                || at(s, beg) == b'_')
                        {
                            beg += 1;
                        }
                        if at(s, beg) != b'(' {
                            return Err(init_error(format!(
                                "[filter.c] Could not parse the expression: {}",
                                shown()
                            )));
                        }
                        let mut name = vec![b'"'];
                        name.extend_from_slice(&s[tmp..beg]);
                        name.push(b'"');
                        let n = name.len();
                        out.push(init.init1(&name, n)?);
                        nargs += 1;
                    }
                    let mut end = beg;
                    while at(s, end) != 0 && at(s, end) != b')' {
                        end += 1;
                    }
                    if at(s, end) == 0 {
                        return Err(init_error(format!(
                            "[filter.c] Could not parse the expression: {}",
                            shown()
                        )));
                    }
                    if end < beg + 1 {
                        // filter.c copies a negative length here and then
                        // reads a NULL string.
                        return Err(unsupported("filter function without an argument list"));
                    }
                    for arg in parse_tag_list(&s[beg + 1..end]) {
                        nargs += 1;
                        let n = arg.len();
                        out.push(init.init1(&arg, n)?);
                    }
                    out.push(Token {
                        tok_type: ret,
                        nargs,
                        hdr_id: -1,
                        pass_site: -1,
                        threshold: -1.0,
                        ..Token::default()
                    });
                    tmp = end + 1;
                    continue;
                } else {
                    while ops
                        .last()
                        .is_some_and(|t| OP_PREC[ret as usize] < OP_PREC[t.tok_type as usize])
                    {
                        out.push(ops.pop().unwrap());
                    }
                }
                ops.push(Token {
                    tok_type: ret,
                    ..Token::default()
                });
            } else if len == 0 {
                if at(s, tmp) != 0 && !is_space(at(s, tmp)) {
                    return Err(init_error(format!(
                        "Could not parse the expression: [{}]",
                        shown()
                    )));
                }
                break;
            } else {
                let n = if at(s, tmp + len - 1) == b',' {
                    len - 1
                } else {
                    len
                };
                let tok = init.init1(&s[tmp..], n)?;
                out.push(tok);
                tmp += len;
            }
            last_op = ret;
        }
        while let Some(tok) = ops.pop() {
            if tok.tok_type == TOK_LFT || tok.tok_type == TOK_RGT {
                return Err(init_error(format!(
                    "Could not parse the expression: [{}]",
                    shown()
                )));
            }
            out.push(tok);
        }
        if out.is_empty() {
            return Err(init_error(format!(
                "[filter.c] Error: could not parse the expression \"{}\"",
                shown()
            )));
        }
        determine_ext_types(&mut init, &mut out, &shown())?;
        postprocess(&init, &mut out, &shown())?;

        let nsamples = if init.max_unpack & UN_FMT != 0 {
            init.hdr.nsamples()
        } else {
            0
        };
        for tok in out.iter_mut() {
            let func = match tok.tok_type {
                TOK_MAX => Some(Func::Max),
                TOK_MIN => Some(Func::Min),
                TOK_AVG => Some(Func::Avg),
                TOK_MEDIAN => Some(Func::Median),
                TOK_STDEV => Some(Func::Stddev),
                TOK_SUM => Some(Func::Sum),
                TOK_ABS => Some(Func::Abs),
                TOK_CNT => Some(Func::Count),
                TOK_LEN => Some(Func::Strlen),
                TOK_PHRED => Some(Func::Phred),
                TOK_BINOM => Some(Func::Binom),
                TOK_FISHER => Some(Func::Fisher),
                TOK_PERLSUB => Some(Func::Perl),
                TOK_SCOUNT => Some(Func::SCount),
                TOK_SMAX => Some(Func::SMax),
                TOK_SMIN => Some(Func::SMin),
                TOK_SAVG => Some(Func::SAvg),
                TOK_SMEDIAN => Some(Func::SMedian),
                TOK_SSTDEV => Some(Func::SStddev),
                TOK_SSUM => Some(Func::SSum),
                _ => None,
            };
            if let Some(f) = func {
                tok.func = Some(f);
                tok.tok_type = TOK_FUNC;
            }
            tok.values.expand0(1);
            if nsamples > 0 {
                tok.pass_samples = vec![0; nsamples as usize];
            }
        }
        let n = out.len();
        Ok(Filter {
            str: s.to_vec(),
            filters: out,
            stack: vec![0; n],
            nsamples,
            tmpi: Vec::new(),
            tmpf: Vec::new(),
            tmps: KStr::default(),
            cached_gt: CachedGt::default(),
            cached_site: CachedSite::default(),
        })
    }
}

/// `determine_ext_types`: the type of each `{}` from the token before it.
fn determine_ext_types(init: &mut Init, out: &mut [Token], shown: &str) -> Result<()> {
    let n = out.len();
    let mut i = 0;
    while i < n {
        if out[i].iext == 0 || out[i].ht_type != -1 {
            i += 1;
            continue;
        }
        if i == 0 || i + 1 == n || out[i - 1].ht_type == -1 {
            break;
        }
        out[i].ht_type = out[i - 1].ht_type;
        i += 1;
    }
    if i != n {
        return Err(init_error(format!(
            "[filter.c] Error: unable to determine the type, use explicit notation: {shown}"
        )));
    }
    for tok in out.iter_mut() {
        let j = tok.iext - 1;
        if j < 0 {
            continue;
        }
        let e = &mut init.ext[j as usize];
        if *e != -1 && *e != tok.ht_type {
            return Err(init_error("[filter.c] FIXME: this should not happen"));
        }
        *e = tok.ht_type;
        if tok.ht_type == HT_STR {
            tok.is_str = true;
        }
    }
    Ok(())
}

/// The walk filter.c makes over the compiled list after parsing: FILTER
/// names to ids, TYPE and genotype keywords, regular expressions, `.` as
/// a missing number, `@file` comparators and the logic operators.
fn postprocess(init: &Init, out: &mut [Token], shown: &str) -> Result<()> {
    let hdr = &init.hdr;
    let n = out.len() as isize;
    let parse_err = || init_error(format!("Could not parse the expression: {shown}"));
    let deep_err = || {
        init_error(format!(
            "[filter.c] Could not parse the expression: {shown}"
        ))
    };
    let mut i: isize = 0;
    while i < n {
        let iu = i as usize;
        if i + 1 < n
            && (out[iu].tok_type == TOK_LT || out[iu].tok_type == TOK_BT)
            && out[iu + 1].tok_type == TOK_EQ
        {
            return Err(init_error(format!(
                "Error parsing the expression: \"{shown}\""
            )));
        }
        if out[iu].hash.is_some() {
            let j = if peek(out, i + 1)?.tok_type == TOK_VAL {
                i + 1
            } else {
                i - 1
            };
            let cmp = peek(out, j)?.comparator;
            if cmp != Some(Cmp::Id) {
                if cmp.is_some() {
                    return Err(init_error(
                        "Error: could not parse the expression with \"@file_name\" syntax (possible todo)",
                    ));
                }
                if (j as usize) < out.len() {
                    out[j as usize].comparator = Some(Cmp::StringHash);
                }
            }
        }
        match out[iu].tok_type {
            TOK_OR | TOK_OR_VEC => out[iu].func = Some(Func::LogicOr),
            TOK_AND | TOK_AND_VEC => out[iu].func = Some(Func::LogicAnd),
            _ => {}
        }
        if out[iu].tok_type == TOK_EQ || out[iu].tok_type == TOK_NE {
            let (mut j, mut k) = (i - 1, i - 2);
            if !peek(out, j)?.is_str {
                k = i - 1;
                j = i - 2;
            }
            let tj = peek(out, j)?;
            if tj.is_str && tj.key.as_deref() == Some(b".") {
                let tk = peek(out, k)?;
                let mut set_missing = false;
                if tk.hdr_id > 0 {
                    let ty = hdr.id2type(tk.hl_type, tk.hdr_id)?;
                    set_missing = ty == HT_INT || ty == HT_REAL;
                } else if tk.tag.is_none() {
                    return Err(init_error("Error: could not parse the expression"));
                } else if tk.tag.as_deref() == Some(b"QUAL") {
                    set_missing = true;
                }
                if set_missing && (j as usize) < out.len() {
                    let t = &mut out[j as usize];
                    t.is_str = false;
                    t.is_missing = true;
                    t.values.set(0, missing());
                }
            }
        }
        if out[iu].tok_type == TOK_LIKE || out[iu].tok_type == TOK_NLIKE {
            let j = i - 1;
            if peek(out, j)?.key.is_none() {
                return Err(init_error(format!(
                    "Could not parse the expression, wrong value for regex operator: {shown}"
                )));
            }
            let t = &mut out[j as usize];
            let key = t.key.as_mut().unwrap();
            let len = key.len();
            let mut icase = false;
            if len > 2 && key[len - 1] == b'i' && key[len - 2] == b'/' && key[len - 3] != b'\\' {
                key.truncate(len - 2);
                icase = true;
            }
            let re = Regex::new(key, icase).map_err(|_| {
                init_error(format!(
                    "Could not compile the regex expression \"{}\": {shown}",
                    String::from_utf8_lossy(key)
                ))
            })?;
            t.regex = Some(Box::new(re));
        }
        if out[iu].is_str
            && out[iu].tok_type == TOK_VAL
            && out[iu].key.as_ref().is_some_and(|k| k.contains(&b','))
        {
            // filter.c prints a note here; it reads the following tokens.
            let a = peek(out, i + 1)?;
            if a.tok_type != TOK_EQ && !(a.is_str && peek(out, i + 2)?.tok_type == TOK_EQ) {
                let _ = a.tok_type == TOK_NE || (a.is_str && peek(out, i + 2)?.tok_type == TOK_NE);
            }
        }
        if out[iu].tok_type != TOK_VAL || out[iu].tag.is_none() {
            i += 1;
            continue;
        }
        if out[iu].setter == Some(Setter::Type) {
            if i + 1 == n || peek(out, i + 1)?.key.is_none() {
                return Err(parse_err());
            }
            let t1 = peek(out, i + 1)?.tok_type;
            let t2 = peek(out, i + 2)?.tok_type;
            let (itok, ival) = if matches!(t1, TOK_EQ | TOK_NE | TOK_LIKE | TOK_NLIKE) {
                (i + 1, i - 1)
            } else if matches!(t2, TOK_EQ | TOK_NE | TOK_LIKE | TOK_NLIKE) {
                (i + 2, i + 1)
            } else {
                return Err(deep_err());
            };
            peek(out, ival)?;
            let tv = &mut out[ival as usize];
            let key = tv.key.clone().unwrap_or_default();
            let th = if case_eq(&key, b"snp") || case_eq(&key, b"snps") {
                2.0
            } else if case_eq(&key, b"indel") || case_eq(&key, b"indels") {
                8.0
            } else if case_eq(&key, b"mnp") || case_eq(&key, b"mnps") {
                4.0
            } else if case_eq(&key, b"other") {
                16.0
            } else if case_eq(&key, b"bnd") {
                32.0
            } else if case_eq(&key, b"overlap") {
                64.0
            } else if case_eq(&key, b"ref") {
                1.0
            } else {
                return Err(init_error(format!(
                    "The type \"{}\" not recognised: {shown}",
                    String::from_utf8_lossy(&key)
                )));
            };
            tv.threshold = th;
            tv.is_str = false;
            tv.is_constant = true;
            tv.tag = tv.key.take();
            if matches!(peek(out, itok)?.tok_type, TOK_LIKE | TOK_NLIKE)
                && (itok as usize) < out.len()
            {
                out[itok as usize].comparator = Some(Cmp::BitAnd);
            }
            i = itok + 1;
            continue;
        }
        if out[iu].tag.as_deref() == Some(b"GT") {
            if i + 1 == n {
                return Err(parse_err());
            }
            let t1 = peek(out, i + 1)?.tok_type;
            let ival = if matches!(t1, TOK_EQ | TOK_NE | TOK_LIKE | TOK_NLIKE) {
                i - 1
            } else if matches!(
                peek(out, i + 2)?.tok_type,
                TOK_EQ | TOK_NE | TOK_LIKE | TOK_NLIKE
            ) {
                i + 1
            } else {
                return Err(deep_err());
            };
            let Some(key) = peek(out, ival)?.key.clone() else {
                return Err(init_error(
                    "Comparison between samples is not supported, sorry!",
                ));
            };
            let lower = key.to_ascii_lowercase();
            let set = |s: Setter, k: Vec<u8>| (Some(s), k);
            let (setter, newkey) = match &lower[..] {
                b"hom" | b"het" | b"hap" => set(Setter::Genotype3, lower.clone()),
                b"mis" | b"ref" | b"alt" => set(Setter::Genotype4, lower.clone()),
                b"rr" => set(Setter::Genotype2, lower.clone()),
                b"ra" | b"ar" => set(Setter::Genotype2, b"ra".to_vec()),
                _ if key == b"aA" || key == b"Aa" => set(Setter::Genotype2, b"aA".to_vec()),
                b"aa" => set(Setter::Genotype2, b"aa".to_vec()),
                b"a" => set(Setter::Genotype2, b"a".to_vec()),
                b"r" => set(Setter::Genotype2, b"r".to_vec()),
                _ => (None, key),
            };
            if let Some(s) = setter {
                out[iu].setter = Some(s);
                out[ival as usize].key = Some(newkey);
            }
            i += 1;
            continue;
        }
        if out[iu].hl_type == HL_FLT {
            if i + 1 == n {
                return Err(parse_err());
            }
            let itok = iu;
            let t1 = peek(out, i + 1)?.tok_type;
            let ival = if t1 == TOK_EQ || t1 == TOK_NE {
                i - 1
            } else if t1 == TOK_LIKE {
                out[iu + 1].tok_type = TOK_IN;
                i - 1
            } else if t1 == TOK_NLIKE {
                out[iu + 1].tok_type = TOK_NOT_IN;
                i - 1
            } else {
                let t2 = peek(out, i + 2)?.tok_type;
                let inside = (iu + 2) < out.len();
                if t2 == TOK_EQ || t2 == TOK_NE {
                    i += 1;
                } else if t2 == TOK_LIKE && inside {
                    out[iu + 2].tok_type = TOK_IN;
                    i += 1;
                } else if t2 == TOK_NLIKE && inside {
                    out[iu + 2].tok_type = TOK_NOT_IN;
                    i += 1;
                } else {
                    return Err(deep_err());
                }
                i
            };
            let tv = peek(out, ival)?;
            if tv.tok_type != TOK_VAL || tv.key.is_none() {
                return Err(init_error(format!(
                    "[filter.c] Could not parse the expression, an unquoted string value perhaps? {shown}"
                )));
            }
            let key = tv.key.clone().unwrap();
            let mut ids = Vec::new();
            let mut has_missing = false;
            let mut bp = 0;
            while bp < key.len() {
                let ep = key[bp..]
                    .iter()
                    .position(|&c| c == b';')
                    .map_or(key.len(), |p| bp + p);
                let name = &key[bp..ep];
                if name == b"." {
                    has_missing = true;
                } else {
                    let id = hdr.id(name);
                    if !hdr.exists(HL_FLT, id) {
                        return Err(init_error(format!(
                            "The filter \"{}\" not present in the VCF header",
                            String::from_utf8_lossy(name)
                        )));
                    }
                    ids.push(id);
                }
                if ep == key.len() {
                    break;
                }
                bp = ep + 1;
            }
            if has_missing && !ids.is_empty() {
                return Err(init_error(format!(
                    "The FILTER expression cannot contain missing value AND filters: \"{}\" ({})",
                    String::from_utf8_lossy(&key),
                    ids.len()
                )));
            }
            let tv = &mut out[ival as usize];
            tv.nvalues = ids.len() as i32;
            tv.values.expand(ids.len());
            for (k, &id) in ids.iter().enumerate() {
                tv.values.set(k as isize, f64::from(id));
            }
            tv.tag = tv.key.take();
            let id = tv.hdr_id;
            out[itok].hdr_id = id;
            i += 1;
            continue;
        }
        i += 1;
    }
    Ok(())
}

#[cfg(test)]
mod shadow_tests {
    use super::{KStr, take_undefined};

    #[test]
    fn nominal_capacity_keeps_only_written_bytes_and_backfills_sparse_points() {
        take_undefined();
        let mut s = KStr::default();
        s.resize(usize::MAX / 8);
        assert!(s.m >= usize::MAX / 8);
        assert!(s.buf.is_empty());
        assert!(s.sparse.is_none());
        s.set(7, b'Z');
        s.set(8, 0);
        assert_eq!(s.get(7), b'Z');
        assert_eq!(s.get(8), 0);
        assert!(take_undefined().is_none());
        assert_eq!(s.get(5), 0);
        assert!(take_undefined().is_some());
        for (i, &c) in b"abcdefg".iter().enumerate() {
            s.set(i, c);
        }
        assert_eq!(s.buf, b"abcdefgZ\0");
        assert_eq!(s.cstr(0), b"abcdefgZ");
        assert!(take_undefined().is_none());
        s.set(20, b'Q');
        s.reallocate(4);
        assert_eq!(s.buf, b"abcd");
        assert!(s.sparse.as_ref().is_none_or(|far| far.is_empty()));
        assert_eq!(s.get(20), 0);
        assert!(take_undefined().is_some());
        s.grow_exact(8);
        assert_eq!(s.get(4), 0);
        assert!(take_undefined().is_some());
    }

    #[test]
    fn overflowing_string_metadata_returns_undefined_without_writing() {
        take_undefined();
        let mut s = KStr {
            l: usize::MAX,
            ..KStr::default()
        };
        s.puts(b"x");
        assert!(take_undefined().is_some());
        s.putc(b'x');
        assert!(take_undefined().is_some());
        s.putc_(b'x');
        assert!(take_undefined().is_some());
        s.set(usize::MAX, b'x');
        assert!(take_undefined().is_some());
        assert!(s.buf.is_empty());
        assert!(s.sparse.is_none());
    }

    #[test]
    fn bulk_string_writes_preserve_suffixes_and_sparse_definedness() {
        for prefix in [0, 1, 7, 64] {
            for start in [0, 1, 7, 64, 80] {
                for len in [0, 1, 8, 65] {
                    for sparse in [false, true] {
                        let mut a = KStr::default();
                        a.resize(200);
                        for i in 0..prefix {
                            a.set(i, i as u8);
                        }
                        if sparse {
                            a.set(70, 0xfe);
                            a.set(71, 0);
                        }
                        a.l = start;
                        let mut b = a.clone();
                        let bytes: Vec<_> = (0..len).map(|i| (i * 37) as u8).collect();
                        a.puts(&bytes);
                        let end = b.l + bytes.len();
                        b.resize(end + 2);
                        for (k, &c) in bytes.iter().enumerate() {
                            b.set(b.l + k, c);
                        }
                        b.l = end;
                        b.set(end, 0);
                        assert_eq!(a.buf, b.buf);
                        assert_eq!(a.l, b.l);
                        assert_eq!(a.m, b.m);
                        assert_eq!(a.sparse, b.sparse);
                    }
                }
            }
        }
    }
}
