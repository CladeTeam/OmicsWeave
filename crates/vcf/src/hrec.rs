//! VCF header text parsed into htslib's header records (`bcf_hrec_t`) as
//! `bcf_hdr_parse` builds them, and formatted back as `bcf_hdr_format`
//! writes them into a BCF header.
//!
//! htslib does not keep the header text: it keeps one record per accepted
//! `##` line, normalised on the way in (keys and values trimmed, quoted
//! values re-quoted, `[...]` values kept whole, trailing junk dropped),
//! drops duplicates and lines it cannot register, inserts the `PASS`
//! filter after the first line and appends `IDX` to every dictionary
//! declaration. The BCF header is those records printed in order, so the
//! dictionary ids and the header bytes come from the same model here.

use std::collections::{BTreeMap, HashMap, HashSet};

use omicsweave_core::{Error, Result};

use crate::cnum::is_space;
use crate::header::{Contig, Header, Key, Number, ValueType, version};

fn invalid(detail: impl Into<String>) -> Error {
    Error::invalid("vcf header", detail)
}

/// `BCF_HL_*`: what a header record declares.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum Kind {
    Filter = 0,
    Info = 1,
    Format = 2,
    Contig = 3,
    Structured = 4,
    Generic = 5,
}

/// One header record (`bcf_hrec_t`). A generic line has `value`; a
/// structured one has `keys`/`vals` (quoted values keep their quotes).
#[derive(Debug, Clone, PartialEq, Eq)]
struct Hrec {
    key: Vec<u8>,
    value: Option<Vec<u8>>,
    keys: Vec<Vec<u8>>,
    vals: Vec<Vec<u8>>,
}

impl Hrec {
    /// `bcf_hrec_set_type`.
    fn kind(&self) -> Kind {
        match &self.key[..] {
            b"contig" => Kind::Contig,
            b"INFO" => Kind::Info,
            b"FILTER" => Kind::Filter,
            b"FORMAT" => Kind::Format,
            _ if !self.keys.is_empty() => Kind::Structured,
            _ => Kind::Generic,
        }
    }

    /// `bcf_hrec_find_key`: the first key equal ignoring ASCII case.
    fn find_key(&self, name: &[u8]) -> Option<usize> {
        self.keys.iter().position(|k| k.eq_ignore_ascii_case(name))
    }

    /// `hrec_add_idx`.
    fn add_idx(&mut self, id: u32) {
        self.keys.push(b"IDX".to_vec());
        self.vals.push(id.to_string().into_bytes());
    }

    /// `_bcf_hrec_format` with `is_bcf` set (IDX kept).
    fn format(&self, out: &mut Vec<u8>) {
        out.extend_from_slice(b"##");
        out.extend_from_slice(&self.key);
        out.push(b'=');
        match &self.value {
            Some(v) => out.extend_from_slice(v),
            None => {
                out.push(b'<');
                for (j, (k, v)) in self.keys.iter().zip(&self.vals).enumerate() {
                    if j > 0 {
                        out.push(b',');
                    }
                    out.extend_from_slice(k);
                    out.push(b'=');
                    out.extend_from_slice(v);
                }
                out.push(b'>');
            }
        }
        out.push(b'\n');
    }
}

enum Line {
    Record(Hrec),
    /// A `##` line htslib cannot parse; it is skipped.
    Malformed,
    /// Not a `##` line.
    NotMeta,
}

/// `is_escaped`: an odd number of backslashes before `at`, not looking
/// before `min`.
fn is_escaped(line: &[u8], min: usize, at: usize) -> bool {
    let mut n = 0;
    let mut s = at;
    while s > min {
        s -= 1;
        if line[s] != b'\\' {
            break;
        }
        n += 1;
    }
    n % 2 == 1
}

/// `bcf_hdr_parse_line` on one line (without its newline).
fn parse_line(line: &[u8]) -> Line {
    let c = |i: usize| line.get(i).copied().unwrap_or(0);
    if c(0) != b'#' || c(1) != b'#' {
        return Line::NotMeta;
    }
    let mut p = 2;
    let mut q = p;
    while c(q) != 0 && c(q) != b'=' {
        q += 1;
    }
    if c(q) != b'=' || q == p {
        return Line::Malformed;
    }
    let mut hrec = Hrec {
        key: line[p..q].to_vec(),
        value: None,
        keys: Vec::new(),
        vals: Vec::new(),
    };
    q += 1;
    p = q;
    if c(p) != b'<' {
        // Generic: everything to the end of the line.
        hrec.value = Some(line[p..].to_vec());
        return Line::Record(hrec);
    }

    let mut nopen = 1u32;
    while c(q) != 0 && nopen > 0 {
        q += 1;
        p = q;
        while c(q) == b' ' {
            p += 1;
            q += 1;
        }
        // ^[A-Za-z_][0-9A-Za-z_.]*$
        if c(q).is_ascii_alphabetic() || c(q) == b'_' {
            q += 1;
            while c(q).is_ascii_alphanumeric() || c(q) == b'_' || c(q) == b'.' {
                q += 1;
            }
        }
        let n = q - p;
        while c(q) == b' ' {
            q += 1;
        }
        if c(q) != b'=' || n == 0 {
            return Line::Malformed;
        }
        hrec.keys.push(line[p..p + n].to_vec());
        q += 1;
        p = q;
        while c(q) == b' ' {
            p += 1;
            q += 1;
        }
        let (mut quoted, ending) = match c(p) {
            b'"' => {
                p += 1;
                (true, b'"')
            }
            b'[' => (true, b']'),
            _ => (false, 0),
        };
        if quoted {
            q += 1;
        }
        while c(q) != 0 {
            if quoted {
                if c(q) == ending && !is_escaped(line, p, q) {
                    break;
                }
            } else {
                if c(q) == b'<' {
                    nopen += 1;
                }
                if c(q) == b'>' {
                    nopen -= 1;
                }
                if nopen == 0 || (c(q) == b',' && nopen == 1) {
                    break;
                }
            }
            q += 1;
        }
        let mut r = q;
        if quoted && ending == b']' {
            if c(q) != ending {
                return Line::Malformed;
            }
            r += 1;
            q += 1;
            quoted = false;
        }
        while r > p && line[r - 1] == b' ' {
            r -= 1;
        }
        let mut val = Vec::with_capacity(r - p + 2);
        if quoted {
            val.push(b'"');
        }
        val.extend_from_slice(&line[p..r]);
        if quoted {
            val.push(b'"');
        }
        hrec.vals.push(val);
        if quoted && c(q) == ending {
            q += 1;
        }
        if c(q) == b'>' {
            nopen = nopen.saturating_sub(1);
            q += 1;
        }
    }
    if nopen != 0 {
        return Line::Malformed;
    }
    // Anything after the closing '>' is dropped.
    Line::Record(hrec)
}

/// C `strtol`/`strtoll` in base 10: the value (clamped to the `long`
/// range) and the bytes read; 0 bytes when there are no digits.
fn strtol(s: &[u8]) -> (i64, usize) {
    let mut i = s.iter().take_while(|&&c| is_space(c)).count();
    let neg = match s.get(i) {
        Some(b'-') => {
            i += 1;
            true
        }
        Some(b'+') => {
            i += 1;
            false
        }
        _ => false,
    };
    let digits = s[i..].iter().take_while(|c| c.is_ascii_digit()).count();
    if digits == 0 {
        return (0, 0);
    }
    let magnitude = s[i..i + digits]
        .iter()
        .fold(0i128, |v, &d| (v * 10 + i128::from(d - b'0')).min(1 << 64));
    let v = if neg { -magnitude } else { magnitude };
    (v.clamp(i64::MIN.into(), i64::MAX.into()) as i64, i + digits)
}

/// `IDX` as htslib reads it, `strtol` into an `int`: `None` unless the
/// whole value is read and lies in `[0, INT_MAX - 1)`.
fn parse_idx(s: &[u8]) -> Option<u32> {
    let (v, used) = strtol(s);
    // With no digits, the end is the start: only an empty value passes.
    if used != s.len() {
        return None;
    }
    let v = v as i32;
    (0..i32::MAX - 1).contains(&v).then_some(v as u32)
}

/// `bcf_idinfo_t` of a FILTER/INFO/FORMAT key: its id and the packed
/// declaration word per line kind (FILTER, INFO, FORMAT).
#[derive(Debug, Clone, Copy)]
struct IdEntry {
    id: u32,
    words: [Option<u32>; 3],
}

/// One of htslib's id dictionaries: ids in use and the next free one
/// (`hdr->n[dict_type]`).
#[derive(Debug, Default)]
struct Slots<T> {
    by_id: BTreeMap<u32, T>,
    n: u64,
}

impl<T> Slots<T> {
    /// `bcf_hdr_set_idx`: the given IDX, or the next id.
    fn assign(&mut self, idx: Option<u32>, value: T) -> Result<u32> {
        let id = match idx {
            None => u32::try_from(self.n).map_err(|_| invalid("too many dictionary ids"))?,
            Some(i) => {
                if u64::from(i) < self.n && self.by_id.contains_key(&i) {
                    return Err(invalid(format!("conflicting IDX={i} lines")));
                }
                i
            }
        };
        self.n = self.n.max(u64::from(id) + 1);
        self.by_id.insert(id, value);
        Ok(id)
    }
}

/// `bcf_hdr_t` while it is being built from the text.
#[derive(Debug, Default)]
struct Builder {
    hrecs: Vec<Hrec>,
    ids: HashMap<Vec<u8>, IdEntry>,
    id_slots: Slots<Vec<u8>>,
    contigs: HashMap<Vec<u8>, u32>,
    contig_slots: Slots<(Vec<u8>, Option<u64>)>,
    /// `aux->gen`: generic lines and structured lines with an ID, as text.
    seen: HashSet<Vec<u8>>,
    /// The value of the last `##fileformat` line added.
    fileformat: Option<Vec<u8>>,
    samples: Vec<String>,
}

impl Builder {
    /// `bcf_hdr_add_hrec`.
    fn add(&mut self, mut hrec: Hrec) -> Result<()> {
        let kind = hrec.kind();
        let dedup = match kind {
            Kind::Contig => {
                if !self.register_contig(&mut hrec)? {
                    return Ok(());
                }
                None
            }
            Kind::Filter | Kind::Info | Kind::Format => {
                if !self.register_key(&mut hrec, kind)? {
                    return Ok(());
                }
                None
            }
            Kind::Generic => {
                let mut s = b"##".to_vec();
                s.extend_from_slice(&hrec.key);
                s.push(b'=');
                s.extend_from_slice(hrec.value.as_deref().unwrap_or_default());
                if self.seen.contains(&s) {
                    return Ok(());
                }
                if hrec.key == b"fileformat" {
                    self.fileformat = hrec.value.clone();
                }
                Some(s)
            }
            Kind::Structured => match hrec.find_key(b"ID") {
                Some(i) => {
                    let mut s = b"##".to_vec();
                    s.extend_from_slice(&hrec.key);
                    s.extend_from_slice(b"=<ID=");
                    s.extend_from_slice(&hrec.vals[i]);
                    s.push(b'>');
                    if self.seen.contains(&s) {
                        return Ok(());
                    }
                    Some(s)
                }
                None => None,
            },
        };
        self.hrecs.push(hrec);
        if let Some(s) = dedup {
            self.seen.insert(s);
        }
        Ok(())
    }

    /// `bcf_hdr_register_hrec_ctg`: whether the record is kept.
    fn register_contig(&mut self, hrec: &mut Hrec) -> Result<bool> {
        let Some(id_i) = hrec.find_key(b"ID") else {
            return Ok(false);
        };
        let len = hrec.find_key(b"length").and_then(|i| {
            let v = &hrec.vals[i];
            let v = v.strip_prefix(b"\"").unwrap_or(v);
            match strtol(v) {
                (_, 0) => None,
                (l, _) => u64::try_from(l).ok(),
            }
        });
        let name = hrec.vals[id_i].clone();
        if self.contigs.contains_key(&name) {
            return Ok(false);
        }
        let idx = match hrec.find_key(b"IDX") {
            Some(i) => match parse_idx(&hrec.vals[i]) {
                Some(v) => Some(v),
                None => return Ok(false),
            },
            None => None,
        };
        let id = self.contig_slots.assign(idx, (name.clone(), len))?;
        self.contigs.insert(name, id);
        if idx.is_none() {
            hrec.add_idx(id);
        }
        Ok(true)
    }

    /// `bcf_hdr_register_hrec_info` for FILTER, INFO and FORMAT lines:
    /// whether the record is kept.
    fn register_key(&mut self, hrec: &mut Hrec, kind: Kind) -> Result<bool> {
        let mut id = None;
        let mut idx = None;
        let (mut ty, mut var, mut num) = (u32::MAX, u32::MAX, -1i32);
        for (k, v) in hrec.keys.iter().zip(&hrec.vals) {
            match &k[..] {
                b"ID" => id = Some(v.clone()),
                b"IDX" => match parse_idx(v) {
                    Some(i) => idx = Some(i),
                    None => return Ok(false),
                },
                b"Type" => {
                    ty = match &v[..] {
                        b"Integer" => 1,
                        b"Float" => 2,
                        b"Flag" => 0,
                        _ => 3,
                    }
                }
                b"Number" => {
                    let format = kind == Kind::Format;
                    match &v[..] {
                        b"A" => var = 2,
                        b"R" => var = 4,
                        b"G" => var = 3,
                        b"." => var = 1,
                        b"P" if format => var = 5,
                        b"LA" if format => var = 6,
                        b"LR" if format => var = 8,
                        b"LG" if format => var = 7,
                        b"M" if format => var = 9,
                        _ => {
                            if let Some(n) = scan_int(v) {
                                num = n;
                                var = 0;
                            }
                        }
                    }
                    if var != 0 {
                        num = 0xfffff;
                    }
                }
                _ => {}
            }
        }
        if matches!(kind, Kind::Info | Kind::Format) {
            if ty == u32::MAX {
                ty = 3;
            }
            if var == u32::MAX {
                var = 1;
            }
            if ty == 0 && (var != 0 || num != 0) {
                var = 0;
                num = 0;
            }
        }
        let word = (num as u32 & 0xfffff) << 12 | (var & 0xf) << 8 | (ty & 0xf) << 4 | kind as u32;
        let Some(name) = id else {
            return Ok(false);
        };
        let slot = kind as usize;
        if let Some(entry) = self.ids.get_mut(&name) {
            if entry.words[slot].is_some() {
                // A user PASS replaces the one htslib inserted after the
                // first line (it checks the second record is a FILTER).
                if name == b"PASS"
                    && kind == Kind::Filter
                    && self.hrecs.len() > 1
                    && self.hrecs[1].kind() == Kind::Filter
                {
                    if idx.is_none() {
                        hrec.add_idx(entry.id);
                    }
                    std::mem::swap(&mut self.hrecs[1], hrec);
                }
                return Ok(false);
            }
            entry.words[slot] = Some(word);
            if idx.is_none() {
                hrec.add_idx(entry.id);
            }
            return Ok(true);
        }
        let id = self.id_slots.assign(idx, name.clone())?;
        let mut words = [None; 3];
        words[slot] = Some(word);
        self.ids.insert(name, IdEntry { id, words });
        if idx.is_none() {
            hrec.add_idx(id);
        }
        Ok(true)
    }

    /// `bcf_hdr_parse_sample_line` on the text from the `#CHROM` line on.
    fn parse_samples(&mut self, rest: &[u8]) -> Result<()> {
        const MANDATORY: &[u8] = b"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO";
        let c = |i: usize| rest.get(i).copied().unwrap_or(0);
        if !rest.starts_with(MANDATORY) {
            return Err(invalid("invalid #CHROM columns"));
        }
        let mut beg = MANDATORY.len();
        if c(beg) == 0 || c(beg) == b'\n' {
            return Ok(());
        }
        if !rest[beg..].starts_with(b"\tFORMAT\t") {
            return Err(invalid("invalid FORMAT column"));
        }
        beg += 8;
        let mut names = HashSet::new();
        while c(beg) != 0 {
            let mut end = beg;
            while c(end) != 0 && c(end) != b'\t' && c(end) != b'\n' {
                end += 1;
            }
            let name = &rest[beg..end];
            if name.iter().all(|&b| is_space(b)) {
                return Err(invalid("empty sample name"));
            }
            let name = std::str::from_utf8(name)
                .map_err(|_| invalid("sample name is not UTF-8"))?
                .to_owned();
            if !names.insert(name.clone()) {
                return Err(invalid(format!("duplicated sample name '{name}'")));
            }
            self.samples.push(name);
            if c(end) == 0 || c(end) == b'\n' {
                break;
            }
            beg = end + 1;
        }
        Ok(())
    }

    /// The BCF header text (`bcf_hdr_format` with `is_bcf`), without the
    /// terminating NUL.
    fn format(&self) -> Vec<u8> {
        let mut out = Vec::new();
        for h in &self.hrecs {
            h.format(&mut out);
        }
        out.extend_from_slice(b"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO");
        if !self.samples.is_empty() {
            out.extend_from_slice(b"\tFORMAT");
            for s in &self.samples {
                out.push(b'\t');
                out.extend_from_slice(s.as_bytes());
            }
        }
        out.push(b'\n');
        out
    }
}

/// `sscanf(s, "%d")`: leading space, a sign and digits, the rest ignored.
fn scan_int(s: &[u8]) -> Option<i32> {
    let (v, used) = strtol(s);
    (used > 0).then_some(v as i32)
}

fn value_type(word: u32) -> ValueType {
    match word >> 4 & 0xf {
        0 => ValueType::Flag,
        1 => ValueType::Integer,
        2 => ValueType::Float,
        _ => ValueType::String,
    }
}

fn number(word: u32) -> Number {
    match word >> 8 & 0xf {
        0 => Number::Fixed(word >> 12),
        2 => Number::A,
        3 => Number::G,
        4 => Number::R,
        5 => Number::P,
        6 => Number::LA,
        7 => Number::LG,
        8 => Number::LR,
        9 => Number::M,
        _ => Number::Var,
    }
}

fn name(bytes: Vec<u8>, what: &str) -> Result<String> {
    String::from_utf8(bytes).map_err(|_| invalid(format!("{what} name is not UTF-8")))
}

/// Parses header text as `bcf_hdr_parse` does. The text ends at its first
/// NUL, as a C string does.
pub(crate) fn parse(text: &[u8]) -> Result<Header> {
    let text = match text.iter().position(|&b| b == 0) {
        Some(end) => &text[..end],
        None => text,
    };
    let line_at = |p: usize| -> (&[u8], usize) {
        match text[p..].iter().position(|&b| b == b'\n') {
            Some(n) => (&text[p..p + n], p + n + 1),
            None => (&text[p..], text.len()),
        }
    };
    let mut b = Builder::default();
    // The first line is added before htslib's PASS filter, then read again
    // with the others (a second copy is dropped unless nothing dedupes it).
    if let Line::Record(h) = parse_line(line_at(0).0) {
        b.add(h)?;
    }
    if let Line::Record(h) = parse_line(b"##FILTER=<ID=PASS,Description=\"All filters passed\">") {
        b.add(h)?;
    }
    let mut p = 0;
    loop {
        let (line, next) = line_at(p);
        match parse_line(line) {
            Line::Record(h) => {
                b.add(h)?;
                p = next;
                continue;
            }
            Line::Malformed => {
                p = next;
                continue;
            }
            Line::NotMeta => {}
        }
        if line.starts_with(b"#CHROM\t") || line.starts_with(b"#CHROM ") {
            break;
        }
        if p + line.len() >= text.len() {
            return Err(invalid("no #CHROM line"));
        }
        p = next;
    }
    b.parse_samples(&text[p..])?;

    let bcf_text = b.format();
    let version = b
        .fileformat
        .as_deref()
        .map_or(crate::header::DEFAULT_VERSION, |v| {
            version(&String::from_utf8_lossy(v))
        });
    let mut contig_ids = HashMap::with_capacity(b.contig_slots.by_id.len());
    let mut key_ids = HashMap::with_capacity(b.id_slots.by_id.len());
    let mut contigs = Vec::new();
    let n = usize::try_from(b.contig_slots.n).map_err(|_| invalid("too many contigs"))?;
    contigs
        .try_reserve(n)
        .map_err(|_| invalid("contig allocation failed"))?;
    contigs.resize(
        n,
        Contig {
            name: String::new(),
            len: None,
        },
    );
    for (id, (n, len)) in std::mem::take(&mut b.contig_slots.by_id) {
        let name = name(n, "contig")?;
        contig_ids.insert(name.clone(), id);
        contigs[id as usize] = Contig { name, len };
    }
    let mut keys: Vec<Option<Key>> = Vec::new();
    let n = usize::try_from(b.id_slots.n).map_err(|_| invalid("too many dictionary ids"))?;
    keys.try_reserve(n)
        .map_err(|_| invalid("dictionary allocation failed"))?;
    keys.resize(n, None);
    for (id, n) in std::mem::take(&mut b.id_slots.by_id) {
        let words = b.ids[&n].words;
        let name = name(n, "dictionary")?;
        key_ids.insert(name.clone(), id);
        let mut key = Key {
            name,
            filter: words[0].is_some(),
            ..Key::default()
        };
        if let Some(w) = words[1] {
            key.info = Some(value_type(w));
            key.info_number = number(w);
        }
        if let Some(w) = words[2] {
            key.format = Some(value_type(w));
            key.format_number = number(w);
        }
        keys[id as usize] = Some(key);
    }
    Ok(Header::from_parts(
        version,
        contigs,
        std::mem::take(&mut b.samples),
        keys,
        contig_ids,
        key_ids,
        bcf_text,
    ))
}

#[cfg(test)]
mod tests {
    use super::*;

    fn record(line: &str) -> Hrec {
        match parse_line(line.as_bytes()) {
            Line::Record(h) => h,
            _ => panic!("not parsed: {line}"),
        }
    }

    #[test]
    fn values_are_trimmed_and_requoted() {
        let h = record(r#"##INFO=< ID = DP ,Number=1, Description=" a \"b\" ",X=[1, 2]>junk"#);
        assert_eq!(h.keys, [&b"ID"[..], b"Number", b"Description", b"X"]);
        assert_eq!(h.vals, [&b"DP"[..], b"1", br#"" a \"b\"""#, b"[1, 2]"]);
        // After a closing quote htslib skips one byte before the next key.
        assert!(matches!(
            parse_line(br#"##INFO=<ID=DP,Description="x"  ,Number=1>"#),
            Line::Malformed
        ));
    }

    #[test]
    fn idx_reads_as_strtol_into_int() {
        assert_eq!(parse_idx(b""), Some(0));
        assert_eq!(parse_idx(b" +7"), Some(7));
        assert_eq!(parse_idx(b"4294967296"), Some(0));
        assert_eq!(parse_idx(b"4294967295"), None);
        assert_eq!(parse_idx(b"7x"), None);
        assert_eq!(parse_idx(b"\"7\""), None);
    }
}
