//! The VCF header: version, contigs, samples and the dictionary of FILTER,
//! INFO and FORMAT keys with their declared types. BCF carries the same text.

use std::collections::HashMap;
use std::sync::Arc;

use omicsweave_core::{Error, Result};

/// Version used when `##fileformat` is absent or unreadable (htslib: 4.2).
pub(crate) const DEFAULT_VERSION: i32 = 4_002_000;
/// VCF 4.4, from which a genotype may carry a leading phase indicator.
pub(crate) const VCF44: i32 = 4_004_000;

/// The value type a key is declared with (`Character` reads as `String`).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ValueType {
    Flag,
    Integer,
    Float,
    String,
}

impl ValueType {
    fn parse(s: &str) -> ValueType {
        match s {
            "Integer" => ValueType::Integer,
            "Float" => ValueType::Float,
            "Flag" => ValueType::Flag,
            _ => ValueType::String,
        }
    }

    /// htslib's `BCF_HT_*` code.
    pub fn code(self) -> u32 {
        match self {
            ValueType::Flag => 0,
            ValueType::Integer => 1,
            ValueType::Float => 2,
            ValueType::String => 3,
        }
    }
}

/// The declared count of a key's values (`Number=`) as htslib keeps it: a
/// fixed count (20 bits wide) or one of the variable kinds. `P`, `LA`,
/// `LG`, `LR` and `M` are recognised for FORMAT keys only.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum Number {
    Fixed(u32),
    /// `.`, also assumed when `Number` is absent or unreadable.
    #[default]
    Var,
    A,
    G,
    R,
    P,
    LA,
    LG,
    LR,
    M,
}

impl Number {
    /// `sscanf("%d")`-style reading of `Number=`; `None` when it fails.
    fn parse(s: &str, format: bool) -> Option<Number> {
        Some(match s {
            "A" => Number::A,
            "R" => Number::R,
            "G" => Number::G,
            "." => Number::Var,
            "P" if format => Number::P,
            "LA" if format => Number::LA,
            "LR" if format => Number::LR,
            "LG" if format => Number::LG,
            "M" if format => Number::M,
            _ => {
                let t = s.trim_start_matches([' ', '\t', '\n', '\r', '\x0b', '\x0c']);
                let (negative, digits) = match t.as_bytes().first() {
                    Some(b'-') => (true, &t[1..]),
                    Some(b'+') => (false, &t[1..]),
                    _ => (false, t),
                };
                let n = digits.bytes().take_while(u8::is_ascii_digit).count();
                if n == 0 {
                    return None;
                }
                let v = digits[..n].bytes().fold(0i64, |v, d| {
                    v.saturating_mul(10).saturating_add(i64::from(d - b'0'))
                });
                let v = if negative { v.saturating_neg() } else { v };
                Number::Fixed(v as i32 as u32 & 0xfffff)
            }
        })
    }

    /// htslib's `BCF_VL_*` code (`bcf_hdr_id2length`).
    pub fn length_code(self) -> u32 {
        match self {
            Number::Fixed(_) => 0,
            Number::Var => 1,
            Number::A => 2,
            Number::G => 3,
            Number::R => 4,
            Number::P => 5,
            Number::LA => 6,
            Number::LG => 7,
            Number::LR => 8,
            Number::M => 9,
        }
    }

    /// The count `bcf_hdr_id2number` reports: 0xfffff for variable kinds.
    pub fn count(self) -> u32 {
        match self {
            Number::Fixed(n) => n & 0xfffff,
            _ => 0xfffff,
        }
    }
}

/// The header line kinds a dictionary key can be declared under
/// (`BCF_HL_FLT`, `BCF_HL_INFO`, `BCF_HL_FMT`).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum LineKind {
    Filter = 0,
    Info = 1,
    Format = 2,
}

/// A key of the shared FILTER/INFO/FORMAT dictionary.
#[derive(Debug, Clone, Default, PartialEq, Eq)]
pub struct Key {
    pub name: String,
    pub info: Option<ValueType>,
    pub format: Option<ValueType>,
    pub filter: bool,
    /// `Number=` of the INFO declaration (meaningful when `info` is set).
    pub info_number: Number,
    /// `Number=` of the FORMAT declaration (meaningful when `format` is set).
    pub format_number: Number,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Contig {
    pub name: String,
    pub len: Option<u64>,
}

/// A FILTER/INFO/FORMAT declaration: (ID, kind, type and number, IDX).
type Declaration = (String, Kind, Option<(ValueType, Number)>, Option<u32>);

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum Kind {
    Info,
    Format,
    Filter,
}

#[derive(Debug, Clone, Default)]
pub struct Header {
    /// `major * 1_000_000 + minor * 1_000`, as htslib numbers versions.
    pub version: i32,
    /// Indexed by contig id ([`crate::Record::rid`]).
    pub contigs: Vec<Contig>,
    pub samples: Vec<String>,
    /// Indexed by dictionary id; gaps where IDX values skip ids.
    pub keys: Vec<Option<Key>>,
    contig_ids: HashMap<String, u32>,
    key_ids: HashMap<String, u32>,
    /// The header text htslib writes into a BCF header, when the header
    /// was parsed by [`Header::from_vcf_text`].
    pub(crate) bcf_text: Option<Arc<[u8]>>,
}

/// Splits `ID=x,Number=1,Description="a, b"` into key-value pairs.
fn structured(body: &str) -> Vec<(String, String)> {
    let mut out = Vec::new();
    let b = body.as_bytes();
    let mut i = 0;
    while i < b.len() {
        let key_start = i;
        while i < b.len() && b[i] != b'=' && b[i] != b',' {
            i += 1;
        }
        let key = body[key_start..i].trim().to_owned();
        let mut value = String::new();
        if i < b.len() && b[i] == b'=' {
            i += 1;
            if i < b.len() && b[i] == b'"' {
                i += 1;
                let start = i;
                while i < b.len() && b[i] != b'"' {
                    i += if b[i] == b'\\' { 2 } else { 1 };
                }
                value = body[start..i.min(b.len())].to_owned();
                while i < b.len() && b[i] != b',' {
                    i += 1;
                }
            } else {
                let start = i;
                while i < b.len() && b[i] != b',' {
                    i += 1;
                }
                value = body[start..i].to_owned();
            }
        }
        out.push((key, value));
        i += 1;
    }
    out
}

pub(crate) fn version(value: &str) -> i32 {
    let parse = || -> Option<i32> {
        let rest = &value[value.find("VCFv")? + 4..];
        let (major, minor) = rest.split_once('.')?;
        let digits = |s: &str| {
            s.chars()
                .take_while(char::is_ascii_digit)
                .collect::<String>()
        };
        let major = digits(major).parse::<i32>().ok()?;
        let minor = digits(minor).parse::<i32>().ok()?;
        if major >= 1000 || minor >= 1000 {
            return None;
        }
        Some(major * 1_000_000 + minor * 1_000)
    };
    parse().unwrap_or(DEFAULT_VERSION)
}

impl Header {
    /// Parses VCF header text (the `##` lines and the `#CHROM` line) as
    /// htslib's `bcf_hdr_parse` does, keeping the header records htslib
    /// would write into a BCF header (see [`crate::bcf_header_bytes`]).
    /// Dictionary ids are assigned as htslib assigns them: `PASS` after
    /// the first line, then declarations in order unless they carry `IDX`.
    /// The text ends at its first NUL. Contig, key and sample names must
    /// be UTF-8; other header bytes are kept as they are.
    pub fn from_vcf_text(text: &[u8]) -> Result<Header> {
        crate::hrec::parse(text)
    }

    pub(crate) fn from_parts(
        version: i32,
        contigs: Vec<Contig>,
        samples: Vec<String>,
        keys: Vec<Option<Key>>,
        contig_ids: HashMap<String, u32>,
        key_ids: HashMap<String, u32>,
        bcf_text: Vec<u8>,
    ) -> Header {
        Header {
            version,
            contigs,
            samples,
            keys,
            contig_ids,
            key_ids,
            bcf_text: Some(bcf_text.into()),
        }
    }

    /// Parses header text: the `##` lines and the `#CHROM` line.
    pub fn parse(text: &str) -> Result<Header> {
        let mut h = Header {
            version: DEFAULT_VERSION,
            ..Header::default()
        };
        let mut declared: Vec<Declaration> = Vec::new();
        let mut contigs: Vec<(Contig, Option<u32>)> = Vec::new();
        let mut columns = None;
        for line in text.split('\n') {
            let line = line.trim_end_matches(['\r', '\0']);
            if let Some(meta) = line.strip_prefix("##") {
                let Some((key, value)) = meta.split_once('=') else {
                    continue;
                };
                if key == "fileformat" {
                    h.version = version(value);
                    continue;
                }
                let Some(body) = value.strip_prefix('<').and_then(|v| v.strip_suffix('>')) else {
                    continue;
                };
                let fields = structured(body);
                let get = |k: &str| fields.iter().find(|(f, _)| f == k).map(|(_, v)| v.as_str());
                let Some(id) = get("ID") else { continue };
                let idx = get("IDX")
                    .map(|v| {
                        v.parse::<u32>()
                            .ok()
                            .filter(|&i| i <= i32::MAX as u32)
                            .ok_or_else(|| Error::invalid("vcf header", "invalid IDX"))
                    })
                    .transpose()?;
                let kind = match key {
                    "contig" => {
                        let len = get("length").and_then(|v| v.parse().ok());
                        contigs.push((
                            Contig {
                                name: id.to_owned(),
                                len,
                            },
                            idx,
                        ));
                        continue;
                    }
                    "INFO" => Kind::Info,
                    "FORMAT" => Kind::Format,
                    "FILTER" => Kind::Filter,
                    _ => continue,
                };
                let ty = (kind != Kind::Filter).then(|| {
                    let ty = ValueType::parse(get("Type").unwrap_or(""));
                    // htslib keeps the last Number; a Flag is always Number=0.
                    let number = fields
                        .iter()
                        .rev()
                        .find(|(f, _)| f == "Number")
                        .and_then(|(_, v)| Number::parse(v, kind == Kind::Format))
                        .unwrap_or(Number::Var);
                    let number = if ty == ValueType::Flag {
                        Number::Fixed(0)
                    } else {
                        number
                    };
                    (ty, number)
                });
                declared.push((id.to_owned(), kind, ty, idx));
            } else if let Some(cols) = line.strip_prefix('#') {
                columns = Some(cols.to_owned());
            }
        }
        let cols = columns.ok_or_else(|| Error::invalid("vcf header", "no #CHROM line"))?;
        let mut fields = cols.split('\t');
        for expected in ["CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO"] {
            if fields.next() != Some(expected) {
                return Err(Error::invalid("vcf header", "invalid #CHROM columns"));
            }
        }
        if let Some(format) = fields.next() {
            if format != "FORMAT" {
                return Err(Error::invalid("vcf header", "invalid FORMAT column"));
            }
        }
        h.samples = cols.split('\t').skip(9).map(str::to_owned).collect();

        // htslib always injects PASS at zero before processing declarations,
        // including headers whose other dictionary entries carry explicit IDX.
        h.declare("PASS", Kind::Filter, None, Some(0));
        for (name, kind, ty, idx) in declared {
            if !h.key_ids.contains_key(&name) {
                let id = idx.map_or(h.keys.len(), |i| i as usize);
                if h.keys.get(id).is_some_and(Option::is_some) {
                    return Err(Error::invalid("vcf header", "conflicting dictionary IDX"));
                }
                if id >= h.keys.len() {
                    h.keys.try_reserve(id + 1 - h.keys.len()).map_err(|_| {
                        Error::invalid("vcf header", "dictionary allocation failed")
                    })?;
                }
            }
            h.declare(&name, kind, ty, idx);
        }
        for (contig, idx) in contigs {
            if h.contig_ids.contains_key(&contig.name) {
                continue;
            }
            let id = idx.unwrap_or(h.contigs.len() as u32);
            if h.contigs
                .get(id as usize)
                .is_some_and(|c| !c.name.is_empty())
            {
                return Err(Error::invalid("vcf header", "conflicting contig IDX"));
            }
            if id as usize >= h.contigs.len() {
                h.contigs
                    .try_reserve(id as usize + 1 - h.contigs.len())
                    .map_err(|_| Error::invalid("vcf header", "contig allocation failed"))?;
            }
            while h.contigs.len() <= id as usize {
                h.contigs.push(Contig {
                    name: String::new(),
                    len: None,
                });
            }
            h.contig_ids.insert(contig.name.clone(), id);
            h.contigs[id as usize] = contig;
        }
        Ok(h)
    }

    fn declare(
        &mut self,
        name: &str,
        kind: Kind,
        ty: Option<(ValueType, Number)>,
        idx: Option<u32>,
    ) -> u32 {
        let id = match self.key_ids.get(name) {
            Some(&id) => id,
            None => {
                let id = idx.unwrap_or(self.keys.len() as u32);
                if self.keys.len() <= id as usize {
                    self.keys.resize(id as usize + 1, None);
                }
                self.keys[id as usize] = Some(Key {
                    name: name.to_owned(),
                    ..Key::default()
                });
                self.key_ids.insert(name.to_owned(), id);
                id
            }
        };
        let key = self.keys[id as usize].as_mut().unwrap();
        match (kind, ty) {
            (Kind::Info, Some((ty, number))) if key.info.is_none() => {
                key.info = Some(ty);
                key.info_number = number;
            }
            (Kind::Format, Some((ty, number))) if key.format.is_none() => {
                key.format = Some(ty);
                key.format_number = number;
            }
            (Kind::Filter, _) => key.filter = true,
            _ => {}
        }
        id
    }

    pub fn contig_id(&self, name: &str) -> Option<u32> {
        self.contig_ids.get(name).copied()
    }

    pub fn key_id(&self, name: &str) -> Option<u32> {
        self.key_ids.get(name).copied()
    }

    pub fn key(&self, id: u32) -> Option<&Key> {
        self.keys.get(id as usize).and_then(Option::as_ref)
    }

    pub fn info_type(&self, name: &str) -> Option<ValueType> {
        self.key(self.key_id(name)?)?.info
    }

    pub fn format_type(&self, name: &str) -> Option<ValueType> {
        self.key(self.key_id(name)?)?.format
    }

    /// A contig met in a record but not declared, added as htslib does.
    pub(crate) fn add_contig(&mut self, name: &str) -> u32 {
        if let Some(id) = self.contig_id(name) {
            return id;
        }
        let id = self.contigs.len() as u32;
        self.contigs.push(Contig {
            name: name.to_owned(),
            len: None,
        });
        self.contig_ids.insert(name.to_owned(), id);
        id
    }

    /// An INFO or FORMAT key met in a record but not declared for its kind,
    /// added as `Type=String` as htslib does; returns its id.
    pub(crate) fn add_undeclared(&mut self, name: &str, format: bool) -> u32 {
        let kind = if format { Kind::Format } else { Kind::Info };
        self.declare(
            name,
            kind,
            Some((ValueType::String, Number::Fixed(1))),
            None,
        )
    }

    /// A FILTER met in a record that is not a dictionary key yet, added as
    /// htslib does; returns its id. A name already in the dictionary (as
    /// any kind) keeps its id unchanged.
    pub(crate) fn add_undeclared_filter(&mut self, name: &str) -> u32 {
        match self.key_id(name) {
            Some(id) => id,
            None => self.declare(name, Kind::Filter, None, None),
        }
    }

    pub fn info_number(&self, name: &str) -> Option<Number> {
        let key = self.key(self.key_id(name)?)?;
        key.info.map(|_| key.info_number)
    }

    pub fn format_number(&self, name: &str) -> Option<Number> {
        let key = self.key(self.key_id(name)?)?;
        key.format.map(|_| key.format_number)
    }

    /// The index of a sample by name.
    pub fn sample_id(&self, name: &str) -> Option<usize> {
        self.samples.iter().position(|s| s == name)
    }

    /// The packed declaration word htslib keeps per key and line kind
    /// (`hdr->id[BCF_DT_ID][id].val->info[kind]`): number << 12 |
    /// length << 8 | type << 4 | kind, or 15 for a kind the key is not
    /// declared as. `None` for an id outside the dictionary.
    pub fn declaration_word(&self, id: u32, kind: LineKind) -> Option<u32> {
        let key = self.key(id)?;
        let word = |ty: Option<ValueType>, number: Number| match ty {
            None => 15,
            Some(ty) => {
                number.count() << 12 | number.length_code() << 8 | ty.code() << 4 | kind as u32
            }
        };
        Some(match kind {
            LineKind::Info => word(key.info, key.info_number),
            LineKind::Format => word(key.format, key.format_number),
            LineKind::Filter if key.filter => 0xfffff << 12 | 0xf << 8 | 0xf << 4,
            LineKind::Filter => 15,
        })
    }

    /// `bcf_hdr_idinfo_exists`.
    pub fn declared_as(&self, id: u32, kind: LineKind) -> bool {
        self.declaration_word(id, kind)
            .is_some_and(|w| w & 0xf != 0xf)
    }

    /// Whether records can carry genotypes: htslib reads GT only when it is
    /// a `String` FORMAT field (declared so, or added undeclared).
    pub fn has_gt_field(&self) -> bool {
        self.format_type("GT") == Some(ValueType::String)
    }
}
