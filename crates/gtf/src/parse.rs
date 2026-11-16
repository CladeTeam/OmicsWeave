//! Parsing a batch of lines into raw columns with batch-local string ids.
//!
//! Batches parse independently on worker threads; [`Raw::append`] then
//! moves each batch's strings into the global pool in batch order, so the
//! result does not depend on the number of threads.

use std::sync::Arc;

use omicsweave_core::{Error, Result};

use crate::model::{AttributePolicy, Dialect, PHASE_NONE, Strand, UNQUOTED};
use crate::strings::{IdMap, NONE, StringPool};

/// Slots of [`Raw::std`]: the attributes the model reads whatever the
/// attribute policy.
pub(crate) const GENE_ID: usize = 0;
pub(crate) const TRANSCRIPT_ID: usize = 1;
pub(crate) const ID: usize = 2;
pub(crate) const GENE_NAME: usize = 3;
pub(crate) const GENE_BIOTYPE: usize = 4;
pub(crate) const TX_NAME: usize = 5;
pub(crate) const TX_BIOTYPE: usize = 6;
const STD: usize = 7;

/// File records as columns, in file order.
#[derive(Default)]
pub(crate) struct Raw {
    pub pool: StringPool,
    pub seqid: Vec<u32>,
    pub source: Vec<u32>,
    pub feature_type: Vec<u32>,
    pub score: Vec<u32>,
    pub start: Vec<u32>,
    pub end: Vec<u32>,
    pub strand: Vec<Strand>,
    pub phase: Vec<u8>,
    pub std: Vec<[u32; STD]>,
    pub parents: Vec<u32>,
    pub parent_offsets: Vec<u32>,
    pub tags: Vec<u32>,
    pub tag_offsets: Vec<u32>,
    /// Distinct (key, value) pairs; keys carry the quoting bit.
    pub pairs: Vec<(u32, u32)>,
    pub pair_ids: IdMap<u64, u32>,
    /// Pair ids of each record, `attr_offsets` delimited.
    pub attrs: Vec<u32>,
    pub attr_offsets: Vec<u32>,
    pub line: Vec<u64>,
    /// (records before the line, text of the line) for `#` lines.
    pub comments: Vec<(u32, u32)>,
    /// A `##FASTA` directive ended the records.
    pub fasta: bool,
}

/// What the parser keeps of column 9.
#[derive(Clone)]
pub(crate) enum Keep {
    All,
    Only(Arc<[Box<[u8]>]>),
    None,
}

impl Keep {
    pub fn new(p: &AttributePolicy) -> Keep {
        match p {
            AttributePolicy::All => Keep::All,
            AttributePolicy::None => Keep::None,
            AttributePolicy::Only(keys) => Keep::Only(
                keys.iter()
                    .map(|k| k.as_bytes().to_vec().into_boxed_slice())
                    .collect(),
            ),
        }
    }

    #[inline]
    fn keeps(&self, key: &[u8]) -> bool {
        match self {
            Keep::All => true,
            Keep::None => false,
            Keep::Only(keys) => keys.iter().any(|k| &k[..] == key),
        }
    }
}

impl Raw {
    pub fn new() -> Raw {
        Raw {
            parent_offsets: vec![0],
            tag_offsets: vec![0],
            attr_offsets: vec![0],
            ..Raw::default()
        }
    }

    pub fn len(&self) -> usize {
        self.start.len()
    }

    /// The id of pair `(key, value)`, adding it if new.
    #[inline]
    pub fn pair(&mut self, kv: (u32, u32)) -> Result<u32> {
        let n = self.pairs.len();
        if n >= NONE as usize {
            return Err(Error::Unsupported("more than 2^32 attribute pairs"));
        }
        let id = *self
            .pair_ids
            .entry(((kv.0 as u64) << 32) | kv.1 as u64)
            .or_insert(n as u32);
        if id as usize == n {
            self.pairs.push(kv);
        }
        Ok(id)
    }

    pub fn parents(&self, r: usize) -> &[u32] {
        &self.parents[self.parent_offsets[r] as usize..self.parent_offsets[r + 1] as usize]
    }

    pub fn tags(&self, r: usize) -> &[u32] {
        &self.tags[self.tag_offsets[r] as usize..self.tag_offsets[r + 1] as usize]
    }

    /// Moves the records of `part` (a later batch) to the end, re-interning
    /// its strings in this pool.
    pub fn append(&mut self, part: Raw) -> Result<()> {
        if self.len().saturating_add(part.len()) > crate::model::ROW_LIMIT as usize {
            return Err(Error::Unsupported("more than 2^30 annotation records"));
        }
        let mut map = Vec::with_capacity(part.pool.len());
        for id in 0..part.pool.len() as u32 {
            map.push(self.pool.intern(part.pool.bytes(id))?);
        }
        let m = |id: u32| if id == NONE { NONE } else { map[id as usize] };
        let base = u32::try_from(self.len())
            .ok()
            .filter(|&b| b < crate::model::ROW_LIMIT)
            .ok_or(Error::Unsupported("more than 2^30 annotation records"))?;
        for (records, text) in part.comments {
            self.comments.push((base + records, m(text)));
        }
        self.seqid.extend(part.seqid.iter().map(|&i| m(i)));
        self.source.extend(part.source.iter().map(|&i| m(i)));
        self.feature_type
            .extend(part.feature_type.iter().map(|&i| m(i)));
        self.score.extend(part.score.iter().map(|&i| m(i)));
        self.start.extend_from_slice(&part.start);
        self.end.extend_from_slice(&part.end);
        self.strand.extend_from_slice(&part.strand);
        self.phase.extend_from_slice(&part.phase);
        self.std.extend(part.std.iter().map(|s| s.map(m)));
        self.line.extend_from_slice(&part.line);
        append_csr(
            &mut self.parents,
            &mut self.parent_offsets,
            &part.parents,
            &part.parent_offsets,
            m,
        )?;
        append_csr(
            &mut self.tags,
            &mut self.tag_offsets,
            &part.tags,
            &part.tag_offsets,
            m,
        )?;
        let mut pmap = Vec::with_capacity(part.pairs.len());
        for &(k, v) in &part.pairs {
            pmap.push(self.pair((m(k & !UNQUOTED) | (k & UNQUOTED), m(v)))?);
        }
        append_csr(
            &mut self.attrs,
            &mut self.attr_offsets,
            &part.attrs,
            &part.attr_offsets,
            |p| pmap[p as usize],
        )?;
        self.fasta |= part.fasta;
        Ok(())
    }
}

fn append_csr(
    values: &mut Vec<u32>,
    offsets: &mut Vec<u32>,
    part_values: &[u32],
    part_offsets: &[u32],
    map: impl Fn(u32) -> u32,
) -> Result<()> {
    let base = values.len();
    if base + part_values.len() > u32::MAX as usize {
        return Err(Error::Unsupported("more than 2^32 attribute values"));
    }
    values.extend(part_values.iter().map(|&v| map(v)));
    offsets.extend(part_offsets[1..].iter().map(|&o| base as u32 + o));
    Ok(())
}

/// The dialect of a file from its first lines: a `##gff-version 3`
/// directive, else the attribute syntax of the first record.
pub(crate) fn sniff_dialect(text: &[u8]) -> Dialect {
    for line in text.split(|&c| c == b'\n') {
        if line.starts_with(b"##gff-version") {
            let v = line[13..].trim_ascii();
            if v.starts_with(b"3") {
                return Dialect::Gff3;
            }
            continue;
        }
        if line.is_empty() || line[0] == b'#' {
            continue;
        }
        let Some(attrs) = line.splitn(9, |&c| c == b'\t').nth(8) else {
            continue;
        };
        let eq = memchr::memchr(b'=', attrs);
        let gtf = memchr::memchr2(b' ', b'"', attrs);
        return match (eq, gtf) {
            (Some(e), Some(g)) if e < g => Dialect::Gff3,
            (Some(_), None) => Dialect::Gff3,
            _ => Dialect::Gtf,
        };
    }
    Dialect::Gtf
}

fn bad(line: u64, detail: impl Into<String>) -> Error {
    Error::Format {
        line,
        detail: detail.into(),
    }
}

/// A 1-based coordinate.
fn coordinate(b: &[u8], line: u64, what: &str) -> Result<u32> {
    if b.is_empty() || b.len() > 10 {
        return Err(bad(
            line,
            format!("bad {what} {:?}", String::from_utf8_lossy(b)),
        ));
    }
    let mut v: u64 = 0;
    for &c in b {
        if !c.is_ascii_digit() {
            return Err(bad(
                line,
                format!("bad {what} {:?}", String::from_utf8_lossy(b)),
            ));
        }
        v = v * 10 + (c - b'0') as u64;
    }
    u32::try_from(v).map_err(|_| bad(line, format!("{what} {v} exceeds 2^32 - 1")))
}

struct Fields<'a> {
    std: [u32; STD],
    rank: [u8; STD],
    pool: &'a mut StringPool,
}

impl Fields<'_> {
    /// Fills `slot` with `value` unless a key of higher `rank` already did.
    #[inline]
    fn set(&mut self, slot: usize, rank: u8, value: &[u8]) -> Result<()> {
        if rank > self.rank[slot] && !value.is_empty() {
            self.std[slot] = self.pool.intern(value)?;
            self.rank[slot] = rank;
        }
        Ok(())
    }
}

/// Parses one batch of lines.
pub(crate) fn parse_batch(
    text: &[u8],
    first_line: u64,
    dialect: Dialect,
    keep: &Keep,
) -> Result<Raw> {
    let mut raw = Raw::new();
    let mut at = 0;
    let mut no = first_line;
    while at < text.len() {
        let end = memchr::memchr(b'\n', &text[at..]).map_or(text.len(), |i| at + i);
        let mut line = &text[at..end];
        at = end + 1;
        let line_no = no;
        no += 1;
        if let [rest @ .., b'\r'] = line {
            line = rest;
        }
        if line.is_empty() {
            continue;
        }
        if line[0] == b'#' || line.starts_with(b"track ") || line.starts_with(b"browser ") {
            if line.starts_with(b"##FASTA") {
                raw.fasta = true;
                break;
            }
            let id = raw.pool.intern(line)?;
            raw.comments.push((raw.len() as u32, id));
            continue;
        }
        parse_record(&mut raw, line, line_no, dialect, keep)?;
    }
    Ok(raw)
}

fn parse_record(raw: &mut Raw, line: &[u8], no: u64, dialect: Dialect, keep: &Keep) -> Result<()> {
    if raw.len() >= crate::model::ROW_LIMIT as usize {
        return Err(Error::Unsupported("more than 2^30 annotation records"));
    }
    let mut cols: [&[u8]; 9] = [b""; 9];
    let mut n = 0;
    let mut rest = line;
    while n < 8 {
        match memchr::memchr(b'\t', rest) {
            Some(i) => {
                cols[n] = &rest[..i];
                rest = &rest[i + 1..];
                n += 1;
            }
            None => {
                cols[n] = rest;
                rest = b"";
                n += 1;
                break;
            }
        }
    }
    if n == 8 && !rest.is_empty() {
        cols[8] = rest;
        n = 9;
    } else if n == 8 && line.ends_with(b"\t") {
        n = 9;
    }
    if n < 8 {
        return Err(bad(no, format!("{n} columns, need 8 or 9")));
    }
    let start = coordinate(cols[3], no, "start")?;
    let end = coordinate(cols[4], no, "end")?;
    if start == 0 {
        return Err(bad(no, "start must be at least 1"));
    }
    if end < start {
        return Err(bad(no, format!("end {end} before start {start}")));
    }
    let strand = match cols[6] {
        [c] => Strand::from_byte(*c),
        _ => None,
    }
    .ok_or_else(|| {
        bad(
            no,
            format!("bad strand {:?}", String::from_utf8_lossy(cols[6])),
        )
    })?;
    let phase = match cols[7] {
        b"." => PHASE_NONE,
        [c @ b'0'..=b'2'] => c - b'0',
        p => {
            return Err(bad(
                no,
                format!("bad phase {:?}", String::from_utf8_lossy(p)),
            ));
        }
    };
    let score = match cols[5] {
        b"." => NONE,
        s => {
            let ok = std::str::from_utf8(s)
                .ok()
                .and_then(|t| t.parse::<f64>().ok())
                .is_some();
            if !ok {
                return Err(bad(
                    no,
                    format!("bad score {:?}", String::from_utf8_lossy(s)),
                ));
            }
            raw.pool.intern(s)?
        }
    };
    let pool = &mut raw.pool;
    raw.seqid.push(pool.intern(cols[0])?);
    raw.source.push(pool.intern(cols[1])?);
    raw.feature_type.push(pool.intern(cols[2])?);
    raw.score.push(score);
    raw.start.push(start - 1);
    raw.end.push(end);
    raw.strand.push(strand);
    raw.phase.push(phase);
    raw.line.push(no);
    let mut f = Fields {
        std: [NONE; STD],
        rank: [0; STD],
        pool: &mut raw.pool,
    };
    let attrs = match cols[8] {
        b"." => b"".as_slice(),
        a => a,
    };
    let (pairs, pair_ids, attrs_out) = (&mut raw.pairs, &mut raw.pair_ids, &mut raw.attrs);
    let mut push_pair =
        |f: &mut Fields<'_>, key: &[u8], value: &[u8], quoted: bool| -> Result<()> {
            if keep.keeps(key) {
                let k = f.pool.intern(key)?;
                let v = f.pool.intern(value)?;
                let kv = (if quoted { k } else { k | UNQUOTED }, v);
                let n = u32::try_from(pairs.len())
                    .ok()
                    .filter(|&n| n != NONE)
                    .ok_or(Error::Unsupported("more than 2^32 attribute pairs"))?;
                let id = *pair_ids
                    .entry(((kv.0 as u64) << 32) | kv.1 as u64)
                    .or_insert(n);
                if id == n {
                    pairs.push(kv);
                }
                attrs_out.push(id);
            }
            Ok(())
        };
    match dialect {
        Dialect::Gtf => {
            let b = attrs;
            let mut i = 0;
            loop {
                while i < b.len() && matches!(b[i], b' ' | b'\t' | b';') {
                    i += 1;
                }
                if i >= b.len() {
                    break;
                }
                let ks = i;
                while i < b.len() && !matches!(b[i], b' ' | b'\t' | b';') {
                    i += 1;
                }
                let key = &b[ks..i];
                while i < b.len() && matches!(b[i], b' ' | b'\t') {
                    i += 1;
                }
                let (value, quoted) = if i < b.len() && b[i] == b'"' {
                    let vs = i + 1;
                    let ve = memchr::memchr(b'"', &b[vs..])
                        .map(|j| vs + j)
                        .ok_or_else(|| bad(no, "unterminated quoted attribute value"))?;
                    i = ve + 1;
                    (&b[vs..ve], true)
                } else {
                    let vs = i;
                    while i < b.len() && b[i] != b';' {
                        i += 1;
                    }
                    (b[vs..i].trim_ascii_end(), false)
                };
                while i < b.len() && b[i] != b';' {
                    i += 1;
                }
                gtf_field(&mut f, key, value, &mut raw.tags)?;
                push_pair(&mut f, key, value, quoted)?;
            }
        }
        Dialect::Gff3 => {
            for piece in attrs.split(|&c| c == b';') {
                let piece = piece.trim_ascii();
                if piece.is_empty() {
                    continue;
                }
                let (key, value) = match memchr::memchr(b'=', piece) {
                    Some(e) => (piece[..e].trim_ascii(), piece[e + 1..].trim_ascii()),
                    None => (piece, b"".as_slice()),
                };
                gff3_field(&mut f, key, value, &mut raw.parents, &mut raw.tags)?;
                push_pair(&mut f, key, value, true)?;
            }
        }
    }
    let std = f.std;
    raw.std.push(std);
    let offset = |n: usize| {
        u32::try_from(n).map_err(|_| Error::Unsupported("more than 2^32 attribute values"))
    };
    raw.parent_offsets.push(offset(raw.parents.len())?);
    raw.tag_offsets.push(offset(raw.tags.len())?);
    raw.attr_offsets.push(offset(raw.attrs.len())?);
    Ok(())
}

fn gtf_field(f: &mut Fields<'_>, key: &[u8], value: &[u8], tags: &mut Vec<u32>) -> Result<()> {
    match key {
        b"gene_id" => f.set(GENE_ID, 1, value),
        b"transcript_id" => f.set(TRANSCRIPT_ID, 1, value),
        b"gene_name" => f.set(GENE_NAME, 3, value),
        b"gene" => f.set(GENE_NAME, 1, value),
        b"gene_type" => f.set(GENE_BIOTYPE, 3, value),
        b"gene_biotype" => f.set(GENE_BIOTYPE, 2, value),
        b"transcript_name" => f.set(TX_NAME, 3, value),
        b"transcript_type" => f.set(TX_BIOTYPE, 3, value),
        b"transcript_biotype" => f.set(TX_BIOTYPE, 2, value),
        b"tag" if !value.is_empty() => {
            tags.push(f.pool.intern(value)?);
            Ok(())
        }
        _ => Ok(()),
    }
}

fn gff3_field(
    f: &mut Fields<'_>,
    key: &[u8],
    value: &[u8],
    parents: &mut Vec<u32>,
    tags: &mut Vec<u32>,
) -> Result<()> {
    match key {
        b"ID" => f.set(ID, 1, value),
        b"gene_id" => f.set(GENE_ID, 1, value),
        b"transcript_id" => f.set(TRANSCRIPT_ID, 1, value),
        b"Name" => {
            f.set(GENE_NAME, 2, value)?;
            f.set(TX_NAME, 2, value)
        }
        b"gene_name" => f.set(GENE_NAME, 3, value),
        b"gene" => f.set(GENE_NAME, 1, value),
        b"transcript_name" => f.set(TX_NAME, 3, value),
        b"biotype" => {
            f.set(GENE_BIOTYPE, 1, value)?;
            f.set(TX_BIOTYPE, 1, value)
        }
        b"gene_biotype" => f.set(GENE_BIOTYPE, 3, value),
        b"gene_type" => f.set(GENE_BIOTYPE, 2, value),
        b"transcript_biotype" => f.set(TX_BIOTYPE, 3, value),
        b"transcript_type" => f.set(TX_BIOTYPE, 2, value),
        b"Parent" => {
            for p in value.split(|&c| c == b',') {
                if !p.is_empty() {
                    parents.push(f.pool.intern(p)?);
                }
            }
            Ok(())
        }
        b"tag" => {
            for t in value.split(|&c| c == b',') {
                if !t.is_empty() {
                    tags.push(f.pool.intern(t)?);
                }
            }
            Ok(())
        }
        _ => Ok(()),
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn gtf_attributes_keep_quoting_and_order() {
        let line = b"chr1\tHAVANA\texon\t11869\t12227\t.\t+\t.\tgene_id \"G1\"; transcript_id \"T1\"; level 2; tag \"basic\"; tag \"CCDS\"; note \"a; b\";\n";
        let raw = parse_batch(line, 1, Dialect::Gtf, &Keep::All).unwrap();
        assert_eq!(raw.len(), 1);
        assert_eq!((raw.start[0], raw.end[0]), (11868, 12227));
        let p = &raw.pool;
        let pairs: Vec<(String, String, bool)> = raw
            .attrs
            .iter()
            .map(|&id| raw.pairs[id as usize])
            .map(|(k, v)| {
                (
                    p.str(k & !UNQUOTED).to_string(),
                    p.str(v).to_string(),
                    k & UNQUOTED == 0,
                )
            })
            .collect();
        assert_eq!(pairs[2], ("level".into(), "2".into(), false));
        assert_eq!(pairs[5], ("note".into(), "a; b".into(), true));
        assert_eq!(raw.tags(0).len(), 2);
        assert_eq!(p.str(raw.std[0][GENE_ID]), "G1");
    }

    #[test]
    fn sniffs_dialects() {
        assert_eq!(sniff_dialect(b"##gff-version 3\n"), Dialect::Gff3);
        assert_eq!(
            sniff_dialect(b"#c\nchr1\ts\tgene\t1\t2\t.\t+\t.\tID=g;Name=x\n"),
            Dialect::Gff3
        );
        assert_eq!(
            sniff_dialect(b"chr1\ts\texon\t1\t2\t.\t+\t.\tgene_id \"a=b\";\n"),
            Dialect::Gtf
        );
    }
}
