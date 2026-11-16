//! Base modifications: the `MM` / `ML` / `MN` tags (SAMtags §1.7), read as
//! htslib's `bam_parse_basemod2` and `bam_next_basemod` read them.
//!
//! - [`ModTags`] locates the tags among a record's optional fields (BAM,
//!   or a CRAM record's fields in their BAM encoding), `MM` before the
//!   older `Mm`, `ML` before `Ml`.
//! - [`ModDecoder`] resolves them against the record's bases into
//!   [`ModificationEntry`] rows: every listed base of every entry and code,
//!   in SEQ order (htslib's order: by position, then entry, then code).
//!   Positions count the canonical bases of the read as sequenced, so a
//!   reverse-strand read is walked on its complement from the right end.
//! - [`reference_positions`] then places the rows on the reference
//!   through the CIGAR.
//!
//! The decoder rejects what htslib rejects (an unknown canonical base or
//! strand, a missing `;`, more listed bases than the read has, too few or
//! too many `ML` values, an `MN` other than the SEQ length, `MM` not of
//! type `Z`, `ML` not `B:C`); a delta that is not an unsigned integer is
//! also an error, where htslib may stop reading the list silently.

use omicsweave_core::{Error, ModificationEntry, Result};

use crate::record::aux_get;

/// `seq_nt16_table` for the bases of an ASCII SEQ.
const fn nt16_table() -> [u8; 256] {
    let mut t = [15u8; 256];
    let codes = b"=ACMGRSVTWYHKDBN";
    let mut i = 0;
    while i < 16 {
        t[codes[i] as usize] = i as u8;
        t[codes[i].to_ascii_lowercase() as usize] = i as u8;
        i += 1;
    }
    t[b'U' as usize] = 8;
    t[b'u' as usize] = 8;
    t[b'0' as usize] = 1;
    t[b'1' as usize] = 2;
    t[b'2' as usize] = 4;
    t[b'3' as usize] = 8;
    t
}

const NT16_CODE: [u8; 256] = nt16_table();

/// `seqi_rc`: the complement of each 4-bit base code.
const NT16_RC: [u8; 16] = [0, 8, 4, 12, 2, 10, 6, 14, 1, 9, 5, 13, 3, 11, 7, 15];

/// htslib keeps at most this many modification types per record.
const MAX_TYPES: usize = 256;

/// The base-modification tags of one record.
#[derive(Debug, Clone, Copy)]
pub struct ModTags<'a> {
    /// The `MM` string, without its NUL.
    mm: &'a [u8],
    /// The `ML` values.
    ml: Option<&'a [u8]>,
    /// `MN` as `bam_aux2i` reads it: 0 for a value that is not an integer.
    mn: Option<i64>,
}

fn malformed(detail: &'static str) -> BadModTags {
    BadModTags(detail)
}

/// Why a record's modification tags cannot be read.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct BadModTags(pub &'static str);

impl BadModTags {
    /// The error of record `name`.
    pub fn error(self, name: &[u8]) -> Error {
        Error::invalid(
            "MM/ML tags",
            format!("{}: {}", String::from_utf8_lossy(name), self.0),
        )
    }
}

impl std::fmt::Display for BadModTags {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.write_str(self.0)
    }
}

/// The value of `tag` (its type byte first) in the first piece holding it.
fn find<'a>(
    pieces: &[&'a [u8]],
    tag: [u8; 2],
) -> std::result::Result<Option<&'a [u8]>, BadModTags> {
    for piece in pieces {
        match aux_get(piece, tag) {
            Ok(Some(at)) => return Ok(Some(&piece[at..])),
            Ok(None) => {}
            Err(()) => return Err(malformed("corrupted optional fields")),
        }
    }
    Ok(None)
}

fn find_either<'a>(
    pieces: &[&'a [u8]],
    tags: [[u8; 2]; 2],
) -> std::result::Result<Option<&'a [u8]>, BadModTags> {
    Ok(match find(pieces, tags[0])? {
        Some(v) => Some(v),
        None => find(pieces, tags[1])?,
    })
}

/// `bam_aux2i`: an integer value, else 0.
fn aux_int(v: &[u8]) -> i64 {
    let le = |n: usize| -> Option<&[u8]> { v.get(1..1 + n) };
    match v[0] {
        b'c' => le(1).map_or(0, |b| b[0] as i8 as i64),
        b'C' => le(1).map_or(0, |b| b[0] as i64),
        b's' => le(2).map_or(0, |b| i16::from_le_bytes([b[0], b[1]]) as i64),
        b'S' => le(2).map_or(0, |b| u16::from_le_bytes([b[0], b[1]]) as i64),
        b'i' => le(4).map_or(0, |b| i32::from_le_bytes(b.try_into().unwrap()) as i64),
        b'I' => le(4).map_or(0, |b| u32::from_le_bytes(b.try_into().unwrap()) as i64),
        _ => 0,
    }
}

impl<'a> ModTags<'a> {
    /// The tags among `pieces`, a record's optional fields in BAM
    /// encoding (in one or more pieces); `None` when it has no `MM` / `Mm`.
    pub fn find(pieces: &[&'a [u8]]) -> std::result::Result<Option<ModTags<'a>>, BadModTags> {
        let Some(mm) = find_either(pieces, [*b"MM", *b"Mm"])? else {
            return Ok(None);
        };
        if mm[0] != b'Z' {
            return Err(malformed("MM is not of type Z"));
        }
        let body = &mm[1..];
        let end = body
            .iter()
            .position(|&c| c == 0)
            .ok_or(malformed("MM is not NUL-terminated"))?;
        let mm = &body[..end];
        let ml = match find_either(pieces, [*b"ML", *b"Ml"])? {
            None => None,
            Some(v) => {
                if v.len() < 6 || v[0] != b'B' || v[1] != b'C' {
                    return Err(malformed("ML is not of type B:C"));
                }
                let n = u32::from_le_bytes(v[2..6].try_into().unwrap()) as usize;
                let end = 6usize
                    .checked_add(n)
                    .ok_or(malformed("ML is shorter than its count"))?;
                Some(
                    v.get(6..end)
                        .ok_or(malformed("ML is shorter than its count"))?,
                )
            }
        };
        let mn = find(pieces, *b"MN")?.map(aux_int);
        Ok(Some(ModTags { mm, ml, mn }))
    }

    /// The `MM` string.
    pub fn mm(&self) -> &'a [u8] {
        self.mm
    }

    /// The `ML` values, when the record has them.
    pub fn ml(&self) -> Option<&'a [u8]> {
        self.ml
    }

    /// `MN` as an integer (0 when it is not one), when present.
    pub fn mn(&self) -> Option<i64> {
        self.mn
    }
}

/// The bases of a record, as stored.
#[derive(Debug, Clone, Copy)]
pub enum ModSeq<'a> {
    /// BAM's packing: two 4-bit codes per byte, high nibble first.
    Packed { bytes: &'a [u8], len: usize },
    /// One ASCII base per position (`=ACMGRSVTWYHKDBN`, either case);
    /// htslib's aliases `U`/`u` for T and `0`..`3` for A/C/G/T are accepted.
    Ascii(&'a [u8]),
}

impl ModSeq<'_> {
    fn len(&self) -> usize {
        match self {
            ModSeq::Packed { len, .. } => *len,
            ModSeq::Ascii(b) => b.len(),
        }
    }

    #[inline]
    unsafe fn code_unchecked(&self, i: usize) -> u8 {
        match self {
            // SAFETY: the caller proves i < len and decode validates the
            // packed byte count before any sequence walk. Both nibbles
            // are in 0..16, including the unused last low nibble.
            ModSeq::Packed { bytes, .. } => unsafe {
                (*bytes.get_unchecked(i / 2) >> (4 * (1 - i % 2))) & 0xf
            },
            // SAFETY: i is within the ASCII sequence and any byte indexes
            // the complete 256-entry conversion table.
            ModSeq::Ascii(b) => unsafe { *NT16_CODE.get_unchecked(*b.get_unchecked(i) as usize) },
        }
    }
}

/// One modification type of a record's `MM` tag: a code of one entry.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct ModType {
    /// As [`ModificationEntry::code`].
    pub code: i32,
    /// As [`ModificationEntry::canonical`].
    pub canonical: u8,
    /// Whether the entry is for the opposite strand (`-`).
    pub minus: bool,
    /// Whether the entry is explicit (`?`).
    pub explicit: bool,
}

/// Decodes `MM` / `ML` into [`ModificationEntry`] rows; holds scratch
/// reused from record to record.
#[derive(Debug, Default, Clone)]
pub struct ModDecoder {
    deltas: Vec<u32>,
    types: Vec<ModType>,
    selected: Vec<(usize, i32)>,
    positions: [Vec<u32>; 4],
}

/// `hts_str2uint(.., 31, ..)` at `s[*i..]`: an optional `+` then at least
/// one digit, below 2^31.
fn uint31(s: &[u8], i: &mut usize) -> Option<u32> {
    if s.get(*i) == Some(&b'+') {
        *i += 1;
    }
    let start = *i;
    let mut n: u64 = 0;
    while let Some(&c) = s.get(*i)
        && c.is_ascii_digit()
    {
        n = n * 10 + (c - b'0') as u64;
        if n > i32::MAX as u64 {
            return None;
        }
        *i += 1;
    }
    (*i > start).then_some(n as u32)
}

impl ModDecoder {
    pub fn new() -> ModDecoder {
        ModDecoder::default()
    }

    /// The modification types of the last decoded record, in `MM` order
    /// (htslib's `bam_mods_recorded` / `bam_mods_queryi`), filtered or not.
    /// Empty after an error, or an untagged [`decode_record`] call.
    pub fn types(&self) -> &[ModType] {
        &self.types
    }

    /// Appends the modifications of a record with `tags` and bases `seq`
    /// (reverse-complemented by the aligner when `reverse`) to `out`, in
    /// SEQ order with `ref_pos` -1, keeping the codes in `codes` (all of
    /// them when empty). On error, `out` is left as it was.
    pub fn decode(
        &mut self,
        tags: &ModTags<'_>,
        seq: ModSeq<'_>,
        reverse: bool,
        codes: &[i32],
        out: &mut Vec<ModificationEntry>,
    ) -> std::result::Result<(), BadModTags> {
        let first = out.len();
        let result = self.decode_into(tags, seq, reverse, codes, out);
        if result.is_err() {
            out.truncate(first);
            self.types.clear();
        }
        result
    }

    fn decode_into(
        &mut self,
        tags: &ModTags<'_>,
        seq: ModSeq<'_>,
        reverse: bool,
        codes: &[i32],
        out: &mut Vec<ModificationEntry>,
    ) -> std::result::Result<(), BadModTags> {
        self.types.clear();
        let l = seq.len();
        if l > i32::MAX as usize {
            return Err(malformed("SEQ is too long for modification positions"));
        }
        if let ModSeq::Packed { bytes, .. } = seq
            && bytes.len() < l.div_ceil(2)
        {
            return Err(malformed("packed SEQ is shorter than its length"));
        }
        if let Some(mn) = tags.mn
            && l != 0
            && mn != l as i64
        {
            return Err(malformed("MN differs from the SEQ length"));
        }
        let mm = tags.mm;
        // Multiple entries share a lazily extended canonical-base prefix.
        // Scratch lengths reset per record while allocations are reused.
        let mut multiple_groups = false;
        for p in &mut self.positions {
            p.clear();
        }
        let mut scanned = 0usize;
        let first = out.len();
        let mut ml_at = 0usize;
        let mut hit_groups = 0usize;
        let mut i = 0;
        while i < mm.len() {
            let (canonical, base) = match mm[i] {
                b'A' => (1, b'A'),
                b'C' => (2, b'C'),
                b'G' => (4, b'G'),
                b'T' | b'U' => (8, b'T'),
                b'N' => (15, b'N'),
                _ => return Err(malformed("unknown canonical base in MM")),
            };
            i += 1;
            let minus = match mm.get(i) {
                Some(b'+') => false,
                Some(b'-') => true,
                _ => return Err(malformed("MM strand is neither + nor -")),
            };
            i += 1;
            let codes_at = i;
            let chebi = if mm.get(i).is_some_and(u8::is_ascii_digit) {
                let v = uint31(mm, &mut i).ok_or(malformed("invalid ChEBI code in MM"))?;
                Some(v)
            } else {
                while mm.get(i).is_some_and(u8::is_ascii_alphabetic) {
                    i += 1;
                }
                if i == mm.len() {
                    return Err(malformed("MM ends without a semicolon"));
                }
                None
            };
            let n_types = if chebi.is_some() { 1 } else { i - codes_at };
            let code_of = |t: usize| -> i32 {
                match chebi {
                    // htslib reads ChEBI 0 as the code of its last digit.
                    Some(0) => b'0' as i32,
                    Some(v) => -(v as i32),
                    None => mm[codes_at + t] as i32,
                }
            };
            let explicit = match mm.get(i) {
                Some(b'.') => {
                    i += 1;
                    false
                }
                Some(b'?') => {
                    i += 1;
                    true
                }
                Some(b',' | b';') => false,
                _ => return Err(malformed("unexpected character in MM")),
            };
            self.deltas.clear();
            loop {
                match mm.get(i) {
                    Some(b';') => {
                        i += 1;
                        break;
                    }
                    Some(b',') => {
                        i += 1;
                        let d = uint31(mm, &mut i).ok_or(malformed("invalid delta in MM"))?;
                        self.deltas.push(d);
                    }
                    None => return Err(malformed("MM ends without a semicolon")),
                    Some(_) => return Err(malformed("invalid delta in MM")),
                }
            }
            multiple_groups |= i < mm.len();
            if self.types.len() + n_types >= MAX_TYPES {
                return Err(malformed("too many modification types"));
            }
            let mut flags = 0u8;
            if minus {
                flags |= ModificationEntry::MINUS;
            }
            if explicit {
                flags |= ModificationEntry::EXPLICIT;
            }
            if tags.ml.is_none() {
                flags |= ModificationEntry::UNSCORED;
            }
            self.selected.clear();
            for t in 0..n_types {
                let code = code_of(t);
                self.types.push(ModType {
                    code,
                    canonical: base,
                    minus,
                    explicit,
                });
                if codes.is_empty() || codes.contains(&code) {
                    self.selected.push((t, code));
                }
            }
            let n_ml = self
                .deltas
                .len()
                .checked_mul(n_types)
                .ok_or(malformed("too many ML values"))?;
            let ml_end = ml_at
                .checked_add(n_ml)
                .ok_or(malformed("too many ML values"))?;
            if let Some(ml) = tags.ml
                && ml_end > ml.len()
            {
                return Err(malformed("too few ML values"));
            }
            // A forward read without bases lists nothing (htslib stops at
            // the end of SEQ before its first base).
            // htslib tracks deltas through the entry's codes: an entry
            // without codes places nothing.
            if self.deltas.is_empty() || n_types == 0 || (l == 0 && !reverse) {
                ml_at = ml_end;
                continue;
            }
            let group_start = out.len();
            let mut emit = |di: usize, j: usize| {
                for &(t, code) in &self.selected {
                    // SAFETY: di indexes deltas, t is below n_types, and
                    // ml_end <= ml.len() was checked for the whole group.
                    let probability = tags.ml.map_or(0, |ml| unsafe {
                        *ml.get_unchecked(ml_at + di * n_types + t)
                    });
                    out.push(ModificationEntry {
                        ref_pos: -1,
                        qpos: j as i32,
                        code,
                        canonical: base,
                        probability,
                        flags,
                        _pad: [0; 5],
                    });
                }
            };
            if multiple_groups || canonical == 15 {
                let needed = self
                    .deltas
                    .iter()
                    .try_fold(0usize, |n, &d| n.checked_add(d as usize)?.checked_add(1))
                    .ok_or(malformed("MM lists more bases than the read has"))?;
                if needed > l {
                    return Err(malformed("MM lists more bases than the read has"));
                }
                let bucket = match canonical {
                    1 => 0,
                    2 => 1,
                    4 => 2,
                    _ => 3,
                };
                if canonical != 15 {
                    while self.positions[bucket].len() < needed && scanned < l {
                        let j = if reverse { l - 1 - scanned } else { scanned };
                        scanned += 1;
                        // SAFETY: j is in 0..l; packed lengths were validated.
                        let mut b = unsafe { seq.code_unchecked(j) };
                        if reverse {
                            // SAFETY: sequence conversion returns a 4-bit code.
                            b = unsafe { *NT16_RC.get_unchecked(b as usize) };
                        }
                        let observed = match b {
                            1 => 0,
                            2 => 1,
                            4 => 2,
                            8 => 3,
                            _ => continue,
                        };
                        self.positions[observed].push(j as u32);
                    }
                    if self.positions[bucket].len() < needed {
                        return Err(malformed("MM lists more bases than the read has"));
                    }
                }
                let positions = &self.positions[bucket];
                let mut next = 0usize;
                for (di, &delta) in self.deltas.iter().enumerate() {
                    // needed is the checked sum of every (delta + 1),
                    // so each occurrence is below needed and fits usize.
                    let occurrence = next + delta as usize;
                    let j = if canonical == 15 {
                        if reverse {
                            l - 1 - occurrence
                        } else {
                            occurrence
                        }
                    } else {
                        // SAFETY: occurrence < needed and the prefix was
                        // extended until positions.len() >= needed.
                        unsafe { *positions.get_unchecked(occurrence) as usize }
                    };
                    emit(di, j);
                    next = occurrence + 1;
                }
            } else {
                let mut di = 0usize;
                let mut skip = self.deltas[0];
                for k in 0..l {
                    let j = if reverse { l - 1 - k } else { k };
                    // SAFETY: j is in 0..l; packed lengths were validated.
                    let mut b = unsafe { seq.code_unchecked(j) };
                    if reverse {
                        // SAFETY: sequence conversion returns a 4-bit code.
                        b = unsafe { *NT16_RC.get_unchecked(b as usize) };
                    }
                    if canonical != 15 && b != canonical {
                        continue;
                    }
                    if skip > 0 {
                        skip -= 1;
                        continue;
                    }
                    emit(di, j);
                    di += 1;
                    if di == self.deltas.len() {
                        break;
                    }
                    skip = self.deltas[di];
                }
                if di < self.deltas.len() {
                    return Err(malformed("MM lists more bases than the read has"));
                }
            }
            if out.len() > group_start {
                hit_groups += 1;
                if reverse {
                    reverse_runs(&mut out[group_start..]);
                }
            }
            ml_at = ml_end;
        }
        if let Some(ml) = tags.ml
            && ml_at != ml.len()
        {
            return Err(malformed("too many ML values"));
        }
        if hit_groups > 1 {
            out[first..].sort_by_key(|e| e.qpos);
        }
        Ok(())
    }
}

/// Reverses the order of positions in rows that came out right to left,
/// keeping the codes of one position in `MM` order.
fn reverse_runs(rows: &mut [ModificationEntry]) {
    rows.reverse();
    let mut s = 0;
    while s < rows.len() {
        let mut e = s + 1;
        while e < rows.len() && rows[e].qpos == rows[s].qpos {
            e += 1;
        }
        rows[s..e].reverse();
        s = e;
    }
}

/// Fills `ref_pos` of `rows` (in SEQ order) for a record placed at `pos`
/// with `cigar` (`(op, len)`). Resets every `ref_pos` to -1 first; bases the
/// CIGAR does not align, and records with negative `pos`, retain -1.
/// Coordinate overflow stops placement; remaining rows retain -1.
pub fn reference_positions(
    cigar: impl Iterator<Item = (u8, u32)>,
    pos: i64,
    rows: &mut [ModificationEntry],
) {
    for row in rows.iter_mut() {
        row.ref_pos = -1;
    }
    if pos < 0 {
        return;
    }
    place_reference_positions(cigar, pos, rows);
}

/// Places rows whose reference positions are already initialized to -1.
fn place_reference_positions(
    cigar: impl Iterator<Item = (u8, u32)>,
    pos: i64,
    rows: &mut [ModificationEntry],
) {
    let mut q: i64 = 0;
    let mut r = pos;
    let mut at = 0;
    for (op, len) in cigar {
        if at == rows.len() {
            return;
        }
        let len = len as i64;
        let on_query = crate::record::consumes_query(op);
        let on_ref = crate::record::consumes_reference(op);
        let Some(q_end) = q.checked_add(if on_query { len } else { 0 }) else {
            return;
        };
        let Some(r_end) = r.checked_add(if on_ref { len } else { 0 }) else {
            return;
        };
        if on_query {
            while at < rows.len() {
                // SAFETY: at is checked against the slice length on every
                // iteration; the borrow ends before at is incremented.
                let row = unsafe { rows.get_unchecked_mut(at) };
                if row.qpos as i64 >= q_end {
                    break;
                }
                if on_ref && row.qpos as i64 >= q {
                    // The offset is nonnegative and below len, and r_end
                    // proved the entire reference-consuming span fits i64.
                    row.ref_pos = r + (row.qpos as i64 - q);
                }
                at += 1;
            }
        }
        q = q_end;
        r = r_end;
    }
}

/// [`ModDecoder::decode`] then [`reference_positions`] for a record, as a
/// region query keeps it; `Ok(0)` rows when it has no `MM` tag.
#[allow(clippy::too_many_arguments)]
pub fn decode_record(
    decoder: &mut ModDecoder,
    aux: &[&[u8]],
    seq: ModSeq<'_>,
    flag: u16,
    pos: i64,
    cigar: impl Iterator<Item = (u8, u32)>,
    name: &[u8],
    codes: &[i32],
    out: &mut Vec<ModificationEntry>,
) -> Result<usize> {
    decoder.types.clear();
    let first = out.len();
    let Some(tags) = ModTags::find(aux).map_err(|e| e.error(name))? else {
        return Ok(0);
    };
    decoder
        .decode(&tags, seq, flag & crate::flag::REVERSE != 0, codes, out)
        .map_err(|e| e.error(name))?;
    if flag & crate::flag::UNMAPPED == 0 && pos >= 0 {
        // decode initializes every appended row to -1.
        place_reference_positions(cigar, pos, &mut out[first..]);
    }
    Ok(out.len() - first)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn decode(mm: &str, ml: Option<&[u8]>, seq: &[u8], reverse: bool) -> Vec<(i32, i32, u8)> {
        let tags = ModTags {
            mm: mm.as_bytes(),
            ml,
            mn: None,
        };
        let mut out = Vec::new();
        ModDecoder::new()
            .decode(&tags, ModSeq::Ascii(seq), reverse, &[], &mut out)
            .unwrap();
        out.iter()
            .map(|e| (e.qpos, e.code, e.probability))
            .collect()
    }

    #[test]
    fn skips_count_canonical_bases() {
        // C at 1, 3, 6, 8: skip one, take the next, skip none, take.
        let got = decode("C+m,1,0;", Some(&[10, 20]), b"ACGCTTCAC", false);
        assert_eq!(got, vec![(3, b'm' as i32, 10), (6, b'm' as i32, 20)]);
    }

    #[test]
    fn ascii_aliases_count_the_same_canonical_bases_as_htslib() {
        assert_eq!(
            decode("T+m,0,0,0;", None, b"0123Uu", false),
            vec![
                (3, b'm' as i32, 0),
                (4, b'm' as i32, 0),
                (5, b'm' as i32, 0)
            ]
        );
        assert_eq!(
            decode("A+a,0;C+m,0;G+h,0;T+t,0,0,0;", None, b"0123Uu", false),
            vec![
                (0, b'a' as i32, 0),
                (1, b'm' as i32, 0),
                (2, b'h' as i32, 0),
                (3, b't' as i32, 0),
                (4, b't' as i32, 0),
                (5, b't' as i32, 0)
            ]
        );
        assert_eq!(
            decode("A+a,0,0;", Some(&[1, 2]), b"Uu", true),
            vec![(0, b'a' as i32, 2), (1, b'a' as i32, 1)]
        );
    }

    #[test]
    fn reference_placement_reinitializes_unaligned_rows() {
        let mut rows: Vec<_> = (0..4)
            .map(|qpos| ModificationEntry {
                qpos,
                ref_pos: 999,
                ..Default::default()
            })
            .collect();
        reference_positions([(4, 1), (0, 2), (1, 1)].into_iter(), 10, &mut rows);
        assert_eq!(
            rows.iter().map(|r| r.ref_pos).collect::<Vec<_>>(),
            [-1, 10, 11, -1]
        );
        reference_positions([(0, 1)].into_iter(), 20, &mut rows);
        assert_eq!(
            rows.iter().map(|r| r.ref_pos).collect::<Vec<_>>(),
            [20, -1, -1, -1]
        );
        reference_positions([(0, 4)].into_iter(), -1, &mut rows);
        assert!(rows.iter().all(|r| r.ref_pos == -1));
    }

    #[test]
    fn failed_or_untagged_records_leave_no_modification_types() {
        let mut decoder = ModDecoder::new();
        let mut rows = Vec::new();
        let tags = ModTags {
            mm: b"C+m,0;",
            ml: None,
            mn: None,
        };
        decoder
            .decode(&tags, ModSeq::Ascii(b"C"), false, &[], &mut rows)
            .unwrap();
        assert_eq!(decoder.types().len(), 1);
        assert_eq!(
            decode_record(
                &mut decoder,
                &[],
                ModSeq::Ascii(b"C"),
                0,
                0,
                [(0, 1)].into_iter(),
                b"untagged",
                &[],
                &mut rows
            )
            .unwrap(),
            0
        );
        assert!(decoder.types().is_empty());
        assert_eq!(rows.len(), 1);
        let bad = ModTags {
            mm: b"C+m,0;A+a,9;",
            ml: None,
            mn: None,
        };
        assert!(
            decoder
                .decode(&bad, ModSeq::Ascii(b"CA"), false, &[], &mut rows)
                .is_err()
        );
        assert!(decoder.types().is_empty());
        assert_eq!(rows.len(), 1);
    }

    #[test]
    fn reverse_reads_count_from_the_right_on_the_complement() {
        // Stored GTTG...: as sequenced CAAC (C at stored 3 and 0).
        let got = decode("C+m,0,0;", Some(&[1, 2]), b"GTTG", true);
        assert_eq!(got, vec![(0, b'm' as i32, 2), (3, b'm' as i32, 1)]);
    }

    #[test]
    fn interleaved_codes_and_chebi() {
        let got = decode("C+mh,0;A+17596,1;", Some(&[5, 6, 7]), b"CAA", false);
        assert_eq!(
            got,
            vec![(0, b'm' as i32, 5), (0, b'h' as i32, 6), (2, -17596, 7)]
        );
    }

    #[test]
    fn indexed_groups_preserve_position_group_and_code_order() {
        let mut decoder = ModDecoder::new();
        // Two C groups share their canonical index; N counts ambiguous
        // bases too. Reusing scratch for a reverse read must rebuild it.
        let tags = ModTags {
            mm: b"C+mh,0,1;N+n,1;C+f,0,0;",
            ml: Some(&[10, 11, 20, 21, 30, 40, 50]),
            mn: None,
        };
        for (seq, reverse, expected) in [
            (
                &b"CNCAC"[..],
                false,
                vec![
                    (0, b'm', 10),
                    (0, b'h', 11),
                    (0, b'f', 40),
                    (1, b'n', 30),
                    (2, b'f', 50),
                    (4, b'm', 20),
                    (4, b'h', 21),
                ],
            ),
            (
                &b"GTGNG"[..],
                true,
                vec![
                    (0, b'm', 20),
                    (0, b'h', 21),
                    (2, b'f', 50),
                    (3, b'n', 30),
                    (4, b'm', 10),
                    (4, b'h', 11),
                    (4, b'f', 40),
                ],
            ),
        ] {
            let mut out = Vec::new();
            decoder
                .decode(&tags, ModSeq::Ascii(seq), reverse, &[], &mut out)
                .unwrap();
            let got: Vec<_> = out
                .iter()
                .map(|m| (m.qpos, m.code as u8, m.probability))
                .collect();
            assert_eq!(got, expected);
            decoder
                .decode(&tags, ModSeq::Ascii(seq), reverse, &[b'h' as i32], &mut out)
                .unwrap();
            assert_eq!(
                out[out.len() - 2..]
                    .iter()
                    .map(|m| m.probability)
                    .collect::<Vec<_>>(),
                if reverse { vec![21, 11] } else { vec![11, 21] }
            );
        }
        // A later canonical type may extend the prefix; an earlier type
        // still resolves against the positions retained before that scan.
        assert_eq!(
            decode("A+a,1;C+m,0;T+t,0;C+h,2;", None, b"CCAACT", false),
            vec![
                (0, b'm' as i32, 0),
                (3, b'a' as i32, 0),
                (4, b'h' as i32, 0),
                (5, b't' as i32, 0)
            ]
        );
    }

    #[test]
    fn packed_lengths_are_checked_before_sequence_access_and_errors_roll_back() {
        let tags = ModTags {
            mm: b"C+m,0;",
            ml: None,
            mn: None,
        };
        let mut decoder = ModDecoder::new();
        let sentinel = ModificationEntry {
            code: 999,
            ..Default::default()
        };
        let mut out = vec![sentinel];
        for seq in [
            ModSeq::Packed { bytes: &[], len: 1 },
            ModSeq::Packed {
                bytes: &[0x22],
                len: 3,
            },
            ModSeq::Packed {
                bytes: &[],
                len: i32::MAX as usize + 1,
            },
        ] {
            assert!(decoder.decode(&tags, seq, false, &[], &mut out).is_err());
            assert_eq!(out, [sentinel]);
        }
        let bad = ModTags {
            mm: b"C+m,0;C+h,3;",
            ml: None,
            mn: None,
        };
        assert!(
            decoder
                .decode(&bad, ModSeq::Ascii(b"CC"), false, &[], &mut out)
                .is_err()
        );
        assert_eq!(out, [sentinel]);
        // An odd packed length never examines the trailing low nibble.
        decoder
            .decode(
                &tags,
                ModSeq::Packed {
                    bytes: &[0x21],
                    len: 1,
                },
                false,
                &[],
                &mut out,
            )
            .unwrap();
        assert_eq!(out.len(), 2);
        assert_eq!(out[1].qpos, 0);
    }

    #[test]
    fn reference_placement_stops_at_coordinate_overflow() {
        let mut rows = [ModificationEntry {
            qpos: 0,
            ref_pos: -1,
            ..Default::default()
        }];
        reference_positions([(0, 2)].into_iter(), i64::MAX, &mut rows);
        assert_eq!(rows[0].ref_pos, -1);
        reference_positions([(2, 1), (0, 1)].into_iter(), i64::MAX - 2, &mut rows);
        assert_eq!(rows[0].ref_pos, i64::MAX - 1);
    }

    #[test]
    fn errors() {
        let bad = |mm: &str, ml: Option<&[u8]>| {
            let tags = ModTags {
                mm: mm.as_bytes(),
                ml,
                mn: None,
            };
            ModDecoder::new()
                .decode(&tags, ModSeq::Ascii(b"CCA"), false, &[], &mut Vec::new())
                .is_err()
        };
        assert!(bad("C+m,2;", None));
        assert!(bad("C+m,0", None));
        assert!(bad("X+m,0;", None));
        assert!(bad("C+m,0;", Some(&[1, 2])));
        assert!(bad("C+m,0,0;", Some(&[1])));
        assert!(!bad("C+m,0,0;", Some(&[1, 2])));
    }
}
