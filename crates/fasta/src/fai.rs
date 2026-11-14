//! The FASTA index (`.fai`): where each sequence starts in the FASTA text and
//! how its lines are laid out, which makes any base range addressable.
//!
//! Offsets are positions in the uncompressed text, for plain and bgzip FASTA
//! alike; for bgzip files, [`crate::Gzi`] maps them to compressed blocks.

use std::collections::HashMap;
use std::fmt::Write as _;
use std::ops::Range;
use std::path::{Path, PathBuf};

use omicsweave_core::{FileSource, ReadAt, sidecar};

use crate::error::{Error, Result, invalid};
use crate::index::build_index;

/// Appended to a FASTA file's path to name its `.fai`.
pub const FAI_SUFFIX: &str = ".fai";

/// The `.fai` path of a FASTA file, as samtools names it (`a.fa` ->
/// `a.fa.fai`, `a.fa.gz` -> `a.fa.gz.fai`).
pub fn fai_path(fasta: impl AsRef<Path>) -> PathBuf {
    let mut name = fasta.as_ref().as_os_str().to_owned();
    name.push(FAI_SUFFIX);
    PathBuf::from(name)
}

#[inline]
fn is_layout_space(c: &u8) -> bool {
    // FASTA line layout follows the six whitespace bytes of the C locale.
    matches!(*c, b' ' | b'\t' | b'\n' | 0x0B | 0x0C | b'\r')
}

/// One sequence of a FASTA file.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct FaiEntry {
    pub name: String,
    /// Bases in the sequence.
    pub len: u32,
    /// Text offset of the first base.
    pub offset: u64,
    /// Bases per full line.
    pub line_bases: u32,
    /// Bytes per full line, terminator included.
    pub line_width: u32,
}

impl FaiEntry {
    fn check(&self, start: u32, end: u32) -> Result<()> {
        if start > end {
            return Err(Error::InvalidRegion { start, end });
        }
        if end > self.len {
            return Err(Error::RegionOutOfRange {
                start,
                end,
                len: self.len,
            });
        }
        if end > start && self.line_bases == 0 {
            return Err(invalid(
                "fai",
                format!("{}: zero bases per line", self.name),
            ));
        }
        if self.line_width < self.line_bases {
            return Err(invalid("fai", "line width below line bases"));
        }
        Ok(())
    }

    /// Text offset of base `pos`.
    #[inline]
    fn position(&self, pos: u32) -> Result<u64> {
        let (pos, bases, width) = (pos as u64, self.line_bases as u64, self.line_width as u64);
        self.offset
            .checked_add(pos / bases * width + pos % bases)
            .ok_or_else(|| invalid("fai", "text offset overflows u64"))
    }

    /// Text bytes holding bases `[start, end)`, line terminators between
    /// them included. Read these and pass them to [`FaiEntry::decode`].
    pub fn text_range(&self, start: u32, end: u32) -> Result<Range<u64>> {
        self.check(start, end)?;
        if start == end {
            return Ok(self.offset..self.offset);
        }
        let last = self
            .position(end - 1)?
            .checked_add(1)
            .ok_or_else(|| invalid("fai", "text range overflows u64"))?;
        Ok(self.position(start)?..last)
    }

    /// Writes bases `[start, start + out.len())` into `out` from `raw`, the
    /// bytes of [`FaiEntry::text_range`] (longer buffers are fine). Bases
    /// are copied as stored, case included.
    pub fn decode(&self, start: u32, raw: &[u8], out: &mut [u8]) -> Result<()> {
        let end = u32::try_from(start as u64 + out.len() as u64).map_err(|_| {
            Error::RegionOutOfRange {
                start,
                end: u32::MAX,
                len: self.len,
            }
        })?;
        let range = self.text_range(start, end)?;
        let need = range.end - range.start;
        if (raw.len() as u64) < need {
            return Err(Error::Truncated {
                what: "FASTA text",
                need,
                have: raw.len() as u64,
            });
        }
        if out.is_empty() {
            return Ok(());
        }
        let bases = self.line_bases as usize;
        let gap = (self.line_width - self.line_bases) as usize;
        let first = (bases - (start % self.line_bases) as usize).min(out.len());
        let rest = out.len() - first;
        if rest <= bases {
            if rest > 0 && !raw[first..first + gap].iter().all(is_layout_space) {
                return Err(invalid("FASTA text", "line layout differs from the index"));
            }
            // SAFETY: text_range checked the metadata and coordinates. A
            // window spanning at most two lines needs first+rest base bytes
            // plus one gap iff rest>0, all covered by the raw length check.
            // first <= out.len(); the two exclusive destinations are disjoint.
            unsafe {
                out.get_unchecked_mut(..first)
                    .copy_from_slice(raw.get_unchecked(..first));
                if rest > 0 {
                    out.get_unchecked_mut(first..)
                        .copy_from_slice(raw.get_unchecked(first + gap..first + gap + rest));
                }
            }
            return Ok(());
        }
        let mut decoder = TextDecoder::new(self, start, end)?;
        decoder.push(raw, out)?;
        decoder.finish(out.len())
    }
}

/// Line-layout state across BGZF boundaries, including a split CRLF gap.
pub(crate) struct TextDecoder<'a> {
    entry: &'a FaiEntry,
    start: u32,
    column: usize,
    pending_gap: usize,
    written: usize,
}

impl<'a> TextDecoder<'a> {
    pub(crate) fn new(entry: &'a FaiEntry, start: u32, end: u32) -> Result<Self> {
        entry.check(start, end)?;
        Ok(Self {
            entry,
            start,
            column: if entry.line_bases == 0 {
                0
            } else {
                (start % entry.line_bases) as usize
            },
            pending_gap: 0,
            written: 0,
        })
    }

    pub(crate) fn written(&self) -> usize {
        self.written
    }

    pub(crate) fn advance_contiguous(&mut self, n: usize) {
        debug_assert_eq!(self.entry.line_width, self.entry.line_bases);
        self.written += n;
    }

    pub(crate) fn push(&mut self, mut raw: &[u8], out: &mut [u8]) -> Result<()> {
        debug_assert!(self.written <= out.len());
        let gap = (self.entry.line_width - self.entry.line_bases) as usize;
        if gap == 0 {
            let n = raw.len().min(out.len() - self.written);
            out[self.written..self.written + n].copy_from_slice(&raw[..n]);
            self.written += n;
            return Ok(());
        }
        while !raw.is_empty() && self.written < out.len() {
            if self.pending_gap > 0 {
                let n = self.pending_gap.min(raw.len());
                if !raw[..n].iter().all(is_layout_space) {
                    return Err(invalid(
                        "FASTA text",
                        format!(
                            "{}: line layout differs from the index near base {}",
                            self.entry.name,
                            self.start as u64 + self.written as u64
                        ),
                    ));
                }
                self.pending_gap -= n;
                raw = &raw[n..];
                if self.pending_gap == 0 {
                    self.column = 0;
                }
                if self.pending_gap != 0 {
                    break;
                }
            }
            let bases = self.entry.line_bases as usize;
            let width = self.entry.line_width as usize;
            if self.column == 0 && bases > 0 {
                while raw.len() >= width && out.len() - self.written > bases {
                    if !raw[bases..width].iter().all(is_layout_space) {
                        return Err(invalid("FASTA text", "line layout differs from the index"));
                    }
                    // SAFETY: width >= bases and raw contains a complete
                    // line. The remaining output is strictly larger than
                    // bases, so both ranges are in bounds. Input and output
                    // cannot alias through these shared/exclusive borrows.
                    unsafe {
                        out.get_unchecked_mut(self.written..self.written + bases)
                            .copy_from_slice(raw.get_unchecked(..bases));
                    }
                    self.written += bases;
                    raw = &raw[width..];
                }
            }
            if raw.is_empty() {
                break;
            }
            let n = (self.entry.line_bases as usize - self.column)
                .min(raw.len())
                .min(out.len() - self.written);
            // SAFETY: n is bounded by both the remaining input and output.
            // The layout was checked at construction; column stays <= bases.
            unsafe { out.get_unchecked_mut(self.written..self.written + n) }
                .copy_from_slice(&raw[..n]);
            self.written += n;
            self.column += n;
            raw = &raw[n..];
            if self.column == self.entry.line_bases as usize {
                self.pending_gap = gap;
            }
        }
        Ok(())
    }

    pub(crate) fn finish(&self, need: usize) -> Result<()> {
        if self.written != need {
            return Err(Error::Truncated {
                what: "FASTA bases",
                need: need as u64,
                have: self.written as u64,
            });
        }
        Ok(())
    }
}

/// A parsed `.fai`. Index metadata is small and kept by the value.
#[derive(Debug, Clone, Default, PartialEq, Eq)]
pub struct Fai {
    entries: Vec<FaiEntry>,
    by_name: HashMap<String, usize>,
}

impl Fai {
    /// Builds an index from entries; a repeated name keeps its first entry,
    /// as htslib does.
    pub fn from_entries(entries: impl IntoIterator<Item = FaiEntry>) -> Fai {
        let mut fai = Fai::default();
        for entry in entries {
            fai.insert(entry);
        }
        fai
    }

    fn insert(&mut self, entry: FaiEntry) {
        if !self.by_name.contains_key(&entry.name) {
            self.insert_missing(entry);
        }
    }

    /// Callers check the name immediately before insertion, with no
    /// intervening mutation of this index.
    pub(crate) fn insert_missing(&mut self, entry: FaiEntry) {
        let previous = self.by_name.insert(entry.name.clone(), self.entries.len());
        debug_assert!(previous.is_none());
        self.entries.push(entry);
    }

    /// Parses `.fai` text: one `name, length, offset, line bases, line width`
    /// row per sequence, tab-separated.
    pub fn parse(text: &[u8]) -> Result<Fai> {
        let text = std::str::from_utf8(text).map_err(|e| invalid("fai", e.to_string()))?;
        let mut fai = Fai::default();
        for (i, line) in text.lines().enumerate() {
            if line.is_empty() {
                continue;
            }
            let mut split = line.split('\t');
            let fields: [Option<&str>; 5] = std::array::from_fn(|_| split.next());
            let extra = split.count();
            if extra == 1 {
                return Err(invalid("fai", "FASTQ indexes are not supported"));
            }
            if extra > 0 || fields.iter().any(Option::is_none) {
                return Err(invalid("fai", format!("row {}: expected 5 fields", i + 1)));
            }
            let fields = fields.map(|field| field.expect("five fields checked"));
            let bad =
                |k: usize| invalid("fai", format!("row {}: bad number {:?}", i + 1, fields[k]));
            let small = |k: usize| -> Result<u32> { fields[k].parse().map_err(|_| bad(k)) };
            let entry = FaiEntry {
                name: fields[0].to_owned(),
                len: small(1)?,
                offset: fields[2].parse().map_err(|_| bad(2))?,
                line_bases: small(3)?,
                line_width: small(4)?,
            };
            if entry.line_width < entry.line_bases {
                return Err(invalid(
                    "fai",
                    format!("row {}: line width below line bases", i + 1),
                ));
            }
            fai.insert(entry);
        }
        Ok(fai)
    }

    /// The `.fai` text, as htslib writes it.
    pub fn to_text(&self) -> String {
        let mut out = String::new();
        for e in &self.entries {
            let _ = writeln!(
                out,
                "{}\t{}\t{}\t{}\t{}",
                e.name, e.len, e.offset, e.line_bases, e.line_width
            );
        }
        out
    }

    /// The `.fai` file bytes ([`Fai::to_text`]).
    pub fn to_bytes(&self) -> Vec<u8> {
        self.to_text().into_bytes()
    }

    /// Builds the index of a plain or BGZF FASTA source, streaming it in
    /// bounded reads; see [`crate::build_index`], which also returns the
    /// `.gzi` of BGZF sources.
    pub fn build<S: ReadAt + ?Sized>(source: &S) -> Result<Fai> {
        Ok(build_index(source)?.0)
    }

    /// Writes [`fai_path`]`(fasta)` atomically.
    pub fn save(&self, fasta: impl AsRef<Path>) -> Result<()> {
        sidecar::write_file(&fai_path(fasta), &self.to_bytes())
    }

    /// Reads and parses [`fai_path`]`(fasta)`.
    pub fn load(fasta: impl AsRef<Path>) -> Result<Fai> {
        Fai::parse(&std::fs::read(fai_path(fasta))?)
    }

    /// Builds the index of the FASTA file at `fasta` and saves it beside it.
    pub fn build_file(fasta: impl AsRef<Path>) -> Result<Fai> {
        let fasta = fasta.as_ref();
        let fai = Fai::build(&FileSource::open(fasta)?)?;
        fai.save(fasta)?;
        Ok(fai)
    }

    pub fn get(&self, name: &str) -> Option<&FaiEntry> {
        self.by_name.get(name).map(|&i| &self.entries[i])
    }

    /// Sequences in file order.
    pub fn entries(&self) -> &[FaiEntry] {
        &self.entries
    }

    pub fn len(&self) -> usize {
        self.entries.len()
    }

    pub fn is_empty(&self) -> bool {
        self.entries.is_empty()
    }
}

#[cfg(test)]
mod text_tests {
    use super::*;

    #[test]
    fn streaming_lines_match_bases_at_every_chunk_boundary() {
        let bases: Vec<_> = (0..512).map(|i| b"ACGTNacgtn"[i % 10]).collect();
        for width in [1usize, 2, 60] {
            for gap in [b"\n".as_slice(), b"\r\n", b" \t \t \t \t \t \t \t \t\n"] {
                let mut text = Vec::new();
                for line in bases.chunks(width) {
                    text.extend(line);
                    text.extend(gap);
                }
                let entry = FaiEntry {
                    name: "s".into(),
                    len: 512,
                    offset: 0,
                    line_bases: width as u32,
                    line_width: (width + gap.len()) as u32,
                };
                for start in 0..2 * width {
                    for len in [0, 1, width - 1, width, width + 1, 2 * width + 1] {
                        let end = start + len;
                        let range = entry.text_range(start as u32, end as u32).unwrap();
                        let raw = &text[range.start as usize..range.end as usize];
                        let mut decoded = vec![0xAA; len + 2];
                        entry
                            .decode(start as u32, raw, &mut decoded[1..len + 1])
                            .unwrap();
                        assert_eq!(&decoded[1..len + 1], &bases[start..end]);
                        assert_eq!((decoded[0], decoded[len + 1]), (0xAA, 0xAA));
                        for chunk_len in [1, 2, 3, 31, 60, 61, 62, 128] {
                            let mut out = vec![0xAA; len + 2];
                            let mut decoder =
                                TextDecoder::new(&entry, start as u32, end as u32).unwrap();
                            for chunk in raw.chunks(chunk_len) {
                                decoder.push(chunk, &mut out[1..len + 1]).unwrap();
                            }
                            decoder.finish(len).unwrap();
                            assert_eq!(&out[1..len + 1], &bases[start..end]);
                            assert_eq!((out[0], out[len + 1]), (0xAA, 0xAA));
                        }
                    }
                }
            }
        }
    }
}
