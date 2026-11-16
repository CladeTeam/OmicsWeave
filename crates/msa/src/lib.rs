//! bigMaf alignments projected onto forward reference coordinates.
//!
//! [`MsaReader::query`] decodes once into reusable staging storage and
//! returns its shape. [`MsaReader::write`] copies those rows into independent
//! caller storage, orienting each row on demand. Coordinates beyond the
//! chromosome end retain the requested length and are filled.

use std::collections::HashMap;
use std::sync::Arc;

use omicsweave_bbi::{BbiReader, Chroms, Kind};
use omicsweave_core::cache::BlockCache;
use omicsweave_core::orient::{complement, reverse_complement};
use omicsweave_core::{
    Error, FileSource, IntervalRead, ReadAt, Region, Result, Strand, check_output,
};

pub const NOT_BIGMAF: &str = "not a bigMaf: a bigMaf is a bigBed of bed3 and a mafBlock field";

/// Row-major output shape. Names are available through [`MsaReader::species`].
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct MsaSize {
    pub species: usize,
    pub len: usize,
}

impl MsaSize {
    pub fn elements(self) -> Result<usize> {
        self.species
            .checked_mul(self.len)
            .filter(|&size| size <= isize::MAX as usize)
            .ok_or_else(|| Error::invalid("MSA shape", "exceeds addressable output size"))
    }
}

/// A dedicated bigMaf reader; caches retain the underlying decoded bigBed blocks.
pub struct MsaReader<S> {
    bbi: BbiReader<S>,
    fixed: bool,
    names: Vec<String>,
    rows: HashMap<Vec<u8>, usize>,
    fill: u8,
    bases: Vec<u8>,
    /// Provenance is needed only when the chosen fill has a complement.
    covered: Vec<bool>,
    keep: Vec<usize>,
    seen: Vec<bool>,
    size: Option<MsaSize>,
}

impl<S: ReadAt> MsaReader<S> {
    pub fn open(source: S, species: Option<Vec<String>>, fill: u8) -> Result<Self> {
        Self::from_bbi(BbiReader::open(source)?, species, fill)
    }

    fn from_bbi(bbi: BbiReader<S>, species: Option<Vec<String>>, fill: u8) -> Result<Self> {
        if bbi.header().kind != Kind::BigBed
            || bbi.header().field_count != 4
            || !bbi.auto_sql()?.is_some_and(|s| maf_schema(&s))
        {
            return Err(Error::invalid("bigMaf", NOT_BIGMAF));
        }
        let fixed = species.is_some();
        let names = species.unwrap_or_default();
        let mut rows = HashMap::new();
        for (i, name) in names.iter().enumerate() {
            if rows.insert(name.as_bytes().to_vec(), i).is_some() {
                return Err(Error::invalid(
                    "bigMaf species",
                    "species are named more than once",
                ));
            }
        }
        Ok(Self {
            bbi,
            fixed,
            names,
            rows,
            fill,
            bases: Vec::new(),
            covered: Vec::new(),
            keep: Vec::new(),
            seen: Vec::new(),
            size: None,
        })
    }

    pub fn chroms(&self) -> &Chroms {
        self.bbi.chroms()
    }

    pub fn species(&self) -> &[String] {
        &self.names
    }

    pub fn set_cache(&mut self, cache: Option<Arc<BlockCache>>) -> Result<()> {
        self.bbi.set_cache(cache)
    }

    pub fn query(&mut self, region: Region<'_>) -> Result<MsaSize> {
        self.size = None;
        if !self.fixed {
            self.names.clear();
            self.rows.clear();
        }
        let len = region
            .end
            .checked_sub(region.start)
            .ok_or(Error::InvalidRegion {
                start: region.start,
                end: region.end,
            })? as usize;
        self.bases.clear();
        self.covered.clear();
        let provenance = complement(self.fill) != self.fill;
        let initial = MsaSize {
            species: self.names.len(),
            len,
        }
        .elements()?;
        resize(&mut self.bases, initial, self.fill)?;
        if provenance {
            resize(&mut self.covered, initial, false)?;
        }
        self.bbi.query(region, false)?;
        let chrom_len = self
            .bbi
            .chroms()
            .id(region.chrom)
            .and_then(|id| self.bbi.chroms().get(id))
            .map(|c| c.len);
        let Self {
            bbi,
            names,
            rows,
            bases,
            covered,
            keep,
            seen,
            fill,
            fixed,
            ..
        } = self;
        bbi.visit_intervals(|start, end, block| {
            let lo = start.max(region.start);
            let hi = end.min(region.end);
            if lo >= hi {
                return Ok(());
            }
            let error = |detail: String| {
                Error::invalid(
                    "bigMaf",
                    format!("{}:{start}-{end}: {detail}", region.chrom),
                )
            };
            // All slices borrow the retained decompressed block. Scratch is
            // reused between records, including the reference column mapping.
            let reference = lines(block)
                .find(|line| is_sequence(line))
                .ok_or_else(|| error("a block without s lines".into()))?;
            let reference_row = sequence(reference).map_err(&error)?;
            let reference = reference_row.text;
            keep.clear();
            keep.extend(
                reference
                    .iter()
                    .enumerate()
                    .filter_map(|(i, &b)| (b != b'-').then_some(i)),
            );
            if keep.len() != (end - start) as usize {
                return Err(error(format!(
                    "the reference row has {} bases, not {}",
                    keep.len(),
                    end - start
                )));
            }
            if reference_row.strand != b'+' {
                return Err(error("the reference row is not on the + strand".into()));
            }
            if reference_row.start != start as u64 || reference_row.size != (end - start) as u64 {
                return Err(error(
                    "reference coordinates differ from the bigBed interval".into(),
                ));
            }
            if chrom_len != u32::try_from(reference_row.src_size).ok() {
                return Err(error(
                    "reference srcSize differs from the bigBed chromosome".into(),
                ));
            }
            let src = reference_row.source.as_bytes();
            let chrom = region.chrom.as_bytes();
            if src != chrom
                && !(src.len() > chrom.len()
                    && src.ends_with(chrom)
                    && src[src.len() - chrom.len() - 1] == b'.')
            {
                return Err(error(
                    "reference source differs from the bigBed chromosome".into(),
                ));
            }
            let columns = &keep[(lo - start) as usize..(hi - start) as usize];
            let first_column = columns[0];
            // keep comes from enumerate() and is strictly increasing. Its
            // endpoints therefore identify a contiguous run in O(1).
            let contiguous = columns[columns.len() - 1] - first_column + 1 == columns.len();
            seen.clear();
            seen.resize(names.len(), false);
            for (i, line) in lines(block).filter(|line| is_sequence(line)).enumerate() {
                let parsed = if i == 0 {
                    reference_row
                } else {
                    sequence(line).map_err(&error)?
                };
                let text = parsed.text;
                if text.len() != reference.len() {
                    return Err(error(format!(
                        "a row has {} columns, not {}",
                        text.len(),
                        reference.len()
                    )));
                }
                if parsed.bases as u64 != parsed.size {
                    return Err(error(format!(
                        "row {} has {} bases, not {}",
                        parsed.source, parsed.bases, parsed.size
                    )));
                }
                let name = parsed.source.split('.').next().unwrap().as_bytes();
                if name.is_empty() {
                    return Err(error("an empty species name".into()));
                }
                let row = if let Some(&row) = rows.get(name) {
                    row
                } else {
                    if *fixed {
                        continue;
                    }
                    let name_text = std::str::from_utf8(name).map_err(|e| error(e.to_string()))?;
                    let row = names.len();
                    names.push(name_text.to_owned());
                    rows.insert(name.to_vec(), row);
                    let size = MsaSize {
                        species: names.len(),
                        len,
                    }
                    .elements()?;
                    resize(bases, size, *fill)?;
                    if provenance {
                        resize(covered, size, false)?;
                    }
                    seen.push(false);
                    row
                };
                if seen[row] {
                    continue;
                }
                seen[row] = true;
                let offset = row * len + (lo - region.start) as usize;
                let output = &mut bases[offset..offset + columns.len()];
                if contiguous {
                    output.copy_from_slice(&text[first_column..first_column + columns.len()]);
                } else {
                    for (out, &column) in output.iter_mut().zip(columns) {
                        // SAFETY: keep contains only indices from reference's
                        // enumerate(), and this row's text length was checked to
                        // equal reference.len(). The decoded block (or cache lease)
                        // remains borrowed by this visitor for the whole copy.
                        *out = unsafe { *text.get_unchecked(column) };
                    }
                }
                if provenance {
                    covered[offset..offset + columns.len()].fill(true);
                }
            }
            Ok(())
        })?;
        let size = MsaSize {
            species: names.len(),
            len,
        };
        self.size = Some(size);
        Ok(size)
    }

    pub fn write(&self, strand: Strand, out: &mut [u8]) -> Result<()> {
        let size = self.size.ok_or(Error::NotQueried)?;
        check_output("MSA output", out.len(), size.elements()?)?;
        out.copy_from_slice(&self.bases);
        if strand.is_reverse() && size.len != 0 {
            if self.covered.is_empty() {
                for row in out.chunks_exact_mut(size.len) {
                    reverse_complement(row);
                }
            } else {
                for (row, covered) in out
                    .chunks_exact_mut(size.len)
                    .zip(self.covered.chunks_exact(size.len))
                {
                    row.reverse();
                    for (base, &known) in row.iter_mut().zip(covered.iter().rev()) {
                        if known {
                            *base = complement(*base);
                        }
                    }
                }
            }
        }
        Ok(())
    }
}

fn lines(mut bytes: &[u8]) -> impl Iterator<Item = &[u8]> {
    std::iter::from_fn(move || {
        if bytes.is_empty() {
            return None;
        }
        let stop = memchr::memchr(b';', bytes).unwrap_or(bytes.len());
        let line = bytes[..stop].trim_ascii();
        bytes = &bytes[(stop + 1).min(bytes.len())..];
        Some(line)
    })
}

fn resize<T: Clone>(values: &mut Vec<T>, len: usize, fill: T) -> Result<()> {
    values
        .try_reserve(len.saturating_sub(values.len()))
        .map_err(|e| Error::invalid("MSA output", e.to_string()))?;
    values.resize(len, fill);
    Ok(())
}

fn is_sequence(line: &[u8]) -> bool {
    line.first() == Some(&b's') && line.get(1).is_some_and(u8::is_ascii_whitespace)
}

#[derive(Clone, Copy)]
struct SequenceRow<'a> {
    source: &'a str,
    start: u64,
    size: u64,
    strand: u8,
    src_size: u64,
    text: &'a [u8],
    bases: usize,
}

fn sequence(line: &[u8]) -> std::result::Result<SequenceRow<'_>, String> {
    let mut fields: [&[u8]; 6] = [&[]; 6];
    let mut rest = line;
    for (index, field) in fields.iter_mut().enumerate() {
        rest = rest.trim_ascii_start();
        if rest.is_empty() {
            return Err(format!("an s line of {index} words, not 7"));
        }
        let end = rest
            .iter()
            .position(u8::is_ascii_whitespace)
            .unwrap_or(rest.len());
        *field = &rest[..end];
        rest = &rest[end..];
    }
    let text = rest.trim_ascii();
    if text.is_empty() {
        return Err("an s line of 6 words, not 7".into());
    }
    // Only the short metadata fields are tokenized. The large text payload
    // is checked and counted in one pass, without UTF-8 conversion.
    let mut bases = 0usize;
    let mut invalid_text = false;
    for &base in text {
        invalid_text |= !base.is_ascii_graphic();
        bases += usize::from(base != b'-');
    }
    if invalid_text {
        // Whitespace inside the payload denotes extra fields. Preserve the
        // malformed-row diagnostic without rescanning valid alignments.
        let count = 6 + text
            .split(u8::is_ascii_whitespace)
            .filter(|w| !w.is_empty())
            .count();
        if count != 7 {
            return Err(format!("an s line of {count} words, not 7"));
        }
    }
    let source = std::str::from_utf8(fields[1]).map_err(|e| format!("source is not UTF-8: {e}"))?;
    let number = |field: &[u8], what: &str| -> std::result::Result<u64, String> {
        std::str::from_utf8(field)
            .ok()
            .and_then(|s| s.parse().ok())
            .ok_or_else(|| format!("invalid {what}"))
    };
    let start = number(fields[2], "start")?;
    let size = number(fields[3], "size")?;
    let src_size = number(fields[5], "srcSize")?;
    if start.checked_add(size).is_none_or(|end| end > src_size) {
        return Err("an s row extends beyond srcSize".into());
    }
    let strand = match fields[4] {
        b"+" => b'+',
        b"-" => b'-',
        _ => return Err("invalid strand".into()),
    };
    if invalid_text {
        return Err("alignment text is not printable ASCII".into());
    }
    Ok(SequenceRow {
        source,
        start,
        size,
        strand,
        src_size,
        text,
        bases,
    })
}

fn maf_schema(sql: &str) -> bool {
    let mut unquoted = String::with_capacity(sql.len());
    let (mut quoted, mut escaped) = (false, false);
    for c in sql.chars() {
        if quoted {
            if escaped {
                escaped = false;
            } else if c == '\\' {
                escaped = true;
            } else if c == '"' {
                quoted = false;
            }
        } else if c == '"' {
            quoted = true;
        } else {
            unquoted.push(c);
        }
    }
    if quoted {
        return false;
    }
    let Some((_, body)) = unquoted.split_once('(') else {
        return false;
    };
    let Some((body, _)) = body.split_once(')') else {
        return false;
    };
    let mut fields = body.split(';').map(str::trim).filter(|s| !s.is_empty());
    for _ in 0..3 {
        if fields.next().is_none() {
            return false;
        }
    }
    let Some(field) = fields.next() else {
        return false;
    };
    let mut words = field.split_whitespace();
    matches!(words.next(), Some("string" | "lstring"))
        && words.next() == Some("mafBlock")
        && words.next().is_none()
        && fields.next().is_none()
}

impl MsaReader<FileSource> {
    pub fn try_clone(&self) -> Result<Self> {
        Ok(Self {
            bbi: self.bbi.try_clone()?,
            fixed: self.fixed,
            names: if self.fixed {
                self.names.clone()
            } else {
                Vec::new()
            },
            rows: if self.fixed {
                self.rows.clone()
            } else {
                HashMap::new()
            },
            fill: self.fill,
            bases: Vec::new(),
            covered: Vec::new(),
            keep: Vec::new(),
            seen: Vec::new(),
            size: None,
        })
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn invalid_utf8_source_is_rejected_before_species_filtering() {
        assert!(
            sequence(b"s \xff.c 0 1 + 10 A")
                .err()
                .unwrap()
                .contains("UTF-8")
        );
    }

    #[test]
    fn sequence_payload_has_one_field_and_counts_gaps() {
        let row = sequence(b"s\tref.c 0 2 + 10  A-C \t").unwrap();
        assert_eq!(row.text, b"A-C");
        assert_eq!(row.bases, 2);
        assert!(
            sequence(b"s ref.c 0 2 + 10 AC extra")
                .err()
                .unwrap()
                .contains("8 words")
        );
        assert!(
            sequence(b"s ref.c 0 2 + 10 A\xff")
                .err()
                .unwrap()
                .contains("ASCII")
        );
    }

    #[test]
    fn output_size_must_fit_an_addressable_byte_slice() {
        assert!(
            MsaSize {
                species: usize::MAX,
                len: 2
            }
            .elements()
            .is_err()
        );
        assert!(
            MsaSize {
                species: 1,
                len: isize::MAX as usize + 1
            }
            .elements()
            .is_err()
        );
        assert_eq!(
            MsaSize {
                species: 0,
                len: 10
            }
            .elements()
            .unwrap(),
            0
        );
    }
}
