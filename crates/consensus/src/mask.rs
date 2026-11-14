//! Masks (`-m` with `--mask-with`): BED intervals written over the
//! reference before any record applies.

use std::collections::HashMap;

use omicsweave_core::{Error, Result};

/// What a mask writes (`--mask-with`).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum MaskWith {
    /// This character over each base (bcftools' default `N`); records that
    /// overlap the mask are left out.
    Char(u8),
    /// Uppercase the bases (`uc`).
    Upper,
    /// Lowercase the bases (`lc`).
    Lower,
}

/// The intervals of a BED file, merged and sorted per chromosome, with what
/// to write over them.
#[derive(Debug, Clone)]
pub struct Mask {
    /// 0-based inclusive `[begin, end]`.
    intervals: HashMap<String, Vec<(i64, i64)>>,
    with: MaskWith,
}

impl Mask {
    /// Parses BED text: chromosome, start and end, whitespace-separated;
    /// blank and `#` lines skipped, further columns ignored.
    pub fn parse(text: &[u8], with: MaskWith) -> Result<Mask> {
        let mut intervals: HashMap<String, Vec<(i64, i64)>> = HashMap::new();
        for (i, line) in text.split(|&c| c == b'\n').enumerate() {
            let bad = |detail: &str| Error::Format {
                line: i as u64 + 1,
                detail: format!("mask: {detail}"),
            };
            let line = std::str::from_utf8(line).map_err(|_| bad("not UTF-8"))?;
            let line = line.trim();
            if line.is_empty() || line.starts_with('#') {
                continue;
            }
            let mut fields = line.split_whitespace();
            let chrom = fields.next().unwrap_or_default();
            let mut coordinate = |what| -> Result<i64> {
                fields
                    .next()
                    .and_then(|f| f.parse().ok())
                    .ok_or_else(|| bad(what))
            };
            let start = coordinate("bad start")?;
            let end = coordinate("bad end")?;
            if start < 0 || end < start {
                return Err(bad("end before start"));
            }
            if end > start {
                intervals
                    .entry(chrom.to_owned())
                    .or_default()
                    .push((start, end - 1));
            }
        }
        for list in intervals.values_mut() {
            list.sort_unstable();
            let mut merged: Vec<(i64, i64)> = Vec::with_capacity(list.len());
            for &(b, e) in list.iter() {
                match merged.last_mut() {
                    Some(last) if b <= last.1 + 1 => last.1 = last.1.max(e),
                    _ => merged.push((b, e)),
                }
            }
            *list = merged;
        }
        Ok(Mask { intervals, with })
    }

    pub fn with(&self) -> MaskWith {
        self.with
    }

    /// Whether records overlapping the mask are left out.
    pub(crate) fn skips_records(&self) -> bool {
        matches!(self.with, MaskWith::Char(_))
    }

    /// Appends mask overlaps for inclusive spans in nondecreasing start
    /// order. Span ends need not be monotone: nested records stay valid.
    pub(crate) fn overlapping_records(
        &self,
        chrom: &str,
        spans: impl IntoIterator<Item = (i64, i64)>,
        out: &mut Vec<bool>,
    ) {
        let mut list = self.intervals.get(chrom).map_or(&[][..], Vec::as_slice);
        out.extend(spans.into_iter().map(|(start, end)| {
            if list.first().is_some_and(|&(_, e)| e < start) {
                // Nearby records usually pass few mask intervals. Exponential
                // search also skips a distant BED prefix in logarithmic time.
                let mut after = 1usize;
                while after < list.len() && list[after].1 < start {
                    after = after.saturating_mul(2);
                }
                let search = after.min(list.len());
                let first = list[..search].partition_point(|&(_, e)| e < start);
                list = &list[first..];
            }
            end >= start && list.first().is_some_and(|&(b, _)| b <= end)
        }));
    }

    /// Writes the mask over `seq`, the reference from `start`.
    pub(crate) fn apply(&self, chrom: &str, start: i64, seq: &mut [u8]) {
        if seq.is_empty() {
            return;
        }
        let Some(list) = self.intervals.get(chrom) else {
            return;
        };
        let end = start + seq.len() as i64 - 1;
        let first = list.partition_point(|&(_, e)| e < start);
        for &(b, e) in list[first..].iter().take_while(|&&(b, _)| b <= end) {
            let span = &mut seq[(b.max(start) - start) as usize..=(e.min(end) - start) as usize];
            match self.with {
                MaskWith::Char(c) => span.fill(c),
                MaskWith::Upper => span.make_ascii_uppercase(),
                MaskWith::Lower => span.make_ascii_lowercase(),
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::{Mask, MaskWith};
    use std::fmt::Write as _;

    #[test]
    fn sorted_record_masks_match_all_pairs_overlap_with_nested_ends() {
        let mut bed = String::new();
        let mut intervals = Vec::new();
        for i in 0..257 {
            let begin = i * 11;
            let end = begin + 2 + i % 5;
            intervals.push((begin, end));
            writeln!(bed, "c\t{begin}\t{end}").unwrap();
        }
        let mask = Mask::parse(bed.as_bytes(), MaskWith::Char(b'N')).unwrap();
        let mut spans = Vec::new();
        for start in -1..2850 {
            for len in [3000, 1, 0, 2, 17] {
                spans.push((start, start + len - 1));
            }
        }
        // A separate sparse batch forces searches across distant BED prefixes.
        for batch in [
            spans,
            vec![(1200, 2500), (1200, 1200), (2800, 2801), (9000, 9001)],
        ] {
            let expected: Vec<_> = batch
                .iter()
                .map(|&(start, end)| {
                    end >= start && intervals.iter().any(|&(b, e)| b <= end && e > start)
                })
                .collect();
            let mut actual = vec![true, false];
            mask.overlapping_records("c", batch.iter().copied(), &mut actual);
            assert_eq!(&actual[..2], &[true, false]);
            assert_eq!(&actual[2..], expected);
            actual.clear();
            mask.overlapping_records("missing", batch.iter().copied(), &mut actual);
            assert_eq!(actual, vec![false; batch.len()]);
        }
    }
}
