//! The `.crai` index: one gzip-compressed text line per slice (per
//! reference of a multi-reference slice), loaded as htslib's
//! `cram_index_load` loads it and searched as `cram_index_query` searches
//! it.

use std::io::Read;
use std::path::{Path, PathBuf};

use omicsweave_core::{Error, Result};

/// One line of a `.crai`.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct CraiEntry {
    /// Reference id; -1 for unplaced reads.
    pub ref_id: i32,
    /// 1-based first position of the slice's alignments.
    pub start: i64,
    pub span: i64,
    /// File offset of the slice's container.
    pub container: i64,
    /// Offset of the slice from the end of the container header.
    pub landmark: i32,
    /// Bytes of the slice.
    pub slice_len: i32,
}

impl CraiEntry {
    /// The line htslib writes for this entry.
    pub fn line(&self) -> String {
        format!(
            "{}\t{}\t{}\t{}\t{}\t{}\n",
            self.ref_id, self.start, self.span, self.container, self.landmark, self.slice_len
        )
    }
}

/// A top-level entry of a reference's nested containment list.
#[derive(Debug, Clone, Copy)]
struct Top {
    start: i32,
    end: i32,
    offset: i64,
}

/// A loaded `.crai`.
#[derive(Debug, Clone, Default)]
pub struct CramIndex {
    entries: Vec<CraiEntry>,
    /// By reference id + 1: entries not contained in an earlier one.
    tops: Vec<Vec<Top>>,
    monotone: Vec<bool>,
}

fn invalid(detail: impl Into<String>) -> Error {
    Error::invalid("CRAM index", detail)
}

/// `kget_int32` / `kget_int64`: blanks, an optional `-`, digits.
fn next_int(text: &[u8], pos: &mut usize) -> Option<i64> {
    let mut p = *pos;
    while p < text.len() && (text[p] == b' ' || text[p] == b'\t') {
        p += 1;
    }
    let mut negative = false;
    if p < text.len() && text[p] == b'-' {
        negative = true;
        p += 1;
    }
    if p >= text.len() || !text[p].is_ascii_digit() {
        return None;
    }
    let mut v: i64 = 0;
    while p < text.len() && text[p].is_ascii_digit() {
        v = v.checked_mul(10)?;
        let digit = (text[p] - b'0') as i64;
        v = if negative {
            v.checked_sub(digit)?
        } else {
            v.checked_add(digit)?
        };
        p += 1;
    }
    *pos = p;
    Some(v)
}

/// The text of index bytes: gzip (one member or several, BGZF included)
/// or plain.
fn text_of(raw: &[u8]) -> Result<Vec<u8>> {
    if !raw.starts_with(&[0x1f, 0x8b]) {
        return Ok(raw.to_vec());
    }
    let mut out = Vec::new();
    flate2::read::MultiGzDecoder::new(raw)
        .read_to_end(&mut out)
        .map_err(|_| Error::Inflate("CRAM index does not inflate"))?;
    Ok(out)
}

impl CramIndex {
    /// Parses a `.crai` (gzip-compressed or plain text) of a file with
    /// `n_refs` references.
    pub fn parse(raw: &[u8], n_refs: usize) -> Result<CramIndex> {
        let text = text_of(raw)?;
        if text.is_empty() {
            return Err(invalid("empty index"));
        }
        if n_refs > i32::MAX as usize {
            return Err(invalid("too many references"));
        }
        let mut idx = CramIndex {
            entries: Vec::new(),
            tops: vec![Vec::new(); n_refs + 1],
            monotone: Vec::new(),
        };
        // The containment stack: (refid, start, end); the bottom is the
        // reference's root.
        let mut stack: Vec<(i32, i32, i32)> = vec![(-1, i32::MIN, i32::MAX)];
        let mut pos = 0;
        while pos < text.len() {
            let mut field = || next_int(&text, &mut pos).ok_or_else(|| invalid("malformed line"));
            let fields = [field()?, field()?, field()?, field()?, field()?, field()?];
            let int32 = |n| i32::try_from(n).map_err(|_| invalid("32-bit field overflow"));
            let ref_id = int32(fields[0])?;
            let start = int32(fields[1])?;
            let span = int32(fields[2])?;
            let container = fields[3];
            let landmark = int32(fields[4])?;
            let slice_len = int32(fields[5])?;
            if start < 0 || span < 0 || container < 0 || landmark < 0 || slice_len < 0 {
                return Err(invalid("negative coordinate, offset or extent"));
            }
            let end = int32(start as i64 + span as i64 - 1)?;
            if ref_id < -1 || ref_id as i64 >= n_refs as i64 {
                return Err(invalid(format!("reference id {ref_id}")));
            }
            let current = stack.last().unwrap().0;
            if ref_id != current {
                if !idx.tops[(ref_id + 1) as usize].is_empty() {
                    return Err(invalid("index is not sorted"));
                }
                stack.clear();
                stack.push((ref_id, i32::MIN, i32::MAX));
            }
            loop {
                let &(r, s, e) = stack.last().unwrap();
                let contains = start >= s && end <= e;
                if contains && !(s == 0 && r == -1) {
                    break;
                }
                stack.pop();
                if stack.is_empty() {
                    return Err(invalid("index entries do not nest"));
                }
            }
            if stack.len() == 1 {
                idx.tops[(ref_id + 1) as usize].push(Top {
                    start,
                    end,
                    offset: container,
                });
            }
            stack.push((ref_id, start, end));
            idx.entries.push(CraiEntry {
                ref_id,
                start: start as i64,
                span: span as i64,
                container,
                landmark,
                slice_len,
            });
            while pos < text.len() && text[pos] != b'\n' {
                pos += 1;
            }
            pos += 1;
        }
        idx.monotone = idx
            .tops
            .iter()
            .map(|tops| {
                tops.windows(2)
                    .all(|w| w[0].start <= w[1].start && w[0].end <= w[1].end)
            })
            .collect();
        Ok(idx)
    }

    /// Reads the `.crai` at `path`.
    pub fn read_file(path: impl AsRef<Path>, n_refs: usize) -> Result<CramIndex> {
        let path = path.as_ref();
        let raw = std::fs::read(path).map_err(|e| Error::Io {
            kind: e.kind(),
            detail: format!("{}: {e}", path.display()),
        })?;
        CramIndex::parse(&raw, n_refs)
    }

    /// Every entry, in file order.
    pub fn entries(&self) -> &[CraiEntry] {
        &self.entries
    }

    /// `cram_index_query(fd, refid, pos, NULL)`: the container offset to
    /// start a query of `ref_id` at 1-based `pos` from; `None` when the
    /// reference has no entries.
    pub fn query(&self, ref_id: i32, pos: i64) -> Option<i64> {
        let reference = usize::try_from(ref_id.checked_add(1)?).ok()?;
        let list = self.tops.get(reference)?;
        if list.is_empty() {
            return None;
        }
        let pos = pos.clamp(i32::MIN as i64, i32::MAX as i64) as i32;
        if self.monotone[reference] {
            return Some(
                list[list
                    .partition_point(|top| top.end < pos)
                    .min(list.len() - 1)]
                .offset,
            );
        }
        let (mut i, mut j) = (0usize, list.len() - 1);
        let mut k = j / 2;
        while k != i {
            if list[k].start >= pos {
                j = k;
            } else {
                i = k;
            }
            k = (j - i) / 2 + i;
        }
        if list[j].start < pos {
            i = j;
        }
        while i > 0 && list[i - 1].end >= pos {
            i -= 1;
        }
        while i + 1 < list.len() && list[i].end < pos {
            i += 1;
        }
        Some(list[i].offset)
    }
}

/// The index htslib finds beside `cram`: `<cram>.crai`, else the name with
/// its extension replaced by `.crai`.
pub fn locate_index(cram: impl AsRef<Path>) -> Option<PathBuf> {
    let cram = cram.as_ref();
    let mut beside = cram.as_os_str().to_owned();
    beside.push(".crai");
    let beside = PathBuf::from(beside);
    if beside.is_file() {
        return Some(beside);
    }
    let replaced = cram.with_extension("crai");
    (cram.extension().is_some() && replaced.is_file()).then_some(replaced)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn monotone_and_irregular_tops_match_the_linear_search() {
        for n in 1usize..100 {
            for mode in 0..4 {
                let tops: Vec<_> = (0..n)
                    .map(|i| Top {
                        start: if mode == 3 {
                            (i * 17 % n) as i32
                        } else {
                            i as i32 * 3
                        },
                        end: match mode {
                            0 => i as i32 * 3 + 2,
                            1 => i as i32 + n as i32,
                            _ => (i * 17 % n) as i32 + n as i32,
                        },
                        offset: i as i64 * 80,
                    })
                    .collect();
                let monotone = tops
                    .windows(2)
                    .all(|w| w[0].start <= w[1].start && w[0].end <= w[1].end);
                let idx = CramIndex {
                    entries: Vec::new(),
                    tops: vec![Vec::new(), tops],
                    monotone: vec![true, monotone],
                };
                let mut linear = idx.clone();
                linear.monotone.fill(false);
                for pos in (-3i64..(n * 4) as i64).chain([i64::MIN, i64::MAX]) {
                    assert_eq!(idx.query(0, pos), linear.query(0, pos));
                }
                assert_eq!(idx.query(i32::MAX, 0), None);
            }
        }
    }

    #[test]
    fn numeric_overflow_and_negative_extents_are_errors() {
        for text in [
            "4294967296\t1\t1\t10\t0\t1\n",
            "0\t4294967296\t1\t10\t0\t1\n",
            "0\t1\t4294967296\t10\t0\t1\n",
            "0\t1\t1\t9223372036854775808\t0\t1\n",
            "0\t1\t1\t-9223372036854775809\t0\t1\n",
            "0\t2147483647\t2\t10\t0\t1\n",
            "0\t1\t-1\t10\t0\t1\n",
            "0\t1\t1\t-10\t0\t1\n",
            "0\t1\t1\t10\t0\t-1\n",
        ] {
            assert!(CramIndex::parse(text.as_bytes(), 1).is_err(), "{text}");
        }
        assert!(CramIndex::parse(b"0\t2147483647\t1\t10\t0\t1\n", 1).is_ok());
        let mut at = 0;
        assert_eq!(next_int(b"-9223372036854775808", &mut at), Some(i64::MIN));
    }

    #[test]
    fn nested_entries_are_skipped_by_the_search() {
        let text =
            b"0\t1\t100\t10\t0\t5\n0\t20\t10\t20\t0\t5\n0\t150\t100\t30\t0\t5\n1\t5\t5\t40\t0\t5\n";
        let idx = CramIndex::parse(text, 2).unwrap();
        assert_eq!(idx.entries().len(), 4);
        assert_eq!(idx.query(0, 1), Some(10));
        assert_eq!(idx.query(0, 50), Some(10));
        assert_eq!(idx.query(0, 120), Some(30));
        assert_eq!(idx.query(0, 1000), Some(30));
        assert_eq!(idx.query(1, 1), Some(40));
        assert_eq!(idx.query(-1, 0), None);
        assert!(
            CramIndex::parse(b"1\t1\t1\t0\t0\t1\n0\t1\t1\t0\t0\t1\n1\t5\t1\t0\t0\t1\n", 2).is_err()
        );
    }
}
