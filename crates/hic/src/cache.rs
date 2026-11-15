//! Immutable contact blocks prepared for repeated rectangular queries.

use std::ops::Range;

use omicsweave_core::cache::BlockKey;
use omicsweave_core::{Error, Result};

use crate::{BlockEntry, Record, decode_block};

/// Resource limits for one cached block. When explicitly installed on a
/// reader these are strict; default thresholds only control cache admission.
#[derive(Debug, Clone, Copy)]
pub struct HicCacheLimits {
    /// Maximum decompressed wire bytes accepted for one block.
    pub max_inflated_bytes: usize,
    /// Maximum retained, nonzero finite records in one block.
    pub max_records: usize,
}

#[cfg(test)]
mod tests {
    use super::*;

    fn literal(version: i32, records: &[(i16, i16, f32)]) -> Vec<u8> {
        let mut bytes = (records.len() as i32).to_le_bytes().to_vec();
        bytes.extend_from_slice(&[0; 8]);
        bytes.push(1);
        if version == 9 {
            bytes.extend_from_slice(&[0, 0]);
        }
        bytes.push(1);
        bytes.extend_from_slice(&(records.len() as i16).to_le_bytes());
        for &(x, y, value) in records {
            bytes.extend_from_slice(&y.to_le_bytes());
            bytes.extend_from_slice(&1i16.to_le_bytes());
            bytes.extend_from_slice(&x.to_le_bytes());
            bytes.extend_from_slice(&value.to_le_bytes());
        }
        bytes
    }

    #[test]
    fn prepared_queries_match_literal_stored_order_with_duplicates_and_mirroring() {
        let records = [
            (9, 2, 7.0),
            (2, 9, 1.0),
            (3, 3, 4.0),
            (9, 2, -3.0),
            (2, 2, f32::NAN),
            (-1, 3, 7.0),
            (3, 2, 0.0),
            (2, 9, 2.0),
            (2, 2, f32::INFINITY),
        ];
        for version in [8, 9] {
            let bytes = literal(version, &records);
            let estimate = PreparedBlock::estimated_heap_bytes(version, &bytes).unwrap();
            let mut understated = bytes.clone();
            understated[..4].copy_from_slice(&0i32.to_le_bytes());
            assert_eq!(
                PreparedBlock::estimated_heap_bytes(version, &understated),
                Some(estimate)
            );
            let block = PreparedBlock::decode(version, &bytes, 100)
                .unwrap()
                .unwrap();
            for intra in [false, true] {
                for x in [0..3, 2..10, 8..12, 2..2] {
                    for y in [0..4, 2..10, 8..12, 3..3] {
                        let mut expected = Vec::new();
                        for &(a, b, value) in &records {
                            let (Ok(a), Ok(b)) = (u32::try_from(a), u32::try_from(b)) else {
                                continue;
                            };
                            if !value.is_finite() || value == 0.0 {
                                continue;
                            }
                            if x.contains(&a) && y.contains(&b) {
                                expected.push((a, b, value));
                            }
                            if intra && a != b && x.contains(&b) && y.contains(&a) {
                                expected.push((b, a, value));
                            }
                        }
                        let mut actual = Vec::new();
                        block.visit(&x, &y, intra, &mut Vec::new(), |a, b, value| {
                            actual.push((a, b, value))
                        });
                        assert_eq!(actual, expected);
                    }
                }
            }
            assert!(PreparedBlock::decode(version, &bytes[..bytes.len() - 1], 100).is_err());
            assert!(PreparedBlock::decode(version, &bytes, 1).unwrap().is_none());
        }
    }
}

impl Default for HicCacheLimits {
    fn default() -> Self {
        Self {
            max_inflated_bytes: 64 << 20,
            max_records: 4 << 20,
        }
    }
}

impl HicCacheLimits {
    pub(crate) fn check(self) -> Result<()> {
        if self.max_inflated_bytes == 0 || self.max_records == 0 {
            return Err(Error::invalid("hic cache", "block limits must be positive"));
        }
        Ok(())
    }
}

/// One source-bound block dependency; obtain it with `HicReader::cache_plan`.
#[derive(Debug, Clone)]
pub struct HicBlockRequest {
    pub(crate) key: BlockKey,
    pub(crate) entry: BlockEntry,
}

impl HicBlockRequest {
    pub fn key(&self) -> BlockKey {
        self.key
    }

    pub fn entry(&self) -> BlockEntry {
        self.entry
    }
}

#[derive(Debug)]
struct Row {
    y: u32,
    range: Range<usize>,
}

/// An immutable validated block, owned by an external cache lease.
///
/// Its index selects intersecting rows and columns without decoding or
/// traversing unrelated records. Selected records are emitted in their wire
/// order, preserving duplicate-cell assignment and cis-mirroring behavior.
pub struct PreparedBlock {
    wire_bytes: usize,
    records: Vec<Record>,
    order: Vec<usize>,
    rows: Vec<Row>,
}

impl PreparedBlock {
    /// Bounds preparation allocations from row counts or dense dimensions.
    /// The declared top-level record count is not trusted: independent files
    /// can contain more decoded records than that header field advertises.
    pub(crate) fn estimated_heap_bytes(version: i32, bytes: &[u8]) -> Option<usize> {
        fn integer(bytes: &[u8], at: &mut usize, width: usize) -> Option<usize> {
            let end = at.checked_add(width)?;
            let field = bytes.get(*at..end)?;
            *at = end;
            usize::try_from(if width == 2 {
                i16::from_le_bytes(field.try_into().ok()?) as i32
            } else {
                i32::from_le_bytes(field.try_into().ok()?)
            })
            .ok()
        }
        let float = *bytes.get(12)? != 0;
        let (xw, yw, mut at) = if version == 9 {
            (
                if *bytes.get(13)? != 0 { 4 } else { 2 },
                if *bytes.get(14)? != 0 { 4 } else { 2 },
                15,
            )
        } else {
            (2, 2, 13)
        };
        let kind = *bytes.get(at)?;
        at += 1;
        let value_width = if float { 4 } else { 2 };
        let (records, rows) = match kind {
            1 => {
                let rows = integer(bytes, &mut at, yw)?;
                let mut records = 0usize;
                for _ in 0..rows {
                    at = at.checked_add(yw)?;
                    let count = integer(bytes, &mut at, xw)?;
                    records = records.checked_add(count)?;
                    at = at.checked_add(count.checked_mul(xw + value_width)?)?;
                    bytes.get(..at)?;
                }
                (records, rows.min(records))
            }
            2 => {
                let points = integer(bytes, &mut at, 4)?;
                let width = integer(bytes, &mut at, 2)?;
                if points != 0 && width == 0 {
                    return None;
                }
                bytes.get(at..at.checked_add(points.checked_mul(value_width)?)?)?;
                (
                    points,
                    if points == 0 {
                        0
                    } else {
                        points.div_ceil(width)
                    },
                )
            }
            _ => return None,
        };
        let capacity = |n: usize| {
            if n == 0 {
                Some(0)
            } else {
                n.checked_next_power_of_two().map(|n| n.max(4))
            }
        };
        capacity(records)?
            .checked_mul(std::mem::size_of::<Record>())?
            .checked_add(records.checked_mul(std::mem::size_of::<usize>())?)?
            .checked_add(capacity(rows)?.checked_mul(std::mem::size_of::<Row>())?)
    }

    pub(crate) fn decode(version: i32, bytes: &[u8], max_records: usize) -> Result<Option<Self>> {
        let mut records = Vec::new();
        let mut overflow = false;
        let mut allocation_failed = false;
        decode_block(version, bytes, |x, y, value| {
            if allocation_failed {
                return;
            }
            let (Ok(x), Ok(y)) = (u32::try_from(x), u32::try_from(y)) else {
                return;
            };
            if !value.is_finite() || value == 0.0 {
                return;
            }
            if records.len() == max_records {
                overflow = true;
                return;
            }
            if records.len() == records.capacity() && records.try_reserve(1).is_err() {
                allocation_failed = true;
                return;
            }
            records.push(Record { x, y, value });
        })?;
        if allocation_failed {
            return Err(Error::invalid("hic cache", "record allocation failed"));
        }
        if overflow {
            return Ok(None);
        }
        let mut order = Vec::new();
        order
            .try_reserve_exact(records.len())
            .map_err(|_| Error::invalid("hic cache", "index allocation failed"))?;
        order.extend(0..records.len());
        if !records
            .windows(2)
            .all(|w| (w[0].y, w[0].x) <= (w[1].y, w[1].x))
        {
            order.sort_unstable_by_key(|&i| (records[i].y, records[i].x, i));
        }
        let mut rows = Vec::new();
        let mut first = 0;
        while first < order.len() {
            let y = records[order[first]].y;
            let mut end = first + 1;
            while end < order.len() && records[order[end]].y == y {
                end += 1;
            }
            rows.try_reserve(1)
                .map_err(|_| Error::invalid("hic cache", "row allocation failed"))?;
            rows.push(Row {
                y,
                range: first..end,
            });
            first = end;
        }
        Ok(Some(Self {
            wire_bytes: bytes.len(),
            records,
            order,
            rows,
        }))
    }

    /// Retained allocations, including the row/column query index.
    pub fn allocated_bytes(&self) -> usize {
        self.records
            .capacity()
            .saturating_mul(std::mem::size_of::<Record>())
            .saturating_add(
                self.order
                    .capacity()
                    .saturating_mul(std::mem::size_of::<usize>()),
            )
            .saturating_add(
                self.rows
                    .capacity()
                    .saturating_mul(std::mem::size_of::<Row>()),
            )
    }

    pub fn len(&self) -> usize {
        self.records.len()
    }

    pub fn is_empty(&self) -> bool {
        self.records.is_empty()
    }

    pub(crate) fn within(&self, limits: HicCacheLimits) -> bool {
        self.wire_bytes <= limits.max_inflated_bytes && self.records.len() <= limits.max_records
    }

    fn select(&self, x: &Range<u32>, y: &Range<u32>, selected: &mut Vec<usize>) {
        let first = self.rows.partition_point(|r| r.y < y.start);
        for row in self.rows[first..].iter().take_while(|r| r.y < y.end) {
            let entries = &self.order[row.range.clone()];
            let lo = entries.partition_point(|&i| self.records[i].x < x.start);
            let hi = entries.partition_point(|&i| self.records[i].x < x.end);
            selected.extend_from_slice(&entries[lo..hi]);
        }
    }

    pub(crate) fn visit(
        &self,
        x: &Range<u32>,
        y: &Range<u32>,
        intra: bool,
        selected: &mut Vec<usize>,
        mut emit: impl FnMut(u32, u32, f32),
    ) -> usize {
        selected.clear();
        self.select(x, y, selected);
        if intra {
            self.select(y, x, selected);
        }
        if !selected.is_sorted() {
            selected.sort_unstable();
        }
        selected.dedup();
        for &i in selected.iter() {
            let r = self.records[i];
            if x.contains(&r.x) && y.contains(&r.y) {
                emit(r.x, r.y, r.value);
            }
            if intra && r.x != r.y && x.contains(&r.y) && y.contains(&r.x) {
                emit(r.y, r.x, r.value);
            }
        }
        selected.len()
    }
}
