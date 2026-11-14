//! BED, loaded whole: BED lines vary in length and the format has no index,
//! so it is not read randomly. [`BedTable::load`] parses `.bed` text or
//! gzip `.bed.gz` in one pass and builds a columnar table sorted for region
//! queries; queries then follow the standard
//! [`omicsweave_core::IntervalRead`], with the same output as bigBed.
//!
//! Line rules: empty lines, `#` comments and `track` / `browser` lines are
//! skipped; fields are tab-separated, or whitespace-separated on lines
//! without a tab; a trailing `\r` is ignored. Coordinates are 0-based,
//! half-open `u32`, with `start <= end`. Chromosome names must be UTF-8.
//! A record's annotation is its text
//! after the third field (fields joined by tabs), as in bigBed.

mod gz;

use std::collections::HashMap;
use std::ops::Range;

pub use gz::{is_gzip, text_chunks};
use omicsweave_core::{
    Error, IntervalEntry, IntervalRead, IntervalSize, IntervalsOut, Region, Result, Strings,
};

const ANNOTATION_BLOCK: usize = 64;

/// A loaded BED file in split layout: interval entries grouped by
/// chromosome and sorted by start (stable, so ties keep file order), each
/// chromosome's balanced interval index for overlap search, and the annotations
/// as an Arrow `LargeString` column in file order, referenced by
/// [`IntervalEntry::annotation`].
#[derive(Debug, Clone, Default)]
pub struct BedTable {
    chroms: Vec<String>,
    by_name: HashMap<String, usize>,
    /// Entries of each chromosome, by chromosome id.
    ranges: Vec<Range<usize>>,
    entries: Vec<IntervalEntry>,
    /// Monotone ends, or subtree maxima of an implicit balanced tree.
    max_end: Vec<u32>,
    ordered_ends: Vec<bool>,
    full_spans: Vec<FullSpan>,
    contiguous_annotations: Vec<bool>,
    /// Byte totals at sorted-row block boundaries for noncontiguous monotone groups.
    annotation_blocks: Vec<Vec<usize>>,
    offsets: Vec<i64>,
    data: Vec<u8>,
    /// Reusable positions for nonmonotone tree queries.
    hits: Vec<usize>,
    staged: Option<Staged>,
}

#[derive(Debug, Clone, Copy)]
struct FullSpan {
    min_end: u32,
    annotation_bytes: usize,
}

#[derive(Debug, Clone)]
struct Staged {
    size: IntervalSize,
    entries: StagedEntries,
}

#[derive(Debug, Clone)]
enum StagedEntries {
    /// A monotone-end query produces a contiguous slice without hit storage.
    Span(Range<usize>),
    /// The private chromosome certificate also proves consecutive file rows.
    ContiguousSpan(Range<usize>),
    Hits,
}

fn error(line: u64, detail: impl Into<String>) -> Error {
    Error::Format {
        line,
        detail: detail.into(),
    }
}

fn is_header(line: &[u8]) -> bool {
    [&b"track"[..], b"browser"]
        .iter()
        .any(|w| line.starts_with(w) && line.get(w.len()).is_none_or(|c| matches!(c, b' ' | b'\t')))
}

fn coordinate(field: &[u8], what: &str, line: u64) -> Result<u32> {
    std::str::from_utf8(field)
        .ok()
        .and_then(|s| s.parse().ok())
        .ok_or_else(|| {
            error(
                line,
                format!("bad {what} {:?}", String::from_utf8_lossy(field)),
            )
        })
}

/// Accumulates records while the text streams in.
#[derive(Default)]
struct Builder {
    table: BedTable,
    /// (chromosome id, start, end) per record, in file order.
    records: Vec<(u32, u32, u32)>,
    carry: Vec<u8>,
    line: u64,
    joined: Vec<u8>,
}

impl Builder {
    fn record(&mut self, line: &[u8]) -> Result<()> {
        let n = self.line;
        let line = line.strip_suffix(b"\r").unwrap_or(line);
        if line.is_empty() || line[0] == b'#' || is_header(line) {
            return Ok(());
        }
        let (chrom, start, end, rest): (&[u8], &[u8], &[u8], &[u8]);
        if memchr::memchr(b'\t', line).is_some() {
            let mut fields = line.splitn(4, |&c| c == b'\t');
            chrom = fields.next().unwrap_or_default();
            start = fields
                .next()
                .ok_or_else(|| error(n, "fewer than 3 fields"))?;
            end = fields
                .next()
                .ok_or_else(|| error(n, "fewer than 3 fields"))?;
            rest = fields.next().unwrap_or_default();
        } else {
            let mut fields = line
                .split(|c| c.is_ascii_whitespace())
                .filter(|f| !f.is_empty());
            chrom = fields
                .next()
                .ok_or_else(|| error(n, "fewer than 3 fields"))?;
            start = fields
                .next()
                .ok_or_else(|| error(n, "fewer than 3 fields"))?;
            end = fields
                .next()
                .ok_or_else(|| error(n, "fewer than 3 fields"))?;
            self.joined.clear();
            for (i, f) in fields.enumerate() {
                if i > 0 {
                    self.joined.push(b'\t');
                }
                self.joined.extend_from_slice(f);
            }
            rest = &[];
        }
        if chrom.is_empty() {
            return Err(error(n, "empty chromosome name"));
        }
        let (s, e) = (coordinate(start, "start", n)?, coordinate(end, "end", n)?);
        if s > e {
            return Err(error(n, format!("start {s} after end {e}")));
        }
        let table = &mut self.table;
        let name =
            std::str::from_utf8(chrom).map_err(|_| error(n, "chromosome name is not UTF-8"))?;
        let id = match table.by_name.get(name) {
            Some(&id) => id,
            None => {
                let id = table.chroms.len();
                table.by_name.insert(name.to_string(), id);
                table.chroms.push(name.to_owned());
                id
            }
        };
        let annotation = if rest.is_empty() {
            &self.joined[..]
        } else {
            rest
        };
        table.data.extend_from_slice(annotation);
        table.offsets.push(
            i64::try_from(table.data.len()).map_err(|_| error(n, "annotation data exceeds i64"))?,
        );
        self.joined.clear();
        if self.records.len() >= u32::MAX as usize {
            return Err(error(n, "more than 2^32 - 1 records"));
        }
        self.records.push((id as u32, s, e));
        Ok(())
    }

    fn push(&mut self, chunk: &[u8]) -> Result<()> {
        let mut rest = chunk;
        if !self.carry.is_empty() {
            let Some(p) = memchr::memchr(b'\n', rest) else {
                self.carry.extend_from_slice(rest);
                return Ok(());
            };
            let mut line = std::mem::take(&mut self.carry);
            line.extend_from_slice(&rest[..p]);
            self.line += 1;
            self.record(&line)?;
            line.clear();
            self.carry = line;
            rest = &rest[p + 1..];
        }
        while let Some(p) = memchr::memchr(b'\n', rest) {
            self.line += 1;
            self.record(&rest[..p])?;
            rest = &rest[p + 1..];
        }
        self.carry.extend_from_slice(rest);
        Ok(())
    }

    fn finish(mut self) -> Result<BedTable> {
        if !self.carry.is_empty() {
            let line = std::mem::take(&mut self.carry);
            self.line += 1;
            self.record(&line)?;
        }
        let mut table = self.table;
        // Group by chromosome (counting sort, stable), then sort each group
        // by start (stable).
        let mut counts = vec![0usize; table.chroms.len()];
        for &(c, _, _) in &self.records {
            counts[c as usize] += 1;
        }
        let mut next = Vec::with_capacity(counts.len());
        let mut at = 0;
        for &count in &counts {
            table.ranges.push(at..at + count);
            next.push(at);
            at += count;
        }
        table.entries = vec![IntervalEntry::default(); self.records.len()];
        for (row, &(c, start, end)) in self.records.iter().enumerate() {
            let slot = &mut next[c as usize];
            table.entries[*slot] = IntervalEntry {
                start,
                end,
                annotation: row as u32,
            };
            *slot += 1;
        }
        table.max_end = vec![0; table.entries.len()];
        let mut sort_scratch = Vec::new();
        for range in &table.ranges {
            let group = &mut table.entries[range.clone()];
            // Stable chromosome grouping leaves file row ids strictly increasing.
            // Their endpoints certify contiguity only before sorting. A start
            // inversion must permute these unique ids, invalidating that order.
            let contiguous = group.first().zip(group.last()).is_none_or(|(first, last)| {
                last.annotation as usize - first.annotation as usize + 1 == group.len()
            });
            let unchanged = sort_entries(group, &mut sort_scratch);
            let mut ordered = true;
            let mut previous_end = 0;
            let mut min_end = u32::MAX;
            let maxima = &mut table.max_end[range.clone()];
            for (dst, entry) in maxima.iter_mut().zip(group.iter()) {
                ordered &= previous_end <= entry.end;
                previous_end = entry.end;
                min_end = min_end.min(entry.end);
                *dst = entry.end;
            }
            table.ordered_ends.push(ordered);
            table.contiguous_annotations.push(contiguous && unchanged);
            let mut annotation_blocks = Vec::new();
            if ordered && !(contiguous && unchanged) && group.len() >= ANNOTATION_BLOCK {
                // Only the two partial edge blocks need per-row lengths during
                // a query. All complete blocks use these bounded byte totals.
                annotation_blocks.reserve(group.len() / ANNOTATION_BLOCK + 1);
                annotation_blocks.push(0);
                let mut bytes = 0;
                for block in group.chunks_exact(ANNOTATION_BLOCK) {
                    bytes += annotation_bytes(block, &table.offsets);
                    annotation_blocks.push(bytes);
                }
            }
            table.annotation_blocks.push(annotation_blocks);
            let annotation_bytes = if ordered {
                0
            } else if contiguous && unchanged {
                let first = group[0].annotation as usize;
                (table.offsets[first + group.len()] - table.offsets[first]) as usize
            } else {
                // Every file row occurs once in a chromosome, so this sum
                // cannot exceed the validated total annotation byte length.
                group
                    .iter()
                    .map(|entry| {
                        let row = entry.annotation as usize;
                        (table.offsets[row + 1] - table.offsets[row]) as usize
                    })
                    .sum()
            };
            table.full_spans.push(FullSpan {
                min_end,
                annotation_bytes,
            });
            if !ordered {
                build_max(group, maxima);
            }
        }
        Ok(table)
    }
}

/// Stable LSD radix passes keep file order for equal u32 starts. Small
/// groups use comparison sorting; sorted input needs no temporary storage.
/// Returns whether the input was already sorted and retained its row order.
fn sort_entries(entries: &mut [IntervalEntry], scratch: &mut Vec<IntervalEntry>) -> bool {
    if entries.windows(2).all(|w| w[0].start <= w[1].start) {
        return true;
    }
    if entries.len() < 4096 {
        entries.sort_by_key(|e| e.start);
        return false;
    }
    let first = entries[0].start;
    let changed = entries.iter().fold(0, |bits, e| bits | (e.start ^ first));
    scratch.resize(entries.len(), IntervalEntry::default());
    let mut in_scratch = false;
    for shift in (0..32).step_by(8) {
        if (changed >> shift) & 255 == 0 {
            continue;
        }
        if in_scratch {
            radix_pass(&scratch[..entries.len()], entries, shift);
        } else {
            radix_pass(entries, &mut scratch[..entries.len()], shift);
        }
        in_scratch = !in_scratch;
    }
    if in_scratch {
        entries.copy_from_slice(&scratch[..entries.len()]);
    }
    false
}

fn radix_pass(src: &[IntervalEntry], dst: &mut [IntervalEntry], shift: u32) {
    assert_eq!(src.len(), dst.len());
    let mut next = [0usize; 256];
    for entry in src {
        next[((entry.start >> shift) & 255) as usize] += 1;
    }
    let mut at = 0;
    for count in &mut next {
        let n = *count;
        *count = at;
        at += n;
    }
    for &entry in src {
        let slot = &mut next[((entry.start >> shift) & 255) as usize];
        // SAFETY: the counting pass assigns disjoint bucket ranges whose
        // lengths sum to src.len() == dst.len(). Each source item consumes
        // one slot in its bucket, so the slot is in bounds before increment.
        // Rust's disjoint shared/mutable borrows exclude aliasing; IntervalEntry
        // is Copy and the destination slice is initialized and aligned.
        unsafe {
            *dst.get_unchecked_mut(*slot) = entry;
        }
        *slot += 1;
    }
}

impl BedTable {
    /// Loads a whole `.bed` or gzip `.bed.gz` (detected by its magic bytes;
    /// any number of members, bgzip included), inflating in chunks.
    pub fn load(raw: &[u8]) -> Result<BedTable> {
        let mut builder = Builder::default();
        builder.table.offsets.push(0);
        text_chunks(raw, |chunk| builder.push(chunk))?;
        builder.finish()
    }

    /// Chromosome names, in order of first appearance.
    pub fn chroms(&self) -> &[String] {
        &self.chroms
    }

    /// Records in the file.
    pub fn len(&self) -> usize {
        self.entries.len()
    }

    pub fn is_empty(&self) -> bool {
        self.entries.is_empty()
    }

    /// The sorted entries of `chrom`; empty if it has none.
    pub fn entries(&self, chrom: &str) -> &[IntervalEntry] {
        self.by_name
            .get(chrom)
            .map_or(&[], |&id| &self.entries[self.ranges[id].clone()])
    }

    /// All annotations, by file row ([`IntervalEntry::annotation`]).
    pub fn annotations(&self) -> Strings<'_> {
        Strings {
            offsets: &self.offsets,
            data: &self.data,
        }
    }

    fn stage_span(&self, span: Range<usize>, chrom_id: usize) -> Staged {
        let entries = &self.entries[span.clone()];
        let Some(first) = entries.first() else {
            return Staged {
                size: IntervalSize::default(),
                entries: StagedEntries::Span(span),
            };
        };
        let mut bulk = false;
        let annotation_bytes = if entries.len() == 1 || self.contiguous_annotations[chrom_id] {
            let lo = first.annotation as usize;
            // SAFETY: load bounds each annotation id and stores one extra
            // offset. A singleton uses id+1; otherwise the private chromosome
            // flag certifies every adjacent annotation id is consecutive, so
            // lo+entries.len() is the last id+1 and cannot overflow usize.
            // Offsets are immutable, monotone and bounded by data.len().
            let bytes = unsafe {
                (*self.offsets.get_unchecked(lo + entries.len()) - *self.offsets.get_unchecked(lo))
                    as usize
            };
            bulk = entries.len() > 1;
            bytes
        } else {
            let blocks = &self.annotation_blocks[chrom_id];
            let lo = span.start - self.ranges[chrom_id].start;
            let hi = span.end - self.ranges[chrom_id].start;
            let first_block = lo.div_ceil(ANNOTATION_BLOCK);
            let last_block = hi / ANNOTATION_BLOCK;
            if !blocks.is_empty() && first_block < last_block {
                // Both boundaries lie inside this query span. Thus these
                // products cannot overflow even for a maximal record count.
                let first = first_block * ANNOTATION_BLOCK - lo;
                let last = last_block * ANNOTATION_BLOCK - lo;
                blocks[last_block] - blocks[first_block]
                    + annotation_bytes(&entries[..first], &self.offsets)
                    + annotation_bytes(&entries[last..], &self.offsets)
            } else {
                annotation_bytes(entries, &self.offsets)
            }
        };
        Staged {
            size: IntervalSize {
                entries: entries.len(),
                annotation_bytes,
            },
            entries: if bulk {
                StagedEntries::ContiguousSpan(span)
            } else {
                StagedEntries::Span(span)
            },
        }
    }

    #[inline(never)]
    fn write_contiguous(&self, out: IntervalsOut<'_>, entries: &[IntervalEntry]) {
        // The staged certificate guarantees consecutive file rows, including
        // any empty annotations. All slices retain checked safe boundaries.
        let first = entries[0].annotation as usize;
        let offsets = &self.offsets[first..first + entries.len() + 1];
        let base = offsets[0];
        let data = &self.data[base as usize..offsets[entries.len()] as usize];
        out.annotations.data.copy_from_slice(data);
        out.annotations.offsets[0] = 0;
        for (row, ((dst, entry), (dst_offset, src_offset))) in out
            .entries
            .iter_mut()
            .zip(entries)
            .zip(out.annotations.offsets[1..].iter_mut().zip(&offsets[1..]))
            .enumerate()
        {
            *dst = IntervalEntry {
                annotation: row as u32,
                ..*entry
            };
            *dst_offset = *src_offset - base;
        }
    }

    fn write_entries(&self, out: IntervalsOut<'_>, entries: impl Iterator<Item = IntervalEntry>) {
        out.annotations.offsets[0] = 0;
        let mut byte = 0;
        for (row, e) in entries.enumerate() {
            // SAFETY: write supplies exactly the entries staged by a successful
            // query, with load-validated annotation ids and monotone offsets
            // inside data. Query sums each annotation length once, bounded by
            // data.len() <= i64::MAX, and write checks every destination size.
            // The record count fits u32. Exclusive output borrows cannot alias
            // this immutable table or each other; copy_from_slice writes only
            // initialized u8 data and all offsets/entries are aligned values.
            unsafe {
                let id = e.annotation as usize;
                let lo = *self.offsets.get_unchecked(id) as usize;
                let hi = *self.offsets.get_unchecked(id + 1) as usize;
                let note = self.data.get_unchecked(lo..hi);
                let end = byte + note.len();
                out.annotations
                    .data
                    .get_unchecked_mut(byte..end)
                    .copy_from_slice(note);
                *out.entries.get_unchecked_mut(row) = IntervalEntry {
                    annotation: row as u32,
                    ..e
                };
                *out.annotations.offsets.get_unchecked_mut(row + 1) = end as i64;
                byte = end;
            }
        }
    }
}

fn annotation_bytes(entries: &[IntervalEntry], offsets: &[i64]) -> usize {
    entries
        .iter()
        .map(|entry| {
            let id = entry.annotation as usize;
            (offsets[id + 1] - offsets[id]) as usize
        })
        .sum()
}

fn build_max(entries: &[IntervalEntry], maxima: &mut [u32]) -> u32 {
    if entries.is_empty() {
        return 0;
    }
    let mid = entries.len() / 2;
    let (left, rest) = maxima.split_at_mut(mid);
    let (root, right) = rest.split_at_mut(1);
    let max = entries[mid]
        .end
        .max(build_max(&entries[..mid], left))
        .max(build_max(&entries[mid + 1..], right));
    root[0] = max;
    max
}

/// Median roots bound recursion to log2(N); subtree maxima prune nested
/// intervals without a long enclosing interval forcing a prefix scan.
fn visit_tree(
    entries: &[IntervalEntry],
    maxima: &[u32],
    base: usize,
    region: &Region<'_>,
    emit: &mut impl FnMut(usize),
) {
    if entries.is_empty() {
        return;
    }
    let mid = entries.len() / 2;
    // SAFETY: Both slices are equally sized views of a tree built by
    // build_max; mid is in range for every nonempty subtree.
    let (entry, max) = unsafe { (*entries.get_unchecked(mid), *maxima.get_unchecked(mid)) };
    if max <= region.start || entries[0].start >= region.end {
        return;
    }
    if entries.len() <= 16 {
        for (i, entry) in entries.iter().enumerate() {
            if entry.start >= region.end {
                break;
            }
            if entry.end > region.start {
                emit(base + i);
            }
        }
        return;
    }
    visit_tree(&entries[..mid], &maxima[..mid], base, region, emit);
    if entry.start >= region.end {
        return;
    }
    if entry.end > region.start {
        emit(base + mid);
    }
    visit_tree(
        &entries[mid + 1..],
        &maxima[mid + 1..],
        base + mid + 1,
        region,
        emit,
    );
}

impl IntervalRead for BedTable {
    fn query(&mut self, region: Region<'_>) -> Result<IntervalSize> {
        self.staged = None;
        self.hits.clear();
        region.check()?;
        let id = match self.by_name.get(region.chrom) {
            Some(&id) if !region.is_empty() => id,
            _ => {
                let size = IntervalSize::default();
                self.staged = Some(Staged {
                    size,
                    entries: StagedEntries::Span(0..0),
                });
                return Ok(size);
            }
        };
        let range = self.ranges[id].clone();
        let group = &self.entries[range.clone()];
        let max_end = &self.max_end[range.clone()];
        if self.ordered_ends[id] {
            let lo = max_end.partition_point(|&end| end <= region.start);
            let hi = group.partition_point(|e| e.start < region.end);
            let span = range.start + lo..range.start + hi;
            let staged = self.stage_span(span, id);
            let size = staged.size;
            self.staged = Some(staged);
            return Ok(size);
        }
        let summary = self.full_spans[id];
        if region.start < summary.min_end
            && group.last().is_some_and(|last| last.start < region.end)
        {
            // The minimum end and maximum start prove that every record
            // overlaps. Annotation contiguity still needs its independent
            // file-order certificate; a reordered group uses ordinary output.
            let size = IntervalSize {
                entries: group.len(),
                annotation_bytes: summary.annotation_bytes,
            };
            self.staged = Some(Staged {
                size,
                entries: if self.contiguous_annotations[id] {
                    StagedEntries::ContiguousSpan(range)
                } else {
                    StagedEntries::Span(range)
                },
            });
            return Ok(size);
        }
        let hits = &mut self.hits;
        let mut size = IntervalSize::default();
        visit_tree(group, max_end, range.start, &region, &mut |i| {
            size.entries += 1;
            // SAFETY: visit_tree emits each matching table index once; load
            // assigns an in-range annotation id and monotone i64 offsets.
            // Summed lengths are bounded by this table's data length. No
            // annotation bytes need to be accessed to compute their size.
            unsafe {
                let id = self.entries.get_unchecked(i).annotation as usize;
                size.annotation_bytes += (*self.offsets.get_unchecked(id + 1)
                    - *self.offsets.get_unchecked(id))
                    as usize;
            }
            hits.push(i);
        });
        self.staged = Some(Staged {
            size,
            entries: StagedEntries::Hits,
        });
        Ok(size)
    }

    fn write(&self, out: IntervalsOut<'_>) -> Result<()> {
        let staged = self.staged.as_ref().ok_or(Error::NotQueried)?;
        out.check(staged.size)?;
        match &staged.entries {
            StagedEntries::Span(span) => {
                self.write_entries(out, self.entries[span.clone()].iter().copied());
            }
            StagedEntries::ContiguousSpan(span) => {
                self.write_contiguous(out, &self.entries[span.clone()]);
            }
            StagedEntries::Hits => self.write_entries(
                out,
                self.hits.iter().map(|&i| {
                    // SAFETY: successful query records only valid indices emitted
                    // by visit_tree, and private table entries cannot change.
                    unsafe { *self.entries.get_unchecked(i) }
                }),
            ),
        }
        Ok(())
    }
}
