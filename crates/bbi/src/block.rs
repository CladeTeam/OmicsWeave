//! Views over decompressed data blocks. Each view borrows the caller's bytes
//! and decodes items on access; nothing is copied or allocated.

use crate::bytes::{f32_unchecked, require, u16_at, u32_at, u32_unchecked, u64_unchecked};
use crate::chroms::ChromId;
use crate::error::{Error, Result, invalid};
use crate::search::{Order, OrderScan, Point};
use std::iter::FusedIterator;
use std::ops::Range;

pub const WIG_HEADER_LEN: usize = 24;
pub const ZOOM_RECORD_LEN: usize = 32;
const BED_FIXED_LEN: usize = 12;

#[inline]
fn bed_text_len(buf: &[u8]) -> Option<usize> {
    // Keep short BED fields inline; use SIMD search for longer text.
    if buf.len() < 8 {
        return buf.iter().position(|&b| b == 0);
    }
    // SAFETY: The length check guarantees the complete unaligned word.
    let word = unsafe { u64_unchecked(buf, 0) };
    let zeros = word.wrapping_sub(0x0101_0101_0101_0101) & !word & 0x8080_8080_8080_8080;
    if zeros != 0 {
        // The least significant marked byte is the first zero in the LE
        // word; borrow propagation can only mark subsequent bytes.
        return Some(zeros.trailing_zeros() as usize / 8);
    }
    let tail = &buf[8..];
    if tail.first() == Some(&0) {
        return Some(8);
    }
    memchr::memchr(0, tail).map(|at| at + 8)
}

#[inline]
fn overlaps(start: u32, end: u32, query_start: u32, query_end: u32) -> bool {
    query_start < query_end && end > query_start && start < query_end
}

/// How a bigWig section stores its items.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum WigKind {
    /// (start, end, value) per item.
    BedGraph,
    /// (start, value) per item; every item covers `span` bases.
    VarStep,
    /// value per item; item `i` starts at `start + i * step` and covers `span`.
    FixedStep,
}

impl WigKind {
    fn item_len(self) -> usize {
        match self {
            WigKind::BedGraph => 12,
            WigKind::VarStep => 8,
            WigKind::FixedStep => 4,
        }
    }
}

/// One bigWig data block: a section header and its items, all on one
/// chromosome.
#[derive(Debug, Clone, Copy)]
pub struct WigBlock<'a> {
    buf: &'a [u8],
    chrom: ChromId,
    start: u32,
    end: u32,
    step: u32,
    span: u32,
    kind: WigKind,
    count: u16,
    order: Order,
}

impl<'a> WigBlock<'a> {
    /// Parses decompressed block bytes.
    pub fn parse(buf: &'a [u8]) -> Result<WigBlock<'a>> {
        require(buf, WIG_HEADER_LEN, "wig section")?;
        let kind = match buf[20] {
            1 => WigKind::BedGraph,
            2 => WigKind::VarStep,
            3 => WigKind::FixedStep,
            t => return Err(invalid("wig section", format!("unknown type {t}"))),
        };
        let block = WigBlock {
            buf,
            chrom: u32_at(buf, 0),
            start: u32_at(buf, 4),
            end: u32_at(buf, 8),
            step: u32_at(buf, 12),
            span: u32_at(buf, 16),
            kind,
            count: u16_at(buf, 22),
            order: Order::UNKNOWN,
        };
        require(
            buf,
            WIG_HEADER_LEN + block.count as usize * kind.item_len(),
            "wig section items",
        )?;
        if kind == WigKind::FixedStep && block.count > 0 {
            let last = block.start as u64
                + (block.count as u64 - 1) * block.step as u64
                + block.span as u64;
            if last > u32::MAX as u64 {
                return Err(invalid("wig section", "fixedStep items run past 2^32"));
            }
        }
        if kind == WigKind::VarStep {
            let limit = u32::MAX - block.span;
            let items = &buf[WIG_HEADER_LEN..WIG_HEADER_LEN + block.len() * kind.item_len()];
            if items.chunks_exact(8).any(|item| u32_at(item, 0) > limit) {
                return Err(invalid("wig section", "variableStep items run past 2^32"));
            }
        }
        Ok(block)
    }

    /// Prepares this borrowed view for repeated range queries in O(N) time,
    /// without allocating. Sorted records use binary search; unsorted records
    /// retain stored-order semantics. Keep this view with the cached block.
    pub fn prepare(mut self) -> Self {
        if self.kind != WigKind::FixedStep {
            self.order = Order::scan(self.len(), |i| {
                // SAFETY: scan only requests indices below the parsed count.
                let (s, e, _) = unsafe { self.item_unchecked(i) };
                ((self.chrom, s), (self.chrom, e))
            });
        }
        self
    }

    fn candidates(&self, start: u32, end: u64) -> Range<usize> {
        if start as u64 >= end {
            return 0..0;
        }
        if self.kind == WigKind::FixedStep {
            if self.is_empty() {
                return 0..0;
            }
            if self.step == 0 {
                return if self.start as u64 >= end
                    || self.start as u64 + self.span as u64 <= start as u64
                {
                    0..0
                } else {
                    0..self.len()
                };
            }
            let first_end = self.start as u64 + self.span as u64;
            let first = if start as u64 >= first_end {
                (((start as u64 - first_end) / self.step as u64 + 1).min(self.len() as u64))
                    as usize
            } else {
                0
            };
            let last = if end <= self.start as u64 {
                0
            } else {
                (((end - 1 - self.start as u64) / self.step as u64 + 1).min(self.len() as u64))
                    as usize
            };
            return first.min(last)..last;
        }
        self.order.candidates(
            self.len(),
            (self.chrom, start),
            (self.chrom, end.min(u32::MAX as u64) as u32),
            |i| {
                // SAFETY: candidates only requests indices below the parsed count.
                let (s, e, _) = unsafe { self.item_unchecked(i) };
                ((self.chrom, s), (self.chrom, e))
            },
        )
    }

    pub fn chrom(&self) -> ChromId {
        self.chrom
    }

    pub fn kind(&self) -> WigKind {
        self.kind
    }

    /// The section's declared span on its chromosome.
    pub fn range(&self) -> (u32, u32) {
        (self.start, self.end)
    }

    pub fn len(&self) -> usize {
        self.count as usize
    }

    pub fn is_empty(&self) -> bool {
        self.count == 0
    }

    /// Item `i` as (start, end, value).
    ///
    /// # Panics
    /// Panics if `i >= self.len()`.
    #[inline]
    pub fn item(&self, i: usize) -> (u32, u32, f32) {
        assert!(i < self.len(), "wig item index out of bounds");
        // SAFETY: The assertion bounds i; parse validated all item bytes and
        // the step coordinate arithmetic. Fields are private and immutable.
        unsafe { self.item_unchecked(i) }
    }

    /// # Safety
    /// `i < self.len()`.
    #[inline]
    unsafe fn item_unchecked(&self, i: usize) -> (u32, u32, f32) {
        // SAFETY: Each specialization matches the parsed kind and the caller
        // bounds i to the item array validated by parse.
        unsafe {
            match self.kind {
                WigKind::BedGraph => self.item_of_kind::<1>(i),
                WigKind::VarStep => self.item_of_kind::<2>(i),
                WigKind::FixedStep => self.item_of_kind::<3>(i),
            }
        }
    }

    /// # Safety
    /// `i < self.len()` and KIND corresponds to self.kind (1, 2 or 3).
    #[inline]
    unsafe fn item_of_kind<const KIND: u8>(&self, i: usize) -> (u32, u32, f32) {
        let width = match KIND {
            1 => 12,
            2 => 8,
            _ => 4,
        };
        let at = WIG_HEADER_LEN + i * width;
        let buf = self.buf;
        // SAFETY: parse validated count * item_len bytes, and the caller
        // guarantees i < count and the matching layout. Reads are unaligned.
        unsafe {
            match KIND {
                1 => (
                    u32_unchecked(buf, at),
                    u32_unchecked(buf, at + 4),
                    f32_unchecked(buf, at + 8),
                ),
                2 => {
                    let start = u32_unchecked(buf, at);
                    (start, start + self.span, f32_unchecked(buf, at + 4))
                }
                _ => {
                    let start = self.start + i as u32 * self.step;
                    (start, start + self.span, f32_unchecked(buf, at))
                }
            }
        }
    }

    /// Every item, in stored order.
    pub fn intervals(
        &self,
    ) -> impl ExactSizeIterator<Item = (u32, u32, f32)> + DoubleEndedIterator + FusedIterator + use<'a>
    {
        WigIntervals {
            block: *self,
            range: 0..self.len(),
        }
    }

    /// Items on `chrom` overlapping `[start, end)`, unclipped.
    pub fn overlapping(
        &self,
        chrom: ChromId,
        start: u32,
        end: u32,
    ) -> impl DoubleEndedIterator<Item = (u32, u32, f32)> + FusedIterator + use<'a> {
        let range = if chrom == self.chrom {
            self.candidates(start, end as u64)
        } else {
            0..0
        };
        WigIntervals {
            block: *self,
            range,
        }
        .filter(move |&(s, e, _)| overlaps(s, e, start, end))
    }

    /// Writes the value of every base of `chrom` from `start` covered by this
    /// block into `values[base - start]`. Uncovered positions are left as the
    /// caller initialised them; later items overwrite earlier ones. Returns
    /// the number of writes, counting positions covered by multiple items
    /// once per item.
    ///
    /// # Panics
    /// Panics if the logical write count cannot fit usize.
    pub fn fill(&self, chrom: ChromId, start: u32, values: &mut [f32]) -> usize {
        if chrom != self.chrom || values.is_empty() {
            return 0;
        }
        match self.kind {
            WigKind::BedGraph => {
                let range = self.candidates(start, start as u64 + values.len() as u64);
                self.fill_selected::<1>(start, values, range)
            }
            WigKind::VarStep => {
                let range = self.candidates(start, start as u64 + values.len() as u64);
                self.fill_selected::<2>(start, values, range)
            }
            WigKind::FixedStep => {
                if self.span == 0 || self.is_empty() {
                    return 0;
                }
                if self.span > self.step {
                    return self.fill_overlapping_fixed(start, values);
                }
                let end = start as u64 + values.len() as u64;
                let range = self.candidates(start, end);
                let (first, last) = (range.start, range.end);
                if first >= last {
                    return 0;
                }
                #[cfg(target_endian = "little")]
                if self.step == 1 && self.span == 1 {
                    let count = last - first;
                    let dst = (self.start as u64 + first as u64 - start as u64) as usize;
                    // SAFETY: With step = span = 1, [first, last) contains
                    // only whole bases within the output window. parse
                    // validated all count * 4 source bytes. Shared input and
                    // exclusive output cannot alias; u8 needs no alignment,
                    // every f32 bit pattern is valid, and both are LE here.
                    unsafe {
                        std::ptr::copy_nonoverlapping(
                            self.buf.as_ptr().add(WIG_HEADER_LEN + first * 4),
                            values.as_mut_ptr().add(dst).cast::<u8>(),
                            count * 4,
                        );
                    }
                    return count;
                }
                self.fill_kind::<3, false>(start, values, first, last)
            }
        }
    }

    /// Owning readers reuse scratch only when heavily overlapping records
    /// justify its O(L) initialization. Fixed-step and monotone intervals
    /// already paint each output base at most once without scratch.
    pub(crate) fn fill_reusing(
        &self,
        chrom: ChromId,
        start: u32,
        values: &mut [f32],
        scratch: &mut Vec<usize>,
    ) -> Result<usize> {
        if self.kind == WigKind::FixedStep
            || self.len() < 64
            || chrom != self.chrom
            || values.is_empty()
        {
            return Ok(self.fill(chrom, start, values));
        }
        let block = self.prepare();
        if block.order.starts && block.order.ends {
            return Ok(block.fill(chrom, start, values));
        }
        let end = start as u64 + values.len() as u64;
        let range = block.candidates(start, end);
        let threshold = (values.len() as u64)
            .saturating_mul(8)
            .saturating_add((range.len() as u64).saturating_mul(8));
        let mut work = 0u64;
        for i in range {
            // SAFETY: candidates bounds i to the validated immutable block.
            let (s, e, _) = unsafe { block.item_unchecked(i) };
            work += (e as u64).min(end).saturating_sub(s.max(start) as u64);
            if work > threshold {
                let need = values
                    .len()
                    .checked_add(1)
                    .and_then(|n| n.checked_mul(3))
                    .ok_or_else(|| invalid("wig fill scratch", "slot count overflows usize"))?;
                if scratch.len() < need {
                    scratch
                        .try_reserve(need - scratch.len())
                        .map_err(|e| invalid("wig fill scratch", e.to_string()))?;
                    scratch.resize(need, 0);
                }
                return block.fill_with_scratch(chrom, start, values, scratch);
            }
        }
        Ok(block.fill(chrom, start, values))
    }

    fn fill_overlapping_fixed(&self, start: u32, values: &mut [f32]) -> usize {
        let n = self.len() as u64;
        let origin = self.start as u64;
        let step = self.step as u64;
        let span = self.span as u64;
        let lo = (start as u64).max(origin);
        let hi = (start as u64 + values.len() as u64).min(origin + (n - 1) * step + span);
        if lo >= hi {
            return 0;
        }
        let first = (lo - origin).checked_div(step).unwrap_or(n - 1).min(n - 1);
        let last = (hi - 1 - origin)
            .checked_div(step)
            .unwrap_or(n - 1)
            .min(n - 1);
        for i in first..=last {
            let s = origin + i * step;
            // The next item overwrites this item's suffix; equal spans make
            // ends monotone, so each output position only needs one write.
            let e = if i + 1 == n { s + span } else { s + step };
            let dst_lo = (s.max(lo) - start as u64) as usize;
            let dst_hi = (e.min(hi) - start as u64) as usize;
            // SAFETY: first..=last is within the parsed fixedStep count;
            // parse validated the complete four-byte value of each item.
            let v = unsafe { f32_unchecked(self.buf, WIG_HEADER_LEN + i as usize * 4) };
            // SAFETY: These effective intervals partition [lo,hi), clipped
            // to the caller's output. In particular dst_lo < dst_hi <= len.
            unsafe { values.get_unchecked_mut(dst_lo..dst_hi) }.fill(v);
        }
        // Preserve the API's logical write count, including overwritten bases,
        // without visiting every overlapping item. Integrate the arithmetic
        // progression of interval starts using triangular sums.
        let writes = if step == 0 {
            (hi - lo) * n
        } else {
            let ramp = |x: u64| {
                if x == 0 {
                    return 0;
                }
                let k = ((x - 1) / step + 1).min(n);
                k * x - step * k * (k - 1) / 2
            };
            let prefix =
                |t: u64| ramp(t.saturating_sub(origin)) - ramp(t.saturating_sub(origin + span));
            prefix(hi) - prefix(lo)
        };
        usize::try_from(writes).expect("number of writes overflows usize")
    }

    /// Fills arbitrary overlapping items without repeatedly writing a base.
    /// `scratch` needs 3*(values.len()+1) slots and is reset by this call.
    /// Cost is O((N + L) * alpha(L)) with union by size and path compression;
    /// prefer fill for sparse or non-overlapping input. Later items still win,
    /// and the return value includes logically overwritten positions.
    pub fn fill_with_scratch(
        &self,
        chrom: ChromId,
        start: u32,
        values: &mut [f32],
        scratch: &mut [usize],
    ) -> Result<usize> {
        if chrom != self.chrom || values.is_empty() {
            return Ok(0);
        }
        let slots = values.len() + 1;
        let need = slots
            .checked_mul(3)
            .ok_or_else(|| invalid("wig fill scratch", "slot count overflows usize"))?;
        if scratch.len() < need {
            return Err(invalid("wig fill scratch", "need 3*(values.len()+1) slots"));
        }
        let (parents, rest) = scratch[..need].split_at_mut(slots);
        let (sizes, ends) = rest.split_at_mut(slots);
        for i in 0..slots {
            parents[i] = i;
            sizes[i] = 1;
            ends[i] = i;
        }
        let end = start as u64 + values.len() as u64;
        let mut written = 0u64;
        for i in self.candidates(start, end).rev() {
            // SAFETY: candidates is bounded by the parsed item count.
            let (s, e, v) = unsafe { self.item_unchecked(i) };
            let lo = s.max(start) as u64;
            let hi = (e as u64).min(end);
            if lo >= hi {
                continue;
            }
            let (lo, hi) = ((lo - start as u64) as usize, (hi - start as u64) as usize);
            written += (hi - lo) as u64;
            // SAFETY: parents is an identity forest over 0..=values.len();
            // unions below always link different roots within that range.
            let mut root = unsafe { find_root(parents, lo) };
            let mut at = ends[root];
            while at < hi {
                // SAFETY: at < hi <= values.len().
                unsafe {
                    *values.get_unchecked_mut(at) = v;
                }
                // SAFETY: at+1 <= values.len(); all parents remain in range
                // and unions of different roots preserve the acyclic forest.
                let mut other = unsafe { find_root(parents, at + 1) };
                // Components consist of consecutive painted bases followed
                // by one unpainted base. Join adjacent components after painting.
                if sizes[root] < sizes[other] {
                    std::mem::swap(&mut root, &mut other);
                }
                parents[other] = root;
                sizes[root] += sizes[other];
                ends[root] = ends[root].max(ends[other]);
                at = ends[root];
            }
        }
        usize::try_from(written)
            .map_err(|_| invalid("wig fill", "number of writes overflows usize"))
    }

    #[inline]
    fn fill_selected<const KIND: u8>(
        &self,
        start: u32,
        values: &mut [f32],
        range: Range<usize>,
    ) -> usize {
        if self.order.starts && self.order.ends && self.order.overlapping {
            self.fill_kind::<KIND, true>(start, values, range.start, range.end)
        } else {
            self.fill_kind::<KIND, false>(start, values, range.start, range.end)
        }
    }

    // KIND is selected by fill; first and last are bounded by self.len().
    // MERGE requires certified monotone starts and ends.
    #[inline]
    fn fill_kind<const KIND: u8, const MERGE: bool>(
        &self,
        start: u32,
        values: &mut [f32],
        first: usize,
        last: usize,
    ) -> usize {
        let end = start as u64 + values.len() as u64;
        // A section has at most u16::MAX intervals of at most u32::MAX bases.
        let mut written = 0u64;
        for i in first..last {
            // SAFETY: fill selects KIND for the parsed layout and bounds
            // the range to the item count validated by parse.
            let (s, e, v) = unsafe { self.item_of_kind::<KIND>(i) };
            let lo = s.max(start) as u64;
            let hi = (e as u64).min(end);
            if lo >= hi {
                continue;
            }
            written += hi - lo;
            let hi = if MERGE && i + 1 < self.len() {
                // SAFETY: i+1 < len and KIND is the parsed layout.
                let (next, _, _) = unsafe { self.item_of_kind::<KIND>(i + 1) };
                hi.min(next as u64)
            } else {
                hi
            };
            if lo >= hi {
                continue;
            }
            let (lo, hi) = ((lo - start as u64) as usize, (hi - start as u64) as usize);
            // SAFETY: lo < hi after clipping; lo >= start and
            // hi <= start + values.len(), so this is within values.
            unsafe { values.get_unchecked_mut(lo..hi) }.fill(v);
        }
        usize::try_from(written).expect("number of writes overflows usize")
    }
}

#[derive(Debug, Clone)]
struct WigIntervals<'a> {
    block: WigBlock<'a>,
    range: Range<usize>,
}

impl WigIntervals<'_> {
    fn fold_kind<const KIND: u8, B, F>(self, mut value: B, mut f: F) -> B
    where
        F: FnMut(B, (u32, u32, f32)) -> B,
    {
        for i in self.range {
            // SAFETY: Construction bounds the range to the parsed count;
            // fold dispatches KIND once to match the immutable block layout.
            value = f(value, unsafe { self.block.item_of_kind::<KIND>(i) });
        }
        value
    }
}

impl Iterator for WigIntervals<'_> {
    type Item = (u32, u32, f32);
    fn next(&mut self) -> Option<Self::Item> {
        self.nth(0)
    }
    fn nth(&mut self, n: usize) -> Option<Self::Item> {
        let i = self.range.nth(n)?;
        // SAFETY: All range mutations preserve its bounds within block.len().
        Some(unsafe { self.block.item_unchecked(i) })
    }
    fn size_hint(&self) -> (usize, Option<usize>) {
        self.range.size_hint()
    }
    fn count(self) -> usize {
        self.range.len()
    }
    fn last(mut self) -> Option<Self::Item> {
        self.next_back()
    }
    fn fold<B, F>(self, init: B, f: F) -> B
    where
        F: FnMut(B, Self::Item) -> B,
    {
        match self.block.kind {
            WigKind::BedGraph => self.fold_kind::<1, B, F>(init, f),
            WigKind::VarStep => self.fold_kind::<2, B, F>(init, f),
            WigKind::FixedStep => self.fold_kind::<3, B, F>(init, f),
        }
    }
}

impl DoubleEndedIterator for WigIntervals<'_> {
    fn next_back(&mut self) -> Option<Self::Item> {
        self.nth_back(0)
    }
    fn nth_back(&mut self, n: usize) -> Option<Self::Item> {
        let i = self.range.nth_back(n)?;
        // SAFETY: Both ends of the range stay within the parsed item count.
        Some(unsafe { self.block.item_unchecked(i) })
    }
}

impl ExactSizeIterator for WigIntervals<'_> {
    fn len(&self) -> usize {
        self.range.len()
    }
}
impl FusedIterator for WigIntervals<'_> {}

/// # Safety
/// at and every parent must index parents, which must be an acyclic forest
/// whose roots point to themselves.
#[inline]
unsafe fn find_root(parents: &mut [usize], mut at: usize) -> usize {
    // SAFETY: The caller maintains an acyclic forest with in-range parents;
    // path halving preserves both properties, including root self-links.
    unsafe {
        loop {
            let parent = *parents.get_unchecked(at);
            if parent == at {
                return at;
            }
            let grandparent = *parents.get_unchecked(parent);
            *parents.get_unchecked_mut(at) = grandparent;
            at = grandparent;
        }
    }
}

/// One bigBed record. `rest` holds the tab-separated fields after the first
/// three, without the NUL terminator.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct BedEntry<'a> {
    pub chrom: ChromId,
    pub start: u32,
    pub end: u32,
    pub rest: &'a [u8],
}

impl BedEntry<'_> {
    pub fn overlaps(&self, chrom: ChromId, start: u32, end: u32) -> bool {
        self.chrom == chrom && overlaps(self.start, self.end, start, end)
    }
}

/// One bigBed data block: variable-length records, sorted by (chrom, start),
/// possibly spanning chromosomes.
#[derive(Debug, Clone, Copy)]
pub struct BedBlock<'a> {
    buf: &'a [u8],
}

impl<'a> BedBlock<'a> {
    /// Wraps decompressed block bytes; records are validated as they are
    /// read.
    pub fn parse(buf: &'a [u8]) -> Result<BedBlock<'a>> {
        Ok(BedBlock { buf })
    }

    pub fn entries(&self) -> BedEntries<'a> {
        BedEntries {
            buf: self.buf,
            at: 0,
        }
    }

    /// Validates the records, records their offsets and certifies their
    /// order in one pass over the encoded bytes. `offsets` is cleared first;
    /// reusing one vector across blocks avoids allocating after warm-up. The
    /// returned view borrows the offsets and the original bytes; item access
    /// no longer scans NUL-terminated fields.
    pub fn prepare<'w>(&self, offsets: &'w mut Vec<u32>) -> Result<PreparedBedBlock<'w>>
    where
        'a: 'w,
    {
        if self.buf.len() > u32::MAX as usize {
            return Err(invalid("bigBed block", "longer than 4 GiB"));
        }
        offsets.clear();
        let mut order = OrderScan::new();
        let mut entries = self.entries();
        while entries.at < self.buf.len() {
            let at = entries.at;
            let entry = entries.next().expect("nonempty remainder")?;
            offsets.push(at as u32);
            order.push((entry.chrom, entry.start), (entry.chrom, entry.end));
        }
        Ok(PreparedBedBlock {
            buf: self.buf,
            offsets,
            order: order.finish(),
        })
    }
}

/// A validated bigBed block with record offsets in caller-owned storage.
#[derive(Debug, Clone, Copy)]
pub struct PreparedBedBlock<'a> {
    buf: &'a [u8],
    offsets: &'a [u32],
    order: Order,
}

impl<'a> PreparedBedBlock<'a> {
    pub fn len(&self) -> usize {
        self.offsets.len()
    }
    pub fn is_empty(&self) -> bool {
        self.offsets.is_empty()
    }

    #[inline]
    fn bounds(&self, i: usize) -> (Point, Point) {
        let at = self.offsets[i] as usize;
        // SAFETY: prepare validated every stored offset and its 12-byte
        // header. The caller-owned offsets are immutably borrowed by this view.
        unsafe {
            let chrom = u32_unchecked(self.buf, at);
            (
                (chrom, u32_unchecked(self.buf, at + 4)),
                (chrom, u32_unchecked(self.buf, at + 8)),
            )
        }
    }

    /// O(1) record access, including its borrowed text.
    ///
    /// # Panics
    /// Panics when i >= self.len().
    #[inline]
    pub fn item(&self, i: usize) -> BedEntry<'a> {
        assert!(i < self.len(), "bigBed record index out of bounds");
        // SAFETY: The assertion bounds i to offsets validated by prepare.
        unsafe { self.item_unchecked(i) }
    }

    /// # Safety
    /// i < self.len().
    #[inline]
    unsafe fn item_unchecked(&self, i: usize) -> BedEntry<'a> {
        // SAFETY: i is bounded by the caller. prepare validated a complete
        // record between each offset and the next (or buf.len()), ending in
        // one NUL byte. Both input and offset storage remain immutable.
        unsafe {
            let at = *self.offsets.get_unchecked(i) as usize;
            let end = if i + 1 < self.len() {
                *self.offsets.get_unchecked(i + 1) as usize
            } else {
                self.buf.len()
            };
            BedEntry {
                chrom: u32_unchecked(self.buf, at),
                start: u32_unchecked(self.buf, at + 4),
                end: u32_unchecked(self.buf, at + 8),
                rest: self.buf.get_unchecked(at + BED_FIXED_LEN..end - 1),
            }
        }
    }

    pub fn entries(
        &self,
    ) -> impl ExactSizeIterator<Item = BedEntry<'a>> + DoubleEndedIterator + FusedIterator + use<'a>
    {
        let block = *self;
        (0..self.len()).map(move |i| {
            // SAFETY: The range is bounded by the validated record count.
            unsafe { block.item_unchecked(i) }
        })
    }

    /// Matching entries in stored order. Sorted starts and ends give
    /// O(log N + K) access; overlapping ends use a conservative max-span bound.
    /// Unsorted records use a linear scan, preserving their original order.
    pub fn overlapping(
        &self,
        chrom: ChromId,
        start: u32,
        end: u32,
    ) -> impl DoubleEndedIterator<Item = BedEntry<'a>> + FusedIterator + use<'a> {
        let range = self
            .order
            .candidates(self.len(), (chrom, start), (chrom, end), |i| self.bounds(i));
        let block = *self;
        range.filter_map(move |i| {
            let (s, e) = block.bounds(i);
            if s.0 == chrom && overlaps(s.1, e.1, start, end) {
                // SAFETY: candidates bounds i to the validated record count.
                Some(unsafe { block.item_unchecked(i) })
            } else {
                None
            }
        })
    }
}

/// Iterator over [`BedBlock::entries`]. Yields an error once and stops on a
/// malformed record.
#[derive(Debug, Clone)]
pub struct BedEntries<'a> {
    buf: &'a [u8],
    at: usize,
}

impl<'a> Iterator for BedEntries<'a> {
    type Item = Result<BedEntry<'a>>;

    fn size_hint(&self) -> (usize, Option<usize>) {
        let remaining = self.buf.len() - self.at;
        (
            0,
            Some(if remaining == 0 {
                0
            } else {
                remaining / 13 + 1
            }),
        )
    }

    fn next(&mut self) -> Option<Self::Item> {
        let buf = self.buf;
        if self.at >= buf.len() {
            return None;
        }
        let at = self.at;
        if buf.len() - at < BED_FIXED_LEN {
            self.at = buf.len();
            return Some(Err(Error::Truncated {
                what: "bigBed record",
                need: BED_FIXED_LEN as u64,
                have: (buf.len() - at) as u64,
            }));
        }
        let text = at + BED_FIXED_LEN;
        let Some(nul) = bed_text_len(&buf[text..]) else {
            self.at = buf.len();
            return Some(Err(invalid("bigBed record", "missing NUL terminator")));
        };
        self.at = text + nul + 1;
        // SAFETY: The remaining-length check validates the 12-byte header;
        // position found a NUL inside buf[text..], bounding the text slice.
        Some(Ok(unsafe {
            BedEntry {
                chrom: u32_unchecked(buf, at),
                start: u32_unchecked(buf, at + 4),
                end: u32_unchecked(buf, at + 8),
                rest: buf.get_unchecked(text..text + nul),
            }
        }))
    }
}

impl FusedIterator for BedEntries<'_> {}

#[cfg(test)]
mod fill_tests {
    use super::*;
    use proptest::prelude::*;

    proptest! {
        #[test]
        fn adaptive_fill_matches_stored_order(
            records in prop::collection::vec((0u32..128, 0u32..1024, any::<u32>()), 64..256),
            start in 0u32..128, len in 0usize..256, sorted in any::<bool>(),
        ) {
            let mut records: Vec<_> = records.into_iter().map(|(s, width, v)| (s, s + width, v)).collect();
            if sorted { records.sort_by_key(|e| e.0); }
            let mut buf = Vec::new();
            for n in [0u32, 0, 2048, 0, 0] { buf.extend(n.to_le_bytes()); }
            buf.extend([1, 0]);
            buf.extend((records.len() as u16).to_le_bytes());
            for &(s, e, v) in &records {
                buf.extend(s.to_le_bytes());
                buf.extend(e.to_le_bytes());
                buf.extend(v.to_le_bytes());
            }
            let block = WigBlock::parse(&buf).unwrap();
            let mut actual = vec![f32::from_bits(0xDEAD_BEEF); len + 2];
            let mut expected = actual.clone();
            let mut count = 0;
            for &(s, e, bits) in &records {
                for (i, value) in expected[1..len + 1].iter_mut().enumerate() {
                    let pos = start as u64 + i as u64;
                    if s as u64 <= pos && pos < e as u64 {
                        *value = f32::from_bits(bits);
                        count += 1;
                    }
                }
            }
            let mut scratch = Vec::new();
            let got = block.fill_reusing(0, start, &mut actual[1..len + 1], &mut scratch).unwrap();
            prop_assert_eq!(got, count);
            prop_assert!(actual.iter().zip(&expected).all(|(a, b)| a.to_bits() == b.to_bits()));
            // Reuse capacity with a different destination length and window.
            let mut small = [0.0; 17];
            let mut reference = small;
            for &(s, e, bits) in &records {
                for (i, value) in reference.iter_mut().enumerate() {
                    let pos = start as u64 + 19 + i as u64;
                    if s as u64 <= pos && pos < e as u64 {
                        *value = f32::from_bits(bits);
                    }
                }
            }
            block.fill_reusing(0, start + 19, &mut small, &mut scratch).unwrap();
            prop_assert!(small.iter().zip(reference).all(|(a, b)| a.to_bits() == b.to_bits()));
        }
    }
}

/// One zoom-level summary over `[start, end)` of `chrom`.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Summary {
    pub chrom: ChromId,
    pub start: u32,
    pub end: u32,
    /// Covered bases in the window.
    pub valid: u32,
    pub min: f32,
    pub max: f32,
    pub sum: f32,
    pub sum_squares: f32,
}

/// One zoom-level data block: fixed-size summary records.
#[derive(Debug, Clone, Copy)]
pub struct ZoomBlock<'a> {
    buf: &'a [u8],
    order: Order,
}

impl<'a> ZoomBlock<'a> {
    pub fn parse(buf: &'a [u8]) -> Result<ZoomBlock<'a>> {
        if buf.len() % ZOOM_RECORD_LEN != 0 {
            return Err(invalid(
                "zoom block",
                format!("{} bytes is not a multiple of {ZOOM_RECORD_LEN}", buf.len()),
            ));
        }
        Ok(ZoomBlock {
            buf,
            order: Order::UNKNOWN,
        })
    }

    /// Certifies record ordering in O(N) time for repeated range queries.
    /// No data is copied or allocated; unordered records remain queryable.
    pub fn prepare(mut self) -> Self {
        self.order = Order::scan(self.len(), |i| self.bounds(i));
        self
    }

    #[inline]
    fn bounds(&self, i: usize) -> (Point, Point) {
        debug_assert!(i < self.len());
        let at = i * ZOOM_RECORD_LEN;
        // SAFETY: Internal callers bound i to the complete records validated
        // by parse; all three coordinate fields are inside the record.
        unsafe {
            let chrom = u32_unchecked(self.buf, at);
            (
                (chrom, u32_unchecked(self.buf, at + 4)),
                (chrom, u32_unchecked(self.buf, at + 8)),
            )
        }
    }

    /// Matching records in stored order. A prepared, ordered block uses
    /// binary search before decoding the matching records' numeric fields.
    pub fn overlapping(
        &self,
        chrom: ChromId,
        start: u32,
        end: u32,
    ) -> impl DoubleEndedIterator<Item = Summary> + FusedIterator + use<'a> {
        let range = self
            .order
            .candidates(self.len(), (chrom, start), (chrom, end), |i| self.bounds(i));
        let block = *self;
        range.filter_map(move |i| {
            let (s, e) = block.bounds(i);
            if s.0 == chrom && overlaps(s.1, e.1, start, end) {
                // SAFETY: candidates bounds i to the parsed record count.
                Some(unsafe { block.record_unchecked(i) })
            } else {
                None
            }
        })
    }

    pub fn len(&self) -> usize {
        self.buf.len() / ZOOM_RECORD_LEN
    }

    pub fn is_empty(&self) -> bool {
        self.buf.is_empty()
    }

    /// # Panics
    /// Panics if `i >= self.len()`.
    #[inline]
    pub fn record(&self, i: usize) -> Summary {
        assert!(i < self.len(), "zoom record index out of bounds");
        // SAFETY: The assertion bounds i to the complete records validated
        // by parse. The buffer cannot change through this shared view.
        unsafe { self.record_unchecked(i) }
    }

    /// # Safety
    /// `i < self.len()`.
    #[inline]
    unsafe fn record_unchecked(&self, i: usize) -> Summary {
        let at = i * ZOOM_RECORD_LEN;
        let buf = self.buf;
        // SAFETY: The caller bounds i to the complete 32-byte records
        // validated by parse. All fields lie inside that record.
        unsafe {
            Summary {
                chrom: u32_unchecked(buf, at),
                start: u32_unchecked(buf, at + 4),
                end: u32_unchecked(buf, at + 8),
                valid: u32_unchecked(buf, at + 12),
                min: f32_unchecked(buf, at + 16),
                max: f32_unchecked(buf, at + 20),
                sum: f32_unchecked(buf, at + 24),
                sum_squares: f32_unchecked(buf, at + 28),
            }
        }
    }

    pub fn records(
        &self,
    ) -> impl ExactSizeIterator<Item = Summary> + DoubleEndedIterator + FusedIterator + use<'a>
    {
        let block = *self;
        (0..block.len()).map(move |i| {
            // SAFETY: The range bounds i to the complete record count.
            unsafe { block.record_unchecked(i) }
        })
    }
}
