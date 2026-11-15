//! [`IntervalIndex`]: per-seqid sorted interval arrays with an implicit
//! augmented interval tree (the cgranges layout), and [`SweepCursor`] for
//! queries that move along a chromosome.
//!
//! Items of a seqid are stored sorted by (start, end, id) in four parallel
//! arrays: start, end, id, and `max`, the largest end in the implicit
//! binary subtree rooted at each position (node `i` sits at the level of
//! its number of trailing one bits). The tree is the array itself, so
//! building is one sort and one bottom-up pass, and queries touch only
//! contiguous memory.
//!
//! Intervals are 0-based and half-open. An item `[s, e)` overlaps a query
//! `[qs, qe)` when `s < qe && qs < e`; an empty query overlaps nothing, and
//! an empty item overlaps a query that strictly contains its position.

use std::ops::Range;

/// Overlap search over intervals grouped by seqid. Ids are what the
/// builder was given (row numbers of the table indexed).
#[derive(Debug, Clone, Default)]
pub struct IntervalIndex {
    offsets: Vec<u32>,
    start: Vec<u32>,
    end: Vec<u32>,
    max: Vec<u32>,
    id: Vec<u32>,
    /// Level of each seqid's root, `-1` for none.
    root: Vec<i8>,
    /// Longest item of each seqid, for bounded nearest-neighbour scans.
    max_len: Vec<u32>,
}

/// Linear scans replace tree descent below subtrees of this level.
const SMALL: i32 = 3;

impl IntervalIndex {
    /// Indexes rows `0..seqid.len()`: row `i` is `seqid[i]`,
    /// `[start[i], end[i])`. Seqids are below `n_seqids`.
    pub fn build(n_seqids: usize, seqid: &[u32], start: &[u32], end: &[u32]) -> IntervalIndex {
        assert_eq!(
            seqid.len(),
            start.len(),
            "interval columns differ in length"
        );
        assert_eq!(seqid.len(), end.len(), "interval columns differ in length");
        assert!(
            seqid.len() <= u32::MAX as usize,
            "too many indexed intervals"
        );
        IntervalIndex::build_with(
            n_seqids,
            (0..seqid.len() as u32).map(|i| {
                let r = i as usize;
                (seqid[r], start[r], end[r], i)
            }),
        )
    }

    /// Indexes `(seqid, start, end, id)` items.
    pub fn build_with(
        n_seqids: usize,
        items: impl IntoIterator<Item = (u32, u32, u32, u32)>,
    ) -> IntervalIndex {
        let mut items: Vec<(u32, u32, u32, u32)> = items.into_iter().collect();
        assert!(
            items.len() <= u32::MAX as usize,
            "too many indexed intervals"
        );
        assert!(n_seqids <= u32::MAX as usize, "too many indexed seqids");
        items.sort_unstable();
        let n = items.len();
        let mut ix = IntervalIndex {
            offsets: Vec::with_capacity(n_seqids + 1),
            start: Vec::with_capacity(n),
            end: Vec::with_capacity(n),
            max: vec![0; n],
            id: Vec::with_capacity(n),
            root: vec![-1; n_seqids],
            max_len: vec![0; n_seqids],
        };
        ix.offsets.push(0);
        let mut at = 0;
        for s in 0..n_seqids as u32 {
            while at < n && items[at].0 == s {
                let (_, st, en, id) = items[at];
                ix.start.push(st);
                ix.end.push(en);
                ix.id.push(id);
                assert!(st <= en, "reversed interval in interval index");
                ix.max_len[s as usize] = ix.max_len[s as usize].max(en - st);
                at += 1;
            }
            ix.offsets.push(at as u32);
        }
        assert!(at == n, "seqid out of range in interval index");
        for s in 0..n_seqids {
            let r = ix.range(s as u32);
            ix.root[s] = index_tree(&ix.start[r.clone()], &ix.end[r.clone()], &mut ix.max[r]);
        }
        ix
    }

    /// Items indexed.
    pub fn len(&self) -> usize {
        self.start.len()
    }

    pub fn is_empty(&self) -> bool {
        self.start.is_empty()
    }

    /// Seqids the index was built for.
    pub fn seqids(&self) -> usize {
        self.root.len()
    }

    fn range(&self, seqid: u32) -> Range<usize> {
        let s = seqid as usize;
        if s >= self.root.len() {
            return 0..0;
        }
        self.offsets[s] as usize..self.offsets[s + 1] as usize
    }

    /// The items of `seqid` sorted by start: (starts, ends, ids).
    pub fn items(&self, seqid: u32) -> (&[u32], &[u32], &[u32]) {
        let r = self.range(seqid);
        (&self.start[r.clone()], &self.end[r.clone()], &self.id[r])
    }

    /// Calls `f` with the array position of every item overlapping
    /// `[start, end)`, in ascending position (start order).
    fn each_position(&self, seqid: u32, start: u32, end: u32, mut f: impl FnMut(usize)) {
        let r = self.range(seqid);
        let n = r.len();
        if n == 0 || start >= end {
            return;
        }
        let off = r.start;
        let s = &self.start[r.clone()];
        let e = &self.end[r.clone()];
        let m = &self.max[r];
        let k = self.root[seqid as usize] as i32;
        let mut stack = [(0usize, 0i32, false); 64];
        let mut t = 1;
        stack[0] = ((1usize << k) - 1, k, false);
        while t > 0 {
            t -= 1;
            let (x, k, w) = stack[t];
            if k <= SMALL {
                let i0 = x >> k << k;
                let i1 = (i0 + (1 << (k + 1)) - 1).min(n);
                let mut i = i0;
                while i < i1 && s[i] < end {
                    if start < e[i] {
                        f(off + i);
                    }
                    i += 1;
                }
            } else if !w {
                let y = x - (1 << (k - 1));
                stack[t] = (x, k, true);
                t += 1;
                if y >= n || m[y] > start {
                    stack[t] = (y, k - 1, false);
                    t += 1;
                }
            } else if x < n && s[x] < end {
                if start < e[x] {
                    f(off + x);
                }
                stack[t] = (x + (1 << (k - 1)), k - 1, false);
                t += 1;
            }
        }
    }

    /// Calls `f` with the id of every item overlapping `[start, end)` on
    /// `seqid`, in start order. Unknown seqids have no items.
    #[inline]
    pub fn for_each_overlap(&self, seqid: u32, start: u32, end: u32, mut f: impl FnMut(u32)) {
        self.each_position(seqid, start, end, |p| f(self.id[p]));
    }

    /// Appends the ids of the items overlapping `[start, end)`, in start
    /// order.
    pub fn overlaps(&self, seqid: u32, start: u32, end: u32, out: &mut Vec<u32>) {
        self.for_each_overlap(seqid, start, end, |id| out.push(id));
    }

    /// Items overlapping `[start, end)`.
    pub fn count_overlaps(&self, seqid: u32, start: u32, end: u32) -> usize {
        let mut n = 0;
        self.each_position(seqid, start, end, |_| n += 1);
        n
    }

    /// Appends the ids of the items lying within `[start, end)`.
    pub fn contained_in(&self, seqid: u32, start: u32, end: u32, out: &mut Vec<u32>) {
        self.each_position(seqid, start, end, |p| {
            if self.start[p] >= start && self.end[p] <= end {
                out.push(self.id[p]);
            }
        });
    }

    /// Appends the ids of the items covering all of `[start, end)` (for an
    /// empty query, the items covering position `start`).
    pub fn containing(&self, seqid: u32, start: u32, end: u32, out: &mut Vec<u32>) {
        let qe = end.max(start.saturating_add(1));
        self.each_position(seqid, start, qe, |p| {
            if self.start[p] <= start && self.end[p] >= qe {
                out.push(self.id[p]);
            }
        });
    }

    /// Appends the ids of the items nearest to `[start, end)` and returns
    /// their distance: `0` for overlapping items (an empty query counts as
    /// position `start`), else the number of bases between the query and
    /// the closest items on either side, all ties included. `None` when
    /// the seqid has no items.
    pub fn nearest(&self, seqid: u32, start: u32, end: u32, out: &mut Vec<u32>) -> Option<u32> {
        let r = self.range(seqid);
        if r.is_empty() {
            return None;
        }
        let qe = end.max(start.saturating_add(1));
        let before = out.len();
        self.overlaps(seqid, start, qe, out);
        if out.len() > before {
            return Some(0);
        }
        let s = &self.start[r.clone()];
        let ids = &self.id[r.clone()];
        let e = &self.end[r.clone()];
        // With no overlap, every item starting before the query end ends
        // at or before the query start.
        let p = s.partition_point(|&x| x < qe);
        let down = (p < s.len()).then(|| s[p] - qe);
        let (up_end, local_start) = if p == 0 {
            (None, None)
        } else {
            let (end, local_start) = self.nearest_left(seqid, p);
            (Some(end), local_start)
        };
        let up = up_end.map(|b| start - b);
        let best = match (up, down) {
            (Some(u), Some(d)) => u.min(d),
            (Some(u), None) => u,
            (None, Some(d)) => d,
            (None, None) => return None,
        };
        if up == Some(best) {
            let b = up_end.unwrap();
            if let Some(local_start) = local_start {
                for i in local_start..p {
                    if e[i] == b {
                        out.push(ids[i]);
                    }
                }
            } else {
                self.each_ending_at(seqid, p, b, |i| out.push(self.id[i]));
            }
        }
        if down == Some(best) {
            let mut i = p;
            while i < s.len() && s[i] == s[p] {
                out.push(ids[i]);
                i += 1;
            }
        }
        Some(best)
    }

    /// A nearby short interval often proves the maximum immediately.
    /// Bound that scan so long intervals cannot turn a query linear.
    fn nearest_left(&self, seqid: u32, prefix: usize) -> (u32, Option<usize>) {
        let r = self.range(seqid);
        let starts = &self.start[r.clone()];
        let ends = &self.end[r];
        let max_len = self.max_len[seqid as usize];
        let mut best = 0;
        let mut i = prefix;
        for _ in 0..8 {
            if i == 0 {
                return (best, Some(0));
            }
            i -= 1;
            if starts[i].saturating_add(max_len) < best {
                return (best, Some(i + 1));
            }
            best = best.max(ends[i]);
        }
        if i == 0 {
            (best, Some(0))
        } else {
            (self.prefix_max_end(seqid, prefix), None)
        }
    }

    /// Maximum end in a nonempty prefix, following its boundary down the
    /// tree and consuming complete left subtrees in O(log n).
    fn prefix_max_end(&self, seqid: u32, prefix: usize) -> u32 {
        let r = self.range(seqid);
        let ends = &self.end[r.clone()];
        let max = &self.max[r];
        let mut k = self.root[seqid as usize] as u32;
        let mut x = (1usize << k) - 1;
        let mut best = 0;
        loop {
            if x < prefix {
                best = best.max(ends[x]);
                if k == 0 {
                    return best;
                }
                let step = 1usize << (k - 1);
                best = best.max(max[x - step]);
                x += step;
            } else {
                if k == 0 {
                    return best;
                }
                x -= 1usize << (k - 1);
            }
            k -= 1;
        }
    }

    /// Visits all ties at the prefix's maximum end in start order. A
    /// subtree whose maximum is smaller cannot contain a tie; virtual
    /// nodes beyond the final row descend only into their left subtree.
    fn each_ending_at(&self, seqid: u32, prefix: usize, end: u32, mut f: impl FnMut(usize)) {
        let r = self.range(seqid);
        let ends = &self.end[r.clone()];
        let max = &self.max[r.clone()];
        let k = self.root[seqid as usize] as u32;
        let mut stack = [(0usize, 0u32, false); 64];
        stack[0] = ((1usize << k) - 1, k, false);
        let mut t = 1;
        while t != 0 {
            t -= 1;
            let (x, k, visit) = stack[t];
            if x >= prefix {
                if k != 0 {
                    stack[t] = (x - (1usize << (k - 1)), k - 1, false);
                    t += 1;
                }
                continue;
            }
            if visit || k == 0 {
                if ends[x] == end {
                    f(r.start + x);
                }
                if k != 0 {
                    stack[t] = (x + (1usize << (k - 1)), k - 1, false);
                    t += 1;
                }
            } else if max[x] >= end {
                stack[t] = (x, k, true);
                stack[t + 1] = (x - (1usize << (k - 1)), k - 1, false);
                t += 2;
            }
        }
    }

    /// Overlaps of many queries `(seqid, start, end)` as one CSR result:
    /// query `q`'s ids are `ids[offsets[q]..offsets[q + 1]]`. Both buffers
    /// are cleared first.
    pub fn overlaps_batch(
        &self,
        queries: impl IntoIterator<Item = (u32, u32, u32)>,
        offsets: &mut Vec<u64>,
        ids: &mut Vec<u32>,
    ) {
        offsets.clear();
        ids.clear();
        offsets.push(0);
        for (seqid, start, end) in queries {
            self.overlaps(seqid, start, end, ids);
            offsets.push(ids.len() as u64);
        }
    }

    /// A cursor over `seqid` for queries in ascending start order.
    pub fn cursor(&self, seqid: u32) -> SweepCursor<'_> {
        let mut c = SweepCursor {
            index: self,
            hi: 0,
            next: 0,
            active: Vec::new(),
            min_end: u32::MAX,
            floor: 0,
            seqid,
        };
        c.seek(seqid);
        c
    }

    /// Heap bytes held.
    pub fn heap_bytes(&self) -> usize {
        4 * (self.offsets.capacity()
            + self.start.capacity()
            + self.end.capacity()
            + self.max.capacity()
            + self.id.capacity()
            + self.max_len.capacity())
            + self.root.capacity()
    }
}

/// Fills `max` for the implicit tree over `start` / `end` (sorted by
/// start); returns the root level, `-1` when empty.
fn index_tree(start: &[u32], end: &[u32], max: &mut [u32]) -> i8 {
    let n = start.len();
    if n == 0 {
        return -1;
    }
    debug_assert!(start.windows(2).all(|w| w[0] <= w[1]));
    let mut last_i = 0usize;
    let mut last = 0u32;
    let mut i = 0;
    while i < n {
        last_i = i;
        last = end[i];
        max[i] = end[i];
        i += 2;
    }
    let mut k = 1;
    while (1usize << k) <= n {
        let x = 1usize << (k - 1);
        let i0 = (x << 1) - 1;
        let step = x << 2;
        let mut i = i0;
        while i < n {
            let el = max[i - x];
            let er = if i + x < n { max[i + x] } else { last };
            max[i] = end[i].max(el).max(er);
            i += step;
        }
        last_i = if (last_i >> k) & 1 == 1 {
            last_i - x
        } else {
            last_i + x
        };
        if last_i < n && max[last_i] > last {
            last = max[last_i];
        }
        k += 1;
    }
    (k - 1) as i8
}

/// Overlap queries along one seqid whose starts do not decrease, as when
/// scanning a coordinate-sorted BAM.
///
/// The cursor keeps the *active* items: admitted (start below the end of
/// some query so far) and not yet expired (end above the *floor*). Each
/// item is admitted and expired once, so moving the cursor is amortized
/// O(1) per item; a query costs O(active) and allocates nothing once the
/// active buffer has grown to its peak.
///
/// - [`SweepCursor::query`]`(start, end, out)` raises the floor to `start`
///   and appends the overlapping ids.
/// - For several queries per step (the aligned blocks of one spliced read,
///   whose starts are not monotonic across reads), call
///   [`SweepCursor::advance`]`(read_start)` once, then
///   [`SweepCursor::overlaps`] for each block: any query with start at or
///   above the floor is answered exactly.
/// - A query or advance below the floor is still answered exactly: the
///   cursor reseeds itself from the tree (O(log n + k)).
///
/// Results come in start order; ids are those of the index.
#[derive(Debug, Clone)]
pub struct SweepCursor<'a> {
    index: &'a IntervalIndex,
    /// End of the seqid's array range.
    hi: usize,
    /// Next array position not yet admitted.
    next: usize,
    /// Array positions of active items, ascending.
    active: Vec<u32>,
    /// Smallest end among active items.
    min_end: u32,
    floor: u32,
    seqid: u32,
}

impl SweepCursor<'_> {
    /// Moves to the start of `seqid` (keeping the buffer allocation).
    pub fn seek(&mut self, seqid: u32) {
        let r = self.index.range(seqid);
        self.seqid = seqid;
        self.next = r.start;
        self.hi = r.end;
        self.active.clear();
        self.min_end = u32::MAX;
        self.floor = 0;
    }

    pub fn seqid(&self) -> u32 {
        self.seqid
    }

    /// The lowest start future queries may have without reseeding.
    pub fn floor(&self) -> u32 {
        self.floor
    }

    /// Items currently admitted and not expired.
    pub fn active_len(&self) -> usize {
        self.active.len()
    }

    /// Declares that later queries start at or after `floor`, expiring
    /// items that end at or before it. A lower floor than the current one
    /// reseeds the cursor.
    #[inline]
    pub fn advance(&mut self, floor: u32) {
        if floor < self.floor {
            self.reseed(floor);
            return;
        }
        self.floor = floor;
        if floor >= self.min_end {
            let end = &self.index.end;
            let mut min_end = u32::MAX;
            self.active.retain(|&p| {
                let e = end[p as usize];
                if e > floor {
                    min_end = min_end.min(e);
                    true
                } else {
                    false
                }
            });
            self.min_end = min_end;
        }
    }

    fn reseed(&mut self, floor: u32) {
        let ix = self.index;
        let r = ix.range(self.seqid);
        let s = &ix.start[r.clone()];
        self.active.clear();
        self.min_end = u32::MAX;
        self.floor = floor;
        self.next = r.start + s.partition_point(|&x| x < floor);
        let (active, min_end) = (&mut self.active, &mut self.min_end);
        let (starts, ends) = (&ix.start, &ix.end);
        ix.each_position(self.seqid, floor, floor.saturating_add(1), |p| {
            if starts[p] < floor {
                active.push(p as u32);
                *min_end = (*min_end).min(ends[p]);
            }
        });
    }

    #[inline]
    fn admit(&mut self, end: u32) {
        let ix = self.index;
        let (starts, ends) = (&ix.start, &ix.end);
        while self.next < self.hi && starts[self.next] < end {
            let e = ends[self.next];
            if e > self.floor {
                self.active.push(self.next as u32);
                self.min_end = self.min_end.min(e);
            }
            self.next += 1;
        }
    }

    /// Calls `f` with the id of every item overlapping `[start, end)`, in
    /// start order, without moving the floor (unless `start` is below it).
    #[inline]
    pub fn for_each_overlap(&mut self, start: u32, end: u32, mut f: impl FnMut(u32)) {
        if start < self.floor {
            self.reseed(start);
        }
        if start >= end {
            return;
        }
        self.admit(end);
        let ix = self.index;
        for &p in &self.active {
            let p = p as usize;
            // SAFETY: active positions come only from `admit`'s bounded
            // range or `each_position` in `reseed`. The index owns equal
            // length columns, exposes no mutable columns, and remains
            // borrowed for the cursor's entire lifetime.
            let (s, e, id) = unsafe {
                (
                    *ix.start.get_unchecked(p),
                    *ix.end.get_unchecked(p),
                    *ix.id.get_unchecked(p),
                )
            };
            if s < end && e > start {
                f(id);
            }
        }
    }

    /// Appends the ids overlapping `[start, end)` without moving the floor.
    #[inline]
    pub fn overlaps(&mut self, start: u32, end: u32, out: &mut Vec<u32>) {
        self.for_each_overlap(start, end, |id| out.push(id));
    }

    /// [`SweepCursor::advance`] to `start`, then appends the ids
    /// overlapping `[start, end)`.
    #[inline]
    pub fn query(&mut self, start: u32, end: u32, out: &mut Vec<u32>) {
        self.advance(start);
        self.overlaps(start, end, out);
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn small_tree_matches_scan() {
        let items: Vec<(u32, u32)> = (0..300u32)
            .map(|i| (i * 7 % 101, i * 7 % 101 + i % 37 + 1))
            .collect();
        let seq = vec![0u32; items.len()];
        let st: Vec<u32> = items.iter().map(|x| x.0).collect();
        let en: Vec<u32> = items.iter().map(|x| x.1).collect();
        let ix = IntervalIndex::build(1, &seq, &st, &en);
        for qs in 0..150 {
            for ql in [0, 1, 5, 40] {
                let mut got = Vec::new();
                ix.overlaps(0, qs, qs + ql, &mut got);
                got.sort_unstable();
                let want: Vec<u32> = (0..items.len() as u32)
                    .filter(|&i| st[i as usize] < qs + ql && qs < en[i as usize] && ql > 0)
                    .collect();
                assert_eq!(got, want);
            }
        }
    }

    #[test]
    fn prefix_max_and_all_ties_match_scans_for_partial_trees() {
        for n in [
            1, 2, 3, 4, 5, 7, 8, 9, 15, 16, 17, 31, 32, 33, 63, 64, 65, 127, 128, 129, 255, 256,
            257, 1023, 1024, 1025,
        ] {
            let index =
                IntervalIndex::build_with(1, (0..n).map(|i| (0, i, i + i * 7919 % 1000, i)));
            for prefix in 1..=n as usize {
                let want = *index.end[..prefix].iter().max().unwrap();
                assert_eq!(index.prefix_max_end(0, prefix), want);
                let mut got = Vec::new();
                index.each_ending_at(0, prefix, want, |p| got.push(p));
                let positions: Vec<usize> = (0..prefix).filter(|&p| index.end[p] == want).collect();
                assert_eq!(got, positions);
            }
        }
    }
}
