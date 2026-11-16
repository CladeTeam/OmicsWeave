//! [`BinningIndex`]: the BAI / CSI index of a BGZF file of sorted records,
//! built as htslib's `hts_idx_push` / `hts_idx_finish` build it.
//!
//! Records are pushed while they are written, with [`Logical`] offsets; the
//! offsets become virtual offsets in [`BinningIndex::finish`], before the
//! steps that compare them (merging small bins, joining chunks that share a
//! block), so the index holds the same bins, chunks and linear offsets as
//! one htslib builds from the finished file. Bins are written in ascending
//! order where htslib follows its hash table; readers do not depend on it.

use omicsweave_core::exec::Inline;
use omicsweave_core::{Error, Result};
use std::collections::BTreeMap;

use crate::writer::{Blocks, Logical, Writer};

const MIN_MARKER_DIST: u64 = 0x10000;
const UNSET: u64 = u64::MAX;
const NO_BIN: u32 = u32::MAX;

/// Which index to build.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum IndexFormat {
    /// `.bai`: 16 KiB bins over 6 levels, positions below 2^29.
    Bai,
    /// `.csi` with `min_shift` and `levels`.
    Csi { min_shift: u32, levels: u32 },
}

/// Most levels whose bin numbers (and the pseudo-bin after them) fit the
/// `u32` bin field.
const MAX_LEVELS: u32 = 10;
/// Most bits a position may take: `min_shift + 3 * levels`, so that every
/// `i64` coordinate shift and the index's `maxpos` stay representable.
const MAX_BITS: u32 = 63;

fn check_params(min_shift: u32, levels: u32) -> Result<()> {
    let bits = levels.checked_mul(3).and_then(|b| b.checked_add(min_shift));
    if levels > MAX_LEVELS || bits.is_none_or(|b| b > MAX_BITS) {
        return Err(Error::invalid(
            "index",
            format!(
                "CSI min_shift {min_shift} with {levels} levels: need at most {MAX_LEVELS} levels and min_shift + 3 * levels <= {MAX_BITS}"
            ),
        ));
    }
    Ok(())
}

/// The CSI levels htslib picks (`hts_adjust_csi_settings` from 0 levels) so
/// that `max_len` fits with bins of `2^min_shift` bases: more levels up to
/// 9, then larger bins. Returns `(min_shift, levels)`; fails when no
/// supported setting holds `max_len`.
pub fn try_csi_levels(max_len: u64, min_shift: u32) -> Result<(u32, u32)> {
    const HTS_MAX_LEVELS: u32 = 9;
    check_params(min_shift, 0)?;
    let need = max_len as u128 + 256;
    let maxpos = |shift: u32, levels: u32| 1u128 << (shift + levels * 3).min(127);
    let (shift, levels) = if need <= maxpos(min_shift, HTS_MAX_LEVELS) {
        let mut levels = 0;
        while need > maxpos(min_shift, levels) {
            levels += 1;
        }
        (min_shift, levels)
    } else {
        let mut shift = min_shift;
        while need > maxpos(shift, HTS_MAX_LEVELS) && shift <= MAX_BITS {
            shift += 1;
        }
        (shift, HTS_MAX_LEVELS)
    };
    check_params(shift, levels)?;
    Ok((shift, levels))
}

/// `hts_reg2bin`: the smallest bin holding `[beg, end)`, for parameters
/// already checked.
fn reg2bin(beg: i64, end: i64, min_shift: u32, levels: u32) -> u32 {
    let end = end - 1;
    let mut s = min_shift;
    let mut t = bin_first(levels) as i64;
    let mut l = levels;
    while l > 0 {
        if beg >> s == end >> s {
            return (t + (beg >> s)) as u32;
        }
        l -= 1;
        s += 3;
        t -= 1i64 << (l * 3);
    }
    0
}

/// `hts_reg2bin` for an in-range interval, or `(-1, 0)` for a record
/// without coordinates, with its parameters and arithmetic checked.
pub fn try_reg2bin(beg: i64, end: i64, min_shift: u32, levels: u32) -> Result<u32> {
    check_params(min_shift, levels)?;
    let maxpos = 1i128 << (min_shift + levels * 3);
    if beg < -1 || end < 0 || end < beg || beg as i128 > maxpos || end as i128 > maxpos {
        return Err(Error::invalid(
            "index",
            format!("region {beg}..{end} out of range for this index"),
        ));
    }
    Ok(reg2bin(beg, end, min_shift, levels))
}

fn bin_first(level: u32) -> u32 {
    ((1u32 << (level * 3)) - 1) / 7
}

fn bin_level(bin: u32) -> u32 {
    let (mut l, mut b) = (0, bin);
    while b != 0 {
        l += 1;
        b = (b - 1) >> 3;
    }
    l
}

fn bin_parent(bin: u32) -> u32 {
    (bin - 1) >> 3
}

fn bin_bot(bin: u32, levels: u32) -> u64 {
    let l = bin_level(bin);
    ((bin - bin_first(l)) as u64) << ((levels - l) * 3)
}

#[derive(Default)]
struct Bin {
    loff: u64,
    chunks: Vec<(u64, u64)>,
}

/// Ordinary indexes use a compact direct table. CSI can represent billions
/// of bottom-level windows without serializing that table, so unusually
/// wide spans switch to first-offset runs instead of allocating by span.
enum Linear {
    Dense(Vec<u64>),
    Sparse { len: usize, runs: Vec<(usize, u64)> },
}

impl Default for Linear {
    fn default() -> Self {
        Self::Dense(Vec::new())
    }
}

impl Linear {
    const DENSE_LIMIT: usize = 1 << 16;

    fn run(runs: &mut Vec<(usize, u64)>, end: usize, offset: u64) {
        if let Some(last) = runs.last_mut()
            && last.1 == offset
        {
            last.0 = end;
        } else {
            runs.push((end, offset));
        }
    }

    fn extend(&mut self, lo: usize, hi: usize, offset: u64, csi: bool) {
        if csi
            && hi >= Self::DENSE_LIMIT
            && let Self::Dense(dense) = self
        {
            let mut runs = Vec::new();
            for (at, &value) in dense.iter().enumerate() {
                Self::run(&mut runs, at + 1, value);
            }
            *self = Self::Sparse {
                len: dense.len(),
                runs,
            };
        }
        match self {
            Self::Dense(dense) => {
                let old = dense.len();
                if old <= hi {
                    dense.resize(hi + 1, UNSET);
                    let from = lo.max(old);
                    if from <= hi {
                        dense[from..=hi].fill(offset);
                    }
                }
            }
            Self::Sparse { len, runs } => {
                if *len <= hi {
                    if lo > hi {
                        // A boundary-empty interval extends the table but
                        // covers no window and supplies no first offset.
                        Self::run(runs, hi + 1, UNSET);
                    } else {
                        if lo > *len {
                            Self::run(runs, lo, UNSET);
                        }
                        Self::run(runs, hi + 1, offset);
                    }
                    *len = hi + 1;
                }
            }
        }
    }

    fn resolve(&mut self, resolve: impl Fn(u64) -> u64) {
        match self {
            Self::Dense(dense) => {
                for offset in dense.iter_mut() {
                    if *offset != UNSET {
                        *offset = resolve(*offset);
                    }
                }
                for i in (0..dense.len().saturating_sub(1)).rev() {
                    if dense[i] == UNSET {
                        dense[i] = dense[i + 1];
                    }
                }
            }
            Self::Sparse { runs, .. } => {
                let mut next = UNSET;
                for (_, offset) in runs.iter_mut().rev() {
                    if *offset == UNSET {
                        *offset = next;
                    } else {
                        *offset = resolve(*offset);
                        next = *offset;
                    }
                }
            }
        }
    }

    fn get(&self, at: usize) -> Option<u64> {
        match self {
            Self::Dense(dense) => dense.get(at).copied(),
            Self::Sparse { len, runs } => {
                if at >= *len {
                    None
                } else {
                    Some(runs[runs.partition_point(|&(end, _)| end <= at)].1)
                }
            }
        }
    }

    fn dense(&self) -> &[u64] {
        match self {
            Self::Dense(dense) => dense,
            Self::Sparse { .. } => unreachable!("BAI retains its serialized linear table"),
        }
    }
}

#[derive(Default)]
struct RefIndex {
    present: bool,
    bins: BTreeMap<u32, Bin>,
    linear: Linear,
    /// Pseudo-bin: (first offset, end offset, mapped, unmapped).
    meta: Option<(u64, u64, u64, u64)>,
}

/// An index built while records are written; see the module docs.
pub struct BinningIndex {
    format: IndexFormat,
    min_shift: u32,
    levels: u32,
    n_bins: u32,
    refs: Vec<RefIndex>,
    n_no_coor: u64,
    save_tid: i32,
    last_tid: i32,
    save_bin: u32,
    last_bin: u32,
    save_off: u64,
    last_off: u64,
    off_beg: u64,
    last_coor: i64,
    n_mapped: u64,
    n_unmapped: u64,
    /// Set by the first failure (including invalid parameters): the index
    /// then rejects every further record and its finish.
    failed: Option<Error>,
}

impl BinningIndex {
    /// An index of `n_refs` references whose first record starts at
    /// `start` (after the header), with parameters checked.
    pub fn try_new(format: IndexFormat, n_refs: usize, start: Logical) -> Result<BinningIndex> {
        let ix = BinningIndex::new(format, n_refs, start);
        match ix.failed {
            Some(e) => Err(e),
            None => Ok(ix),
        }
    }

    /// [`BinningIndex::try_new`] whose parameter error, if any, is reported
    /// by every later [`BinningIndex::push`] and [`BinningIndex::finish`].
    pub fn new(format: IndexFormat, n_refs: usize, start: Logical) -> BinningIndex {
        let (min_shift, levels) = match format {
            IndexFormat::Bai => (14, 5),
            IndexFormat::Csi { min_shift, levels } => (min_shift, levels),
        };
        let failed = check_params(min_shift, levels).err();
        let n_bins = if failed.is_none() {
            (((1u64 << (3 * levels + 3)) - 1) / 7) as u32
        } else {
            0
        };
        BinningIndex {
            format,
            min_shift,
            levels,
            n_bins,
            refs: (0..n_refs).map(|_| RefIndex::default()).collect(),
            n_no_coor: 0,
            save_tid: -1,
            last_tid: -1,
            save_bin: NO_BIN,
            last_bin: NO_BIN,
            save_off: start.0,
            last_off: start.0,
            off_beg: start.0,
            last_coor: u32::MAX as i64,
            n_mapped: 0,
            n_unmapped: 0,
            failed,
        }
    }

    fn meta_bin(&self) -> u32 {
        self.n_bins + 1
    }

    /// `hts_idx_amend_last`: the next record starts at `at` (after a block
    /// was cut in front of it).
    pub fn amend_last(&mut self, at: Logical) {
        if at.0 < self.last_off {
            self.failed
                .get_or_insert_with(|| Error::invalid("index", "record offsets must not decrease"));
            return;
        }
        self.last_off = at.0;
    }

    fn insert(&mut self, tid: i32, bin: u32, beg: u64, end: u64) {
        self.refs[tid as usize]
            .bins
            .entry(bin)
            .or_default()
            .chunks
            .push((beg, end));
    }

    fn insert_meta(&mut self, tid: i32, off_end: u64) {
        self.refs[tid as usize].meta =
            Some((self.off_beg, off_end, self.n_mapped, self.n_unmapped));
    }

    /// `hts_idx_push`: a record of reference `tid` (any negative value: no
    /// coordinates) over `[beg, end)`, ending at `end_at`. After an error
    /// the index rejects everything that follows.
    pub fn push(
        &mut self,
        tid: i32,
        beg: i64,
        end: i64,
        end_at: Logical,
        mapped: bool,
    ) -> Result<()> {
        if let Some(e) = &self.failed {
            return Err(e.clone());
        }
        let r = self.push_checked(tid.max(-1), beg, end, end_at, mapped);
        if let Err(e) = &r {
            self.failed = Some(e.clone());
        }
        r
    }

    fn push_checked(
        &mut self,
        tid: i32,
        beg: i64,
        end: i64,
        end_at: Logical,
        mapped: bool,
    ) -> Result<()> {
        if end_at.0 < self.last_off {
            return Err(Error::invalid("index", "record offsets must not decrease"));
        }
        let (mut beg, mut end) = if tid < 0 { (-1, 0) } else { (beg, end) };
        let maxpos = 1i128 << (self.min_shift + self.levels * 3);
        if tid >= 0 && (beg as i128 > maxpos || end as i128 > maxpos) {
            return Err(Error::invalid(
                "index",
                format!(
                    "region {beg}..{end} exceeds what the index can hold ({maxpos}); use a CSI index with more levels"
                ),
            ));
        }
        if tid >= 0 && tid as usize >= self.refs.len() {
            return Err(Error::invalid(
                "index",
                format!("reference {tid} out of range"),
            ));
        }
        if self.last_tid != tid || (self.last_tid >= 0 && tid < 0) {
            if tid >= 0 && self.n_no_coor > 0 {
                return Err(Error::invalid(
                    "index",
                    "records without coordinates must come last",
                ));
            }
            if tid >= 0 && self.refs[tid as usize].present {
                return Err(Error::invalid(
                    "index",
                    "records of a reference are not contiguous (input is not sorted)",
                ));
            }
            self.last_tid = tid;
            self.last_bin = NO_BIN;
        } else if tid >= 0 && self.last_coor > beg {
            return Err(Error::invalid(
                "index",
                format!(
                    "unsorted positions on reference #{}: {} followed by {}",
                    i64::from(tid) + 1,
                    i128::from(self.last_coor) + 1,
                    i128::from(beg) + 1
                ),
            ));
        }
        if end < beg {
            return Err(Error::invalid(
                "index",
                format!("record end {end} before its start {}", i128::from(beg) + 1),
            ));
        }
        if tid >= 0 {
            let r = &mut self.refs[tid as usize];
            r.present = true;
            if beg < 0 {
                beg = 0;
            }
            if end <= 0 {
                end = 1;
            }
            let (lo, hi) = (
                (beg >> self.min_shift) as usize,
                ((end - 1) >> self.min_shift) as usize,
            );
            // Within a reference `beg` never decreases, so every window
            // below the current length that this record covers was set by
            // an earlier record; only windows past the length can be new.
            r.linear.extend(
                lo,
                hi,
                self.last_off,
                matches!(self.format, IndexFormat::Csi { .. }),
            );
        } else {
            self.n_no_coor += 1;
        }
        let bin = reg2bin(beg, end, self.min_shift, self.levels);
        if self.last_bin != bin {
            if self.save_bin != NO_BIN {
                self.insert(self.save_tid, self.save_bin, self.save_off, self.last_off);
            }
            if self.last_bin == NO_BIN && self.save_bin != NO_BIN {
                let off_end = self.last_off;
                self.insert_meta(self.save_tid, off_end);
                self.n_mapped = 0;
                self.n_unmapped = 0;
                self.off_beg = off_end;
            }
            self.save_off = self.last_off;
            self.save_bin = bin;
            self.last_bin = bin;
            self.save_tid = tid;
        }
        if mapped {
            self.n_mapped += 1;
        } else {
            self.n_unmapped += 1;
        }
        self.last_off = end_at.0;
        self.last_coor = beg;
        Ok(())
    }

    /// `hts_idx_finish` at `end` (the position after the last record), with
    /// positions resolved through `blocks`; returns the index file's bytes
    /// (BGZF-compressed for CSI, as htslib writes it).
    pub fn finish(mut self, end: Logical, blocks: &Blocks) -> Result<Vec<u8>> {
        if let Some(e) = self.failed {
            return Err(e);
        }
        if end.0 < self.last_off {
            return Err(Error::invalid("index", "end offset before the last record"));
        }
        if self.save_tid >= 0 {
            self.insert(self.save_tid, self.save_bin, self.save_off, end.0);
            self.insert_meta(self.save_tid, end.0);
        }
        let resolve = |v: u64| blocks.resolve(Logical(v));
        let (n_bins, levels) = (self.n_bins, self.levels);
        for r in &mut self.refs {
            r.linear.resolve(resolve);
            for bin in r.bins.values_mut() {
                for c in &mut bin.chunks {
                    *c = (resolve(c.0), resolve(c.1));
                }
            }
            if let Some(m) = &mut r.meta {
                m.0 = resolve(m.0);
                m.1 = resolve(m.1);
            }
            for (&key, bin) in r.bins.iter_mut() {
                let bot = bin_bot(key, levels) as usize;
                bin.loff = if key < n_bins {
                    r.linear.get(bot).unwrap_or(0)
                } else {
                    0
                };
            }
            compress_binning(r, levels, n_bins);
        }
        let mut body = Vec::new();
        match self.format {
            IndexFormat::Bai => body.extend_from_slice(b"BAI\x01"),
            IndexFormat::Csi { .. } => {
                body.extend_from_slice(b"CSI\x01");
                body.extend_from_slice(&(self.min_shift as i32).to_le_bytes());
                body.extend_from_slice(&(self.levels as i32).to_le_bytes());
                body.extend_from_slice(&0u32.to_le_bytes());
            }
        }
        let csi = matches!(self.format, IndexFormat::Csi { .. });
        body.extend_from_slice(&(self.refs.len() as i32).to_le_bytes());
        let meta_bin = self.meta_bin();
        for r in &self.refs {
            let n = r.bins.len() + r.meta.is_some() as usize;
            body.extend_from_slice(&(n as i32).to_le_bytes());
            for (&key, bin) in &r.bins {
                body.extend_from_slice(&key.to_le_bytes());
                if csi {
                    body.extend_from_slice(&bin.loff.to_le_bytes());
                }
                body.extend_from_slice(&(bin.chunks.len() as i32).to_le_bytes());
                for &(u, v) in &bin.chunks {
                    body.extend_from_slice(&u.to_le_bytes());
                    body.extend_from_slice(&v.to_le_bytes());
                }
            }
            if let Some((beg, end, mapped, unmapped)) = r.meta {
                body.extend_from_slice(&meta_bin.to_le_bytes());
                if csi {
                    body.extend_from_slice(&0u64.to_le_bytes());
                }
                body.extend_from_slice(&2i32.to_le_bytes());
                for v in [beg, end, mapped, unmapped] {
                    body.extend_from_slice(&v.to_le_bytes());
                }
            }
            if !csi {
                let linear = r.linear.dense();
                body.extend_from_slice(&(linear.len() as i32).to_le_bytes());
                for &o in linear {
                    body.extend_from_slice(&o.to_le_bytes());
                }
            }
        }
        body.extend_from_slice(&self.n_no_coor.to_le_bytes());
        if !csi {
            return Ok(body);
        }
        let mut w = Writer::new(Vec::new(), crate::DEFAULT_LEVEL, Box::new(Inline::new()))?;
        w.write(&body)?;
        Ok(w.finish()?.out)
    }
}

fn compress_binning(r: &mut RefIndex, levels: u32, n_bins: u32) {
    for l in (1..=levels).rev() {
        let lo = bin_first(l);
        let hi = (((1u64 << (3 * (l + 1))) - 1) / 7).min(n_bins as u64) as u32;
        let keys: Vec<u32> = r.bins.range(lo..hi).map(|(&k, _)| k).collect();
        for key in keys {
            let parent = bin_parent(key);
            let bin = r.bins.get_mut(&key).expect("key listed");
            if l < levels && bin.chunks.len() > 1 {
                bin.chunks.sort_unstable_by_key(|c| c.0);
            }
            let span = (bin.chunks[bin.chunks.len() - 1].1 >> 16) - (bin.chunks[0].0 >> 16);
            if span < MIN_MARKER_DIST && r.bins.contains_key(&parent) {
                let moved = r.bins.remove(&key).expect("key listed").chunks;
                r.bins
                    .get_mut(&parent)
                    .expect("parent present")
                    .chunks
                    .extend(moved);
            }
        }
    }
    if let Some(root) = r.bins.get_mut(&0) {
        root.chunks.sort_unstable_by_key(|c| c.0);
    }
    for bin in r.bins.values_mut() {
        let c = &mut bin.chunks;
        let mut m = 0;
        for l in 1..c.len() {
            if c[m].1 >> 16 >= c[l].0 >> 16 {
                if c[m].1 < c[l].1 {
                    c[m].1 = c[l].1;
                }
            } else {
                m += 1;
                c[m] = c[l];
            }
        }
        c.truncate(m + 1);
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn reg2bin_matches_htslib_layout() {
        // BAI: bins of the bottom level start at 4681.
        assert_eq!(reg2bin(0, 1, 14, 5), 4681);
        assert_eq!(reg2bin(16384, 16385, 14, 5), 4682);
        assert_eq!(reg2bin(0, 16385, 14, 5), 585);
        assert_eq!(reg2bin(0, 1 << 29, 14, 5), 0);
        // Records without coordinates land in 4680.
        assert_eq!(reg2bin(-1, 0, 14, 5), 4680);
    }

    #[test]
    fn csi_levels_cover_the_longest_reference() {
        assert_eq!(try_csi_levels(248_956_422, 14).unwrap(), (14, 5));
        assert_eq!(try_csi_levels(1 << 31, 14).unwrap(), (14, 6));
        assert_eq!(try_csi_levels(0, 14).unwrap(), (14, 0));
        assert_eq!(try_csi_levels(1 << 50, 14).unwrap(), (24, 9));
    }

    #[test]
    fn sparse_first_offsets_match_dense_sorted_records_and_gap_resolution() {
        for seed in 0..64u64 {
            let mut sparse = Linear::Sparse {
                len: 0,
                runs: Vec::new(),
            };
            let mut dense = Vec::new();
            let mut start = 0;
            let mut rng = seed + 1;
            for record in 0..200 {
                rng = rng.wrapping_mul(6364136223846793005).wrapping_add(1);
                start += (rng >> 59) as usize;
                let end = start + ((rng >> 32) as usize % 300);
                let offset = record * 17 + 1;
                if dense.len() <= end {
                    dense.resize(end + 1, UNSET);
                }
                for value in &mut dense[start..=end] {
                    if *value == UNSET {
                        *value = offset;
                    }
                }
                sparse.extend(start, end, offset, true);
            }
            for (i, &value) in dense.iter().enumerate() {
                assert_eq!(sparse.get(i), Some(value));
            }
            let resolve = |value| value * 3 + 11;
            sparse.resolve(resolve);
            for value in &mut dense {
                if *value != UNSET {
                    *value = resolve(*value);
                }
            }
            for i in (0..dense.len().saturating_sub(1)).rev() {
                if dense[i] == UNSET {
                    dense[i] = dense[i + 1];
                }
            }
            for (i, &value) in dense.iter().enumerate() {
                assert_eq!(sparse.get(i), Some(value));
            }
            assert_eq!(sparse.get(dense.len()), None);
        }
    }

    #[test]
    fn csi_wide_spans_use_bounded_runs_while_bai_remains_dense() {
        let mut csi = Linear::default();
        csi.extend(10, 20, 7, true);
        csi.extend(1_000_000_000, 1_000_000_100, 9, true);
        let Linear::Sparse { len, runs } = &csi else {
            panic!("expected sparse CSI windows")
        };
        assert_eq!(*len, 1_000_000_101);
        assert_eq!(runs.len(), 4);
        assert!(runs.capacity() <= 8);
        assert_eq!(csi.get(10), Some(7));
        assert_eq!(csi.get(500_000_000), Some(UNSET));
        csi.resolve(|value| value);
        assert_eq!(csi.get(0), Some(7));
        assert_eq!(csi.get(500_000_000), Some(9));

        let mut bai = Linear::default();
        bai.extend(0, Linear::DENSE_LIMIT + 10, 7, false);
        assert!(matches!(bai, Linear::Dense(_)));
    }

    #[test]
    fn sparse_empty_windows_and_unset_tails_match_dense_offsets() {
        let limit = Linear::DENSE_LIMIT;
        let cases = [
            vec![(limit + 2, limit + 1, 7)],
            vec![(limit, limit - 1, 5), (limit + 2, limit + 1, 7)],
            vec![(limit + 2, limit + 1, 7), (limit + 10, limit + 9, 11)],
            vec![(3, 7, 1), (limit + 2, limit + 1, 7)],
            vec![(0, limit + 20, 1), (limit + 8, limit + 7, 7)],
            vec![(0, limit + 1, 1), (limit + 8, limit + 7, 7)],
            vec![
                (0, limit + 1, 1),
                (limit + 8, limit + 7, 7),
                (limit + 9, limit + 10, 11),
            ],
            vec![(limit + 2, limit + 1, 7), (limit + 2, limit + 2, 11)],
        ];
        for records in cases {
            let mut actual = Linear::default();
            let mut expected = Vec::new();
            for (lo, hi, offset) in records {
                expected.resize(expected.len().max(hi + 1), UNSET);
                if lo <= hi {
                    for value in &mut expected[lo..=hi] {
                        if *value == UNSET {
                            *value = offset;
                        }
                    }
                }
                actual.extend(lo, hi, offset, true);
            }
            for (i, &offset) in expected.iter().enumerate() {
                assert_eq!(actual.get(i), Some(offset));
            }
            actual.resolve(|offset| offset * 3 + 1);
            for offset in &mut expected {
                if *offset != UNSET {
                    *offset = *offset * 3 + 1;
                }
            }
            for i in (0..expected.len().saturating_sub(1)).rev() {
                if expected[i] == UNSET {
                    expected[i] = expected[i + 1];
                }
            }
            for (i, &offset) in expected.iter().enumerate() {
                assert_eq!(actual.get(i), Some(offset));
            }
            let Linear::Sparse { runs, .. } = actual else {
                panic!("expected sparse")
            };
            assert!(runs.windows(2).all(|w| w[0].0 < w[1].0));
        }
    }
}
