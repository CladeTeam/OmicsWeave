//! BAI and CSI indexes as htslib loads them (`idx_read`, `idx_read_core`,
//! `update_loff`), the file ranges of a region as `hts_itr_query` computes
//! them, and the walk over those ranges that `hts_itr_next` makes. BAM and
//! BCF share them: only the record layout differs.

use std::collections::HashMap;
use std::io::Read;
use std::path::{Path, PathBuf};

use omicsweave_core::exec::Inline;
use omicsweave_core::{Error, FileSource, ReadAt, Result};

use crate::seek::SeekReader;

/// Which index a file is.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum IndexKind {
    Bai,
    Csi,
}

impl IndexKind {
    /// `bai` or `csi`.
    pub fn name(self) -> &'static str {
        match self {
            IndexKind::Bai => "bai",
            IndexKind::Csi => "csi",
        }
    }
}

/// A range of virtual offsets `[start, end)`: compressed block address in
/// the high 48 bits, offset in the block's text in the low 16.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct Chunk {
    pub start: u64,
    pub end: u64,
}

#[derive(Debug, Clone, Default)]
struct Bin {
    loff: u64,
    /// In file order: the pseudo-bin's second pair holds counts.
    chunks: Vec<Chunk>,
}

#[derive(Debug, Clone, Default)]
struct RefIndex {
    bins: HashMap<u32, Bin>,
    /// BAI only: the smallest offset of a record in each `2^min_shift`
    /// window.
    linear: Vec<u64>,
    /// Allows a prefix search without assuming untrusted offsets are sorted.
    linear_nondecreasing: bool,
}

impl RefIndex {
    fn preceding_linear(&self, at: usize, ceiling: u64) -> Option<usize> {
        let prefix = &self.linear[..at];
        if self.linear_nondecreasing {
            prefix
                .partition_point(|&offset| offset < ceiling)
                .checked_sub(1)
        } else {
            prefix.iter().rposition(|&offset| offset < ceiling)
        }
    }
}

/// A loaded BAI or CSI index.
#[derive(Debug, Clone)]
pub struct BinIndex {
    kind: IndexKind,
    min_shift: u32,
    levels: u32,
    n_bins: u32,
    aux: Vec<u8>,
    refs: Vec<RefIndex>,
    n_no_coor: Option<u64>,
}

fn bin_first(level: u32) -> u32 {
    ((1u32 << (level * 3)) - 1) / 7
}

#[cfg(test)]
fn bin_parent(bin: u32) -> u32 {
    (bin - 1) >> 3
}

fn bin_level(bin: u32) -> u32 {
    let (mut l, mut b) = (0, bin);
    while b != 0 {
        l += 1;
        b = (b - 1) >> 3;
    }
    l
}

/// Little-endian fields from the front of the index bytes.
struct Cursor<'a> {
    bytes: &'a [u8],
    at: usize,
}

impl<'a> Cursor<'a> {
    fn take(&mut self, n: usize) -> Result<&'a [u8]> {
        let end = self
            .at
            .checked_add(n)
            .filter(|&e| e <= self.bytes.len())
            .ok_or(Error::Truncated {
                what: "index",
                need: (self.at as u64).saturating_add(n as u64),
                have: self.bytes.len() as u64,
            })?;
        let s = &self.bytes[self.at..end];
        self.at = end;
        Ok(s)
    }

    fn i32(&mut self) -> Result<i32> {
        Ok(i32::from_le_bytes(self.take(4)?.try_into().unwrap()))
    }

    fn u32(&mut self) -> Result<u32> {
        Ok(u32::from_le_bytes(self.take(4)?.try_into().unwrap()))
    }

    fn u64(&mut self) -> Result<u64> {
        Ok(u64::from_le_bytes(self.take(8)?.try_into().unwrap()))
    }

    fn count(&mut self, what: &str) -> Result<usize> {
        let n = self.i32()?;
        if n < 0 {
            return Err(Error::invalid(
                "index",
                format!("negative {what} count {n}"),
            ));
        }
        Ok(n as usize)
    }
}

/// The text of an index file: BGZF (or gzip) inflated, anything else as is.
fn index_text(raw: &[u8]) -> Result<Vec<u8>> {
    if !raw.starts_with(&[0x1f, 0x8b]) {
        return Ok(raw.to_vec());
    }
    let mut text = Vec::new();
    let decoded = if has_bgzf_extra(raw) {
        crate::Reader::new(raw, Box::new(Inline::new())).read_to_end(&mut text)
    } else {
        flate2::read::MultiGzDecoder::new(raw).read_to_end(&mut text)
    };
    decoded.map_err(
        |e| match e.into_inner().and_then(|e| e.downcast::<Error>().ok()) {
            Some(e) => *e,
            None => Error::Inflate("index does not inflate"),
        },
    )?;
    Ok(text)
}

/// Recognize the reserved subfield even if its payload is malformed, so an
/// invalid BGZF frame cannot bypass its checks through the generic decoder.
fn has_bgzf_extra(raw: &[u8]) -> bool {
    if raw.len() < 12 || raw[3] & 4 == 0 {
        return false;
    }
    let end = (12 + u16::from_le_bytes([raw[10], raw[11]]) as usize).min(raw.len());
    let mut at = 12;
    while at + 4 <= end {
        if &raw[at..at + 2] == b"BC" {
            return true;
        }
        let len = u16::from_le_bytes([raw[at + 2], raw[at + 3]]) as usize;
        if len > end - at - 4 {
            break;
        }
        at += 4 + len;
    }
    false
}

impl BinIndex {
    /// Parses an index file's bytes: a BAI (plain, as `samtools index`
    /// writes it, or BGZF) or a CSI (BGZF, or plain).
    pub fn parse(raw: &[u8]) -> Result<BinIndex> {
        let text = index_text(raw)?;
        let mut c = Cursor {
            bytes: &text,
            at: 0,
        };
        let magic = c.take(4)?;
        let (kind, min_shift, levels, aux) = match magic {
            b"BAI\x01" => (IndexKind::Bai, 14, 5, Vec::new()),
            b"CSI\x01" => {
                let min_shift = c.i32()?;
                let levels = c.i32()?;
                let l_aux = c.u32()? as usize;
                let aux = c.take(l_aux)?.to_vec();
                let bits = levels.checked_mul(3).and_then(|b| b.checked_add(min_shift));
                if min_shift < 0 || !(0..=10).contains(&levels) || bits.is_none_or(|b| b > 63) {
                    return Err(Error::invalid(
                        "index",
                        format!("CSI min_shift {min_shift} with {levels} levels"),
                    ));
                }
                (IndexKind::Csi, min_shift as u32, levels as u32, aux)
            }
            _ => {
                return Err(Error::BadMagic {
                    what: "BAM index",
                    found: u32::from_le_bytes(magic.try_into().unwrap()),
                });
            }
        };
        let n_bins = (((1u64 << (3 * levels + 3)) - 1) / 7) as u32;
        let n_ref = c.count("reference")?;
        let mut refs = Vec::with_capacity(n_ref.min(1 << 16));
        for _ in 0..n_ref {
            let mut r = RefIndex::default();
            let n = c.count("bin")?;
            for _ in 0..n {
                let key = c.u32()?;
                let loff = if kind == IndexKind::Csi { c.u64()? } else { 0 };
                let n_chunk = c.count("chunk")?;
                let mut chunks = Vec::with_capacity(n_chunk.min(1 << 16));
                for _ in 0..n_chunk {
                    let start = c.u64()?;
                    let end = c.u64()?;
                    chunks.push(Chunk { start, end });
                }
                if r.bins.insert(key, Bin { loff, chunks }).is_some() {
                    return Err(Error::invalid("index", format!("duplicate bin {key}")));
                }
            }
            if kind == IndexKind::Bai {
                let n = c.count("linear index")?;
                let raw = c.take(
                    n.checked_mul(8)
                        .ok_or_else(|| Error::invalid("index", "linear index size overflows"))?,
                )?;
                r.linear = raw
                    .chunks_exact(8)
                    .map(|b| u64::from_le_bytes(b.try_into().unwrap()))
                    .collect();
                fill_linear(&mut r.linear);
                update_loff(&mut r, levels, n_bins);
                r.linear_nondecreasing = r.linear.windows(2).all(|w| w[0] <= w[1]);
            }
            refs.push(r);
        }
        let n_no_coor = (text.len() - c.at >= 8).then(|| c.u64()).transpose()?;
        Ok(BinIndex {
            kind,
            min_shift,
            levels,
            n_bins,
            aux,
            refs,
            n_no_coor,
        })
    }

    /// Reads and parses the index file at `path`.
    pub fn read_file(path: impl AsRef<Path>) -> Result<BinIndex> {
        let source = FileSource::open(path)?;
        let len =
            usize::try_from(source.len()).map_err(|_| Error::invalid("index", "file too large"))?;
        let mut raw = vec![0; len];
        source.read_exact_at(&mut raw, 0)?;
        BinIndex::parse(&raw)
    }

    pub fn kind(&self) -> IndexKind {
        self.kind
    }

    /// Bases of the smallest bins are `2^min_shift`.
    pub fn min_shift(&self) -> u32 {
        self.min_shift
    }

    pub fn levels(&self) -> u32 {
        self.levels
    }

    /// The CSI auxiliary data (empty for BAI).
    pub fn aux(&self) -> &[u8] {
        &self.aux
    }

    /// References the index covers.
    pub fn n_refs(&self) -> usize {
        self.refs.len()
    }

    /// Records without coordinates, when the index records it.
    pub fn n_no_coor(&self) -> Option<u64> {
        self.n_no_coor
    }

    fn meta_bin(&self) -> u32 {
        self.n_bins + 1
    }

    /// Mapped and unmapped records of reference `tid`, from its pseudo-bin.
    pub fn counts(&self, tid: usize) -> Option<(u64, u64)> {
        let meta = self.refs.get(tid)?.bins.get(&self.meta_bin())?;
        (meta.chunks.len() >= 2).then(|| (meta.chunks[1].start, meta.chunks[1].end))
    }

    fn maxpos(&self) -> i64 {
        // A 63-bit CSI geometry covers all coordinates expressible by this
        // signed API; shifting signed one by 63 would make its bound negative.
        (1u64 << (self.min_shift + 3 * self.levels)).min(i64::MAX as u64) as i64
    }

    /// The bins of `r` that `[beg, end)` overlaps (`reg2bins`).
    fn overlapping_bins(&self, r: &RefIndex, beg: i64, end: i64) -> Vec<u32> {
        let levels = self.levels;
        let top_shift = self.min_shift + 3 * levels;
        let end = end.min(self.maxpos());
        if beg >= end {
            return Vec::new();
        }
        let end1 = end - 1;
        let count: i64 = (0..=levels)
            .map(|l| {
                let s = top_shift - 3 * l;
                (end1 >> s) - (beg >> s) + 1
            })
            .sum();
        let mut out = Vec::new();
        if (count as u64) < r.bins.len() as u64 {
            for l in 0..=levels {
                let s = top_shift - 3 * l;
                let t = bin_first(l) as i64;
                for b in t + (beg >> s)..=t + (end1 >> s) {
                    if r.bins.contains_key(&(b as u32)) {
                        out.push(b as u32);
                    }
                }
            }
        } else {
            for &bin in r.bins.keys() {
                let level = bin_level(bin);
                if level > levels {
                    continue;
                }
                let first = bin_first(level) as i64;
                let s = top_shift - 3 * level;
                if first + (beg >> s) <= bin as i64 && bin as i64 <= first + (end1 >> s) {
                    out.push(bin);
                }
            }
        }
        out
    }

    /// The file ranges that hold every record of reference `tid` overlapping
    /// `[beg, end)`, sorted, as `hts_itr_query` computes them: chunks of the
    /// overlapping bins above the linear-index bound, clipped and merged.
    /// Empty for a reference the index does not cover.
    pub fn chunks(&self, tid: usize, beg: i64, end: i64) -> Vec<Chunk> {
        let Some(r) = self.refs.get(tid) else {
            return Vec::new();
        };
        let beg = beg.max(0);
        if end <= beg || r.bins.is_empty() {
            return Vec::new();
        }
        let unmapped = match r.bins.get(&self.meta_bin()) {
            Some(m) if m.chunks.len() > 1 => m.chunks[1].end,
            _ => 1,
        };
        let maxpos = self.maxpos();
        if beg >= maxpos {
            return Vec::new();
        }
        let rel_off = (beg >> self.min_shift) as usize;
        // min_off: the loff of the nearest existing bin at or left of `beg`
        // in the bottom level, climbing up.
        let mut bin = bin_first(self.levels) as i64 + rel_off as i64;
        let mut found = None;
        while bin > 0 {
            if let Some(b) = r.bins.get(&(bin as u32)) {
                found = Some(b);
                break;
            }
            let parent = (bin - 1) >> 3;
            let first = (parent << 3) + 1;
            bin = if bin > first { bin - 1 } else { parent };
        }
        if bin <= 0 {
            found = r.bins.get(&0);
        }
        let mut min_off = found.map_or(0, |b| b.loff);
        if !r.linear.is_empty() && rel_off < r.linear.len() {
            min_off = min_off.max(r.linear[rel_off]);
            if unmapped != 0 {
                // Unmapped reads are not in the linear index: look back for
                // a smaller offset.
                let back = r.preceding_linear(rel_off, min_off);
                if let Some(l) = back {
                    min_off = r.linear[l];
                }
                if let Some(b) = found
                    && (min_off < b.loff || back.is_none())
                {
                    min_off = b.loff;
                }
            }
        } else if unmapped != 0
            && let Some(b) = found
        {
            min_off = b.loff;
        }
        // max_off: the first chunk of the nearest existing bin right of
        // `end`.
        let mut max_off = u64::MAX;
        if end <= maxpos {
            let mut bin = bin_first(self.levels) as i64 + ((end - 1) >> self.min_shift) + 1;
            if bin >= self.n_bins as i64 {
                bin = 0;
            }
            loop {
                while bin % 8 == 1 {
                    bin = (bin - 1) >> 3;
                }
                if bin == 0 {
                    break;
                }
                if let Some(b) = r.bins.get(&(bin as u32))
                    && let Some(c) = b.chunks.first()
                {
                    max_off = c.start;
                    break;
                }
                bin += 1;
            }
        }
        let mut off: Vec<Chunk> = Vec::new();
        for key in self.overlapping_bins(r, beg, end) {
            for c in &r.bins[&key].chunks {
                if c.end > min_off && c.start < max_off {
                    off.push(Chunk {
                        start: c.start.max(min_off),
                        end: c.end.min(max_off),
                    });
                }
            }
        }
        if off.is_empty() {
            return off;
        }
        off.sort_unstable();
        // Drop chunks contained in the one before.
        let mut l = 0;
        for i in 1..off.len() {
            if off[l].end < off[i].end {
                l += 1;
                off[l] = off[i];
            }
        }
        off.truncate(l + 1);
        // Trim overlaps left by merging bins at indexing time.
        for i in 1..off.len() {
            if off[i - 1].end >= off[i].start {
                off[i - 1].end = off[i].start;
            }
        }
        // Join chunks that meet in one block.
        let mut l = 0;
        for i in 1..off.len() {
            if off[l].end >> 16 == off[i].start >> 16 {
                off[l].end = off[i].end;
            } else {
                l += 1;
                off[l] = off[i];
            }
        }
        off.truncate(l + 1);
        off
    }
}

/// The index of the file at `path` where htslib looks for it:
/// `<path>.csi`, the name with its extension replaced by `.csi`, then the
/// same for `.bai`. `None` when there is none.
pub fn locate_index(path: impl AsRef<Path>) -> Option<PathBuf> {
    let path = path.as_ref();
    [".csi", ".bai"].into_iter().find_map(|ext| {
        let mut appended = path.as_os_str().to_owned();
        appended.push(ext);
        let replaced = path.with_extension(&ext[1..]);
        [PathBuf::from(appended), replaced]
            .into_iter()
            .find(|p| p.is_file())
    })
}

/// The chunks of one region query read in order, as `hts_itr_next` reads
/// them: a chunk ends once the read position reaches its end, and the
/// reader seeks only when the next chunk does not start where the last one
/// ended.
#[derive(Debug, Clone)]
pub struct ChunkWalk {
    chunks: Vec<Chunk>,
    at: Option<usize>,
}

impl ChunkWalk {
    pub fn new(chunks: Vec<Chunk>) -> ChunkWalk {
        ChunkWalk { chunks, at: None }
    }

    /// Positions `reader` for the next record of the query; `false` once
    /// every chunk is read.
    pub fn next<S: ReadAt + ?Sized>(
        &mut self,
        reader: &mut SeekReader,
        source: &S,
    ) -> Result<bool> {
        if self.at.is_some_and(|i| reader.tell() < self.chunks[i].end) {
            return Ok(true);
        }
        let next = self.at.map_or(0, |i| i + 1);
        let Some(chunk) = self.chunks.get(next) else {
            return Ok(false);
        };
        let adjacent = self.at.is_some_and(|i| self.chunks[i].end == chunk.start);
        if !adjacent {
            reader.seek(source, chunk.start, chunk.end)?;
        }
        self.at = Some(next);
        Ok(true)
    }
}

/// `idx_read_core`: zero entries before the first nonzero one stay; later
/// zeros take the next window's offset (old samtools left them).
fn fill_linear(linear: &mut [u64]) {
    let k = linear.iter().position(|&o| o != 0).unwrap_or(linear.len());
    for j in (k + 1..linear.len()).rev() {
        if linear[j - 1] == 0 {
            linear[j - 1] = linear[j];
        }
    }
}

/// `update_loff`: unset windows take the next one's offset, and each bin's
/// `loff` is its first window's offset.
fn update_loff(r: &mut RefIndex, levels: u32, n_bins: u32) {
    for l in (0..r.linear.len().saturating_sub(1)).rev() {
        if r.linear[l] == u64::MAX {
            r.linear[l] = r.linear[l + 1];
        }
    }
    for (&key, bin) in r.bins.iter_mut() {
        bin.loff = if key < n_bins {
            let l = bin_level(key);
            let bot = ((key - bin_first(l)) as u64) << ((levels - l) * 3);
            usize::try_from(bot)
                .ok()
                .and_then(|b| r.linear.get(b))
                .copied()
                .unwrap_or(0)
        } else {
            0
        };
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn bin_arithmetic() {
        assert_eq!(bin_first(5), 4681);
        assert_eq!(bin_level(4681), 5);
        assert_eq!(bin_level(0), 0);
        assert_eq!(bin_parent(4681), 585);
        assert_eq!(bin_level(37450), 6);
    }

    #[test]
    fn linear_fill_keeps_leading_zeros() {
        let mut l = vec![0, 0, 5, 0, 0, 9, 0];
        fill_linear(&mut l);
        assert_eq!(l, vec![0, 0, 5, 9, 9, 9, 0]);
    }

    #[test]
    fn preceding_linear_matches_reverse_scan_for_all_prefixes_and_thresholds() {
        for n in 0..=6 {
            for mut code in 0..4usize.pow(n) {
                let linear: Vec<_> = (0..n)
                    .map(|_| {
                        let value = [0, 1, 7, u64::MAX][code % 4];
                        code /= 4;
                        value
                    })
                    .collect();
                let r = RefIndex {
                    linear_nondecreasing: linear.windows(2).all(|w| w[0] <= w[1]),
                    linear,
                    ..RefIndex::default()
                };
                for at in 0..=r.linear.len() {
                    for ceiling in [0, 1, 2, 7, 8, u64::MAX] {
                        let expected = (0..at).rev().find(|&i| r.linear[i] < ceiling);
                        assert_eq!(r.preceding_linear(at, ceiling), expected);
                    }
                }
            }
        }
    }
}
