//! External sorting: records gather in a buffer of bounded size; a full
//! buffer is sorted (in parallel parts) and written to a run file; the
//! final buffer and the runs are merged back in key order.
//!
//! Keys must be unique (callers add a sequence number when the input order
//! must break ties), so the order — and every output built from it — does
//! not depend on where the buffer happened to spill. A record type whose
//! equal keys stand for one record (counts of one cell) folds them instead
//! ([`Spill::ABSORBS`]): run files then hold each key once, and the merged
//! stream may still yield a key once per run, for the caller to fold.

use std::cmp::Reverse;
use std::collections::BinaryHeap;
use std::fs::File;
use std::io::{BufReader, BufWriter, Read, Write};
use std::path::{Path, PathBuf};

use omicsweave_core::Result;
use omicsweave_core::exec::{Ordered, Pool};

use super::input::with_path;
use super::output::TempDir;

/// A record that can be sorted out of core.
pub trait Spill: Sized + Send + 'static {
    type Key: Ord + Copy + Send + Sync + 'static;
    fn key(&self) -> Self::Key;
    /// Bytes the record holds in memory, for the buffer's accounting.
    fn mem_size(&self) -> usize;
    fn encode(&self, out: &mut Vec<u8>);
    /// The next record of a run, `None` at its end.
    fn decode(input: &mut impl Read) -> Result<Option<Self>>;
    /// Whether records with equal keys are folded together
    /// ([`Spill::absorb`]) where the sorter writes them out.
    const ABSORBS: bool = false;
    /// Folds `next`, whose key equals this record's, into this record. The
    /// fold must be associative and commutative (exact integer sums, say),
    /// so the result does not depend on where the buffer spilled.
    fn absorb(&mut self, _next: Self) {}
}

/// Writes the records `next` yields, in key order, through `write`,
/// folding equal keys of a type that absorbs them.
fn write_folded<T: Spill>(
    mut next: impl FnMut() -> Result<Option<T>>,
    mut write: impl FnMut(&T) -> Result<()>,
) -> Result<()> {
    if !T::ABSORBS {
        while let Some(r) = next()? {
            write(&r)?;
        }
        return Ok(());
    }
    let mut pending: Option<T> = None;
    while let Some(r) = next()? {
        match &mut pending {
            Some(p) if p.key() == r.key() => p.absorb(r),
            _ => {
                if let Some(p) = pending.replace(r) {
                    write(&p)?;
                }
            }
        }
    }
    if let Some(p) = pending {
        write(&p)?;
    }
    Ok(())
}

/// Reads exactly `buf.len()` bytes, or none at a clean end.
pub fn read_or_end(input: &mut impl Read, buf: &mut [u8]) -> Result<bool> {
    let mut got = 0;
    while got < buf.len() {
        let n = match input.read(&mut buf[got..]) {
            Err(e) if e.kind() == std::io::ErrorKind::Interrupted => continue,
            result => result?,
        };
        if n == 0 {
            if got == 0 {
                return Ok(false);
            }
            return Err(omicsweave_core::Error::invalid(
                "spill",
                "truncated run file",
            ));
        }
        got += n;
    }
    Ok(true)
}

/// Records per buffered segment; segments are the units sorted in
/// parallel.
const SEGMENT: usize = 1 << 20;
/// Keys sampled per key range to choose the splitters.
const OVERSAMPLE: usize = 16;
/// Leave room for input/output files and other pipeline stages.
const MAX_FAN_IN: usize = 64;

fn merge_sizing(limit: u64) -> (usize, usize) {
    let bytes = (limit / 2).max(2);
    let fan_in = (bytes / (8 << 10)).clamp(2, MAX_FAN_IN as u64) as usize;
    let reader_buf = (bytes / fan_in as u64).clamp(1, 4 << 20) as usize;
    (fan_in, reader_buf)
}

/// A generation's run names are determined by its count and merge pass;
/// retaining one allocated path per spill would grow with the whole input.
#[derive(Default)]
struct RunSet {
    count: usize,
    pass: Option<usize>,
    /// A final singleton can keep its name without copying or renaming it.
    carry: Option<PathBuf>,
}

impl RunSet {
    fn path(&self, dir: &Path, index: usize) -> PathBuf {
        if self.count.checked_sub(1) == Some(index)
            && let Some(path) = &self.carry
        {
            return path.clone();
        }
        match self.pass {
            None => dir.join(format!("run-{index}")),
            Some(pass) => dir.join(format!("merge-{pass}-{index}")),
        }
    }
}

pub struct Sorter<'p, T: Spill> {
    pool: &'p Pool,
    segments: Vec<Vec<T>>,
    current: Vec<T>,
    len: usize,
    bytes: u64,
    limit: u64,
    parent: PathBuf,
    dir: Option<TempDir>,
    runs: RunSet,
    spilled: u64,
}

impl<'p, T: Spill> Sorter<'p, T> {
    /// A sorter holding up to `limit` bytes before it spills into a
    /// directory created inside `parent`.
    pub fn new(pool: &'p Pool, limit: u64, parent: &Path) -> Sorter<'p, T> {
        Sorter {
            pool,
            segments: Vec::new(),
            current: Vec::new(),
            len: 0,
            bytes: 0,
            limit: limit.max(1),
            parent: parent.to_owned(),
            dir: None,
            runs: RunSet::default(),
            spilled: 0,
        }
    }

    fn account(&mut self, r: &T) {
        self.bytes = self
            .bytes
            .saturating_add(r.mem_size() as u64)
            .saturating_add(std::mem::size_of::<T>() as u64);
        self.len += 1;
    }

    /// Adds a batch of records as one segment.
    pub fn extend(&mut self, records: Vec<T>) -> Result<()> {
        if records.is_empty() {
            return Ok(());
        }
        for r in &records {
            self.account(r);
        }
        if records.len() < SEGMENT / 16 {
            self.current.extend(records);
            if self.current.len() >= SEGMENT {
                self.segments.push(std::mem::take(&mut self.current));
            }
        } else {
            self.segments.push(records);
        }
        self.check()
    }

    fn check(&mut self) -> Result<()> {
        if self.bytes >= self.limit {
            self.spill()?;
        }
        Ok(())
    }

    /// Bytes written to initial run files so far. Intermediate merge I/O
    /// performed by `finish` is not included.
    pub fn spilled(&self) -> u64 {
        self.spilled
    }

    /// Sorts the buffer as a sample sort: splitters drawn from the keys cut
    /// the key space into ranges; every segment is sorted and cut at the
    /// splitters on a worker; each range's pieces are merged later, also on
    /// a worker. Returns the ranges' pieces in key order.
    fn sorted_ranges(&mut self) -> Result<Vec<Vec<Vec<T>>>> {
        let mut segs = std::mem::take(&mut self.segments);
        if !self.current.is_empty() {
            segs.push(std::mem::take(&mut self.current));
        }
        let n = self.len;
        self.len = 0;
        self.bytes = 0;
        if n == 0 {
            return Ok(Vec::new());
        }
        let threads = self.pool.threads();
        let ranges = if n < 1 << 16 {
            1
        } else {
            threads.saturating_mul(8).min(n >> 14).max(1)
        };
        let splitters = {
            let want = ranges * OVERSAMPLE;
            let mut sample = Vec::with_capacity(want);
            for (j, seg) in segs.iter().enumerate() {
                let len = seg.len() as u128;
                let take = ((len * want as u128).div_ceil(n as u128)).min(len) as usize;
                // Each segment's samples start at an offset of its own: a
                // segment already in key order (a sorted batch) then still
                // contributes keys from across the key space rather than
                // the same few quantiles every other one does.
                let offset = spread(j as u64) as u128 % (len / take as u128).max(1);
                for i in 0..take {
                    sample.push(
                        seg[((offset + i as u128 * len / take as u128) % len) as usize].key(),
                    );
                }
            }
            sample.sort_unstable();
            let mut out: Vec<T::Key> = (1..ranges)
                .map(|j| sample[(j as u128 * sample.len() as u128 / ranges as u128) as usize])
                .collect();
            out.dedup();
            std::sync::Arc::new(out)
        };
        let n_ranges = splitters.len() + 1;
        let mut q = self.pool.ordered::<Vec<Vec<T>>>(threads.saturating_mul(2));
        let mut by_range: Vec<Vec<Vec<T>>> = (0..n_ranges).map(|_| Vec::new()).collect();
        let mut take = |pieces: Vec<Vec<T>>| {
            for (r, p) in pieces.into_iter().enumerate() {
                if !p.is_empty() {
                    by_range[r].push(p);
                }
            }
        };
        for mut seg in segs {
            while q.full() {
                take(q.next(true)?.expect("pending"));
            }
            let splitters = std::sync::Arc::clone(&splitters);
            q.submit(Box::new(move || {
                seg.sort_unstable_by_key(|r| r.key());
                let cuts: Vec<usize> = splitters
                    .iter()
                    .map(|s| seg.partition_point(|r| r.key() < *s))
                    .collect();
                let mut pieces = Vec::with_capacity(cuts.len() + 1);
                for &at in cuts.iter().rev() {
                    pieces.push(seg.split_off(at));
                }
                pieces.push(seg);
                pieces.reverse();
                Ok(pieces)
            }))?;
        }
        while let Some(p) = q.next(true)? {
            take(p);
        }
        Ok(by_range)
    }

    fn spill(&mut self) -> Result<()> {
        let ranges = self.sorted_ranges()?;
        if ranges.is_empty() {
            return Ok(());
        }
        if self.dir.is_none() {
            self.dir = Some(TempDir::new(&self.parent)?);
        }
        let path = self.runs.path(
            self.dir.as_ref().expect("created above").path(),
            self.runs.count,
        );
        let file = File::create(&path).map_err(|e| with_path(e, &path))?;
        let mut w = BufWriter::with_capacity(1 << 20, file);
        let mut buf = Vec::new();
        let mut merged = Ranges::new(self.pool, ranges);
        let mut spilled = self.spilled;
        write_folded(
            || merged.next_record(),
            |r| {
                buf.clear();
                r.encode(&mut buf);
                spilled = spilled.saturating_add(buf.len() as u64);
                w.write_all(&buf)?;
                Ok(())
            },
        )?;
        self.spilled = spilled;
        w.flush()?;
        self.runs.count += 1;
        Ok(())
    }

    /// Every record in key order.
    pub fn finish(mut self) -> Result<Merge<T>> {
        let (fan_in, reader_buf) = merge_sizing(self.limit);
        self.compact_runs(fan_in, reader_buf)?;
        let ranges = Ranges::new(self.pool, self.sorted_ranges()?);
        if self.runs.count == 0 {
            return Ok(Merge::Ranges(ranges));
        }
        let mut sources = Vec::with_capacity(self.runs.count + 1);
        for i in 0..self.runs.count {
            let path = self
                .runs
                .path(self.dir.as_ref().expect("runs have a directory").path(), i);
            let f = File::open(&path).map_err(|e| with_path(e, &path))?;
            sources.push(Source::Run(BufReader::with_capacity(reader_buf, f)));
        }
        sources.push(Source::Ranges(ranges));
        let mut heap = Heap::new(sources)?;
        heap.dir = self.dir.take();
        Ok(Merge::Heap(heap))
    }

    /// Bound open run files and reader buffers independently of run count.
    fn compact_runs(&mut self, fan_in: usize, reader_buf: usize) -> Result<()> {
        let mut pass = 0;
        while self.runs.count > fan_in {
            let runs = std::mem::take(&mut self.runs);
            self.runs.pass = Some(pass);
            let dir = self.dir.as_ref().expect("runs have a directory").path();
            for begin in (0..runs.count).step_by(fan_in) {
                let end = begin.saturating_add(fan_in).min(runs.count);
                let mut paths: Vec<_> = (begin..end).map(|i| runs.path(dir, i)).collect();
                let path = self.runs.path(dir, self.runs.count);
                if paths.len() == 1 {
                    self.runs.carry = paths.pop();
                    self.runs.count += 1;
                    continue;
                }
                let sources = paths
                    .iter()
                    .map(|p| {
                        let file = File::open(p).map_err(|e| with_path(e, p))?;
                        Ok(Source::Run(BufReader::with_capacity(reader_buf, file)))
                    })
                    .collect::<Result<Vec<_>>>()?;
                let mut merged = Heap::<T>::new(sources)?;
                let file = File::create(&path).map_err(|e| with_path(e, &path))?;
                let mut out = BufWriter::with_capacity(reader_buf, file);
                let mut buf = Vec::new();
                write_folded(
                    || merged.next_record(),
                    |record| {
                        buf.clear();
                        record.encode(&mut buf);
                        out.write_all(&buf)?;
                        Ok(())
                    },
                )?;
                out.flush()?;
                drop(out);
                drop(merged);
                for p in &paths {
                    std::fs::remove_file(p).map_err(|e| with_path(e, p))?;
                }
                self.runs.count += 1;
            }
            pass += 1;
        }
        Ok(())
    }
}

/// A well-mixed function of `i` (splitmix64's finalizer).
fn spread(i: u64) -> u64 {
    let mut z = i.wrapping_add(0x9E37_79B9_7F4A_7C15);
    z = (z ^ (z >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
    z = (z ^ (z >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
    z ^ (z >> 31)
}

/// Merges sorted pieces into one sorted vector.
fn merge_pieces<T: Spill>(mut pieces: Vec<Vec<T>>) -> Vec<T> {
    if pieces.len() == 1 {
        return pieces.pop().expect("one piece");
    }
    let total = pieces.iter().map(Vec::len).sum();
    let mut previous = None;
    let disjoint = pieces.iter().filter(|p| !p.is_empty()).all(|piece| {
        let first = piece.first().expect("nonempty piece").key();
        let ordered = previous.is_none_or(|last| last < first);
        previous = Some(piece.last().expect("nonempty piece").key());
        ordered
    });
    if disjoint {
        // Input segments often already follow key order. Their certified
        // disjoint ranges need no per-record heap comparisons or reversal.
        let mut pieces = pieces.into_iter();
        let mut out = pieces.next().unwrap_or_default();
        out.reserve(total - out.len());
        for mut piece in pieces {
            out.append(&mut piece);
        }
        return out;
    }
    let mut out = Vec::with_capacity(total);
    for p in &mut pieces {
        p.reverse();
    }
    let mut heap: BinaryHeap<Reverse<(T::Key, usize)>> = pieces
        .iter()
        .enumerate()
        .filter_map(|(i, p)| p.last().map(|r| Reverse((r.key(), i))))
        .collect();
    while let Some(Reverse((_, i))) = heap.pop() {
        out.push(pieces[i].pop().expect("heap entry has a record"));
        if let Some(r) = pieces[i].last() {
            heap.push(Reverse((r.key(), i)));
        }
    }
    out
}

/// The key ranges of a sample sort, merged on workers a few ranges ahead
/// of the reader and handed out in order.
pub struct Ranges<T: Spill> {
    todo: std::collections::VecDeque<Vec<Vec<T>>>,
    queue: omicsweave_core::exec::PoolQueue<Vec<T>>,
    current: std::vec::IntoIter<T>,
}

impl<T: Spill> Ranges<T> {
    fn new(pool: &Pool, ranges: Vec<Vec<Vec<T>>>) -> Ranges<T> {
        Ranges {
            todo: ranges.into(),
            queue: pool.ordered(pool.threads()),
            current: Vec::new().into_iter(),
        }
    }

    fn next_record(&mut self) -> Result<Option<T>> {
        loop {
            if let Some(r) = self.current.next() {
                return Ok(Some(r));
            }
            while !self.queue.full() {
                let Some(pieces) = self.todo.pop_front() else {
                    break;
                };
                self.queue
                    .submit(Box::new(move || Ok(merge_pieces(pieces))))?;
            }
            match self.queue.next(true)? {
                Some(v) => self.current = v.into_iter(),
                None => return Ok(None),
            }
        }
    }
}

enum Source<T: Spill> {
    Run(BufReader<File>),
    Ranges(Ranges<T>),
}

impl<T: Spill> Source<T> {
    fn next_record(&mut self) -> Result<Option<T>> {
        match self {
            Source::Run(r) => T::decode(r),
            Source::Ranges(r) => r.next_record(),
        }
    }
}

/// Run files and the final buffer merged by key.
pub struct Heap<T: Spill> {
    sources: Vec<Source<T>>,
    heads: Vec<Option<T>>,
    heap: BinaryHeap<Reverse<(T::Key, usize)>>,
    /// Kept alive until the runs are read.
    dir: Option<TempDir>,
}

impl<T: Spill> Heap<T> {
    fn new(sources: Vec<Source<T>>) -> Result<Heap<T>> {
        let mut m = Heap {
            heads: (0..sources.len()).map(|_| None).collect(),
            sources,
            heap: BinaryHeap::new(),
            dir: None,
        };
        for i in 0..m.sources.len() {
            m.advance(i)?;
        }
        Ok(m)
    }

    fn advance(&mut self, i: usize) -> Result<()> {
        if let Some(rec) = self.sources[i].next_record()? {
            self.heap.push(Reverse((rec.key(), i)));
            self.heads[i] = Some(rec);
        }
        Ok(())
    }

    fn next_record(&mut self) -> Result<Option<T>> {
        let Some(mut head) = self.heap.peek_mut() else {
            return Ok(None);
        };
        let i = head.0.1;
        let record = self.heads[i].take().expect("heap entry has a head");
        match self.sources[i].next_record()? {
            Some(next) => {
                *head = Reverse((next.key(), i));
                self.heads[i] = Some(next);
            }
            None => {
                std::collections::binary_heap::PeekMut::pop(head);
            }
        }
        Ok(Some(record))
    }
}

/// Records in key order.
pub enum Merge<T: Spill> {
    Ranges(Ranges<T>),
    Heap(Heap<T>),
}

impl<T: Spill> Merge<T> {
    pub fn next_record(&mut self) -> Result<Option<T>> {
        match self {
            Merge::Ranges(r) => r.next_record(),
            Merge::Heap(h) => h.next_record(),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::sync::atomic::{AtomicUsize, Ordering};

    static MERGE_COMPARISONS: AtomicUsize = AtomicUsize::new(0);

    #[derive(Clone, Copy, Debug, PartialEq, Eq)]
    struct CountedKey(u64);

    impl Ord for CountedKey {
        fn cmp(&self, other: &Self) -> std::cmp::Ordering {
            MERGE_COMPARISONS.fetch_add(1, Ordering::Relaxed);
            self.0.cmp(&other.0)
        }
    }

    impl PartialOrd for CountedKey {
        fn partial_cmp(&self, other: &Self) -> Option<std::cmp::Ordering> {
            Some(self.cmp(other))
        }
    }

    struct Counted(Rec);

    impl Spill for Counted {
        type Key = CountedKey;
        fn key(&self) -> Self::Key {
            CountedKey(self.0.0)
        }
        fn mem_size(&self) -> usize {
            self.0.mem_size()
        }
        fn encode(&self, out: &mut Vec<u8>) {
            self.0.encode(out)
        }
        fn decode(input: &mut impl Read) -> Result<Option<Self>> {
            Rec::decode(input).map(|record| record.map(Self))
        }
    }

    #[derive(Debug, PartialEq, Eq)]
    struct Rec(u64, Vec<u8>);

    impl Spill for Rec {
        type Key = u64;
        fn key(&self) -> u64 {
            self.0
        }
        fn mem_size(&self) -> usize {
            self.1.len()
        }
        fn encode(&self, out: &mut Vec<u8>) {
            out.extend_from_slice(&self.0.to_le_bytes());
            out.extend_from_slice(&(self.1.len() as u32).to_le_bytes());
            out.extend_from_slice(&self.1);
        }
        fn decode(input: &mut impl Read) -> Result<Option<Rec>> {
            let mut h = [0u8; 12];
            if !read_or_end(input, &mut h)? {
                return Ok(None);
            }
            let k = u64::from_le_bytes(h[..8].try_into().unwrap());
            let n = u32::from_le_bytes(h[8..].try_into().unwrap()) as usize;
            let mut b = vec![0u8; n];
            input.read_exact(&mut b)?;
            Ok(Some(Rec(k, b)))
        }
    }

    #[test]
    fn ordered_and_overlapping_pieces_match_sort_oracle() {
        for stride in [1, 17] {
            let mut pieces = vec![Vec::new()];
            let mut expected = Vec::new();
            for i in 0..8 {
                let mut piece = Vec::new();
                for j in 0..257 {
                    let k = if stride == 1 { i * 257 + j } else { j * 8 + i };
                    piece.push(Rec(k, vec![k as u8; (k % 7) as usize]));
                    expected.push(Rec(k, vec![k as u8; (k % 7) as usize]));
                }
                pieces.push(piece);
            }
            pieces.push(Vec::new());
            expected.sort_unstable_by_key(|r| r.0);
            assert_eq!(merge_pieces(pieces), expected);
        }
    }

    #[test]
    fn ordered_piece_merge_compares_ranges_instead_of_every_record() {
        let pieces = (0..16)
            .map(|i| {
                (0..4096)
                    .map(|j| Counted(Rec(i * 4096 + j, vec![])))
                    .collect()
            })
            .collect();
        MERGE_COMPARISONS.store(0, Ordering::Relaxed);
        let merged = merge_pieces(pieces);
        let comparisons = MERGE_COMPARISONS.load(Ordering::Relaxed);
        eprintln!(
            "ordered piece merge: records={}, pieces=16, comparisons={comparisons}",
            merged.len()
        );
        assert!(
            comparisons < 64,
            "per-record merge comparisons: {comparisons}"
        );
        assert_eq!(merged.len(), 16 * 4096);
        assert!(merged.iter().enumerate().all(|(i, r)| r.0.0 == i as u64));
    }

    #[test]
    fn merge_resources_are_bounded_by_budget_and_fan_in() {
        for limit in [0, 1, 3_000, 1 << 20, u64::MAX] {
            let (fan_in, reader_buf) = merge_sizing(limit);
            assert!((2..=MAX_FAN_IN).contains(&fan_in));
            assert!(reader_buf >= 1);
            assert!(fan_in as u64 * reader_buf as u64 <= (limit / 2).max(2));
        }
    }

    #[test]
    fn multi_pass_merge_matches_sorted_oracle_and_cleans_up() {
        let pool = Pool::new(2);
        let dir = tempfile::tempdir().unwrap();
        // A two-way merge forces several complete generations of runs.
        let mut sorter = Sorter::new(&pool, 1, dir.path());
        let mut expected = Vec::new();
        for i in (0..37).rev() {
            let key = i * 17 % 37;
            expected.push(key);
            sorter.extend(vec![Rec(key, vec![key as u8; 5])]).unwrap();
        }
        expected.sort_unstable();
        let mut merged = sorter.finish().unwrap();
        let Merge::Heap(heap) = &merged else {
            panic!("expected disk merge")
        };
        assert!(heap.sources.len() <= 3);
        let mut actual = Vec::new();
        while let Some(record) = merged.next_record().unwrap() {
            assert_eq!(record.1, vec![record.0 as u8; 5]);
            actual.push(record.0);
        }
        assert_eq!(actual, expected);
        drop(merged);
        assert_eq!(std::fs::read_dir(dir.path()).unwrap().count(), 0);
    }

    #[test]
    fn truncated_intermediate_run_and_early_drop_clean_all_runs() {
        let pool = Pool::new(1);
        let dir = tempfile::tempdir().unwrap();
        let mut sorter = Sorter::new(&pool, 1, dir.path());
        for key in 0..5 {
            sorter.extend(vec![Rec(key, Vec::new())]).unwrap();
        }
        let path = sorter.runs.path(sorter.dir.as_ref().unwrap().path(), 0);
        std::fs::write(path, [0]).unwrap();
        assert!(sorter.finish().is_err());
        assert_eq!(std::fs::read_dir(dir.path()).unwrap().count(), 0);
        let mut sorter = Sorter::new(&pool, 1, dir.path());
        for key in 0..5 {
            sorter.extend(vec![Rec(key, Vec::new())]).unwrap();
        }
        let mut merged = sorter.finish().unwrap();
        assert!(merged.next_record().unwrap().is_some());
        drop(merged);
        assert_eq!(std::fs::read_dir(dir.path()).unwrap().count(), 0);
    }

    #[test]
    fn sorts_the_same_with_and_without_spilling() {
        let pool = Pool::new(3);
        let dir = tempfile::tempdir().unwrap();
        let keys: Vec<u64> = (0..200_000u64)
            .map(|i| i.wrapping_mul(2_654_435_761) % 1_000_003 * 100_000 + i)
            .collect();
        let mut expect = keys.clone();
        expect.sort();
        for limit in [1 << 40, 50_000, 3_000] {
            let mut s = Sorter::new(&pool, limit, dir.path());
            for chunk in keys.chunks(777) {
                let recs = chunk
                    .iter()
                    .map(|&k| Rec(k, vec![k as u8; (k % 5) as usize]))
                    .collect();
                s.extend(recs).unwrap();
            }
            let spilled = s.spilled() > 0;
            assert_eq!(spilled, limit < 1 << 40);
            let mut m = s.finish().unwrap();
            let mut got = Vec::new();
            while let Some(r) = m.next_record().unwrap() {
                assert_eq!(r.1, vec![r.0 as u8; (r.0 % 5) as usize]);
                got.push(r.0);
            }
            assert_eq!(got, expect);
        }
        assert_eq!(std::fs::read_dir(dir.path()).unwrap().count(), 0);
    }
}
