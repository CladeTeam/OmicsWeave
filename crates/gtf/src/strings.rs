//! Interned strings: one byte arena, one offset per distinct string, and an
//! open-addressing table of string ids for lookup by content.

use std::collections::HashMap;
use std::hash::{BuildHasherDefault, Hasher};

use omicsweave_core::{Error, Result};

/// The id of "no string" in every `u32` string column.
pub const NONE: u32 = u32::MAX;

/// Distinct byte strings stored back to back in one arena.
///
/// String `i` is `data[offsets[i]..offsets[i + 1]]`; ids are dense, in
/// order of first interning. Lookup by content goes through an
/// open-addressing table of ids, so the pool holds no per-string heap
/// objects. Annotation text is ASCII or UTF-8; [`StringPool::str`] returns
/// `""` for a string that is not valid UTF-8 (use [`StringPool::bytes`]).
#[derive(Debug, Clone)]
pub struct StringPool {
    data: Vec<u8>,
    offsets: Vec<u32>,
    /// Ids by hash slot; `NONE` marks an empty slot.
    table: Vec<u32>,
}

impl Default for StringPool {
    fn default() -> StringPool {
        StringPool::new()
    }
}

/// 64-bit hash of a byte string (word-at-a-time multiply-xor with a final
/// avalanche), good enough for an open-addressing table.
#[inline]
pub(crate) fn hash_bytes(b: &[u8]) -> u64 {
    const K: u64 = 0x9e37_79b9_7f4a_7c15;
    let mut h = (b.len() as u64).wrapping_mul(K);
    let mut chunks = b.chunks_exact(8);
    for c in &mut chunks {
        let w = u64::from_le_bytes(c.try_into().unwrap());
        h = (h.rotate_left(23) ^ w).wrapping_mul(K);
    }
    let rest = chunks.remainder();
    if !rest.is_empty() {
        let mut w = [0u8; 8];
        w[..rest.len()].copy_from_slice(rest);
        h = (h.rotate_left(23) ^ u64::from_le_bytes(w)).wrapping_mul(K);
    }
    h ^= h >> 33;
    h = h.wrapping_mul(0xff51_afd7_ed55_8ccd);
    h ^= h >> 33;
    h
}

impl StringPool {
    pub fn new() -> StringPool {
        StringPool {
            data: Vec::new(),
            offsets: vec![0],
            table: vec![NONE; 16],
        }
    }

    /// Distinct strings held.
    pub fn len(&self) -> usize {
        self.offsets.len() - 1
    }

    pub fn is_empty(&self) -> bool {
        self.len() == 0
    }

    /// The bytes of string `id`.
    #[inline]
    pub fn bytes(&self, id: u32) -> &[u8] {
        let i = id as usize;
        &self.data[self.offsets[i] as usize..self.offsets[i + 1] as usize]
    }

    /// String `id` as text; `""` for [`NONE`] or non-UTF-8 bytes.
    #[inline]
    pub fn str(&self, id: u32) -> &str {
        if id == NONE {
            return "";
        }
        std::str::from_utf8(self.bytes(id)).unwrap_or("")
    }

    /// String `id`, or `None` for [`NONE`].
    #[inline]
    pub fn get_str(&self, id: u32) -> Option<&str> {
        (id != NONE).then(|| self.str(id))
    }

    /// The id of `s`, if interned.
    pub fn get(&self, s: &[u8]) -> Option<u32> {
        let mask = self.table.len() - 1;
        let mut slot = hash_bytes(s) as usize & mask;
        loop {
            let id = self.table[slot];
            if id == NONE {
                return None;
            }
            if self.bytes(id) == s {
                return Some(id);
            }
            slot = (slot + 1) & mask;
        }
    }

    /// The id of `s`, interning it first if new.
    pub fn intern(&mut self, s: &[u8]) -> Result<u32> {
        let h = hash_bytes(s);
        let mask = self.table.len() - 1;
        let mut slot = h as usize & mask;
        loop {
            let id = self.table[slot];
            if id == NONE {
                break;
            }
            if self.bytes(id) == s {
                return Ok(id);
            }
            slot = (slot + 1) & mask;
        }
        let id = self.len();
        let end = self.data.len() + s.len();
        if id >= NONE as usize - 1 || end > u32::MAX as usize {
            return Err(Error::Unsupported("annotation strings exceed 4 GiB"));
        }
        self.data.extend_from_slice(s);
        self.offsets.push(end as u32);
        self.table[slot] = id as u32;
        if (id + 1) * 2 > self.table.len() {
            self.grow();
        }
        Ok(id as u32)
    }

    fn grow(&mut self) {
        let cap = self.table.len() * 2;
        let mut table = vec![NONE; cap];
        let mask = cap - 1;
        for id in 0..self.len() as u32 {
            let mut slot = hash_bytes(self.bytes(id)) as usize & mask;
            while table[slot] != NONE {
                slot = (slot + 1) & mask;
            }
            table[slot] = id;
        }
        self.table = table;
    }

    /// Heap bytes held.
    pub fn heap_bytes(&self) -> usize {
        self.data.capacity() + 4 * (self.offsets.capacity() + self.table.capacity())
    }

    pub(crate) fn shrink(&mut self) {
        self.data.shrink_to_fit();
        self.offsets.shrink_to_fit();
    }
}

/// A multiply-shift hasher for integer keys (string ids, packed pairs).
#[derive(Default, Clone, Copy)]
pub struct IdHasher(u64);

impl Hasher for IdHasher {
    #[inline]
    fn finish(&self) -> u64 {
        self.0
    }

    fn write(&mut self, bytes: &[u8]) {
        for &b in bytes {
            self.write_u64(b as u64);
        }
    }

    #[inline]
    fn write_u32(&mut self, n: u32) {
        self.write_u64(n as u64);
    }

    #[inline]
    fn write_u64(&mut self, n: u64) {
        let h = (self.0 ^ n).wrapping_mul(0x9e37_79b9_7f4a_7c15);
        self.0 = h ^ (h >> 29);
    }
}

pub type IdMap<K, V> = HashMap<K, V, BuildHasherDefault<IdHasher>>;

pub fn id_map<K, V>(capacity: usize) -> IdMap<K, V> {
    HashMap::with_capacity_and_hasher(capacity, BuildHasherDefault::default())
}

/// Rows by string id, for keys that may repeat (a gene id on two
/// chromosomes, a shared gene name): one run per key, with its rows in
/// ascending order. Building takes expected linear time.
#[derive(Debug, Clone, Default)]
pub struct IdIndex {
    runs: IdMap<u32, (u32, u32)>,
    rows: Vec<u32>,
}

impl IdIndex {
    /// Indexes `keys[row]` for every row whose key is not [`NONE`].
    pub fn build(keys: &[u32]) -> IdIndex {
        assert!(keys.len() <= u32::MAX as usize, "too many indexed rows");
        let count = keys.iter().filter(|&&key| key != NONE).count();
        let mut runs: IdMap<u32, (u32, u32)> = id_map(count);
        for &key in keys {
            if key != NONE {
                runs.entry(key).or_default().1 += 1;
            }
        }
        runs.shrink_to_fit();
        let mut total = 0;
        for (at, n) in runs.values_mut() {
            *at = total;
            total += *n;
        }
        let mut rows = vec![0; total as usize];
        for (row, key) in keys.iter().enumerate() {
            if let Some((at, _)) = runs.get_mut(key) {
                rows[*at as usize] = row as u32;
                *at += 1;
            }
        }
        // The insertion cursor ends at the run's end; retain its start
        // for lookup, without allocating a second table of cursors.
        for (at, n) in runs.values_mut() {
            *at -= *n;
        }
        IdIndex { runs, rows }
    }

    /// The rows with key `id`, ascending.
    pub fn rows(&self, id: u32) -> &[u32] {
        match self.runs.get(&id) {
            Some(&(at, n)) => &self.rows[at as usize..(at + n) as usize],
            None => &[],
        }
    }

    pub fn heap_bytes(&self) -> usize {
        4 * self.rows.capacity() + 16 * self.runs.capacity()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn interns_and_finds_by_content() {
        let mut p = StringPool::new();
        let words: Vec<String> = (0..5000).map(|i| format!("w{}", i * 7919 % 3001)).collect();
        let ids: Vec<u32> = words
            .iter()
            .map(|w| p.intern(w.as_bytes()).unwrap())
            .collect();
        for (w, &id) in words.iter().zip(&ids) {
            assert_eq!(p.str(id), w);
            assert_eq!(p.get(w.as_bytes()), Some(id));
        }
        assert_eq!(p.len(), 3001);
        assert_eq!(p.get(b"absent"), None);
        assert_eq!(p.intern(b"").unwrap(), p.intern(b"").unwrap());
    }

    #[test]
    fn id_runs_preserve_rows_with_missing_and_repeated_keys() {
        for n in [0, 1, 2, 15, 257, 10_000] {
            let keys: Vec<u32> = (0..n)
                .map(|i| if i % 7 == 0 { NONE } else { i * 7919 % 113 })
                .collect();
            let index = IdIndex::build(&keys);
            for key in (0..115).chain([NONE]) {
                let want: Vec<u32> = keys
                    .iter()
                    .enumerate()
                    .filter_map(|(r, &k)| (k == key && k != NONE).then_some(r as u32))
                    .collect();
                assert_eq!(index.rows(key), want);
            }
        }
    }
}
