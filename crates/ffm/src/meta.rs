//! Per-motif metadata: IDs with an in-memory lookup table, and the index of
//! each motif's block and stream offsets.

use xxhash_rust::xxh3::xxh3_64;

use crate::block::RecordOffsets;
use crate::error::{Result, invalid};
use crate::layout::DirEntry;
use crate::schema::forbidden;

const EMPTY: u32 = u32::MAX;

/// Motif IDs in motif order, with an open-addressing table built in memory
/// (load factor at most 1/2) for expected O(1) lookup.
pub(crate) struct IdTable {
    text: String,
    offsets: Vec<u32>,
    slots: Vec<u32>,
}

impl IdTable {
    pub fn new() -> IdTable {
        IdTable {
            text: String::new(),
            offsets: vec![0],
            slots: vec![EMPTY; 2],
        }
    }

    /// Rebuilds the table from stored IDs, checking each and their
    /// uniqueness.
    pub fn from_stored(text: String, offsets: Vec<u32>, limit: usize) -> Result<IdTable> {
        const WHAT: &str = "FFM IDs";
        let n = offsets.len().saturating_sub(1);
        if offsets.first() != Some(&0) || offsets.last().map(|&e| e as usize) != Some(text.len()) {
            return Err(invalid(WHAT, "offsets do not span the ID bytes"));
        }
        let mut slots = crate::options::reserve(slot_count(n), limit, WHAT)?;
        slots.resize(slot_count(n), EMPTY);
        let mut table = IdTable {
            text,
            offsets,
            slots,
        };
        for i in 0..n {
            let (a, b) = (table.offsets[i] as usize, table.offsets[i + 1] as usize);
            if a >= b || !table.text.is_char_boundary(b) {
                return Err(invalid(WHAT, format!("ID {i} is empty or not UTF-8")));
            }
            check_id(&table.text[a..b])?;
            if !table.insert(i as u32) {
                return Err(invalid(WHAT, format!("ID {:?} repeats", &table.text[a..b])));
            }
        }
        Ok(table)
    }

    pub fn len(&self) -> usize {
        self.offsets.len() - 1
    }

    pub fn get(&self, i: usize) -> Option<&str> {
        let a = *self.offsets.get(i)? as usize;
        let b = *self.offsets.get(i + 1)? as usize;
        Some(&self.text[a..b])
    }

    pub fn find(&self, id: &str) -> Option<usize> {
        let mask = self.slots.len() - 1;
        let mut slot = xxh3_64(id.as_bytes()) as usize & mask;
        loop {
            let i = self.slots[slot];
            if i == EMPTY {
                return None;
            }
            if self.get(i as usize) == Some(id) {
                return Some(i as usize);
            }
            slot = (slot + 1) & mask;
        }
    }

    pub fn check_push(&self, id: &str) -> Result<()> {
        check_id(id)?;
        if self.find(id).is_some() {
            return Err(invalid("FFM IDs", format!("ID {id:?} repeats")));
        }
        let stored_len = self
            .offsets
            .len()
            .checked_add(1)
            .and_then(|n| n.checked_mul(4))
            .and_then(|n| n.checked_add(self.text.len()))
            .and_then(|n| n.checked_add(id.len()));
        if self.len() >= crate::MAX_MOTIFS || stored_len.is_none_or(|n| n > u32::MAX as usize) {
            return Err(invalid(
                "FFM IDs",
                "IDs stream or motif count exceeds format limits",
            ));
        }
        Ok(())
    }

    /// Appends an ID; fails if it is malformed or already present.
    pub fn push(&mut self, id: &str) -> Result<()> {
        self.check_push(id)?;
        let end = self.text.len() + id.len();
        self.text.push_str(id);
        self.offsets.push(end as u32);
        let n = self.len();
        if slot_count(n) > self.slots.len() {
            self.slots = vec![EMPTY; slot_count(n)];
            for i in 0..n {
                self.insert(i as u32);
            }
        } else {
            self.insert(n as u32 - 1);
        }
        Ok(())
    }

    /// Inserts motif `i`; false when an equal ID is already there.
    fn insert(&mut self, i: u32) -> bool {
        let id = self.get(i as usize).expect("stored ID");
        let mask = self.slots.len() - 1;
        let mut slot = xxh3_64(id.as_bytes()) as usize & mask;
        loop {
            let j = self.slots[slot];
            if j == EMPTY {
                self.slots[slot] = i;
                return true;
            }
            if self.get(j as usize) == Some(id) {
                return false;
            }
            slot = (slot + 1) & mask;
        }
    }

    /// The decoded IDs stream: offsets, then the ID bytes.
    pub fn stored(&self) -> Vec<u8> {
        let mut out = Vec::with_capacity(self.offsets.len() * 4 + self.text.len());
        for &o in &self.offsets {
            out.extend_from_slice(&o.to_le_bytes());
        }
        out.extend_from_slice(self.text.as_bytes());
        out
    }

    pub fn parse_stored(bytes: &[u8], motif_count: usize, limit: usize) -> Result<IdTable> {
        const WHAT: &str = "FFM IDs";
        let head = motif_count
            .checked_add(1)
            .and_then(|n| n.checked_mul(4))
            .ok_or_else(|| invalid(WHAT, "offset array too long"))?;
        if bytes.len() < head {
            return Err(invalid(WHAT, "shorter than its offsets"));
        }
        let offsets = bytes[..head]
            .chunks_exact(4)
            .map(|c| u32::from_le_bytes(c.try_into().expect("4")))
            .collect();
        let text = String::from_utf8(bytes[head..].to_vec())
            .map_err(|_| invalid(WHAT, "IDs are not UTF-8"))?;
        IdTable::from_stored(text, offsets, limit)
    }

    pub fn heap_bytes(&self) -> usize {
        self.text.capacity() + 4 * (self.offsets.capacity() + self.slots.capacity())
    }
}

fn slot_count(n: usize) -> usize {
    (2 * n).next_power_of_two().max(2)
}

fn check_id(id: &str) -> Result<()> {
    if id.is_empty() || id.bytes().any(forbidden) {
        return Err(invalid(
            "FFM IDs",
            format!("ID {id:?} is empty or contains a tab, newline or NUL"),
        ));
    }
    Ok(())
}

/// Each motif's block and starts within that block's streams, stored by
/// column: blocks, structure starts, value starts, annotation starts.
pub(crate) struct MotifIndex {
    n: usize,
    columns: Vec<u32>,
    /// First motif of each block.
    block_first: Vec<u32>,
}

impl MotifIndex {
    /// Checks the stored columns against the directory.
    pub fn new(columns: Vec<u32>, directory: &[DirEntry], annotated: bool) -> Result<MotifIndex> {
        const WHAT: &str = "FFM index";
        if !columns.len().is_multiple_of(4) {
            return Err(invalid(WHAT, "column lengths differ"));
        }
        let n = columns.len() / 4;
        let mut block_first = Vec::with_capacity(directory.len());
        let mut i = 0usize;
        for (b, entry) in directory.iter().enumerate() {
            block_first.push(i as u32);
            let mut prev: Option<RecordOffsets> = None;
            for _ in 0..entry.record_count {
                if i >= n || columns[i] as usize != b {
                    return Err(invalid(WHAT, format!("motif {i} is not in block {b}")));
                }
                let st = offsets(&columns, n, i);
                let ok = match prev {
                    None => st == RecordOffsets::default(),
                    Some(p) => {
                        st.structure > p.structure
                            && st.values > p.values
                            && (if annotated {
                                st.annotation > p.annotation
                            } else {
                                st.annotation == 0
                            })
                    }
                };
                if !ok {
                    return Err(invalid(
                        WHAT,
                        format!("motif {i} starts are not increasing within block {b}"),
                    ));
                }
                prev = Some(st);
                i += 1;
            }
        }
        if i != n {
            return Err(invalid(
                WHAT,
                "index and directory hold different motif counts",
            ));
        }
        Ok(MotifIndex {
            n,
            columns,
            block_first,
        })
    }

    pub fn block(&self, i: usize) -> usize {
        self.columns[i] as usize
    }

    pub fn offsets(&self, i: usize) -> RecordOffsets {
        offsets(&self.columns, self.n, i)
    }

    pub fn block_first(&self, b: usize) -> usize {
        self.block_first[b] as usize
    }

    /// RecordOffsets of every motif of block `b`.
    pub fn block_offsets(&self, b: usize, record_count: u32) -> Vec<RecordOffsets> {
        let first = self.block_first(b);
        (first..first + record_count as usize)
            .map(|i| self.offsets(i))
            .collect()
    }
}

fn offsets(columns: &[u32], n: usize, i: usize) -> RecordOffsets {
    RecordOffsets {
        structure: columns[n + i],
        values: columns[2 * n + i],
        annotation: columns[3 * n + i],
    }
}

/// Builds the stored columns from per-motif rows.
pub(crate) fn index_columns(blocks: &[u32], rows: &[RecordOffsets]) -> Vec<u32> {
    let mut out = Vec::with_capacity(rows.len() * 4);
    out.extend_from_slice(blocks);
    out.extend(rows.iter().map(|r| r.structure));
    out.extend(rows.iter().map(|r| r.values));
    out.extend(rows.iter().map(|r| r.annotation));
    out
}
