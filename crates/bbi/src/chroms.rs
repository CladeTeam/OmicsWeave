//! The chromosome B+ tree: name -> (id, length). The names are lightweight
//! metadata, so the parsed table owns them.

use std::collections::{BTreeMap, HashMap};

use crate::bytes::{locate, require, u16_at, u32_at, u64_at};
use crate::error::{Error, Result, invalid};

pub const CHROM_TREE_MAGIC: u32 = 0x78CA_8C91;
const TREE_HEADER_LEN: usize = 32;
const NODE_HEADER_LEN: usize = 4;
/// Values hold the chromosome id and length.
const VALUE_LEN: u32 = 8;
/// Bounds the traversal so a cyclic tree fails instead of looping.
const MAX_DEPTH: usize = 32;

/// Dense chromosome index used by the data and the R-tree.
pub type ChromId = u32;

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Chrom {
    pub name: String,
    pub len: u32,
}

#[derive(Debug, Clone, Default, PartialEq, Eq)]
pub struct Chroms {
    by_id: Vec<Chrom>,
    by_name: HashMap<String, ChromId>,
}

impl Chroms {
    /// Parses the tree in `buf`, which holds the bytes of
    /// [`crate::Layout::chrom_tree`] read from file offset `base`.
    pub fn parse(buf: &[u8], base: u64) -> Result<Chroms> {
        require(buf, TREE_HEADER_LEN, "chrom tree")?;
        locate(base, buf, base, buf.len(), "chrom tree")?;
        let magic = u32_at(buf, 0);
        if magic != CHROM_TREE_MAGIC {
            return Err(Error::BadMagic {
                what: "chrom tree",
                found: magic,
            });
        }
        let key_size = u32_at(buf, 8);
        let value_size = u32_at(buf, 12);
        let count = u64_at(buf, 16);
        if key_size == 0 {
            return Err(invalid("chrom tree", "zero key size"));
        }
        if value_size != VALUE_LEN {
            return Err(invalid(
                "chrom tree",
                format!("value size {value_size}, expected {VALUE_LEN}"),
            ));
        }
        // Each item takes at least key + value bytes, which bounds a
        // plausible count by the buffer before allocating for it.
        if count > buf.len() as u64 / (key_size as u64 + VALUE_LEN as u64) {
            return Err(invalid("chrom tree", format!("{count} items cannot fit")));
        }

        let mut slots: Vec<Option<Chrom>> = vec![None; count as usize];
        let mut seen = 0u64;
        let key_size = key_size as usize;
        let mut reader = TreeReader {
            buf,
            base,
            key_size,
            remaining_nodes: (buf.len() / NODE_HEADER_LEN).max(MAX_DEPTH),
            occupied: BTreeMap::new(),
        };
        reader.visit(base + TREE_HEADER_LEN as u64, 0, &mut |key, id, len| {
            let name_len = memchr::memchr(0, key).unwrap_or(key_size);
            let name = std::str::from_utf8(&key[..name_len])
                .map_err(|e| invalid("chrom name", e.to_string()))?;
            let slot = slots
                .get_mut(id as usize)
                .ok_or_else(|| invalid("chrom tree", format!("id {id} >= item count {count}")))?;
            if slot.is_some() {
                return Err(invalid("chrom tree", format!("duplicate id {id}")));
            }
            *slot = Some(Chrom {
                name: name.to_owned(),
                len,
            });
            seen += 1;
            Ok(())
        })?;
        if seen != count {
            return Err(invalid(
                "chrom tree",
                format!("{seen} leaves for {count} items"),
            ));
        }

        let by_id: Vec<Chrom> = slots
            .into_iter()
            .map(|c| c.expect("every id seen once"))
            .collect();
        let mut by_name = HashMap::with_capacity(by_id.len());
        for (id, chrom) in by_id.iter().enumerate() {
            if by_name.insert(chrom.name.clone(), id as ChromId).is_some() {
                return Err(invalid(
                    "chrom tree",
                    format!("duplicate name {:?}", chrom.name),
                ));
            }
        }
        Ok(Chroms { by_id, by_name })
    }

    pub fn len(&self) -> usize {
        self.by_id.len()
    }

    pub fn is_empty(&self) -> bool {
        self.by_id.is_empty()
    }

    pub fn id(&self, name: &str) -> Option<ChromId> {
        self.by_name.get(name).copied()
    }

    pub fn get(&self, id: ChromId) -> Option<&Chrom> {
        self.by_id.get(id as usize)
    }

    /// Chromosomes in id order.
    pub fn iter(&self) -> impl ExactSizeIterator<Item = (ChromId, &Chrom)> {
        self.by_id
            .iter()
            .enumerate()
            .map(|(id, c)| (id as ChromId, c))
    }
}

// Recursion uses O(depth) state. Disjoint node ranges bound record work by
// encoded bytes; the O(nodes) map also rejects reused empty subtrees.
struct TreeReader<'a> {
    buf: &'a [u8],
    base: u64,
    key_size: usize,
    remaining_nodes: usize,
    occupied: BTreeMap<u64, u64>,
}

impl TreeReader<'_> {
    fn visit(
        &mut self,
        offset: u64,
        depth: usize,
        leaf_item: &mut impl FnMut(&[u8], ChromId, u32) -> Result<()>,
    ) -> Result<()> {
        if depth == MAX_DEPTH {
            return Err(Error::TooDeep(MAX_DEPTH));
        }
        if self.remaining_nodes == 0 {
            return Err(invalid("chrom tree", "too many node references"));
        }
        self.remaining_nodes -= 1;
        let at = locate(
            self.base,
            self.buf,
            offset,
            NODE_HEADER_LEN,
            "chrom tree node",
        )?;
        let leaf = self.buf[at] != 0;
        let n = u16_at(self.buf, at + 2) as usize;
        let item = self
            .key_size
            .checked_add(8)
            .ok_or_else(|| invalid("chrom tree node", "item size overflows usize"))?;
        let bytes = n
            .checked_mul(item)
            .ok_or_else(|| invalid("chrom tree node", "node size overflows usize"))?;
        let items = locate(
            self.base,
            self.buf,
            offset + NODE_HEADER_LEN as u64,
            bytes,
            "chrom tree node",
        )?;
        if offset < self.base + TREE_HEADER_LEN as u64 {
            return Err(invalid("chrom tree", "node overlaps tree header"));
        }
        let end = self.base + (items + bytes) as u64;
        if self
            .occupied
            .range(..=offset)
            .next_back()
            .is_some_and(|(_, &previous_end)| previous_end > offset)
            || self
                .occupied
                .range(offset..)
                .next()
                .is_some_and(|(&start, _)| start < end)
        {
            return Err(invalid("chrom tree", "repeated or overlapping node ranges"));
        }
        if depth > 0 || !leaf {
            self.occupied.insert(offset, end);
        }
        for i in 0..n {
            let at = items + i * item;
            if leaf {
                leaf_item(
                    &self.buf[at..at + self.key_size],
                    u32_at(self.buf, at + self.key_size),
                    u32_at(self.buf, at + self.key_size + 4),
                )?;
            } else {
                let child = u64_at(self.buf, at + self.key_size);
                self.visit(child, depth + 1, leaf_item)?;
            }
        }
        Ok(())
    }
}
