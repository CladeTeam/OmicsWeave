//! The R-tree that maps a genomic range to the data blocks overlapping it.
//!
//! The index bytes stay in the caller's buffer; queries walk the nodes in
//! place and allocate nothing.

use std::collections::BTreeMap;
use std::iter::FusedIterator;
use std::mem::MaybeUninit;
use std::ops::Range;

use crate::bytes::{locate, require, u16_at, u32_at, u32_unchecked, u64_at, u64_unchecked};
use crate::chroms::ChromId;
use crate::error::{Error, Result, invalid};
use crate::search::{Order, OrderScan, Point};

pub const INDEX_MAGIC: u32 = 0x2468_ACE0;
const INDEX_HEADER_LEN: usize = 48;
const NODE_HEADER_LEN: usize = 4;
const LEAF_ITEM_LEN: usize = 32;
const BRANCH_ITEM_LEN: usize = 24;
/// Far deeper than any index a file can hold; reaching it means a cycle.
pub const MAX_DEPTH: usize = 32;

/// A data block in the file: read `size` bytes at `offset`, inflate if the
/// file is compressed, then parse.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub struct BlockRef {
    pub offset: u64,
    pub size: u64,
}

impl BlockRef {
    /// # Panics
    /// Panics if a caller-constructed block's end overflows u64. Blocks
    /// returned by [`Index::query`] have already passed this check.
    pub fn range(&self) -> Range<u64> {
        let end = self
            .offset
            .checked_add(self.size)
            .expect("block range overflows u64");
        self.offset..end
    }
}

/// One leaf item of an R-tree: a data block and the (chromosome, base) span
/// the index records for it, `start` inclusive, `end` exclusive.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub struct Leaf {
    pub start: (ChromId, u32),
    pub end: (ChromId, u32),
    pub block: BlockRef,
}

/// Caller-owned storage for one prepared R-tree node. Fields are private so
/// query bounds can only come from validating the original index bytes.
#[derive(Debug, Clone, Copy, Default)]
pub struct NodeInfo {
    offset: u64,
    items: usize,
    count: u16,
    leaf: bool,
    order: Order,
}

/// The fields of the R-tree header.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct IndexHeader {
    /// Maximum children per node.
    pub block_size: u32,
    /// Leaf items, i.e. indexed data blocks.
    pub item_count: u64,
    pub start: (ChromId, u32),
    pub end: (ChromId, u32),
    /// End of the data region this tree indexes.
    pub end_file_offset: u64,
    pub items_per_slot: u32,
}

/// An R-tree whose bytes, starting with its header, were read from file
/// offset `base` into a caller buffer.
#[derive(Debug, Clone, Copy)]
pub struct Index<'a> {
    buf: &'a [u8],
    base: u64,
    header: IndexHeader,
    nodes: &'a [NodeInfo],
}

impl<'a> Index<'a> {
    /// Parses checked header metadata in O(1). Node ranges are checked lazily
    /// during query, or once up front by prepare for repeated queries.
    pub fn new(buf: &'a [u8], base: u64) -> Result<Index<'a>> {
        require(buf, INDEX_HEADER_LEN + NODE_HEADER_LEN, "index")?;
        locate(base, buf, base, buf.len(), "index")?;
        let magic = u32_at(buf, 0);
        if magic != INDEX_MAGIC {
            return Err(Error::BadMagic {
                what: "index",
                found: magic,
            });
        }
        let header = IndexHeader {
            block_size: u32_at(buf, 4),
            item_count: u64_at(buf, 8),
            start: (u32_at(buf, 16), u32_at(buf, 20)),
            end: (u32_at(buf, 24), u32_at(buf, 28)),
            end_file_offset: u64_at(buf, 32),
            items_per_slot: u32_at(buf, 40),
        };
        Ok(Index {
            buf,
            base,
            header,
            nodes: &[],
        })
    }

    pub fn header(&self) -> &IndexHeader {
        &self.header
    }

    /// Reborrows prepared storage held by an owning composition such as BbiReader.
    ///
    /// # Safety
    /// nodes must be the used slots produced by a successful prepare for exactly self.buf
    /// and self.base. Neither bytes nor metadata may have changed since that
    /// validation, and they must stay immutable for the returned view's lifetime.
    pub(crate) unsafe fn with_prepared_nodes(mut self, nodes: &'a [NodeInfo]) -> Self {
        self.nodes = nodes;
        self
    }

    /// Prepares all nodes once for repeated queries. `nodes` is cleared and
    /// filled in a single walk; reuse it to avoid allocating. Each node's
    /// items are visited once, certifying their order and their file ranges
    /// together. Cost: O(items + nodes * log(nodes)); the returned view
    /// borrows `nodes` and the caller's index bytes. Monotone bounds allow
    /// binary search within nodes; arbitrary overlaps use conservative
    /// bounds, preserving all hits and their file order.
    pub fn prepare<'w>(&self, nodes: &'w mut Vec<NodeInfo>) -> Result<Index<'w>>
    where
        'a: 'w,
    {
        nodes.clear();
        self.walk_nodes(&mut |mut node| {
            let item = if node.leaf {
                LEAF_ITEM_LEN
            } else {
                BRANCH_ITEM_LEN
            };
            let mut order = OrderScan::new();
            for i in 0..node.count as usize {
                let at = node.items + i * item;
                order.push(
                    (u32_at(self.buf, at), u32_at(self.buf, at + 4)),
                    (u32_at(self.buf, at + 8), u32_at(self.buf, at + 12)),
                );
                if node.leaf
                    && u64_at(self.buf, at + 16)
                        .checked_add(u64_at(self.buf, at + 24))
                        .is_none()
                {
                    return Err(invalid("data block", "file range overflows u64"));
                }
            }
            node.order = order.finish();
            nodes.push(node);
            Ok(())
        })?;
        nodes.sort_unstable_by_key(|node| node.offset);
        Ok(Index {
            buf: self.buf,
            base: self.base,
            header: self.header,
            nodes,
        })
    }

    fn read_node(&self, offset: u64) -> Result<NodeInfo> {
        let at = locate(self.base, self.buf, offset, NODE_HEADER_LEN, "index node")?;
        let leaf = self.buf[at] != 0;
        let count = u16_at(self.buf, at + 2);
        let item = if leaf { LEAF_ITEM_LEN } else { BRANCH_ITEM_LEN };
        let items = locate(
            self.base,
            self.buf,
            offset + NODE_HEADER_LEN as u64,
            count as usize * item,
            "index node",
        )?;
        Ok(NodeInfo {
            offset,
            items,
            count,
            leaf,
            order: Order::UNKNOWN,
        })
    }

    fn walk_nodes(&self, visitor: &mut impl FnMut(NodeInfo) -> Result<()>) -> Result<()> {
        let mut remaining = (self.buf.len() / NODE_HEADER_LEN).max(MAX_DEPTH);
        let mut occupied = BTreeMap::new();
        let mut ancestors = [0; MAX_DEPTH];
        self.walk_node(
            self.base + INDEX_HEADER_LEN as u64,
            0,
            &mut remaining,
            &mut occupied,
            &mut ancestors,
            visitor,
        )
    }

    fn walk_node(
        &self,
        offset: u64,
        depth: usize,
        remaining: &mut usize,
        occupied: &mut BTreeMap<u64, u64>,
        ancestors: &mut [u64; MAX_DEPTH],
        visitor: &mut impl FnMut(NodeInfo) -> Result<()>,
    ) -> Result<()> {
        if depth == MAX_DEPTH {
            return Err(Error::TooDeep(MAX_DEPTH));
        }
        if *remaining == 0 {
            return Err(invalid("index", "too many node references"));
        }
        *remaining -= 1;
        let node = self.read_node(offset)?;
        if offset < self.base + INDEX_HEADER_LEN as u64 {
            return Err(invalid("index", "node overlaps index header"));
        }
        // read_node validated this complete byte range. Disjoint ranges
        // bound total record validation by the index's actual byte length.
        let item = if node.leaf {
            LEAF_ITEM_LEN
        } else {
            BRANCH_ITEM_LEN
        };
        let end = self.base + (node.items + node.count as usize * item) as u64;
        if let Some((&start, &previous_end)) = occupied.range(..=offset).next_back() {
            if start == offset {
                return Err(if ancestors[..depth].contains(&offset) {
                    Error::TooDeep(MAX_DEPTH)
                } else {
                    invalid("index", "repeated node reference")
                });
            }
            if previous_end > offset {
                return Err(invalid("index", "overlapping node ranges"));
            }
        }
        if occupied
            .range(offset..)
            .next()
            .is_some_and(|(&start, _)| start < end)
        {
            return Err(invalid("index", "overlapping node ranges"));
        }
        if depth > 0 || !node.leaf {
            occupied.insert(offset, end);
        }
        ancestors[depth] = offset;
        visitor(node)?;
        if !node.leaf {
            for i in 0..node.count as usize {
                let child = u64_at(self.buf, node.items + i * BRANCH_ITEM_LEN + 16);
                self.walk_node(child, depth + 1, remaining, occupied, ancestors, visitor)?;
            }
        }
        Ok(())
    }

    /// Every leaf item into `out` (cleared first), in tree order: each data
    /// block with the span the index gives it. Validates the whole tree.
    pub fn leaves(&self, out: &mut Vec<Leaf>) -> Result<()> {
        out.clear();
        self.walk_nodes(&mut |node| {
            if !node.leaf {
                return Ok(());
            }
            for i in 0..node.count as usize {
                let at = node.items + i * LEAF_ITEM_LEN;
                let block = BlockRef {
                    offset: u64_at(self.buf, at + 16),
                    size: u64_at(self.buf, at + 24),
                };
                if block.offset.checked_add(block.size).is_none() {
                    return Err(invalid("data block", "file range overflows u64"));
                }
                out.push(Leaf {
                    start: (u32_at(self.buf, at), u32_at(self.buf, at + 4)),
                    end: (u32_at(self.buf, at + 8), u32_at(self.buf, at + 12)),
                    block,
                });
            }
            Ok(())
        })
    }

    /// Data blocks overlapping `[start, end)` on `chrom`, in file order.
    ///
    /// A block overlaps when its (chromosome, base) span intersects the query
    /// in the index's lexicographic order, so blocks spanning chromosomes
    /// (as bigBed blocks may) are found from either side.
    pub fn query(&self, chrom: ChromId, start: u32, end: u32) -> Query<'a> {
        let mut query = Query {
            index: *self,
            lo: (chrom, start),
            hi: (chrom, end),
            stack: [MaybeUninit::uninit(); MAX_DEPTH],
            depth: 0,
            remaining_bytes: self.buf.len() - INDEX_HEADER_LEN,
            done: start >= end,
            pending: None,
        };
        if !query.done
            && let Err(e) = query.push(self.base + INDEX_HEADER_LEN as u64)
        {
            query.done = true;
            query.pending = Some(e);
        }
        query
    }
}

#[derive(Debug, Clone, Copy)]
struct Frame {
    /// Buffer position of the node's first item.
    items: usize,
    count: u16,
    next: u16,
    leaf: bool,
    ordered_starts: bool,
}

/// Iterator over the blocks of [`Index::query`]. Yields an error once and
/// stops if the index is malformed.
#[derive(Debug)]
pub struct Query<'a> {
    index: Index<'a>,
    lo: Point,
    hi: Point,
    stack: [MaybeUninit<Frame>; MAX_DEPTH],
    depth: usize,
    remaining_bytes: usize,
    done: bool,
    pending: Option<Error>,
}

impl Query<'_> {
    // Keep raw validation out of prepared queries' node-selection path.
    #[inline(never)]
    fn read_raw_node(&mut self, offset: u64) -> Result<NodeInfo> {
        // Active frames are initialized by push and remain initialized until
        // popped. Detect cycles before charging the byte budget so they have
        // the same explicit depth error as preparation.
        for frame in &self.stack[..self.depth] {
            // SAFETY: push initializes every slot below depth; Frame is Copy.
            let items = unsafe { frame.assume_init_ref() }.items;
            // read_node validated items as the position immediately after
            // the node header, including its absolute file range.
            if self.index.base + items as u64 - NODE_HEADER_LEN as u64 == offset {
                return Err(Error::TooDeep(MAX_DEPTH));
            }
        }
        let node = self.index.read_node(offset)?;
        // A valid tree's nodes occupy disjoint encoded ranges. Charging all
        // record bytes, including repeated references, bounds raw traversal
        // by actual input bytes rather than node_count * maximum fanout.
        let encoded = NODE_HEADER_LEN
            + node.count as usize
                * if node.leaf {
                    LEAF_ITEM_LEN
                } else {
                    BRANCH_ITEM_LEN
                };
        self.remaining_bytes = self
            .remaining_bytes
            .checked_sub(encoded)
            .ok_or_else(|| invalid("index", "node references exceed encoded byte budget"))?;
        Ok(node)
    }

    fn push(&mut self, offset: u64) -> Result<()> {
        if self.depth == MAX_DEPTH {
            return Err(Error::TooDeep(MAX_DEPTH));
        }
        let node = if self.index.nodes.is_empty() {
            self.read_raw_node(offset)?
        } else if self.depth == 0 {
            // prepare requires the root at base + INDEX_HEADER_LEN and
            // rejects children overlapping it or the header. Sorting by
            // offset therefore puts that unique root first. Only query's
            // initial push has depth zero; descendants keep using lookup.
            self.index.nodes[0]
        } else {
            let i = self
                .index
                .nodes
                .binary_search_by_key(&offset, |node| node.offset)
                .map_err(|_| invalid("index workspace", "node not prepared"))?;
            self.index.nodes[i]
        };
        let range = node
            .order
            .candidates(node.count as usize, self.lo, self.hi, |i| {
                let at = node.items
                    + i * if node.leaf {
                        LEAF_ITEM_LEN
                    } else {
                        BRANCH_ITEM_LEN
                    };
                // SAFETY: read_node validated all node items before this metadata
                // was stored. candidates bounds i to the immutable node count.
                unsafe {
                    (
                        (
                            u32_unchecked(self.index.buf, at),
                            u32_unchecked(self.index.buf, at + 4),
                        ),
                        (
                            u32_unchecked(self.index.buf, at + 8),
                            u32_unchecked(self.index.buf, at + 12),
                        ),
                    )
                }
            });
        self.stack[self.depth].write(Frame {
            items: node.items,
            count: range.end as u16,
            next: range.start as u16,
            leaf: node.leaf,
            ordered_starts: node.order.starts,
        });
        self.depth += 1;
        Ok(())
    }

    fn step(&mut self) -> Result<Option<BlockRef>> {
        let buf = self.index.buf;
        while self.depth > 0 {
            // SAFETY: depth <= MAX_DEPTH and push writes each frame before
            // increasing depth. Pop only decreases depth; no uninitialized
            // slot is read. Frame is Copy and owns no resources to drop.
            let frame = unsafe {
                self.stack
                    .get_unchecked_mut(self.depth - 1)
                    .assume_init_mut()
            };
            if frame.next == frame.count {
                self.depth -= 1;
                continue;
            }
            let item = if frame.leaf {
                LEAF_ITEM_LEN
            } else {
                BRANCH_ITEM_LEN
            };
            let at = frame.items + frame.next as usize * item;
            frame.next += 1;
            // SAFETY: push validated count * item bytes for this frame.
            // next was less than count, so this complete item is in buf;
            // the coordinate fields fit both branch and leaf layouts.
            let (first, last) = unsafe {
                (
                    (u32_unchecked(buf, at), u32_unchecked(buf, at + 4)),
                    (u32_unchecked(buf, at + 8), u32_unchecked(buf, at + 12)),
                )
            };
            if first >= self.hi {
                // Only certified ordering permits stopping at the first miss.
                if frame.ordered_starts {
                    frame.next = frame.count;
                }
                continue;
            }
            if last <= self.lo {
                continue;
            }
            // SAFETY: The complete item was validated by push; both layouts
            // contain the eight-byte offset at byte 16.
            let offset = unsafe { u64_unchecked(buf, at + 16) };
            if frame.leaf {
                // SAFETY: push used LEAF_ITEM_LEN (32) for this frame.
                let size = unsafe { u64_unchecked(buf, at + 24) };
                if offset.checked_add(size).is_none() {
                    return Err(invalid("data block", "file range overflows u64"));
                }
                return Ok(Some(BlockRef { offset, size }));
            }
            self.push(offset)?;
        }
        Ok(None)
    }
}

impl FusedIterator for Query<'_> {}

impl Iterator for Query<'_> {
    type Item = Result<BlockRef>;

    fn next(&mut self) -> Option<Self::Item> {
        if let Some(e) = self.pending.take() {
            return Some(Err(e));
        }
        if self.done {
            return None;
        }
        match self.step() {
            Ok(Some(block)) => Some(Ok(block)),
            Ok(None) => {
                self.done = true;
                None
            }
            Err(e) => {
                self.done = true;
                Some(Err(e))
            }
        }
    }
}
