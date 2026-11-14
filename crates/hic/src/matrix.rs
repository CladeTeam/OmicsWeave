//! A chromosome pair's matrix at one resolution: its block grid and block
//! index, and which blocks a tile needs.

use std::ops::Range;

use omicsweave_core::{ContactUnit, Error, ReadAt, Result};

use crate::cursor::Cursor;
use crate::meta::Depth;

/// One stored block: read `size` bytes at `position`, then inflate.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub struct BlockEntry {
    pub number: u32,
    pub position: u64,
    pub size: u32,
}

impl BlockEntry {
    pub fn range(&self) -> Range<u64> {
        self.position..self.position + self.size as u64
    }
}

/// How a map is cut into blocks: square blocks of `block_bins` bins in a
/// grid `columns` blocks wide; within a chromosome from v9, blocks are laid
/// along the diagonal (`columns` positions) and by depth away from it.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Grid {
    pub block_bins: u32,
    pub columns: u32,
}

/// Bins per block side the Juicer tools aim for.
const BLOCK_SIZE: u64 = 1000;

/// Most block columns, ⌊√(2³¹ − 1)⌋ − 1, as the Juicer tools cap them:
/// rows × columns block numbers stay below 2³¹.
const MAX_COLUMNS: u64 = 46_339;

impl Grid {
    /// The grid the Juicer tools choose for a map whose longer chromosome
    /// spans `extent` units (base pairs or fragments).
    pub fn juicer(unit: ContactUnit, bin_size: u32, extent: u64, intra: bool) -> Grid {
        let bins = extent / bin_size as u64 + 1;
        // Fine base-pair maps get fewer, wider columns.
        let cutoff = match unit {
            ContactUnit::Frag => 0,
            ContactUnit::Bp if intra => 500,
            ContactUnit::Bp => 5000,
        };
        let mut columns = bins / BLOCK_SIZE + 1;
        if (bin_size as u64) < cutoff {
            columns = bins * bin_size as u64 / (BLOCK_SIZE * cutoff) + 1;
        }
        let columns = columns.min(MAX_COLUMNS);
        Grid {
            block_bins: (bins / columns + 1) as u32,
            columns: columns as u32,
        }
    }

    /// A grid of square blocks of `block_bins` bins per side for a map whose
    /// longer chromosome spans `extent` units at `bin_size`. Fails when the
    /// blocks are so small that block numbers would pass the format's
    /// signed 32-bit range (more than √(2³¹) columns).
    pub fn with_block_bins(bin_size: u32, extent: u64, block_bins: u32) -> Result<Grid> {
        if block_bins == 0 {
            return Err(Error::invalid("hic block grid", "zero bins per block"));
        }
        let bins = extent / bin_size as u64 + 1;
        let columns = bins / block_bins as u64 + 1;
        if columns > MAX_COLUMNS {
            return Err(Error::invalid(
                "hic block grid",
                format!("{block_bins}-bin blocks need {columns} columns, at most {MAX_COLUMNS}"),
            ));
        }
        Ok(Grid {
            block_bins,
            columns: columns as u32,
        })
    }

    /// The block holding the record at (`x`, `y`). Within a chromosome
    /// (`intra`) from v9, blocks follow `depth`; otherwise the square grid.
    pub fn block_number(&self, version: i32, depth: Depth, intra: bool, x: u32, y: u32) -> u64 {
        let bb = self.block_bins as i64;
        let columns = self.columns as u64;
        let (x, y) = (x as i64, y as i64);
        if version > 8 && intra {
            let pad = ((x + y) / 2 / bb) as u64;
            depth.of(x - y, bb) as u64 * columns + pad
        } else {
            (y / bb) as u64 * columns + (x / bb) as u64
        }
    }
}

/// A chromosome pair's matrix at one resolution. `chr1 <= chr2`; records'
/// x bins lie on `chr1`, y bins on `chr2`.
#[derive(Debug, Clone, PartialEq)]
pub struct ZoomData {
    pub chr1: u32,
    pub chr2: u32,
    pub unit: ContactUnit,
    pub bin_size: u32,
    /// Sum of the whole map (both halves of an intra-chromosomal one).
    pub sum_counts: f32,
    pub block_bin_count: u32,
    pub block_column_count: u32,
    /// Sorted by block number.
    pub blocks: Vec<BlockEntry>,
}

/// Inclusive ranges of block rows and columns (in grid or depth/pad space).
#[derive(Debug, Clone, Copy)]
struct Rect {
    rows: (i64, i64),
    cols: (i64, i64),
}

impl Rect {
    fn ranges(&self, columns: i64) -> u64 {
        if self.cols.1 - self.cols.0 + 1 >= columns {
            1
        } else {
            (self.rows.1 - self.rows.0 + 1) as u64
        }
    }

    fn contains(&self, number: u32, columns: i64) -> bool {
        let n = number as i64;
        let first = -((self.cols.1 - n).div_euclid(columns));
        let last = (n - self.cols.0).div_euclid(columns);
        first.max(self.rows.0) <= last.min(self.rows.1)
    }
}

impl ZoomData {
    /// Reads the matrix record at `position` and keeps the resolution
    /// `bin_size` in `unit`; `None` if the matrix has none.
    pub(crate) fn parse<S: ReadAt + ?Sized>(
        c: &mut Cursor<'_, S>,
        position: u64,
        unit: ContactUnit,
        bin_size: u32,
    ) -> Result<Option<ZoomData>> {
        let mut found = None;
        let mut first = true;
        ZoomData::parse_record(
            c,
            position,
            |u, b| {
                let hit = first && u == unit && b as i64 == bin_size as i64;
                first &= !hit;
                hit
            },
            &mut |z| found = Some(z),
        )?;
        Ok(found)
    }

    /// Reads the matrix record at `position` and hands every resolution
    /// `keep` accepts (by unit and bin size) to `out`, in stored order.
    /// Resolutions in units this crate does not know are skipped.
    pub(crate) fn parse_record<S: ReadAt + ?Sized>(
        c: &mut Cursor<'_, S>,
        position: u64,
        mut keep: impl FnMut(ContactUnit, i32) -> bool,
        out: &mut dyn FnMut(ZoomData),
    ) -> Result<()> {
        c.seek(position);
        let (Ok(chr1), Ok(chr2)) = (u32::try_from(c.i32()?), u32::try_from(c.i32()?)) else {
            return Err(Error::invalid("matrix", "negative chromosome index"));
        };
        for _ in 0..c.count("matrix resolutions")? {
            let this_unit = ContactUnit::parse(&c.string()?);
            c.i32()?; // resolution index
            let sum_counts = c.f32()?;
            c.bytes(12)?; // occupied cell count, 5th and 95th percentiles
            let size = c.i32()?;
            let block_bin_count = c.i32()?;
            let block_column_count = c.i32()?;
            let n = c.count("blocks")?;
            let Some(unit) = this_unit.filter(|&u| keep(u, size)) else {
                c.skip(n as u64 * 16)?;
                continue;
            };
            let (Ok(block_bin_count), Ok(block_column_count), Ok(bin_size)) = (
                u32::try_from(block_bin_count),
                u32::try_from(block_column_count),
                u32::try_from(size),
            ) else {
                return Err(Error::invalid("matrix", "negative block grid or bin size"));
            };
            if block_bin_count == 0 || block_column_count == 0 {
                return Err(Error::invalid("matrix", "zero block grid dimension"));
            }
            let mut blocks = (0..n)
                .map(|_| {
                    let number = c.i32()?;
                    let position = c.offset("block position")?;
                    let size = c.i32()?;
                    match (u32::try_from(number), u32::try_from(size)) {
                        (Ok(number), Ok(size)) => Ok(BlockEntry {
                            number,
                            position,
                            size,
                        }),
                        _ => Err(Error::invalid(
                            "block index",
                            "negative block number or size",
                        )),
                    }
                })
                .collect::<Result<Vec<_>>>()?;
            blocks.sort_unstable_by_key(|b| b.number);
            if blocks.windows(2).any(|w| w[0].number == w[1].number) {
                return Err(Error::invalid("block index", "duplicate block number"));
            }
            out(ZoomData {
                chr1,
                chr2,
                unit,
                bin_size,
                sum_counts,
                block_bin_count,
                block_column_count,
                blocks,
            });
        }
        Ok(())
    }

    pub fn grid(&self) -> Grid {
        Grid {
            block_bins: self.block_bin_count,
            columns: self.block_column_count,
        }
    }

    pub fn block(&self, number: u32) -> Option<&BlockEntry> {
        self.blocks
            .binary_search_by_key(&number, |b| b.number)
            .ok()
            .map(|i| &self.blocks[i])
    }

    /// Stored blocks that may hold records with x bins in `x` and y bins in
    /// `y` (and, for an intra-chromosomal map, the mirrored records), in
    /// block-number order: the candidates the Juicer tools read.
    pub fn blocks_for(
        &self,
        version: i32,
        depth: Depth,
        x: Range<u32>,
        y: Range<u32>,
    ) -> Vec<BlockEntry> {
        let mut out = Vec::new();
        self.blocks_for_into(version, depth, x, y, &mut out);
        out
    }

    /// [`ZoomData::blocks_for`] into caller-owned storage. Costs at most the
    /// size of the block index, however large the tile.
    pub fn blocks_for_into(
        &self,
        version: i32,
        depth: Depth,
        x: Range<u32>,
        y: Range<u32>,
        out: &mut Vec<BlockEntry>,
    ) {
        out.clear();
        if x.is_empty()
            || y.is_empty()
            || self.blocks.is_empty()
            || self.block_bin_count == 0
            || self.block_column_count == 0
        {
            return;
        }
        if matches!(depth, Depth::Log(0 | 1) | Depth::Linear(0)) {
            return;
        }
        let intra = self.chr1 == self.chr2;
        let (x1, x2) = (x.start as i64, x.end as i64 - 1);
        let (y1, y2) = (y.start as i64, y.end as i64 - 1);
        let bb = self.block_bin_count as i64;
        let columns = self.block_column_count as i64;
        let mut rects = [None, None];
        if version > 8 && intra {
            let (a, b) = (depth.of(x1 - y2, bb), depth.of(x2 - y1, bb));
            let crosses = x1 <= y2 && y1 <= x2;
            rects[0] = Some(Rect {
                rows: (if crosses { 0 } else { a.min(b) }, a.max(b) + 1),
                cols: ((x1 + y1) / 2 / bb, (x2 + y2) / 2 / bb + 1),
            });
        } else {
            let (cols, rows) = ((x1 / bb, (x2 + 1) / bb), (y1 / bb, (y2 + 1) / bb));
            rects[0] = Some(Rect { rows, cols });
            if intra {
                rects[1] = Some(Rect {
                    rows: cols,
                    cols: rows,
                });
            }
        }
        let rects = rects.iter().flatten();
        let scan = |out: &mut Vec<BlockEntry>| {
            out.clear();
            out.extend(
                self.blocks
                    .iter()
                    .filter(|b| rects.clone().any(|r| r.contains(b.number, columns))),
            );
        };
        // A wide rectangle's overlapping row ranges form one interval;
        // estimate the binary searches actually needed after that merge.
        let ranges: u64 = rects.clone().map(|r| r.ranges(columns)).sum();
        let seek_cost = self.blocks.len().ilog2() as u64 + 1;
        if ranges.saturating_mul(seek_cost) >= self.blocks.len() as u64 {
            scan(out);
            return;
        }
        let mirrored = rects.clone().count() > 1;
        let sort_budget = self.blocks.len() / seek_cost as usize;
        for r in rects.clone() {
            // Wide column bounds create overlapping numerical row ranges;
            // their union is one interval, so no duplicate copying is needed.
            let merged = r.cols.1 - r.cols.0 + 1 >= columns;
            let last_row = if merged { r.rows.0 } else { r.rows.1 };
            for row in r.rows.0..=last_row {
                let Some(base) = row.checked_mul(columns) else {
                    continue;
                };
                let Some(lo) = base.checked_add(r.cols.0) else {
                    continue;
                };
                let hi_base = if merged {
                    r.rows.1.saturating_mul(columns)
                } else {
                    base
                };
                let hi = hi_base.saturating_add(r.cols.1);
                let first = self.blocks.partition_point(|b| (b.number as i64) < lo);
                let last = first + self.blocks[first..].partition_point(|b| b.number as i64 <= hi);
                // Bound merging/sorting work by a full index scan as well.
                if mirrored && out.len() + last - first > sort_budget {
                    scan(out);
                    return;
                }
                out.extend_from_slice(&self.blocks[first..last]);
            }
        }
        if mirrored {
            out.sort_unstable_by_key(|b| b.number);
            out.dedup_by_key(|b| b.number);
        }
    }
}
