//! [`BbiSidecar`]: every data block of a bigWig or bigBed file, base level
//! and zoom levels, with its span, file range, decompressed size and item
//! count. BBI indexes record only a block's compressed size; the sidecar
//! adds what a loader must know to allocate before inflating.
//!
//! # Body layout (version [`BBI_SIDECAR_VERSION`])
//!
//! Inside the container of [`omicsweave_core::sidecar`], little-endian:
//!
//! | Offset | Type | Field |
//! |---|---|---|
//! | 0 | `u8` | kind: 0 bigWig, 1 bigBed |
//! | 1 | `u8` | reserved (0) |
//! | 2 | `u16` | levels `L`: the base level, then each zoom level |
//! | 4 | `u32` | the header's `uncompressBufSize` (0: uncompressed blocks) |
//! | 8 | `L × 16` | per level: `u32` reduction (0 for the base level), `u32` reserved, `u64` block count |
//! | | `Σ blocks × 40` | every level's blocks, level by level |
//!
//! A block record:
//!
//! | Offset | Type | Field |
//! |---|---|---|
//! | 0 | `u32` × 4 | start chromosome, start base, end chromosome, end base (the index's span) |
//! | 16 | `u64` | file offset |
//! | 24 | `u32` | bytes in the file |
//! | 28 | `u32` | bytes once inflated (the stored size when uncompressed) |
//! | 32 | `u32` | items: wig items, bigBed records or zoom summaries |
//! | 36 | `u32` | reserved (0) |
//!
//! Blocks keep the index's leaf order; a block's position in its level is
//! its ordinal.

use std::io::Write;
use std::ops::Range;
use std::path::Path;

use omicsweave_core::sidecar::{
    self, BodyReader, BodyWriter, SidecarKind, SourceIdentity, sidecar_path,
};
use omicsweave_core::{FileSource, ReadAt};

use crate::block::{BedBlock, WigBlock, ZoomBlock};
use crate::chroms::ChromId;
use crate::error::{Error, Result, invalid};
use crate::header::{HEADER_LEN, Header, Kind, ZoomLevel};
use crate::index::{Index, Leaf};
use crate::inflate::Inflater;
use crate::layout::{Layout, Level};

/// Body version this crate writes and reads.
pub const BBI_SIDECAR_VERSION: u32 = 1;

const LEVEL_LEN: usize = 16;
const BLOCK_LEN: usize = 40;
/// Contiguous blocks read together while building.
const READ_BATCH: u64 = 8 << 20;
/// A modest declaration can size the first real block without repeatedly
/// decompressing it; larger declarations use the inflater's geometric growth.
const MAX_INITIAL_INFLATE_HINT: usize = 256 << 10;

/// One data block.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub struct BlockInfo {
    /// The (chromosome, base) span the index records, end exclusive.
    pub start: (ChromId, u32),
    pub end: (ChromId, u32),
    pub offset: u64,
    /// Bytes in the file.
    pub size: u32,
    /// Bytes once inflated; equal to `size` in an uncompressed file.
    pub inflated: u32,
    /// Wig items, bigBed records or zoom summaries in the block.
    pub items: u32,
}

impl BlockInfo {
    /// The block's bytes in the file.
    ///
    /// # Panics
    /// Panics if a caller-constructed block's end overflows u64. Parsed
    /// sidecars and blocks built from a file have already passed this check.
    pub fn range(&self) -> Range<u64> {
        let end = self
            .offset
            .checked_add(self.size as u64)
            .expect("block range overflows u64");
        self.offset..end
    }
}

/// The blocks of one level, in index order.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct LevelBlocks {
    /// Bases per zoom summary; 0 for the base level.
    pub reduction: u32,
    pub blocks: Vec<BlockInfo>,
}

/// Every block of a BBI file.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct BbiSidecar {
    pub kind: Kind,
    /// The header's `uncompressBufSize`; 0 when blocks are not compressed.
    pub uncompress_buf_size: u32,
    /// The base level, then each zoom level in the file's order.
    pub levels: Vec<LevelBlocks>,
}

fn read<S: ReadAt + ?Sized>(source: &S, range: Range<u64>, buf: &mut Vec<u8>) -> Result<()> {
    if range.start > range.end || range.end > source.len() {
        return Err(Error::OutOfRange {
            what: "BBI read",
            offset: range.start,
            lo: 0,
            hi: source.len(),
        });
    }
    let len = usize::try_from(range.end - range.start)
        .map_err(|_| invalid("BBI read", "range exceeds usize"))?;
    buf.clear();
    buf.try_reserve(len)
        .map_err(|e| invalid("BBI buffer", e.to_string()))?;
    buf.resize(len, 0);
    source.read_exact_at(buf, range.start)
}

fn narrow(v: usize, what: &'static str) -> Result<u32> {
    u32::try_from(v).map_err(|_| invalid(what, "exceeds 2^32 - 1"))
}

impl BbiSidecar {
    /// Blocks are compressed.
    pub fn compressed(&self) -> bool {
        self.uncompress_buf_size != 0
    }

    /// The blocks of `level`.
    pub fn level(&self, level: Level) -> Option<&LevelBlocks> {
        match level {
            Level::Base => self.levels.first(),
            Level::Zoom(i) => i.checked_add(1).and_then(|i| self.levels.get(i)),
        }
    }

    /// Reads every index of `source` and inflates every block once to
    /// measure it. Contiguous blocks are read in batches.
    pub fn build<S: ReadAt + ?Sized>(source: &S) -> Result<BbiSidecar> {
        let mut buf = Vec::new();
        read(source, 0..HEADER_LEN as u64, &mut buf)?;
        let header = Header::parse(&buf)?;
        read(source, header.zoom_headers(), &mut buf)?;
        let zooms = ZoomLevel::parse_all(&header, &buf)?;
        let layout = Layout::new(&header, &zooms, source.len())?;

        let mut inflater = Inflater::new();
        let mut inflated = Vec::new();
        let mut leaves = Vec::new();
        let mut raw = Vec::new();
        let mut levels = Vec::with_capacity(1 + zooms.len());
        for i in 0..=zooms.len() {
            let (level, reduction) = match i {
                0 => (Level::Base, 0),
                _ => (Level::Zoom(i - 1), zooms[i - 1].reduction),
            };
            let range = layout.index(level).expect("every level has an index");
            read(source, range.clone(), &mut buf)?;
            Index::new(&buf, range.start)?.leaves(&mut leaves)?;
            let mut blocks = Vec::with_capacity(leaves.len());
            let mut k = 0;
            while k < leaves.len() {
                let run = batch(&leaves, k, source.len())?;
                read(source, run.file.clone(), &mut raw)?;
                for leaf in &leaves[run.leaves.clone()] {
                    let lo = (leaf.block.offset - run.file.start) as usize;
                    let bytes = &raw[lo..lo + leaf.block.size as usize];
                    let data = if header.compressed() {
                        if inflated.is_empty()
                            && header.block_buffer_len() <= MAX_INITIAL_INFLATE_HINT
                        {
                            inflated.resize(header.block_buffer_len(), 0);
                        }
                        let n = inflater.inflate_growing(
                            bytes,
                            &mut inflated,
                            0,
                            header.block_buffer_len(),
                        )?;
                        &inflated[..n]
                    } else {
                        bytes
                    };
                    let items = match (level, header.kind) {
                        (Level::Zoom(_), _) => ZoomBlock::parse(data)?.len(),
                        (Level::Base, Kind::BigWig) => WigBlock::parse(data)?.len(),
                        (Level::Base, Kind::BigBed) => {
                            let mut n = 0;
                            for entry in BedBlock::parse(data)?.entries() {
                                entry?;
                                n += 1;
                            }
                            n
                        }
                    };
                    blocks.push(BlockInfo {
                        start: leaf.start,
                        end: leaf.end,
                        offset: leaf.block.offset,
                        size: narrow(leaf.block.size as usize, "BBI block size")?,
                        inflated: narrow(data.len(), "BBI block size")?,
                        items: narrow(items, "BBI block items")?,
                    });
                }
                k = run.leaves.end;
            }
            levels.push(LevelBlocks { reduction, blocks });
        }
        Ok(BbiSidecar {
            kind: header.kind,
            uncompress_buf_size: header.uncompress_buf_size,
            levels,
        })
    }

    /// The body bytes, without the container.
    ///
    /// # Panics
    /// Panics if the number of levels, including the base level, exceeds
    /// 65535, the maximum representable in sidecar version 1.
    pub fn body(&self) -> Vec<u8> {
        let levels = u16::try_from(self.levels.len())
            .expect("level count exceeds sidecar version 1's u16 limit");
        let blocks: usize = self.levels.iter().map(|l| l.blocks.len()).sum();
        let mut w = BodyWriter(Vec::with_capacity(
            8 + self.levels.len() * LEVEL_LEN + blocks * BLOCK_LEN,
        ));
        w.u8(match self.kind {
            Kind::BigWig => 0,
            Kind::BigBed => 1,
        });
        w.u8(0);
        w.u16(levels);
        w.u32(self.uncompress_buf_size);
        for level in &self.levels {
            w.u32(level.reduction);
            w.u32(0);
            w.u64(level.blocks.len() as u64);
        }
        for b in self.levels.iter().flat_map(|l| &l.blocks) {
            w.u32(b.start.0);
            w.u32(b.start.1);
            w.u32(b.end.0);
            w.u32(b.end.1);
            w.u64(b.offset);
            w.u32(b.size);
            w.u32(b.inflated);
            w.u32(b.items);
            w.u32(0);
        }
        w.0
    }

    /// The whole sidecar, bound to the source with identity `source`.
    ///
    /// # Panics
    /// Panics when [`Self::body`] cannot represent the number of levels.
    pub fn to_bytes(&self, source: SourceIdentity) -> Vec<u8> {
        sidecar::seal(SidecarKind::Bbi, BBI_SIDECAR_VERSION, source, &self.body())
    }

    /// Writes the sidecar of `source` to `out`.
    pub fn write_to<S: ReadAt + ?Sized>(&self, source: &S, mut out: impl Write) -> Result<()> {
        out.write_all(&self.to_bytes(SourceIdentity::of(source)?))?;
        Ok(())
    }

    /// Parses a sidecar, checking that it belongs to the source with
    /// identity `source` ([`Error::StaleSidecar`] otherwise).
    pub fn parse(bytes: &[u8], source: SourceIdentity) -> Result<BbiSidecar> {
        let (version, body) = sidecar::unseal(bytes, SidecarKind::Bbi, source)?;
        if version != BBI_SIDECAR_VERSION {
            return Err(invalid(
                "BBI sidecar",
                format!("body version {version}, expected {BBI_SIDECAR_VERSION}"),
            ));
        }
        let mut r = BodyReader::new(body, "BBI sidecar");
        let kind = match r.u8()? {
            0 => Kind::BigWig,
            1 => Kind::BigBed,
            k => return Err(invalid("BBI sidecar", format!("unknown kind {k}"))),
        };
        r.u8()?;
        let n = r.u16()? as usize;
        if n == 0 {
            return Err(invalid("BBI sidecar", "no base level"));
        }
        let uncompress_buf_size = r.u32()?;
        if n > r.remaining() / LEVEL_LEN {
            return Err(invalid(
                "BBI sidecar",
                format!("{n} level headers cannot fit"),
            ));
        }
        let mut heads = Vec::with_capacity(n);
        for _ in 0..n {
            let reduction = r.u32()?;
            r.u32()?;
            heads.push((reduction, r.count(0)?));
        }
        let mut levels = Vec::with_capacity(n);
        for (reduction, count) in heads {
            if count > r.remaining() / BLOCK_LEN {
                return Err(invalid("BBI sidecar", format!("{count} blocks cannot fit")));
            }
            let mut blocks = Vec::with_capacity(count);
            for _ in 0..count {
                let start = (r.u32()?, r.u32()?);
                let end = (r.u32()?, r.u32()?);
                let offset = r.u64()?;
                let size = r.u32()?;
                let inflated = r.u32()?;
                let items = r.u32()?;
                r.u32()?;
                if offset
                    .checked_add(size as u64)
                    .is_none_or(|end| end > source.len)
                {
                    return Err(Error::OutOfRange {
                        what: "BBI sidecar block",
                        offset,
                        lo: 0,
                        hi: source.len,
                    });
                }
                blocks.push(BlockInfo {
                    start,
                    end,
                    offset,
                    size,
                    inflated,
                    items,
                });
            }
            levels.push(LevelBlocks { reduction, blocks });
        }
        r.finish()?;
        Ok(BbiSidecar {
            kind,
            uncompress_buf_size,
            levels,
        })
    }

    /// Writes the sidecar of the file at `path` to `<path>.sidecar`.
    pub fn save(&self, path: impl AsRef<Path>) -> Result<()> {
        let path = path.as_ref();
        let identity = SourceIdentity::of(&FileSource::open(path)?)?;
        sidecar::write_file(&sidecar_path(path), &self.to_bytes(identity))
    }

    /// Reads `<path>.sidecar`, checking that it belongs to the file at
    /// `path`.
    pub fn load(path: impl AsRef<Path>) -> Result<BbiSidecar> {
        let path = path.as_ref();
        let identity = SourceIdentity::of(&FileSource::open(path)?)?;
        BbiSidecar::parse(&std::fs::read(sidecar_path(path))?, identity)
    }

    /// Builds the sidecar of the file at `path` and saves it beside it.
    pub fn build_file(path: impl AsRef<Path>) -> Result<BbiSidecar> {
        let path = path.as_ref();
        let source = FileSource::open(path)?;
        let sidecar = BbiSidecar::build(&source)?;
        let identity = SourceIdentity::of(&source)?;
        sidecar::write_file(&sidecar_path(path), &sidecar.to_bytes(identity))?;
        Ok(sidecar)
    }
}

struct Run {
    leaves: Range<usize>,
    file: Range<u64>,
}

/// Leaves from `k` whose blocks follow each other in the file, up to
/// [`READ_BATCH`] bytes (always at least leaf `k`).
fn batch(leaves: &[Leaf], k: usize, file_len: u64) -> Result<Run> {
    let end_of = |leaf: &Leaf| -> Result<u64> {
        let b = leaf.block;
        b.offset
            .checked_add(b.size)
            .filter(|&e| e <= file_len)
            .ok_or(Error::OutOfRange {
                what: "BBI data block",
                offset: b.offset,
                lo: 0,
                hi: file_len,
            })
    };
    let start = leaves[k].block.offset;
    let mut end = end_of(&leaves[k])?;
    let mut next = k + 1;
    while let Some(leaf) = leaves.get(next) {
        if leaf.block.offset != end {
            break;
        }
        let e = end_of(leaf)?;
        if e - start > READ_BATCH {
            break;
        }
        end = e;
        next += 1;
    }
    Ok(Run {
        leaves: k..next,
        file: start..end,
    })
}
