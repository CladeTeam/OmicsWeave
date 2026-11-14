//! [`HicSidecar`]: every block of a `.hic` file with its decompressed size
//! and decoded record count. A `.hic` block index records only where a
//! block is and how many compressed bytes it takes; its zlib stream does not
//! record the inflated length, and the block encoding (list of rows or
//! dense, short or wide fields) is chosen per block, so neither size is
//! known before inflating. The sidecar measures each block once.
//!
//! # Body layout (version [`HIC_SIDECAR_VERSION`])
//!
//! Inside the container of [`omicsweave_core::sidecar`], little-endian:
//!
//! | Offset | Type | Field |
//! |---|---|---|
//! | 0 | `u32` | `.hic` version (8 or 9) |
//! | 4 | `u32` | reserved (0) |
//! | 8 | `u64` | matrices `M` |
//! | 16 | `M × 40` | matrices |
//! | | `Σ blocks × 24` | every matrix's blocks, matrix by matrix |
//!
//! A matrix record (sorted by chromosome pair, unit, bin size):
//!
//! | Offset | Type | Field |
//! |---|---|---|
//! | 0 | `u32` × 2 | chromosome indexes `chr1 <= chr2` |
//! | 8 | `u32` | unit: 0 BP, 1 FRAG |
//! | 12 | `u32` | bin size |
//! | 16 | `u32` × 2 | block bin count, block column count |
//! | 24 | `u64` | index of its first block in the block table |
//! | 32 | `u64` | its block count |
//!
//! A block record (sorted by block number within its matrix):
//!
//! | Offset | Type | Field |
//! |---|---|---|
//! | 0 | `u32` | block number |
//! | 4 | `u32` | compressed bytes |
//! | 8 | `u64` | file position |
//! | 16 | `u32` | bytes once inflated |
//! | 20 | `u32` | records [`crate::decode_block`] yields (12 bytes each as [`crate::Record`]) |
//!
//! A block's index in the block table is its ordinal.

use std::io::Write;
use std::ops::Range;
use std::path::Path;

use omicsweave_core::sidecar::{
    self, BodyReader, BodyWriter, SidecarKind, SourceIdentity, sidecar_path,
};
use omicsweave_core::{ContactUnit, Error, FileSource, ReadAt, Result};

use crate::block::{Inflater, decode_block};
use crate::cursor::Cursor;
use crate::matrix::ZoomData;
use crate::meta::{Footer, Header, MAX_VERSION, MIN_VERSION};

/// Body version this crate writes and reads.
pub const HIC_SIDECAR_VERSION: u32 = 1;

const MATRIX_LEN: usize = 40;
const BLOCK_LEN: usize = 24;

/// One stored block.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub struct BlockSize {
    pub number: u32,
    pub position: u64,
    /// Compressed bytes in the file.
    pub size: u32,
    /// Bytes once inflated.
    pub inflated: u32,
    /// Records the block decodes to.
    pub records: u32,
}

impl BlockSize {
    pub fn range(&self) -> Range<u64> {
        self.position..self.position + self.size as u64
    }
}

/// One chromosome pair's matrix at one resolution.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct MatrixBlocks {
    pub chr1: u32,
    pub chr2: u32,
    pub unit: ContactUnit,
    pub bin_size: u32,
    pub block_bin_count: u32,
    pub block_column_count: u32,
    /// Sorted by block number.
    pub blocks: Vec<BlockSize>,
}

impl MatrixBlocks {
    fn key(&self) -> (u32, u32, ContactUnit, u32) {
        (self.chr1, self.chr2, self.unit, self.bin_size)
    }

    pub fn block(&self, number: u32) -> Option<&BlockSize> {
        self.blocks
            .binary_search_by_key(&number, |b| b.number)
            .ok()
            .map(|i| &self.blocks[i])
    }
}

/// Every block of a `.hic` file, or of the matrices chosen when building.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct HicSidecar {
    pub version: i32,
    /// Sorted by (chr1, chr2, unit, bin size).
    pub matrices: Vec<MatrixBlocks>,
}

fn unit_code(u: ContactUnit) -> u32 {
    match u {
        ContactUnit::Bp => 0,
        ContactUnit::Frag => 1,
    }
}

impl HicSidecar {
    /// The matrix of a chromosome pair (`chr1 <= chr2`) at a resolution.
    pub fn matrix(
        &self,
        chr1: u32,
        chr2: u32,
        unit: ContactUnit,
        bin_size: u32,
    ) -> Option<&MatrixBlocks> {
        let key = (chr1, chr2, unit, bin_size);
        self.matrices
            .binary_search_by_key(&key, MatrixBlocks::key)
            .ok()
            .map(|i| &self.matrices[i])
    }

    /// Inflates and decodes every block of every matrix of `source` once.
    pub fn build<S: ReadAt + ?Sized>(source: &S) -> Result<HicSidecar> {
        HicSidecar::build_filtered(source, |_, _, _, _| true)
    }

    /// Like [`HicSidecar::build`], for the matrices `keep` accepts by
    /// chromosome pair, unit and bin size.
    pub fn build_filtered<S: ReadAt + ?Sized>(
        source: &S,
        mut keep: impl FnMut(u32, u32, ContactUnit, u32) -> bool,
    ) -> Result<HicSidecar> {
        let mut cursor = Cursor::new(source, 0);
        let header = Header::parse(&mut cursor)?;
        let footer = Footer::parse(&mut cursor, &header, source.len())?;
        let mut pairs: Vec<(&(u32, u32), &u64)> = footer.matrices.iter().collect();
        pairs.sort_unstable();
        let mut zooms = Vec::new();
        for (&(c1, c2), &position) in pairs {
            let mut found = Vec::new();
            ZoomData::parse_record(
                &mut Cursor::new(source, position),
                position,
                |unit, bin| u32::try_from(bin).is_ok_and(|bin| keep(c1, c2, unit, bin)),
                &mut |z| found.push(z),
            )?;
            if found.iter().any(|z| (z.chr1, z.chr2) != (c1, c2)) {
                return Err(Error::invalid(
                    "hic matrix",
                    "chromosome pair disagrees with master index",
                ));
            }
            zooms.extend(found);
        }
        let mut inflater = Inflater::new();
        let mut raw = Vec::new();
        let mut block = Vec::new();
        let mut matrices = Vec::with_capacity(zooms.len());
        for z in zooms {
            let mut blocks = Vec::with_capacity(z.blocks.len());
            let mut first = 0;
            while first < z.blocks.len() {
                let start = z.blocks[first].position;
                let end_of = |b: &crate::matrix::BlockEntry| {
                    b.position
                        .checked_add(b.size as u64)
                        .filter(|&end| end <= source.len())
                        .ok_or(Error::OutOfRange {
                            what: "hic block",
                            offset: b.position,
                            lo: 0,
                            hi: source.len(),
                        })
                };
                let mut end = end_of(&z.blocks[first])?;
                let mut last = first + 1;
                while last < z.blocks.len() && z.blocks[last].position == end {
                    let next = end_of(&z.blocks[last])?;
                    if next - start > 1 << 20 {
                        break;
                    }
                    end = next;
                    last += 1;
                }
                let len = usize::try_from(end - start)
                    .map_err(|_| Error::invalid("hic block", "read exceeds usize"))?;
                raw.resize(len, 0);
                source.read_exact_at(&mut raw, start)?;
                for b in &z.blocks[first..last] {
                    let lo = (b.position - start) as usize;
                    inflater.inflate(&raw[lo..lo + b.size as usize], &mut block)?;
                    let mut records = 0u64;
                    decode_block(header.version, &block, |_, _, _| records += 1)?;
                    let narrow = |v: u64| {
                        u32::try_from(v)
                            .map_err(|_| Error::invalid("hic block", "exceeds 2^32 - 1"))
                    };
                    blocks.push(BlockSize {
                        number: b.number,
                        position: b.position,
                        size: b.size,
                        inflated: narrow(block.len() as u64)?,
                        records: narrow(records)?,
                    });
                }
                first = last;
            }
            matrices.push(MatrixBlocks {
                chr1: z.chr1,
                chr2: z.chr2,
                unit: z.unit,
                bin_size: z.bin_size,
                block_bin_count: z.block_bin_count,
                block_column_count: z.block_column_count,
                blocks,
            });
        }
        matrices.sort_by_key(MatrixBlocks::key);
        if matrices.windows(2).any(|w| w[0].key() == w[1].key()) {
            return Err(Error::invalid("hic matrix", "resolution stored twice"));
        }
        Ok(HicSidecar {
            version: header.version,
            matrices,
        })
    }

    /// The body bytes, without the container.
    pub fn body(&self) -> Vec<u8> {
        let blocks: usize = self.matrices.iter().map(|m| m.blocks.len()).sum();
        let mut w = BodyWriter(Vec::with_capacity(
            16 + self.matrices.len() * MATRIX_LEN + blocks * BLOCK_LEN,
        ));
        w.u32(self.version as u32);
        w.u32(0);
        w.u64(self.matrices.len() as u64);
        let mut first = 0u64;
        for m in &self.matrices {
            w.u32(m.chr1);
            w.u32(m.chr2);
            w.u32(unit_code(m.unit));
            w.u32(m.bin_size);
            w.u32(m.block_bin_count);
            w.u32(m.block_column_count);
            w.u64(first);
            w.u64(m.blocks.len() as u64);
            first += m.blocks.len() as u64;
        }
        for b in self.matrices.iter().flat_map(|m| &m.blocks) {
            w.u32(b.number);
            w.u32(b.size);
            w.u64(b.position);
            w.u32(b.inflated);
            w.u32(b.records);
        }
        w.0
    }

    /// The whole sidecar, bound to the source with identity `source`.
    pub fn to_bytes(&self, source: SourceIdentity) -> Vec<u8> {
        sidecar::seal(SidecarKind::Hic, HIC_SIDECAR_VERSION, source, &self.body())
    }

    /// Writes the sidecar of `source` to `out`.
    pub fn write_to<S: ReadAt + ?Sized>(&self, source: &S, mut out: impl Write) -> Result<()> {
        out.write_all(&self.to_bytes(SourceIdentity::of(source)?))?;
        Ok(())
    }

    /// Parses a sidecar, checking that it belongs to the source with
    /// identity `source` ([`Error::StaleSidecar`] otherwise).
    pub fn parse(bytes: &[u8], source: SourceIdentity) -> Result<HicSidecar> {
        let (version, body) = sidecar::unseal(bytes, SidecarKind::Hic, source)?;
        if version != HIC_SIDECAR_VERSION {
            return Err(Error::invalid(
                "hic sidecar",
                format!("body version {version}, expected {HIC_SIDECAR_VERSION}"),
            ));
        }
        let bad = |detail: &str| Error::invalid("hic sidecar", detail.to_owned());
        let mut r = BodyReader::new(body, "hic sidecar");
        let hic_version = r.u32()? as i32;
        if !(MIN_VERSION..=MAX_VERSION).contains(&hic_version) {
            return Err(Error::Unsupported(".hic versions other than 8 and 9"));
        }
        r.u32()?;
        let n = r.count(MATRIX_LEN)?;
        let mut heads = Vec::with_capacity(n);
        let mut next = 0u64;
        for _ in 0..n {
            let chr1 = r.u32()?;
            let chr2 = r.u32()?;
            let unit = match r.u32()? {
                0 => ContactUnit::Bp,
                1 => ContactUnit::Frag,
                _ => return Err(bad("unknown unit")),
            };
            let bin_size = r.u32()?;
            let block_bin_count = r.u32()?;
            let block_column_count = r.u32()?;
            if chr1 > chr2 || bin_size == 0 || block_bin_count == 0 || block_column_count == 0 {
                return Err(bad("invalid chromosome pair, bin size or block grid"));
            }
            let first = r.u64()?;
            let count = r.u64()?;
            if first != next {
                return Err(bad("block ranges are not consecutive"));
            }
            next = first
                .checked_add(count)
                .ok_or_else(|| bad("block count overflows"))?;
            heads.push((
                chr1,
                chr2,
                unit,
                bin_size,
                block_bin_count,
                block_column_count,
                count,
            ));
        }
        if next.checked_mul(BLOCK_LEN as u64) != Some(r.remaining() as u64) {
            return Err(bad("block table length disagrees with the matrices"));
        }
        let mut matrices = Vec::with_capacity(n);
        for (chr1, chr2, unit, bin_size, block_bin_count, block_column_count, count) in heads {
            let mut blocks = Vec::with_capacity(count as usize);
            for _ in 0..count {
                let number = r.u32()?;
                let size = r.u32()?;
                let position = r.u64()?;
                let inflated = r.u32()?;
                let records = r.u32()?;
                if blocks
                    .last()
                    .is_some_and(|b: &BlockSize| b.number >= number)
                {
                    return Err(bad("blocks are not sorted with unique numbers"));
                }
                if position
                    .checked_add(size as u64)
                    .is_none_or(|end| end > source.len)
                {
                    return Err(bad("block range outside source"));
                }
                blocks.push(BlockSize {
                    number,
                    position,
                    size,
                    inflated,
                    records,
                });
            }
            matrices.push(MatrixBlocks {
                chr1,
                chr2,
                unit,
                bin_size,
                block_bin_count,
                block_column_count,
                blocks,
            });
        }
        r.finish()?;
        if matrices.windows(2).any(|w| w[0].key() >= w[1].key()) {
            return Err(bad("matrices are not sorted"));
        }
        Ok(HicSidecar {
            version: hic_version,
            matrices,
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
    pub fn load(path: impl AsRef<Path>) -> Result<HicSidecar> {
        let path = path.as_ref();
        let identity = SourceIdentity::of(&FileSource::open(path)?)?;
        HicSidecar::parse(&std::fs::read(sidecar_path(path))?, identity)
    }

    /// Builds the sidecar of the file at `path` and saves it beside it.
    pub fn build_file(path: impl AsRef<Path>) -> Result<HicSidecar> {
        let path = path.as_ref();
        let source = FileSource::open(path)?;
        let sidecar = HicSidecar::build(&source)?;
        let identity = SourceIdentity::of(&source)?;
        sidecar::write_file(&sidecar_path(path), &sidecar.to_bytes(identity))?;
        Ok(sidecar)
    }
}
