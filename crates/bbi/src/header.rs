//! The fixed header, the zoom-level table, the total summary and the autoSql
//! string: the metadata at the start of every BBI file.

use std::ops::Range;

use crate::bytes::{f64_at, require, u16_at, u32_at, u64_at};
use crate::error::{Error, Result, invalid};

pub const BIGWIG_MAGIC: u32 = 0x888F_FC26;
pub const BIGBED_MAGIC: u32 = 0x8789_F2EB;

/// Bytes of the fixed header at file offset 0.
pub const HEADER_LEN: usize = 64;
/// Bytes of one zoom-level header; the table follows the fixed header.
pub const ZOOM_HEADER_LEN: usize = 24;
/// Bytes of the total summary at [`Header::total_summary_offset`].
pub const TOTAL_SUMMARY_LEN: usize = 40;

/// Which BBI format a file holds; both share the container layout and differ
/// in the content of their data blocks.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum Kind {
    BigWig,
    BigBed,
}

/// The fixed 64-byte file header. Offsets are absolute file offsets; zero
/// means the section is absent.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Header {
    pub kind: Kind,
    pub version: u16,
    pub zoom_levels: u16,
    pub chrom_tree_offset: u64,
    pub full_data_offset: u64,
    pub full_index_offset: u64,
    /// Total BED fields per bigBed record (0 in bigWig).
    pub field_count: u16,
    /// Leading fields that follow the standard BED definitions.
    pub defined_field_count: u16,
    pub auto_sql_offset: u64,
    pub total_summary_offset: u64,
    /// Upper bound of a decompressed data block; zero means blocks are stored
    /// uncompressed.
    pub uncompress_buf_size: u32,
    pub extension_offset: u64,
}

impl Header {
    /// Parses the first [`HEADER_LEN`] bytes of a file.
    pub fn parse(buf: &[u8]) -> Result<Header> {
        require(buf, HEADER_LEN, "header")?;
        let magic = u32_at(buf, 0);
        let kind = match magic {
            BIGWIG_MAGIC => Kind::BigWig,
            BIGBED_MAGIC => Kind::BigBed,
            m if m.swap_bytes() == BIGWIG_MAGIC || m.swap_bytes() == BIGBED_MAGIC => {
                return Err(Error::BigEndian);
            }
            found => {
                return Err(Error::BadMagic {
                    what: "header",
                    found,
                });
            }
        };
        Ok(Header {
            kind,
            version: u16_at(buf, 4),
            zoom_levels: u16_at(buf, 6),
            chrom_tree_offset: u64_at(buf, 8),
            full_data_offset: u64_at(buf, 16),
            full_index_offset: u64_at(buf, 24),
            field_count: u16_at(buf, 32),
            defined_field_count: u16_at(buf, 34),
            auto_sql_offset: u64_at(buf, 36),
            total_summary_offset: u64_at(buf, 44),
            uncompress_buf_size: u32_at(buf, 52),
            extension_offset: u64_at(buf, 56),
        })
    }

    /// File range of the zoom-level table, which directly follows the header.
    pub fn zoom_headers(&self) -> Range<u64> {
        let start = HEADER_LEN as u64;
        start..start + self.zoom_levels as u64 * ZOOM_HEADER_LEN as u64
    }

    /// Data blocks are zlib streams that must be inflated before parsing.
    pub fn compressed(&self) -> bool {
        self.uncompress_buf_size != 0
    }

    /// Bytes an inflate output buffer needs to hold any decompressed block of
    /// this file (base or zoom). Zero for uncompressed files, whose blocks are
    /// parsed as read.
    pub fn block_buffer_len(&self) -> usize {
        self.uncompress_buf_size as usize
    }
}

/// One zoom level: summaries over windows of `reduction` bases, with their own
/// data blocks and R-tree index.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct ZoomLevel {
    pub reduction: u32,
    pub data_offset: u64,
    pub index_offset: u64,
}

impl ZoomLevel {
    /// Parses the table at [`Header::zoom_headers`]; `buf` holds exactly
    /// that range.
    pub fn parse_all(header: &Header, buf: &[u8]) -> Result<Vec<ZoomLevel>> {
        let need = header.zoom_levels as usize * ZOOM_HEADER_LEN;
        require(buf, need, "zoom headers")?;
        Ok((0..header.zoom_levels as usize)
            .map(|i| {
                let at = i * ZOOM_HEADER_LEN;
                ZoomLevel {
                    reduction: u32_at(buf, at),
                    data_offset: u64_at(buf, at + 8),
                    index_offset: u64_at(buf, at + 16),
                }
            })
            .collect())
    }
}

/// Whole-file statistics over every covered base.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct TotalSummary {
    pub bases_covered: u64,
    pub min: f64,
    pub max: f64,
    pub sum: f64,
    pub sum_squares: f64,
}

impl TotalSummary {
    /// Parses the [`TOTAL_SUMMARY_LEN`] bytes at
    /// [`Header::total_summary_offset`].
    pub fn parse(buf: &[u8]) -> Result<TotalSummary> {
        require(buf, TOTAL_SUMMARY_LEN, "total summary")?;
        Ok(TotalSummary {
            bases_covered: u64_at(buf, 0),
            min: f64_at(buf, 8),
            max: f64_at(buf, 16),
            sum: f64_at(buf, 24),
            sum_squares: f64_at(buf, 32),
        })
    }
}

/// The NUL-terminated autoSql schema of a bigBed file, from the bytes of
/// [`crate::Layout::auto_sql`].
pub fn parse_auto_sql(buf: &[u8]) -> Result<&str> {
    let end = memchr::memchr(0, buf).ok_or_else(|| invalid("autoSql", "missing NUL terminator"))?;
    std::str::from_utf8(&buf[..end]).map_err(|e| invalid("autoSql", e.to_string()))
}
