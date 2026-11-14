//! Byte ranges of the variable-length sections, so the caller knows what to
//! read and how much memory to provide before parsing them.
//!
//! The file stores where sections start but not where they end. A section is
//! taken to end where the next known section starts (or at the end of the
//! file); parsers then bound-check every access within that range.

use std::ops::Range;

use crate::error::{Error, Result};
use crate::header::{HEADER_LEN, Header, TOTAL_SUMMARY_LEN, ZoomLevel};

/// Which R-tree a query runs against.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum Level {
    /// The full-resolution data.
    Base,
    /// The zoom level at this position of the zoom table.
    Zoom(usize),
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Layout {
    file_len: u64,
    /// Every known section start, sorted and deduplicated, ending with
    /// `file_len`.
    starts: Vec<u64>,
    chrom_tree: Range<u64>,
    data: Range<u64>,
    base_index: Range<u64>,
    auto_sql: Option<Range<u64>>,
    total_summary: Option<Range<u64>>,
    zoom_indexes: Vec<Range<u64>>,
}

impl Layout {
    pub fn new(header: &Header, zooms: &[ZoomLevel], file_len: u64) -> Result<Layout> {
        if zooms.len() != header.zoom_levels as usize {
            return Err(crate::error::invalid(
                "layout",
                format!(
                    "{} zoom headers for {} levels",
                    zooms.len(),
                    header.zoom_levels
                ),
            ));
        }
        let table_end = header.zoom_headers().end;
        let check = |what: &'static str, offset: u64| -> Result<u64> {
            if offset < table_end || offset >= file_len {
                return Err(Error::OutOfRange {
                    what,
                    offset,
                    lo: table_end,
                    hi: file_len,
                });
            }
            Ok(offset)
        };

        let mut starts = Vec::with_capacity(7 + 2 * zooms.len());
        starts.extend([
            check("chrom tree", header.chrom_tree_offset)?,
            check("data", header.full_data_offset)?,
            check("index", header.full_index_offset)?,
            file_len,
        ]);
        for (what, offset) in [
            ("autoSql", header.auto_sql_offset),
            ("total summary", header.total_summary_offset),
            ("extension", header.extension_offset),
        ] {
            if offset != 0 {
                starts.push(check(what, offset)?);
            }
        }
        for zoom in zooms {
            starts.push(check("zoom data", zoom.data_offset)?);
            starts.push(check("zoom index", zoom.index_offset)?);
        }
        starts.sort_unstable();
        starts.dedup();

        let section = |start: u64| -> Range<u64> {
            let i = starts.partition_point(|&s| s <= start);
            start..starts[i.min(starts.len() - 1)]
        };

        if !(header.chrom_tree_offset < header.full_data_offset
            && header.full_data_offset < header.full_index_offset)
        {
            return Err(crate::error::invalid(
                "layout",
                "expected chrom tree < data < index offsets",
            ));
        }

        let total_summary = match header.total_summary_offset {
            0 => None,
            start => {
                let limit = section(start).end;
                let Some(end) = start
                    .checked_add(TOTAL_SUMMARY_LEN as u64)
                    .filter(|&end| end <= limit)
                else {
                    return Err(Error::OutOfRange {
                        what: "total summary",
                        offset: start,
                        lo: HEADER_LEN as u64,
                        hi: limit,
                    });
                };
                Some(start..end)
            }
        };

        Ok(Layout {
            file_len,
            chrom_tree: section(header.chrom_tree_offset),
            data: header.full_data_offset..header.full_index_offset,
            base_index: section(header.full_index_offset),
            auto_sql: (header.auto_sql_offset != 0).then(|| section(header.auto_sql_offset)),
            total_summary,
            zoom_indexes: zooms.iter().map(|z| section(z.index_offset)).collect(),
            starts,
        })
    }

    pub fn file_len(&self) -> u64 {
        self.file_len
    }

    /// The chromosome B+ tree.
    pub fn chrom_tree(&self) -> Range<u64> {
        self.chrom_tree.clone()
    }

    /// The full-resolution data blocks (with their leading record count).
    pub fn data(&self) -> Range<u64> {
        self.data.clone()
    }

    /// The R-tree of a level, header included.
    pub fn index(&self, level: Level) -> Option<Range<u64>> {
        match level {
            Level::Base => Some(self.base_index.clone()),
            Level::Zoom(i) => self.zoom_indexes.get(i).cloned(),
        }
    }

    /// The autoSql string, NUL terminator included.
    pub fn auto_sql(&self) -> Option<Range<u64>> {
        self.auto_sql.clone()
    }

    pub fn total_summary(&self) -> Option<Range<u64>> {
        self.total_summary.clone()
    }

    /// Every known section start, sorted, ending with the file length.
    pub fn section_starts(&self) -> &[u64] {
        &self.starts
    }
}
