//! The contact family (Hi-C): tiles of raw contact counts from a binned,
//! symmetric contact map.

use std::ops::Range;

use crate::error::{Error, Result};
use crate::region::Region;

/// What a bin coordinate counts.
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq, Hash, PartialOrd, Ord)]
pub enum ContactUnit {
    /// Base pairs: a region is a base range.
    #[default]
    Bp,
    /// Restriction fragments: a region is a range of fragment indices
    /// (fragment `i` runs from restriction site `i - 1` to site `i`).
    Frag,
}

impl ContactUnit {
    /// The unit's name as `.hic` files spell it.
    pub fn as_str(self) -> &'static str {
        match self {
            ContactUnit::Bp => "BP",
            ContactUnit::Frag => "FRAG",
        }
    }

    pub fn parse(s: &str) -> Option<ContactUnit> {
        match s {
            "BP" => Some(ContactUnit::Bp),
            "FRAG" => Some(ContactUnit::Frag),
            _ => None,
        }
    }
}

/// A tile of a contact map: the bins of `rows` against the bins of `cols`
/// at `resolution` units per bin. The two regions may lie on different
/// chromosomes. A region covers every bin it overlaps, so the tile has
/// [`ContactQuery::shape`] cells.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct ContactQuery<'a> {
    pub rows: Region<'a>,
    pub cols: Region<'a>,
    pub resolution: u32,
    /// The unit of the regions and of `resolution`.
    pub unit: ContactUnit,
}

impl<'a> ContactQuery<'a> {
    /// A base-pair tile.
    pub fn new(rows: Region<'a>, cols: Region<'a>, resolution: u32) -> ContactQuery<'a> {
        ContactQuery {
            rows,
            cols,
            resolution,
            unit: ContactUnit::Bp,
        }
    }

    pub fn with_unit(self, unit: ContactUnit) -> ContactQuery<'a> {
        ContactQuery { unit, ..self }
    }

    /// Bins overlapped by `region` at `resolution`. Empty or reversed
    /// regions yield an empty range at their start bin.
    ///
    /// # Panics
    ///
    /// Panics if `resolution` is zero. Use [`ContactQuery::check`] to
    /// validate an untrusted query before computing its bins or shape.
    pub fn bins(region: &Region<'_>, resolution: u32) -> Range<u32> {
        assert!(resolution != 0, "contact resolution must be positive");
        if region.is_empty() {
            let at = region.start / resolution;
            return at..at;
        }
        region.start / resolution..region.end.div_ceil(resolution)
    }

    pub fn row_bins(&self) -> Range<u32> {
        Self::bins(&self.rows, self.resolution)
    }

    pub fn col_bins(&self) -> Range<u32> {
        Self::bins(&self.cols, self.resolution)
    }

    /// (rows, columns) of the tile.
    pub fn shape(&self) -> (usize, usize) {
        (self.row_bins().len(), self.col_bins().len())
    }

    /// Fails unless both regions are valid and the resolution positive.
    pub fn check(&self) -> Result<()> {
        self.rows.check()?;
        self.cols.check()?;
        if self.resolution == 0 {
            return Err(Error::invalid("contact query", "zero resolution"));
        }
        Ok(())
    }
}

/// One cell of a contact tile: the start (in the query's unit) of its row
/// bin and of its column bin, and its contact count.
///
/// numpy: `np.dtype([("row", "<u4"), ("col", "<u4"), ("value", "<f4")])`.
#[repr(C)]
#[derive(Debug, Clone, Copy, Default, PartialEq)]
pub struct ContactEntry {
    pub row: u32,
    pub col: u32,
    pub value: f32,
}

const _: () =
    assert!(std::mem::size_of::<ContactEntry>() == 12 && std::mem::align_of::<ContactEntry>() == 4);

/// Contact maps (`.hic`): raw contact counts, without normalization. Tiles
/// are given whole in the orientation asked for, whichever half of a
/// symmetric map the source stores. Unknown chromosomes are errors; regions
/// may run past a chromosome's end, where there are no contacts.
pub trait ContactRead {
    /// Length of `chrom` in base pairs; [`Error::UnknownChrom`] if absent.
    fn chrom_len(&self, chrom: &str) -> Result<u32>;

    /// Bin sizes available in `unit`.
    fn resolutions(&self, unit: ContactUnit) -> &[u32];

    /// Writes the tile row-major into `out` (`rows × cols` of
    /// [`ContactQuery::shape`]), 0 where there is no contact.
    fn matrix(&mut self, query: ContactQuery<'_>, out: &mut [f32]) -> Result<()>;

    /// Decodes the tile's cells that hold a contact once and keeps those with
    /// a nonzero count in the reader until the next query, ordered by (row,
    /// col). Returns how many there are.
    fn query_entries(&mut self, query: ContactQuery<'_>) -> Result<usize>;

    /// Copies the cells of the last [`ContactRead::query_entries`] into
    /// `out`, whose length must equal that count. Only copies.
    fn write_entries(&self, out: &mut [ContactEntry]) -> Result<()>;
}
