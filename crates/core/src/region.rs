use crate::error::{Error, Result};

/// The access key of every genomic source: bases `[start, end)` of a
/// chromosome, 0-based, half-open.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub struct Region<'a> {
    pub chrom: &'a str,
    pub start: u32,
    pub end: u32,
}

impl<'a> Region<'a> {
    pub fn new(chrom: &'a str, start: u32, end: u32) -> Region<'a> {
        Region { chrom, start, end }
    }

    /// Bases covered; zero for an invalid region.
    pub fn len(&self) -> u32 {
        self.end.saturating_sub(self.start)
    }

    pub fn is_empty(&self) -> bool {
        self.len() == 0
    }

    /// Fails unless `start <= end`.
    pub fn check(&self) -> Result<()> {
        if self.start > self.end {
            return Err(Error::InvalidRegion {
                start: self.start,
                end: self.end,
            });
        }
        Ok(())
    }

    /// Fails unless the region is valid and lies within a sequence of `len`
    /// bases.
    pub fn check_within(&self, len: u32) -> Result<()> {
        self.check()?;
        if self.end > len {
            return Err(Error::RegionOutOfRange {
                start: self.start,
                end: self.end,
                len,
            });
        }
        Ok(())
    }

    /// `[start, end)` intersects this region. Empty records count when they
    /// lie strictly inside, as in UCSC and libBigWig queries.
    pub fn overlaps(&self, start: u32, end: u32) -> bool {
        self.start < self.end && end > self.start && start < self.end
    }
}
