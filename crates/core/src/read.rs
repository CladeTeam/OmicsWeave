//! One access interface per data family. Readers of different formats in a
//! family are interchangeable: same key, same output layout.
//!
//! Fixed-size outputs are decoded straight into the caller's buffer.
//! Variable-size outputs are decoded once by a query, held by the reader
//! (outputs themselves are always the caller's memory), sized for the
//! caller, then copied out; nothing is decoded twice.

use crate::Result;
use crate::bins::BinSpec;
use crate::layout::{IntervalSize, IntervalsOut, TrackEntry};
use crate::locus::Locus;
use crate::region::Region;

/// Sequences (`.fa`, `.fa.gz`, `.2bit`).
pub trait SequenceRead {
    /// Length of `chrom`; [`crate::Error::UnknownChrom`] if absent.
    fn chrom_len(&self, chrom: &str) -> Result<u32>;

    /// Writes the bases of `locus` as ASCII into `out`
    /// (`out.len() == locus.len()`), 5' to 3' of its strand: a reverse
    /// locus reads the reverse complement of its concatenated spans. Every
    /// span must lie within the sequence. Each stored block the spans
    /// overlap is read and decoded once.
    fn sequence(&mut self, locus: &Locus<'_>, out: &mut [u8]) -> Result<()>;
}

/// Signal tracks (bigWig). Unknown chromosomes are errors; keys may run
/// past a chromosome's end, where there is no data; an empty key yields
/// no values, bins or entries.
pub trait TrackRead {
    /// Length of `chrom`; [`crate::Error::UnknownChrom`] if absent.
    fn chrom_len(&self, chrom: &str) -> Result<u32>;

    /// Writes one value per position of `locus` into `out`
    /// (`out.len() == locus.len()`), 5' to 3' of its strand; positions
    /// without data read NaN. Each data block the spans overlap is decoded
    /// once.
    fn values(&mut self, locus: &Locus<'_>, out: &mut [f32]) -> Result<()>;

    /// Writes the statistics of `spec` over the bins of `locus` into `out`:
    /// one row of [`BinSpec::count`] bins per statistic, in `spec.stats`
    /// order, bins 5' to 3' (`out.len() == spec.output_len(locus.len())`).
    /// See [`crate::BinSpec`] for the bin layout and statistics. The data
    /// are decoded once for all bins and statistics.
    fn bins(&mut self, locus: &Locus<'_>, spec: &BinSpec<'_>, out: &mut [f32]) -> Result<()>;

    /// Decodes the intervals with data in `region` once, clipped to it, in
    /// position order, and holds them in the reader until the next query.
    /// Returns how many there are, for the caller to allocate.
    fn query_entries(&mut self, region: Region<'_>) -> Result<usize>;

    /// Copies the intervals of the last [`TrackRead::query_entries`] into
    /// `out`, whose length must equal that count. Only copies: nothing is
    /// read or decoded again.
    fn write_entries(&self, out: &mut [TrackEntry]) -> Result<()>;
}

/// Interval features (bigBed). Entries are the features overlapping
/// the region, unclipped, ordered by start; an empty region, or a chromosome
/// without features (including one the source never names), yields none.
pub trait IntervalRead {
    /// Decodes the features of `region` once and holds them in the reader
    /// until the next query. Returns the output sizes, for the caller to
    /// allocate.
    fn query(&mut self, region: Region<'_>) -> Result<IntervalSize>;

    /// Copies the features of the last [`IntervalRead::query`] into `out`,
    /// sized exactly by it. Annotations are each feature's fields after the
    /// third, tab-separated. Only copies: nothing is read or decoded again.
    fn write(&self, out: IntervalsOut<'_>) -> Result<()>;
}
