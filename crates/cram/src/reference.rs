//! Reference sequences for decoding, from any [`SequenceRead`].

use omicsweave_core::{Error, Locus, Region, Result, SequenceRead};

use crate::header::Reference;
use crate::slice::References;

/// A [`SequenceRead`] that has no sequences, for files read without a
/// reference: queries that need reference bases fail.
#[derive(Debug, Clone, Copy, Default)]
pub struct NoReference;

impl SequenceRead for NoReference {
    fn chrom_len(&self, chrom: &str) -> Result<u32> {
        Err(Error::UnknownChrom(chrom.to_owned()))
    }

    fn sequence(&mut self, locus: &Locus<'_>, _out: &mut [u8]) -> Result<()> {
        Err(Error::UnknownChrom(locus.chrom.to_owned()))
    }
}

/// The references of a CRAM header read from a [`SequenceRead`] by name,
/// uppercased as htslib loads them.
pub(crate) struct ByName<'a, R: ?Sized> {
    pub reader: &'a mut R,
    pub refs: &'a [Reference],
}

impl<R: SequenceRead + ?Sized> References for ByName<'_, R> {
    fn ref_len(&mut self, tid: usize) -> Result<i64> {
        Ok(self.reader.chrom_len(&self.refs[tid].name)? as i64)
    }

    fn bases(&mut self, tid: usize, start: i64, end: i64, out: &mut Vec<u8>) -> Result<()> {
        out.clear();
        if end < start {
            return Ok(());
        }
        let name = &self.refs[tid].name;
        let (lo, hi) = (start - 1, end);
        let (lo, hi) = (
            u32::try_from(lo)
                .map_err(|_| Error::invalid("CRAM reference", "position out of range"))?,
            u32::try_from(hi)
                .map_err(|_| Error::invalid("CRAM reference", "position out of range"))?,
        );
        out.resize((hi - lo) as usize, 0);
        self.reader
            .sequence(&Locus::from(Region::new(name, lo, hi)), out)?;
        out.make_ascii_uppercase();
        Ok(())
    }
}
