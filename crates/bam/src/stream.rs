//! BAM records read through a [`SeekReader`].

use std::ops::Range;

use omicsweave_bgzf::SeekReader;
use omicsweave_core::{Error, ReadAt, Result};

use crate::header::Take;

/// The next BAM record after its `block_size`, as a [`SeekReader::slice`]
/// range; `None` at a clean end of the file.
pub(crate) fn next_record<S: ReadAt + ?Sized>(
    stream: &mut SeekReader,
    source: &S,
) -> Result<Option<Range<usize>>> {
    let have = stream.fill(source, 4)?;
    if have == 0 {
        return Ok(None);
    }
    if have < 4 {
        return Err(Error::Truncated {
            what: "BAM record",
            need: 4,
            have: have as u64,
        });
    }
    let size = i32::from_le_bytes(stream.peek()[..4].try_into().unwrap());
    if size < 32 {
        return Err(Error::invalid(
            "BAM record",
            format!("block_size {size} below the 32 bytes of fixed fields"),
        ));
    }
    let total = stream.take(source, 4 + size as usize, "BAM record")?;
    Ok(Some(total.start + 4..total.end))
}

/// A [`SeekReader`] with its source, for reading a header.
pub(crate) struct Reading<'a, S: ?Sized> {
    pub stream: &'a mut SeekReader,
    pub source: &'a S,
}

impl<S: ReadAt + ?Sized> Take for Reading<'_, S> {
    fn take(&mut self, n: usize, what: &'static str) -> Result<&[u8]> {
        let range = self.stream.take(self.source, n, what)?;
        Ok(self.stream.slice(range))
    }
}
