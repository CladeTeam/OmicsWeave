//! Standard output layouts. Every reader of a family writes the same layout
//! into caller-allocated memory, so outputs map directly onto numpy arrays
//! (and Arrow buffers) without conversion.

use crate::error::{Error, Result, check_output};

/// One track interval, clipped to the queried region.
///
/// numpy: `np.dtype([("start", "<u4"), ("end", "<u4"), ("value", "<f4")])`.
#[repr(C)]
#[derive(Debug, Clone, Copy, Default, PartialEq)]
pub struct TrackEntry {
    pub start: u32,
    pub end: u32,
    pub value: f32,
}

/// One interval feature, as stored (not clipped: the feature and its
/// annotation describe the whole interval). `annotation` is the row of its
/// annotation in the accompanying [`StringColumn`].
///
/// numpy: `np.dtype([("start", "<u4"), ("end", "<u4"), ("annotation", "<u4")])`.
#[repr(C)]
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq, Hash)]
pub struct IntervalEntry {
    pub start: u32,
    pub end: u32,
    pub annotation: u32,
}

const _: () =
    assert!(std::mem::size_of::<TrackEntry>() == 12 && std::mem::align_of::<TrackEntry>() == 4);
const _: () = assert!(
    std::mem::size_of::<IntervalEntry>() == 12 && std::mem::align_of::<IntervalEntry>() == 4
);

/// A caller-allocated string column in Arrow `LargeString` layout: string
/// `i` is `data[offsets[i]..offsets[i + 1]]`; `offsets` has one more element
/// than there are strings and starts at 0.
pub struct StringColumn<'a> {
    pub offsets: &'a mut [i64],
    pub data: &'a mut [u8],
}

/// Read access to a string column in the same layout.
#[derive(Debug, Clone, Copy)]
pub struct Strings<'a> {
    pub offsets: &'a [i64],
    pub data: &'a [u8],
}

impl<'a> Strings<'a> {
    pub fn len(&self) -> usize {
        self.offsets.len().saturating_sub(1)
    }

    pub fn is_empty(&self) -> bool {
        self.len() == 0
    }

    /// The bytes of row `i`, without UTF-8 validation.
    ///
    /// # Panics
    ///
    /// Panics if the row is absent or its offsets are negative, reversed,
    /// not representable as `usize`, or outside `data`.
    pub fn get(&self, i: usize) -> &'a [u8] {
        let next = i.checked_add(1).expect("string row index overflows usize");
        let start = usize::try_from(self.offsets[i]).expect("invalid string start offset");
        let end = usize::try_from(self.offsets[next]).expect("invalid string end offset");
        &self.data[start..end]
    }
}

/// Sizes of an interval query's output, to allocate before filling.
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub struct IntervalSize {
    pub entries: usize,
    /// Total bytes of all annotations.
    pub annotation_bytes: usize,
}

/// Caller-allocated output of an interval query, sized by [`IntervalSize`].
pub struct IntervalsOut<'a> {
    pub entries: &'a mut [IntervalEntry],
    pub annotations: StringColumn<'a>,
}

impl IntervalsOut<'_> {
    /// Fails unless the buffers match `size` exactly and annotation indices
    /// and byte offsets fit their `u32` and `i64` output fields.
    pub fn check(&self, size: IntervalSize) -> Result<()> {
        if u32::try_from(size.entries.saturating_sub(1)).is_err() {
            return Err(Error::invalid(
                "interval size",
                "annotation row index exceeds u32",
            ));
        }
        if i64::try_from(size.annotation_bytes).is_err() {
            return Err(Error::invalid(
                "interval size",
                "annotation byte offset exceeds i64",
            ));
        }
        check_output("interval entries", self.entries.len(), size.entries)?;
        check_output(
            "annotation offsets",
            self.annotations.offsets.len(),
            size.entries
                .checked_add(1)
                .ok_or_else(|| Error::invalid("interval size", "offset count overflows usize"))?,
        )?;
        check_output(
            "annotation bytes",
            self.annotations.data.len(),
            size.annotation_bytes,
        )
    }
}

/// Appends entries and annotations to an [`IntervalsOut`] in order, so every
/// reader fills it the same way.
pub struct IntervalWriter<'o, 'a> {
    out: &'o mut IntervalsOut<'a>,
    row: usize,
    byte: usize,
}

impl<'o, 'a> IntervalWriter<'o, 'a> {
    /// Starts filling `out`, which must have been checked against its size.
    pub fn new(out: &'o mut IntervalsOut<'a>) -> IntervalWriter<'o, 'a> {
        if let Some(first) = out.annotations.offsets.first_mut() {
            *first = 0;
        }
        IntervalWriter {
            out,
            row: 0,
            byte: 0,
        }
    }

    /// Writes the next interval; its annotation becomes the next string row.
    /// Panics if the output is full or its annotation index exceeds `u32`,
    /// which a correct checked size prevents.
    pub fn push(&mut self, start: u32, end: u32, annotation: &[u8]) {
        let row = self.row;
        let annotation_row = u32::try_from(row).expect("annotation row exceeds u32");
        let byte = self.byte + annotation.len();
        self.out.entries[row] = IntervalEntry {
            start,
            end,
            annotation: annotation_row,
        };
        self.out.annotations.data[self.byte..byte].copy_from_slice(annotation);
        self.out.annotations.offsets[row + 1] = byte as i64;
        self.row += 1;
        self.byte = byte;
    }

    /// Entries written so far.
    pub fn rows(&self) -> usize {
        self.row
    }
}
