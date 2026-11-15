//! Oriented keys of dense queries: one or more spans of a chromosome read
//! as one sequence, 5' to 3' on a strand.

use crate::error::{Error, Result};
use crate::region::Region;

/// Strand of a key or a record, as the strand column of BED, GTF and GFF3
/// spells it.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, PartialOrd, Ord, Default)]
#[repr(u8)]
pub enum Strand {
    /// `+`
    Forward = 0,
    /// `-`
    Reverse = 1,
    /// `.`: not stranded.
    #[default]
    Unstranded = 2,
    /// `?`: stranded, strand unknown.
    Unknown = 3,
}

impl Strand {
    pub fn from_byte(b: u8) -> Option<Strand> {
        Some(match b {
            b'+' => Strand::Forward,
            b'-' => Strand::Reverse,
            b'.' => Strand::Unstranded,
            b'?' => Strand::Unknown,
            _ => return None,
        })
    }

    pub fn as_char(self) -> char {
        match self {
            Strand::Forward => '+',
            Strand::Reverse => '-',
            Strand::Unstranded => '.',
            Strand::Unknown => '?',
        }
    }

    /// `1`, `-1`, or `0` when not `+` / `-`.
    pub fn sign(self) -> i8 {
        match self {
            Strand::Forward => 1,
            Strand::Reverse => -1,
            _ => 0,
        }
    }

    /// Whether the 5' end is the high coordinate. Only `-` is; keys and
    /// records without a known strand read as `+`.
    #[inline]
    pub fn is_reverse(self) -> bool {
        self == Strand::Reverse
    }
}

/// Bases `[start, end)` of a chromosome, 0-based, half-open.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, PartialOrd, Ord, Default)]
#[repr(C)]
pub struct Span {
    pub start: u32,
    pub end: u32,
}

impl Span {
    pub fn new(start: u32, end: u32) -> Span {
        Span { start, end }
    }

    pub fn len(&self) -> u32 {
        self.end.saturating_sub(self.start)
    }

    pub fn is_empty(&self) -> bool {
        self.len() == 0
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum Spans<'a> {
    One(Span),
    Many(&'a [Span]),
}

/// The key of every dense query (bases, per-base signal, signal bins): the
/// spans of one chromosome, in ascending order and disjoint, concatenated
/// into one sequence of [`Locus::len`] positions and read 5' to 3' of
/// `strand`. A [`Region`] is the one-span case; a spliced transcript is a
/// locus of its exons.
///
/// Position `i` of the output is base `i` of the forward concatenation, or
/// base `len - 1 - i` of it when the strand is [`Strand::Reverse`]; bases
/// of a reverse locus are also complemented.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Locus<'a> {
    pub chrom: &'a str,
    spans: Spans<'a>,
    pub strand: Strand,
    len: u32,
}

impl<'a> Locus<'a> {
    /// A locus of `spans`, which must be ascending and disjoint (they may
    /// touch). Empty spans are allowed and contribute nothing.
    pub fn new(chrom: &'a str, spans: &'a [Span], strand: Strand) -> Result<Locus<'a>> {
        let mut len: u32 = 0;
        let mut floor = 0;
        for (i, span) in spans.iter().enumerate() {
            if span.start > span.end {
                return Err(Error::InvalidRegion {
                    start: span.start,
                    end: span.end,
                });
            }
            if i > 0 && span.start < floor {
                return Err(Error::invalid(
                    "locus",
                    format!(
                        "span {}..{} starts before the end {floor} of the span before it",
                        span.start, span.end
                    ),
                ));
            }
            floor = span.end;
            len += span.len();
        }
        let spans = match spans {
            [one] => Spans::One(*one),
            _ => Spans::Many(spans),
        };
        Ok(Locus {
            chrom,
            spans,
            strand,
            len,
        })
    }

    /// The one-span locus of `region` on `strand`.
    pub fn region(region: Region<'a>, strand: Strand) -> Result<Locus<'a>> {
        region.check()?;
        Ok(Locus {
            chrom: region.chrom,
            spans: Spans::One(Span::new(region.start, region.end)),
            strand,
            len: region.len(),
        })
    }

    pub fn spans(&self) -> &[Span] {
        match &self.spans {
            Spans::One(span) => std::slice::from_ref(span),
            Spans::Many(spans) => spans,
        }
    }

    /// Positions of the output: the sum of the span lengths.
    pub fn len(&self) -> u32 {
        self.len
    }

    pub fn is_empty(&self) -> bool {
        self.len == 0
    }

    pub fn is_reverse(&self) -> bool {
        self.strand.is_reverse()
    }

    /// The smallest region holding every span; empty at the first span's
    /// start (or 0) when the locus is empty.
    pub fn hull(&self) -> Region<'a> {
        let spans = self.spans();
        match (spans.first(), spans.last()) {
            (Some(first), Some(last)) => Region::new(self.chrom, first.start, last.end),
            _ => Region::new(self.chrom, 0, 0),
        }
    }

    /// The non-empty spans as regions, each with its offset in the forward
    /// concatenation.
    pub fn pieces(&self) -> impl Iterator<Item = (u32, Region<'a>)> + '_ {
        let chrom = self.chrom;
        let mut at = 0;
        self.spans().iter().filter(|s| !s.is_empty()).map(move |s| {
            let piece = (at, Region::new(chrom, s.start, s.end));
            at += s.len();
            piece
        })
    }

    /// Fails unless every span lies within a sequence of `len` bases.
    pub fn check_within(&self, len: u32) -> Result<()> {
        self.check()?;
        if let Some(last) = self.spans().last()
            && last.end > len
        {
            return Err(Error::RegionOutOfRange {
                start: last.start,
                end: last.end,
                len,
            });
        }
        Ok(())
    }
}

impl<'a> From<Region<'a>> for Locus<'a> {
    /// The unstranded one-span locus. A reversed region keeps its bounds and
    /// fails when a reader checks it.
    fn from(region: Region<'a>) -> Locus<'a> {
        Locus {
            chrom: region.chrom,
            spans: Spans::One(Span::new(region.start, region.end)),
            strand: Strand::Unstranded,
            len: region.len(),
        }
    }
}

impl Locus<'_> {
    /// Fails if a span is reversed. Loci built with [`Locus::new`] and
    /// [`Locus::region`] always pass; one converted from a [`Region`] may
    /// not.
    pub fn check(&self) -> Result<()> {
        // Many is created only by `new`, which validates its immutable
        // borrowed spans once. Only the unchecked Region conversion can
        // introduce a reversed one-span key.
        if let Spans::One(span) = self.spans
            && span.start > span.end
        {
            return Err(Error::InvalidRegion {
                start: span.start,
                end: span.end,
            });
        }
        Ok(())
    }
}
