//! SAM text encoded as uncompressed BAM bytes, byte for byte as htslib
//! encodes it (`sam_hdr_read` / `sam_read1` from SAM, then `sam_hdr_write`
//! / `sam_write1` to BAM).
//!
//! The work splits the way a parallel SAM → BAM conversion does:
//!
//! 1. [`SamHeader::read`] consumes the header lines from decompressed SAM
//!    text; [`bam_header_bytes`] gives the BAM header.
//! 2. The caller cuts the remaining text into lines and hands batches to
//!    workers. [`BamEncoder`] is shared (`Send + Sync`); each worker owns a
//!    [`BamScratch`] and appends each line's complete record (with its
//!    `block_size`) to a buffer of its own with
//!    [`BamEncoder::encode_line`], which also returns the [`BamSpan`] a BAI
//!    or CSI index records.
//! 3. The header bytes followed by the records, in line order, are the
//!    uncompressed BAM stream; BGZF framing, compression and indexing
//!    happen elsewhere.
//!
//! ```
//! use omicsweave_sam::{BamEncoder, BamScratch, SamHeader, bam_header_bytes};
//!
//! let sam = b"@SQ\tSN:chr1\tLN:1000\nr1\t0\tchr1\t11\t60\t4M\t*\t0\t0\tACGT\tIIII\n";
//! let mut input = &sam[..];
//! let header = SamHeader::read(&mut input).unwrap();
//! let mut bam = bam_header_bytes(&header);
//! let encoder = BamEncoder::new(&header);
//! let mut scratch = BamScratch::new();
//! // Lines without their "\n"; the text after the last "\n" is not a line.
//! let body = input.strip_suffix(b"\n").unwrap_or(input);
//! for (i, line) in body.split(|&c| c == b'\n').enumerate() {
//!     let line_no = header.lines() + 1 + i as u64;
//!     let span = encoder.encode_line(&mut scratch, line, line_no, &mut bam).unwrap();
//!     assert_eq!((span.tid, span.beg, span.end, span.mapped), (0, 10, 14, true));
//! }
//! ```
//!
//! Deliberate departures from htslib, and the line-ending contract, are
//! listed on [`BamEncoder`] and [`SamHeader`].

mod cnum;
mod header;
mod record;

pub use header::{Reference, SamHeader, bam_header_bytes};
pub use record::{BamEncoder, BamScratch, BamSpan};
