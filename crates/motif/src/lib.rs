//! Motif matching over [FFM](omicsweave_ffm) models.
//!
//! A collection of motifs is compiled once under an explicit [`Background`]
//! and [`Scoring`] policy into a [`Program`]; a scan then encodes each
//! stretch of sequence once and updates every motif's [`Reducer`]s from it.
//! Compilation runs in steps whose products are public (see
//! `docs/motif.md`):
//!
//! 1. estimate — normalized, smoothed factor probabilities;
//! 2. orient — the reverse strand's factors;
//! 3. contrast — additive log-odds terms against the background
//!    ([`ScoreModel`]);
//! 4. quantize — integer terms ([`IntModel`]);
//! 5. null — the exact integer score distribution ([`NullDist`]) and
//!    p-value/score conversions ([`Variant`]);
//! 6. lower — lookup blocks over packed k-mers;
//! 7. link — the reducers' needs ([`Need`]) become the kernel's mode.
//!
//! Sequences are given in memory as ASCII bytes ([`Strings`] rows). The
//! framework reports exact single-window statistics under the declared
//! background; multiple-testing corrections and enrichment statistics are
//! left to the caller.
//!
//! ```
//! use omicsweave_ffm::{Alphabet, FfmReader, FfmWriter, Motif, Schema, WriterOptions};
//! use omicsweave_motif::{Background, Hits, Program, ScanOptions, Scoring, Strands};
//!
//! let rows = vec![
//!     vec![0.9, 0.05, 0.03, 0.02],
//!     vec![0.02, 0.9, 0.05, 0.03],
//!     vec![0.03, 0.02, 0.9, 0.05],
//!     vec![0.05, 0.03, 0.02, 0.9],
//!     vec![0.9, 0.05, 0.03, 0.02],
//!     vec![0.02, 0.9, 0.05, 0.03],
//! ];
//! let mut w = FfmWriter::new(Vec::new(), Alphabet::dna(), Schema::default(), WriterOptions::default())?;
//! w.push("M1", &Motif::pfm(&rows, Some(20.0)), &[])?;
//! let reader = FfmReader::open(w.finish()?)?;
//! let motif = reader.motif(0)?;
//!
//! let alphabet = Alphabet::dna();
//! let program = Program::compile(
//!     &[motif.view()],
//!     &alphabet,
//!     &Background::uniform(&alphabet),
//!     &Scoring::fimo(),
//!     Strands::Both,
//!     Some(1),
//! )?;
//! let hits = program.scan_sequence(b"TTTTACGTACTTTT", Hits::pvalue(1e-3), &ScanOptions::default())?;
//! assert_eq!(hits.len(), 1);
//! assert_eq!(hits[0].start, 4);
//! # Ok::<(), omicsweave_core::Error>(())
//! ```

mod background;
mod factor;
pub mod fimo;
mod lower;
pub mod mast;
pub mod mcast;
mod model;
mod null;
pub mod profile;
mod program;
mod reduce;
mod scan;
mod score;
mod scoring;
pub mod stats;
pub mod substitution;

pub use background::Background;
pub use model::{MotifModel, StrandModel, Threshold, Variant};
pub use null::NullDist;
pub use program::{Part, Program};
pub use reduce::{
    BestEntry, CountMatrix, Demand, Dense, DenseScores, HitEntry, Hits, Match, Need, Reducer,
    ScanInfo, SequenceBest, SequenceChunk, SequenceCounts, Window, Windows,
};
pub use scan::ScanOptions;
pub use score::{FORBIDDEN, IntModel, IntTerm, ScoreModel, Term};
pub use scoring::{
    Ambiguity, AmbiguityRule, BackgroundRule, ColumnPrior, LogScale, Quantizer, Rounding, Scoring,
    Strands,
};

pub use omicsweave_core::{Strand, Strings};
