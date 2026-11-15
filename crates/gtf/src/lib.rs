//! Gene annotation (GTF 2.2 and GFF3) as a compact, columnar, indexed
//! model: the reusable annotation layer for quantification, feature
//! counting and annotation-aware sampling.
//!
//! # Reading
//!
//! [`Annotation::read`] reads GENCODE, Ensembl, RefSeq and UCSC GTF and
//! GFF3 files, plain, gzip or BGZF (told by content). The text is cut into
//! line batches parsed in parallel on an [`omicsweave_core::exec`] pool and
//! merged in order, so the result is the same for any thread count.
//! [`Annotation::parse`] does the same for bytes in memory.
//!
//! ```no_run
//! use omicsweave_gtf::{Annotation, ReadOptions};
//! # fn main() -> omicsweave_core::Result<()> {
//! let ann = Annotation::read("gencode.v44.annotation.gtf.gz", &ReadOptions::default())?;
//! println!("{} genes, {} problems", ann.genes().len(), ann.report().problems.len());
//! # Ok(())
//! # }
//! ```
//!
//! Normalization:
//! - Coordinates are 0-based, half-open `u32` (`start - 1`, `end` of the
//!   file); strand is a [`Strand`], phase `0..=2` or [`PHASE_NONE`]; the
//!   score is kept as written ([`NONE`] for `.`).
//! - The hierarchy is gene → transcript → feature. GTF groups records by
//!   (seqid, `gene_id`) and (seqid, `transcript_id`); GFF3 follows `ID` /
//!   `Parent`, recognizing transcripts by their exon-level children
//!   (any type name: `mRNA`, `lnc_RNA`, `V_gene_segment`, ...). Genes and
//!   transcripts without a record of their own are synthesized from their
//!   children (span = min / max of the children).
//! - Gene and transcript ids, names, biotypes and transcript tags are read
//!   from the usual keys of each source (`gene_type` / `gene_biotype` /
//!   `biotype`, `gene_name` / `gene` / `Name`, `tag`); known tags are also
//!   bits ([`Tags`]).
//! - [`ReadOptions::aliases`] renames seqids at read;
//!   [`Annotation::match_seqids`] maps another file's names (a BAM header)
//!   to seqids, with the built-in `chr1` ↔ `1`, `chrM` ↔ `MT` rules.
//! - Problems (orphans, duplicate ids, strand conflicts, features outside
//!   their parent, ...) are collected in [`Annotation::report`];
//!   [`ReadOptions::strict`] makes the first one an error. Unsorted input
//!   is fine.
//!
//! # Memory model
//!
//! No per-record heap objects: [`Genes`], [`Transcripts`], [`Features`]
//! and [`Others`] are struct-of-arrays; every string (seqids, sources,
//! types, ids, names, biotypes, attribute keys and values) is interned once
//! in a [`StringPool`] (one byte arena + offsets) and referenced by `u32`.
//! A gene's transcripts are one contiguous row range
//! ([`Genes::transcripts`]), a transcript's features another
//! ([`Transcripts::features`]), sorted by position; genes are sorted by
//! (seqid, start, end), so gene rows are also genomic order. Attributes
//! are one flat (key, value) table by file record, kept according to
//! [`AttributePolicy`]: `All` (lossless, original order, needed to write
//! the GTF back), `Only(keys)` or `None`.
//!
//! # Lookups and queries
//!
//! - By id or name: [`Annotation::gene`], [`Annotation::genes_by_id`],
//!   [`Annotation::genes_by_name`], [`Annotation::transcript`], ...
//! - By position: [`Annotation::interval_index`]`(`[`Level`]`)` is an
//!   [`IntervalIndex`] over genes, transcripts or exons (ids are table
//!   rows): `overlaps`, `contained_in`, `containing`, `nearest`,
//!   `overlaps_batch`.
//! - Along a chromosome: [`SweepCursor`] answers queries with
//!   non-decreasing starts (a sorted BAM) in amortized O(1) per step,
//!   writing ids into a caller buffer:
//!
//! ```no_run
//! # use omicsweave_gtf::{Annotation, Level};
//! # fn f(ann: &Annotation, seqid: u32, reads: &[(u32, Vec<(u32, u32)>)]) {
//! let mut exons = ann.cursor(Level::Exon, seqid);
//! let mut hits = Vec::new();
//! for (read_start, blocks) in reads {
//!     exons.advance(*read_start); // read starts never decrease
//!     hits.clear();
//!     for &(s, e) in blocks {
//!         exons.overlaps(s, e, &mut hits); // feature rows of exons
//!     }
//! }
//! # }
//! ```
//!
//! # Derived views
//!
//! Computed on demand from the tables, the parameterless ones cached:
//! [`Annotation::exon_blocks`], [`Annotation::introns`],
//! [`Annotation::exonic_regions`] (union exons per gene),
//! [`Annotation::junctions`] (donor / acceptor, strand, genes and
//! transcripts per junction), [`Annotation::exon_bins`] (DEXSeq-like
//! disjoint bins), [`Annotation::transcript_lengths`],
//! [`Annotation::cds_lengths`]. [`Annotation::derive`] builds them for a
//! [`Selection`], clipped to chromosome lengths, and adds count-feature
//! windows: TSS / TES, promoters and TES windows with configurable
//! upstream / downstream ([`Flank`]), and gene bodies extended upstream
//! (ArchR / Signac gene activity).
//!
//! Interval views are [`IntervalSet`]s: (seqid, start, end, strand, owner,
//! name) columns with their own index and sweep cursor. External features
//! (peaks, fixed bins) use the same type through [`IntervalSetBuilder`].
//!
//! # Subsets
//!
//! [`Annotation::select`] returns a [`Selection`] (gene and transcript row
//! lists, no copy) by biotype, tags, gene ids / names or seqids.
//!
//! # Writing
//!
//! [`Annotation::write_gtf`] writes the original records back (lossless
//! for a GTF read with all attributes) or the model in canonical GTF
//! ([`GtfLayout`]); [`Annotation::write_bed12`] writes transcripts.

mod alias;
mod annotation;
mod build;
mod derived;
mod index;
mod input;
mod model;
mod parse;
mod read;
mod report;
mod select;
mod set;
mod strings;
mod write;

pub use alias::{Aliases, builtin_alternatives, match_names};
pub use annotation::Annotation;
pub use derived::{Derive, ExonBins, Flank, Junctions, TranscriptPart};
pub use index::{IntervalIndex, SweepCursor};
pub use input::{Compression, detect_compression};
pub use model::{
    AttributePolicy, Attributes, Dialect, FeatureKind, Features, Genes, Others, PHASE_NONE, Row,
    Strand, Tags, Transcripts,
};
pub use read::ReadOptions;
pub use report::{Problem, ProblemKind, Report};
pub use select::{Select, Selection};
pub use set::{Csr, IntervalSet, IntervalSetBuilder, Level};
pub use strings::{IdIndex, NONE, StringPool};
pub use write::GtfLayout;
