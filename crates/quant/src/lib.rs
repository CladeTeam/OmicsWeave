//! Quantification of aligned reads: one pass over coordinate-sorted BAMs
//! into coverage tracks, count matrices (with velocity and labeling
//! layers), single-cell ATAC fragments, splice junctions, metabolic
//! labeling statistics, DNA methylation, Hi-C contacts and allele counts,
//! under a budget of threads and memory.
//!
//! The input is an already aligned BAM; nothing is realigned. Unlike the
//! format crates and `omicsweave-convert`, this layer analyses: it filters,
//! assigns reads to features and counts.
//!
//! # Pipeline
//!
//! - [`source`]: each input's header is checked (`@HD SO:coordinate`);
//!   BGZF blocks inflate on the worker pool and the records are cut into
//!   batches of consecutive records, their coordinate order checked on the
//!   way (an unsorted BAM is an error that names the fix).
//! - [`prelude`]: per record, once and only on demand — FLAG / MAPQ / NH
//!   filters, optional fields, library strandedness, aligned blocks,
//!   junctions, fragment spans, Tn5 cut sites, aligned segments, the
//!   reference bases under the read and the read's bases against them
//!   ([`ReadView`]) — and mate pairing within a batch and across batches.
//! - [`reference`]: the task's reference genome (bgzip FASTA or 2bit) and
//!   known variants (BCF), read by each batch through a window over its
//!   records; SNV panels for allele counts ([`SitePanel`]).
//! - [`assign`]: feature spaces over [`omicsweave_gtf::IntervalSet`]s
//!   (genes, gene windows, bigBed peaks, genome bins) and per-batch sweep
//!   cursors; a record is assigned to a space once and the result shared.
//! - [`accum`]: the [`Accumulator`]s, one per product: batches are
//!   observed on workers and merged in record order on one thread, so the
//!   output does not depend on the thread count; totals spill within
//!   their share of the budget.
//! - [`model`]: products as data ([`Track`], [`CountMatrix`],
//!   [`FragmentTable`], [`JunctionTable`], [`Histogram`], ...), apart from
//!   the [`sink`]s that write them (bigWig, 10x directories with CellPack /
//!   Matrix Market / CSV matrices, `fragments.tsv.gz` + `.tbi`,
//!   `junctions.tab`, `summary.json`).
//!
//! [`plan`] checks a [`Task`] against its [`Budget`] before any record is
//! read ([`MemoryPlan`]: fixed tables, the batches in flight, the buffers
//! the accumulators share); [`run`] makes the one pass and writes every
//! product, each file under a temporary name renamed into place at the
//! end.
//!
//! # Output layout
//!
//! Under the output directory: `<coverage>.bw` (or `.forward.bw` /
//! `.reverse.bw`); `<counts>/` (`barcodes.tsv.gz`, `features.tsv.gz`,
//! `matrix.cpk` | `matrix.mtx.gz` | `matrix.csv.gz`, further layers
//! `spliced.cpk`, `unspliced.cpk`, `ambiguous.cpk`, `labeled.cpk`, and
//! `counts.tsv` for bulk tables); `atac/fragments.tsv.gz` (+ `.tbi`),
//! `atac/barcode_metrics.csv`, `atac/peaks/`, `atac/tiles/`,
//! `atac/gene_activity/`; `junctions/junctions.tab` (with a matrix, a 10x
//! directory beside it); `labeling/conversions.tsv.gz`,
//! `labeling/mismatches.tsv`; `methylation/cytosines.tsv.gz` (+ `.tbi`),
//! `methylation/beta.bw`, `methylation/coverage.bw`,
//! `methylation/mbias.tsv`, `methylation/<matrix>/` (`methylated.cpk`,
//! `total.cpk`); `hic/contacts.pairs.gz`, `hic/contacts.hic`,
//! `hic/cell_metrics.csv`; `alleles/` (`ref.cpk`, `alt.cpk`, `other.cpk`),
//! `alleles/variants.tsv.gz`, `alleles/genes/` (`hap1.cpk`, `hap2.cpk`,
//! `both.cpk`); and `summary.json`.

pub mod accum;
pub mod assign;
pub mod model;
pub mod prelude;
pub mod reference;
pub mod sink;
pub mod source;
pub mod spool;
pub mod task;

pub use accum::alleles::{AlleleGenes, Alleles};
pub use accum::atac::Atac;
pub use accum::counts::{CountUnit, Counts, Dedup, FeatureSource, Labeled, RowKey};
pub use accum::coverage::{Coverage, CoverageUnit, Scale, StrandSplit};
pub use accum::hic::Contacts;
pub use accum::labeling::{Conversion, Labeling};
pub use accum::methylation::{
    CHG, CHH, CPG, CallSource, MethylMatrix, Methylation, contexts_from_names,
};
pub use accum::splice::{JunctionMatrix, Junctions};
pub use accum::{Accumulator, MemoryNeed, Needs, PlanInfo};
pub use assign::{FeatureSpace, Query, SpaceId};
pub use model::{
    BarcodeTable, Contact, ContactTable, CountMatrix, Cytosine, FeatureList, FeatureTable,
    FragmentTable, Histogram, Junction, JunctionTable, Mbias, MismatchTable, Site, Track,
    VariantTable,
};
pub use omicsweave_convert::{Budget, MemoryPlan, OnProgress, Progress, human, parse_size};
pub use prelude::{AlignedBase, Filter, ReadView, Strandedness};
pub use reference::{Reference, SitePanel, Variants};
pub use sink::json::Json;
pub use source::{Input, SORT_HINT};
pub use task::{MatrixFormat, Plan, Product, Report, Task, plan, quantify, run};
