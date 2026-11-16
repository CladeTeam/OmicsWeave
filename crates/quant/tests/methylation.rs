//! DNA methylation against direct computation on bisulfite read pairs drawn
//! from a generated genome (directional library, both original strands,
//! molecules methylated per cytosine at random, sequencing errors, `N`s,
//! indels and clips, overlapping mates, filtered records, known SNPs) and on
//! long reads with `MM` / `ML` tags: calls from the reference, from `MD`
//! tags and from Bismark `XM` tags; the per-cytosine table (with the tabix
//! index htslib builds for it), the beta and coverage tracks, M-bias,
//! region × barcode and × file matrices; the same files at any thread count
//! and with spilled buffers.

mod common;
mod support;

use std::collections::{BTreeMap, BTreeSet, HashMap};
use std::path::{Path, PathBuf};

use common::{
    Read, Rng, assert_tabix_as_htslib, bigwig_values, complement, genome, gunzip, md_tag,
    read_10x_layer, write_bam, write_bcf, write_fasta,
};
use omicsweave_quant::{
    Budget, CHG, CHH, CPG, CallSource, FeatureSource, MethylMatrix, Methylation, Product, RowKey,
    Task, plan, quantify,
};

#[path = "cases/methylation/fixtures.rs"]
mod fixtures;
use fixtures::*;
#[path = "cases/methylation/oracle.rs"]
mod oracle;
use oracle::*;
#[path = "cases/methylation/output.rs"]
mod output;
use output::*;
#[path = "cases/methylation/execution.rs"]
mod execution;
#[path = "cases/methylation/matrices.rs"]
mod matrices;
#[path = "cases/methylation/modification_tags.rs"]
mod modification_tags;
#[path = "cases/methylation/reference_calls.rs"]
mod reference_calls;
#[path = "cases/methylation/tag_calls.rs"]
mod tag_calls;
