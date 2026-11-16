//! Tasks as a whole: one pass for every product, output independent of the
//! thread count, unsorted inputs and small budgets refused before any
//! output appears.

mod common;

use std::path::Path;
use std::sync::Arc;

use common::{Rng, Spec, files, genes, gtf, synth, write_bam};
use omicsweave_gtf::{Annotation, ReadOptions};
use omicsweave_quant::{
    Atac, Budget, Counts, Coverage, CoverageUnit, Product, StrandSplit, Strandedness, Task, plan,
    quantify,
};

#[path = "cases/task/boundaries.rs"]
mod boundaries;
#[path = "cases/task/execution.rs"]
mod execution;
#[path = "cases/task/validation.rs"]
mod validation;
