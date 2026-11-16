//! Splicing products against direct computation on generated reads: the
//! junction table (STAR `SJ.out.tab` columns, motifs from a bgzip FASTA
//! and a 2bit genome), the junction × barcode matrix, and the velocity
//! layers of gene matrices; outputs independent of threads and spilling.

mod common;
mod support;

use std::collections::{BTreeMap, BTreeSet, HashMap};
use std::path::{Path, PathBuf};
use std::sync::Arc;

use common::{
    Gene, Read, Rng, Spec, genes, genome, gtf, gunzip, read_10x, read_10x_layer, reverse_for,
    synth, write_bam, write_fasta, write_twobit,
};
use omicsweave_gtf::{Annotation, ReadOptions};
use omicsweave_quant::{
    Budget, Counts, JunctionMatrix, Junctions, MatrixFormat, Product, Task, quantify,
};

#[path = "cases/splicing/fixtures.rs"]
mod fixtures;
use fixtures::*;
#[path = "cases/splicing/oracle.rs"]
mod oracle;
use oracle::*;
#[path = "cases/splicing/execution.rs"]
mod execution;
#[path = "cases/splicing/junctions.rs"]
mod junctions;
#[path = "cases/splicing/velocity.rs"]
mod velocity;
