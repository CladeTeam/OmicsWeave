//! Accumulators: the state machines one pass over the BAM updates.
//!
//! Each [`Accumulator`] keeps a **batch** state, updated by worker threads
//! for a run of consecutive records ([`Accumulator::observe`], then
//! [`Accumulator::seal`]), and a **total**, into which the batches are
//! merged on one thread in record order ([`Accumulator::merge`]). The
//! total may spill (sort buffers, spools) within the share of the memory
//! budget its [`RunCtx`] grants; [`Accumulator::finish`] turns it into a
//! format-free product (see [`crate::model`]), which
//! [`Accumulator::deliver`] hands to the sinks. Because batches are merged
//! in order and every fold is associative, results do not depend on the
//! number of threads or on where batches were cut.
//!
//! Adding a product means writing one accumulator: say what it
//! [`Needs`], what memory it [`estimate`](Accumulator::estimate)s, how a
//! batch observes a record, how batches merge, and what it delivers.

pub mod alleles;
pub mod atac;
pub mod counts;
pub mod coverage;
pub mod hic;
pub mod labeling;
pub mod methylation;
pub mod splice;
pub mod velocity;

use std::any::Any;
use std::path::{Path, PathBuf};

use omicsweave_convert::stage::Outputs;
use omicsweave_core::Result;
use omicsweave_core::exec::Pool;

use crate::prelude::ReadView;
use crate::sink::json::Json;
use crate::source::Input;

/// What an accumulator reads from each record, so the pipeline computes
/// only what some accumulator uses. Combined with `|`.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub struct Needs(pub u32);

impl Needs {
    pub const NONE: Needs = Needs(0);
    /// CIGAR-derived blocks, junctions, ends.
    pub const CIGAR: Needs = Needs(1);
    /// Optional fields (`CB`, `UB`, `RG`, `NH`, `MD`, ...).
    pub const TAGS: Needs = Needs(1 << 1);
    /// The mate of each pair, paired across batches.
    pub const MATES: Needs = Needs(1 << 2);
    /// Feature assignment through a registered space.
    pub const ASSIGN: Needs = Needs(1 << 3);
    /// Read bases.
    pub const SEQUENCE: Needs = Needs(1 << 4);
    /// Base qualities.
    pub const QUALITY: Needs = Needs(1 << 5);
    /// Reference bases under each alignment (from a FASTA or 2bit).
    pub const REFERENCE: Needs = Needs(1 << 6);
    /// Known variant sites to mask (from a BCF).
    pub const VARIANTS: Needs = Needs(1 << 7);

    pub fn contains(self, other: Needs) -> bool {
        self.0 & other.0 == other.0
    }

    /// Names of the needs set, for reports.
    pub fn names(self) -> Vec<&'static str> {
        const NAMES: [&str; 8] = [
            "cigar",
            "tags",
            "mates",
            "assign",
            "sequence",
            "quality",
            "reference",
            "variants",
        ];
        (0..8)
            .filter(|i| self.0 >> i & 1 == 1)
            .map(|i| NAMES[i])
            .collect()
    }
}

impl std::ops::BitOr for Needs {
    type Output = Needs;
    fn bitor(self, o: Needs) -> Needs {
        Needs(self.0 | o.0)
    }
}

impl std::ops::BitOrAssign for Needs {
    fn bitor_assign(&mut self, o: Needs) {
        self.0 |= o.0;
    }
}

/// What a plan knows before reading records: the inputs' headers and the
/// sizes the pipeline uses.
#[derive(Debug, Clone)]
pub struct PlanInfo {
    pub threads: usize,
    /// Each input: references and compressed size.
    pub inputs: Vec<(Vec<(String, u32)>, u64)>,
    /// Record bytes per batch.
    pub batch_bytes: u64,
    /// Batches in flight.
    pub batches: usize,
}

impl PlanInfo {
    /// Compressed bytes of every input.
    pub fn input_bytes(&self) -> u64 {
        self.inputs.iter().map(|i| i.1).sum()
    }

    /// The longest reference of any input.
    pub fn longest_reference(&self) -> u32 {
        self.inputs
            .iter()
            .flat_map(|i| i.0.iter().map(|r| r.1))
            .max()
            .unwrap_or(0)
    }
}

/// The memory an accumulator needs.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub struct MemoryNeed {
    /// Held for the whole run and not spillable (tables, dictionaries,
    /// pending mates expected at the budget's scale).
    pub fixed: u64,
    /// Per batch in flight.
    pub per_batch: u64,
    /// The least spillable buffer it can work with.
    pub min_buffer: u64,
    /// Relative share of the spillable buffer it should get (0: none).
    pub weight: u32,
}

/// What a merging total works with: the worker pool, its share of the
/// buffer memory, where to spill, and the inputs.
pub struct RunCtx<'p> {
    pub pool: &'p Pool,
    /// Bytes of spillable buffer granted to this accumulator.
    pub buffer: u64,
    /// The directory spill directories are created in.
    pub tmp: PathBuf,
    pub inputs: &'p [Input],
}

/// Where products go: the output directory, the files being written, and
/// how matrices are stored.
pub struct Delivery<'a> {
    pub outputs: &'a mut Outputs<&'static str>,
    pub dir: &'a Path,
    pub pool: &'a Pool,
    pub matrix: crate::task::MatrixFormat,
    /// Where sinks put temporary files.
    pub tmp: &'a Path,
    /// Buffer bytes sinks may use (dense tables, compression in flight).
    pub buffer: u64,
}

impl Delivery<'_> {
    /// `dir/rel`, creating its parent directories.
    pub fn path(&self, rel: &str) -> Result<PathBuf> {
        let p = self.dir.join(rel);
        if let Some(parent) = p.parent() {
            std::fs::create_dir_all(parent)?;
        }
        Ok(p)
    }
}

/// What the batch states of one batch know.
#[derive(Debug, Clone, Copy)]
pub struct BatchInfo {
    pub input: u32,
}

/// A product's state machine; see the module documentation.
pub trait Accumulator: Send + Sync + 'static {
    /// State of one batch, built on a worker.
    type Batch: Send + 'static;
    /// Merged state, on the merging thread for the whole run (delivered on
    /// a thread of its own).
    type Total<'p>: Send;
    /// What [`Accumulator::finish`] makes of the total.
    type Product;

    /// The product's name: its output directory or file stem.
    fn name(&self) -> &str;
    fn needs(&self) -> Needs;
    fn estimate(&self, info: &PlanInfo) -> MemoryNeed;
    fn batch(&self, info: BatchInfo) -> Self::Batch;
    fn observe(&self, b: &mut Self::Batch, r: &ReadView<'_>) -> Result<()>;
    /// Ends a batch, still on its worker (sort, collapse, pair locally).
    fn seal(&self, _b: &mut Self::Batch) -> Result<()> {
        Ok(())
    }
    fn total<'p>(&'p self, rc: &RunCtx<'p>) -> Result<Self::Total<'p>>;
    /// Merges the next batch, in record order. `frontier` is the coordinate
    /// key ([`crate::prelude::coord_key`]) of the batch's last record: no
    /// later record of the input starts before it.
    fn merge<'p>(
        &'p self,
        t: &mut Self::Total<'p>,
        b: Self::Batch,
        frontier: u64,
        rc: &RunCtx<'p>,
    ) -> Result<()>;
    /// An input has ended (every batch of it is merged).
    fn end_input<'p>(
        &'p self,
        _t: &mut Self::Total<'p>,
        _input: u32,
        _rc: &RunCtx<'p>,
    ) -> Result<()> {
        Ok(())
    }
    fn finish<'p>(&'p self, t: Self::Total<'p>, rc: &RunCtx<'p>) -> Result<Self::Product>;
    /// Writes the product through the sinks; returns its summary.
    fn deliver(&self, p: Self::Product, d: &mut Delivery<'_>) -> Result<Json>;
}

/// An accumulator with its types erased, so a task holds any mix.
pub(crate) trait Stage: Send + Sync {
    fn name(&self) -> &str;
    fn needs(&self) -> Needs;
    fn estimate(&self, info: &PlanInfo) -> MemoryNeed;
    fn batch(&self, info: BatchInfo) -> Box<dyn Any + Send>;
    fn observe(&self, b: &mut (dyn Any + Send), r: &ReadView<'_>) -> Result<()>;
    fn seal(&self, b: &mut (dyn Any + Send)) -> Result<()>;
    fn start<'p>(&'p self, rc: RunCtx<'p>) -> Result<Box<dyn Merging + 'p>>;
}

/// A running total of an erased accumulator.
pub(crate) trait Merging: Send {
    fn merge(&mut self, b: Box<dyn Any + Send>, frontier: u64) -> Result<()>;
    fn end_input(&mut self, input: u32) -> Result<()>;
    fn finish(self: Box<Self>, d: &mut Delivery<'_>) -> Result<Json>;
}

pub(crate) struct Erased<A>(pub A);

struct Running<'p, A: Accumulator> {
    acc: &'p A,
    total: A::Total<'p>,
    rc: RunCtx<'p>,
}

impl<A: Accumulator> Stage for Erased<A> {
    fn name(&self) -> &str {
        self.0.name()
    }

    fn needs(&self) -> Needs {
        self.0.needs()
    }

    fn estimate(&self, info: &PlanInfo) -> MemoryNeed {
        self.0.estimate(info)
    }

    fn batch(&self, info: BatchInfo) -> Box<dyn Any + Send> {
        Box::new(self.0.batch(info))
    }

    #[inline]
    fn observe(&self, b: &mut (dyn Any + Send), r: &ReadView<'_>) -> Result<()> {
        let b = b.downcast_mut::<A::Batch>().expect("batch of this stage");
        self.0.observe(b, r)
    }

    fn seal(&self, b: &mut (dyn Any + Send)) -> Result<()> {
        self.0
            .seal(b.downcast_mut::<A::Batch>().expect("batch of this stage"))
    }

    fn start<'p>(&'p self, rc: RunCtx<'p>) -> Result<Box<dyn Merging + 'p>> {
        let total = self.0.total(&rc)?;
        Ok(Box::new(Running {
            acc: &self.0,
            total,
            rc,
        }))
    }
}

impl<'p, A: Accumulator> Merging for Running<'p, A> {
    fn merge(&mut self, b: Box<dyn Any + Send>, frontier: u64) -> Result<()> {
        let b = *b.downcast::<A::Batch>().expect("batch of this stage");
        self.acc.merge(&mut self.total, b, frontier, &self.rc)
    }

    fn end_input(&mut self, input: u32) -> Result<()> {
        self.acc.end_input(&mut self.total, input, &self.rc)
    }

    fn finish(self: Box<Self>, d: &mut Delivery<'_>) -> Result<Json> {
        let Running { acc, total, rc } = *self;
        let product = acc.finish(total, &rc)?;
        acc.deliver(product, d)
    }
}
