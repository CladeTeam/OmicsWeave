//! A quantification task: inputs, products and options; checking it
//! against a budget ([`plan`]) and running it ([`run`]).

use std::any::Any;
use std::cell::RefCell;
use std::collections::HashMap;
use std::io::Read;
use std::path::{Path, PathBuf};
use std::sync::Arc;
use std::time::{Duration, Instant};

use omicsweave_bam::Record;
use omicsweave_convert::stage::{Outputs, Sizing, sizing, text_batch};
use omicsweave_convert::{Budget, MemoryPlan, OnProgress, Progress, human};
use omicsweave_core::exec::{Ordered, Pool};
use omicsweave_core::{Error, Result};
use omicsweave_gtf::{Aliases, Annotation, Flank, Selection};

use crate::accum::alleles::{AlleleStage, Alleles};
use crate::accum::atac::{Atac, FragmentStage};
use crate::accum::counts::{
    CountStage, CountUnit, Counts, Dedup, FeatureSource, Fold, Reduce, RowKey, Rows, Velocity,
};
use crate::accum::coverage::Coverage;
use crate::accum::hic::{ContactStage, Contacts};
use crate::accum::labeling::Labeling;
use crate::accum::methylation::{MatrixStage, MethylStage, Methylation};
use crate::accum::splice::{JunctionStage, Junctions};
use crate::accum::velocity::GeneModels;
use crate::accum::{
    BatchInfo, Delivery, Erased, MemoryNeed, Merging, Needs, PlanInfo, RunCtx, Stage,
};
use crate::assign::{BatchSpaces, FeatureSpace, SpaceId, Spaces};
use crate::prelude::ReadView;
use crate::reference::{BatchCtx, MaskWindow, RefWindow, Reference, SitePanel, Variants};
use crate::sink::json::Json;
use crate::source::{BatchReader, Input};

/// How count matrices are stored.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum MatrixFormat {
    /// CellPack `.cpk`.
    #[default]
    Cpk,
    /// Matrix Market `.mtx.gz`.
    Mtx,
    /// Dense `.csv.gz`.
    Csv,
}

impl MatrixFormat {
    pub fn from_name(name: &str) -> Result<MatrixFormat> {
        Ok(match name {
            "cpk" => MatrixFormat::Cpk,
            "mtx" => MatrixFormat::Mtx,
            "csv" => MatrixFormat::Csv,
            _ => {
                return Err(Error::invalid(
                    "matrix",
                    format!("unknown matrix format {name:?}: use cpk, mtx or csv"),
                ));
            }
        })
    }
}

/// A product of a task.
#[derive(Debug, Clone)]
pub enum Product {
    Coverage(Coverage),
    Counts(Counts),
    Atac(Atac),
    Junctions(Junctions),
    Labeling(Labeling),
    Methylation(Methylation),
    Contacts(Contacts),
    Alleles(Alleles),
}

/// Input BAMs, the products of one pass over them, and where they go.
#[derive(Debug, Clone)]
pub struct Task {
    /// Coordinate-sorted BAMs, read one after the other.
    pub inputs: Vec<PathBuf>,
    pub out_dir: PathBuf,
    pub annotation: Option<Arc<Annotation>>,
    /// Chromosome aliases for matching BAM references to features.
    pub aliases: Option<Aliases>,
    pub products: Vec<Product>,
    pub matrix: MatrixFormat,
    /// The genome the BAMs were aligned to (bgzip or plain FASTA with its
    /// `.fai` / `.gzi`, or 2bit): reference bases for labeling (else the
    /// `MD` tags) and junction motifs.
    pub reference: Option<PathBuf>,
    /// Known variant sites (BGZF BCF with its CSI), masked when counting
    /// conversions and methylation calls.
    pub variants: Option<PathBuf>,
}

impl Task {
    pub fn new(inputs: Vec<PathBuf>, out_dir: impl Into<PathBuf>) -> Task {
        Task {
            inputs,
            out_dir: out_dir.into(),
            annotation: None,
            aliases: None,
            products: Vec::new(),
            matrix: MatrixFormat::Cpk,
            reference: None,
            variants: None,
        }
    }
}

/// What the stages share with the batch jobs.
pub(crate) struct Shared {
    pub stages: Vec<Box<dyn Stage>>,
    pub spaces: Spaces,
    pub n_refs: Vec<usize>,
    /// The reference and variants, when some stage reads them per record.
    pub reference: Option<Arc<Reference>>,
    pub variants: Option<Arc<Variants>>,
}

/// A task checked against its budget, with its memory division.
pub struct Plan {
    pub task: Task,
    pub budget: Budget,
    pub memory: MemoryPlan,
    /// Everything the accumulators read from records.
    pub needs: Needs,
    /// Each accumulator's name and the buffer bytes granted to it.
    pub stages: Vec<(String, u64)>,
    pub inputs: Vec<Input>,
    pub(crate) shared: Arc<Shared>,
    batch_bytes: usize,
    batches: usize,
    inflate: Sizing,
}

impl std::fmt::Debug for Plan {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("Plan")
            .field("budget", &self.budget)
            .field("memory", &self.memory)
            .field("stages", &self.stages)
            .finish()
    }
}

/// What a finished run did.
#[derive(Debug, Clone)]
pub struct Report {
    /// Files written with their kind (`bigwig`, `cpk`, `fragments`, ...).
    pub files: Vec<(&'static str, PathBuf)>,
    /// Records read, all inputs.
    pub records: u64,
    pub bytes_read: u64,
    pub elapsed: Duration,
    /// Of which the pass over the inputs (the rest writes the products).
    pub scan: Duration,
    /// What `summary.json` holds.
    pub summary: Json,
}

/// What the scanner holds besides the blocks in flight.
const SCAN_TEXT: u64 = 8 << 20;
const READ_BUF: u64 = 1 << 20;

fn read_list(path: &Path) -> Result<Vec<Vec<u8>>> {
    let raw = std::fs::read(path).map_err(|e| Error::Io {
        kind: e.kind(),
        detail: format!("{}: {e}", path.display()),
    })?;
    let text = if raw.starts_with(&[0x1f, 0x8b]) {
        let mut out = Vec::new();
        flate2::read::MultiGzDecoder::new(&raw[..]).read_to_end(&mut out)?;
        out
    } else {
        raw
    };
    Ok(text
        .split(|&c| c == b'\n')
        .map(|l| l.strip_suffix(b"\r").unwrap_or(l))
        .filter(|l| !l.is_empty())
        .map(|l| l.split(|&c| c == b'\t').next().unwrap_or(l).to_vec())
        .collect())
}

fn index_of(names: &[Vec<u8>]) -> HashMap<Vec<u8>, u32> {
    let mut m = HashMap::with_capacity(names.len());
    for (i, n) in names.iter().enumerate() {
        m.entry(n.clone()).or_insert(i as u32);
    }
    m
}

/// Builds what a task's products need: spaces and resolved stages.
struct Resolver<'a> {
    task: &'a Task,
    inputs: &'a [Input],
    spaces: Spaces,
    lengths: Option<Vec<Option<u32>>>,
    reference: Option<Arc<Reference>>,
}

impl Resolver<'_> {
    fn annotation(&self, what: &str) -> Result<&Arc<Annotation>> {
        self.task.annotation.as_ref().ok_or_else(|| {
            Error::invalid("quant", format!("{what} needs a gene annotation (--gtf)"))
        })
    }

    fn lengths(&mut self) -> Result<Vec<Option<u32>>> {
        if self.lengths.is_none() {
            let ann = self.annotation("a gene window")?;
            let chroms: Vec<(&str, u32)> = self.inputs[0]
                .refs
                .iter()
                .map(|r| (r.0.as_str(), r.1))
                .collect();
            self.lengths = Some(ann.seqid_lengths(&chroms, self.task.aliases.as_ref()));
        }
        Ok(self.lengths.clone().expect("set above"))
    }

    fn space(&mut self, c: &Counts) -> Result<SpaceId> {
        self.space_of(&c.features, c.subset.as_ref())
    }

    fn space_of(
        &mut self,
        features: &FeatureSource,
        subset: Option<&Arc<Selection>>,
    ) -> Result<SpaceId> {
        if subset.is_some()
            && !matches!(
                features,
                FeatureSource::GeneExons
                    | FeatureSource::GeneSpans
                    | FeatureSource::Promoters(_)
                    | FeatureSource::GeneActivity(_)
            )
        {
            return Err(Error::invalid("quant", "a gene subset needs gene features"));
        }
        let selection = subset.map(Arc::as_ref);
        // One space per subset object: a shared subset builds its space once.
        let sel_key = format!("{:?}", subset.map(Arc::as_ptr));
        match features {
            FeatureSource::GeneExons => {
                let ann = Arc::clone(self.annotation("gene counts")?);
                self.spaces.register(format!("gene-exons {sel_key}"), || {
                    if let Some(s) = selection {
                        s.check(&ann)?;
                    }
                    Ok(FeatureSpace::gene_exons(&ann, selection))
                })
            }
            FeatureSource::GeneSpans => {
                let ann = Arc::clone(self.annotation("gene counts")?);
                self.spaces.register(format!("gene-spans {sel_key}"), || {
                    if let Some(s) = selection {
                        s.check(&ann)?;
                    }
                    Ok(FeatureSpace::gene_spans(&ann, selection))
                })
            }
            FeatureSource::Promoters(f) | FeatureSource::GeneActivity(f) => {
                let promoter = matches!(features, FeatureSource::Promoters(_));
                let ann = Arc::clone(self.annotation("gene windows")?);
                let lengths = self.lengths()?;
                let f = *f;
                self.spaces.register(
                    format!(
                        "gene-windows {promoter} {} {} {sel_key}",
                        f.upstream, f.downstream
                    ),
                    || {
                        if let Some(s) = selection {
                            s.check(&ann)?;
                        }
                        Ok(FeatureSpace::gene_windows(
                            &ann, selection, &lengths, promoter, f,
                        ))
                    },
                )
            }
            FeatureSource::BigBed(p) => self.peaks(p),
            FeatureSource::Bins(size) => {
                let refs = self.inputs[0].refs.clone();
                let size = *size;
                self.spaces
                    .register(format!("bins {size}"), || FeatureSpace::bins(&refs, size))
            }
            FeatureSource::Space(s) => {
                let s = Arc::clone(s);
                self.spaces
                    .register(format!("space {:p}", Arc::as_ptr(&s)), || Ok((*s).clone()))
            }
        }
    }

    fn peaks(&mut self, p: &Path) -> Result<SpaceId> {
        let p = p.to_owned();
        self.spaces.register(format!("bigbed {}", p.display()), || {
            FeatureSpace::big_bed(&p, "Peaks")
        })
    }

    fn rows(&self, key: &RowKey, whitelist: Option<&Path>) -> Result<(Rows, u64)> {
        Ok(match (key, whitelist) {
            (RowKey::Tag(t), None) => (Rows::Tag(*t), 0),
            (RowKey::Tag(t), Some(w)) => {
                let names = read_list(w)?;
                if names.is_empty() {
                    return Err(Error::invalid(
                        "whitelist",
                        format!("{} lists no barcodes", w.display()),
                    ));
                }
                let bytes = names.iter().map(|n| 3 * n.len() as u64 + 96).sum();
                (
                    Rows::Listed(*t, Arc::new(index_of(&names)), Arc::new(names)),
                    bytes,
                )
            }
            (RowKey::ReadGroup, w) => {
                let mut names: Vec<Vec<u8>> = match w {
                    Some(w) => read_list(w)?,
                    None => Vec::new(),
                };
                if w.is_none() {
                    for i in self.inputs {
                        for g in &i.read_groups {
                            if !names.iter().any(|n| n == g.as_bytes()) {
                                names.push(g.as_bytes().to_vec());
                            }
                        }
                    }
                }
                if names.is_empty() {
                    return Err(Error::invalid(
                        "quant",
                        "rows by read group, but no @RG line in the headers",
                    ));
                }
                (
                    Rows::ReadGroup(Arc::new(index_of(&names)), Arc::new(names)),
                    0,
                )
            }
            (RowKey::File, _) => (
                Rows::File(Arc::new(
                    self.inputs
                        .iter()
                        .map(|i| i.path.display().to_string().into_bytes())
                        .collect(),
                )),
                0,
            ),
        })
    }

    fn counts(&mut self, c: &Counts) -> Result<CountStage> {
        let space = self.space(c)?;
        let (rows, table_bytes) = self.rows(&c.rows, c.whitelist.as_deref())?;
        let reduce = match c.dedup {
            Dedup::None => Reduce::Sum,
            Dedup::Umi => Reduce::Distinct,
            Dedup::Position => Reduce::MaxPerSub,
        };
        let layered = c.velocity || c.labeled.is_some();
        if layered {
            let genes = matches!(
                c.features,
                FeatureSource::GeneExons | FeatureSource::GeneSpans
            );
            let units = matches!(c.unit, CountUnit::Read | CountUnit::Pair);
            if !genes || !units || c.dedup == Dedup::Position {
                return Err(Error::invalid(
                    "quant",
                    format!(
                        "{}: velocity and labeled layers are counted over genes, by reads, \
                         pairs or UMIs",
                        c.name
                    ),
                ));
            }
        }
        let mut layers = vec!["counts".to_string()];
        let mut fold = Fold::plain(reduce);
        let velocity = if c.velocity {
            fold.spliced = Some(layers.len());
            layers.extend(["spliced", "unspliced", "ambiguous"].map(String::from));
            let mut spans = c.clone();
            spans.features = FeatureSource::GeneSpans;
            let vspace = self.space(&spans)?;
            let ann = self.annotation("velocity")?;
            let genes: Vec<u32> = match &c.subset {
                Some(s) => s.genes.clone(),
                None => (0..ann.genes().len() as u32).collect(),
            };
            Some(Velocity {
                space: vspace,
                models: Arc::new(GeneModels::new(ann, &genes)),
            })
        } else {
            None
        };
        if c.labeled.is_some() {
            fold.labeled = Some(layers.len());
            layers.push("labeled".into());
        }
        Ok(CountStage {
            spec: c.clone(),
            space,
            features: Arc::clone(&self.spaces.get(space).features),
            rows,
            fold,
            layers,
            velocity,
            table_bytes,
        })
    }

    fn junctions(&mut self, j: &Junctions) -> Result<JunctionStage> {
        self.same_references(&format!("{}: junctions", j.name))?;
        let first: Vec<&str> = self.inputs[0].refs.iter().map(|r| r.0.as_str()).collect();
        let (rows, table_bytes, reduce) = match &j.matrix {
            Some(m) => {
                if m.dedup == Dedup::Position {
                    return Err(Error::invalid(
                        "quant",
                        format!("{}: junction matrices count reads, pairs or UMIs", j.name),
                    ));
                }
                let (rows, bytes) = self.rows(&m.rows, m.whitelist.as_deref())?;
                let reduce = match m.dedup {
                    Dedup::Umi => Reduce::Distinct,
                    _ => Reduce::Sum,
                };
                (Some(rows), bytes, reduce)
            }
            None => (None, 0, Reduce::Sum),
        };
        let annotation = self.task.annotation.as_ref().map(|ann| {
            let names: Vec<&str> = ann.seqid_names().collect();
            let seqids = omicsweave_gtf::match_names(&names, &first, self.task.aliases.as_ref());
            (Arc::clone(ann), seqids)
        });
        Ok(JunctionStage {
            spec: j.clone(),
            rows,
            reduce,
            annotation,
            reference: self.reference.clone(),
            table_bytes,
        })
    }

    /// Every input's references, which must be those of the first in the
    /// same order.
    fn same_references(&self, what: &str) -> Result<()> {
        let first = &self.inputs[0].refs;
        if self.inputs.iter().any(|i| &i.refs != first) {
            return Err(Error::invalid(
                "quant",
                format!(
                    "{what} of several BAMs need the same references in the same order (including lengths)"
                ),
            ));
        }
        Ok(())
    }

    fn methylation(&mut self, m: &Methylation) -> Result<MethylStage> {
        if self.inputs.len() > 1 && (m.table || m.bigwig) {
            return Err(Error::invalid(
                "quant",
                format!(
                    "{}: the per-cytosine table and tracks are computed for one BAM at a time; \
                     turn them off to count matrices over several BAMs",
                    m.name
                ),
            ));
        }
        if m.contexts == 0 || m.matrices.iter().any(|x| x.contexts == 0) {
            return Err(Error::invalid(
                "quant",
                format!("{}: no cytosine context to count", m.name),
            ));
        }
        if !(m.min_probability > 0.0 && m.min_probability <= 1.0) {
            return Err(Error::invalid(
                "quant",
                format!("{}: min_probability must be in (0, 1]", m.name),
            ));
        }
        let mut matrices = Vec::new();
        let mut table_bytes = 0;
        let mut names = std::collections::HashSet::new();
        for x in &m.matrices {
            if x.name.is_empty() || !names.insert(x.name.clone()) {
                return Err(Error::invalid(
                    "quant",
                    format!("{}: methylation matrices need distinct names", m.name),
                ));
            }
            let space = self.space_of(&x.features, x.subset.as_ref())?;
            let (rows, bytes) = self.rows(&x.rows, x.whitelist.as_deref())?;
            table_bytes += bytes;
            matrices.push(MatrixStage {
                spec: x.clone(),
                space,
                features: Arc::clone(&self.spaces.get(space).features),
                rows,
            });
        }
        if matrices.len() > u8::MAX as usize {
            return Err(Error::invalid("quant", "too many methylation matrices"));
        }
        Ok(MethylStage {
            spec: m.clone(),
            matrices,
            table_bytes,
        })
    }

    fn contacts(&mut self, c: &Contacts) -> Result<ContactStage> {
        self.same_references(&format!("{}: contacts", c.name))?;
        if c.resolutions.contains(&0) {
            return Err(Error::invalid(
                "quant",
                format!("{}: a .hic resolution must be at least 1 bp", c.name),
            ));
        }
        let (rows, table_bytes) = match c.barcode_tag {
            Some(t) => {
                let (r, b) = self.rows(&RowKey::Tag(t), c.whitelist.as_deref())?;
                (Some(r), b)
            }
            None => (None, 0),
        };
        Ok(ContactStage {
            spec: c.clone(),
            rows,
            table_bytes,
        })
    }

    fn alleles(&mut self, a: &Alleles) -> Result<AlleleStage> {
        let phasing = a
            .genes
            .as_ref()
            .map(|g| g.sample.clone().unwrap_or_default());
        let panel = SitePanel::open(
            &a.sites,
            self.inputs,
            self.task.aliases.as_ref(),
            phasing.as_deref(),
        )?;
        let (rows, table_bytes) = self.rows(&a.rows, a.whitelist.as_deref())?;
        let genes = match &a.genes {
            Some(_) => {
                let space = self.space_of(&FeatureSource::GeneExons, None)?;
                Some((space, Arc::clone(&self.spaces.get(space).features)))
            }
            None => None,
        };
        Ok(AlleleStage {
            spec: a.clone(),
            panel: Arc::new(panel),
            rows,
            genes,
            table_bytes,
        })
    }

    fn atac(&mut self, a: &Atac) -> Result<(FragmentStage, Vec<CountStage>)> {
        let has_ann = self.task.annotation.is_some();
        let peaks = match &a.peaks {
            Some(p) => Some(self.peaks(p)?),
            None => None,
        };
        let tss = if has_ann {
            let ann = Arc::clone(self.annotation("TSS metrics")?);
            let lengths = self.lengths()?;
            let centre = self.spaces.register("tss centre".into(), || {
                Ok(FeatureSpace::tss_windows(
                    &ann,
                    &lengths,
                    Flank::new(50, 51),
                ))
            })?;
            let flank = self
                .spaces
                .register("tss flanks".into(), || tss_flanks(&ann, &lengths))?;
            Some((centre, flank))
        } else {
            None
        };
        let (rows, table_bytes) = self.rows(&RowKey::Tag(a.barcode_tag), a.whitelist.as_deref())?;
        let mut counts = Vec::new();
        for c in a.matrices(has_ann) {
            counts.push(self.counts(&c)?);
        }
        Ok((
            FragmentStage {
                spec: a.clone(),
                rows,
                peaks,
                tss,
                table_bytes,
            },
            counts,
        ))
    }
}

/// The two 100 bp windows 1900-2000 bp from each gene TSS (ArchR's TSS
/// enrichment flanks).
fn tss_flanks(ann: &Annotation, lengths: &[Option<u32>]) -> Result<FeatureSpace> {
    use omicsweave_gtf::{IntervalSetBuilder, NONE, Strand};
    let tss = ann.derive().clip(lengths).tss(omicsweave_gtf::Level::Gene);
    let names: Vec<String> = ann.seqid_names().map(str::to_owned).collect();
    let mut b = IntervalSetBuilder::new().chroms(names.iter().map(String::as_str))?;
    let mut t = crate::model::FeatureTable::new("TSS flank");
    for r in 0..tss.len() {
        let p = tss.start[r];
        let chrom = &names[tss.seqid[r] as usize];
        let len = lengths[tss.seqid[r] as usize].unwrap_or(u32::MAX);
        let windows = [
            (p.checked_sub(2000), p.checked_sub(1900)),
            (p.checked_add(1901), p.checked_add(2001)),
        ];
        for (s, e) in windows {
            if let (Some(s), Some(e)) = (s, e) {
                let e = e.min(len);
                if s < e {
                    b.push(chrom, s, e, Strand::Unstranded, NONE, None)?;
                    t.push("flank", "flank", chrom, s, e, '.', (e - s) as u64);
                }
            }
        }
    }
    let set = b.build();
    let n = set.len() as u32;
    FeatureSpace::from_set(set, (0..n).collect(), t, names)
}

/// Checks `task` against `budget` and divides the memory; reads the
/// headers of the inputs and builds the feature spaces, nothing else.
pub fn plan(task: Task, budget: Budget) -> Result<Plan> {
    if budget.threads == 0 {
        return Err(Error::invalid("budget", "threads must be at least one"));
    }
    if task.inputs.is_empty() {
        return Err(Error::invalid("quant", "no input BAM"));
    }
    if task.products.is_empty() {
        return Err(Error::invalid(
            "quant",
            "nothing to compute: ask for coverage, counts, ATAC, junction, labeling, \
             methylation, contact or allele products",
        ));
    }
    if task.out_dir.is_file() {
        return Err(Error::invalid(
            "output",
            format!("{} is a file, not a directory", task.out_dir.display()),
        ));
    }
    let inputs = task
        .inputs
        .iter()
        .map(|p| Input::open(p))
        .collect::<Result<Vec<_>>>()?;
    let single = inputs.len() == 1;
    let reference = match &task.reference {
        Some(p) => Some(Arc::new(Reference::open(
            p,
            &inputs,
            task.aliases.as_ref(),
        )?)),
        None => None,
    };
    let variants = match &task.variants {
        Some(p) => Some(Arc::new(Variants::open(p, &inputs, task.aliases.as_ref())?)),
        None => None,
    };
    let mut r = Resolver {
        task: &task,
        inputs: &inputs,
        spaces: Spaces::default(),
        lengths: None,
        reference: reference.clone(),
    };
    let mut stages: Vec<Box<dyn Stage>> = Vec::new();
    for p in &task.products {
        match p {
            Product::Coverage(c) => {
                if !single {
                    return Err(Error::invalid(
                        "quant",
                        "coverage is computed for one BAM at a time",
                    ));
                }
                stages.push(Box::new(Erased(c.clone())));
            }
            Product::Counts(c) => stages.push(Box::new(Erased(r.counts(c)?))),
            Product::Atac(a) => {
                if !single {
                    return Err(Error::invalid(
                        "quant",
                        "ATAC products are computed for one BAM at a time",
                    ));
                }
                let (f, counts) = r.atac(a)?;
                stages.push(Box::new(Erased(f)));
                for c in counts {
                    stages.push(Box::new(Erased(c)));
                }
            }
            Product::Junctions(j) => stages.push(Box::new(Erased(r.junctions(j)?))),
            Product::Labeling(l) => stages.push(Box::new(Erased(l.clone()))),
            Product::Methylation(m) => stages.push(Box::new(Erased(r.methylation(m)?))),
            Product::Contacts(c) => stages.push(Box::new(Erased(r.contacts(c)?))),
            Product::Alleles(a) => stages.push(Box::new(Erased(r.alleles(a)?))),
        }
    }
    let mut seen = std::collections::HashSet::new();
    for s in &stages {
        if !seen.insert(s.name().to_owned()) {
            return Err(Error::invalid(
                "quant",
                format!("two products are called {:?}; name them apart", s.name()),
            ));
        }
    }
    let mut spaces = r.spaces;
    let refs: Vec<Vec<String>> = inputs
        .iter()
        .map(|i| i.refs.iter().map(|r| r.0.clone()).collect())
        .collect();
    spaces.match_inputs(&refs, task.aliases.as_ref());

    let threads = budget.threads;
    let batch_bytes = text_batch(&budget);
    let inflate = sizing(
        threads,
        2 * omicsweave_bgzf::MAX_BLOCK_LEN as u64,
        budget.memory / 8,
    );
    let probe = PlanInfo {
        threads,
        inputs: inputs.iter().map(|i| (i.refs.clone(), i.len)).collect(),
        batch_bytes: batch_bytes as u64,
        batches: 1,
    };
    let needs_of: Vec<MemoryNeed> = stages.iter().map(|s| s.estimate(&probe)).collect();
    let mut needs = Needs::NONE;
    for s in &stages {
        needs |= s.needs();
    }
    // Each batch in flight holds a window of reference bases (and masked
    // positions) about its records.
    let windows = (needs.contains(Needs::REFERENCE) && reference.is_some()) as u64 * (2 << 20)
        + (needs.contains(Needs::VARIANTS) && variants.is_some()) as u64 * (1 << 20);
    let per_batch: u64 =
        batch_bytes as u64 + windows + needs_of.iter().map(|n| n.per_batch).sum::<u64>();
    let batches = sizing(threads, per_batch, budget.memory / 4);
    let fixed = needs_of.iter().map(|n| n.fixed).sum::<u64>()
        + spaces.heap_bytes()
        + task
            .annotation
            .as_ref()
            .map_or(0, |a| a.heap_bytes() as u64);
    let pipeline = batches.per_queue + inflate.per_queue + SCAN_TEXT + READ_BUF;
    let min_buffer: u64 = needs_of.iter().map(|n| n.min_buffer).sum();
    let need = fixed.saturating_add(pipeline).saturating_add(min_buffer);
    if need > budget.memory {
        return Err(Error::invalid(
            "budget",
            format!(
                "this task needs at least {} of memory (fixed {}, pipeline {}, buffers {}); \
                 the budget is {}",
                human(need),
                human(fixed),
                human(pipeline),
                human(min_buffer),
                human(budget.memory)
            ),
        ));
    }
    let buffer = budget.memory - fixed - pipeline;
    let spare = buffer - min_buffer;
    let weights: u64 = needs_of.iter().map(|n| n.weight as u64).sum::<u64>().max(1);
    let shares: Vec<u64> = needs_of
        .iter()
        .map(|n| n.min_buffer + spare / weights * n.weight as u64)
        .collect();
    let names: Vec<(String, u64)> = stages
        .iter()
        .zip(&shares)
        .map(|(s, &b)| (s.name().to_owned(), b))
        .collect();
    let n_refs = inputs.iter().map(|i| i.refs.len()).collect();
    Ok(Plan {
        budget,
        memory: MemoryPlan {
            fixed,
            pipeline,
            buffer,
        },
        needs,
        stages: names,
        inputs,
        shared: Arc::new(Shared {
            stages,
            spaces,
            n_refs,
            reference: reference.filter(|_| needs.contains(Needs::REFERENCE)),
            variants: variants.filter(|_| needs.contains(Needs::VARIANTS)),
        }),
        batch_bytes,
        batches: batches.limit,
        inflate,
        task,
    })
}

/// A batch's states, sealed, with its frontier.
type Sealed = (Vec<Box<dyn Any + Send>>, u64);

fn process(shared: &Shared, batch: crate::source::Batch) -> Result<Sealed> {
    let input = batch.input as usize;
    let ctx = BatchCtx {
        spaces: BatchSpaces::new(&shared.spaces, input),
        reference: shared
            .reference
            .as_deref()
            .map(|r| RefCell::new(RefWindow::new(r, input))),
        variants: shared
            .variants
            .as_deref()
            .map(|v| RefCell::new(MaskWindow::new(v, input))),
    };
    let info = BatchInfo { input: batch.input };
    let mut states: Vec<Box<dyn Any + Send>> =
        shared.stages.iter().map(|s| s.batch(info)).collect();
    let n_refs = shared.n_refs[input];
    let data = &batch.data;
    let mut at = 0usize;
    let mut ordinal = batch.first;
    while at + 4 <= data.len() {
        let len = u32::from_le_bytes(data[at..at + 4].try_into().unwrap()) as usize;
        let rec = Record::parse(&data[at + 4..at + 4 + len], n_refs)?;
        let view = ReadView::new(batch.input, ordinal, rec, &ctx);
        for (s, st) in shared.stages.iter().zip(states.iter_mut()) {
            s.observe(st.as_mut(), &view)?;
        }
        at += 4 + len;
        ordinal += 1;
    }
    for (s, st) in shared.stages.iter().zip(states.iter_mut()) {
        s.seal(st.as_mut())?;
    }
    Ok((states, batch.frontier))
}

fn merge_all(m: &mut [Box<dyn Merging + '_>], sealed: Sealed) -> Result<()> {
    let (states, frontier) = sealed;
    for (mg, st) in m.iter_mut().zip(states) {
        mg.merge(st, frontier)?;
    }
    Ok(())
}

/// Runs a plan: one pass over each input, then every product is written;
/// files appear only once all of them are.
pub fn run(plan: &Plan, progress: OnProgress<'_>) -> Result<Report> {
    let start = Instant::now();
    let out_dir = &plan.task.out_dir;
    std::fs::create_dir_all(out_dir)?;
    let tmp = plan
        .budget
        .tmp_dir
        .clone()
        .unwrap_or_else(|| out_dir.clone());
    let pool = Pool::new(plan.budget.threads);
    let total: u64 = plan.inputs.iter().map(|i| i.len).sum();
    let report = |phase: &'static str, read: u64| {
        if let Some(f) = progress {
            f(Progress {
                phase,
                bytes_read: read,
                bytes_total: total,
            });
        }
    };
    let shared = &plan.shared;
    let mut merging: Vec<Box<dyn Merging + '_>> = Vec::new();
    for (s, (_, share)) in shared.stages.iter().zip(&plan.stages) {
        merging.push(s.start(RunCtx {
            pool: &pool,
            buffer: *share,
            tmp: tmp.clone(),
            inputs: &plan.inputs,
        })?);
    }
    let mut records = 0u64;
    let mut read_before = 0u64;
    let mut per_input = Vec::new();
    let mut last_report = Instant::now();
    for (i, input) in plan.inputs.iter().enumerate() {
        let mut reader = BatchReader::open(
            input,
            i as u32,
            Box::new(pool.ordered(plan.inflate.limit)),
            plan.batch_bytes,
        )?;
        let mut q = pool.ordered::<Sealed>(plan.batches);
        let mut n = 0u64;
        while let Some(batch) = reader.next_batch()? {
            n += batch.records as u64;
            while q.full() {
                merge_all(&mut merging, q.next(true)?.expect("pending"))?;
            }
            let sh = Arc::clone(shared);
            q.submit(Box::new(move || process(&sh, batch)))?;
            while let Some(done) = q.next(false)? {
                merge_all(&mut merging, done)?;
            }
            if last_report.elapsed() > Duration::from_millis(200) {
                report("read", read_before + reader.bytes_read());
                last_report = Instant::now();
            }
        }
        while let Some(done) = q.next(true)? {
            merge_all(&mut merging, done)?;
        }
        for m in &mut merging {
            m.end_input(i as u32)?;
        }
        read_before += reader.bytes_read();
        records += n;
        per_input.push(Json::obj([
            ("path", Json::from(input.path.display().to_string())),
            ("placed_records", Json::from(n)),
        ]));
    }
    let scan = start.elapsed();
    report("write", read_before);
    // Products are written side by side, each to files of its own; they
    // are published together once every one succeeded.
    let share = plan.memory.buffer / merging.len().max(1) as u64;
    let delivered: Vec<Result<(Json, Outputs<&'static str>)>> = std::thread::scope(|s| {
        let handles: Vec<_> = merging
            .into_iter()
            .map(|m| {
                let (pool, tmp) = (&pool, &tmp);
                s.spawn(move || {
                    let mut outputs = Outputs::new();
                    let mut d = Delivery {
                        outputs: &mut outputs,
                        dir: out_dir,
                        pool,
                        matrix: plan.task.matrix,
                        tmp,
                        buffer: share,
                    };
                    let summary = m.finish(&mut d)?;
                    Ok((summary, outputs))
                })
            })
            .collect();
        handles
            .into_iter()
            .map(|h| {
                h.join()
                    .unwrap_or_else(|_| Err(Error::invalid("quant", "writing a product panicked")))
            })
            .collect()
    });
    let mut products = Vec::new();
    let mut written = Vec::new();
    for d in delivered {
        let (summary, outputs) = d?;
        products.push(summary);
        written.push(outputs);
    }
    let mut outputs: Outputs<&'static str> = Outputs::new();
    let summary = Json::obj([
        ("omicsweave", Json::from(env!("CARGO_PKG_VERSION"))),
        ("inputs", Json::Arr(per_input)),
        (
            "needs",
            Json::Arr(plan.needs.names().into_iter().map(Json::from).collect()),
        ),
        ("products", Json::Arr(products)),
    ]);
    outputs.write(
        &out_dir.join("summary.json"),
        "json",
        summary.pretty().as_bytes(),
    )?;
    let mut files = Vec::new();
    for o in written {
        files.extend(o.commit()?);
    }
    files.extend(outputs.commit()?);
    Ok(Report {
        files,
        records,
        bytes_read: read_before,
        elapsed: start.elapsed(),
        scan,
        summary,
    })
}

/// [`plan`] then [`run`].
pub fn quantify(task: Task, budget: Budget, progress: OnProgress<'_>) -> Result<Report> {
    run(&plan(task, budget)?, progress)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::accum::Accumulator;

    #[test]
    fn products_share_one_space_per_gene_subset() {
        let ann = Annotation::parse(
            b"c\t.\texon\t1\t5\t.\t+\t.\tgene_id \"a\"; transcript_id \"a1\";\n\
              c\t.\texon\t20\t25\t.\t+\t.\tgene_id \"b\"; transcript_id \"b1\";\n",
            &omicsweave_gtf::ReadOptions {
                threads: 1,
                ..Default::default()
            },
        )
        .unwrap();
        let mut task = Task::new(vec![], "unused");
        task.annotation = Some(Arc::new(ann));
        let mut resolver = Resolver {
            task: &task,
            inputs: &[],
            spaces: Spaces::default(),
            lengths: None,
            reference: None,
        };
        let subset = Arc::new(Selection {
            genes: vec![0],
            transcripts: vec![0],
        });
        let first = resolver
            .space_of(&FeatureSource::GeneExons, Some(&subset))
            .unwrap();
        let again = resolver
            .space_of(&FeatureSource::GeneExons, Some(&Arc::clone(&subset)))
            .unwrap();
        assert_eq!(first, again);
        assert_eq!(resolver.spaces.get(first).features.len(), 1);
        let all = resolver.space_of(&FeatureSource::GeneExons, None).unwrap();
        assert_ne!(first, all);
        assert_eq!(resolver.spaces.get(all).features.len(), 2);
        assert!(
            resolver
                .space_of(&FeatureSource::Bins(10), Some(&subset))
                .is_err()
        );
    }

    #[test]
    fn fragment_only_whitelist_is_included_in_fixed_memory() {
        let dir = tempfile::tempdir().unwrap();
        let whitelist = dir.path().join("barcodes");
        let mut names = String::new();
        for i in 0..10_000 {
            names.push_str(&format!("B{i:08}\n"));
        }
        std::fs::write(&whitelist, names).unwrap();
        let task = Task::new(vec![], "unused");
        let mut resolver = Resolver {
            task: &task,
            inputs: &[],
            spaces: Spaces::default(),
            lengths: None,
            reference: None,
        };
        let spec = Atac {
            whitelist: Some(whitelist),
            tiles: None,
            gene_activity: None,
            ..Atac::default()
        };
        let (stage, counts) = resolver.atac(&spec).unwrap();
        assert!(counts.is_empty());
        let memory = stage.estimate(&PlanInfo {
            threads: 1,
            inputs: vec![],
            batch_bytes: 4096,
            batches: 1,
        });
        assert!(memory.fixed > 1_000_000, "fixed={}", memory.fixed);
        assert!(memory.fixed >= 10_000 * size_of::<[u64; 6]>() as u64);
    }

    #[test]
    fn combined_coordinate_outputs_check_names_lengths_and_order() {
        let task = Task::new(vec![], "unused");
        let dictionary = vec![("chr1".to_string(), 100), ("chr2".to_string(), 200)];
        for (other, valid) in [
            (dictionary.clone(), true),
            (vec![("chr1".into(), 101), ("chr2".into(), 200)], false),
            (vec![("chr2".into(), 200), ("chr1".into(), 100)], false),
            (vec![("chr1".into(), 100)], false),
        ] {
            let inputs = [
                Input {
                    path: PathBuf::from("one"),
                    refs: dictionary.clone(),
                    read_groups: vec![],
                    len: 0,
                },
                Input {
                    path: PathBuf::from("two"),
                    refs: other,
                    read_groups: vec![],
                    len: 0,
                },
            ];
            let mut resolver = Resolver {
                task: &task,
                inputs: &inputs,
                spaces: Spaces::default(),
                lengths: None,
                reference: None,
            };
            assert_eq!(resolver.same_references("contacts").is_ok(), valid);
            assert_eq!(resolver.junctions(&Junctions::default()).is_ok(), valid);
        }
    }
}
