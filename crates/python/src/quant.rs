//! Quantification. The product presets are this crate's: [`defaults`]
//! hands them to Python as plain values, and products come back as the
//! package's dataclasses, read attribute by attribute; a missing or
//! mistyped attribute is an error, never a silent default. A task is
//! planned once ([`QuantPlan`]) and the plan is what runs.

use std::collections::HashMap;
use std::path::PathBuf;
use std::sync::Arc;

use numpy::PyReadonlyArray1;
use omicsweave_core::Result;
use omicsweave_gtf::{Aliases, Flank, Selection};
use omicsweave_quant::{
    AlleleGenes, Alleles, Atac, Budget, CHG, CHH, CPG, CallSource, Contacts, Conversion, CountUnit,
    Counts, Coverage, CoverageUnit, Dedup, FeatureSource, FeatureSpace, Filter, JunctionMatrix,
    Junctions, Labeled, Labeling, MatrixFormat, MethylMatrix, Methylation, Plan, Product, Progress,
    Report, RowKey, Scale, StrandSplit, Strandedness, Task, contexts_from_names, plan, run,
};
use pyo3::IntoPyObjectExt;
use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;
use pyo3::types::{PyDict, PyList, PyTuple};

use crate::annotation::{GeneAnnotation, IntervalSet};
use crate::error::OrRaise;
use crate::shared::{PlanRunState, path_str};

// ---------------------------------------------------------------------------
// Python objects to products

fn attr<'py, T: FromPyObjectOwned<'py>>(o: &Bound<'py, PyAny>, name: &str) -> PyResult<T> {
    o.getattr(name)?.extract::<T>().map_err(Into::into)
}

fn invalid(detail: String) -> PyErr {
    PyValueError::new_err(detail)
}

fn tag(s: &str) -> PyResult<[u8; 2]> {
    match s.as_bytes() {
        &[a, b] => Ok([a, b]),
        _ => Err(invalid(format!("a tag has two characters, not {s:?}"))),
    }
}

fn tag_name(t: [u8; 2]) -> String {
    String::from_utf8_lossy(&t).into_owned()
}

fn strandedness(s: &str) -> PyResult<Strandedness> {
    Ok(match s {
        "unstranded" => Strandedness::Unstranded,
        "forward" => Strandedness::Forward,
        "reverse" => Strandedness::Reverse,
        _ => {
            return Err(invalid(format!(
                "strand must be unstranded, forward or reverse, not {s:?}"
            )));
        }
    })
}

fn strandedness_name(s: Strandedness) -> &'static str {
    match s {
        Strandedness::Unstranded => "unstranded",
        Strandedness::Forward => "forward",
        Strandedness::Reverse => "reverse",
    }
}

fn strand_of(o: &Bound<'_, PyAny>) -> PyResult<Strandedness> {
    strandedness(&attr::<String>(o, "strand")?)
}

fn filter(o: &Bound<'_, PyAny>) -> PyResult<Filter> {
    let f = o.getattr("filter")?;
    Ok(Filter {
        require: attr(&f, "require")?,
        exclude: attr(&f, "exclude")?,
        min_mapq: attr(&f, "min_mapq")?,
        unique: attr(&f, "unique")?,
    })
}

fn row_key(s: &str) -> PyResult<RowKey> {
    Ok(match s {
        "@RG" => RowKey::ReadGroup,
        "@file" => RowKey::File,
        t => RowKey::Tag(tag(t)?),
    })
}

fn rows(o: &Bound<'_, PyAny>) -> PyResult<RowKey> {
    row_key(&attr::<String>(o, "rows")?)
}

fn dedup(o: &Bound<'_, PyAny>) -> PyResult<Dedup> {
    Ok(match attr::<String>(o, "dedup")?.as_str() {
        "none" => Dedup::None,
        "umi" => Dedup::Umi,
        "position" => Dedup::Position,
        d => {
            return Err(invalid(format!(
                "dedup must be none, umi or position, not {d:?}"
            )));
        }
    })
}

fn flank(o: &Bound<'_, PyAny>) -> PyResult<Flank> {
    Ok(Flank::new(attr(o, "upstream")?, attr(o, "downstream")?))
}

fn interval_set(o: &Bound<'_, PyAny>) -> PyResult<Arc<omicsweave_gtf::IntervalSet>> {
    let native = match o.cast::<IntervalSet>() {
        Ok(set) => set.clone(),
        Err(_) => o.getattr("_native")?.cast_into::<IntervalSet>()?,
    };
    native.get().shared()
}

fn features(o: &Bound<'_, PyAny>) -> PyResult<FeatureSource> {
    let f = o.getattr("features")?;
    Ok(match attr::<String>(&f, "kind")?.as_str() {
        "gene-exons" => FeatureSource::GeneExons,
        "gene-spans" => FeatureSource::GeneSpans,
        "promoters" => FeatureSource::Promoters(flank(&f)?),
        "gene-activity" => FeatureSource::GeneActivity(flank(&f)?),
        "peaks" => FeatureSource::BigBed(attr(&f, "path")?),
        "bins" => FeatureSource::Bins(attr(&f, "size")?),
        "regions" => {
            let set = interval_set(&f.getattr("intervals")?)?;
            let kind: String = attr(&f, "feature_type")?;
            let space = o
                .py()
                .detach(|| FeatureSpace::of_intervals((*set).clone(), &kind))
                .or_raise()?;
            FeatureSource::Space(Arc::new(space))
        }
        k => return Err(invalid(format!("unknown features {k:?}"))),
    })
}

/// The `subset` of a product: an `annotation.Selection` (its gene and
/// transcript rows), checked against the task's annotation when planned.
type Subsets = HashMap<usize, (Py<PyAny>, Arc<Selection>)>;

fn subset(o: &Bound<'_, PyAny>, subsets: &mut Subsets) -> PyResult<Option<Arc<Selection>>> {
    let s = o.getattr("subset")?;
    if s.is_none() {
        return Ok(None);
    }
    let key = s.as_ptr() as usize;
    if let Some((_, rows)) = subsets.get(&key) {
        return Ok(Some(Arc::clone(rows)));
    }
    let rows = |name: &str| -> PyResult<Vec<u32>> {
        let a = s.getattr(name)?;
        let a: PyReadonlyArray1<'_, u32> = a.extract().map_err(|_| {
            PyValueError::new_err(format!("subset.{name} must be a uint32 array of rows"))
        })?;
        Ok(a.as_slice()?.to_vec())
    };
    let rows = Arc::new(Selection {
        genes: rows("genes")?,
        transcripts: rows("transcripts")?,
    });
    // Retain the Python owner so its address cannot be reused while this task is parsed.
    subsets.insert(key, (s.unbind(), Arc::clone(&rows)));
    Ok(Some(rows))
}

fn conversion(o: &Bound<'_, PyAny>) -> PyResult<Conversion> {
    Conversion::parse(&attr::<String>(o, "conversion")?, attr(o, "min_quality")?).or_raise()
}

fn contexts(o: &Bound<'_, PyAny>) -> PyResult<u8> {
    contexts_from_names(&attr::<Vec<String>>(o, "contexts")?).or_raise()
}

fn coverage(o: &Bound<'_, PyAny>) -> PyResult<Coverage> {
    let deletions: bool = attr(o, "deletions")?;
    let max_fragment: u32 = attr(o, "max_fragment")?;
    let unit = match attr::<String>(o, "unit")?.as_str() {
        "bases" => CoverageUnit::Bases { deletions },
        "span" => CoverageUnit::Span,
        "five-prime" => CoverageUnit::FivePrime,
        "three-prime" => CoverageUnit::ThreePrime,
        "fragment" => CoverageUnit::Fragment {
            max_len: max_fragment,
        },
        u => return Err(invalid(format!("unknown coverage unit {u:?}"))),
    };
    if deletions && !matches!(unit, CoverageUnit::Bases { .. }) {
        return Err(invalid("deletions apply to the bases unit".into()));
    }
    if max_fragment != 0 && !matches!(unit, CoverageUnit::Fragment { .. }) {
        return Err(invalid("max_fragment applies to the fragment unit".into()));
    }
    let strand = strand_of(o)?;
    let strands = match attr::<String>(o, "strands")?.as_str() {
        "merged" if strand == Strandedness::Unstranded => StrandSplit::Merged,
        "merged" => {
            return Err(invalid(
                "strand applies to split strands; merged coverage has none".into(),
            ));
        }
        "split" => StrandSplit::Split(strand),
        s => {
            return Err(invalid(format!(
                "strands must be merged or split, not {s:?}"
            )));
        }
    };
    let scale = o.getattr("scale")?;
    let scale = match scale.extract::<String>() {
        Ok(s) if s == "none" => Scale::None,
        Ok(s) if s == "cpm" => Scale::Cpm,
        Ok(s) => {
            return Err(invalid(format!(
                "scale must be none, cpm or a number, not {s:?}"
            )));
        }
        Err(_) => Scale::Factor(scale.extract::<f64>()?),
    };
    Ok(Coverage {
        name: attr(o, "name")?,
        unit,
        strands,
        filter: filter(o)?,
        scale,
        sidecar: attr(o, "sidecar")?,
    })
}

fn labeled(o: &Bound<'_, PyAny>) -> PyResult<Option<Labeled>> {
    let l = o.getattr("labeled")?;
    if l.is_none() {
        return Ok(None);
    }
    Ok(Some(Labeled {
        conversion: conversion(&l)?,
        min_conversions: attr(&l, "min_conversions")?,
    }))
}

fn counts(o: &Bound<'_, PyAny>, subsets: &mut Subsets) -> PyResult<Counts> {
    let unit = match attr::<String>(o, "unit")?.as_str() {
        "read" => CountUnit::Read,
        "pair" => CountUnit::Pair,
        "fragment" => CountUnit::Fragment { tn5: false },
        "tn5-fragment" => CountUnit::Fragment { tn5: true },
        "cut-site" => CountUnit::CutSite,
        u => return Err(invalid(format!("unknown count unit {u:?}"))),
    };
    Ok(Counts {
        name: attr(o, "name")?,
        features: features(o)?,
        subset: subset(o, subsets)?,
        unit,
        rows: rows(o)?,
        whitelist: attr(o, "whitelist")?,
        umi_tag: tag(&attr::<String>(o, "umi_tag")?)?,
        dedup: dedup(o)?,
        strand: strand_of(o)?,
        multi_overlap: attr(o, "multi_overlap")?,
        filter: filter(o)?,
        max_fragment: attr(o, "max_fragment")?,
        table: attr(o, "table")?,
        velocity: attr(o, "velocity")?,
        labeled: labeled(o)?,
    })
}

fn atac(o: &Bound<'_, PyAny>) -> PyResult<Atac> {
    let activity = o.getattr("gene_activity")?;
    Ok(Atac {
        name: attr(o, "name")?,
        barcode_tag: tag(&attr::<String>(o, "barcode_tag")?)?,
        whitelist: attr(o, "whitelist")?,
        filter: filter(o)?,
        max_fragment: attr(o, "max_fragment")?,
        fragments: attr(o, "fragments")?,
        peaks: attr(o, "peaks")?,
        tiles: attr(o, "tiles")?,
        gene_activity: if activity.is_none() {
            None
        } else {
            Some(flank(&activity)?)
        },
    })
}

fn junction_matrix(o: &Bound<'_, PyAny>) -> PyResult<JunctionMatrix> {
    let dedup = dedup(o)?;
    if dedup == Dedup::Position {
        return Err(invalid(
            "junction matrices deduplicate by none or umi, not position".into(),
        ));
    }
    Ok(JunctionMatrix {
        rows: rows(o)?,
        whitelist: attr(o, "whitelist")?,
        umi_tag: tag(&attr::<String>(o, "umi_tag")?)?,
        dedup,
        filter: filter(o)?,
    })
}

fn junctions(o: &Bound<'_, PyAny>) -> PyResult<Junctions> {
    let matrix = o.getattr("matrix")?;
    Ok(Junctions {
        name: attr(o, "name")?,
        filter: filter(o)?,
        pairs: attr(o, "pairs")?,
        matrix: if matrix.is_none() {
            None
        } else {
            Some(junction_matrix(&matrix)?)
        },
    })
}

fn labeling(o: &Bound<'_, PyAny>) -> PyResult<Labeling> {
    Ok(Labeling {
        name: attr(o, "name")?,
        filter: filter(o)?,
        strand: strand_of(o)?,
        conversion: conversion(o)?,
    })
}

fn methyl_matrix(o: &Bound<'_, PyAny>, subsets: &mut Subsets) -> PyResult<MethylMatrix> {
    Ok(MethylMatrix {
        name: attr(o, "name")?,
        features: features(o)?,
        subset: subset(o, subsets)?,
        rows: rows(o)?,
        whitelist: attr(o, "whitelist")?,
        contexts: contexts(o)?,
    })
}

fn methylation(o: &Bound<'_, PyAny>, subsets: &mut Subsets) -> PyResult<Methylation> {
    let matrices = o
        .getattr("matrices")?
        .try_iter()?
        .map(|m| methyl_matrix(&m?, subsets))
        .collect::<PyResult<Vec<_>>>()?;
    Ok(Methylation {
        name: attr(o, "name")?,
        source: CallSource::from_name(&attr::<String>(o, "source")?).or_raise()?,
        contexts: contexts(o)?,
        merge_cpg: attr(o, "merge_cpg")?,
        filter: filter(o)?,
        min_quality: attr(o, "min_quality")?,
        min_probability: attr(o, "min_probability")?,
        overlap: attr(o, "overlap")?,
        table: attr(o, "table")?,
        bigwig: attr(o, "bigwig")?,
        matrices,
    })
}

fn contacts(o: &Bound<'_, PyAny>) -> PyResult<Contacts> {
    Ok(Contacts {
        name: attr(o, "name")?,
        filter: filter(o)?,
        min_distance: attr(o, "min_distance")?,
        self_ligation: attr(o, "self_ligation")?,
        dedup: attr(o, "dedup")?,
        pairs: attr(o, "pairs")?,
        resolutions: attr(o, "resolutions")?,
        genome: attr(o, "genome")?,
        barcode_tag: attr::<Option<String>>(o, "barcode_tag")?
            .map(|t| tag(&t))
            .transpose()?,
        whitelist: attr(o, "whitelist")?,
        sidecar: attr(o, "sidecar")?,
    })
}

fn alleles(o: &Bound<'_, PyAny>) -> PyResult<Alleles> {
    let umi = match attr::<String>(o, "unit")?.as_str() {
        "molecule" => true,
        "fragment" => false,
        u => {
            return Err(invalid(format!(
                "allele units are molecule or fragment, not {u:?}"
            )));
        }
    };
    let genes = o.getattr("genes")?;
    Ok(Alleles {
        name: attr(o, "name")?,
        sites: attr(o, "sites")?,
        rows: rows(o)?,
        whitelist: attr(o, "whitelist")?,
        umi,
        umi_tag: tag(&attr::<String>(o, "umi_tag")?)?,
        filter: filter(o)?,
        min_quality: attr(o, "min_quality")?,
        table: attr(o, "table")?,
        genes: if genes.is_none() {
            None
        } else {
            Some(AlleleGenes {
                sample: attr(&genes, "sample")?,
                strand: strand_of(&genes)?,
            })
        },
    })
}

fn product(o: &Bound<'_, PyAny>, subsets: &mut Subsets) -> PyResult<Product> {
    Ok(match attr::<String>(o, "kind")?.as_str() {
        "coverage" => Product::Coverage(coverage(o)?),
        "counts" => Product::Counts(counts(o, subsets)?),
        "atac" => Product::Atac(atac(o)?),
        "junctions" => Product::Junctions(junctions(o)?),
        "labeling" => Product::Labeling(labeling(o)?),
        "methylation" => Product::Methylation(methylation(o, subsets)?),
        "contacts" => Product::Contacts(contacts(o)?),
        "alleles" => Product::Alleles(alleles(o)?),
        k => return Err(invalid(format!("unknown product kind {k:?}"))),
    })
}

/// The task a `quant.Task` describes.
fn task(t: &Bound<'_, PyAny>) -> PyResult<Task> {
    let mut task = Task::new(attr(t, "inputs")?, attr::<PathBuf>(t, "out_dir")?);
    let annotation = t.getattr("annotation")?;
    if !annotation.is_none() {
        let native = annotation
            .getattr("_native")?
            .cast_into::<GeneAnnotation>()?;
        task.annotation = Some(native.get().shared()?);
    }
    let aliases: Option<Vec<(String, String)>> = {
        let a = t.getattr("aliases")?;
        if a.is_none() {
            None
        } else {
            Some(
                a.call_method0("items")?
                    .try_iter()?
                    .map(|kv| kv?.extract::<(String, String)>())
                    .collect::<PyResult<_>>()?,
            )
        }
    };
    task.aliases = aliases.map(Aliases::from_pairs);
    task.matrix = MatrixFormat::from_name(&attr::<String>(t, "matrix")?).or_raise()?;
    task.reference = attr(t, "reference")?;
    task.variants = attr(t, "variants")?;
    let mut subsets = Subsets::new();
    for p in t.getattr("products")?.try_iter()? {
        task.products.push(product(&p?, &mut subsets)?);
    }
    Ok(task)
}

/// The budget of a run: `None` threads are every core the process may use,
/// `None` memory is the default budget.
pub fn budget(
    threads: Option<usize>,
    memory: Option<u64>,
    tmp_dir: Option<PathBuf>,
) -> PyResult<Budget> {
    let mut b = Budget::default();
    if let Some(t) = threads {
        if t == 0 {
            return Err(invalid("threads must be at least one".into()));
        }
        b.threads = t;
    }
    if let Some(m) = memory {
        b.memory = m;
    }
    b.tmp_dir = tmp_dir;
    Ok(b)
}

// ---------------------------------------------------------------------------
// Presets to Python values

type Dict<'py> = Bound<'py, PyDict>;

fn dict<'py>(py: Python<'py>, items: &[(&str, Bound<'py, PyAny>)]) -> PyResult<Dict<'py>> {
    let d = PyDict::new(py);
    for (k, v) in items {
        d.set_item(*k, v)?;
    }
    Ok(d)
}

fn any<'py, T: IntoPyObject<'py>>(py: Python<'py>, v: T) -> PyResult<Bound<'py, PyAny>> {
    v.into_bound_py_any(py)
}

fn filter_py<'py>(py: Python<'py>, f: &Filter) -> PyResult<Bound<'py, PyAny>> {
    Ok(dict(
        py,
        &[
            ("require", any(py, f.require)?),
            ("exclude", any(py, f.exclude)?),
            ("min_mapq", any(py, f.min_mapq)?),
            ("unique", any(py, f.unique)?),
        ],
    )?
    .into_any())
}

fn row_key_name(r: &RowKey) -> String {
    match r {
        RowKey::Tag(t) => tag_name(*t),
        RowKey::ReadGroup => "@RG".into(),
        RowKey::File => "@file".into(),
    }
}

fn dedup_name(d: Dedup) -> &'static str {
    match d {
        Dedup::None => "none",
        Dedup::Umi => "umi",
        Dedup::Position => "position",
    }
}

fn flank_py<'py>(py: Python<'py>, kind: &str, f: Flank) -> PyResult<Bound<'py, PyAny>> {
    Ok(dict(
        py,
        &[
            ("kind", any(py, kind)?),
            ("upstream", any(py, f.upstream)?),
            ("downstream", any(py, f.downstream)?),
        ],
    )?
    .into_any())
}

fn contexts_py<'py>(py: Python<'py>, bits: u8) -> PyResult<Bound<'py, PyAny>> {
    let names: Vec<&str> = [(CPG, "CpG"), (CHG, "CHG"), (CHH, "CHH")]
        .into_iter()
        .filter(|(b, _)| bits & b != 0)
        .map(|(_, n)| n)
        .collect();
    Ok(PyTuple::new(py, names)?.into_any())
}

fn coverage_py<'py>(py: Python<'py>, c: &Coverage) -> PyResult<Dict<'py>> {
    let (unit, deletions, max_fragment) = match c.unit {
        CoverageUnit::Bases { deletions } => ("bases", deletions, 0),
        CoverageUnit::Span => ("span", false, 0),
        CoverageUnit::FivePrime => ("five-prime", false, 0),
        CoverageUnit::ThreePrime => ("three-prime", false, 0),
        CoverageUnit::Fragment { max_len } => ("fragment", false, max_len),
    };
    let (strands, strand) = match c.strands {
        StrandSplit::Merged => ("merged", Strandedness::Unstranded),
        StrandSplit::Split(s) => ("split", s),
    };
    let scale = match c.scale {
        Scale::None => any(py, "none")?,
        Scale::Cpm => any(py, "cpm")?,
        Scale::Factor(f) => any(py, f)?,
    };
    dict(
        py,
        &[
            ("name", any(py, &c.name)?),
            ("unit", any(py, unit)?),
            ("deletions", any(py, deletions)?),
            ("max_fragment", any(py, max_fragment)?),
            ("strands", any(py, strands)?),
            ("strand", any(py, strandedness_name(strand))?),
            ("filter", filter_py(py, &c.filter)?),
            ("scale", scale),
            ("sidecar", any(py, c.sidecar)?),
        ],
    )
}

fn labeled_py<'py>(py: Python<'py>, l: &Labeled) -> PyResult<Bound<'py, PyAny>> {
    Ok(dict(
        py,
        &[
            ("conversion", any(py, l.conversion.name())?),
            ("min_quality", any(py, l.conversion.min_quality)?),
            ("min_conversions", any(py, l.min_conversions)?),
        ],
    )?
    .into_any())
}

/// A count matrix's settings; its features are the caller's.
fn counts_py<'py>(py: Python<'py>, c: &Counts) -> PyResult<Dict<'py>> {
    let unit = match c.unit {
        CountUnit::Read => "read",
        CountUnit::Pair => "pair",
        CountUnit::Fragment { tn5: false } => "fragment",
        CountUnit::Fragment { tn5: true } => "tn5-fragment",
        CountUnit::CutSite => "cut-site",
    };
    let features = match c.features {
        FeatureSource::GeneExons => Some("gene-exons"),
        FeatureSource::GeneSpans => Some("gene-spans"),
        _ => None,
    };
    let d = dict(
        py,
        &[
            ("name", any(py, &c.name)?),
            ("unit", any(py, unit)?),
            ("dedup", any(py, dedup_name(c.dedup))?),
            ("rows", any(py, row_key_name(&c.rows))?),
            ("umi_tag", any(py, tag_name(c.umi_tag))?),
            ("strand", any(py, strandedness_name(c.strand))?),
            ("multi_overlap", any(py, c.multi_overlap)?),
            ("filter", filter_py(py, &c.filter)?),
            ("max_fragment", any(py, c.max_fragment)?),
            ("table", any(py, c.table)?),
            ("velocity", any(py, c.velocity)?),
        ],
    )?;
    if let Some(kind) = features {
        d.set_item("features", dict(py, &[("kind", any(py, kind)?)])?)?;
    }
    Ok(d)
}

fn atac_py<'py>(py: Python<'py>, a: &Atac) -> PyResult<Dict<'py>> {
    dict(
        py,
        &[
            ("name", any(py, &a.name)?),
            ("barcode_tag", any(py, tag_name(a.barcode_tag))?),
            ("filter", filter_py(py, &a.filter)?),
            ("max_fragment", any(py, a.max_fragment)?),
            ("fragments", any(py, a.fragments)?),
            ("tiles", any(py, a.tiles)?),
            (
                "gene_activity",
                match a.gene_activity {
                    Some(f) => flank_py(py, "gene-activity", f)?,
                    None => py.None().into_bound(py),
                },
            ),
        ],
    )
}

fn junction_matrix_py<'py>(py: Python<'py>, m: &JunctionMatrix) -> PyResult<Dict<'py>> {
    dict(
        py,
        &[
            ("rows", any(py, row_key_name(&m.rows))?),
            ("umi_tag", any(py, tag_name(m.umi_tag))?),
            ("dedup", any(py, dedup_name(m.dedup))?),
            ("filter", filter_py(py, &m.filter)?),
        ],
    )
}

fn methylation_py<'py>(py: Python<'py>, m: &Methylation) -> PyResult<Dict<'py>> {
    dict(
        py,
        &[
            ("name", any(py, &m.name)?),
            ("source", any(py, m.source.name())?),
            ("contexts", contexts_py(py, m.contexts)?),
            ("merge_cpg", any(py, m.merge_cpg)?),
            ("filter", filter_py(py, &m.filter)?),
            ("min_quality", any(py, m.min_quality)?),
            ("min_probability", any(py, m.min_probability)?),
            ("overlap", any(py, m.overlap)?),
            ("table", any(py, m.table)?),
            ("bigwig", any(py, m.bigwig)?),
        ],
    )
}

fn contacts_py<'py>(py: Python<'py>, c: &Contacts) -> PyResult<Dict<'py>> {
    dict(
        py,
        &[
            ("name", any(py, &c.name)?),
            ("filter", filter_py(py, &c.filter)?),
            ("min_distance", any(py, c.min_distance)?),
            ("self_ligation", any(py, c.self_ligation)?),
            ("dedup", any(py, c.dedup)?),
            ("pairs", any(py, c.pairs)?),
            ("resolutions", PyTuple::new(py, &c.resolutions)?.into_any()),
            ("genome", any(py, &c.genome)?),
            ("barcode_tag", any(py, c.barcode_tag.map(tag_name))?),
            ("sidecar", any(py, c.sidecar)?),
        ],
    )
}

fn alleles_py<'py>(py: Python<'py>, a: &Alleles) -> PyResult<Dict<'py>> {
    dict(
        py,
        &[
            ("name", any(py, &a.name)?),
            ("rows", any(py, row_key_name(&a.rows))?),
            (
                "unit",
                any(py, if a.umi { "molecule" } else { "fragment" })?,
            ),
            ("umi_tag", any(py, tag_name(a.umi_tag))?),
            ("filter", filter_py(py, &a.filter)?),
            ("min_quality", any(py, a.min_quality)?),
            ("table", any(py, a.table)?),
        ],
    )
}

/// The defaults and presets of every product, as the keyword arguments of
/// the package's dataclasses (nested values as dicts; features as dicts
/// with a `kind`). Only what the preset fixes is given: a preset's
/// features are the caller's unless it names gene features.
#[pyfunction]
pub fn quant_defaults(py: Python<'_>) -> PyResult<Bound<'_, PyDict>> {
    let atac = Atac::default();
    let no_sites = PathBuf::new();
    let out = PyDict::new(py);
    out.set_item("filter", filter_py(py, &Filter::new())?)?;
    out.set_item("coverage", coverage_py(py, &Coverage::default())?)?;
    out.set_item(
        "counts",
        counts_py(py, &Counts::new("counts", FeatureSource::Bins(1)))?,
    )?;
    out.set_item(
        "counts.single_cell",
        counts_py(py, &Counts::rna_single_cell())?,
    )?;
    out.set_item("counts.bulk", counts_py(py, &Counts::bulk_rna())?)?;
    out.set_item(
        "counts.atac",
        counts_py(
            py,
            &Counts::atac("atac", FeatureSource::Bins(1), atac.filter),
        )?,
    )?;
    out.set_item(
        "counts.chip",
        counts_py(py, &Counts::chip(FeatureSource::Bins(1)))?,
    )?;
    out.set_item("labeled", labeled_py(py, &Labeled::default())?)?;
    out.set_item("atac", atac_py(py, &atac)?)?;
    out.set_item(
        "junction_matrix.single_cell",
        junction_matrix_py(py, &JunctionMatrix::single_cell())?,
    )?;
    out.set_item(
        "junction_matrix.bulk",
        junction_matrix_py(py, &JunctionMatrix::bulk())?,
    )?;
    let j = Junctions::default();
    out.set_item(
        "junctions",
        dict(
            py,
            &[
                ("name", any(py, &j.name)?),
                ("filter", filter_py(py, &j.filter)?),
                ("pairs", any(py, j.pairs)?),
            ],
        )?,
    )?;
    let l = Labeling::default();
    out.set_item(
        "labeling",
        dict(
            py,
            &[
                ("name", any(py, &l.name)?),
                ("filter", filter_py(py, &l.filter)?),
                ("strand", any(py, strandedness_name(l.strand))?),
                ("conversion", any(py, l.conversion.name())?),
                ("min_quality", any(py, l.conversion.min_quality)?),
            ],
        )?,
    )?;
    let mm = MethylMatrix::new("matrix", FeatureSource::Bins(1));
    out.set_item(
        "methyl_matrix",
        dict(
            py,
            &[
                ("rows", any(py, row_key_name(&mm.rows))?),
                ("contexts", contexts_py(py, mm.contexts)?),
            ],
        )?,
    )?;
    out.set_item("methylation", methylation_py(py, &Methylation::default())?)?;
    out.set_item("contacts", contacts_py(py, &Contacts::default())?)?;
    out.set_item(
        "alleles.single_cell",
        alleles_py(py, &Alleles::single_cell(no_sites.clone()))?,
    )?;
    out.set_item("alleles.bulk", alleles_py(py, &Alleles::bulk(no_sites))?)?;
    Ok(out)
}

// ---------------------------------------------------------------------------
// Plans and runs

/// A task checked against its budget, ready to run (more than once).
#[pyclass(frozen, module = "omicsweave._native")]
pub struct QuantPlan {
    plan: Plan,
    runs: PlanRunState,
}

#[pymethods]
impl QuantPlan {
    /// Plans `task` (a `quant.Task`), reading only the BAM headers and
    /// building the feature spaces; raises when the smallest working set
    /// does not fit the budget.
    #[new]
    #[pyo3(signature = (task, *, threads=None, memory=None, tmp_dir=None))]
    fn new(
        py: Python<'_>,
        task: &Bound<'_, PyAny>,
        threads: Option<usize>,
        memory: Option<u64>,
        tmp_dir: Option<PathBuf>,
    ) -> PyResult<QuantPlan> {
        let t = self::task(task)?;
        let b = budget(threads, memory, tmp_dir)?;
        let p = py.detach(|| plan(t, b)).or_raise()?;
        Ok(QuantPlan {
            plan: p,
            runs: PlanRunState::new(),
        })
    }

    /// (threads, memory bytes, spill directory or None).
    #[getter]
    fn budget(&self) -> (usize, u64, Option<String>) {
        let p = &self.plan;
        (
            p.budget.threads,
            p.budget.memory,
            p.budget.tmp_dir.as_deref().map(path_str),
        )
    }

    /// (fixed, pipeline, buffer) bytes.
    #[getter]
    fn memory(&self) -> (u64, u64, u64) {
        let p = &self.plan;
        (p.memory.fixed, p.memory.pipeline, p.memory.buffer)
    }

    #[getter]
    fn needs(&self) -> Vec<&'static str> {
        self.plan.needs.names()
    }

    #[getter]
    fn stages(&self) -> Vec<(String, u64)> {
        self.plan.stages.clone()
    }

    /// (path, [(reference, length)], read groups, compressed bytes) per input.
    #[getter]
    fn inputs<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyList>> {
        let p = &self.plan;
        let rows = p
            .inputs
            .iter()
            .map(|i| {
                (
                    path_str(&i.path),
                    i.refs.clone(),
                    i.read_groups.clone(),
                    i.len,
                )
            })
            .collect::<Vec<_>>();
        PyList::new(py, rows)
    }

    /// Runs the plan: one pass over each input, then every product is
    /// written. Returns ((kind, path) of every file, records, bytes read,
    /// seconds, seconds of the pass, summary JSON).
    #[pyo3(signature = (progress=None))]
    fn run(&self, py: Python<'_>, progress: Option<Py<PyAny>>) -> PyResult<Quantified> {
        let _run = self.runs.begin()?;
        let report = py
            .detach(|| -> Result<Report> {
                let p = &self.plan;
                match &progress {
                    None => run(p, None),
                    Some(cb) => {
                        // Called on this thread between batches, with the
                        // GIL taken only for the call; an exception is
                        // reported as unraisable and the run goes on.
                        let f = |p: Progress| {
                            Python::attach(|py| {
                                if let Err(e) = cb.call1(py, (p.phase, p.bytes_read, p.bytes_total))
                                {
                                    e.write_unraisable(py, Some(cb.bind(py)));
                                }
                            })
                        };
                        run(p, Some(&f))
                    }
                }
            })
            .or_raise()?;
        Ok((
            report
                .files
                .iter()
                .map(|(k, p)| (*k, path_str(p)))
                .collect(),
            report.records,
            report.bytes_read,
            report.elapsed.as_secs_f64(),
            report.scan.as_secs_f64(),
            report.summary.pretty(),
        ))
    }
}

/// ((kind, path) of every file written, records, bytes read, seconds,
/// seconds of the pass, summary JSON).
type Quantified = (Vec<(&'static str, String)>, u64, u64, f64, f64, String);
