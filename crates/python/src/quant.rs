//! Quantification: product descriptions in (dicts of primitives), the task
//! built, planned and run by `omicsweave_quant` with the GIL released.

use std::path::PathBuf;

use omicsweave_core::Result;
use omicsweave_gtf::{Aliases, Flank};
use omicsweave_quant::{
    AlleleGenes, Alleles, Atac, Budget, CallSource, Contacts, Conversion, CountUnit, Counts,
    Coverage, CoverageUnit, Dedup, FeatureSource, Filter, JunctionMatrix, Junctions, Labeled,
    Labeling, MatrixFormat, MethylMatrix, Methylation, Plan, Product, Progress, Report, RowKey,
    Scale, StrandSplit, Strandedness, Task, contexts_from_names, parse_size, plan, run,
};
use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;
use pyo3::types::PyDict;

use crate::annotation::GeneAnnotation;
use crate::error::OrRaise;
use crate::shared::path_str;

fn get<'py, T: FromPyObjectOwned<'py>>(d: &Bound<'py, PyDict>, key: &str) -> PyResult<Option<T>> {
    match d.get_item(key)? {
        Some(v) if !v.is_none() => Ok(Some(v.extract::<T>().map_err(Into::into)?)),
        _ => Ok(None),
    }
}

fn need<'py, T: FromPyObjectOwned<'py>>(d: &Bound<'py, PyDict>, key: &str) -> PyResult<T> {
    get(d, key)?.ok_or_else(|| PyValueError::new_err(format!("product option {key:?} is missing")))
}

fn tag(s: &str) -> PyResult<[u8; 2]> {
    let b = s.as_bytes();
    if b.len() != 2 {
        return Err(PyValueError::new_err(format!(
            "a tag has two characters, not {s:?}"
        )));
    }
    Ok([b[0], b[1]])
}

fn library(d: &Bound<'_, PyDict>) -> PyResult<Strandedness> {
    let code: u8 = get(d, "library")?.unwrap_or(0);
    Strandedness::from_code(code)
        .ok_or_else(|| PyValueError::new_err(format!("library must be 0, 1 or 2, not {code}")))
}

fn filter(d: &Bound<'_, PyDict>) -> PyResult<Filter> {
    filter_at(d, "filter")
}

fn filter_at(d: &Bound<'_, PyDict>, key: &str) -> PyResult<Filter> {
    let (require, exclude, min_mapq, unique): (u16, u16, u8, bool) = need(d, key)?;
    Ok(Filter {
        require,
        exclude,
        min_mapq,
        unique,
    })
}

fn coverage(d: &Bound<'_, PyDict>) -> PyResult<Coverage> {
    let max_len: u32 = get(d, "max_fragment")?.unwrap_or(0);
    let unit = match need::<String>(d, "unit")?.as_str() {
        "bases" => CoverageUnit::Bases { deletions: false },
        "bases+deletions" => CoverageUnit::Bases { deletions: true },
        "span" => CoverageUnit::Span,
        "five-prime" => CoverageUnit::FivePrime,
        "three-prime" => CoverageUnit::ThreePrime,
        "fragment" => CoverageUnit::Fragment { max_len },
        u => {
            return Err(PyValueError::new_err(format!(
                "unknown coverage unit {u:?}"
            )));
        }
    };
    let strands = match need::<String>(d, "strands")?.as_str() {
        "merged" => StrandSplit::Merged,
        "split" => StrandSplit::Split(library(d)?),
        s => {
            return Err(PyValueError::new_err(format!(
                "strands must be merged or split, not {s:?}"
            )));
        }
    };
    let scale = match d.get_item("scale")? {
        None => Scale::None,
        Some(v) if v.is_none() => Scale::None,
        Some(v) => match v.extract::<String>() {
            Ok(s) if s == "none" => Scale::None,
            Ok(s) if s == "cpm" => Scale::Cpm,
            Ok(s) => {
                return Err(PyValueError::new_err(format!(
                    "scale must be none, cpm or a number, not {s:?}"
                )));
            }
            Err(_) => Scale::Factor(v.extract::<f64>()?),
        },
    };
    Ok(Coverage {
        name: need(d, "name")?,
        unit,
        strands,
        filter: filter(d)?,
        scale,
        sidecar: get(d, "sidecar")?.unwrap_or(false),
    })
}

fn flank(d: &Bound<'_, PyDict>) -> PyResult<Flank> {
    Ok(Flank::new(
        get(d, "upstream")?.unwrap_or(0),
        get(d, "downstream")?.unwrap_or(0),
    ))
}

fn features(d: &Bound<'_, PyDict>) -> PyResult<FeatureSource> {
    Ok(match need::<String>(d, "features")?.as_str() {
        "genes" => FeatureSource::GeneExons,
        "gene-full" => FeatureSource::GeneSpans,
        "promoters" => FeatureSource::Promoters(flank(d)?),
        "gene-activity" => FeatureSource::GeneActivity(flank(d)?),
        "peaks" => FeatureSource::BigBed(need::<PathBuf>(d, "peaks")?),
        "bins" => FeatureSource::Bins(need(d, "bin_size")?),
        f => return Err(PyValueError::new_err(format!("unknown features {f:?}"))),
    })
}

fn counts(d: &Bound<'_, PyDict>) -> PyResult<Counts> {
    let features = features(d)?;
    let (unit, molecule) = match need::<String>(d, "unit")?.as_str() {
        "read" => (CountUnit::Read, false),
        "molecule" => (CountUnit::Read, true),
        "pair" => (CountUnit::Pair, false),
        "fragment" => (CountUnit::Fragment { tn5: false }, false),
        "tn5-fragment" => (CountUnit::Fragment { tn5: true }, false),
        "cut-site" => (CountUnit::CutSite, false),
        u => return Err(PyValueError::new_err(format!("unknown count unit {u:?}"))),
    };
    let rows = row_key(&need::<String>(d, "rows")?)?;
    let dedup = match get::<String>(d, "dedup")?.as_deref() {
        None | Some("none") if !molecule => Dedup::None,
        None | Some("umi") => Dedup::Umi,
        Some("position") => Dedup::Position,
        Some("none") => {
            return Err(PyValueError::new_err("molecules are counted by UMI"));
        }
        Some(x) => return Err(PyValueError::new_err(format!("unknown dedup {x:?}"))),
    };
    Ok(Counts {
        name: need(d, "name")?,
        features,
        biotypes: get(d, "biotypes")?,
        unit,
        rows,
        whitelist: get(d, "whitelist")?,
        umi_tag: tag(&get::<String>(d, "umi_tag")?.unwrap_or_else(|| "UB".into()))?,
        dedup,
        strand: library(d)?,
        multi_overlap: get(d, "multi_overlap")?.unwrap_or(false),
        filter: filter(d)?,
        max_fragment: get(d, "max_fragment")?.unwrap_or(0),
        table: get(d, "table")?.unwrap_or(false),
        velocity: get(d, "velocity")?.unwrap_or(false),
        labeled: match get::<(String, u8, u32)>(d, "labeled")? {
            Some((c, q, n)) => Some(Labeled {
                conversion: Conversion::parse(&c, q).or_raise()?,
                min_conversions: n,
            }),
            None => None,
        },
    })
}

fn row_key(s: &str) -> PyResult<RowKey> {
    Ok(match s {
        "@RG" => RowKey::ReadGroup,
        "@file" => RowKey::File,
        t => RowKey::Tag(tag(t)?),
    })
}

fn junctions(d: &Bound<'_, PyDict>) -> PyResult<Junctions> {
    let matrix = match get::<String>(d, "matrix_rows")? {
        None => None,
        Some(rows) => Some(JunctionMatrix {
            rows: row_key(&rows)?,
            whitelist: get(d, "matrix_whitelist")?,
            umi_tag: tag(&get::<String>(d, "matrix_umi_tag")?.unwrap_or_else(|| "UB".into()))?,
            dedup: match get::<String>(d, "matrix_dedup")?.as_deref() {
                None | Some("none") => Dedup::None,
                Some("umi") => Dedup::Umi,
                Some(x) => {
                    return Err(PyValueError::new_err(format!(
                        "junction matrices deduplicate by none or umi, not {x:?}"
                    )));
                }
            },
            filter: filter_at(d, "matrix_filter")?,
        }),
    };
    Ok(Junctions {
        name: need(d, "name")?,
        filter: filter(d)?,
        pairs: get(d, "pairs")?.unwrap_or(true),
        matrix,
    })
}

fn labeling(d: &Bound<'_, PyDict>) -> PyResult<Labeling> {
    Ok(Labeling {
        name: need(d, "name")?,
        filter: filter(d)?,
        strand: library(d)?,
        conversion: Conversion::parse(
            &get::<String>(d, "conversion")?.unwrap_or_else(|| "T>C".into()),
            get(d, "min_quality")?.unwrap_or(20),
        )
        .or_raise()?,
    })
}

fn atac(d: &Bound<'_, PyDict>) -> PyResult<Atac> {
    let activity: Option<(u32, u32)> = get(d, "gene_activity")?;
    Ok(Atac {
        name: need(d, "name")?,
        barcode_tag: tag(&need::<String>(d, "barcode_tag")?)?,
        whitelist: get(d, "whitelist")?,
        filter: filter(d)?,
        max_fragment: get(d, "max_fragment")?.unwrap_or(0),
        fragments: get(d, "fragments")?.unwrap_or(true),
        peaks: get(d, "peaks")?,
        tiles: get(d, "tiles")?,
        gene_activity: activity.map(|(u, w)| Flank::new(u, w)),
    })
}

fn dicts<'py>(d: &Bound<'py, PyDict>, key: &str) -> PyResult<Vec<Bound<'py, PyDict>>> {
    match d.get_item(key)? {
        Some(v) if !v.is_none() => Ok(v.extract::<Vec<Bound<'py, PyDict>>>()?),
        _ => Ok(Vec::new()),
    }
}

fn contexts(d: &Bound<'_, PyDict>) -> PyResult<u8> {
    let names: Vec<String> = get(d, "contexts")?.unwrap_or_else(|| vec!["CpG".into()]);
    contexts_from_names(&names).or_raise()
}

fn methylation(d: &Bound<'_, PyDict>) -> PyResult<Methylation> {
    let mut matrices = Vec::new();
    for m in dicts(d, "matrices")? {
        matrices.push(MethylMatrix {
            name: need(&m, "name")?,
            features: features(&m)?,
            biotypes: get(&m, "biotypes")?,
            rows: row_key(&need::<String>(&m, "rows")?)?,
            whitelist: get(&m, "whitelist")?,
            contexts: contexts(&m)?,
        });
    }
    Ok(Methylation {
        name: need(d, "name")?,
        source: CallSource::from_name(&need::<String>(d, "source")?).or_raise()?,
        contexts: contexts(d)?,
        merge_cpg: get(d, "merge_cpg")?.unwrap_or(false),
        filter: filter(d)?,
        min_quality: get(d, "min_quality")?.unwrap_or(5),
        min_probability: get(d, "min_probability")?.unwrap_or(0.5),
        overlap: get(d, "overlap")?.unwrap_or(true),
        table: get(d, "table")?.unwrap_or(true),
        bigwig: get(d, "bigwig")?.unwrap_or(true),
        matrices,
    })
}

fn contacts(d: &Bound<'_, PyDict>) -> PyResult<Contacts> {
    Ok(Contacts {
        name: need(d, "name")?,
        filter: filter(d)?,
        min_distance: get(d, "min_distance")?.unwrap_or(0),
        self_ligation: get(d, "self_ligation")?.unwrap_or(0),
        dedup: get(d, "dedup")?.unwrap_or(true),
        pairs: get(d, "pairs")?.unwrap_or(true),
        resolutions: get(d, "resolutions")?.unwrap_or_default(),
        genome: get(d, "genome")?.unwrap_or_else(|| "unknown".into()),
        barcode_tag: match get::<String>(d, "barcode_tag")? {
            Some(t) => Some(tag(&t)?),
            None => None,
        },
        whitelist: get(d, "whitelist")?,
        sidecar: get(d, "sidecar")?.unwrap_or(false),
    })
}

fn alleles(d: &Bound<'_, PyDict>) -> PyResult<Alleles> {
    let umi = match need::<String>(d, "unit")?.as_str() {
        "molecule" => true,
        "fragment" => false,
        u => {
            return Err(PyValueError::new_err(format!(
                "allele units are molecule or fragment, not {u:?}"
            )));
        }
    };
    let genes = match d.get_item("genes")? {
        Some(v) if !v.is_none() => {
            let (sample, lib): (Option<String>, u8) = v.extract()?;
            Some(AlleleGenes {
                sample,
                strand: Strandedness::from_code(lib).ok_or_else(|| {
                    PyValueError::new_err(format!("library must be 0, 1 or 2, not {lib}"))
                })?,
            })
        }
        _ => None,
    };
    Ok(Alleles {
        name: need(d, "name")?,
        sites: need(d, "sites")?,
        rows: row_key(&need::<String>(d, "rows")?)?,
        whitelist: get(d, "whitelist")?,
        umi,
        umi_tag: tag(&get::<String>(d, "umi_tag")?.unwrap_or_else(|| "UB".into()))?,
        filter: filter(d)?,
        min_quality: get(d, "min_quality")?.unwrap_or(20),
        table: get(d, "table")?.unwrap_or(true),
        genes,
    })
}

fn product(d: &Bound<'_, PyDict>) -> PyResult<Product> {
    Ok(match need::<String>(d, "kind")?.as_str() {
        "coverage" => Product::Coverage(coverage(d)?),
        "counts" => Product::Counts(counts(d)?),
        "atac" => Product::Atac(atac(d)?),
        "junctions" => Product::Junctions(junctions(d)?),
        "labeling" => Product::Labeling(labeling(d)?),
        "methylation" => Product::Methylation(methylation(d)?),
        "hic" => Product::Contacts(contacts(d)?),
        "alleles" => Product::Alleles(alleles(d)?),
        k => return Err(PyValueError::new_err(format!("unknown product kind {k:?}"))),
    })
}

#[allow(clippy::too_many_arguments)]
fn task(
    inputs: Vec<PathBuf>,
    out_dir: PathBuf,
    products: Vec<Bound<'_, PyDict>>,
    annotation: Option<PyRef<'_, GeneAnnotation>>,
    aliases: Option<Vec<(String, String)>>,
    matrix: &str,
    reference: Option<PathBuf>,
    variants: Option<PathBuf>,
) -> PyResult<Task> {
    let mut t = Task::new(inputs, out_dir);
    t.annotation = annotation.map(|a| a.shared()).transpose()?;
    t.aliases = aliases.map(Aliases::from_pairs);
    t.matrix = MatrixFormat::from_name(matrix).or_raise()?;
    t.reference = reference;
    t.variants = variants;
    for p in &products {
        t.products.push(product(p)?);
    }
    Ok(t)
}

fn budget(threads: Option<usize>, memory: &str, tmp_dir: Option<PathBuf>) -> Result<Budget> {
    let mut b = Budget::default();
    if let Some(t) = threads {
        b.threads = t;
    }
    b.memory = parse_size(memory)?;
    b.tmp_dir = tmp_dir;
    Ok(b)
}

/// (threads, memory budget, fixed, pipeline, buffer, needs, [(stage, buffer)]).
type Planned = (
    usize,
    u64,
    u64,
    u64,
    u64,
    Vec<&'static str>,
    Vec<(String, u64)>,
);

fn planned(p: &Plan) -> Planned {
    (
        p.budget.threads,
        p.budget.memory,
        p.memory.fixed,
        p.memory.pipeline,
        p.memory.buffer,
        p.needs.names(),
        p.stages.clone(),
    )
}

/// Checks a quantification against its budget, reading only the BAM
/// headers (and building the feature spaces).
#[pyfunction]
#[pyo3(signature = (inputs, out_dir, products, *, annotation=None, aliases=None, matrix="cpk",
                    reference=None, variants=None, threads=None, memory="16G"))]
#[allow(clippy::too_many_arguments)]
pub fn plan_quant(
    py: Python<'_>,
    inputs: Vec<PathBuf>,
    out_dir: PathBuf,
    products: Vec<Bound<'_, PyDict>>,
    annotation: Option<PyRef<'_, GeneAnnotation>>,
    aliases: Option<Vec<(String, String)>>,
    matrix: &str,
    reference: Option<PathBuf>,
    variants: Option<PathBuf>,
    threads: Option<usize>,
    memory: &str,
) -> PyResult<Planned> {
    let t = task(
        inputs, out_dir, products, annotation, aliases, matrix, reference, variants,
    )?;
    let memory = memory.to_owned();
    py.detach(|| -> Result<Planned> {
        let p = plan(t, budget(threads, &memory, None)?)?;
        Ok(planned(&p))
    })
    .or_raise()
}

/// ((kind, path) of every file written, records, bytes read, seconds,
/// summary JSON).
type Quantified = (Vec<(&'static str, String)>, u64, u64, f64, String);

/// Runs a quantification: one pass over the inputs, every product written.
#[pyfunction]
#[pyo3(signature = (inputs, out_dir, products, *, annotation=None, aliases=None, matrix="cpk",
                    reference=None, variants=None, threads=None, memory="16G", tmp_dir=None,
                    progress=None))]
#[allow(clippy::too_many_arguments)]
pub fn quantify(
    py: Python<'_>,
    inputs: Vec<PathBuf>,
    out_dir: PathBuf,
    products: Vec<Bound<'_, PyDict>>,
    annotation: Option<PyRef<'_, GeneAnnotation>>,
    aliases: Option<Vec<(String, String)>>,
    matrix: &str,
    reference: Option<PathBuf>,
    variants: Option<PathBuf>,
    threads: Option<usize>,
    memory: &str,
    tmp_dir: Option<PathBuf>,
    progress: Option<Py<PyAny>>,
) -> PyResult<Quantified> {
    let t = task(
        inputs, out_dir, products, annotation, aliases, matrix, reference, variants,
    )?;
    let memory = memory.to_owned();
    let report = py
        .detach(|| -> Result<Report> {
            let p = plan(t, budget(threads, &memory, tmp_dir)?)?;
            match &progress {
                None => run(&p, None),
                Some(cb) => {
                    // Called on this thread between batches, with the GIL
                    // taken only for the call; an exception is reported as
                    // unraisable and the run goes on.
                    let f = |p: Progress| {
                        Python::attach(|py| {
                            if let Err(e) = cb.call1(py, (p.phase, p.bytes_read, p.bytes_total)) {
                                e.write_unraisable(py, Some(cb.bind(py)));
                            }
                        })
                    };
                    run(&p, Some(&f))
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
        report.summary.pretty(),
    ))
}
