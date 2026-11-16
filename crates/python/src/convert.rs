//! Format conversion: a conversion is configured and planned once
//! ([`ConvertPlan`]), and the plan is what runs, with the GIL released.

use std::path::PathBuf;

use omicsweave_bbi::WriteOptions;
use omicsweave_convert::{
    Conversion, Format, Index, Options, Plan, Progress, Report, detect_sources, human, parse_size,
    plan, run,
};
use omicsweave_core::Result;
use pyo3::prelude::*;

use crate::error::OrRaise;
use crate::quant::budget;
use crate::shared::{PlanRunState, path_str};
use crate::track::write_options;

/// The conversion of `src` (and `more_src`, for conversions that read
/// several sources) into `dst`, and the source's format. A format not named
/// is told by content for the sources (every file, which must agree)
/// and by name for the destination.
fn choose(
    src: PathBuf,
    more_src: Vec<PathBuf>,
    dst: PathBuf,
    src_format: Option<&str>,
    dst_format: Option<&str>,
) -> Result<(Conversion, Format)> {
    let from = match src_format {
        Some(name) => Format::from_name(name)?,
        None if more_src.is_empty() && src.exists() => Format::detect(&src)?,
        None => {
            let all: Vec<PathBuf> = std::iter::once(src.clone())
                .chain(more_src.iter().cloned())
                .collect();
            detect_sources(&all)?
        }
    };
    let to = match dst_format {
        Some(name) => Format::from_name(name)?,
        None => Format::from_extension(&dst)?,
    };
    let mut conversion = Conversion::between(from, to, src, dst)?;
    conversion.add_sources(more_src)?;
    Ok((conversion, from))
}

/// `(items_per_slot, block_size, compress, level, zoom_levels, zoom_reductions)`.
type BbiOptions = (Option<u32>, u32, bool, i32, usize, Option<Vec<u32>>);

fn bbi_options(o: BbiOptions) -> WriteOptions {
    let (items_per_slot, block_size, compress, level, zoom_levels, zoom_reductions) = o;
    write_options(
        items_per_slot,
        block_size,
        compress,
        level,
        zoom_levels,
        zoom_reductions,
    )
}

/// A conversion configured and checked against its budget, ready to run.
#[pyclass(frozen, module = "omicsweave._native")]
pub struct ConvertPlan {
    plan: Plan,
    runs: PlanRunState,
    source: Format,
}

/// (format, path) of every file written with the destination first,
/// records, bytes read, bytes spilled, seconds, warnings.
type Converted = (Vec<(&'static str, String)>, u64, u64, u64, f64, Vec<String>);

#[pymethods]
impl ConvertPlan {
    /// Configures the conversion of `src` into `dst` (an option it does not
    /// have is an error) and checks it against the budget without reading
    /// any data.
    #[new]
    #[pyo3(signature = (src, dst, *, more_src=Vec::new(), src_format=None, dst_format=None, threads=None,
                        memory=None, tmp_dir=None, level=None, index=None, min_shift=None,
                        chrom_sizes=None, resolutions=None, value=None, r#mod=None,
                        strand=None, ref_species=None, sidecar=false, add_contigs=false, sort=false,
                        reference=None, verify_checksums=true, bbi=None))]
    #[allow(clippy::too_many_arguments)]
    fn new(
        py: Python<'_>,
        src: PathBuf,
        dst: PathBuf,
        more_src: Vec<PathBuf>,
        src_format: Option<String>,
        dst_format: Option<String>,
        threads: Option<usize>,
        memory: Option<u64>,
        tmp_dir: Option<PathBuf>,
        level: Option<i32>,
        index: Option<String>,
        min_shift: Option<u32>,
        chrom_sizes: Option<PathBuf>,
        resolutions: Option<Vec<u32>>,
        value: Option<String>,
        r#mod: Option<String>,
        strand: Option<String>,
        ref_species: Option<String>,
        sidecar: bool,
        add_contigs: bool,
        sort: bool,
        reference: Option<PathBuf>,
        verify_checksums: bool,
        bbi: Option<BbiOptions>,
    ) -> PyResult<ConvertPlan> {
        let budget = budget(threads, memory, tmp_dir)?;
        let (p, source) = py
            .detach(|| -> Result<(Plan, Format)> {
                let (mut conversion, source) = choose(
                    src,
                    more_src,
                    dst,
                    src_format.as_deref(),
                    dst_format.as_deref(),
                )?;
                conversion.configure(&Options {
                    level,
                    index: index.as_deref().map(Index::from_name).transpose()?,
                    chrom_sizes,
                    resolutions,
                    value,
                    mod_code: r#mod,
                    strand,
                    ref_species,
                    sidecar,
                    add_contigs,
                    sort,
                    reference,
                    min_shift,
                    verify_checksums,
                    bbi: bbi.map(bbi_options),
                })?;
                Ok((plan(conversion, budget)?, source))
            })
            .or_raise()?;
        Ok(ConvertPlan {
            plan: p,
            runs: PlanRunState::new(),
            source,
        })
    }

    /// (source format, target format).
    #[getter]
    fn formats(&self) -> (&'static str, &'static str) {
        let p = &self.plan;
        (self.source.name(), p.conversion.target().name())
    }

    /// (source path, target path).
    #[getter]
    fn paths(&self) -> (String, String) {
        let p = &self.plan;
        (path_str(p.conversion.src()), path_str(p.conversion.dst()))
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

    /// Runs the conversion; the target and its companions appear only once
    /// the whole conversion succeeded.
    #[pyo3(signature = (progress=None))]
    fn run(&self, py: Python<'_>, progress: Option<Py<PyAny>>) -> PyResult<Converted> {
        let _run = self.runs.begin()?;
        let report = py
            .detach(|| -> Result<Report> {
                let p = &self.plan;
                match &progress {
                    None => run(p, None),
                    Some(cb) => {
                        // Called on this thread between batches; the GIL is
                        // taken only for the call. An exception from the
                        // callback is reported as unraisable and the
                        // conversion goes on.
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
                .map(|(f, p)| (f.name(), path_str(p)))
                .collect(),
            report.records,
            report.bytes_read,
            report.spilled,
            report.elapsed.as_secs_f64(),
            report.warnings,
        ))
    }
}

/// Every (source, destination) pair of format names there is a conversion for.
#[pyfunction]
pub fn conversions() -> Vec<(&'static str, &'static str)> {
    Conversion::routes()
        .map(|(from, to)| (from.name(), to.name()))
        .collect()
}

/// `text` (`512M`, `8G`, `1.5GiB`, a byte count) in bytes.
#[pyfunction(name = "parse_size")]
pub fn parse_size_py(text: &str) -> PyResult<u64> {
    parse_size(text).or_raise()
}

/// `bytes` for people: `512 B`, `1.5 GiB`.
#[pyfunction(name = "human_size")]
pub fn human_size(bytes: u64) -> String {
    human(bytes)
}
