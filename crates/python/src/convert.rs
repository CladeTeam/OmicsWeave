//! Format conversion: arguments in, `omicsweave_convert` does the rest with
//! the GIL released.

use std::path::PathBuf;

use omicsweave_convert::{
    Budget, Conversion, Format, Index, Options, Progress, Report, convert as run, human,
    parse_size, plan,
};
use omicsweave_core::Result;
use pyo3::prelude::*;

use crate::error::OrRaise;
use crate::shared::path_str;

/// The conversion of `src` into `dst`, and the source's format. A format
/// not named is told by content for the source and by name for the
/// destination.
fn choose(
    src: PathBuf,
    dst: PathBuf,
    src_format: Option<&str>,
    dst_format: Option<&str>,
) -> Result<(Conversion, Format)> {
    let from = match src_format {
        Some(name) => Format::from_name(name)?,
        None => Format::detect(&src)?,
    };
    let to = match dst_format {
        Some(name) => Format::from_name(name)?,
        None => Format::from_extension(&dst)?,
    };
    Ok((Conversion::between(from, to, src, dst)?, from))
}

fn budget(threads: Option<usize>, memory: &str, tmp_dir: Option<PathBuf>) -> Result<Budget> {
    let mut budget = Budget::default();
    if let Some(t) = threads {
        budget.threads = t;
    }
    budget.memory = parse_size(memory)?;
    budget.tmp_dir = tmp_dir;
    Ok(budget)
}

/// (source format, (format, path) of every file written, records, bytes
/// read, bytes spilled, seconds).
type Converted = (
    &'static str,
    Vec<(&'static str, String)>,
    u64,
    u64,
    u64,
    f64,
);

/// Converts `src` into `dst`; returns (source format, (format, path) of
/// every file written with the destination first, records, bytes read,
/// bytes spilled, seconds).
#[pyfunction]
#[pyo3(signature = (src, dst, *, src_format=None, dst_format=None, threads=None, memory="4G",
                    tmp_dir=None, level=None, index=None, chrom_sizes=None, resolutions=None,
                    sidecar=true, add_contigs=false, progress=None))]
#[allow(clippy::too_many_arguments)]
pub fn convert(
    py: Python<'_>,
    src: PathBuf,
    dst: PathBuf,
    src_format: Option<String>,
    dst_format: Option<String>,
    threads: Option<usize>,
    memory: &str,
    tmp_dir: Option<PathBuf>,
    level: Option<i32>,
    index: Option<String>,
    chrom_sizes: Option<PathBuf>,
    resolutions: Option<Vec<u32>>,
    sidecar: bool,
    add_contigs: bool,
    progress: Option<Py<PyAny>>,
) -> PyResult<Converted> {
    let memory = memory.to_owned();
    let (from, report) = py
        .detach(|| -> Result<(Format, Report)> {
            let (mut conversion, from) =
                choose(src, dst, src_format.as_deref(), dst_format.as_deref())?;
            conversion.configure(&Options {
                level,
                index: index.as_deref().map(Index::from_name).transpose()?,
                chrom_sizes,
                resolutions,
                sidecar,
                add_contigs,
            })?;
            let budget = budget(threads, &memory, tmp_dir)?;
            let report = match &progress {
                None => run(conversion, budget, None)?,
                Some(cb) => {
                    // Called on this thread between batches; the GIL is
                    // taken only for the call. An exception from the
                    // callback is reported as unraisable and the conversion
                    // goes on.
                    let f = |p: Progress| {
                        Python::attach(|py| {
                            if let Err(e) = cb.call1(py, (p.phase, p.bytes_read, p.bytes_total)) {
                                e.write_unraisable(py, Some(cb.bind(py)));
                            }
                        })
                    };
                    run(conversion, budget, Some(&f))?
                }
            };
            Ok((from, report))
        })
        .or_raise()?;
    Ok((
        from.name(),
        report
            .files
            .iter()
            .map(|(f, p)| (f.name(), path_str(p)))
            .collect(),
        report.records,
        report.bytes_read,
        report.spilled,
        report.elapsed.as_secs_f64(),
    ))
}

/// (source format, destination format, threads, memory budget, fixed,
/// pipeline, buffer bytes).
type Planned = (&'static str, &'static str, usize, u64, u64, u64, u64);

/// Checks the conversion of `src` into `dst` against its budget without
/// reading any data; returns (source format, destination format, threads,
/// memory budget, fixed, pipeline, buffer bytes).
#[pyfunction]
#[pyo3(signature = (src, dst, *, src_format=None, dst_format=None, threads=None, memory="4G"))]
pub fn plan_convert(
    py: Python<'_>,
    src: PathBuf,
    dst: PathBuf,
    src_format: Option<String>,
    dst_format: Option<String>,
    threads: Option<usize>,
    memory: &str,
) -> PyResult<Planned> {
    let memory = memory.to_owned();
    py.detach(|| -> Result<Planned> {
        let (conversion, from) = choose(src, dst, src_format.as_deref(), dst_format.as_deref())?;
        let to = conversion.target();
        let p = plan(conversion, budget(threads, &memory, None)?)?;
        Ok((
            from.name(),
            to.name(),
            p.budget.threads,
            p.budget.memory,
            p.memory.fixed,
            p.memory.pipeline,
            p.memory.buffer,
        ))
    })
    .or_raise()
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
