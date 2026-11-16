//! `omicsweave._native`: the extension behind the `omicsweave` package.
//!
//! The Python package defines the public types (keys, results, options and
//! errors); this module takes and returns only primitives and numpy arrays.
//! Each reader or writer keeps its Rust state behind a mutex that is taken
//! only with the GIL released, so decoding runs in parallel across threads
//! and Python never waits on a lock while holding the GIL.

mod alignment;
mod annotation;
mod cache;
mod consensus;
mod contact;
mod convert;
mod error;
mod index;
mod interval;
mod quant;
mod sequence;
mod shared;
mod track;
mod variant;

use pyo3::prelude::*;
use pyo3::types::PyDict;

#[pymodule]
fn _native(m: &Bound<'_, PyModule>) -> PyResult<()> {
    let py = m.py();
    shared::initialize_process_guard()?;
    let os = py.import("os")?;
    if os.hasattr("register_at_fork")? {
        let hook = wrap_pyfunction!(shared::_after_fork_child, m)?;
        let options = PyDict::new(py);
        options.set_item("after_in_child", hook)?;
        os.getattr("register_at_fork")?.call((), Some(&options))?;
    }
    m.add("__version__", env!("CARGO_PKG_VERSION"))?;
    m.add_class::<cache::PyBlockCache>()?;
    m.add_class::<cache::CacheKey>()?;
    m.add_class::<cache::CacheNamespace>()?;
    m.add_class::<sequence::SequenceReader>()?;
    m.add_class::<sequence::FastaWriter>()?;
    m.add_class::<track::TrackReader>()?;
    m.add_class::<track::BigWigWriter>()?;
    m.add_class::<interval::IntervalReader>()?;
    m.add_class::<interval::BigBedWriter>()?;
    m.add_class::<contact::ContactReader>()?;
    m.add_class::<contact::HicWriter>()?;
    m.add_class::<contact::HicBuilder>()?;
    m.add_class::<consensus::ConsensusBuilder>()?;
    m.add_class::<alignment::AlignmentReader>()?;
    m.add_class::<variant::VariantReader>()?;
    m.add_class::<annotation::GeneAnnotation>()?;
    m.add_class::<annotation::IntervalSet>()?;
    m.add_function(wrap_pyfunction!(index::build_fasta_index, m)?)?;
    m.add_function(wrap_pyfunction!(alignment::build_bam_index, m)?)?;
    m.add_function(wrap_pyfunction!(alignment::build_cram_index, m)?)?;
    m.add_function(wrap_pyfunction!(index::build_sidecar, m)?)?;
    m.add_function(wrap_pyfunction!(index::detect, m)?)?;
    m.add_function(wrap_pyfunction!(convert::convert, m)?)?;
    m.add_function(wrap_pyfunction!(convert::plan_convert, m)?)?;
    m.add_function(wrap_pyfunction!(convert::conversions, m)?)?;
    m.add_function(wrap_pyfunction!(convert::parse_size_py, m)?)?;
    m.add_function(wrap_pyfunction!(convert::human_size, m)?)?;
    m.add_function(wrap_pyfunction!(quant::quantify, m)?)?;
    m.add_function(wrap_pyfunction!(quant::plan_quant, m)?)?;
    Ok(())
}
