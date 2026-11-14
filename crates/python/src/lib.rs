//! `omicsweave._native`: the extension behind the `omicsweave` package.
//!
//! The Python package defines the public types (keys, results, options and
//! errors); this module takes and returns only primitives and numpy arrays.
//! Each reader or writer keeps its Rust state behind a mutex that is taken
//! only with the GIL released, so decoding runs in parallel across threads
//! and Python never waits on a lock while holding the GIL.

mod consensus;
mod contact;
mod convert;
mod error;
mod index;
mod interval;
mod sequence;
mod shared;
mod track;

use pyo3::prelude::*;

#[pymodule]
fn _native(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add("__version__", env!("CARGO_PKG_VERSION"))?;
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
    m.add_function(wrap_pyfunction!(index::build_fasta_index, m)?)?;
    m.add_function(wrap_pyfunction!(index::build_sidecar, m)?)?;
    m.add_function(wrap_pyfunction!(index::detect, m)?)?;
    m.add_function(wrap_pyfunction!(convert::convert, m)?)?;
    m.add_function(wrap_pyfunction!(convert::plan_convert, m)?)?;
    m.add_function(wrap_pyfunction!(convert::conversions, m)?)?;
    m.add_function(wrap_pyfunction!(convert::parse_size_py, m)?)?;
    m.add_function(wrap_pyfunction!(convert::human_size, m)?)?;
    Ok(())
}
