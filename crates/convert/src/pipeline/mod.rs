//! The stages every conversion is built from: reading and inflating the
//! source, sorting out of core, encoding text lines into BGZF records, and
//! publishing the destination files.

pub(crate) mod hic;
pub(crate) mod input;
pub(crate) mod output;
pub(crate) mod records;
pub(crate) mod sort;
