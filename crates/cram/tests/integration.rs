//! CRAM behavior against independent htslib fixtures and record expectations.

mod support;

#[path = "cases/cache.rs"]
mod cache;
#[path = "cases/format.rs"]
mod format;
#[path = "cases/index.rs"]
mod index;
#[path = "cases/modifications.rs"]
mod modifications;
#[path = "cases/queries.rs"]
mod queries;
#[path = "cases/reference.rs"]
mod reference;
