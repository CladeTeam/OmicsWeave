//! Public BAM behavior checked against independent htslib decoding and indexes.
mod support;

#[path = "cases/cache.rs"]
mod cache;
#[path = "cases/column_capacity.rs"]
mod column_capacity;
#[path = "cases/columns.rs"]
mod columns;
#[path = "cases/errors.rs"]
mod errors;
#[path = "cases/header.rs"]
mod header;
#[path = "cases/index.rs"]
mod index;
#[path = "cases/modifications/errors.rs"]
mod modification_errors;
#[path = "cases/modifications/records.rs"]
mod modification_records;
#[path = "cases/modifications/regions.rs"]
mod modification_regions;
#[path = "cases/queries.rs"]
mod queries;
#[path = "cases/records.rs"]
mod records;
