//! Public contracts shared by all format readers.

mod support;

#[path = "cases/bins.rs"]
mod bins;
#[path = "cases/cache.rs"]
mod cache;
#[path = "cases/coordinates.rs"]
mod coordinates;
#[path = "cases/layouts.rs"]
mod layouts;
#[path = "cases/orientation.rs"]
mod orientation;
#[path = "cases/sidecar/mod.rs"]
mod sidecar;
#[path = "cases/sources.rs"]
mod sources;

#[path = "cases/alignment_layout.rs"]
mod alignment_layout;
