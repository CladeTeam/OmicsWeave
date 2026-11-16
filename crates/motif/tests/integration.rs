//! Default deterministic motif behavior and independent scan/null oracles.
#[path = "support/rng.rs"]
mod rng;
mod support;

#[path = "cases/background.rs"]
mod background;
#[path = "cases/composition.rs"]
mod composition;
#[path = "cases/fimo.rs"]
mod fimo;
#[path = "cases/glam2.rs"]
mod glam2;
#[path = "cases/mast.rs"]
mod mast;
#[path = "cases/mcast.rs"]
mod mcast;
#[path = "cases/null_distribution.rs"]
mod null_distribution;
#[path = "cases/reverse_strand.rs"]
mod reverse_strand;
#[path = "cases/scanning.rs"]
mod scanning;
#[path = "cases/scoring.rs"]
mod scoring;
