#![allow(dead_code)]
pub mod fixtures;
pub mod io;
pub mod oracle;
pub mod straw;
pub mod wire;

pub use fixtures::{Rng, Spec, random_spec};
pub use straw::straw;
pub use wire::write;
