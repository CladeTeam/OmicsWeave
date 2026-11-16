//! Fixtures, independent oracles and output readers shared by quant tests.
//! Also imported by the reference-window unit tests.
#![allow(dead_code, unused_imports)]

mod alignment;
mod fixtures;
mod oracle;
mod output;
mod support;

pub use alignment::*;
pub use fixtures::*;
pub use oracle::*;
pub use output::*;
pub use support::*;
