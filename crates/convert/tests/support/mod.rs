//! Shared conversion inputs and independent reference-tool adapters.
#![allow(dead_code)]

pub mod assertions;
pub mod rng;
pub mod storage;

use omicsweave_convert::Budget;

/// Budgets that exercise the inline and the threaded paths, roomy and
/// tight.
pub fn budgets() -> Vec<Budget> {
    vec![
        Budget::new(1, 1 << 30),
        Budget::new(6, 1 << 30),
        Budget::new(3, 6 << 20),
    ]
}
