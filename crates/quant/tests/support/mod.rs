//! Execution axes shared by the domain targets; expected biological counts
//! stay with each domain's independent oracle.
#![allow(dead_code)]

use std::path::{Path, PathBuf};

use omicsweave_quant::{Budget, MatrixFormat, Task, plan, quantify, run};

use crate::common::{files, read_10x_layer};

pub const STORAGE_FORMATS: [MatrixFormat; 3] =
    [MatrixFormat::Cpk, MatrixFormat::Mtx, MatrixFormat::Csv];

fn assert_bytes(got: &[u8], want: &[u8], context: &str) {
    if got != want {
        let first = got
            .iter()
            .zip(want)
            .position(|(a, b)| a != b)
            .unwrap_or(got.len().min(want.len()));
        panic!(
            "{context}: {} vs {} bytes; first difference at byte {first}",
            got.len(),
            want.len()
        );
    }
}

pub fn assert_product_files(got: &Path, want: &Path, context: &str) {
    let products = |dir: &Path| {
        files(dir)
            .into_iter()
            .filter(|(name, _)| name != "summary.json")
            .collect::<Vec<_>>()
    };
    let (got, want) = (products(got), products(want));
    assert_eq!(
        got.iter().map(|(n, _)| n).collect::<Vec<_>>(),
        want.iter().map(|(n, _)| n).collect::<Vec<_>>(),
        "{context}: delivered paths"
    );
    for ((path, a), (_, b)) in got.iter().zip(&want) {
        assert_bytes(a, b, &format!("{context}: {path}"));
    }
}

/// Separate worker parallelism from forced spilling, so a mismatch names
/// the changed axis. Every domain retains its own seed and spill threshold.
pub fn assert_execution_invariant(
    dir: &Path,
    task: impl Fn(&Path) -> Task,
    stage_bytes: u64,
    minimum_spilling_stages: usize,
) -> PathBuf {
    let baseline = dir.join("serial-roomy");
    let reference = quantify(task(&baseline), Budget::new(1, 1 << 30), None)
        .unwrap_or_else(|e| panic!("serial-roomy: {e}"));
    let expected_summary = summary_without_spills(&baseline);
    for (name, spill) in [("parallel-roomy", false), ("parallel-spill", true)] {
        let out = dir.join(name);
        let mut checked = plan(task(&out), Budget::new(6, 512 << 20))
            .unwrap_or_else(|e| panic!("{name}: planning: {e}"));
        if spill {
            for (_, bytes) in &mut checked.stages {
                *bytes = stage_bytes;
            }
        }
        let report = run(&checked, None).unwrap_or_else(|e| panic!("{name}: execution: {e}"));
        assert_eq!(report.records, reference.records, "{name}: scanned records");
        assert_eq!(
            report.bytes_read, reference.bytes_read,
            "{name}: scanned bytes"
        );
        if spill {
            assert_spills(&out, minimum_spilling_stages, name);
        }
        assert_product_files(&out, &baseline, name);
        assert_eq!(
            summary_without_spills(&out),
            expected_summary,
            "{name}: biological summary fields"
        );
    }
    baseline
}

pub fn assert_spills(dir: &Path, minimum: usize, context: &str) {
    let summary = std::fs::read_to_string(dir.join("summary.json")).unwrap();
    let counters: Vec<u64> = summary
        .lines()
        .filter_map(|line| {
            line.trim().strip_prefix("\"spilled_bytes\": ").map(|n| {
                n.trim_end_matches(',')
                    .parse()
                    .expect("numeric spill counter")
            })
        })
        .collect();
    assert!(
        counters.iter().filter(|&&n| n > 0).count() >= minimum,
        "{context}: need {minimum} spilling stages; counters={counters:?}\n{summary}"
    );
}

fn summary_without_spills(dir: &Path) -> String {
    std::fs::read_to_string(dir.join("summary.json"))
        .unwrap()
        .lines()
        .filter(|line| !line.trim().starts_with("\"spilled_bytes\":"))
        .collect::<Vec<_>>()
        .join("\n")
}

/// Compare decoded layers and axis labels across actual delivery formats.
/// The baseline itself is checked against a domain oracle in its caller.
pub fn assert_storage_invariant(
    baseline: &Path,
    task: impl Fn(&Path) -> Task,
    layers: &[(&str, &str)],
) {
    let original_format = task(baseline).matrix;
    for format in STORAGE_FORMATS {
        if format == original_format {
            continue;
        }
        let name = baseline.file_name().unwrap().to_string_lossy();
        let out = baseline.with_file_name(format!("{name}-{format:?}"));
        let mut alternate = task(&out);
        alternate.matrix = format;
        quantify(alternate, Budget::new(3, 512 << 20), None)
            .unwrap_or_else(|e| panic!("storage {format:?}: {e}"));
        for &(product, layer) in layers {
            let (got, want) = (out.join(product), baseline.join(product));
            assert_eq!(
                read_10x_layer(&got, layer),
                read_10x_layer(&want, layer),
                "storage {format:?}: {product}/{layer}"
            );
            for labels in ["barcodes.tsv.gz", "features.tsv.gz"] {
                assert_bytes(
                    &std::fs::read(got.join(labels)).unwrap(),
                    &std::fs::read(want.join(labels)).unwrap(),
                    &format!("storage {format:?}: {product}/{labels}"),
                );
            }
        }
    }
}
