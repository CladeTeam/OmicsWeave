//! Checked plans, repeatable execution and public progress reporting.

use std::sync::Mutex;

use omicsweave_convert::{Budget, Conversion, FastaToBgzipFasta, Format, plan, run};

use crate::support::storage::inflate;

#[test]
fn planning_rejects_invalid_budgets_and_same_paths_without_changing_files() {
    let dir = tempfile::tempdir().unwrap();
    let src = dir.path().join("input.fa");
    let dst = dir.path().join("output.fa.gz");
    let input = b">s\nACGT\n";
    std::fs::write(&src, input).unwrap();
    std::fs::write(&dst, b"previous output").unwrap();
    let mut zero_threads = Budget::new(1, 32 << 20);
    zero_threads.threads = 0;
    for budget in [zero_threads, Budget::new(1, 0)] {
        let conversion = Conversion::FastaToBgzipFasta(FastaToBgzipFasta::new(&src, &dst));
        assert!(plan(conversion, budget.clone()).is_err(), "{budget:?}");
        assert_eq!(
            std::fs::read(&dst).unwrap(),
            b"previous output",
            "{budget:?}"
        );
        assert_eq!(std::fs::read_dir(dir.path()).unwrap().count(), 2);
    }
    let conversion = Conversion::FastaToBgzipFasta(FastaToBgzipFasta::new(&src, &src));
    assert!(plan(conversion, super::budget()).is_err());
    assert_eq!(std::fs::read(&src).unwrap(), input);
}

#[test]
fn a_checked_plan_runs_repeatedly_and_reports_input_progress_and_all_files() {
    let dir = tempfile::tempdir().unwrap();
    let src = dir.path().join("input.fa");
    let dst = dir.path().join("output.fa.gz");
    let input = b">a\nACGT\n>b\nNacgT\n";
    std::fs::write(&src, input).unwrap();
    let checked = plan(
        Conversion::FastaToBgzipFasta(FastaToBgzipFasta::new(&src, &dst)),
        super::budget(),
    )
    .unwrap();
    assert!(!dst.exists(), "planning must not publish output");
    let mut first = None;
    for attempt in 0..2 {
        let progress = Mutex::new(Vec::new());
        let observe = |event| progress.lock().unwrap().push(event);
        let report = run(&checked, Some(&observe)).unwrap();
        assert_eq!(report.records, 2);
        assert_eq!(report.bytes_read, input.len() as u64);
        assert_eq!(report.spilled, 0);
        assert!(report.warnings.is_empty());
        assert_eq!(
            report.files,
            vec![
                (Format::BgzipFasta, dst.clone()),
                (Format::Fai, dir.path().join("output.fa.gz.fai")),
                (Format::Gzi, dir.path().join("output.fa.gz.gzi")),
            ]
        );
        assert_eq!(inflate(&dst), input, "attempt {attempt}");
        let bytes = std::fs::read(&dst).unwrap();
        if let Some(expected) = &first {
            assert_eq!(&bytes, expected, "attempt {attempt}");
        } else {
            first = Some(bytes);
        }
        let events = progress.into_inner().unwrap();
        assert!(!events.is_empty());
        for event in events {
            assert_eq!(event.phase, "compress");
            assert_eq!(event.bytes_total, input.len() as u64);
            assert!(event.bytes_read > 0 && event.bytes_read <= event.bytes_total);
        }
        assert_eq!(std::fs::read_dir(dir.path()).unwrap().count(), 4);
    }
}
