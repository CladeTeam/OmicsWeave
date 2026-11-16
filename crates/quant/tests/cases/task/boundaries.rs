use super::*;

fn coverage_task(dir: &Path, out: &Path, records: bool) -> Task {
    let sam = dir.join("small.sam");
    let mut text = "@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:chr1\tLN:100\n".to_owned();
    if records {
        text.push_str("one\t0\tchr1\t11\t60\t3M\t*\t0\t0\tACG\tIII\n");
    }
    std::fs::write(&sam, text).unwrap();
    let bam = dir.join("small.bam");
    common::bam_of(&sam, &bam);
    let mut task = Task::new(vec![bam], out);
    task.products.push(Product::Coverage(Coverage::default()));
    task
}

fn planning_error(task: Task, budget: Budget, expected: &str) {
    let out = task.out_dir.clone();
    let error = plan(task, budget).unwrap_err().to_string();
    assert!(
        error.contains(expected),
        "expected {expected:?}; got {error}"
    );
    assert!(!out.exists(), "planning created {}", out.display());
}

#[test]
fn invalid_task_shapes_fail_without_creating_output() {
    let dir = tempfile::tempdir().unwrap();
    let out = dir.path().join("out");
    let task = coverage_task(dir.path(), &out, true);
    let mut zero_threads = Budget::new(1, 256 << 20);
    zero_threads.threads = 0;
    planning_error(task.clone(), zero_threads, "threads");
    let mut invalid = task.clone();
    invalid.inputs.clear();
    planning_error(invalid, Budget::new(1, 256 << 20), "no input BAM");
    let mut invalid = task.clone();
    invalid.products.clear();
    planning_error(invalid, Budget::new(1, 256 << 20), "nothing to compute");
    let mut invalid = task.clone();
    invalid
        .products
        .push(Product::Coverage(Coverage::default()));
    planning_error(invalid, Budget::new(1, 256 << 20), "two products");
    let mut invalid = task.clone();
    invalid.inputs.push(task.inputs[0].clone());
    planning_error(invalid, Budget::new(1, 256 << 20), "one BAM");
    let mut invalid = task.clone();
    invalid.products = vec![Product::Counts(Counts::bulk_rna())];
    planning_error(invalid, Budget::new(1, 256 << 20), "annotation");
    let mut invalid = task.clone();
    invalid.inputs = vec![dir.path().join("missing.bam")];
    planning_error(invalid, Budget::new(1, 256 << 20), "missing.bam");

    std::fs::write(&out, b"existing user file").unwrap();
    let error = plan(task, Budget::new(1, 256 << 20))
        .unwrap_err()
        .to_string();
    assert!(error.contains("not a directory"), "{error}");
    assert_eq!(std::fs::read(&out).unwrap(), b"existing user file");
}

#[test]
fn memory_plan_accounts_for_the_budget_before_creating_output() {
    let dir = tempfile::tempdir().unwrap();
    let out = dir.path().join("out");
    let task = coverage_task(dir.path(), &out, true);
    for threads in [1, 6] {
        let budget = Budget::new(threads, 128 << 20);
        let checked = plan(task.clone(), budget.clone()).unwrap();
        assert_eq!(checked.memory.total(), budget.memory, "{threads} threads");
        assert!(checked.memory.pipeline > 0);
        assert!(checked.memory.buffer > 0);
        assert_eq!(checked.stages.len(), 1);
        assert!(checked.stages[0].1 <= checked.memory.buffer);
        assert!(!out.exists(), "planning has no output side effects");
    }
    planning_error(task, Budget::new(1, 0), "needs at least");
}

#[test]
fn empty_sorted_input_delivers_a_readable_zero_coverage_track() {
    let dir = tempfile::tempdir().unwrap();
    let out = dir.path().join("out");
    let task = coverage_task(dir.path(), &out, false);
    let report = quantify(task, Budget::new(2, 128 << 20), None).unwrap();
    assert_eq!(report.records, 0);
    assert!(report.files.iter().all(|(_, p)| p.is_file()));
    let values = common::bigwig_values(&out.join("coverage.bw"), &[("chr1".into(), 100)]);
    assert_eq!(values, vec![vec![0.0; 100]]);
    assert!(out.join("summary.json").is_file());
}

#[test]
fn failed_delivery_preserves_existing_files_and_can_be_retried() {
    let dir = tempfile::tempdir().unwrap();
    let out = dir.path().join("out");
    std::fs::create_dir(&out).unwrap();
    let old = out.join("coverage.bw");
    let blocked = out.join("blocked");
    std::fs::write(&old, b"previous coverage").unwrap();
    std::fs::write(&blocked, b"user-owned path").unwrap();
    let mut task = coverage_task(dir.path(), &out, true);
    task.products.push(Product::Coverage(Coverage {
        name: "blocked/second".into(),
        ..Coverage::default()
    }));
    let error = quantify(task.clone(), Budget::new(2, 128 << 20), None).unwrap_err();
    assert!(!error.to_string().is_empty());
    assert_eq!(std::fs::read(&old).unwrap(), b"previous coverage");
    assert_eq!(std::fs::read(&blocked).unwrap(), b"user-owned path");
    assert_eq!(files(&out).len(), 2, "failed delivery left files behind");

    std::fs::remove_file(&blocked).unwrap();
    std::fs::create_dir(&blocked).unwrap();
    let report = quantify(task, Budget::new(2, 128 << 20), None).unwrap();
    assert_eq!(report.records, 1);
    assert!(out.join("blocked/second.bw").is_file());
    assert!(out.join("summary.json").is_file());
    let values = common::bigwig_values(&old, &[("chr1".into(), 100)]);
    assert_eq!(values[0].iter().filter(|&&x| x == 1.0).count(), 3);
}
