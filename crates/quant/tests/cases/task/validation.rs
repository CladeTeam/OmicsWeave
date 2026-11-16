use super::*;

#[test]
fn unsorted_inputs_and_small_budgets_are_refused() {
    let dir = tempfile::tempdir().unwrap();
    let spec = Spec {
        fragments: 500,
        ..Spec::default()
    };
    let mut reads = synth(&mut Rng::new(1), &spec);
    // Declared sorted, but two records swapped.
    let n = reads.len();
    reads.swap(n / 3, n / 2);
    let bam = write_bam(dir.path(), "swapped", &spec.refs, &reads, 1);
    let out = dir.path().join("o");
    let mut t = Task::new(vec![bam], &out);
    t.products.push(Product::Coverage(Coverage::default()));
    let err = quantify(t.clone(), Budget::new(2, 256 << 20), None)
        .unwrap_err()
        .to_string();
    assert!(err.contains("not sorted by coordinate"), "{err}");
    assert!(err.contains("ow convert"), "{err}");
    assert!(!out.join("coverage.bw").exists());
    assert!(
        std::fs::read_dir(&out)
            .map(|d| d.count() == 0)
            .unwrap_or(true),
        "nothing is left behind"
    );

    // A header without SO:coordinate is refused at plan time.
    let sam = dir.path().join("u.sam");
    std::fs::write(
        &sam,
        "@HD\tVN:1.6\tSO:unsorted\n@SQ\tSN:c\tLN:100\nr\t0\tc\t5\t60\t10M\t*\t0\t0\tACGTACGTAC\t*\n",
    )
    .unwrap();
    let bam = dir.path().join("u.bam");
    common::bam_of(&sam, &bam);
    t.inputs = vec![bam];
    let err = plan(t.clone(), Budget::new(1, 256 << 20))
        .unwrap_err()
        .to_string();
    assert!(
        err.contains("SO:unsorted") && err.contains("--sort"),
        "{err}"
    );

    // A budget below the working set fails before reading.
    t.inputs = vec![dir.path().join("swapped.bam")];
    let err = plan(t, Budget::new(4, 4 << 20)).unwrap_err().to_string();
    assert!(err.contains("needs at least"), "{err}");
}
