use super::*;

#[test]
fn several_input_matrices_keep_file_rows_and_output_restrictions() {
    let dir = tempfile::tempdir().unwrap();
    let mut inputs = Vec::new();
    for name in ["one", "two"] {
        let sam = dir.path().join(format!("{name}.sam"));
        std::fs::write(&sam,"@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:chr1\tLN:1000\nr\t0\tchr1\t101\t60\t3M\t*\t0\t0\tCGC\tIII\tXM:Z:Z.z\tXG:Z:CT\n").unwrap();
        let bam = dir.path().join(format!("{name}.bam"));
        common::bam_of(&sam, &bam);
        inputs.push(bam);
    }
    let out = dir.path().join("out");
    let mut task = Task::new(inputs, &out);
    task.matrix = omicsweave_quant::MatrixFormat::Mtx;
    let spec = Methylation {
        source: CallSource::Bismark,
        table: false,
        bigwig: false,
        overlap: false,
        matrices: vec![MethylMatrix::new("bins", FeatureSource::Bins(1000))],
        ..Methylation::default()
    };
    task.products = vec![Product::Methylation(spec.clone())];
    quantify(task.clone(), Budget::new(2, 128 << 20), None).unwrap();
    let path = out.join("methylation/bins");
    for (layer, want) in [("methylated", 1), ("total", 2)] {
        let matrix = read_10x_layer(&path, layer);
        assert_eq!(matrix.len(), 2);
        assert!(matrix.values().all(|&n| n == want));
        assert_eq!(
            matrix.keys().map(|(r, _)| r).collect::<BTreeSet<_>>().len(),
            2
        );
    }
    let mut forbidden = spec;
    forbidden.bigwig = true;
    task.products = vec![Product::Methylation(forbidden)];
    assert!(
        plan(task, Budget::new(1, 128 << 20))
            .unwrap_err()
            .to_string()
            .contains("one BAM")
    );
}
