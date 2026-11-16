use super::*;

#[test]
fn coverage_counts_and_atac_share_a_single_pass_at_any_thread_count() {
    let dir = tempfile::tempdir().unwrap();
    let spec = Spec::default();
    let mut rng = Rng::new(31);
    let reads = synth(&mut rng, &spec);
    let gs = genes(&mut rng, &spec.refs, 8);
    let bam = write_bam(dir.path(), "in", &spec.refs, &reads, 31);
    let g = dir.path().join("g.gtf");
    std::fs::write(&g, gtf(&spec.refs, &gs)).unwrap();
    let ann = Arc::new(Annotation::read(&g, &ReadOptions::default()).unwrap());
    let task = |out: &Path| {
        let mut t = Task::new(vec![bam.clone()], out);
        t.annotation = Some(Arc::clone(&ann));
        t.products = vec![
            Product::Coverage(Coverage::default()),
            Product::Coverage(Coverage {
                name: "frag".into(),
                unit: CoverageUnit::Fragment { max_len: 600 },
                strands: StrandSplit::Split(Strandedness::Reverse),
                ..Coverage::default()
            }),
            Product::Counts(Counts::rna_single_cell()),
            Product::Atac(Atac::default()),
        ];
        t
    };
    let mut reference = None;
    for threads in [1, 2, 7] {
        let out = dir.path().join(format!("t{threads}"));
        let report = quantify(task(&out), Budget::new(threads, 256 << 20), None).unwrap();
        assert_eq!(
            report.records,
            reads.len() as u64,
            "{threads} threads: records read once"
        );
        assert!(report.files.iter().any(|(k, _)| *k == "cpk"));
        let all = files(&out);
        let reported: std::collections::BTreeSet<_> = report
            .files
            .iter()
            .map(|(_, p)| p.strip_prefix(&out).unwrap().display().to_string())
            .collect();
        assert_eq!(
            reported.len(),
            report.files.len(),
            "{threads} threads: duplicate report paths"
        );
        assert_eq!(
            reported,
            all.iter()
                .map(|(p, _)| p.clone())
                .collect::<std::collections::BTreeSet<_>>(),
            "{threads} threads: report names every delivered file"
        );
        assert_eq!(
            report.summary.pretty().as_bytes(),
            std::fs::read(out.join("summary.json")).unwrap(),
            "{threads} threads: report summary"
        );
        match &reference {
            None => reference = Some(all),
            Some(r) => assert_eq!(&all, r, "{threads} threads"),
        }
    }
    let names: Vec<String> = reference.unwrap().into_iter().map(|f| f.0).collect();
    for want in [
        "coverage.bw",
        "frag.forward.bw",
        "frag.reverse.bw",
        "genes/matrix.cpk",
        "atac/fragments.tsv.gz",
        "atac/fragments.tsv.gz.tbi",
        "atac/barcode_metrics.csv",
        "atac/tiles/matrix.cpk",
        "atac/gene_activity/features.tsv.gz",
        "summary.json",
    ] {
        assert!(names.iter().any(|n| n == want), "missing {want}: {names:?}");
    }
}
