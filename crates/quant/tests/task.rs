//! Tasks as a whole: one pass for every product, output independent of the
//! thread count, unsorted inputs and small budgets refused before any
//! output appears.

mod common;

use std::path::Path;
use std::sync::Arc;

use common::{Rng, Spec, genes, gtf, synth, write_bam};
use omicsweave_gtf::{Annotation, ReadOptions};
use omicsweave_quant::{
    Atac, Budget, Counts, Coverage, CoverageUnit, Product, StrandSplit, Strandedness, Task, plan,
    quantify,
};

fn files(dir: &Path) -> Vec<(String, Vec<u8>)> {
    let mut out = Vec::new();
    let mut stack = vec![dir.to_owned()];
    while let Some(d) = stack.pop() {
        for e in std::fs::read_dir(&d).unwrap() {
            let p = e.unwrap().path();
            if p.is_dir() {
                stack.push(p);
            } else {
                out.push((
                    p.strip_prefix(dir).unwrap().display().to_string(),
                    std::fs::read(&p).unwrap(),
                ));
            }
        }
    }
    out.sort();
    out
}

#[test]
fn every_product_in_one_pass_is_independent_of_threads() {
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
        assert!(report.records > 0);
        assert!(report.files.iter().any(|(k, _)| *k == "cpk"));
        let all = files(&out);
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
