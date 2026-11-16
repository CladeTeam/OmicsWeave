//! ATAC and methylation deliveries name the indices they actually write.

mod common;

use omicsweave_quant::{Atac, Budget, CallSource, Methylation, Product, Task, quantify};

#[test]
fn sparse_quantification_reports_tbi_or_csi_for_both_tables() {
    let dir = tempfile::tempdir().unwrap();
    for (len, pos, kind) in [
        (1000u32, 100u32, "tbi"),
        (1 << 29, (1 << 29) - 500, "tbi"),
        ((1 << 29) + 1, (1 << 29) - 500, "csi"),
        (830_000_000, 800_000_000, "csi"),
    ] {
        let sam = dir.path().join("reads.sam");
        let text = format!(
            "@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:long\tLN:{len}\n\
             pair\t99\tlong\t{}\t60\t1M\t=\t{}\t201\tC\tI\tCB:Z:BC\tXM:Z:Z\tXG:Z:CT\n\
             pair\t147\tlong\t{}\t60\t1M\t=\t{}\t-201\tC\tI\tCB:Z:BC\tXM:Z:Z\tXG:Z:CT\n",
            pos + 1,
            pos + 201,
            pos + 201,
            pos + 1,
        );
        std::fs::write(&sam, text).unwrap();
        let bam = dir.path().join("reads.bam");
        common::bam_of(&sam, &bam);
        let out = dir.path().join(format!("{kind}-{len}"));
        let mut task = Task::new(vec![bam], &out);
        task.products = vec![
            Product::Atac(Atac {
                tiles: None,
                gene_activity: None,
                ..Atac::default()
            }),
            Product::Methylation(Methylation {
                source: CallSource::Bismark,
                bigwig: false,
                ..Methylation::default()
            }),
        ];
        let report = quantify(task, Budget::new(2, 128 << 20), None).unwrap();
        let summary = std::fs::read_to_string(out.join("summary.json")).unwrap();
        for table in ["atac/fragments.tsv.gz", "methylation/cytosines.tsv.gz"] {
            let rel = format!("{table}.{kind}");
            assert!(summary.contains(&format!("\"{rel}\"")), "{summary}");
            assert!(report.files.contains(&(kind, out.join(&rel))));
            assert!(out.join(&rel).is_file(), "{len}: {rel}");
            let index = common::gunzip_bytes(&out.join(&rel));
            assert_eq!(
                &index[..4],
                if kind == "tbi" {
                    b"TBI\x01"
                } else {
                    b"CSI\x01"
                },
                "{len}: {rel} index encoding"
            );
            let other = if kind == "tbi" { "csi" } else { "tbi" };
            assert!(!summary.contains(&format!("{table}.{other}")));
            assert!(!common::gunzip(&out.join(table)).is_empty());
        }
    }
}
