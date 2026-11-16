use super::*;

#[test]
fn reference_calls_match_the_reads() {
    let f = fixture(71);
    let out = f.dir.path().join("ref");
    quantify(
        task(
            &f,
            &out,
            product(CPG | CHG | CHH, false, CallSource::Reference),
            true,
            true,
        ),
        Budget::new(3, 256 << 20),
        None,
    )
    .unwrap();
    let cyt = cytosines(&f, Mode::Genome, true, CPG | CHG | CHH);
    assert!(cyt.len() > 5000);
    let want = lines(&cyt, false);
    let table = out.join("methylation/cytosines.tsv.gz");
    assert_eq!(gunzip(&table), text(&f, &want));
    assert_tabix_as_htslib(&table, f.dir.path());
    // The tracks hold each line's beta and coverage.
    let beta = bigwig_values(&out.join("methylation/beta.bw"), &f.refs);
    let cov = bigwig_values(&out.join("methylation/coverage.bw"), &f.refs);
    let mut want_beta: Vec<Vec<f32>> = f.refs.iter().map(|r| vec![0.0; r.1 as usize]).collect();
    let mut want_cov = want_beta.clone();
    for &(t, s, _, _, _, m, u) in &want {
        want_cov[t][s as usize] = (m + u) as f32;
        want_beta[t][s as usize] = (m as f64 / (m + u) as f64) as f32;
    }
    assert_eq!(cov, want_cov);
    assert_eq!(beta, want_beta);
    assert_eq!(
        std::fs::read_to_string(out.join("methylation/mbias.tsv")).unwrap(),
        mbias(&f, Mode::Genome, true)
    );
    let (meth, total) = matrix(&f, Mode::Genome, true, 2000, CPG, |r| {
        r.tag("CB").map(str::to_owned)
    });
    assert!(!meth.is_empty());
    assert_eq!(
        read_10x_layer(&out.join("methylation/bins"), "methylated"),
        meth
    );
    assert_eq!(
        read_10x_layer(&out.join("methylation/bins"), "total"),
        total
    );
    let bam = f.bam.display().to_string();
    let (meth, total) = matrix(&f, Mode::Genome, true, 5000, CPG | CHH, |_| {
        Some(bam.clone())
    });
    assert_eq!(
        read_10x_layer(&out.join("methylation/files"), "methylated"),
        meth
    );
    assert_eq!(
        read_10x_layer(&out.join("methylation/files"), "total"),
        total
    );

    support::assert_storage_invariant(
        &out,
        |target| {
            task(
                &f,
                target,
                product(CPG | CHG | CHH, false, CallSource::Reference),
                true,
                true,
            )
        },
        &[
            ("methylation/bins", "methylated"),
            ("methylation/bins", "total"),
            ("methylation/files", "methylated"),
            ("methylation/files", "total"),
        ],
    );

    let summary = std::fs::read_to_string(out.join("summary.json")).unwrap();
    let json: Vec<&str> = summary.lines().collect();
    let methylated: u64 = cyt.values().map(|v| v.0 as u64).sum();
    assert!(summary.contains("\"mbias\""));
    assert!(
        json.iter()
            .any(|l| l.trim() == format!("\"cytosines\": {},", cyt.len())),
        "{summary}"
    );
    assert!(methylated > 0);
    let dropped = json
        .iter()
        .find_map(|l| l.trim().strip_prefix("\"overlap_calls_dropped\": "))
        .unwrap();
    assert_ne!(dropped, "0,", "overlapping mates count once");

    // Without the mask, the SNP positions are called too.
    let open = f.dir.path().join("open");
    quantify(
        task(
            &f,
            &open,
            product(CPG, true, CallSource::Reference),
            true,
            false,
        ),
        Budget::new(2, 256 << 20),
        None,
    )
    .unwrap();
    let want = lines(&cytosines(&f, Mode::Genome, false, CPG), true);
    assert!(want.iter().any(|l| l.4 == '.'));
    assert_eq!(
        gunzip(&open.join("methylation/cytosines.tsv.gz")),
        text(&f, &want)
    );
}
