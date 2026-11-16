use super::*;

fn lines(p: &Path) -> Vec<String> {
    std::fs::read_to_string(p)
        .unwrap()
        .lines()
        .map(str::to_owned)
        .collect()
}

#[test]
fn junction_table_matches_star_columns() {
    let f = fixture(41);
    let want = expected_table(&f, true);
    assert!(want.iter().any(|l| l.split('\t').nth(5) == Some("1")));
    for m in 1..=6 {
        assert!(
            want.iter()
                .any(|l| l.split('\t').nth(4) == Some(&m.to_string())),
            "motif {m} seen"
        );
    }
    let mut outputs = Vec::new();
    for (i, reference) in [&f.fasta, &f.twobit].into_iter().enumerate() {
        let out = f.dir.path().join(format!("j{i}"));
        let t = task(
            &f,
            &out,
            vec![Product::Junctions(Junctions::default())],
            Some(reference),
        );
        quantify(t, Budget::new(1 + 4 * i, 256 << 20), None).unwrap();
        let got = lines(&out.join("junctions/junctions.tab"));
        assert_eq!(got, want, "reference {}", reference.display());
        outputs.push(std::fs::read(out.join("junctions/junctions.tab")).unwrap());
    }
    assert_eq!(outputs[0], outputs[1]);
    // Without a reference: no motifs, strands from the annotation.
    let out = f.dir.path().join("nomotif");
    let t = task(
        &f,
        &out,
        vec![Product::Junctions(Junctions::default())],
        None,
    );
    quantify(t, Budget::new(3, 256 << 20), None).unwrap();
    assert_eq!(
        lines(&out.join("junctions/junctions.tab")),
        expected_table(&f, false)
    );
}

#[test]
fn junction_matrix_counts_umis_per_barcode() {
    let f = fixture(42);
    let out = f.dir.path().join("m");
    let j = Junctions {
        matrix: Some(JunctionMatrix::single_cell()),
        ..Junctions::default()
    };
    let mut t = task(&f, &out, vec![Product::Junctions(j)], Some(&f.fasta));
    t.matrix = MatrixFormat::Mtx;
    quantify(t.clone(), Budget::new(4, 256 << 20), None).unwrap();

    let matrix_filter = |r: &Read| r.flag & 0x704 == 0 && r.nh() == 1;
    // Barcodes in order of first appearance among reads counted in the
    // matrix; units take the row and UMI of their first such read.
    let mut rows: Vec<String> = Vec::new();
    for r in &f.reads {
        if table_filter(r)
            && matrix_filter(r)
            && let Some(cb) = r.tag("CB")
            && !rows.iter().any(|x| x == cb)
        {
            rows.push(cb.to_owned());
        }
    }
    let mut seen = BTreeSet::new();
    for u in units(&f.reads, table_filter) {
        let Some(first) = u
            .iter()
            .find(|r| matrix_filter(r) && r.tag("CB").is_some() && r.tag("UB").is_some())
        else {
            continue;
        };
        let (cb, ub) = (first.tag("CB").unwrap(), first.tag("UB").unwrap());
        for (t, s, e) in unit_junctions(&u).into_keys() {
            seen.insert((
                cb.to_owned(),
                format!("{}:{s}-{e}", f.spec.refs[t].0),
                ub.to_owned(),
            ));
        }
    }
    let mut want = BTreeMap::new();
    for (cb, j, _) in seen {
        *want.entry((cb, j)).or_insert(0u32) += 1;
    }
    let dir = out.join("junctions");
    assert_eq!(read_10x(&dir), want);
    let barcodes: Vec<String> = gunzip(&dir.join("barcodes.tsv.gz"))
        .lines()
        .map(str::to_owned)
        .collect();
    assert_eq!(barcodes, rows);
    // One column per junction of the table, in its order.
    let table = lines(&dir.join("junctions.tab"));
    let features = gunzip(&dir.join("features.tsv.gz"));
    assert_eq!(features.lines().count(), table.len());
    for (l, t) in features.lines().zip(&table) {
        let fv: Vec<&str> = l.split('\t').collect();
        let tv: Vec<&str> = t.split('\t').collect();
        assert_eq!(fv[2], "Splice Junction");
        assert_eq!(fv[3], tv[0]);
        assert_eq!(
            fv[4].parse::<u32>().unwrap() + 1,
            tv[1].parse::<u32>().unwrap()
        );
        assert_eq!(fv[5], tv[2]);
    }
    support::assert_storage_invariant(
        &out,
        |target| {
            let mut alternate = t.clone();
            alternate.out_dir = target.to_owned();
            alternate
        },
        &[("junctions", "matrix")],
    );
}
