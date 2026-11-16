use super::*;

#[test]
fn velocity_layers_follow_velocyto_rules() {
    let f = fixture(43);
    let out = f.dir.path().join("v");
    let mut sc = Counts::rna_single_cell();
    sc.velocity = true;
    let mut bulk = Counts::bulk_rna();
    bulk.name = "bulk".into();
    bulk.velocity = true;
    bulk.table = false;
    let t = task(
        &f,
        &out,
        vec![Product::Counts(sc), Product::Counts(bulk)],
        None,
    );
    quantify(t.clone(), Budget::new(5, 256 << 20), None).unwrap();
    for (name, single_cell) in [("genes", true), ("bulk", false)] {
        let want = expected_velocity(&f, single_cell);
        assert!(
            want.iter().all(|l| !l.is_empty()),
            "{name}: every class seen"
        );
        for (i, layer) in ["spliced", "unspliced", "ambiguous"].iter().enumerate() {
            assert_eq!(
                read_10x_layer(&out.join(name), layer),
                want[i],
                "{name}/{layer}"
            );
        }
    }
    support::assert_storage_invariant(
        &out,
        |target| {
            let mut alternate = t.clone();
            alternate.out_dir = target.to_owned();
            alternate
        },
        &[
            ("genes", "spliced"),
            ("genes", "unspliced"),
            ("genes", "ambiguous"),
            ("bulk", "spliced"),
            ("bulk", "unspliced"),
            ("bulk", "ambiguous"),
        ],
    );
}
