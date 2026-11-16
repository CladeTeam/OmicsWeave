use super::*;

#[test]
fn md_and_bismark_tags_give_their_calls() {
    let f = fixture(72);
    let md = f.dir.path().join("md");
    quantify(
        task(
            &f,
            &md,
            product(CPG | CHG | CHH, false, CallSource::Reference),
            false,
            false,
        ),
        Budget::new(4, 256 << 20),
        None,
    )
    .unwrap();
    let want = lines(&cytosines(&f, Mode::Md, false, CPG | CHG | CHH), false);
    assert_eq!(
        gunzip(&md.join("methylation/cytosines.tsv.gz")),
        text(&f, &want)
    );

    // Auto takes XM where it is, with or without a reference.
    for (name, reference) in [("xm", false), ("xm-ref", true)] {
        let out = f.dir.path().join(name);
        quantify(
            task(
                &f,
                &out,
                product(CPG | CHH, true, CallSource::Auto),
                reference,
                true,
            ),
            Budget::new(2, 256 << 20),
            None,
        )
        .unwrap();
        let want = lines(&cytosines(&f, Mode::Xm, true, CPG | CHH), true);
        assert_eq!(
            gunzip(&out.join("methylation/cytosines.tsv.gz")),
            text(&f, &want)
        );
        assert_eq!(
            std::fs::read_to_string(out.join("methylation/mbias.tsv")).unwrap(),
            mbias(&f, Mode::Xm, true)
        );
    }
}
