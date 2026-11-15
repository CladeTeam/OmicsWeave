mod cfilter;

use cfilter::{Outcome, Run};

#[test]
fn comma_vector_relations_and_borrowed_lookup_match_bcftools() {
    let values = |prefix: &str| {
        (0..64)
            .map(|i| format!("{prefix}{i:03}"))
            .collect::<Vec<_>>()
            .join(",")
    };
    let left = values("left");
    let right = values("right");
    let same = vec!["same"; 64].join(",");
    let missing = vec!["."; 128].join(",");
    let empties = vec![".,,word"; 48].join(",");
    let late = format!("{right},left063");
    // Fully initialized FORMAT storage makes this a deterministic C oracle;
    // width-growth holes are checked separately through the safe error path.
    let rows = [
        (&late, &late),
        (&left, &right),
        (&left, &late),
        (&same, &same),
        (&missing, &missing),
        (&empties, &right),
        (&right, &left),
    ];
    let mut text = "##fileformat=VCFv4.4\n##contig=<ID=chr1>\n\
        ##INFO=<ID=VA,Number=.,Type=String>\n\
        ##INFO=<ID=VB,Number=.,Type=String>\n\
        ##FORMAT=<ID=ST,Number=.,Type=String>\n\
        #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tA\tB\n"
        .to_owned();
    for (i, (a, b)) in rows.into_iter().enumerate() {
        text.push_str(&format!(
            "chr1\t{}\tleft063;other\tA\tC\t30\tPASS\tVA={a};VB={b}\tST\t{a}\t{b}\n",
            i + 1
        ));
    }
    let dir = tempfile::tempdir().unwrap();
    let path = dir.path().join("vectors.vcf");
    std::fs::write(&path, &text).unwrap();
    let lookup = dir.path().join("lookup.txt");
    std::fs::write(&lookup, b"left063\nsame\nword\n").unwrap();
    let mut expressions = vec![
        "INFO/VA=INFO/VB".to_owned(),
        "INFO/VA!=INFO/VB".to_owned(),
        "FMT/ST=INFO/VB".to_owned(),
        "FMT/ST!=INFO/VB".to_owned(),
        "INFO/VA~\"left[0-9]*\"".to_owned(),
        "FMT/ST!~\"right\"".to_owned(),
        "ID=\"left063\"".to_owned(),
        "ID!=\"left063\"".to_owned(),
    ];
    for field in ["INFO/VA", "INFO/VA[0]", "FMT/ST", "FMT/ST[0:0]"] {
        for op in ["=", "!="] {
            expressions.push(format!("{field}{op}@{}", lookup.display()));
        }
    }
    for expression in expressions {
        let oracle = cfilter::c_filter(&path, &expression);
        let Run::Records(outcomes) = &oracle else {
            panic!("{expression}: {oracle:?}");
        };
        assert!(
            outcomes.iter().all(|v| matches!(v, Outcome::Pass(_))),
            "{expression}: {oracle:?}"
        );
        let actual = cfilter::ours(text.as_bytes(), &expression);
        assert_eq!(actual, oracle, "{expression}");
    }
}

#[test]
fn growing_format_shadows_reject_unwritten_reads_and_recover() {
    use omicsweave_core::Error;
    use omicsweave_vcf::{Filter, Record, VariantReader};

    let values = |prefix: &str| {
        (0..64)
            .map(|i| format!("{prefix}{i:03}"))
            .collect::<Vec<_>>()
            .join(",")
    };
    let left = values("left");
    let right = values("right");
    let late = format!("{right},left063");
    let mut text = "##fileformat=VCFv4.4\n##contig=<ID=chr1>\n\
        ##INFO=<ID=VB,Number=.,Type=String>\n\
        ##FORMAT=<ID=ST,Number=.,Type=String>\n\
        #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tA\tB\n"
        .to_owned();
    for (i, (a, b)) in [(&left, &right), (&left, &late), (&right, &right)]
        .into_iter()
        .enumerate()
    {
        text.push_str(&format!(
            "chr1\t{}\t.\tA\tC\t30\tPASS\tVB={b}\tST\t{a}\t{b}\n",
            i + 1
        ));
    }
    let mut reader = VariantReader::new(text.as_bytes()).unwrap();
    let mut record = Record::default();
    let mut filter = Filter::new(reader.header(), "FMT/ST=INFO/VB").unwrap();
    assert!(reader.read(&mut record).unwrap());
    assert!(filter.test(reader.header(), &record).unwrap());
    assert!(reader.read(&mut record).unwrap());
    assert!(matches!(
        filter.test(reader.header(), &record),
        Err(Error::Unsupported(_))
    ));
    assert!(reader.read(&mut record).unwrap());
    assert!(filter.test(reader.header(), &record).unwrap());
}
