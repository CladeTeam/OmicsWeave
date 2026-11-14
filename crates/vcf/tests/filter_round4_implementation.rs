use omicsweave_vcf::{Filter, Record, VariantReader};

#[test]
fn literal_regex_has_bounded_stack_and_linear_search_on_long_subjects() {
    let prefix = "a".repeat(65_536);
    let text = format!(
        "##fileformat=VCFv4.4\n##contig=<ID=chr1>\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\nchr1\t1\t{prefix}b\tA\tC\t.\t.\t.\n"
    );
    let mut reader = VariantReader::new(text.as_bytes()).unwrap();
    let mut rec = Record::default();
    assert!(reader.read(&mut rec).unwrap());
    for (expr, expected) in [
        (format!("ID~\"{}b\"", "a".repeat(4096)), true),
        (format!("ID~\"{}c\"", "a".repeat(4096)), false),
        (format!("ID~\"^{prefix}b$\""), true),
        ("ID~\"^b\"".into(), false),
    ] {
        let mut filter = Filter::new(reader.header(), &expr).unwrap();
        assert_eq!(filter.test(reader.header(), &rec).unwrap(), expected);
    }
}

#[test]
fn fisher_filters_can_evaluate_in_parallel_without_shared_signgam() {
    let text = "##fileformat=VCFv4.4\n##contig=<ID=chr1>\n##INFO=<ID=DP4,Number=4,Type=Integer>\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\nchr1\t1\t.\tA\tC\t.\t.\tDP4=10,20,15,25\n";
    let mut reader = VariantReader::new(text.as_bytes()).unwrap();
    let mut rec = Record::default();
    assert!(reader.read(&mut rec).unwrap());
    let h = reader.header();
    std::thread::scope(|scope| {
        for _ in 0..8 {
            let rec = &rec;
            scope.spawn(move || {
                let mut filter = Filter::new(h, "fisher(INFO/DP4)>=0").unwrap();
                for _ in 0..100 {
                    assert!(filter.test(h, rec).unwrap());
                }
            });
        }
    });
}

#[test]
fn invalid_public_record_and_extreme_alt_index_return_errors_and_recover() {
    let text = "##fileformat=VCFv4.4\n##contig=<ID=chr1>\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\nchr1\t1\t.\tA\tC\t.\t.\t.\n";
    let mut reader = VariantReader::new(text.as_bytes()).unwrap();
    let mut rec = Record::default();
    assert!(reader.read(&mut rec).unwrap());
    let h = reader.header();
    let mut alt = Filter::new(h, "ALT[2147483647]=\"A\"").unwrap();
    assert!(matches!(
        alt.test(h, &rec),
        Err(omicsweave_core::Error::Unsupported(_))
    ));
    for expr in ["REF=\"A\"", "ILEN=0", "POS=1"] {
        let mut filter = Filter::new(h, expr).unwrap();
        assert!(filter.test(h, &rec).unwrap());
        let invalid = if expr == "POS=1" {
            let mut invalid = rec.clone();
            invalid.pos = i64::MAX;
            invalid
        } else {
            Record::default()
        };
        assert!(matches!(
            filter.test(h, &invalid),
            Err(omicsweave_core::Error::Unsupported(_))
        ));
        assert!(filter.test(h, &rec).unwrap());
    }
}

#[test]
fn native_backref_size_errors_and_clone_lifetimes_preserve_reuse() {
    let header =
        "##fileformat=VCFv4.4\n##contig=<ID=chr1>\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n";
    let text = format!(
        "{header}chr1\t1\taa\tA\tC\t.\t.\t.\nchr1\t2\t{}\tA\tC\t.\t.\t.\n",
        "a".repeat(1_000_001)
    );
    let mut reader = VariantReader::new(text.as_bytes()).unwrap();
    let mut short = Record::default();
    let mut long = Record::default();
    assert!(reader.read(&mut short).unwrap());
    assert!(reader.read(&mut long).unwrap());
    let h = reader.header();
    assert!(matches!(
        Filter::new(h, "ID~\"^\\(\\(a\\)\\{32767\\}\\)\\{32767\\}\\1$\""),
        Err(omicsweave_core::Error::Unsupported(_))
    ));
    let mut filter = Filter::new(h, "ID~\"^\\(a*\\)\\1$\"").unwrap();
    assert!(filter.test(h, &short).unwrap());
    assert!(matches!(
        filter.test(h, &long),
        Err(omicsweave_core::Error::Unsupported(_))
    ));
    assert!(filter.test(h, &short).unwrap());
    let mut shared = filter.clone();
    drop(filter);
    for _ in 0..32 {
        let mut copy = shared.clone();
        assert!(copy.test(h, &short).unwrap());
    }
    assert!(shared.test(h, &short).unwrap());
}
