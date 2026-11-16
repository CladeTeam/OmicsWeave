use crate::filter_oracle;

use filter_oracle::{Outcome, Run};
use omicsweave_vcf::{Filter, Record, VariantReader};

#[test]
fn reused_literal_tokens_match_defined_bcftools_results() {
    let long = "a".repeat(256);
    let mut text = "##fileformat=VCFv4.4\n##contig=<ID=chr1>\n\
        ##INFO=<ID=TXT,Number=.,Type=String>\n\
        ##FORMAT=<ID=ST,Number=1,Type=String>\n\
        #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tA\tB\n"
        .to_owned();
    for (i, (txt, a, b)) in [
        ("alpha", "alpha", "beta"),
        ("beta", "beta", "alpha"),
        (".", ".", "."),
        ("alpha,beta", "alpha", "alpha"),
        (long.as_str(), long.as_str(), "beta"),
        ("beta", "beta", "beta"),
    ]
    .into_iter()
    .enumerate()
    {
        text.push_str(&format!(
            "chr1\t{}\t.\tA\tC\t30\tPASS\tTXT={txt}\tST\t{a}\t{b}\n",
            i + 1
        ));
    }
    let dir = tempfile::tempdir().unwrap();
    let path = dir.path().join("literals.vcf");
    std::fs::write(&path, &text).unwrap();
    for expr in [
        "INFO/TXT=\"alpha\"".to_owned(),
        "INFO/TXT!=\"alpha\"".to_owned(),
        "INFO/TXT=\".\"".to_owned(),
        format!("INFO/TXT=\"{long}\""),
        "FMT/ST=\"alpha\"".to_owned(),
        "FMT/ST!=\"beta\"".to_owned(),
        "INFO/TXT=\"alpha\" || INFO/TXT=\"beta\"".to_owned(),
        "QUAL>20 && INFO/TXT!=\"beta\"".to_owned(),
    ] {
        let expected = filter_oracle::c_filter(&path, &expr);
        let Run::Records(ref outcomes) = expected else {
            panic!("defined oracle required for {expr}: {expected:?}");
        };
        assert_eq!(outcomes.len(), 6, "{expr}: {expected:?}");
        assert!(outcomes.iter().all(|x| matches!(x, Outcome::Pass(_))));
        assert_eq!(
            filter_oracle::ours(text.as_bytes(), &expr),
            expected,
            "{expr}"
        );

        let mut reader = VariantReader::new(text.as_bytes()).unwrap();
        let mut records = Vec::new();
        let mut rec = Record::default();
        while reader.read(&mut rec).unwrap() {
            records.push(rec.clone());
        }
        let mut filter = Filter::new(reader.header(), &expr).unwrap();
        for _ in 0..3 {
            for (record, outcome) in records.iter().zip(outcomes).rev() {
                let Outcome::Pass(wanted) = outcome else {
                    unreachable!();
                };
                assert_eq!(
                    filter.test(reader.header(), record).unwrap(),
                    *wanted,
                    "{expr}"
                );
            }
        }
    }
}
