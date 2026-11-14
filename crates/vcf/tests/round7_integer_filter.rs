mod cfilter;
mod round7_fixture;

use cfilter::{Outcome, Run};

#[test]
fn encoded_integer_fields_feed_defined_bcftools_comparisons() {
    let text = round7_fixture::fixture(
        17,
        17,
        4,
        20,
        &[".", "-120", "127", "128", "32768", "2147483647"],
        "4.4",
    );
    let dir = tempfile::tempdir().unwrap();
    let path = dir.path().join("integers.vcf");
    std::fs::write(&path, &text).unwrap();
    for expr in [
        "INFO/I0>=128",
        "FMT/I0>127",
        "FMT/I1==32768",
        "FMT/I2==\".\"",
        "GT=\"het\"",
    ] {
        let expected = cfilter::c_filter(&path, expr);
        assert!(
            matches!(&expected, Run::Records(v) if v.len() == 20 && v.iter().all(|r| matches!(r, Outcome::Pass(_)))),
            "defined oracle required: {expr}: {expected:?}"
        );
        assert_eq!(cfilter::ours(text.as_bytes(), expr), expected, "{expr}");
    }
}
