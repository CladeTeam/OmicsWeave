mod cfilter;

use cfilter::{Outcome, Run};

#[test]
fn ordered_info_vectors_match_bcftools_with_missing_nan_and_rounding() {
    let mut text = "##fileformat=VCFv4.4\n##contig=<ID=chr1>\n\
        ##INFO=<ID=X,Number=.,Type=Float>\n\
        ##INFO=<ID=Y,Number=.,Type=Float>\n\
        #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
        .to_owned();
    let cases = [
        ("1", "2"),
        ("2", "1"),
        (".", "1"),
        ("1", "."),
        (".", "."),
        ("nan", "1"),
        ("1", "nan"),
        ("nan", "nan"),
        ("inf", "-inf"),
        ("-inf", "inf"),
        ("16777216", "16777217"),
        ("-16777216", "-16777217"),
        ("0", "-0"),
    ];
    for (i, (x, y)) in cases.iter().enumerate() {
        for width in [15, 16, 17, 64] {
            let values = |value: &str| vec![value; width].join(",");
            text.push_str(&format!(
                "chr1\t{}\t.\tA\tC\t.\t.\tX={};Y={}\n",
                i * 100 + width,
                values(x),
                values(y)
            ));
        }
    }
    // Late extrema ensure both true and false outcomes differ from looking
    // only at the first pair; missingness is independent of the numeric range.
    for (i, (x, y)) in [("100,.", "101,0"), ("1,.", "-1,nan")].iter().enumerate() {
        text.push_str(&format!(
            "chr1\t{}\t.\tA\tC\t.\t.\tX={}{x};Y={}{y}\n",
            2000 + i,
            "50,".repeat(32),
            "50,".repeat(32)
        ));
    }
    let dir = tempfile::tempdir().unwrap();
    let path = dir.path().join("vectors.vcf");
    std::fs::write(&path, &text).unwrap();
    for op in ["<", "<=", ">", ">=", "==", "!="] {
        for (left, right) in [
            ("INFO/X", "INFO/Y"),
            ("INFO/Y", "INFO/X"),
            ("INFO/X+0.25", "INFO/Y+0.5"),
            ("INFO/X*0", "INFO/Y*0"),
        ] {
            let expr = format!("{left}{op}{right}");
            let expected = cfilter::c_filter(&path, &expr);
            assert!(
                matches!(&expected, Run::Records(v) if v.len() == cases.len()*4+2 && v.iter().all(|r| matches!(r, Outcome::Pass(_)))),
                "defined oracle required: {expr}: {expected:?}"
            );
            assert_eq!(cfilter::ours(text.as_bytes(), &expr), expected, "{expr}");
        }
    }
}
