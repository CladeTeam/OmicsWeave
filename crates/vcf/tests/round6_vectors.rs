mod cfilter;

use cfilter::{Outcome, Run};

#[test]
fn wide_equality_vectors_match_bcftools() {
    let mut text = "##fileformat=VCFv4.4\n##contig=<ID=chr1>\n\
        ##INFO=<ID=X,Number=.,Type=Float>\n\
        ##INFO=<ID=Y,Number=.,Type=Float>\n\
        #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
        .to_owned();
    let mut records = 0;
    for width in [4, 63, 64, 65, 256] {
        let uniform = |x: &str| vec![x.to_owned(); width];
        let varied = |parity: usize| {
            (0..width)
                .map(|i| (i * 2 + parity).to_string())
                .collect::<Vec<_>>()
        };
        let mut cases = vec![(varied(0), varied(1))];
        let mut late = varied(0);
        late[width - 1] = (width * 2 - 1).to_string();
        cases.push((late, varied(1)));
        let mut row = varied(0);
        row[0] = (width * 2 - 1).to_string();
        cases.push((row, varied(1)));
        for (x, y) in [
            ("1", "2"),
            ("1", "1"),
            (".", "1"),
            (".", "."),
            ("nan", "nan"),
            ("nan", "1"),
            ("inf", "inf"),
            ("-inf", "-inf"),
            ("16777216", "16777217"),
            ("-16777216", "-16777217"),
            ("0", "-0"),
        ] {
            cases.push((uniform(x), uniform(y)));
        }
        let mut mixed = varied(0);
        mixed[width - 1] = ".".to_owned();
        cases.push((mixed.clone(), uniform(".")));
        mixed[width - 2] = "nan".to_owned();
        cases.push((mixed, varied(1)));
        for (x, y) in cases {
            records += 1;
            text.push_str(&format!(
                "chr1\t{records}\t.\tA\tC\t.\t.\tX={};Y={}\n",
                x.join(","),
                y.join(",")
            ));
        }
    }
    let dir = tempfile::tempdir().unwrap();
    let path = dir.path().join("equality.vcf");
    std::fs::write(&path, &text).unwrap();
    for op in ["==", "!="] {
        for (left, right) in [
            ("INFO/X", "INFO/Y"),
            ("INFO/Y", "INFO/X"),
            ("INFO/X+0.25", "INFO/Y+0.5"),
            ("INFO/X*-1e30", "INFO/Y*-1e30"),
            ("INFO/X*0", "INFO/Y*0"),
        ] {
            let expr = format!("{left}{op}{right}");
            let expected = cfilter::c_filter(&path, &expr);
            assert!(
                matches!(&expected, Run::Records(v) if v.len() == records && v.iter().all(|r| matches!(r, Outcome::Pass(_)))),
                "defined oracle required: {expr}: {expected:?}"
            );
            assert_eq!(cfilter::ours(text.as_bytes(), &expr), expected, "{expr}");
        }
    }
}
