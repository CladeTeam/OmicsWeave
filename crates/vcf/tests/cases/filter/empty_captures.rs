//! Empty repeated captures must participate in following BRE back references.

use crate::filter_oracle;

use omicsweave_vcf::{Filter, Record, VariantReader};

const HEADER: &str =
    "##fileformat=VCFv4.4\n##contig=<ID=chr1>\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n";

#[test]
fn empty_capture_quantifiers_and_backrefs_match_libc_and_bcftools() {
    let patterns = [
        "^b*\\(\\)\\?\\1\\.\\?",
        "^\\(\\)\\?\\1$",
        "^\\(\\)*\\1$",
        "^\\(\\)\\{0\\}\\1$",
        "^\\(\\)\\{1\\}\\1$",
        "^\\(\\)\\{2,4\\}\\1$",
        "^\\(a*\\)\\?\\1$",
        "^\\(a*\\)*\\1$",
        "^\\(\\(a*\\)\\|\\(b*\\)\\)\\?\\2\\3$",
        "^\\(\\(b*\\)a\\|\\(a*\\)\\)\\2\\3$",
        "^\\(\\(a*\\)\\|\\(b*\\)\\)\\?c\\2\\3$",
        "^\\(\\(a\\)\\|\\(a\\)\\)\\?\\2\\3$",
        "^\\(\\(a*\\)\\|\\(b*\\)\\)\\{0\\}\\2\\3$",
        "^\\(\\(\\b\\)\\|\\(a*\\)\\)\\?\\2\\3$",
        "^\\(\\($\\)\\|\\(a*\\)\\)\\?x\\2\\3$",
        "\\(\\(a*\\)\\|\\(b*\\)\\)\\?\\2\\3c",
        "\\(\\(b*\\)a\\|\\(a*\\)\\)\\2\\3",
    ];
    let subjects = [
        "X9.[c_,,", "a", "aa", "aaa", "aaaa", "A", "b", "bb", "b.", ".", "c", "x", "qqc",
    ];
    assert!(crate::regex_oracle::matches(patterns[0], subjects[0], false).unwrap());
    let mut text = HEADER.to_owned();
    for (i, subject) in subjects.iter().enumerate() {
        text.push_str(&format!("chr1\t{}\t{subject}\tA\tC\t.\t.\t.\n", i + 1));
    }
    let dir = tempfile::tempdir().unwrap();
    let path = dir.path().join("empty-captures.vcf");
    std::fs::write(&path, &text).unwrap();
    for pattern in patterns {
        for icase in [false, true] {
            let expression = format!("ID~\"{pattern}{}\"", if icase { "/i" } else { "" });
            let mut reader = VariantReader::new(text.as_bytes()).unwrap();
            let mut filter = Filter::new(reader.header(), &expression).unwrap();
            // A default record's empty ID permits an empty subject without
            // translating the VCF missing-ID representation ('.') into it.
            assert_eq!(
                filter.test(reader.header(), &Record::default()).unwrap(),
                crate::regex_oracle::matches(pattern, "", icase).unwrap(),
                "empty subject: {pattern:?}, icase={icase}"
            );
            let mut record = Record::default();
            for subject in subjects {
                assert!(reader.read(&mut record).unwrap());
                assert_eq!(
                    filter.test(reader.header(), &record).unwrap(),
                    crate::regex_oracle::matches(pattern, subject, icase).unwrap(),
                    "pattern={pattern:?}, subject={subject:?}, icase={icase}"
                );
            }
            let expected = filter_oracle::c_filter(&path, &expression);
            assert!(
                matches!(&expected, filter_oracle::Run::Records(rows) if rows.len() == subjects.len() && rows.iter().all(|row| matches!(row, filter_oracle::Outcome::Pass(_)))),
                "defined oracle required: {expression}: {expected:?}"
            );
            assert_eq!(
                filter_oracle::ours(text.as_bytes(), &expression),
                expected,
                "bcftools: {expression}"
            );
        }
    }
}
