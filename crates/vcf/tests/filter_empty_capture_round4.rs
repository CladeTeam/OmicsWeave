//! Empty repeated captures must participate in following BRE back references.

mod cfilter;

use std::ffi::CString;
use std::mem::MaybeUninit;

use omicsweave_vcf::{Filter, Record, VariantReader};

const HEADER: &str =
    "##fileformat=VCFv4.4\n##contig=<ID=chr1>\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n";

fn libc_match(pattern: &str, subject: &str, icase: bool) -> bool {
    let pattern = CString::new(pattern).unwrap();
    let subject = CString::new(subject).unwrap();
    let mut compiled = MaybeUninit::<libc::regex_t>::uninit();
    let flags = libc::REG_NOSUB | if icase { libc::REG_ICASE } else { 0 };
    // SAFETY: both C strings are live and NUL terminated. Only successful
    // regcomp output is used and freed, once; REG_NOSUB with zero nmatch
    // needs no output array. No reference to the allocation survives regfree.
    unsafe {
        assert_eq!(
            libc::regcomp(compiled.as_mut_ptr(), pattern.as_ptr(), flags),
            0
        );
        let mut compiled = compiled.assume_init();
        let result = libc::regexec(&compiled, subject.as_ptr(), 0, std::ptr::null_mut(), 0);
        libc::regfree(&mut compiled);
        match result {
            0 => true,
            libc::REG_NOMATCH => false,
            other => panic!("regexec failed: {other}"),
        }
    }
}

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
    assert!(libc_match(patterns[0], subjects[0], false));
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
                libc_match(pattern, "", icase),
                "empty subject: {pattern:?}, icase={icase}"
            );
            let mut record = Record::default();
            for subject in subjects {
                assert!(reader.read(&mut record).unwrap());
                assert_eq!(
                    filter.test(reader.header(), &record).unwrap(),
                    libc_match(pattern, subject, icase),
                    "pattern={pattern:?}, subject={subject:?}, icase={icase}"
                );
            }
            assert_eq!(
                cfilter::ours(text.as_bytes(), &expression),
                cfilter::c_filter(&path, &expression),
                "bcftools: {expression}"
            );
        }
    }
}
