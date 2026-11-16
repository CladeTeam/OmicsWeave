//! Independent C-locale BRE oracle and reusable filter state checks.

use omicsweave_core::Error;
use omicsweave_vcf::{Filter, Record, VariantReader};

const HEADER: &str = "##fileformat=VCFv4.4\n##contig=<ID=chr1>\n\
    ##INFO=<ID=TXT,Number=1,Type=String,Description=\"subject\">\n\
    #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n";

fn record(subject: &str) -> (VariantReader<'static>, Record) {
    let text = format!("{HEADER}chr1\t1\t.\tA\tC\t1\t.\tTXT={subject}\n");
    let mut reader = VariantReader::new(std::io::Cursor::new(text.into_bytes())).unwrap();
    let mut record = Record::default();
    assert!(reader.read(&mut record).unwrap());
    (reader, record)
}

fn expression(pattern: &str, icase: bool) -> String {
    format!("INFO/TXT~\"{pattern}{}\"", if icase { "/i" } else { "" })
}

#[test]
fn bre_groups_backrefs_empty_repeats_and_icase_match_libc() {
    let patterns = [
        "a",
        "^a*$",
        "^a\\{0,3\\}$",
        "^a\\{2,\\}$",
        "^a\\?a$",
        "^\\(a\\|aa\\)*b$",
        "^\\(a*\\)\\1$",
        "^\\(a\\|b\\)\\1$",
        "^\\(a*\\)*$",
        "^\\(\\)\\{2,4\\}$",
        "^\\(a\\|\\)b$",
        "[[:lower:]]",
        "[^a-z]",
        "[a-c]",
        "\\ba\\b",
        "\\<a\\>",
        "\\`a\\'",
        "\\a",
        "\\A",
        "[a-]",
        "[]a]",
        "[[=a=]]",
        "[[.a.]]",
        "[[:bogus:]]",
        "\\(a",
        "a\\{3,2\\}",
        "\\1",
        "^aa*b$",
        "\\`aa*b$",
        "^a\\|b",
        "\\(^a\\|b\\)",
        "\\(^a\\|b\\)c",
    ];
    let subjects = [
        "a", "A", "aa", "aaa", "aaaa", "ab", "aab", "b", "!", "a_a", "xxb", "xxbc", "ba",
    ];
    for pattern in patterns {
        for icase in [false, true] {
            for subject in subjects {
                let (reader, rec) = record(subject);
                let expected = crate::regex_oracle::matches(pattern, subject, icase);
                let actual = Filter::new(reader.header(), &expression(pattern, icase));
                match (expected, actual) {
                    (Ok(expected), Ok(mut filter)) => assert_eq!(
                        filter.test(reader.header(), &rec).unwrap(),
                        expected,
                        "pattern={pattern:?}, subject={subject:?}, icase={icase}",
                    ),
                    (Err(()), Err(Error::Invalid { .. })) => {}
                    (expected, actual) => panic!(
                        "compile disagreement: pattern={pattern:?}, subject={subject:?}, \
                         icase={icase}, libc={expected:?}, ours={actual:?}",
                    ),
                }
            }
        }
    }
}

#[test]
fn long_atomic_repeat_uses_bounded_stack_and_keeps_exact_results() {
    let subject = "a".repeat(65_536);
    let (reader, rec) = record(&subject);
    for (pattern, expected) in [("^a*$", true), ("^a*b$", false), ("a*", true)] {
        let mut filter = Filter::new(reader.header(), &expression(pattern, false)).unwrap();
        assert_eq!(filter.test(reader.header(), &rec).unwrap(), expected);
    }
}

#[test]
fn filter_error_does_not_poison_next_call_on_the_same_record_storage() {
    let text = format!(
        "{HEADER}chr1\t1\t.\tA\tC\t2\t.\tTXT=aa\n\
         chr1\t2\t.\tA\tC\t1\t.\tTXT=b\n\
         chr1\t3\t.\tA\tC\t4\t.\tTXT=aa\n",
    );
    let mut reader = VariantReader::new(text.as_bytes()).unwrap();
    let mut filter = Filter::new(reader.header(), "10%(QUAL-1)>0").unwrap();
    let mut rec = Record::default();
    assert!(reader.read(&mut rec).unwrap());
    assert!(!filter.test(reader.header(), &rec).unwrap());
    assert!(reader.read(&mut rec).unwrap());
    assert!(matches!(
        filter.test(reader.header(), &rec),
        Err(Error::Unsupported(_))
    ));
    assert!(reader.read(&mut rec).unwrap());
    assert!(filter.test(reader.header(), &rec).unwrap());
}

#[test]
fn pathological_backtracking_reports_resource_error_then_recovers() {
    let (reader, rec) = record(&"a".repeat(80));
    let (short_reader, short) = record("b");
    let expr = "INFO/TXT~\"^\\(a*\\)*b$\"";
    let mut filter = Filter::new(reader.header(), expr).unwrap();
    let result = filter.test(reader.header(), &rec);
    assert!(
        matches!(
            result,
            Err(Error::Unsupported(_)) | Err(Error::Invalid { .. })
        ),
        "{expr}: {result:?}"
    );
    assert!(filter.test(short_reader.header(), &short).unwrap());
}

#[test]
fn info_negative_regex_follows_bcftools_nonmissing_shortcut() {
    let (reader, long) = record(&"a".repeat(80));
    let (short_reader, short) = record("b");
    let expr = "INFO/TXT!~\"^\\(a*\\)*b$\"";
    let mut filter = Filter::new(reader.header(), expr).unwrap();
    // filter.c::_regex_vector_strings negates its nonmissing precheck and
    // returns true before regexec, so neither subject enters the matcher.
    assert!(filter.test(reader.header(), &long).unwrap());
    assert!(filter.test(short_reader.header(), &short).unwrap());
}

#[test]
fn groups_beyond_nine_preserve_single_digit_backref_semantics() {
    let prefix = "\\(a\\)\\(b\\)\\(c\\)\\(d\\)\\(e\\)\\(f\\)\\(g\\)\\(h\\)\\(i\\)\\(j\\)";
    for (pattern, subjects) in [
        (
            format!("^{prefix}\\(k\\)\\1\\9$"),
            ["abcdefghijkai", "abcdefghijkaj"],
        ),
        (format!("^{prefix}\\10$"), ["abcdefghija0", "abcdefghijj"]),
    ] {
        for subject in subjects {
            let expected = crate::regex_oracle::matches(&pattern, subject, false).unwrap();
            let (reader, rec) = record(subject);
            let mut filter = Filter::new(reader.header(), &expression(&pattern, false)).unwrap();
            assert_eq!(filter.test(reader.header(), &rec).unwrap(), expected);
        }
    }
}

#[test]
fn deep_bre_structure_returns_error_before_recursive_stack_exhaustion() {
    let pattern = format!("{}a{}", "\\(".repeat(4096), "\\)".repeat(4096));
    let (reader, _) = record("a");
    assert!(matches!(
        Filter::new(reader.header(), &expression(&pattern, false)),
        Err(Error::Unsupported(_)) | Err(Error::Invalid { .. })
    ));
}

#[test]
fn anchored_matcher_includes_start_zero_for_empty_subject() {
    let (reader, _) = record("a");
    let rec = Record::default();
    for pattern in ["^a*$", "^a*b$", "\\`a*$", "^\\(a\\|\\)*$"] {
        let expected = crate::regex_oracle::matches(pattern, "", false).unwrap();
        let expr = format!("ID~\"{pattern}\"");
        let mut filter = Filter::new(reader.header(), &expr).unwrap();
        assert_eq!(filter.test(reader.header(), &rec).unwrap(), expected);
    }
}

#[test]
fn empty_repetitions_preserve_participating_captures_for_backrefs() {
    let quantifiers = [
        "\\?",
        "*",
        "\\+",
        "\\{0\\}",
        "\\{1\\}",
        "\\{0,2\\}",
        "\\{2,4\\}",
    ];
    let subjects = ["", "x", "a", "aa", "aaa", "b", "ab", "bb", ".", "b."];
    for quantifier in quantifiers {
        for pattern in [
            format!("^b*\\(\\){quantifier}\\1\\.\\?"),
            format!("^\\(a*\\){quantifier}\\1a$"),
            format!("^\\(a\\|\\){quantifier}\\1b$"),
            format!("^\\(\\(a*\\)\\|\\(b*\\)\\){quantifier}\\2\\3$"),
            format!("^\\(\\(\\)\\?\\){quantifier}\\2$"),
        ] {
            for subject in subjects {
                let expected = crate::regex_oracle::matches(&pattern, subject, false).unwrap();
                let text = format!("{HEADER}chr1\t1\t{subject}\tA\tC\t1\t.\t.\n");
                let mut reader = VariantReader::new(text.as_bytes()).unwrap();
                let mut rec = Record::default();
                if !subject.is_empty() {
                    assert!(reader.read(&mut rec).unwrap());
                }
                let expr = format!("ID~\"{pattern}\"");
                let mut filter = Filter::new(reader.header(), &expr).unwrap();
                assert_eq!(
                    filter.test(reader.header(), &rec).unwrap(),
                    expected,
                    "pattern={pattern:?}, subject={subject:?}",
                );
            }
        }
    }
}

#[test]
fn nullable_backrefs_depend_on_reachable_capture_nodes() {
    let patterns = [
        "^\\(a\\(b*\\)\\|c\\)\\2$",
        "^\\(\\(a\\)\\|\\(b*\\)\\)\\?\\2\\3$",
        "^\\(\\(a*\\)\\|\\(b*\\)\\)\\{0\\}\\2\\3$",
        "^\\(x\\(a*\\)\\|y\\(b*\\)\\)\\?\\2\\3$",
        "^\\(\\(a*\\)\\|\\(b*\\)\\)\\?x\\2\\3$",
        "^\\(\\(^a*\\)\\|\\(b*\\)\\)\\?\\2\\3$",
        "^\\(\\(a\\)\\|\\(b\\)\\)*\\2\\3$",
        "^\\(\\(a*\\)\\|\\(b*\\)\\)\\2\\3$",
        "^\\(\\)\\{0\\}\\1$",
        "^\\(a\\)\\?\\1$",
        "^\\(\\(b*\\)a\\|\\(a*\\)\\)\\2\\3$",
        "^\\(a\\(\\)\\|\\(a\\)\\)\\2\\3$",
        "^\\(\\(a*\\)x\\|b\\)\\2$",
        "^\\(\\(\\b\\)\\|\\(a*\\)\\)\\?\\2\\3$",
        "^\\(\\(\\B\\)\\|\\(a*\\)\\)\\?\\2\\3$",
        "^\\(\\($\\)\\|\\(a*\\)\\)\\?x\\2\\3$",
        "\\(\\(a*\\)\\|\\(b*\\)\\)\\?c\\2\\3$",
        "^\\(a\\(\\)\\|\\(a\\)\\)\\2\\(a\\|\\3\\)$",
        "^\\(a\\(\\)\\|\\(a\\)\\)\\(\\2a\\|\\3\\)$",
        "^\\(\\(a\\)b*\\)*\\2$",
        "^\\(\\(a\\)\\|b\\)*\\2$",
        "^\\(aa\\(\\)\\|\\(aa\\)\\)\\2\\3$",
        "^\\(a\\(\\)a\\|\\(aa\\)\\)\\2\\3$",
    ];
    for pattern in patterns {
        for subject in [
            "", "a", "aa", "aaa", "aaaa", "b", "bb", "ab", "aba", "abab", "c", "x", "y", "qqc",
        ] {
            let expected = crate::regex_oracle::matches(pattern, subject, false).unwrap();
            let text = format!("{HEADER}chr1\t1\t{subject}\tA\tC\t1\t.\t.\n");
            let mut reader = VariantReader::new(text.as_bytes()).unwrap();
            let mut rec = Record::default();
            if !subject.is_empty() {
                assert!(reader.read(&mut rec).unwrap());
            }
            let expr = format!("ID~\"{pattern}\"");
            let mut filter = Filter::new(reader.header(), &expr).unwrap();
            assert_eq!(
                filter.test(reader.header(), &rec).unwrap(),
                expected,
                "pattern={pattern:?}, subject={subject:?}",
            );
        }
    }
}

#[test]
fn cloned_backref_filters_use_c_locale_and_restore_each_callers_locale() {
    struct ThreadLocale {
        active: libc::locale_t,
        previous: libc::locale_t,
    }
    impl ThreadLocale {
        fn utf8() -> Self {
            // SAFETY: a static terminated locale name and valid category mask;
            // the new object belongs to this thread and is restored/freed by
            // Drop after all borrowed filter calls have returned.
            unsafe {
                let active =
                    libc::newlocale(libc::LC_ALL_MASK, c"C.utf8".as_ptr(), std::ptr::null_mut());
                assert!(
                    !active.is_null(),
                    "C.utf8 locale is required for this check"
                );
                let previous = libc::uselocale(active);
                assert!(!previous.is_null());
                Self { active, previous }
            }
        }

        fn assert_active(&self) {
            // SAFETY: a null locale queries this thread without modifying it.
            assert_eq!(
                unsafe { libc::uselocale(std::ptr::null_mut()) },
                self.active
            );
        }
    }
    impl Drop for ThreadLocale {
        fn drop(&mut self) {
            // SAFETY: previous is a live locale or LC_GLOBAL_LOCALE, both valid
            // for uselocale. Only the owned active locale is freed, after it is
            // no longer selected by this thread or referenced by any call.
            unsafe {
                libc::uselocale(self.previous);
                libc::freelocale(self.active);
            }
        }
    }

    let pattern = "^\\([[:alpha:]]\\|[[:digit:]]\\)\\1$";
    let expr = format!("ID~\"{pattern}\"");
    let text = format!("{HEADER}chr1\t1\téé\tA\tC\t1\t.\t.\n");
    let mut reader = VariantReader::new(text.as_bytes()).unwrap();
    let mut rec = Record::default();
    assert!(reader.read(&mut rec).unwrap());
    let header = reader.header();
    let shared = Filter::new(header, &expr).unwrap();
    std::thread::scope(|scope| {
        for _ in 0..4 {
            let mut cloned = shared.clone();
            let expr = &expr;
            let rec = &rec;
            scope.spawn(move || {
                let locale = ThreadLocale::utf8();
                assert!(crate::regex_oracle::matches(pattern, "éé", false).unwrap());
                let mut compiled_here = Filter::new(header, expr).unwrap();
                locale.assert_active();
                for _ in 0..16 {
                    assert!(!cloned.test(header, rec).unwrap());
                    locale.assert_active();
                    assert!(!compiled_here.test(header, rec).unwrap());
                    locale.assert_active();
                }
            });
        }
    });
}

#[test]
fn short_straight_backrefs_match_libc_with_unique_and_empty_captures() {
    let patterns = [
        "^\\(a\\)\\1$",
        "^\\(\\)\\1a$",
        "^\\(\\(a\\)\\)\\1\\2$",
        "^\\(\\)\\(a\\)\\1\\2$",
        "^\\(a\\)b\\1$",
        "\\(a\\)\\1",
        "\\`\\(a\\)\\1\\'",
        "^\\(\\b\\)\\1$",
        "^\\(\\B\\)\\1$",
        "^\\(a\\)\\{0\\}\\1$",
    ];
    for pattern in patterns {
        for subject in ["", "a", "aa", "aaa", "aba", "AA", "Aa", "xxaa", "b"] {
            for icase in [false, true] {
                let expected = crate::regex_oracle::matches(pattern, subject, icase).unwrap();
                let text = format!("{HEADER}chr1\t1\t{subject}\tA\tC\t1\t.\t.\n");
                let mut reader = VariantReader::new(text.as_bytes()).unwrap();
                let mut rec = Record::default();
                if !subject.is_empty() {
                    assert!(reader.read(&mut rec).unwrap());
                }
                let expr = format!("ID~\"{pattern}{}\"", if icase { "/i" } else { "" });
                let mut filter = Filter::new(reader.header(), &expr).unwrap();
                assert_eq!(
                    filter.test(reader.header(), &rec).unwrap(),
                    expected,
                    "pattern={pattern:?}, subject={subject:?}, icase={icase}"
                );
            }
        }
    }
    for padding in [26, 27] {
        let pattern = format!("^\\(a\\)\\1{}$", "b".repeat(padding));
        let subject = format!("aa{}", "b".repeat(padding));
        let expected = crate::regex_oracle::matches(&pattern, &subject, false).unwrap();
        let text = format!("{HEADER}chr1\t1\t{subject}\tA\tC\t1\t.\t.\n");
        let mut reader = VariantReader::new(text.as_bytes()).unwrap();
        let mut rec = Record::default();
        assert!(reader.read(&mut rec).unwrap());
        let mut filter = Filter::new(reader.header(), &format!("ID~\"{pattern}\"")).unwrap();
        assert_eq!(filter.test(reader.header(), &rec).unwrap(), expected);
    }
}

#[test]
fn format_strlen_rejects_uninitialized_gaps_and_recovers_after_overwrite() {
    let mut text = "##fileformat=VCFv4.4\n##contig=<ID=chr1>\n\
        ##INFO=<ID=AF,Number=1,Type=Float>\n\
        ##FORMAT=<ID=GT,Number=1,Type=String>\n\
        ##FORMAT=<ID=LEN,Number=1,Type=String>\n\
        #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS\n\
        chr1\t1\t.\tA\tC\t1\t.\tAF=5.266\tGT\t0/1\n\
        chr1\t2\t.\tA\tC\t1\t.\tAF=5.266\tGT:LEN\t0/1:.\n\
        chr1\t3\t.\tA\tC\t1\t.\tAF=5.266\tLEN\tXX-1Y\n"
        .to_owned();
    text.push_str(&format!(
        "chr1\t4\t.\tA\tC\t1\t.\tAF=5.266\tLEN\t{}\n",
        "Z".repeat(64),
    ));
    text.push_str("chr1\t5\t.\tA\tC\t1\t.\tAF=5.266\tLEN\tY-b1\n");
    let mut reader = VariantReader::new(text.as_bytes()).unwrap();
    let mut filter = Filter::new(reader.header(), "STRLEN(FMT/LEN)>AF").unwrap();
    let mut rec = Record::default();
    for _ in 0..2 {
        assert!(reader.read(&mut rec).unwrap());
        assert!(!filter.test(reader.header(), &rec).unwrap());
    }
    assert!(reader.read(&mut rec).unwrap());
    assert!(matches!(
        filter.test(reader.header(), &rec),
        Err(Error::Unsupported(_))
    ));
    for _ in 0..2 {
        assert!(reader.read(&mut rec).unwrap());
        assert!(filter.test(reader.header(), &rec).unwrap());
    }
}

#[test]
fn format_strlen_keeps_initialized_stale_suffixes() {
    let text = "##fileformat=VCFv4.4\n##contig=<ID=chr1>\n\
        ##INFO=<ID=AF,Number=1,Type=Float>\n\
        ##FORMAT=<ID=LEN,Number=1,Type=String>\n\
        #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS\n\
        chr1\t1\t.\tA\tC\t1\t.\tAF=5.266\tLEN\tABCDEFGH\n\
        chr1\t2\t.\tA\tC\t1\t.\tAF=5.266\tLEN\tY-b1\n";
    let mut reader = VariantReader::new(text.as_bytes()).unwrap();
    let mut filter = Filter::new(reader.header(), "STRLEN(FMT/LEN)>AF").unwrap();
    let mut rec = Record::default();
    for _ in 0..2 {
        assert!(reader.read(&mut rec).unwrap());
        assert!(filter.test(reader.header(), &rec).unwrap());
    }
}

#[test]
fn missing_format_strlen_uses_virtual_capacity_and_recovers_from_int_narrowing() {
    let mut text = "##fileformat=VCFv4.4\n##contig=<ID=chr1>\n\
        ##INFO=<ID=AF,Number=1,Type=Float>\n\
        ##FORMAT=<ID=LEN,Number=1,Type=String>\n\
        #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS\n"
        .to_owned();
    for i in 0..256 {
        text.push_str(&format!("chr1\t{}\t.\tA\tC\t1\t.\tAF=5.266\t.\t.\n", i + 1,));
    }
    text.push_str("chr1\t257\t.\tA\tC\t1\t.\tAF=5.266\tLEN\tABCDEFGH\n");
    let mut reader = VariantReader::new(text.as_bytes()).unwrap();
    let mut filter = Filter::new(reader.header(), "STRLEN(FMT/LEN)>AF").unwrap();
    let mut rec = Record::default();
    let mut unsupported = 0;
    for _ in 0..256 {
        assert!(reader.read(&mut rec).unwrap());
        match filter.test(reader.header(), &rec) {
            Ok(false) => {}
            Err(Error::Unsupported(_)) => unsupported += 1,
            result => panic!("unexpected missing-field result: {result:?}"),
        }
    }
    assert!(unsupported > 0);
    assert!(reader.read(&mut rec).unwrap());
    assert!(filter.test(reader.header(), &rec).unwrap());
}
