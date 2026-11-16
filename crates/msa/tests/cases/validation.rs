use omicsweave_core::{Error, Region, Strand};
use omicsweave_msa::MsaReader;

use crate::support::*;

#[test]
fn malformed_schema_and_alignment_rows_are_rejected() {
    for (sql, fields) in [(None, 4), (Some(SQL), 3), (Some("table bed"), 4)] {
        assert!(
            MsaReader::open(bigmaf(&[], sql, fields), None, b'-')
                .err()
                .unwrap()
                .to_string()
                .contains("not a bigMaf")
        );
    }
    for (block, message) in [
        ("s ref.c 0 2 + 100 AC extra;", "8 words, not 7"),
        ("a score=0;e ref.c 0 2 + 100 I;", "without s lines"),
        ("s ref.c 0 3 + 100 A-C;", "2 bases, not 3"),
        (
            "s ref.c 0 3 + 100 ACG;s skipped.c 0 2 + 100 AC;",
            "2 columns, not 3",
        ),
    ] {
        let mut r = reader(&[(0, 3, block.into())], Some(vec!["ref".into()]), b'-');
        assert!(
            r.query(Region::new("c", 0, 3))
                .unwrap_err()
                .to_string()
                .contains(message)
        );
        assert!(matches!(
            r.write(Strand::Forward, &mut [0; 3]),
            Err(Error::NotQueried)
        ));
    }
    assert!(
        MsaReader::open(
            bigmaf(&[], Some(SQL), 4),
            Some(vec!["ref".into(), "ref".into()]),
            b'-'
        )
        .is_err()
    );
}

#[test]
fn whitespace_and_row_metadata_are_checked_even_for_unselected_species() {
    let mut valid = reader(
        &[(
            2,
            5,
            "  s\tref.c\t2 3 + 100 ACG ; s\tdog.c 0 2 - 90 A-C;".into(),
        )],
        None,
        b'.',
    );
    assert_eq!(fetch(&mut valid, 2, 5, Strand::Forward), b"ACGA-C");
    for (start, end, block, message) in [
        (2, 5, "s ref.c 3 3 + 100 ACG;", "reference coordinates"),
        (
            2,
            5,
            "s ref.c 2 3 - 100 ACG;",
            "reference row is not on the + strand",
        ),
        (2, 5, "s ref.c 2 3 + 99 ACG;", "srcSize differs"),
        (2, 5, "s ref.other 2 3 + 100 ACG;", "source differs"),
        (
            2,
            5,
            "s ref.c 2 3 + 100 ACG;s skipped.x 0 3 + 2 ACG;",
            "beyond srcSize",
        ),
        (
            2,
            5,
            "s ref.c 2 3 + 100 ACG;s skipped.x 0 3 ? 90 ACG;",
            "invalid strand",
        ),
        (
            2,
            5,
            "s ref.c 2 3 + 100 ACG;s skipped.x 0 3 + 90 A-C;",
            "2 bases, not 3",
        ),
        (
            2,
            5,
            "s ref.c 2 3 + 100 ACG;s skipped.x 0 3 + 90 éA;",
            "printable ASCII",
        ),
    ] {
        let mut r = reader(
            &[(start, end, block.into())],
            Some(vec!["ref".into()]),
            b'.',
        );
        assert!(
            r.query(Region::new("c", start, end))
                .unwrap_err()
                .to_string()
                .contains(message),
            "{block}"
        );
        assert!(matches!(
            r.write(Strand::Forward, &mut [0; 3]),
            Err(Error::NotQueried)
        ));
    }
    let false_sql = "table bed\n\"mafBlock is only a description\"\n(string chrom;uint chromStart;uint chromEnd;string name;)";
    assert!(MsaReader::open(bigmaf(&[], Some(false_sql), 4), None, b'-').is_err());
}

#[test]
fn autosql_requires_the_fourth_field_declaration_not_quoted_descriptions() {
    for field in ["string mafBlock", "lstring mafBlock"] {
        let sql = format!(
            "table bigMaf\n\"escaped \\\" ; ()\"\n(string chrom;\"x;\" uint chromStart;uint chromEnd;{field};\"mafBlock\")"
        );
        let r = MsaReader::open(bigmaf(&[], Some(&sql), 4), None, b'.').unwrap();
        assert_eq!(r.chroms().len(), 1);
        assert_eq!(r.chroms().get(0).unwrap().name, "c");
        assert_eq!(r.chroms().get(0).unwrap().len, 100);
        assert!(r.species().is_empty());
    }
    for sql in [
        "table t (string chrom;uint start;uint end;string name;) \"mafBlock\"",
        "table t (string chrom;uint start;uint end;uint mafBlock;)",
        "table t (string chrom;uint start;uint end;string mafBlock extra;)",
        "table t (string chrom;uint start;uint end;string mafBlock;string extra;)",
        "table t (string chrom;uint start;uint end;string mafBlock;) \"unterminated",
        "table t (string chrom;uint start;uint end;)",
    ] {
        let error = MsaReader::open(bigmaf(&[], Some(sql), 4), None, b'.')
            .err()
            .unwrap();
        assert!(
            error.to_string().contains(omicsweave_msa::NOT_BIGMAF),
            "{sql}: {error}"
        );
    }
}

#[test]
fn incomplete_and_invalid_sequence_fields_report_the_block_coordinates() {
    for (row, detail) in [
        ("s ref.c 2 3 + 100", "6 words, not 7"),
        ("s ref.c nope 3 + 100 ACG", "invalid start"),
        ("s ref.c 2 nope + 100 ACG", "invalid size"),
        ("s ref.c 2 3 + nope ACG", "invalid srcSize"),
        (
            "s ref.c 18446744073709551615 3 + 18446744073709551615 ACG",
            "beyond srcSize",
        ),
        ("s .c 2 3 + 100 ACG", "empty species"),
        (
            "s ref.c 2 3 + 100 ACG;s other.x 0 3 + 100 AC\u{7f}",
            "printable ASCII",
        ),
    ] {
        let mut r = reader(&[(2, 5, row.into())], Some(vec!["ref".into()]), b'.');
        let error = r.query(Region::new("c", 2, 5)).unwrap_err().to_string();
        assert!(error.contains("c:2-5"), "{row:?}: {error}");
        assert!(
            error.contains(detail),
            "{row:?}: expected {detail:?}, got {error}"
        );
    }
}

#[test]
fn truncated_bbi_input_is_rejected_without_a_valid_reader() {
    let bytes = bigmaf(&[(2, 5, "s ref.c 2 3 + 100 ACG;".into())], Some(SQL), 4);
    for length in [0, 4, 32, 63] {
        assert!(
            MsaReader::open(bytes[..length].to_vec(), None, b'.').is_err(),
            "truncated length={length}"
        );
    }
}

#[test]
fn shape_multiplication_respects_addressable_slice_limits_and_zero_dimensions() {
    use omicsweave_msa::MsaSize;
    for (species, len) in [
        (usize::MAX, 2),
        (2, usize::MAX),
        (1, isize::MAX as usize + 1),
    ] {
        assert!(
            MsaSize { species, len }.elements().is_err(),
            "species={species}, len={len}"
        );
    }
    for (species, len, expected) in [
        (0, usize::MAX, 0),
        (usize::MAX, 0, 0),
        (2, 3, 6),
        (1, isize::MAX as usize, isize::MAX as usize),
    ] {
        assert_eq!(
            MsaSize { species, len }.elements().unwrap(),
            expected,
            "species={species}, len={len}"
        );
    }
}
