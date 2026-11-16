use crate::support::{Fixture, Rng, cg_line, compare, htslib, long_cigar_sam, ours, run};
use proptest::prelude::*;

proptest! {
    #![proptest_config(ProptestConfig { cases: 24, failure_persistence: Some(Box::new(proptest::test_runner::FileFailurePersistence::Direct(concat!(env!("CARGO_MANIFEST_DIR"), "/tests/cases/cigar.proptest-regressions")))), ..ProptestConfig::default() })]

    #[test]
    fn long_cigars_match_htslib(seed in any::<u64>()) {
        run(&long_cigar_sam(&mut Rng(seed)));
    }
}

/// More than 65535 operations spanning 2^28 or more reference bases fit
/// neither BAM's CIGAR nor its placeholder: htslib fails to write the
/// record, and so does the encoder.
#[test]
fn long_cigar_beyond_placeholder_is_rejected() {
    let ops = "4096D1M".repeat(65_536);
    let sam = format!(
        "@SQ\tSN:chr1\tLN:2000000000\nr\t0\tchr1\t1\t60\t{ops}\t*\t0\t0\t{}\t*\n",
        "A".repeat(65_536)
    );
    let fx = Fixture {
        text: sam.into_bytes(),
        unknown: vec![false],
    };
    let want = htslib(&fx.text, false);
    assert!(want.records[0].is_none());
    let got = ours(&fx.text);
    assert!(got.records[0].is_none());
}

#[test]
fn mixed_cigar_lengths_and_long_placeholders_match_htslib() {
    let mut text = b"@SQ\tSN:chr1\tLN:2147483647\n".to_vec();
    // Each group consumes 25 query and 25 reference bases. H, P and B
    // consume neither; I/S and D/N exercise the two lengths separately.
    for repeats in [1, 17, 6554] {
        for flag in [0, 4] {
            let query = 25 * repeats;
            let line = format!(
                "q\t{flag}\tchr1\t123\t60\t{}\t*\t0\t0\t{}\t*\tNM:i:0\n",
                "1M2I3D4N5S6H7P8=9X10B".repeat(repeats),
                "A".repeat(query),
            );
            text.extend_from_slice(line.as_bytes());
        }
    }
    let got = ours(&text);
    let expected = htslib(&text, true);
    assert_eq!(got.header, expected.header);
    assert!(got.records.iter().all(Option::is_some));
    assert_eq!(got.records, expected.records);
}
#[test]
fn first_cg_tag_controls_migration() {
    let mut text = b"@SQ\tSN:chr1\tLN:2147483647\n".to_vec();
    for (cigar, aux) in [
        ("4S", "CG:Z:shadow\tCG:B:I,64"),
        ("4S", "CG:B:C,64\tCG:B:I,64"),
        ("4S", "CG:B:I,64\tCG:B:I,68"),
        ("4S", "CG:B:i,64\tCG:Z:remaining"),
        ("4S1N", "CG:B:I,64\tCG:B:I,32,32"),
        (
            "4S",
            "AA:Z:before\tAA:B:c,255\tBB:d:0x1.8p1\tCG:B:I,64\tAA:Z:after",
        ),
        ("4S", "CG:B:I\tCG:B:I,64"),
        ("4M", "CG:B:I,68"),
    ] {
        text.extend_from_slice(
            format!("q\t0\tchr1\t1\t60\t{cigar}\t*\t0\t0\tACGT\tIIII\t{aux}\n").as_bytes(),
        );
    }
    compare(&text, false);
}

#[test]
fn long_cg_migration_preserves_duplicate_tags_and_aux_order() {
    for n in [65_535usize, 65_536] {
        let mut text = b"@SQ\tSN:chr1\tLN:2147483647\n".to_vec();
        text.extend_from_slice(format!("q\t0\tchr1\t1\t60\t{n}S\t*\t0\t0\t").as_bytes());
        text.extend(std::iter::repeat_n(b'A', n));
        text.extend_from_slice(b"\t*\tAA:Z:before\tCG:B:I");
        for _ in 0..n {
            text.extend_from_slice(b",16");
        }
        text.extend_from_slice(b"\tCG:Z:after\tAA:i:7\n");
        compare(&text, false);
    }
}

#[test]
fn cg_compaction_matches_htslib_across_inline_boundary() {
    let mut text = b"@SQ\tSN:chr1\tLN:2147483647\n".to_vec();
    for (fake, real) in [
        (1, 1),
        (1, 2),
        (1, 65_535),
        (1, 65_536),
        (16, 16),
        (16, 65_536),
        (65_535, 65_536),
        (65_536, 65_536),
    ] {
        for (before, after) in [(0, 0), (16_384, 0), (0, 16_384), (16_384, 16_384)] {
            text.extend(cg_line(fake, real, before, after));
            text.push(b'\n');
        }
    }
    let got = ours(&text);
    let want = htslib(&text, true);
    assert_eq!(got.header, want.header);
    assert!(got.records.iter().all(Option::is_some));
    assert_eq!(got.records, want.records);
}
