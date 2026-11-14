mod common;

use common::{htslib, ours};
use omicsweave_sam::{BamEncoder, BamScratch, SamHeader};

fn compare(text: &[u8]) {
    let want = htslib(text, false);
    let got = ours(text);
    assert_eq!(got.header, want.header);
    assert_eq!(got.records, want.records);
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
    compare(&text);
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
        compare(&text);
    }
}

#[test]
fn limits_and_failures_preserve_output_and_reusable_scratch() {
    let mut input = &b"@SQ\tSN:chr1\tLN:2147483647\n"[..];
    let header = SamHeader::read(&mut input).unwrap();
    let encoder = BamEncoder::new(&header);
    let valid = b"q\t0\tchr1\t1\t60\t1M\t*\t0\t0\tA\tI";
    let invalid: &[&[u8]] = &[
        b"q\t0\tchr1\t1\t60\t268435456N1M\t*\t0\t0\tA\tI",
        b"q\t0\tchr1\t2147483649\t60\t1M\t*\t0\t0\tA\tI",
        b"q\t0\tchr1\t1\t60\t1M\t*\t0\t0\tA\t ",
        b"q\t0\tchr1\t1\t60\t1S\t*\t0\t0\tA\tI\tCG:B:I,32",
        b"q\t0\tchr1\t1\t60\t1S\t*\t0\t0\tA\tI\tCG:B:I,16\tXX:B:I,4294967296",
    ];
    let mut scratch = BamScratch::new();
    let mut out = b"prefix".to_vec();
    for &line in invalid {
        let before = out.clone();
        assert!(
            encoder
                .encode_line(&mut scratch, line, 2, &mut out)
                .is_err()
        );
        assert_eq!(out, before);
        encoder
            .encode_line(&mut scratch, valid, 3, &mut out)
            .unwrap();
    }
    let mut text = header.raw().to_vec();
    for &line in invalid {
        text.extend_from_slice(line);
        text.push(b'\n');
        text.extend_from_slice(valid);
        text.push(b'\n');
    }
    text.extend_from_slice(b"q\t0\tchr1\t2147483648\t60\t268435455N1M\t*\t0\t0\tA\tI\n");
    compare(&text);
}
