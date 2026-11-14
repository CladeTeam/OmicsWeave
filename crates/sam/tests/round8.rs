mod common;

use common::{htslib, ours};
use omicsweave_sam::{BamEncoder, BamScratch, SamHeader};

fn cg_line(fake_count: usize, real_count: usize, before: usize, after: usize) -> Vec<u8> {
    let mut line = b"q\t0\tchr1\t1\t60\t8S".to_vec();
    for _ in 1..fake_count {
        line.extend_from_slice(b"0N");
    }
    line.extend_from_slice(b"\t*\t0\t0\tACGTACGT\tIIIIIIII\tAA:Z:");
    line.extend(std::iter::repeat_n(b'a', before));
    line.extend_from_slice(b"\tCG:B:I,128");
    for _ in 1..real_count {
        line.extend_from_slice(b",2");
    }
    line.extend_from_slice(b"\tCG:Z:shadow\tZZ:Z:");
    line.extend(std::iter::repeat_n(b'z', after));
    line
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

#[test]
fn failed_cg_compaction_restores_output_and_scratch() {
    let mut input = &b"@SQ\tSN:chr1\tLN:2147483647\n"[..];
    let encoder = BamEncoder::new(&SamHeader::read(&mut input).unwrap());
    let valid = cg_line(16, 65_536, 16_384, 16_384);
    let mut invalid = valid.clone();
    let at = invalid
        .windows(12)
        .position(|s| s == b"CG:B:I,128,2")
        .unwrap();
    invalid[at + 8] = b'4';
    invalid[at + 9] = b'4'; // 144 encodes 9M for an eight-base query.
    let mut scratch = BamScratch::new();
    let mut out = b"prefix".to_vec();
    for _ in 0..3 {
        let before = out.clone();
        assert!(
            encoder
                .encode_line(&mut scratch, &invalid, 2, &mut out)
                .is_err()
        );
        assert_eq!(out, before);
        encoder
            .encode_line(&mut scratch, &valid, 3, &mut out)
            .unwrap();
        let mut expected = b"@SQ\tSN:chr1\tLN:2147483647\n".to_vec();
        expected.extend_from_slice(&valid);
        expected.push(b'\n');
        let oracle = htslib(&expected, true);
        assert_eq!(
            &out[before.len()..],
            oracle.records[0].as_ref().unwrap().0.as_slice()
        );
    }
}

#[test]
fn aux_array_limits_and_incomplete_values_match_htslib() {
    let mut text = b"@SQ\tSN:chr1\tLN:2147483647\n".to_vec();
    for aux in [
        "AA:B:c,-129,128",
        "AA:B:C,-1,255",
        "AA:B:s,-32769,32768",
        "AA:B:S,-1,65535",
        "AA:B:i,-2147483648,2147483647",
        "AA:B:I,0,4294967295",
        "AA:B:i,-2147483649",
        "AA:B:I,4294967296",
        "AA:B:I,-1,4294967295",
        "AA:B:I,",
        "AA:B:I,1,",
        "AA:B:A",
        "AA:B:",
        "AA:i:",
        "AA:f:",
        "AA:A:",
        "AA:Z:",
        "AA:H:",
        "AA:H:0",
        "AA:B:I,16\tCG:B:I,128\tCG:B:I,144",
    ] {
        text.extend_from_slice(b"q\t0\tchr1\t1\t60\t8S\t*\t0\t0\tACGTACGT\tIIIIIIII\t");
        text.extend_from_slice(aux.as_bytes());
        text.push(b'\n');
    }
    let got = ours(&text);
    let want = htslib(&text, true);
    assert_eq!(got.header, want.header);
    assert_eq!(got.records, want.records);
}
