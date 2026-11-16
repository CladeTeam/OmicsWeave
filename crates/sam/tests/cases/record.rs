use crate::support::{Rng, htslib, random_sam, run};
use omicsweave_sam::{BamEncoder, BamScratch, SamHeader};
use proptest::prelude::*;

proptest! {
    #![proptest_config(ProptestConfig { cases: 2000, failure_persistence: Some(Box::new(proptest::test_runner::FileFailurePersistence::Direct(concat!(env!("CARGO_MANIFEST_DIR"), "/tests/cases/record.proptest-regressions")))), ..ProptestConfig::default() })]

    #[test]
    fn random_sam_matches_htslib(seed in any::<u64>()) {
        run(&random_sam(&mut Rng(seed)));
    }
}

/// An RNAME htslib would only warn about and write as unplaced is an
/// error, so the name is never silently dropped.
#[test]
fn unknown_reference_is_an_error() {
    let sam = b"@SQ\tSN:chr1\tLN:100\nr\t0\tchr2\t5\t60\t4M\t*\t0\t0\tACGT\t*\n";
    let mut input = &sam[..];
    let header = SamHeader::read(&mut input).unwrap();
    let encoder = BamEncoder::new(&header);
    let mut out = Vec::new();
    let line = input.strip_suffix(b"\n").unwrap();
    let err = encoder
        .encode_line(&mut BamScratch::new(), line, 2, &mut out)
        .unwrap_err();
    assert!(matches!(
        err,
        omicsweave_core::Error::Format { line: 2, .. }
    ));
    assert!(out.is_empty());
    // htslib writes the record, as unplaced and unmapped.
    let want = htslib(sam, true);
    let (_, span) = want.records[0].as_ref().unwrap();
    assert_eq!((span.tid, span.mapped), (-1, false));
}

#[test]
fn record_bytes_and_span_match_an_independent_bam_fixture() {
    let mut input = b"@SQ\tSN:chr1\tLN:100\n".as_slice();
    let encoder = BamEncoder::new(&SamHeader::read(&mut input).unwrap());
    let mut expected = Vec::new();
    // BAM fixed fields: block size, tid, pos, bin/mapq/name length,
    // flag/CIGAR count, sequence length, mate tid/pos, template length.
    for value in [47u32, 0, 1, 0x1249_3c02, 1, 3, u32::MAX, u32::MAX, 0] {
        expected.extend_from_slice(&value.to_le_bytes());
    }
    expected.extend_from_slice(&[
        b'r', 0, 0x30, 0, 0, 0, 0x12, 0xf0, 40, 0, 2, b'N', b'M', b'C', 1,
    ]);
    let mut output = b"prefix".to_vec();
    let span = encoder
        .encode_line(
            &mut BamScratch::new(),
            b"r\t0\tchr1\t2\t60\t3M\t*\t0\t0\tACN\tI!#\tNM:i:1",
            71,
            &mut output,
        )
        .unwrap();
    assert_eq!(output, [b"prefix".as_slice(), expected.as_slice()].concat());
    assert_eq!(
        span,
        omicsweave_sam::BamSpan {
            tid: 0,
            beg: 1,
            end: 4,
            mapped: true
        }
    );
}

#[test]
fn record_terminator_and_nul_have_explicit_public_error_semantics() {
    let mut input = b"@SQ\tSN:chr1\tLN:100\n".as_slice();
    let encoder = BamEncoder::new(&SamHeader::read(&mut input).unwrap());
    let line = b"r\t0\tchr1\t1\t60\t1M\t*\t0\t0\tA\tI";
    let mut scratch = BamScratch::new();
    let mut plain = Vec::new();
    let span = encoder
        .encode_line(&mut scratch, line, 1, &mut plain)
        .unwrap();
    let mut cr = line.to_vec();
    cr.push(b'\r');
    let mut output = Vec::new();
    assert_eq!(
        encoder
            .encode_line(&mut scratch, &cr, 2, &mut output)
            .unwrap(),
        span
    );
    assert_eq!(output, plain);
    for invalid in [
        b"r\t0\tchr1\t1\t60\t1M\t*\t0\t0\tA\tI\0suffix".as_slice(),
        b"too few fields",
    ] {
        let before = output.clone();
        let error = encoder
            .encode_line(&mut scratch, invalid, 918, &mut output)
            .unwrap_err();
        assert!(
            matches!(error, omicsweave_core::Error::Format { line: 918, ref detail } if !detail.is_empty()),
            "{error:?}"
        );
        assert_eq!(output, before);
        encoder
            .encode_line(&mut scratch, line, 919, &mut output)
            .unwrap();
        assert_eq!(&output[before.len()..], plain);
    }
}

#[test]
fn unknown_mate_reference_is_rejected_without_silently_unplacing_the_mate() {
    let sam = b"@SQ\tSN:chr1\tLN:100\nr\t1\tchr1\t5\t60\t1M\tmissing\t9\t0\tA\tI\n";
    let mut input = sam.as_slice();
    let encoder = BamEncoder::new(&SamHeader::read(&mut input).unwrap());
    let mut output = b"prefix".to_vec();
    let error = encoder
        .encode_line(
            &mut BamScratch::new(),
            input.strip_suffix(b"\n").unwrap(),
            2,
            &mut output,
        )
        .unwrap_err();
    assert!(
        matches!(error, omicsweave_core::Error::Format { line: 2, ref detail } if detail.contains("missing"))
    );
    assert_eq!(output, b"prefix");
    let oracle = htslib(sam, true);
    let record = &oracle.records[0].as_ref().unwrap().0;
    assert_eq!(i32::from_le_bytes(record[24..28].try_into().unwrap()), -1);
}
