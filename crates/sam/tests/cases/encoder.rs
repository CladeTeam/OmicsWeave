use crate::support::{Rng, cg_line, compare, data_lines, htslib};
use omicsweave_sam::{BamEncoder, BamScratch, SamHeader, bam_header_bytes};

/// Workers sharing one encoder, each with its own scratch and buffer,
/// produce the serial stream once their buffers are joined in order.
#[test]
fn encoder_is_shared_across_workers() {
    fn shared<T: Send + Sync>(_: &T) {}
    let mut rng = Rng(7);
    let mut text = b"@SQ\tSN:chr1\tLN:248956422\n@SQ\tSN:chr2\tLN:242193529\n".to_vec();
    for i in 0..4000u64 {
        let pos = rng.range(1, 200_000_000);
        let chrom = if i % 3 == 0 { "chr2" } else { "chr1" };
        text.extend_from_slice(
            format!(
                "r{i}\t{}\t{chrom}\t{pos}\t{}\t50M2I48M\t=\t{}\t300\t{}\t{}\tNM:i:{}\tXB:B:f,1.5,{}\n",
                rng.pick(&[0, 16, 99, 147]),
                rng.range(0, 60),
                pos + 250,
                "ACGT".repeat(25),
                "I".repeat(100),
                rng.range(0, 300),
                rng.range(0, 1000)
            )
            .as_bytes(),
        );
    }
    let mut input = &text[..];
    let header = SamHeader::read(&mut input).unwrap();
    let encoder = BamEncoder::new(&header);
    shared(&encoder);
    let lines = data_lines(input);
    let batches: Vec<Vec<u8>> = std::thread::scope(|s| {
        let workers: Vec<_> = lines
            .chunks(500)
            .enumerate()
            .map(|(b, batch)| {
                let encoder = &encoder;
                s.spawn(move || {
                    let mut scratch = BamScratch::new();
                    let mut out = Vec::new();
                    for (i, line) in batch.iter().enumerate() {
                        let line_no = header_lines(b, i);
                        encoder
                            .encode_line(&mut scratch, line, line_no, &mut out)
                            .unwrap();
                    }
                    out
                })
            })
            .collect();
        workers.into_iter().map(|w| w.join().unwrap()).collect()
    });
    let mut stream = bam_header_bytes(&header);
    batches.iter().for_each(|b| stream.extend_from_slice(b));

    let want = htslib(&text, true);
    let mut serial = want.header.unwrap();
    for r in want.records {
        serial.extend_from_slice(&r.unwrap().0);
    }
    assert!(stream == serial);
}

fn header_lines(batch: usize, i: usize) -> u64 {
    3 + (batch * 500 + i) as u64
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
    let mut valid_text = header.raw().to_vec();
    valid_text.extend_from_slice(valid);
    valid_text.push(b'\n');
    let oracle = htslib(&valid_text, true);
    let expected = oracle.records[0].as_ref().unwrap();
    let mut scratch = BamScratch::new();
    let mut out = b"prefix".to_vec();
    for &line in invalid {
        let before = out.clone();
        let error = encoder
            .encode_line(&mut scratch, line, 2, &mut out)
            .unwrap_err();
        assert!(
            matches!(error, omicsweave_core::Error::Format { line: 2, ref detail } if !detail.is_empty()),
            "input {line:?}: {error:?}"
        );
        assert_eq!(out, before, "input {line:?}");
        let span = encoder
            .encode_line(&mut scratch, valid, 3, &mut out)
            .unwrap();
        assert_eq!(span, expected.1);
        assert_eq!(
            &out[before.len()..],
            expected.0.as_slice(),
            "recovery after {line:?}"
        );
    }
    let mut text = header.raw().to_vec();
    for &line in invalid {
        text.extend_from_slice(line);
        text.push(b'\n');
        text.extend_from_slice(valid);
        text.push(b'\n');
    }
    text.extend_from_slice(b"q\t0\tchr1\t2147483648\t60\t268435455N1M\t*\t0\t0\tA\tI\n");
    compare(&text, false);
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
fn encoder_owns_shared_names_after_header_is_dropped() {
    let encoder = {
        let mut input = &b"@SQ\tSN:chr1\tLN:1000\tAN:alt1\n"[..];
        let header = SamHeader::read(&mut input).unwrap();
        BamEncoder::new(&header).clone()
    };
    let mut output = Vec::new();
    let span = encoder
        .encode_line(
            &mut BamScratch::new(),
            b"r\t0\talt1\t2\t60\t1M\tchr1\t1\t0\tA\tI",
            1,
            &mut output,
        )
        .unwrap();
    assert_eq!((span.tid, span.beg, span.end), (0, 1, 2));
}
