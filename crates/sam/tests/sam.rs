//! Random SAM texts encoded by htslib (`sam_read1` → `sam_write1` to an
//! uncompressed BAM) and by this crate: the BAM header bytes, every
//! record's bytes and index span, and which lines are rejected must agree.

mod common;

use common::{Encoded, Fixture, Rng, data_lines, htslib, long_cigar_sam, ours, random_sam};
use omicsweave_sam::{BamEncoder, BamScratch, SamHeader, bam_header_bytes};
use proptest::prelude::*;

fn show(line: &[u8]) -> String {
    String::from_utf8_lossy(line).chars().take(300).collect()
}

/// Compares two encodings line by line; lines naming references the header
/// lacks must be rejected here whatever htslib does with them.
fn same(fx: &Fixture, got: &Encoded, want: &Encoded) {
    let text = String::from_utf8_lossy(&fx.text);
    assert_eq!(
        got.header.is_some(),
        want.header.is_some(),
        "header accepted?\n{text}"
    );
    let Some(want_header) = &want.header else {
        return;
    };
    assert_eq!(got.header.as_ref(), Some(want_header), "BAM header\n{text}");
    let mut input = &fx.text[..];
    SamHeader::read(&mut input).unwrap();
    let lines = data_lines(input);
    assert_eq!(got.records.len(), lines.len());
    assert_eq!(want.records.len(), lines.len(), "htslib record count");
    for (i, line) in lines.iter().enumerate() {
        if fx.unknown[i] {
            assert!(
                got.records[i].is_none(),
                "unknown reference accepted: {}",
                show(line)
            );
            continue;
        }
        assert_eq!(
            got.records[i].as_ref().map(|r| &r.1),
            want.records[i].as_ref().map(|r| &r.1),
            "span of line {i}: {}",
            show(line)
        );
        assert!(
            got.records[i] == want.records[i],
            "record of line {i}: {}\nours   {:?}\nhtslib {:?}",
            show(line),
            got.records[i].as_ref().map(|r| &r.0),
            want.records[i].as_ref().map(|r| &r.0),
        );
    }
}

fn run(fx: &Fixture) {
    same(fx, &ours(&fx.text), &htslib(&fx.text, true));
}

proptest! {
    #![proptest_config(ProptestConfig { cases: 2000, ..ProptestConfig::default() })]

    #[test]
    fn random_sam_matches_htslib(seed in any::<u64>()) {
        run(&random_sam(&mut Rng(seed)));
    }
}

proptest! {
    #![proptest_config(ProptestConfig { cases: 24, ..ProptestConfig::default() })]

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
