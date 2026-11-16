use crate::support::{htslib, ours};
use omicsweave_sam::SamHeader;

#[test]
fn late_hd_and_large_headers_match_htslib() {
    for hd_at in [0, 1, 513, 1024] {
        for ending in ["\n", "\r\n"] {
            let mut text = Vec::new();
            for i in 0..1024 {
                if i == hd_at {
                    text.extend_from_slice(format!("@HD\tVN:1.6{ending}").as_bytes());
                }
                text.extend_from_slice(
                    format!("@SQ\tSN:chr{i}\tLN:10000\tAN:alt{i}{ending}").as_bytes(),
                );
                if i % 97 == 0 {
                    text.extend_from_slice(format!("@CO{ending}").as_bytes());
                }
            }
            if hd_at == 1024 {
                text.extend_from_slice(format!("@HD\tVN:1.6{ending}").as_bytes());
            }
            text.extend_from_slice(format!("@HD\tVN:1.4{ending}").as_bytes());
            text.extend_from_slice(b"r\t0\talt1023\t1\t60\t1M\talt0\t2\t0\tA\tI\n");
            let got = ours(&text);
            let expected = htslib(&text, true);
            assert_eq!(got.header, expected.header);
            assert_eq!(got.records, expected.records);
        }
    }
}

#[test]
fn buffered_header_consumption_normalization_and_alias_precedence() {
    use std::io::{BufReader, Read};

    let raw = b"@SQ\tSN:chr1\tLN:4294967296\tAN:shared,alt\r\n@CO\r\n@HD\tVN:1.6\r\n@SQ\tSN:shared\tLN:20\tAN:alt\r\n";
    let body = b"r\t4\t*\t0\t0\t*\t*\t0\t0\t*\t*\n@not_a_header\n";
    let mut sam = raw.to_vec();
    sam.extend_from_slice(body);
    let expected_text = b"@HD\tVN:1.6\n@SQ\tSN:chr1\tLN:4294967296\tAN:shared,alt\n@CO\t\n@SQ\tSN:shared\tLN:20\tAN:alt\n";
    for capacity in [1, 2, 3, 7, 16, 4096] {
        let mut input = BufReader::with_capacity(capacity, sam.as_slice());
        let header = SamHeader::read(&mut input).unwrap();
        assert_eq!(header.raw(), raw, "capacity {capacity}");
        assert_eq!(header.text(), expected_text, "capacity {capacity}");
        assert_eq!(header.lines(), 4);
        assert_eq!(header.references()[0].name, b"chr1");
        assert_eq!(header.references()[0].len, 4_294_967_296);
        assert_eq!(header.references()[1].name, b"shared");
        assert_eq!(header.references()[1].len, 20);
        assert_eq!(header.tid(b"chr1"), Some(0));
        assert_eq!(header.tid(b"alt"), Some(0));
        assert_eq!(header.tid(b"shared"), Some(1));
        assert_eq!(header.tid(b"missing"), None);
        let mut remaining = Vec::new();
        input.read_to_end(&mut remaining).unwrap();
        assert_eq!(remaining, body, "capacity {capacity}");

        let mut expected_bam = b"BAM\x01".to_vec();
        expected_bam.extend_from_slice(&(expected_text.len() as u32).to_le_bytes());
        expected_bam.extend_from_slice(expected_text);
        expected_bam.extend_from_slice(
            b"\x02\0\0\0\x05\0\0\0chr1\0\xff\xff\xff\xff\x07\0\0\0shared\0\x14\0\0\0",
        );
        assert_eq!(omicsweave_sam::bam_header_bytes(&header), expected_bam);
    }
}

#[test]
fn absent_and_unterminated_headers_preserve_the_input_contract() {
    for body in [b"".as_slice(), b"r\t4\t*\t0\t0\t*\t*\t0\t0\t*\t*"] {
        let mut input = body;
        let header = SamHeader::read(&mut input).unwrap();
        assert_eq!(input, body);
        assert!(header.raw().is_empty());
        assert!(header.text().is_empty());
        assert!(header.references().is_empty());
        assert_eq!(header.lines(), 0);
        assert_eq!(
            omicsweave_sam::bam_header_bytes(&header),
            b"BAM\x01\0\0\0\0\0\0\0\0"
        );
    }
    let mut input = b"@SQ\tSN:chr1\tLN:1".as_slice();
    let header = SamHeader::read(&mut input).unwrap();
    assert!(input.is_empty());
    assert_eq!(header.raw(), b"@SQ\tSN:chr1\tLN:1");
    assert_eq!(header.text(), b"@SQ\tSN:chr1\tLN:1\n");
    assert_eq!(header.lines(), 1);
}

#[test]
fn malformed_headers_report_the_physical_line() {
    for (bad, line, detail) in [
        (b"@\n".as_slice(), 2, "two character"),
        (b"@XX\tID:x\n", 2, "line type"),
        (b"@SQ\tSNchr1\tLN:2\n", 2, "key:value"),
        (b"@SQ\tLN:2\n", 2, "no SN"),
        (b"@SQ\tSN:chr1\n", 2, "LN"),
        (b"@SQ\tSN:chr1\tLN:0\n", 2, "LN"),
        (b"@SQ\tSN:chr1\tLN:-1\n", 2, "LN"),
        (b"@SQ\tSN:chr1\tLN:2\tLN:3\n", 2, "different values"),
        (b"@SQ\tSN:chr1\tLN:2\n@SQ\tSN:chr1\tLN:3\n", 3, "duplicate"),
        (b"@RG\tSM:x\n", 2, "ID"),
        (b"@PG\tPN:x\n", 2, "ID"),
    ] {
        let mut sam = b"@CO\tvalid\n".to_vec();
        sam.extend_from_slice(bad);
        let error = SamHeader::read(&mut sam.as_slice()).unwrap_err();
        assert!(
            matches!(error, omicsweave_core::Error::Format { line: got, detail: ref text } if got == line && text.contains(detail)),
            "input {sam:?}: {error:?}"
        );
    }
}

#[test]
fn nul_in_header_is_rejected_as_a_format_error() {
    let mut input = b"@SQ\tSN:chr1\tLN:100\0ignored\n".as_slice();
    let error = SamHeader::read(&mut input).unwrap_err();
    assert!(
        matches!(error, omicsweave_core::Error::Format { line: 1, ref detail } if detail.contains("NUL"))
    );
}

#[test]
fn header_io_failure_retains_its_kind() {
    use std::io::{self, BufReader, Read};

    struct FailingRead;
    impl Read for FailingRead {
        fn read(&mut self, _: &mut [u8]) -> io::Result<usize> {
            Err(io::Error::new(
                io::ErrorKind::BrokenPipe,
                "injected header read failure",
            ))
        }
    }
    let reader = b"@CO\tvalid\n".as_slice().chain(FailingRead);
    let error = SamHeader::read(&mut BufReader::with_capacity(2, reader)).unwrap_err();
    assert!(
        matches!(error, omicsweave_core::Error::Io { kind: io::ErrorKind::BrokenPipe, ref detail } if detail.contains("injected header read failure"))
    );
}
