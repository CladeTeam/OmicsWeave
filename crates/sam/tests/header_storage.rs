mod common;

use common::{htslib, ours};
use omicsweave_sam::{BamEncoder, BamScratch, SamHeader};

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
