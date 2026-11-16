use super::fixtures::Fixture;
use super::reference::Encoded;
use super::{data_lines, htslib, ours};
use omicsweave_sam::SamHeader;

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

pub fn run(fx: &Fixture) {
    same(fx, &ours(&fx.text), &htslib(&fx.text, true));
}
