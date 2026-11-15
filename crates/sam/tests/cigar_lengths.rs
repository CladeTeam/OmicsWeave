mod common;

use common::{htslib, ours};

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
