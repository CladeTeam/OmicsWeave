use crate::support::writer_fixture::{bed, wig};
use omicsweave_bbi::BbiSidecar;
#[test]
fn disabled_derived_zooms_match_explicit_empty_zooms_byte_for_byte() {
    for n in [0, 1, 9_999, 10_000, 20_000] {
        for (bytes, expected) in [
            (wig(n, false, 0), wig(n, true, 0)),
            (bed(n, false), bed(n, true)),
        ] {
            assert_eq!(bytes, expected);
            assert_eq!(BbiSidecar::build(&bytes.0).unwrap(), bytes.1);
        }
    }
}
