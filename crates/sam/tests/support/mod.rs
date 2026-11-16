mod assertions;
mod encoder;
mod fixtures;
mod htslib;
mod reference;
mod rng;

pub use assertions::run;
pub use encoder::{data_lines, ours};
pub use fixtures::{Fixture, cg_line, long_cigar_sam, random_sam};
pub use htslib::htslib;
pub use rng::Rng;

/// Compare complete bytes and spans while retaining the input for diagnostics.
pub fn compare(text: &[u8], exact: bool) {
    let want = htslib(text, exact);
    let got = ours(text);
    assert_eq!(
        got.header,
        want.header,
        "header for {}",
        String::from_utf8_lossy(text)
    );
    assert_eq!(
        got.records,
        want.records,
        "records for {}",
        String::from_utf8_lossy(text)
    );
}
