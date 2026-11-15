//! Raw binding to the vendored featureCounts.

use std::ffi::{c_char, c_int};

// featureCounts reads compressed input through zlib; referencing the crate
// makes rustc link it.
use libz_sys as _;

unsafe extern "C" {
    /// `featureCounts argv[1..]` (`argv[0]` is the program name). Returns
    /// its exit status. Not reentrant: getopt and featureCounts' own state
    /// are process-global, so callers serialise calls and reset `optind`.
    pub fn ow_featurecounts_main(argc: c_int, argv: *mut *mut c_char) -> c_int;
}
