//! Raw binding to vendored samtools commands, run in process through
//! `src/shim.c`.

use std::ffi::{c_char, c_int};

// samtools calls htslib and zlib; referencing the crates makes rustc link
// them.
use htslib_sys as _;
use libz_sys as _;

unsafe extern "C" {
    /// `samtools <argv[0]> argv[1..]` (`argv[0]` is the command name, such
    /// as `sort`). Returns its exit status, or -1 for a command not built
    /// in. Not reentrant: callers serialise calls.
    pub fn ow_samtools(argc: c_int, argv: *mut *mut c_char) -> c_int;
}
