//! Raw bindings to the vendored `bcftools consensus` and its -i/-e filter
//! engine (`filter.c`), run in process through `src/shim.c`.

#![allow(non_camel_case_types)]

use std::ffi::{c_char, c_int};

pub use htslib_sys::{bcf_hdr_t, bcf1_t};

#[repr(C)]
pub struct filter_t {
    _opaque: [u8; 0],
}

// bcftools calls htslib; referencing the crate makes rustc link it.
use htslib_sys as _;

unsafe extern "C" {
    /// `bcftools consensus` with `argv` (`argv[0]` is the command name).
    /// Returns its exit status, or -1 when bcftools reports an error. Not
    /// reentrant: callers serialise calls.
    pub fn ow_bcftools_consensus(argc: c_int, argv: *mut *mut c_char) -> c_int;

    /// `filter_init`: exits the process on an invalid expression.
    pub fn filter_init(hdr: *mut bcf_hdr_t, str_: *const c_char) -> *mut filter_t;
    /// `filter_test`; `samples` may be NULL. Exits the process on an error.
    pub fn filter_test(filter: *mut filter_t, rec: *mut bcf1_t, samples: *mut *const u8) -> c_int;
    pub fn filter_destroy(filter: *mut filter_t);

    /// `filter_init`, returning NULL where bcftools reports an error. Not
    /// reentrant. Other failures (assertions, `exit` outside `error()`,
    /// crashes) still end the process.
    pub fn ow_filter_init(hdr: *mut bcf_hdr_t, str_: *const c_char) -> *mut filter_t;
    /// `filter_test(filter, rec, NULL)`: 1 or 0, or -1 where bcftools
    /// reports an error. Not reentrant; other failures end the process.
    pub fn ow_filter_test(filter: *mut filter_t, rec: *mut bcf1_t) -> c_int;
}
