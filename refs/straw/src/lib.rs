//! Raw bindings to the vendored straw (Aiden Lab's `.hic` reader) through the
//! C ABI of `src/shim.cpp`.

#![allow(non_camel_case_types)]

use std::ffi::c_char;

// straw calls zlib; referencing the crate makes rustc link it.
use libz_sys as _;

/// Records of one query: `bin_x[i]`, `bin_y[i]` are bin starts in base pairs.
#[repr(C)]
pub struct straw_result {
    pub n: i64,
    pub bin_x: *mut i32,
    pub bin_y: *mut i32,
    pub counts: *mut f32,
}

unsafe extern "C" {
    /// `straw(matrix_type, norm, file, loc1, loc2, unit, bin_size)`; locations
    /// are `chrom:start:end` with inclusive base-pair bounds on bin starts.
    /// NULL on an exception. Release with [`straw_free`].
    pub fn straw_query(
        matrix_type: *const c_char,
        norm: *const c_char,
        file: *const c_char,
        loc1: *const c_char,
        loc2: *const c_char,
        unit: *const c_char,
        bin_size: i32,
    ) -> *mut straw_result;
    pub fn straw_free(r: *mut straw_result);
}
