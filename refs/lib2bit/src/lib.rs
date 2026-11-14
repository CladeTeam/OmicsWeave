//! Raw bindings to the vendored lib2bit (`vendor/lib2bit/2bit.h`).

use std::ffi::{c_char, c_int};

#[repr(C)]
pub struct TwoBit {
    _opaque: [u8; 0],
}

unsafe extern "C" {
    /// `storeMasked != 0` keeps soft-mask blocks so sequences come back with
    /// lowercase masked bases.
    pub fn twobitOpen(fname: *mut c_char, storeMasked: c_int) -> *mut TwoBit;
    pub fn twobitClose(tb: *mut TwoBit);
    pub fn twobitChromLen(tb: *mut TwoBit, chrom: *mut c_char) -> u32;
    /// Half-open `[start, end)`; the result is malloc'd and NUL-terminated.
    pub fn twobitSequence(tb: *mut TwoBit, chrom: *mut c_char, start: u32, end: u32)
    -> *mut c_char;
}
