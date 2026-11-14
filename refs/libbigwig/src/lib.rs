//! Raw bindings to the vendored libBigWig, built with `-DNOCURL`.
//!
//! Declarations mirror `vendor/libBigWig/bigWig.h` for the subset of the API
//! that tests use. Structs whose fields are read from Rust are mirrored with
//! `#[repr(C)]`; the rest stay opaque.

#![allow(non_camel_case_types, non_snake_case)]

use std::ffi::{c_char, c_int, c_void};

// libBigWig calls zlib; referencing the crate makes rustc link it.
use libz_sys as _;

/// `CURLcode (*)(CURL*)`; with `NOCURL`, `CURLcode` is `int` and `CURL` is `void`.
pub type CurlCallback = Option<unsafe extern "C" fn(*mut c_void) -> c_int>;

#[repr(C)]
pub struct URL_t {
    _opaque: [u8; 0],
}

#[repr(C)]
pub struct bwRTree_t {
    _opaque: [u8; 0],
}

#[repr(C)]
pub struct bwWriteBuffer_t {
    _opaque: [u8; 0],
}

#[repr(C)]
pub struct bwZoomHdr_t {
    pub level: *mut u32,
    pub dataOffset: *mut u64,
    pub indexOffset: *mut u64,
    pub idx: *mut *mut bwRTree_t,
}

#[repr(C)]
pub struct bigWigHdr_t {
    pub version: u16,
    pub nLevels: u16,
    pub ctOffset: u64,
    pub dataOffset: u64,
    pub indexOffset: u64,
    pub fieldCount: u16,
    pub definedFieldCount: u16,
    pub sqlOffset: u64,
    pub summaryOffset: u64,
    pub bufSize: u32,
    pub extensionOffset: u64,
    pub zoomHdrs: *mut bwZoomHdr_t,
    pub nBasesCovered: u64,
    pub minVal: f64,
    pub maxVal: f64,
    pub sumData: f64,
    pub sumSquared: f64,
}

#[repr(C)]
pub struct chromList_t {
    pub nKeys: i64,
    pub chrom: *mut *mut c_char,
    pub len: *mut u32,
}

#[repr(C)]
pub struct bigWigFile_t {
    pub URL: *mut URL_t,
    pub hdr: *mut bigWigHdr_t,
    pub cl: *mut chromList_t,
    pub idx: *mut bwRTree_t,
    pub writeBuffer: *mut bwWriteBuffer_t,
    pub isWrite: c_int,
    /// 0 bigWig, 1 bigBed.
    pub r#type: c_int,
}

#[repr(C)]
pub struct bwOverlappingIntervals_t {
    pub l: u32,
    pub m: u32,
    pub start: *mut u32,
    pub end: *mut u32,
    pub value: *mut f32,
}

#[repr(C)]
pub struct bbOverlappingEntries_t {
    pub l: u32,
    pub m: u32,
    pub start: *mut u32,
    pub end: *mut u32,
    pub str: *mut *mut c_char,
}

unsafe extern "C" {
    pub fn bbIsBigBed(fname: *const c_char, callBack: CurlCallback) -> c_int;
    pub fn bwOpen(
        fname: *const c_char,
        callBack: CurlCallback,
        mode: *const c_char,
    ) -> *mut bigWigFile_t;
    pub fn bbOpen(fname: *const c_char, callBack: CurlCallback) -> *mut bigWigFile_t;
    pub fn bwClose(fp: *mut bigWigFile_t);
    pub fn bbGetSQL(fp: *mut bigWigFile_t) -> *mut c_char;

    pub fn bwGetValues(
        fp: *mut bigWigFile_t,
        chrom: *const c_char,
        start: u32,
        end: u32,
        includeNA: c_int,
    ) -> *mut bwOverlappingIntervals_t;
    pub fn bwGetOverlappingIntervals(
        fp: *mut bigWigFile_t,
        chrom: *const c_char,
        start: u32,
        end: u32,
    ) -> *mut bwOverlappingIntervals_t;
    pub fn bbGetOverlappingEntries(
        fp: *mut bigWigFile_t,
        chrom: *const c_char,
        start: u32,
        end: u32,
        withString: c_int,
    ) -> *mut bbOverlappingEntries_t;
    pub fn bwDestroyOverlappingIntervals(o: *mut bwOverlappingIntervals_t);
    pub fn bbDestroyOverlappingEntries(o: *mut bbOverlappingEntries_t);

    pub fn bwCreateHdr(fp: *mut bigWigFile_t, maxZooms: i32) -> c_int;
    pub fn bwCreateChromList(
        chroms: *const *const c_char,
        lengths: *const u32,
        n: i64,
    ) -> *mut chromList_t;
    pub fn bwWriteHdr(bw: *mut bigWigFile_t) -> c_int;
    pub fn bwAddIntervals(
        fp: *mut bigWigFile_t,
        chrom: *const *const c_char,
        start: *const u32,
        end: *const u32,
        values: *const f32,
        n: u32,
    ) -> c_int;
    pub fn bwAddIntervalSpans(
        fp: *mut bigWigFile_t,
        chrom: *const c_char,
        start: *const u32,
        span: u32,
        values: *const f32,
        n: u32,
    ) -> c_int;
    pub fn bwAddIntervalSpanSteps(
        fp: *mut bigWigFile_t,
        chrom: *const c_char,
        start: u32,
        span: u32,
        step: u32,
        values: *const f32,
        n: u32,
    ) -> c_int;
}
