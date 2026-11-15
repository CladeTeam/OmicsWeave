//! Raw bindings to the decoders of the vendored htscodecs, the entropy
//! codecs of CRAM 3.0 and 3.1 blocks: rANS 4x8 (`rANS_static.h`), rANS
//! Nx16 (`rANS_static4x16.h`), the adaptive arithmetic coder
//! (`arith_dynamic.h`), fqzcomp qualities (`fqzcomp_qual.h`) and the name
//! tokeniser (`tokenise_name3.h`).
//!
//! Functions that return a buffer allocate it with `malloc`; the caller
//! frees it with the C `free`, not a Rust allocator. The `_to` decoders
//! write into a caller buffer of `*out_size` bytes, or allocate with
//! `malloc` when that buffer is NULL. Their successful return is the
//! output pointer and `*out_size` is the decoded byte count. A failed
//! call may modify output bytes and length fields; only a non-NULL
//! return establishes a decoded result.
//!
//! These are raw, unsafe interfaces. They do not validate Rust slice
//! lengths, allocation limits or disjointness. Input pointers must be
//! non-NULL even for zero lengths: the C implementations form an end
//! pointer before handling empty input. The mutable C input prototypes
//! do not transfer ownership of the input.

#![allow(non_camel_case_types)]

use std::ffi::{c_char, c_int, c_uchar, c_uint};

// arith_dynamic calls libbz2; referencing the crate makes rustc link it.
use bzip2_sys as _;

unsafe extern "C" {
    /// A static NUL-terminated version string; do not free or modify it.
    pub fn htscodecs_version() -> *const c_char;

    /// rANS 4x8 (CRAM 3.0 method 4): returns a `malloc`ed buffer of
    /// `*out_size` bytes, or NULL.
    ///
    /// # Safety
    /// `in_` must be non-NULL and readable for `in_size` bytes throughout
    /// the call; `out_size` must point to a writable `c_uint`, disjoint
    /// from the input. A successful result is owned by the caller and
    /// must be released once with C `free`.
    pub fn rans_uncompress(
        in_: *mut c_uchar,
        in_size: c_uint,
        out_size: *mut c_uint,
    ) -> *mut c_uchar;

    /// rANS Nx16 (CRAM 3.1 method 5) into `out`; NULL on failure.
    ///
    /// # Safety
    /// Input and `out_size` follow [`rans_uncompress`]'s requirements.
    /// A non-NULL `out` must be writable for the initial `*out_size`
    /// bytes and disjoint from the input and length pointer. In that
    /// case success returns `out`, whose ownership stays with the caller.
    /// NULL `out` instead requests a C allocation, freed with C `free`.
    pub fn rans_uncompress_to_4x16(
        in_: *mut c_uchar,
        in_size: c_uint,
        out: *mut c_uchar,
        out_size: *mut c_uint,
    ) -> *mut c_uchar;

    /// Adaptive arithmetic coding (CRAM 3.1 method 6) into `out`; NULL on
    /// failure.
    ///
    /// # Safety
    /// The input, output capacity, disjointness and allocation ownership
    /// requirements are those of [`rans_uncompress_to_4x16`]. `out_len`
    /// is the writable capacity/result-length pointer.
    pub fn arith_uncompress_to(
        in_: *mut c_uchar,
        in_size: c_uint,
        out: *mut c_uchar,
        out_len: *mut c_uint,
    ) -> *mut c_uchar;

    /// fqzcomp qualities (CRAM 3.1 method 7): returns a `malloc`ed buffer
    /// of `*out_size` bytes, or NULL. `lengths` may be NULL.
    ///
    /// # Safety
    /// `in_` must be non-NULL and readable for `in_size` bytes, and
    /// `out_size` must point to a writable `usize` disjoint from it.
    /// `nlengths` must be nonnegative; non-NULL `lengths` must be writable
    /// for that many `c_int`s and disjoint from the other arguments.
    /// Release a successful result exactly once with C `free`.
    pub fn fqz_decompress(
        in_: *mut c_char,
        in_size: usize,
        out_size: *mut usize,
        lengths: *mut c_int,
        nlengths: c_int,
    ) -> *mut c_char;

    /// The name tokeniser (CRAM 3.1 method 8): returns a `malloc`ed buffer
    /// of `*out_len` bytes, names NUL-separated, or NULL.
    ///
    /// # Safety
    /// `in_` must be non-NULL and readable for `sz` bytes, and `out_len`
    /// must point to a writable `u32` disjoint from the input. Release a
    /// successful result exactly once with C `free`.
    pub fn tok3_decode_names(in_: *mut u8, sz: u32, out_len: *mut u32) -> *mut u8;
}
