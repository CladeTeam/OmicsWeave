//! Raw bindings to the vendored htslib, for the subset tests use: the FASTA
//! index (`htslib/faidx.h`), BGZF I/O (`htslib/bgzf.h`), VCF/BCF reading,
//! writing, typed value access, index loading and region iteration
//! (`htslib/vcf.h`), SAM/BAM reading and
//! writing, index loading and region iteration (`htslib/sam.h`) and logging
//! (`htslib/hts_log.h`). Handles stay opaque; record and header fields are
//! read through the accessors of `src/fields.c` and `src/bam.c`.

#![allow(non_camel_case_types)]

use std::ffi::{c_char, c_int, c_void};

// htslib calls zlib, libbz2 and liblzma; referencing the crates makes rustc
// link them.
use bzip2_sys as _;
use libz_sys as _;
use lzma_sys as _;

pub type hts_pos_t = i64;

#[repr(C)]
pub struct faidx_t {
    _opaque: [u8; 0],
}

#[repr(C)]
pub struct BGZF {
    _opaque: [u8; 0],
}

#[repr(C)]
pub struct htsFile {
    _opaque: [u8; 0],
}

#[repr(C)]
pub struct bcf_hdr_t {
    _opaque: [u8; 0],
}

#[repr(C)]
pub struct bcf1_t {
    _opaque: [u8; 0],
}

/// `bcf_unpack` levels.
pub const BCF_UN_STR: c_int = 1;
pub const BCF_UN_ALL: c_int = 15;
/// `bcf_get_info_values` / `bcf_get_format_values` value types.
pub const BCF_HT_FLAG: c_int = 0;
pub const BCF_HT_INT: c_int = 1;
pub const BCF_HT_REAL: c_int = 2;
pub const BCF_HT_STR: c_int = 3;

/// `enum htsLogLevel`
pub const HTS_LOG_OFF: c_int = 0;
pub const HTS_LOG_WARNING: c_int = 3;

unsafe extern "C" {
    pub fn hts_set_log_level(level: c_int);

    pub fn fai_build3(fn_: *const c_char, fnfai: *const c_char, fngzi: *const c_char) -> c_int;
    pub fn fai_load3(
        fn_: *const c_char,
        fnfai: *const c_char,
        fngzi: *const c_char,
        flags: c_int,
    ) -> *mut faidx_t;
    pub fn fai_destroy(fai: *mut faidx_t);
    /// `p_end_i` is inclusive; the result is malloc'd and NUL-terminated.
    pub fn faidx_fetch_seq64(
        fai: *const faidx_t,
        c_name: *const c_char,
        p_beg_i: hts_pos_t,
        p_end_i: hts_pos_t,
        len: *mut hts_pos_t,
    ) -> *mut c_char;
    pub fn faidx_nseq(fai: *const faidx_t) -> c_int;
    pub fn faidx_iseq(fai: *const faidx_t, i: c_int) -> *const c_char;
    pub fn faidx_seq_len64(fai: *const faidx_t, seq: *const c_char) -> hts_pos_t;

    pub fn bgzf_open(path: *const c_char, mode: *const c_char) -> *mut BGZF;
    pub fn bgzf_close(fp: *mut BGZF) -> c_int;
    pub fn bgzf_read(fp: *mut BGZF, data: *mut c_void, length: usize) -> isize;
    pub fn bgzf_write(fp: *mut BGZF, data: *const c_void, length: usize) -> isize;
    pub fn bgzf_index_build_init(fp: *mut BGZF) -> c_int;
    pub fn bgzf_index_dump(fp: *mut BGZF, bname: *const c_char, suffix: *const c_char) -> c_int;

    pub fn hts_open(fn_: *const c_char, mode: *const c_char) -> *mut htsFile;
    pub fn hts_close(fp: *mut htsFile) -> c_int;
    pub fn bcf_hdr_read(fp: *mut htsFile) -> *mut bcf_hdr_t;
    pub fn bcf_hdr_write(fp: *mut htsFile, h: *mut bcf_hdr_t) -> c_int;
    pub fn bcf_hdr_destroy(h: *mut bcf_hdr_t);
    /// An empty header; mode `"r"` adds no lines.
    pub fn bcf_hdr_init(mode: *const c_char) -> *mut bcf_hdr_t;
    /// Parses NUL-terminated VCF header text into `hdr`; 0 on success.
    pub fn bcf_hdr_parse(hdr: *mut bcf_hdr_t, htxt: *mut c_char) -> c_int;
    pub fn bcf_init() -> *mut bcf1_t;
    pub fn bcf_destroy(v: *mut bcf1_t);
    /// 0 on success, -1 at the end, < -1 on a malformed record.
    pub fn bcf_read(fp: *mut htsFile, h: *const bcf_hdr_t, v: *mut bcf1_t) -> c_int;
    pub fn bcf_write(fp: *mut htsFile, h: *mut bcf_hdr_t, v: *mut bcf1_t) -> c_int;
    pub fn bcf_unpack(b: *mut bcf1_t, which: c_int) -> c_int;
    /// Writes a CSI index (`fnidx` NULL: `<fn>.csi`); 0 on success.
    pub fn bcf_index_build3(
        fn_: *const c_char,
        fnidx: *const c_char,
        min_shift: c_int,
        n_threads: c_int,
    ) -> c_int;
    pub fn bcf_get_variant_types(rec: *mut bcf1_t) -> c_int;
    /// Values written to a malloc'd `*dst`; returns their count or < 0.
    pub fn bcf_get_format_values(
        hdr: *const bcf_hdr_t,
        line: *mut bcf1_t,
        tag: *const c_char,
        dst: *mut *mut c_void,
        ndst: *mut c_int,
        type_: c_int,
    ) -> c_int;

    /// Values written to a malloc'd `*dst`; returns their count or < 0.
    pub fn bcf_get_info_values(
        hdr: *const bcf_hdr_t,
        line: *mut bcf1_t,
        tag: *const c_char,
        dst: *mut *mut c_void,
        ndst: *mut c_int,
        type_: c_int,
    ) -> c_int;

    pub fn ow_bcf_index_load(fn_: *const c_char, fnidx: *const c_char) -> *mut hts_idx_t;
    pub fn ow_bcf_itr_queryi(
        idx: *const hts_idx_t,
        tid: c_int,
        beg: i64,
        end: i64,
    ) -> *mut hts_itr_t;
    /// >= 0 on success, -1 at the end, < -1 on error.
    pub fn ow_bcf_itr_next(fp: *mut htsFile, itr: *mut hts_itr_t, r: *mut bcf1_t) -> c_int;
    pub fn ow_bcf_hdr_name2id(h: *const bcf_hdr_t, name: *const c_char) -> c_int;

    pub fn ow_bcf_pos(r: *const bcf1_t) -> i64;
    pub fn ow_bcf_rlen(r: *const bcf1_t) -> i64;
    pub fn ow_bcf_n_allele(r: *const bcf1_t) -> c_int;
    /// Requires `bcf_unpack(r, BCF_UN_STR)`.
    pub fn ow_bcf_allele(r: *const bcf1_t, i: c_int) -> *const c_char;
    pub fn ow_bcf_seqname(h: *const bcf_hdr_t, r: *const bcf1_t) -> *const c_char;
    pub fn ow_bcf_hdr_nsamples(h: *const bcf_hdr_t) -> c_int;
    /// Requires `bcf_unpack(r, BCF_UN_STR)`.
    pub fn ow_bcf_id(r: *const bcf1_t) -> *const c_char;
    pub fn ow_bcf_qual(r: *const bcf1_t) -> f32;
    /// Require `bcf_unpack(r, BCF_UN_FLT)`.
    pub fn ow_bcf_n_flt(r: *const bcf1_t) -> c_int;
    pub fn ow_bcf_flt(r: *const bcf1_t) -> *const c_int;
    pub fn ow_bcf_n_sample(r: *const bcf1_t) -> c_int;
    pub fn ow_bcf_n_info(r: *const bcf1_t) -> c_int;
    pub fn ow_bcf_n_fmt(r: *const bcf1_t) -> c_int;
    /// INFO field `i`; requires `bcf_unpack(r, BCF_UN_INFO)`.
    pub fn ow_bcf_info(
        r: *const bcf1_t,
        i: c_int,
        key: *mut c_int,
        type_: *mut c_int,
        len: *mut c_int,
        p: *mut *const u8,
        nbytes: *mut c_int,
    );
    /// FORMAT field `i`; requires `bcf_unpack(r, BCF_UN_FMT)`.
    pub fn ow_bcf_fmt(
        r: *const bcf1_t,
        i: c_int,
        id: *mut c_int,
        type_: *mut c_int,
        n: *mut c_int,
        p: *mut *const u8,
        nbytes: *mut c_int,
    );
    pub fn ow_bcf_hdr_nkeys(h: *const bcf_hdr_t) -> c_int;
    /// NULL for an unused dictionary id.
    pub fn ow_bcf_hdr_key(h: *const bcf_hdr_t, i: c_int) -> *const c_char;

    pub fn free(p: *mut c_void);
}

// SAM/BAM reading and writing (`htslib/sam.h`).

#[repr(C)]
pub struct sam_hdr_t {
    _opaque: [u8; 0],
}

#[repr(C)]
pub struct bam1_t {
    _opaque: [u8; 0],
}

unsafe extern "C" {
    pub fn sam_hdr_read(fp: *mut htsFile) -> *mut sam_hdr_t;
    pub fn sam_hdr_write(fp: *mut htsFile, h: *const sam_hdr_t) -> c_int;
    pub fn sam_hdr_destroy(h: *mut sam_hdr_t);
    pub fn bam_init1() -> *mut bam1_t;
    pub fn bam_destroy1(b: *mut bam1_t);
    /// >= 0 on success, -1 at the end, < -1 on a malformed record.
    pub fn sam_read1(fp: *mut htsFile, h: *mut sam_hdr_t, b: *mut bam1_t) -> c_int;
    /// Bytes written (BAM: 4 + block_size), or < 0 on error.
    pub fn sam_write1(fp: *mut htsFile, h: *const sam_hdr_t, b: *const bam1_t) -> c_int;
    pub fn bam_endpos(b: *const bam1_t) -> hts_pos_t;
    /// Writes a BAI (`min_shift` 0) or CSI index (`fnidx` NULL: beside the
    /// file); 0 on success.
    pub fn sam_index_build3(
        fn_: *const c_char,
        fnidx: *const c_char,
        min_shift: c_int,
        nthreads: c_int,
    ) -> c_int;
}

#[repr(C)]
pub struct hts_idx_t {
    _opaque: [u8; 0],
}

#[repr(C)]
pub struct hts_itr_t {
    _opaque: [u8; 0],
}

// BAM record fields, index loading and region iteration (`src/bam.c`).
unsafe extern "C" {
    pub fn ow_bam_tid(b: *const bam1_t) -> c_int;
    pub fn ow_bam_pos(b: *const bam1_t) -> i64;
    pub fn ow_bam_flag(b: *const bam1_t) -> c_int;
    pub fn ow_bam_mapq(b: *const bam1_t) -> c_int;
    pub fn ow_bam_mtid(b: *const bam1_t) -> c_int;
    pub fn ow_bam_mpos(b: *const bam1_t) -> i64;
    pub fn ow_bam_isize(b: *const bam1_t) -> i64;
    pub fn ow_bam_l_qseq(b: *const bam1_t) -> c_int;
    pub fn ow_bam_n_cigar(b: *const bam1_t) -> c_int;
    pub fn ow_bam_qname(b: *const bam1_t) -> *const c_char;
    pub fn ow_bam_cigar(b: *const bam1_t) -> *const u32;
    pub fn ow_bam_seq(b: *const bam1_t) -> *const u8;
    pub fn ow_bam_qual(b: *const bam1_t) -> *const u8;
    pub fn ow_bam_aux(b: *const bam1_t) -> *const u8;
    pub fn ow_bam_l_aux(b: *const bam1_t) -> c_int;
    pub fn ow_sam_index_load(
        fp: *mut htsFile,
        fn_: *const c_char,
        fnidx: *const c_char,
    ) -> *mut hts_idx_t;
    pub fn ow_sam_itr_queryi(
        idx: *const hts_idx_t,
        tid: c_int,
        beg: i64,
        end: i64,
    ) -> *mut hts_itr_t;
    /// >= 0 on success, -1 at the end, < -1 on error.
    pub fn ow_sam_itr_next(fp: *mut htsFile, itr: *mut hts_itr_t, b: *mut bam1_t) -> c_int;
    pub fn ow_hts_itr_destroy(itr: *mut hts_itr_t);
    pub fn ow_hts_idx_destroy(idx: *mut hts_idx_t);
    pub fn sam_hdr_name2tid(h: *mut sam_hdr_t, name: *const c_char) -> c_int;
}

// Format options (`htslib/hts.h`): `key=value` strings as samtools'
// `--output-fmt-option` takes them, and the reference of a CRAM.

#[repr(C)]
pub struct hts_opt {
    _opaque: [u8; 0],
}

unsafe extern "C" {
    /// Parses one `key=value` option onto the list `opts`; 0 on success.
    pub fn hts_opt_add(opts: *mut *mut hts_opt, c_arg: *const c_char) -> c_int;
    /// Applies the options to an open file; 0 on success.
    pub fn hts_opt_apply(fp: *mut htsFile, opts: *mut hts_opt) -> c_int;
    pub fn hts_opt_free(opts: *mut hts_opt);
    /// Sets the reference FASTA (with its `.fai`) of a CRAM; 0 on success.
    pub fn hts_set_fai_filename(fp: *mut htsFile, fn_aux: *const c_char) -> c_int;
}
