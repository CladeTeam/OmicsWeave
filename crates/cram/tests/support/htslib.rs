//! Independent htslib CRAM writer, indexer, and record oracle.

use super::path_cstring;
use htslib_sys as h;
use std::ffi::{CStr, CString};
use std::path::Path;

/// Opens `path` for reading with htslib, with `reference` for CRAM.
unsafe fn open_read(path: &Path, reference: Option<&Path>) -> *mut h::htsFile {
    unsafe {
        let fp = h::hts_open(path_cstring(path).as_ptr(), c"r".as_ptr());
        assert!(!fp.is_null(), "open {path:?}");
        if let Some(r) = reference {
            assert_eq!(h::hts_set_fai_filename(fp, path_cstring(r).as_ptr()), 0);
        }
        fp
    }
}

/// The CRAM htslib writes for the SAM at `src` against `reference`, with
/// format options (`version=3.0`, `embed_ref=1`, ...).
pub fn cram_of(src: &Path, out: &Path, reference: Option<&Path>, options: &[&str]) {
    unsafe {
        let r = open_read(src, None);
        let hdr = h::sam_hdr_read(r);
        assert!(!hdr.is_null());
        let w = h::hts_open(path_cstring(out).as_ptr(), c"wc".as_ptr());
        assert!(!w.is_null());
        if let Some(reference) = reference {
            assert_eq!(
                h::hts_set_fai_filename(w, path_cstring(reference).as_ptr()),
                0
            );
        }
        let mut opts: *mut h::hts_opt = std::ptr::null_mut();
        for o in options {
            let o = CString::new(*o).unwrap();
            assert_eq!(h::hts_opt_add(&mut opts, o.as_ptr()), 0, "{o:?}");
        }
        assert_eq!(h::hts_opt_apply(w, opts), 0);
        h::hts_opt_free(opts);
        assert_eq!(h::sam_hdr_write(w, hdr), 0);
        let rec = h::bam_init1();
        loop {
            let status = h::sam_read1(r, hdr, rec);
            if status < 0 {
                assert_eq!(status, -1, "htslib read {src:?}, options={options:?}");
                break;
            }
            assert!(h::sam_write1(w, hdr, rec) >= 0);
        }
        h::bam_destroy1(rec);
        h::sam_hdr_destroy(hdr);
        h::hts_close(r);
        assert_eq!(h::hts_close(w), 0);
    }
}

/// htslib's `.crai` of the CRAM at `file`, written to `out`.
pub fn crai_of(file: &Path, out: &Path) {
    unsafe {
        assert_eq!(
            h::sam_index_build3(
                path_cstring(file).as_ptr(),
                path_cstring(out).as_ptr(),
                0,
                0
            ),
            0
        );
    }
}

/// Every field of a record as htslib holds it.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Rec {
    pub tid: i32,
    pub pos: i64,
    pub end: i64,
    pub flag: u16,
    pub mapq: u8,
    pub mtid: i32,
    pub mpos: i64,
    pub isize: i64,
    pub name: Vec<u8>,
    pub cigar: Vec<u32>,
    pub seq: Vec<u8>,
    pub qual: Vec<u8>,
    pub aux: Vec<u8>,
}

unsafe fn rec_of(b: *const h::bam1_t) -> Rec {
    unsafe {
        let l = h::ow_bam_l_qseq(b) as usize;
        let seq_p = h::ow_bam_seq(b);
        let seq = (0..l)
            .map(|i| {
                let byte = *seq_p.add(i / 2);
                let code = if i % 2 == 0 { byte >> 4 } else { byte & 0xf };
                b"=ACMGRSVTWYHKDBN"[code as usize]
            })
            .collect();
        let n = h::ow_bam_n_cigar(b) as usize;
        Rec {
            tid: h::ow_bam_tid(b),
            pos: h::ow_bam_pos(b),
            end: h::bam_endpos(b),
            flag: h::ow_bam_flag(b) as u16,
            mapq: h::ow_bam_mapq(b) as u8,
            mtid: h::ow_bam_mtid(b),
            mpos: h::ow_bam_mpos(b),
            isize: h::ow_bam_isize(b),
            name: CStr::from_ptr(h::ow_bam_qname(b)).to_bytes().to_vec(),
            cigar: std::slice::from_raw_parts(h::ow_bam_cigar(b), n).to_vec(),
            seq,
            qual: std::slice::from_raw_parts(h::ow_bam_qual(b), l).to_vec(),
            aux: std::slice::from_raw_parts(h::ow_bam_aux(b), h::ow_bam_l_aux(b) as usize).to_vec(),
        }
    }
}

/// An open CRAM with its index, queried through htslib's iterator.
pub struct Hts {
    fp: *mut h::htsFile,
    hdr: *mut h::sam_hdr_t,
    idx: *mut h::hts_idx_t,
}

impl Hts {
    pub fn open(cram: &Path, crai: &Path, reference: Option<&Path>) -> Hts {
        unsafe {
            let fp = open_read(cram, reference);
            let hdr = h::sam_hdr_read(fp);
            assert!(!hdr.is_null());
            let idx =
                h::ow_sam_index_load(fp, path_cstring(cram).as_ptr(), path_cstring(crai).as_ptr());
            assert!(!idx.is_null());
            Hts { fp, hdr, idx }
        }
    }

    pub fn tid(&self, name: &str) -> i32 {
        let n = CString::new(name).unwrap();
        unsafe { h::sam_hdr_name2tid(self.hdr, n.as_ptr()) }
    }

    /// The records of `[beg, end)` of reference `tid`.
    pub fn query(&self, tid: i32, beg: i64, end: i64) -> Vec<Rec> {
        let mut out = Vec::new();
        unsafe {
            let itr = h::ow_sam_itr_queryi(self.idx, tid, beg, end);
            assert!(!itr.is_null());
            let b = h::bam_init1();
            loop {
                let r = h::ow_sam_itr_next(self.fp, itr, b);
                if r < 0 {
                    assert_eq!(r, -1, "htslib iterator failed");
                    break;
                }
                out.push(rec_of(b));
            }
            h::bam_destroy1(b);
            h::ow_hts_itr_destroy(itr);
        }
        out
    }
}

impl Drop for Hts {
    fn drop(&mut self) {
        unsafe {
            h::ow_hts_idx_destroy(self.idx);
            h::sam_hdr_destroy(self.hdr);
            h::hts_close(self.fp);
        }
    }
}

/// Every record of a file in order, as htslib reads it.
pub fn read_all(path: &Path, reference: Option<&Path>) -> Vec<Rec> {
    let mut out = Vec::new();
    unsafe {
        let fp = open_read(path, reference);
        let hdr = h::sam_hdr_read(fp);
        let b = h::bam_init1();
        loop {
            let status = h::sam_read1(fp, hdr, b);
            if status < 0 {
                assert_eq!(status, -1, "htslib read {path:?}");
                break;
            }
            out.push(rec_of(b));
        }
        h::bam_destroy1(b);
        h::sam_hdr_destroy(hdr);
        h::hts_close(fp);
    }
    out
}
