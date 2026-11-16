//! htslib encoding, decoding and indexing of CRAM input.

use htslib_sys as h;
use std::ffi::CString;
use std::path::Path;

pub(super) fn c(p: &Path) -> CString {
    CString::new(p.as_os_str().as_encoded_bytes()).unwrap()
}

/// Opens `path` for reading with htslib, with `reference` for CRAM.
unsafe fn open_read(path: &Path, reference: Option<&Path>) -> *mut h::htsFile {
    unsafe {
        let fp = h::hts_open(c(path).as_ptr(), c"r".as_ptr());
        assert!(!fp.is_null(), "open {path:?}");
        if let Some(r) = reference {
            assert_eq!(h::hts_set_fai_filename(fp, c(r).as_ptr()), 0);
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
        let w = h::hts_open(c(out).as_ptr(), c"wc".as_ptr());
        assert!(!w.is_null());
        if let Some(reference) = reference {
            assert_eq!(h::hts_set_fai_filename(w, c(reference).as_ptr()), 0);
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
        while h::sam_read1(r, hdr, rec) >= 0 {
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
            h::sam_index_build3(c(file).as_ptr(), c(out).as_ptr(), 0, 0),
            0
        );
    }
}

/// The BAM htslib writes for the CRAM at `src` (`samtools view -b`), at
/// compression `level`.
pub fn bam_of_cram(src: &Path, reference: Option<&Path>, out: &Path, level: u32) {
    unsafe {
        let r = open_read(src, reference);
        let hdr = h::sam_hdr_read(r);
        assert!(!hdr.is_null());
        let mode = CString::new(format!("wb{level}")).unwrap();
        let w = h::hts_open(c(out).as_ptr(), mode.as_ptr());
        assert_eq!(h::sam_hdr_write(w, hdr), 0);
        let rec = h::bam_init1();
        loop {
            let ret = h::sam_read1(r, hdr, rec);
            if ret < 0 {
                assert_eq!(ret, -1, "htslib failed to read the CRAM");
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
