//! htslib-produced BGZF BCF and CSI, independent of the Rust index builder.
use std::ffi::{CString, c_int};
use std::path::PathBuf;

fn c(text: &str) -> CString {
    CString::new(text).unwrap()
}

pub struct Fixture {
    pub _dir: tempfile::TempDir,
    pub bcf: PathBuf,
    pub csi: PathBuf,
}

/// The VCF written as BGZF BCF and CSI-indexed with `min_shift` by htslib.
pub fn fixture(text: &str, min_shift: c_int) -> Fixture {
    let dir = tempfile::tempdir().unwrap();
    let vcf = dir.path().join("in.vcf");
    std::fs::write(&vcf, text).unwrap();
    let bcf = dir.path().join("in.bcf");
    // SAFETY: handles come from htslib and are checked before use.
    unsafe {
        let fp = htslib_sys::hts_open(c(vcf.to_str().unwrap()).as_ptr(), c("r").as_ptr());
        assert!(!fp.is_null());
        let hdr = htslib_sys::bcf_hdr_read(fp);
        assert!(!hdr.is_null());
        let out = htslib_sys::hts_open(c(bcf.to_str().unwrap()).as_ptr(), c("wb").as_ptr());
        assert!(!out.is_null());
        assert_eq!(htslib_sys::bcf_hdr_write(out, hdr), 0);
        let rec = htslib_sys::bcf_init();
        while htslib_sys::bcf_read(fp, hdr, rec) == 0 {
            assert_eq!(htslib_sys::bcf_write(out, hdr, rec), 0);
        }
        htslib_sys::bcf_destroy(rec);
        htslib_sys::hts_close(out);
        htslib_sys::bcf_hdr_destroy(hdr);
        htslib_sys::hts_close(fp);
        let built = htslib_sys::bcf_index_build3(
            c(bcf.to_str().unwrap()).as_ptr(),
            std::ptr::null(),
            min_shift,
            1,
        );
        assert_eq!(built, 0);
    }
    let csi = dir.path().join("in.bcf.csi");
    assert!(csi.is_file());
    Fixture {
        _dir: dir,
        bcf,
        csi,
    }
}
