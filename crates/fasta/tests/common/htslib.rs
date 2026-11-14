//! htslib as the reference side for FASTA, bgzip and their indexes, built on
//! the raw bindings of `htslib-sys`.

use std::ffi::{CStr, CString, c_void};
use std::path::Path;
use std::ptr;

use htslib_sys as ffi;

unsafe extern "C" {
    fn free(p: *mut c_void);
}

fn c_path(path: &Path) -> CString {
    CString::new(path.as_os_str().as_encoded_bytes()).unwrap()
}

/// Silences htslib's logging (expected failures would otherwise print).
pub fn quiet() {
    // SAFETY: sets a global log level.
    unsafe { ffi::hts_set_log_level(ffi::HTS_LOG_OFF) };
}

/// `fai_build3(fa, NULL, NULL)`: writes `fa.fai`, and `fa.gzi` for bgzip
/// input. False when htslib rejects the file.
pub fn build_index(fa: &Path) -> bool {
    let fa = c_path(fa);
    // SAFETY: valid path; NULL index names select the defaults.
    unsafe { ffi::fai_build3(fa.as_ptr(), ptr::null(), ptr::null()) == 0 }
}

/// An index loaded from the files next to `fa`.
pub struct Faidx(*mut ffi::faidx_t);

impl Faidx {
    pub fn load(fa: &Path) -> Option<Faidx> {
        let fa = c_path(fa);
        // SAFETY: valid path; flags 0 loads existing indexes only.
        let fai = unsafe { ffi::fai_load3(fa.as_ptr(), ptr::null(), ptr::null(), 0) };
        (!fai.is_null()).then_some(Faidx(fai))
    }

    pub fn names(&self) -> Vec<String> {
        // SAFETY: indices below faidx_nseq name sequences of the index.
        unsafe {
            (0..ffi::faidx_nseq(self.0))
                .map(|i| {
                    CStr::from_ptr(ffi::faidx_iseq(self.0, i))
                        .to_string_lossy()
                        .into_owned()
                })
                .collect()
        }
    }

    /// Bases `[start, end)`, `end > start`.
    pub fn fetch(&self, name: &str, start: u32, end: u32) -> Option<Vec<u8>> {
        let name = CString::new(name).unwrap();
        let mut len = 0;
        // SAFETY: the result is NULL or a malloc'd string of `len` bytes,
        // freed after copying.
        unsafe {
            let p = ffi::faidx_fetch_seq64(
                self.0,
                name.as_ptr(),
                start as i64,
                end as i64 - 1,
                &mut len,
            );
            if p.is_null() || len < 0 {
                return None;
            }
            let out = std::slice::from_raw_parts(p.cast::<u8>(), len as usize).to_vec();
            free(p.cast());
            Some(out)
        }
    }
}

impl Drop for Faidx {
    fn drop(&mut self) {
        // SAFETY: loaded by fai_load3, destroyed once.
        unsafe { ffi::fai_destroy(self.0) };
    }
}

/// Writes `text` as bgzip to `path` with `path.gzi`, as `bgzip -i` does.
pub fn bgzip(text: &[u8], path: &Path) {
    let p = c_path(path);
    // SAFETY: valid path and buffers; the handle is closed before returning.
    unsafe {
        let fp = ffi::bgzf_open(p.as_ptr(), c"w".as_ptr());
        assert!(!fp.is_null(), "bgzf_open");
        assert_eq!(ffi::bgzf_index_build_init(fp), 0);
        assert_eq!(
            ffi::bgzf_write(fp, text.as_ptr().cast(), text.len()),
            text.len() as isize
        );
        assert_eq!(ffi::bgzf_index_dump(fp, p.as_ptr(), c".gzi".as_ptr()), 0);
        assert_eq!(ffi::bgzf_close(fp), 0);
    }
}

/// Decompresses a whole BGZF file.
pub fn bgzf_read_all(path: &Path) -> Vec<u8> {
    let p = c_path(path);
    let mut out = Vec::new();
    let mut buf = vec![0u8; 1 << 16];
    // SAFETY: valid path and buffer; the handle is closed before returning.
    unsafe {
        let fp = ffi::bgzf_open(p.as_ptr(), c"r".as_ptr());
        assert!(!fp.is_null(), "bgzf_open");
        loop {
            let n = ffi::bgzf_read(fp, buf.as_mut_ptr().cast(), buf.len());
            assert!(n >= 0, "bgzf_read");
            if n == 0 {
                break;
            }
            out.extend_from_slice(&buf[..n as usize]);
        }
        assert_eq!(ffi::bgzf_close(fp), 0);
    }
    out
}

/// `bgzip -r`: reads the whole BGZF file at `path`, indexing it on the way,
/// and returns the `.gzi` bytes htslib writes; `None` when htslib fails.
/// Only meaningful for compressed input (`bgzip -r` refuses plain files).
pub fn reindex(path: &Path) -> Option<Vec<u8>> {
    let p = c_path(path);
    let mut name = path.as_os_str().to_owned();
    name.push(".reindex");
    let index = c_path(Path::new(&name));
    let mut buf = vec![0u8; 0xff00];
    // SAFETY: valid paths and buffer; the handle is closed before returning.
    let ok = unsafe {
        let fp = ffi::bgzf_open(p.as_ptr(), c"r".as_ptr());
        if fp.is_null() {
            return None;
        }
        let mut ok = ffi::bgzf_index_build_init(fp) == 0;
        while ok {
            let n = ffi::bgzf_read(fp, buf.as_mut_ptr().cast(), buf.len());
            if n <= 0 {
                ok = n == 0;
                break;
            }
        }
        ok = ok && ffi::bgzf_index_dump(fp, index.as_ptr(), c".gzi".as_ptr()) == 0;
        ffi::bgzf_close(fp) == 0 && ok
    };
    name.push(".gzi");
    let bytes = std::fs::read(&name).ok();
    let _ = std::fs::remove_file(&name);
    if ok { bytes } else { None }
}

/// Bytes of a bgzip file made by concatenating htslib's bgzip output of
/// each part, as `cat a.gz b.gz` does: every part ends in an empty
/// end-of-file member, and an empty part is that member alone.
pub fn bgzip_parts(parts: &[&[u8]], dir: &Path) -> Vec<u8> {
    let mut out = Vec::new();
    for (i, part) in parts.iter().enumerate() {
        let path = dir.join(format!("part{i}.gz"));
        bgzip(part, &path);
        out.extend(std::fs::read(&path).unwrap());
    }
    out
}
