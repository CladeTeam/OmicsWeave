//! lib2bit as the reference side for 2bit, built on the raw bindings of
//! `lib2bit-sys`.

use std::ffi::{CString, c_void};
use std::path::Path;

use lib2bit_sys as ffi;

unsafe extern "C" {
    fn free(p: *mut c_void);
}

pub struct TwoBit(*mut ffi::TwoBit);

impl TwoBit {
    /// With `mask`, soft-masked bases come back lowercase.
    pub fn open(path: &Path, mask: bool) -> Option<TwoBit> {
        let path = CString::new(path.as_os_str().as_encoded_bytes()).unwrap();
        // SAFETY: lib2bit reads the path without keeping it.
        let tb = unsafe { ffi::twobitOpen(path.as_ptr().cast_mut(), mask as i32) };
        (!tb.is_null()).then_some(TwoBit(tb))
    }

    pub fn len(&self, name: &str) -> u32 {
        let name = CString::new(name).unwrap();
        // SAFETY: open handle; the name is only read.
        unsafe { ffi::twobitChromLen(self.0, name.as_ptr().cast_mut()) }
    }

    /// Bases `[start, end)`, `end > start`.
    pub fn sequence(&self, name: &str, start: u32, end: u32) -> Option<Vec<u8>> {
        let name = CString::new(name).unwrap();
        // SAFETY: the result is NULL or a malloc'd string of `end - start`
        // bases, freed after copying.
        unsafe {
            let p = ffi::twobitSequence(self.0, name.as_ptr().cast_mut(), start, end);
            if p.is_null() {
                return None;
            }
            let out = std::slice::from_raw_parts(p.cast::<u8>(), (end - start) as usize).to_vec();
            free(p.cast());
            Some(out)
        }
    }
}

impl Drop for TwoBit {
    fn drop(&mut self) {
        // SAFETY: opened by twobitOpen, closed once.
        unsafe { ffi::twobitClose(self.0) };
    }
}
