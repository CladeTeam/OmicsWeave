//! libBigWig as the reference side of the differential tests: reads the same
//! files as the crate. Built on the raw
//! bindings of `libbigwig-sys`; owns and allocates freely, since it only runs
//! in tests.

use std::ffi::{CStr, CString};
use std::path::Path;

use libbigwig_sys as ffi;

fn c_path(path: &Path) -> CString {
    CString::new(path.as_os_str().as_encoded_bytes()).expect("path without NUL")
}

fn c_str(s: &str) -> CString {
    CString::new(s).expect("name without NUL")
}

/// The header fields libBigWig parsed, plus the file kind.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Header {
    /// 0 bigWig, 1 bigBed.
    pub kind: i32,
    pub version: u16,
    pub zoom_levels: u16,
    pub field_count: u16,
    pub defined_field_count: u16,
    pub buf_size: u32,
    pub chrom_tree: u64,
    pub full_data: u64,
    pub full_index: u64,
    pub auto_sql: u64,
    pub total_summary: u64,
    pub extension: u64,
    pub bases_covered: u64,
    pub min: f64,
    pub max: f64,
    pub sum: f64,
    pub sum_squares: f64,
}

/// A bigWig or bigBed file opened for reading.
pub struct Reader {
    fp: *mut ffi::bigWigFile_t,
}

impl Reader {
    pub fn open(path: &Path) -> Option<Reader> {
        let path = c_path(path);
        // SAFETY: `path` is NUL-terminated; no curl callback is used.
        let fp = unsafe {
            if ffi::bbIsBigBed(path.as_ptr(), None) != 0 {
                ffi::bbOpen(path.as_ptr(), None)
            } else {
                ffi::bwOpen(path.as_ptr(), None, c"r".as_ptr())
            }
        };
        (!fp.is_null()).then_some(Reader { fp })
    }

    fn file(&self) -> &ffi::bigWigFile_t {
        // SAFETY: `fp` is a file opened for reading, live until drop.
        unsafe { &*self.fp }
    }

    fn hdr(&self) -> &ffi::bigWigHdr_t {
        // SAFETY: an opened file always has a parsed header.
        unsafe { &*self.file().hdr }
    }

    pub fn header(&self) -> Header {
        let h = self.hdr();
        Header {
            kind: self.file().r#type,
            version: h.version,
            zoom_levels: h.nLevels,
            field_count: h.fieldCount,
            defined_field_count: h.definedFieldCount,
            buf_size: h.bufSize,
            chrom_tree: h.ctOffset,
            full_data: h.dataOffset,
            full_index: h.indexOffset,
            auto_sql: h.sqlOffset,
            total_summary: h.summaryOffset,
            extension: h.extensionOffset,
            bases_covered: h.nBasesCovered,
            min: h.minVal,
            max: h.maxVal,
            sum: h.sumData,
            sum_squares: h.sumSquared,
        }
    }

    /// (reduction, data offset, index offset) of each zoom level.
    pub fn zooms(&self) -> Vec<(u32, u64, u64)> {
        let h = self.hdr();
        if h.zoomHdrs.is_null() {
            return Vec::new();
        }
        // SAFETY: libBigWig allocates `nLevels` entries in each array.
        unsafe {
            let z = &*h.zoomHdrs;
            (0..h.nLevels as usize)
                .map(|i| (*z.level.add(i), *z.dataOffset.add(i), *z.indexOffset.add(i)))
                .collect()
        }
    }

    /// (name, length) in chromosome id order.
    pub fn chroms(&self) -> Vec<(String, u32)> {
        // SAFETY: the chromosome list holds `nKeys` names and lengths.
        unsafe {
            let cl = &*self.file().cl;
            (0..cl.nKeys as usize)
                .map(|i| {
                    let name = CStr::from_ptr(*cl.chrom.add(i))
                        .to_string_lossy()
                        .into_owned();
                    (name, *cl.len.add(i))
                })
                .collect()
        }
    }

    /// `bwGetValues(.., includeNA = 1)`: one value per base, NaN where
    /// uncovered. `None` when libBigWig fails.
    pub fn values(&mut self, chrom: &str, start: u32, end: u32) -> Option<Vec<f32>> {
        let chrom = c_str(chrom);
        // SAFETY: with includeNA the result holds `end - start` values; it is
        // copied and destroyed here.
        unsafe {
            let o = ffi::bwGetValues(self.fp, chrom.as_ptr(), start, end, 1);
            if o.is_null() {
                return None;
            }
            let out = std::slice::from_raw_parts((*o).value, (*o).l as usize).to_vec();
            ffi::bwDestroyOverlappingIntervals(o);
            Some(out)
        }
    }

    /// `bwGetOverlappingIntervals`: unclipped (start, end, value) items.
    pub fn intervals(&mut self, chrom: &str, start: u32, end: u32) -> Option<Vec<(u32, u32, f32)>> {
        let chrom = c_str(chrom);
        // SAFETY: the result holds `l` items in each array; it is copied and
        // destroyed here.
        unsafe {
            let o = ffi::bwGetOverlappingIntervals(self.fp, chrom.as_ptr(), start, end);
            if o.is_null() {
                return None;
            }
            let o_ref = &*o;
            let out = (0..o_ref.l as usize)
                .map(|i| (*o_ref.start.add(i), *o_ref.end.add(i), *o_ref.value.add(i)))
                .collect();
            ffi::bwDestroyOverlappingIntervals(o);
            Some(out)
        }
    }
}

impl Drop for Reader {
    fn drop(&mut self) {
        // SAFETY: `fp` was opened here and is closed once.
        unsafe { ffi::bwClose(self.fp) };
    }
}
