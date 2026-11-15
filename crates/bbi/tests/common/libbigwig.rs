//! libBigWig as the reference side of the differential tests: reads the same
//! files as the crate, and writes random bigWig inputs. Built on the raw
//! bindings of `libbigwig-sys`; owns and allocates freely, since it only runs
//! in tests.

use std::ffi::{CStr, CString, c_char, c_void};
use std::path::Path;

use libbigwig_sys as ffi;

unsafe extern "C" {
    // Releases the string `bbGetSQL` allocates.
    fn free(p: *mut c_void);
}

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

/// One bigBed entry as libBigWig reports it.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Entry {
    pub start: u32,
    pub end: u32,
    pub rest: Vec<u8>,
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

    /// `bwStatsFromFull` (`exact`) or `bwStats` over `bins` equal bins of
    /// `[start, end)`; `kind` is the `enum bwStatsType` value. `None` when
    /// libBigWig fails.
    pub fn stats(
        &mut self,
        chrom: &str,
        start: u32,
        end: u32,
        bins: u32,
        kind: i32,
        exact: bool,
    ) -> Option<Vec<f64>> {
        let chrom = c_str(chrom);
        // SAFETY: a non-NULL result holds `bins` doubles allocated with
        // malloc; it is copied and freed here.
        unsafe {
            let o = if exact {
                ffi::bwStatsFromFull(self.fp, chrom.as_ptr(), start, end, bins, kind)
            } else {
                ffi::bwStats(self.fp, chrom.as_ptr(), start, end, bins, kind)
            };
            if o.is_null() {
                return None;
            }
            let out = std::slice::from_raw_parts(o, bins as usize).to_vec();
            free(o.cast());
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

    /// `bbGetOverlappingEntries(.., withString = 1)`.
    pub fn entries(&mut self, chrom: &str, start: u32, end: u32) -> Option<Vec<Entry>> {
        let chrom = c_str(chrom);
        // SAFETY: the result holds `l` items in each array, strings
        // NUL-terminated; it is copied and destroyed here.
        unsafe {
            let o = ffi::bbGetOverlappingEntries(self.fp, chrom.as_ptr(), start, end, 1);
            if o.is_null() {
                return None;
            }
            let o_ref = &*o;
            let out = (0..o_ref.l as usize)
                .map(|i| Entry {
                    start: *o_ref.start.add(i),
                    end: *o_ref.end.add(i),
                    rest: CStr::from_ptr(*o_ref.str.add(i)).to_bytes().to_vec(),
                })
                .collect();
            ffi::bbDestroyOverlappingEntries(o);
            Some(out)
        }
    }

    /// `bbGetSQL`, without the NUL terminator.
    pub fn auto_sql(&mut self) -> Option<String> {
        // SAFETY: NULL or a malloc'd NUL-terminated string, freed after copying.
        unsafe {
            let p = ffi::bbGetSQL(self.fp);
            if p.is_null() {
                return None;
            }
            let s = CStr::from_ptr(p).to_string_lossy().into_owned();
            free(p.cast());
            Some(s)
        }
    }
}

impl Drop for Reader {
    fn drop(&mut self) {
        // SAFETY: `fp` was opened here and is closed once.
        unsafe { ffi::bwClose(self.fp) };
    }
}

/// A bigWig file being written; dropping it finalizes the file (index and
/// zoom levels).
pub struct Writer {
    fp: *mut ffi::bigWigFile_t,
    // libBigWig keeps pointers to the chromosome names it was created with.
    _names: Vec<CString>,
}

impl Writer {
    pub fn create(path: &Path, max_zooms: i32, chroms: &[(&str, u32)]) -> Option<Writer> {
        let path = c_path(path);
        let names: Vec<CString> = chroms.iter().map(|(n, _)| c_str(n)).collect();
        let ptrs: Vec<*const c_char> = names.iter().map(|n| n.as_ptr()).collect();
        let lens: Vec<u32> = chroms.iter().map(|&(_, l)| l).collect();
        // SAFETY: `ptrs` and `lens` hold `chroms.len()` elements; `names`
        // outlives the file. On failure the file is closed before returning.
        unsafe {
            let fp = ffi::bwOpen(path.as_ptr(), None, c"w".as_ptr());
            if fp.is_null() {
                return None;
            }
            let ok = ffi::bwCreateHdr(fp, max_zooms) == 0
                && {
                    (*fp).cl =
                        ffi::bwCreateChromList(ptrs.as_ptr(), lens.as_ptr(), chroms.len() as i64);
                    !(*fp).cl.is_null()
                }
                && ffi::bwWriteHdr(fp) == 0;
            if !ok {
                ffi::bwClose(fp);
                return None;
            }
            Some(Writer { fp, _names: names })
        }
    }

    /// One bedGraph block per call; items sorted, non-overlapping.
    pub fn add_bedgraph(&mut self, chrom: &str, items: &[(u32, u32, f32)]) -> bool {
        let name = c_str(chrom);
        let chroms = vec![name.as_ptr(); items.len()];
        let starts: Vec<u32> = items.iter().map(|i| i.0).collect();
        let ends: Vec<u32> = items.iter().map(|i| i.1).collect();
        let values: Vec<f32> = items.iter().map(|i| i.2).collect();
        // SAFETY: every array holds `items.len()` elements.
        unsafe {
            ffi::bwAddIntervals(
                self.fp,
                chroms.as_ptr(),
                starts.as_ptr(),
                ends.as_ptr(),
                values.as_ptr(),
                items.len() as u32,
            ) == 0
        }
    }

    pub fn add_varstep(&mut self, chrom: &str, span: u32, items: &[(u32, f32)]) -> bool {
        let name = c_str(chrom);
        let starts: Vec<u32> = items.iter().map(|i| i.0).collect();
        let values: Vec<f32> = items.iter().map(|i| i.1).collect();
        // SAFETY: both arrays hold `items.len()` elements.
        unsafe {
            ffi::bwAddIntervalSpans(
                self.fp,
                name.as_ptr(),
                starts.as_ptr(),
                span,
                values.as_ptr(),
                items.len() as u32,
            ) == 0
        }
    }

    pub fn add_fixedstep(
        &mut self,
        chrom: &str,
        start: u32,
        span: u32,
        step: u32,
        values: &[f32],
    ) -> bool {
        let name = c_str(chrom);
        // SAFETY: `values` holds the count passed.
        unsafe {
            ffi::bwAddIntervalSpanSteps(
                self.fp,
                name.as_ptr(),
                start,
                span,
                step,
                values.as_ptr(),
                values.len() as u32,
            ) == 0
        }
    }
}

impl Drop for Writer {
    fn drop(&mut self) {
        // SAFETY: `fp` was created here and is closed once.
        unsafe { ffi::bwClose(self.fp) };
    }
}
