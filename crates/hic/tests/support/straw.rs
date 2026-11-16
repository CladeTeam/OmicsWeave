//! External reference reader for standard v8/v9 layouts.
use omicsweave_core::ContactUnit;
use std::ffi::CString;
use std::path::Path;
use std::sync::Mutex;

static REFERENCE: Mutex<()> = Mutex::new(());

/// straw's raw observed records (`observed`, `NONE`) for
/// `chrom1:start1:end1` x `chrom2:start2:end2` (inclusive bin-start bounds
/// in `unit`): (x, y, count), x on the lower-indexed chromosome.
pub fn straw(
    path: &Path,
    loc1: &str,
    loc2: &str,
    unit: ContactUnit,
    res: u32,
) -> Vec<(u32, u32, f32)> {
    // The C++ reference can carry process-global state; serialize calls within this test binary.
    let _reference = REFERENCE.lock().unwrap();
    let c = |s: &str| CString::new(s).unwrap();
    let (t, n, f, l1, l2, u) = (
        c("observed"),
        c("NONE"),
        c(path.to_str().unwrap()),
        c(loc1),
        c(loc2),
        c(unit.as_str()),
    );
    // SAFETY: valid C strings; the result is released after copying.
    unsafe {
        let r = straw_sys::straw_query(
            t.as_ptr(),
            n.as_ptr(),
            f.as_ptr(),
            l1.as_ptr(),
            l2.as_ptr(),
            u.as_ptr(),
            res as i32,
        );
        assert!(!r.is_null(), "straw failed");
        let n = (*r).n as usize;
        let out = (0..n)
            .map(|i| {
                (
                    *(*r).bin_x.add(i) as u32,
                    *(*r).bin_y.add(i) as u32,
                    *(*r).counts.add(i),
                )
            })
            .collect();
        straw_sys::straw_free(r);
        out
    }
}
