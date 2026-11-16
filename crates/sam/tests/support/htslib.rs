//! Independent C reference: htslib writes BAM, then its BGZF bytes are decoded.

use super::reference::Encoded;
use omicsweave_sam::BamSpan;
use std::ffi::{CString, c_void};
use std::path::Path;

fn le_i32(b: &[u8], at: usize) -> i32 {
    i32::from_le_bytes(b[at..at + 4].try_into().unwrap())
}

/// Reads all of a BGZF file.
fn inflate(path: &Path) -> Vec<u8> {
    let c = |s: &str| CString::new(s).unwrap();
    let mut out = Vec::new();
    // SAFETY: the handle is closed below; reads go into `buf`.
    unsafe {
        let fp = htslib_sys::bgzf_open(c(path.to_str().unwrap()).as_ptr(), c("r").as_ptr());
        assert!(!fp.is_null());
        let mut buf = vec![0u8; 1 << 16];
        loop {
            let n = htslib_sys::bgzf_read(fp, buf.as_mut_ptr() as *mut c_void, buf.len());
            assert!(n >= 0);
            if n == 0 {
                break;
            }
            out.extend_from_slice(&buf[..n as usize]);
        }
        htslib_sys::bgzf_close(fp);
    }
    out
}

/// htslib's encoding: `sam_read1` each record of the SAM text and
/// `sam_write1` it to an uncompressed BAM (`"wbu"`), then cut the inflated
/// stream at the sizes `sam_write1` reports. With `exact`, the stream must
/// hold nothing else (a failed write can leave partial bytes behind).
pub fn htslib(sam: &[u8], exact: bool) -> Encoded {
    let c = |s: &str| CString::new(s).unwrap();
    let dir = tempfile::tempdir().unwrap();
    let input = dir.path().join("in.sam");
    let output = dir.path().join("out.bam");
    std::fs::write(&input, sam).unwrap();
    let lines = sam.iter().filter(|&&b| b == b'\n').count() + 1;
    let mut written: Vec<Option<(usize, i64)>> = Vec::new();
    // SAFETY: handles come from htslib and are released below.
    unsafe {
        let fp = htslib_sys::hts_open(c(input.to_str().unwrap()).as_ptr(), c("r").as_ptr());
        assert!(!fp.is_null());
        let hdr = htslib_sys::sam_hdr_read(fp);
        if hdr.is_null() {
            htslib_sys::hts_close(fp);
            return Encoded {
                header: None,
                records: Vec::new(),
            };
        }
        let out = htslib_sys::hts_open(c(output.to_str().unwrap()).as_ptr(), c("wbu").as_ptr());
        assert!(!out.is_null());
        assert_eq!(htslib_sys::sam_hdr_write(out, hdr), 0);
        let b = htslib_sys::bam_init1();
        loop {
            let r = htslib_sys::sam_read1(fp, hdr, b);
            if r == -1 {
                break;
            }
            assert!(written.len() <= lines, "htslib keeps failing at the end");
            if r < -1 {
                written.push(None);
                continue;
            }
            let w = htslib_sys::sam_write1(out, hdr, b);
            written.push((w >= 0).then(|| (w as usize, htslib_sys::bam_endpos(b))));
        }
        htslib_sys::bam_destroy1(b);
        assert_eq!(htslib_sys::hts_close(out), 0);
        htslib_sys::sam_hdr_destroy(hdr);
        htslib_sys::hts_close(fp);
    }
    let stream = inflate(&output);
    // BAM header: magic, l_text, text, n_ref, then l_name, name, l_ref each.
    let mut at = 8 + le_i32(&stream, 4) as usize;
    let n_ref = le_i32(&stream, at);
    at += 4;
    for _ in 0..n_ref {
        at += 8 + le_i32(&stream, at) as usize;
    }
    let header = stream[..at].to_vec();
    let records = written
        .into_iter()
        .map(|w| {
            let (size, end) = w?;
            let rec = stream[at..at + size].to_vec();
            at += size;
            let flag = u16::from_le_bytes([rec[18], rec[19]]);
            let span = BamSpan {
                tid: le_i32(&rec, 4),
                beg: i64::from(le_i32(&rec, 8)),
                end,
                mapped: flag & 4 == 0,
            };
            Some((rec, span))
        })
        .collect();
    if exact {
        assert_eq!(at, stream.len(), "bytes after the last record");
    }
    Encoded {
        header: Some(header),
        records,
    }
}
