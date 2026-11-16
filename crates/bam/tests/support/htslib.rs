//! Independent htslib writing, indexing and record-query oracle.

use htslib_sys as h;
use htslib_sys as hb;
use std::ffi::{CStr, CString};
use std::path::Path;

fn c(p: &Path) -> CString {
    CString::new(p.to_str().unwrap()).unwrap()
}

/// The BAM htslib writes for the SAM at `src`, compressed at `level`.
pub fn bam_of(src: &Path, out: &Path, level: u32) {
    unsafe {
        let r = h::hts_open(c(src).as_ptr(), c"r".as_ptr());
        assert!(!r.is_null());
        let hdr = h::sam_hdr_read(r);
        assert!(!hdr.is_null());
        let mode = CString::new(format!("wb{level}")).unwrap();
        let w = h::hts_open(c(out).as_ptr(), mode.as_ptr());
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

/// htslib's BAI (`min_shift` 0) or CSI index of the BAM at `file`.
pub fn index(file: &Path, out: &Path, min_shift: i32) {
    unsafe {
        assert_eq!(
            h::sam_index_build3(c(file).as_ptr(), c(out).as_ptr(), min_shift, 0),
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
        let l = hb::ow_bam_l_qseq(b) as usize;
        let seq_p = hb::ow_bam_seq(b);
        let seq = (0..l)
            .map(|i| {
                let byte = *seq_p.add(i / 2);
                let code = if i % 2 == 0 { byte >> 4 } else { byte & 0xf };
                b"=ACMGRSVTWYHKDBN"[code as usize]
            })
            .collect();
        let n = hb::ow_bam_n_cigar(b) as usize;
        Rec {
            tid: hb::ow_bam_tid(b),
            pos: hb::ow_bam_pos(b),
            end: h::bam_endpos(b),
            flag: hb::ow_bam_flag(b) as u16,
            mapq: hb::ow_bam_mapq(b) as u8,
            mtid: hb::ow_bam_mtid(b),
            mpos: hb::ow_bam_mpos(b),
            isize: hb::ow_bam_isize(b),
            name: CStr::from_ptr(hb::ow_bam_qname(b)).to_bytes().to_vec(),
            cigar: std::slice::from_raw_parts(hb::ow_bam_cigar(b), n).to_vec(),
            seq,
            qual: std::slice::from_raw_parts(hb::ow_bam_qual(b), l).to_vec(),
            aux: std::slice::from_raw_parts(hb::ow_bam_aux(b), hb::ow_bam_l_aux(b) as usize)
                .to_vec(),
        }
    }
}

/// An open BAM with its index, queried through htslib's iterator.
pub struct Hts {
    fp: *mut h::htsFile,
    hdr: *mut h::sam_hdr_t,
    idx: *mut hb::hts_idx_t,
}

impl Hts {
    pub fn open(bam: &Path, idx: &Path) -> Hts {
        unsafe {
            let fp = h::hts_open(c(bam).as_ptr(), c"r".as_ptr());
            assert!(!fp.is_null());
            let hdr = h::sam_hdr_read(fp);
            assert!(!hdr.is_null());
            let idx = hb::ow_sam_index_load(fp, c(bam).as_ptr(), c(idx).as_ptr());
            assert!(!idx.is_null());
            Hts { fp, hdr, idx }
        }
    }

    pub fn tid(&self, name: &str) -> i32 {
        let n = CString::new(name).unwrap();
        unsafe { hb::sam_hdr_name2tid(self.hdr, n.as_ptr()) }
    }

    /// The records of `[beg, end)` of reference `tid`.
    pub fn query(&self, tid: i32, beg: i64, end: i64) -> Vec<Rec> {
        let mut out = Vec::new();
        unsafe {
            let itr = hb::ow_sam_itr_queryi(self.idx, tid, beg, end);
            assert!(!itr.is_null());
            let b = h::bam_init1();
            loop {
                let r = hb::ow_sam_itr_next(self.fp, itr, b);
                if r < 0 {
                    assert_eq!(r, -1, "htslib iterator failed");
                    break;
                }
                out.push(rec_of(b));
            }
            h::bam_destroy1(b);
            hb::ow_hts_itr_destroy(itr);
        }
        out
    }
}

impl Drop for Hts {
    fn drop(&mut self) {
        unsafe {
            hb::ow_hts_idx_destroy(self.idx);
            h::sam_hdr_destroy(self.hdr);
            h::hts_close(self.fp);
        }
    }
}

/// Every record of a BAM in file order, as htslib reads it.
pub fn read_all(bam: &Path) -> Vec<Rec> {
    let mut out = Vec::new();
    unsafe {
        let fp = h::hts_open(c(bam).as_ptr(), c"r".as_ptr());
        assert!(!fp.is_null());
        let hdr = h::sam_hdr_read(fp);
        let b = h::bam_init1();
        while h::sam_read1(fp, hdr, b) >= 0 {
            out.push(rec_of(b));
        }
        h::bam_destroy1(b);
        h::sam_hdr_destroy(hdr);
        h::hts_close(fp);
    }
    out
}
