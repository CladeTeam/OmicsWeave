//! htslib as the reference: what it writes for the same text, and the
//! indexes it builds for a file, parsed into a canonical form.

use std::collections::BTreeMap;
use std::ffi::CString;
use std::path::Path;

use htslib_sys as h;

fn c(p: &Path) -> CString {
    CString::new(p.to_str().unwrap()).unwrap()
}

/// Uncompressed BCF htslib writes from the VCF at `src`.
pub fn bcf_of(src: &Path, out: &Path) {
    unsafe {
        h::hts_set_log_level(h::HTS_LOG_OFF);
        let r = h::hts_open(c(src).as_ptr(), c"r".as_ptr());
        assert!(!r.is_null());
        let hdr = h::bcf_hdr_read(r);
        assert!(!hdr.is_null());
        let w = h::hts_open(c(out).as_ptr(), c"wbu".as_ptr());
        assert_eq!(h::bcf_hdr_write(w, hdr), 0);
        let rec = h::bcf_init();
        while h::bcf_read(r, hdr, rec) == 0 {
            assert_eq!(h::bcf_write(w, hdr, rec), 0);
        }
        h::bcf_destroy(rec);
        h::bcf_hdr_destroy(hdr);
        h::hts_close(r);
        assert_eq!(h::hts_close(w), 0);
    }
}

/// Uncompressed BAM htslib writes from the SAM at `src`.
pub fn bam_of(src: &Path, out: &Path) {
    unsafe {
        h::hts_set_log_level(h::HTS_LOG_OFF);
        let r = h::hts_open(c(src).as_ptr(), c"r".as_ptr());
        assert!(!r.is_null());
        let hdr = h::sam_hdr_read(r);
        assert!(!hdr.is_null());
        let w = h::hts_open(c(out).as_ptr(), c"wbu".as_ptr());
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

/// htslib's CSI index of the BCF at `file`, written to `out`.
pub fn bcf_index(file: &Path, out: &Path) {
    unsafe {
        assert_eq!(
            h::bcf_index_build3(c(file).as_ptr(), c(out).as_ptr(), 14, 0),
            0
        );
    }
}

/// htslib's BAI (`min_shift` 0) or CSI index of the BAM at `file`.
pub fn bam_index(file: &Path, out: &Path, min_shift: i32) {
    unsafe {
        assert_eq!(
            h::sam_index_build3(c(file).as_ptr(), c(out).as_ptr(), min_shift, 0),
            0
        );
    }
}

/// One reference of an index: bins (loff, chunks) and the linear index.
#[derive(Debug, PartialEq, Eq, Default)]
pub struct Ref {
    pub bins: BTreeMap<u32, (u64, Vec<(u64, u64)>)>,
    pub linear: Vec<u64>,
}

#[derive(Debug, PartialEq, Eq)]
pub struct Index {
    pub head: Vec<u8>,
    pub refs: Vec<Ref>,
    pub n_no_coor: Option<u64>,
}

/// Parses a BAI or CSI file (bins in any order).
pub fn parse_index(path: &Path) -> Index {
    let raw = std::fs::read(path).unwrap();
    let bytes = if raw.starts_with(b"BAI\x01") {
        raw
    } else {
        super::inflate(path)
    };
    let mut at = 0usize;
    let mut take = |n: usize| {
        let s = bytes.get(at..at + n).map(<[u8]>::to_vec);
        at += n;
        s.unwrap_or_default()
    };
    let u32_ = |b: Vec<u8>| u32::from_le_bytes(b.try_into().unwrap());
    let u64_ = |b: Vec<u8>| u64::from_le_bytes(b.try_into().unwrap());
    let magic = take(4);
    let csi = magic == b"CSI\x01";
    let mut head = magic.clone();
    if csi {
        head.extend(take(8));
        let l_aux = u32_(take(4)) as usize;
        head.extend(take(l_aux));
    }
    let n_ref = u32_(take(4));
    let mut refs = Vec::new();
    for _ in 0..n_ref {
        let mut r = Ref::default();
        for _ in 0..u32_(take(4)) {
            let bin = u32_(take(4));
            let loff = if csi { u64_(take(8)) } else { 0 };
            let n = u32_(take(4));
            let chunks = (0..n).map(|_| (u64_(take(8)), u64_(take(8)))).collect();
            r.bins.insert(bin, (loff, chunks));
        }
        if !csi {
            for _ in 0..u32_(take(4)) {
                r.linear.push(u64_(take(8)));
            }
        }
        refs.push(r);
    }
    let tail = take(8);
    let n_no_coor = (tail.len() == 8).then(|| u64_(tail));
    Index {
        head,
        refs,
        n_no_coor,
    }
}
