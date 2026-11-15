//! Shared helpers of the BAM tests: random sorted SAM, and htslib as the
//! reference for writing, indexing and region queries.
#![allow(dead_code)]

use std::collections::BTreeMap;
use std::ffi::{CStr, CString};
use std::path::Path;

use htslib_sys as h;
use htslib_sys as hb;

/// A small deterministic generator (splitmix64).
pub struct Rng(u64);

impl Rng {
    pub fn new(seed: u64) -> Rng {
        Rng(seed)
    }

    pub fn next(&mut self) -> u64 {
        self.0 = self.0.wrapping_add(0x9E37_79B9_7F4A_7C15);
        let mut z = self.0;
        z = (z ^ (z >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
        z = (z ^ (z >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
        z ^ (z >> 31)
    }

    /// Uniform in `lo..hi`.
    pub fn range(&mut self, lo: u64, hi: u64) -> u64 {
        lo + self.next() % (hi - lo).max(1)
    }

    pub fn chance(&mut self, percent: u64) -> bool {
        self.next() % 100 < percent
    }
}

pub struct SamSpec {
    pub refs: Vec<(&'static str, u64)>,
    pub per_ref: usize,
    /// Mean gap between read starts.
    pub gap: u64,
    pub long_cigars: bool,
}

fn seq(rng: &mut Rng, n: usize) -> String {
    let bases = b"ACGTN";
    (0..n)
        .map(|_| {
            let hi = if rng.chance(2) { 5 } else { 4 };
            bases[rng.range(0, hi) as usize] as char
        })
        .collect()
}

fn qual(rng: &mut Rng, n: usize) -> String {
    (0..n)
        .map(|_| (b'!' + rng.range(0, 41) as u8) as char)
        .collect()
}

/// A CIGAR of query length `l` mixing every operation.
fn cigar(rng: &mut Rng, l: usize) -> String {
    let mut out = String::new();
    let mut left = l;
    if rng.chance(20) && left > 10 {
        let s = rng.range(1, 6) as usize;
        out.push_str(&format!("{s}S"));
        left -= s;
    }
    if rng.chance(5) {
        out.insert_str(0, &format!("{}H", rng.range(1, 30)));
    }
    let mut first = true;
    while left > 0 {
        let m = (rng.range(1, 60) as usize).min(left);
        let op = if rng.chance(30) {
            "="
        } else if rng.chance(10) {
            "X"
        } else {
            "M"
        };
        out.push_str(&format!("{m}{op}"));
        left -= m;
        first = false;
        if left == 0 {
            break;
        }
        match rng.range(0, 10) {
            0 => out.push_str(&format!("{}D", rng.range(1, 40))),
            1 => out.push_str(&format!("{}N", rng.range(50, 20_000))),
            2 if left > 2 => {
                let i = rng.range(1, 3) as usize;
                out.push_str(&format!("{i}I"));
                left -= i;
            }
            3 => out.push_str(&format!("{}P", rng.range(1, 4))),
            _ => {}
        }
    }
    let _ = first;
    out
}

fn rlen(cigar: &str) -> u64 {
    let mut n = 0u64;
    let mut num = 0u64;
    for c in cigar.chars() {
        if let Some(d) = c.to_digit(10) {
            num = num * 10 + d as u64;
        } else {
            if matches!(c, 'M' | 'D' | 'N' | '=' | 'X') {
                n += num;
            }
            num = 0;
        }
    }
    n
}

/// Random coordinate-sorted SAM text over `spec.refs`.
pub fn random_sam(rng: &mut Rng, spec: &SamSpec) -> Vec<u8> {
    let mut s = String::from("@HD\tVN:1.6\tSO:coordinate\n");
    for (name, len) in &spec.refs {
        s.push_str(&format!("@SQ\tSN:{name}\tLN:{len}\n"));
    }
    s.push_str("@RG\tID:rg1\tSM:x\n");
    let mut serial = 0;
    for (name, len) in &spec.refs {
        // Starts spread over the whole reference, some in dense piles.
        let mut starts = Vec::new();
        while starts.len() < spec.per_ref {
            let room = (*len / 4).min(40_000);
            let p = rng.range(1, len - room);
            starts.push(p);
            if rng.chance(10) {
                for _ in 0..rng.range(5, 40) {
                    starts.push(p + rng.range(0, spec.gap));
                }
            }
        }
        starts.sort_unstable();
        for pos in starts {
            serial += 1;
            let l = rng.range(20, 151) as usize;
            let tags = format!("RG:Z:rg1\tNM:i:{}\tXS:A:+", rng.range(0, 30));
            let (q, s_) = (
                if rng.chance(5) {
                    "*".to_string()
                } else {
                    qual(rng, l)
                },
                seq(rng, l),
            );
            if rng.chance(3) {
                // Unmapped, placed beside its mate.
                s.push_str(&format!(
                    "u{serial}\t{}\t{name}\t{pos}\t0\t*\t=\t{pos}\t0\t{s_}\t{q}\t{tags}\n",
                    4 | 1 | 0x80
                ));
                continue;
            }
            let c = if spec.long_cigars && rng.chance(1) {
                // More operations than a BAM CIGAR holds: htslib moves them
                // into a CG tag.
                let n = 33_000;
                let c: String = (0..n).map(|_| "1M1I").collect();
                let seq_l = 2 * n;
                s.push_str(&format!(
                    "long{serial}\t0\t{name}\t{pos}\t60\t{c}\t*\t0\t0\t{}\t*\n",
                    seq(rng, seq_l)
                ));
                continue;
            } else {
                cigar(rng, l)
            };
            let span = rlen(&c).max(1);
            if pos + span >= *len {
                continue;
            }
            let paired = rng.chance(60);
            let reverse = rng.chance(50);
            let mut flag = if reverse { 0x10 } else { 0 };
            for (bit, p) in [(0x100, 3), (0x200, 2), (0x400, 8), (0x800, 2)] {
                if rng.chance(p) {
                    flag |= bit;
                }
            }
            let (rnext, pnext, tlen) = if paired {
                flag |= 1 | if rng.chance(85) { 2 } else { 0 };
                flag |= if rng.chance(50) { 0x40 } else { 0x80 };
                if rng.chance(5) {
                    flag |= 8;
                    ("=".to_string(), pos, 0i64)
                } else {
                    let d = rng.range(0, 600);
                    if rng.chance(50) {
                        let mate = pos + d;
                        ("=".to_string(), mate, (mate + 100 - pos) as i64)
                    } else {
                        let mate = pos.saturating_sub(d).max(1);
                        ("=".to_string(), mate, -((pos + span - mate) as i64))
                    }
                }
            } else {
                ("*".to_string(), 0, 0)
            };
            s.push_str(&format!(
                "r{serial}\t{flag}\t{name}\t{pos}\t{}\t{c}\t{rnext}\t{pnext}\t{tlen}\t{s_}\t{q}\t{tags}\n",
                rng.range(0, 61)
            ));
        }
    }
    for i in 0..rng.range(0, 30) {
        s.push_str(&format!("x{i}\t4\t*\t0\t0\t*\t*\t0\t0\tACGT\t*\n"));
    }
    s.into_bytes()
}

fn c(p: &Path) -> CString {
    CString::new(p.to_str().unwrap()).unwrap()
}

/// The BAM htslib writes for the SAM at `src`, compressed at `level`.
pub fn bam_of(src: &Path, out: &Path, level: u32) {
    unsafe {
        h::hts_set_log_level(h::HTS_LOG_OFF);
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
            h::hts_set_log_level(h::HTS_LOG_OFF);
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
        h::hts_set_log_level(h::HTS_LOG_OFF);
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

/// One reference of an index: bins (loff, chunks) and the linear index.
#[derive(Debug, PartialEq, Eq, Default)]
pub struct IdxRef {
    pub bins: BTreeMap<u32, (u64, Vec<(u64, u64)>)>,
    pub linear: Vec<u64>,
}

#[derive(Debug, PartialEq, Eq)]
pub struct Idx {
    pub head: Vec<u8>,
    pub refs: Vec<IdxRef>,
    pub n_no_coor: Option<u64>,
}

/// The text of a BGZF file.
pub fn inflate(raw: &[u8]) -> Vec<u8> {
    use std::io::Read;
    let mut out = Vec::new();
    omicsweave_bgzf::Reader::new(raw, Box::new(omicsweave_core::exec::Inline::new()))
        .read_to_end(&mut out)
        .unwrap();
    out
}

/// Parses BAI or CSI bytes (bins in any order).
pub fn parse_index(raw: &[u8]) -> Idx {
    let bytes = if raw.starts_with(b"BAI\x01") {
        raw.to_vec()
    } else {
        inflate(raw)
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
        let mut r = IdxRef::default();
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
    Idx {
        head,
        refs,
        n_no_coor,
    }
}

/// Splits a BGZF file into its blocks.
pub fn blocks(raw: &[u8]) -> Vec<&[u8]> {
    let mut out = Vec::new();
    let mut at = 0;
    while at < raw.len() {
        let len = u16::from_le_bytes([raw[at + 16], raw[at + 17]]) as usize + 1;
        out.push(&raw[at..at + len]);
        at += len;
    }
    out
}
