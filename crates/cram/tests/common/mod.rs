//! Shared helpers of the CRAM tests: random references and sorted SAM
//! drawn from them, and htslib as the reference for writing CRAM,
//! indexing it, region queries and conversion to BAM.
#![allow(dead_code)]

use std::ffi::{CStr, CString};
use std::path::{Path, PathBuf};

use htslib_sys as h;
use omicsweave_core::FileSource;
use omicsweave_fasta::{Fai, FastaReader};

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

fn c(p: &Path) -> CString {
    CString::new(p.to_str().unwrap()).unwrap()
}

/// Random reference sequences: mostly ACGT, with N runs and soft-masked
/// (lowercase) stretches.
pub fn reference(rng: &mut Rng, refs: &[(String, u64)]) -> Vec<Vec<u8>> {
    refs.iter()
        .map(|(_, len)| {
            let mut s = Vec::with_capacity(*len as usize);
            while (s.len() as u64) < *len {
                let run = rng.range(1, 400) as usize;
                let kind = rng.range(0, 100);
                for _ in 0..run {
                    let b = b"ACGT"[rng.range(0, 4) as usize];
                    s.push(match kind {
                        0..=1 => b'N',
                        2..=12 => b.to_ascii_lowercase(),
                        _ => b,
                    });
                }
            }
            s.truncate(*len as usize);
            s
        })
        .collect()
}

/// Writes FASTA (60 bases a line) and its `.fai` (htslib's).
pub fn write_fasta(path: &Path, refs: &[(String, u64)], seqs: &[Vec<u8>]) {
    let mut text = Vec::new();
    for ((name, _), s) in refs.iter().zip(seqs) {
        text.extend_from_slice(format!(">{name}\n").as_bytes());
        for line in s.chunks(60) {
            text.extend_from_slice(line);
            text.push(b'\n');
        }
    }
    std::fs::write(path, text).unwrap();
    unsafe {
        h::hts_set_log_level(h::HTS_LOG_OFF);
        assert_eq!(
            h::fai_build3(c(path).as_ptr(), std::ptr::null(), std::ptr::null()),
            0
        );
    }
}

/// Our reader of the FASTA written by [`write_fasta`].
pub fn fasta_reader(path: &Path) -> FastaReader<FileSource> {
    let mut fai = path.as_os_str().to_owned();
    fai.push(".fai");
    let fai = Fai::parse(&std::fs::read(PathBuf::from(fai)).unwrap()).unwrap();
    FastaReader::new(FileSource::open(path).unwrap(), fai, None)
}

pub struct SamSpec {
    pub per_ref: usize,
    /// Mean gap between read starts.
    pub gap: u64,
    pub long_cigars: bool,
    pub unplaced: usize,
}

fn rand_bases(rng: &mut Rng, n: usize) -> Vec<u8> {
    (0..n)
        .map(|_| {
            if rng.chance(1) {
                b'N'
            } else {
                b"ACGT"[rng.range(0, 4) as usize]
            }
        })
        .collect()
}

/// Qualities that wander along the read, each near the mean of the two
/// before it and drifting down, as sequencers' do.
fn qual(rng: &mut Rng, n: usize) -> String {
    let (mut q1, mut q2) = (38i64, 38i64);
    (0..n)
        .map(|i| {
            let q = if rng.chance(2) {
                2
            } else {
                ((q1 + q2) / 2 + rng.range(0, 5) as i64 - 2 - (i / 60) as i64).clamp(2, 41)
            };
            (q2, q1) = (q1, q);
            (b'!' + q as u8) as char
        })
        .collect()
}

/// A read drawn from `seq` at 0-based `pos`: its CIGAR, bases and the
/// reference length it spans.
fn read_at(rng: &mut Rng, seq: &[u8], pos: usize, l: usize) -> (String, Vec<u8>, usize) {
    let mut cigar = String::new();
    let mut bases = Vec::new();
    let mut at = pos;
    if rng.chance(5) {
        cigar.push_str(&format!("{}H", rng.range(1, 20)));
    }
    if rng.chance(15) {
        let s = rng.range(1, 8) as usize;
        cigar.push_str(&format!("{s}S"));
        bases.extend(rand_bases(rng, s));
    }
    let mut left = l.saturating_sub(bases.len()).max(1);
    while left > 0 {
        let m = (rng.range(1, 60) as usize).min(left);
        let op = if rng.chance(20) {
            '='
        } else if rng.chance(5) {
            'X'
        } else {
            'M'
        };
        cigar.push_str(&format!("{m}{op}"));
        for _ in 0..m {
            let r = seq.get(at).copied().unwrap_or(b'N').to_ascii_uppercase();
            bases.push(if rng.chance(3) {
                b"ACGTN"[rng.range(0, 5) as usize]
            } else {
                r
            });
            at += 1;
        }
        left -= m;
        if left == 0 {
            break;
        }
        match rng.range(0, 12) {
            0 => {
                let d = rng.range(1, 30) as usize;
                cigar.push_str(&format!("{d}D"));
                at += d;
            }
            1 => {
                let n = rng.range(50, 3_000) as usize;
                cigar.push_str(&format!("{n}N"));
                at += n;
            }
            2 if left > 3 => {
                let i = rng.range(1, 4) as usize;
                cigar.push_str(&format!("{i}I"));
                bases.extend(rand_bases(rng, i));
                left -= i;
            }
            3 => cigar.push_str(&format!("{}P", rng.range(1, 3))),
            _ => {}
        }
    }
    if rng.chance(10) {
        let s = rng.range(1, 6) as usize;
        cigar.push_str(&format!("{s}S"));
        bases.extend(rand_bases(rng, s));
    }
    if rng.chance(3) {
        cigar.push_str(&format!("{}H", rng.range(1, 20)));
    }
    (cigar, bases, at - pos)
}

fn tags(rng: &mut Rng) -> String {
    let mut t = vec![format!("RG:Z:rg{}", rng.range(1, 3))];
    if rng.chance(30) {
        t.push(format!("NM:i:{}", rng.range(0, 30)));
    }
    if rng.chance(10) {
        t.push("MD:Z:10A5^C3".to_string());
    }
    if rng.chance(50) {
        t.push(format!("XS:A:{}", if rng.chance(50) { '+' } else { '-' }));
    }
    if rng.chance(30) {
        t.push(format!("XI:i:{}", rng.range(0, 100_000) as i64 - 50_000));
    }
    if rng.chance(10) {
        t.push("XF:f:1.5".to_string());
    }
    if rng.chance(10) {
        t.push(format!("XB:B:s,{},-2,3", rng.range(0, 100)));
    }
    if rng.chance(5) {
        t.push("XH:H:1AE301".to_string());
    }
    if rng.chance(20) {
        t.push(format!("XZ:Z:{}", "abc".repeat(rng.range(0, 4) as usize)));
    }
    t.join("\t")
}

/// Random coordinate-sorted SAM over `refs`, reads drawn from `seqs`.
pub fn random_sam(
    rng: &mut Rng,
    refs: &[(String, u64)],
    seqs: &[Vec<u8>],
    spec: &SamSpec,
) -> Vec<u8> {
    let mut s = String::from("@HD\tVN:1.6\tSO:coordinate\n");
    for (name, len) in refs {
        s.push_str(&format!("@SQ\tSN:{name}\tLN:{len}\n"));
    }
    s.push_str("@RG\tID:rg1\tSM:x\n@RG\tID:rg2\tSM:y\n@PG\tID:gen\tPN:test\n@CO\tA comment\n");
    let mut serial = 0;
    for (tid, (name, len)) in refs.iter().enumerate() {
        let len = *len as usize;
        let seq = &seqs[tid];
        let mut lines: Vec<(usize, String)> = Vec::new();
        let mut starts = Vec::new();
        while starts.len() < spec.per_ref {
            let p = rng.range(0, (len as u64).saturating_sub(5_000).max(1)) as usize;
            starts.push(p);
            if rng.chance(10) {
                for _ in 0..rng.range(5, 40) {
                    starts
                        .push((p + rng.range(0, spec.gap) as usize).min(len.saturating_sub(5_000)));
                }
            }
        }
        for pos in starts {
            serial += 1;
            let l = rng.range(20, 151) as usize;
            if rng.chance(3) {
                // Unmapped, placed.
                let b = String::from_utf8(rand_bases(rng, l)).unwrap();
                let q = qual(rng, l);
                lines.push((
                    pos,
                    format!(
                        "u{serial}\t{}\t{name}\t{}\t0\t*\t=\t{}\t0\t{b}\t{q}\t{}",
                        4 | 1 | 0x80,
                        pos + 1,
                        pos + 1,
                        tags(rng)
                    ),
                ));
                continue;
            }
            if spec.long_cigars && rng.chance(1) {
                let n = 33_000;
                let c: String = (0..n).map(|_| "1M1I").collect();
                let b = String::from_utf8(rand_bases(rng, 2 * n)).unwrap();
                if pos + n + 10 < len {
                    lines.push((
                        pos,
                        format!(
                            "long{serial}\t0\t{name}\t{}\t60\t{c}\t*\t0\t0\t{b}\t*",
                            pos + 1
                        ),
                    ));
                }
                continue;
            }
            let (cigar, bases, span) = read_at(rng, seq, pos, l);
            let span = span.max(1);
            if pos + span >= len {
                continue;
            }
            let no_seq = rng.chance(2);
            let (b, q) = if no_seq {
                ("*".to_string(), "*".to_string())
            } else {
                let q = if rng.chance(5) {
                    "*".to_string()
                } else {
                    qual(rng, bases.len())
                };
                (String::from_utf8(bases.clone()).unwrap(), q)
            };
            let mapq = rng.range(0, 61);
            let reverse = rng.chance(50);
            let mut flag = if reverse { 0x10 } else { 0 };
            for (bit, p) in [(0x100, 3), (0x200, 2), (0x400, 8), (0x800, 2)] {
                if rng.chance(p) {
                    flag |= bit;
                }
            }
            if rng.chance(40) && !no_seq {
                // A proper pair: the mate starts downstream, same reference.
                let mpos = pos + rng.range(0, 400) as usize;
                let ml = rng.range(20, 151) as usize;
                let (mcigar, mbases, mspan) = read_at(rng, seq, mpos, ml);
                if mpos + mspan + 1 >= len {
                    continue;
                }
                let end = (pos + span).max(mpos + mspan);
                let tlen = (end - pos) as i64 * if rng.chance(90) { 1 } else { 2 };
                let mrev = rng.chance(50);
                let f1 = flag | 0x1 | 0x2 | 0x40 | if mrev { 0x20 } else { 0 };
                let f2 =
                    0x1 | 0x2 | 0x80 | if mrev { 0x10 } else { 0 } | if reverse { 0x20 } else { 0 };
                let mq = qual(rng, mbases.len());
                lines.push((
                    pos,
                    format!(
                        "p{serial}\t{f1}\t{name}\t{}\t{mapq}\t{cigar}\t=\t{}\t{tlen}\t{b}\t{q}\t{}",
                        pos + 1,
                        mpos + 1,
                        tags(rng)
                    ),
                ));
                lines.push((
                    mpos,
                    format!(
                        "p{serial}\t{f2}\t{name}\t{}\t{}\t{mcigar}\t=\t{}\t{}\t{}\t{mq}\t{}",
                        mpos + 1,
                        rng.range(0, 61),
                        pos + 1,
                        -tlen,
                        String::from_utf8(mbases).unwrap(),
                        tags(rng)
                    ),
                ));
                continue;
            }
            let (rnext, pnext, tlen) = if rng.chance(20) {
                flag |= 1 | if rng.chance(50) { 0x40 } else { 0x80 };
                if rng.chance(20) {
                    flag |= 8;
                }
                let other = if refs.len() > 1 && rng.chance(30) {
                    refs[(tid + 1) % refs.len()].0.clone()
                } else {
                    "=".to_string()
                };
                (
                    other,
                    rng.range(1, len as u64),
                    rng.range(0, 1000) as i64 - 500,
                )
            } else {
                ("*".to_string(), 0, 0)
            };
            lines.push((
                pos,
                format!(
                    "r{serial}\t{flag}\t{name}\t{}\t{mapq}\t{cigar}\t{rnext}\t{pnext}\t{tlen}\t{b}\t{q}\t{}",
                    pos + 1,
                    tags(rng)
                ),
            ));
        }
        lines.sort_by_key(|(p, _)| *p);
        for (_, l) in lines {
            s.push_str(&l);
            s.push('\n');
        }
    }
    for i in 0..spec.unplaced {
        s.push_str(&format!(
            "x{i}\t4\t*\t0\t0\t*\t*\t0\t0\tACGTN\tIIIII\tRG:Z:rg1\n"
        ));
    }
    s.into_bytes()
}

/// Opens `path` for reading with htslib, with `reference` for CRAM.
unsafe fn open_read(path: &Path, reference: Option<&Path>) -> *mut h::htsFile {
    unsafe {
        h::hts_set_log_level(h::HTS_LOG_OFF);
        let fp = h::hts_open(c(path).as_ptr(), c"r".as_ptr());
        assert!(!fp.is_null(), "open {path:?}");
        if let Some(r) = reference {
            assert_eq!(h::hts_set_fai_filename(fp, c(r).as_ptr()), 0);
        }
        fp
    }
}

/// The CRAM htslib writes for the SAM at `src` against `reference`, with
/// format options (`version=3.0`, `embed_ref=1`, ...).
pub fn cram_of(src: &Path, out: &Path, reference: Option<&Path>, options: &[&str]) {
    unsafe {
        let r = open_read(src, None);
        let hdr = h::sam_hdr_read(r);
        assert!(!hdr.is_null());
        let w = h::hts_open(c(out).as_ptr(), c"wc".as_ptr());
        assert!(!w.is_null());
        if let Some(reference) = reference {
            assert_eq!(h::hts_set_fai_filename(w, c(reference).as_ptr()), 0);
        }
        let mut opts: *mut h::hts_opt = std::ptr::null_mut();
        for o in options {
            let o = CString::new(*o).unwrap();
            assert_eq!(h::hts_opt_add(&mut opts, o.as_ptr()), 0, "{o:?}");
        }
        assert_eq!(h::hts_opt_apply(w, opts), 0);
        h::hts_opt_free(opts);
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

/// htslib's `.crai` of the CRAM at `file`, written to `out`.
pub fn crai_of(file: &Path, out: &Path) {
    unsafe {
        h::hts_set_log_level(h::HTS_LOG_OFF);
        assert_eq!(
            h::sam_index_build3(c(file).as_ptr(), c(out).as_ptr(), 0, 0),
            0
        );
    }
}

/// The BAM htslib writes for the CRAM at `src` (`samtools view -b`), at
/// compression `level`.
pub fn bam_of_cram(src: &Path, reference: Option<&Path>, out: &Path, level: u32) {
    unsafe {
        let r = open_read(src, reference);
        let hdr = h::sam_hdr_read(r);
        assert!(!hdr.is_null());
        let mode = CString::new(format!("wb{level}")).unwrap();
        let w = h::hts_open(c(out).as_ptr(), mode.as_ptr());
        assert_eq!(h::sam_hdr_write(w, hdr), 0);
        let rec = h::bam_init1();
        loop {
            let ret = h::sam_read1(r, hdr, rec);
            if ret < 0 {
                assert_eq!(ret, -1, "htslib failed to read the CRAM");
                break;
            }
            assert!(h::sam_write1(w, hdr, rec) >= 0);
        }
        h::bam_destroy1(rec);
        h::sam_hdr_destroy(hdr);
        h::hts_close(r);
        assert_eq!(h::hts_close(w), 0);
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
        let l = h::ow_bam_l_qseq(b) as usize;
        let seq_p = h::ow_bam_seq(b);
        let seq = (0..l)
            .map(|i| {
                let byte = *seq_p.add(i / 2);
                let code = if i % 2 == 0 { byte >> 4 } else { byte & 0xf };
                b"=ACMGRSVTWYHKDBN"[code as usize]
            })
            .collect();
        let n = h::ow_bam_n_cigar(b) as usize;
        Rec {
            tid: h::ow_bam_tid(b),
            pos: h::ow_bam_pos(b),
            end: h::bam_endpos(b),
            flag: h::ow_bam_flag(b) as u16,
            mapq: h::ow_bam_mapq(b) as u8,
            mtid: h::ow_bam_mtid(b),
            mpos: h::ow_bam_mpos(b),
            isize: h::ow_bam_isize(b),
            name: CStr::from_ptr(h::ow_bam_qname(b)).to_bytes().to_vec(),
            cigar: std::slice::from_raw_parts(h::ow_bam_cigar(b), n).to_vec(),
            seq,
            qual: std::slice::from_raw_parts(h::ow_bam_qual(b), l).to_vec(),
            aux: std::slice::from_raw_parts(h::ow_bam_aux(b), h::ow_bam_l_aux(b) as usize).to_vec(),
        }
    }
}

/// An open CRAM with its index, queried through htslib's iterator.
pub struct Hts {
    fp: *mut h::htsFile,
    hdr: *mut h::sam_hdr_t,
    idx: *mut h::hts_idx_t,
}

impl Hts {
    pub fn open(cram: &Path, crai: &Path, reference: Option<&Path>) -> Hts {
        unsafe {
            let fp = open_read(cram, reference);
            let hdr = h::sam_hdr_read(fp);
            assert!(!hdr.is_null());
            let idx = h::ow_sam_index_load(fp, c(cram).as_ptr(), c(crai).as_ptr());
            assert!(!idx.is_null());
            Hts { fp, hdr, idx }
        }
    }

    pub fn tid(&self, name: &str) -> i32 {
        let n = CString::new(name).unwrap();
        unsafe { h::sam_hdr_name2tid(self.hdr, n.as_ptr()) }
    }

    /// The records of `[beg, end)` of reference `tid`.
    pub fn query(&self, tid: i32, beg: i64, end: i64) -> Vec<Rec> {
        let mut out = Vec::new();
        unsafe {
            let itr = h::ow_sam_itr_queryi(self.idx, tid, beg, end);
            assert!(!itr.is_null());
            let b = h::bam_init1();
            loop {
                let r = h::ow_sam_itr_next(self.fp, itr, b);
                if r < 0 {
                    assert_eq!(r, -1, "htslib iterator failed");
                    break;
                }
                out.push(rec_of(b));
            }
            h::bam_destroy1(b);
            h::ow_hts_itr_destroy(itr);
        }
        out
    }
}

impl Drop for Hts {
    fn drop(&mut self) {
        unsafe {
            h::ow_hts_idx_destroy(self.idx);
            h::sam_hdr_destroy(self.hdr);
            h::hts_close(self.fp);
        }
    }
}

/// Every record of a file in order, as htslib reads it.
pub fn read_all(path: &Path, reference: Option<&Path>) -> Vec<Rec> {
    let mut out = Vec::new();
    unsafe {
        let fp = open_read(path, reference);
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

/// The text of a gzip file (one member or several).
pub fn gunzip(raw: &[u8]) -> Vec<u8> {
    use std::io::Read;
    let mut out = Vec::new();
    flate2::read::MultiGzDecoder::new(raw)
        .read_to_end(&mut out)
        .unwrap();
    out
}
