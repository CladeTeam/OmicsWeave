//! Round-four full-field workloads. The exact same source is usable with the
//! frozen baseline. htslib validates raw typed fields and all seven getters
//! before timing; no known baseline-error fixtures are used for comparisons.

#[path = "../tests/common/mod.rs"]
mod common;
#[path = "../tests/round7_fixture/mod.rs"]
mod round7_fixture;

use std::alloc::{GlobalAlloc, Layout, System};
use std::ffi::{CStr, CString, c_void};
use std::hint::black_box;
use std::path::Path;
use std::sync::atomic::{AtomicBool, AtomicU64, Ordering};
use std::time::{Duration, Instant};

use omicsweave_vcf::{
    BcfEncoder, BcfScratch, GetError, Record, VariantReader, bcf_header_bytes, read_vcf_header,
};

struct Counter;
static COUNT: AtomicBool = AtomicBool::new(false);
static ALLOCS: AtomicU64 = AtomicU64::new(0);
static BYTES: AtomicU64 = AtomicU64::new(0);
fn count(n: usize) {
    if COUNT.load(Ordering::Relaxed) {
        ALLOCS.fetch_add(1, Ordering::Relaxed);
        BYTES.fetch_add(n as u64, Ordering::Relaxed);
    }
}
// SAFETY: delegation preserves every System pointer/layout contract. Counters
// contain no pointers and do not read or write the allocated storage.
unsafe impl GlobalAlloc for Counter {
    unsafe fn alloc(&self, l: Layout) -> *mut u8 {
        count(l.size());
        // SAFETY: layout is supplied under GlobalAlloc's caller contract.
        unsafe { System.alloc(l) }
    }
    unsafe fn alloc_zeroed(&self, l: Layout) -> *mut u8 {
        count(l.size());
        // SAFETY: layout is supplied under GlobalAlloc's caller contract.
        unsafe { System.alloc_zeroed(l) }
    }
    unsafe fn dealloc(&self, p: *mut u8, l: Layout) {
        // SAFETY: pointer and matching layout belong to System.
        unsafe { System.dealloc(p, l) }
    }
    unsafe fn realloc(&self, p: *mut u8, l: Layout, n: usize) -> *mut u8 {
        count(n);
        // SAFETY: caller supplies a live System allocation and valid size.
        unsafe { System.realloc(p, l, n) }
    }
}
#[global_allocator]
static ALLOCATOR: Counter = Counter;

#[derive(Default)]
struct Hash(u64);
impl Hash {
    fn word(&mut self, n: u64) {
        self.0 = self.0.wrapping_mul(0x100000001b3).wrapping_add(n);
    }
    fn bytes(&mut self, b: &[u8]) {
        self.word(b.len() as u64);
        for &v in b {
            self.word(u64::from(v));
        }
    }
    fn values(&mut self, result: Result<usize, GetError>, values: impl Iterator<Item = u32>) {
        match result {
            Err(e) => self.word(e.code() as u64),
            Ok(n) => {
                self.word(n as u64);
                for v in values {
                    self.word(u64::from(v));
                }
            }
        }
    }
}

fn digest(bytes: &[u8], getters: bool) -> (usize, u64) {
    let mut reader = VariantReader::new(bytes).unwrap();
    let mut rec = Record::default();
    let mut hash = Hash::default();
    let mut count = 0;
    let (mut ints, mut floats, mut strings) = (Vec::new(), Vec::new(), Vec::new());
    while reader.read(&mut rec).unwrap() {
        count += 1;
        hash.word(rec.rid as u64);
        hash.word(rec.pos as u64);
        hash.word(rec.rlen as u64);
        hash.bytes(rec.id());
        hash.word(u64::from(rec.qual_bits()));
        hash.word(rec.n_sample() as u64);
        hash.word(rec.n_alleles() as u64);
        for a in rec.alleles() {
            hash.bytes(a);
        }
        hash.word(rec.filters().len() as u64);
        for &v in rec.filters() {
            hash.word(v as u64);
        }
        hash.word(rec.info_fields().len() as u64);
        for f in rec.info_fields() {
            hash.word(u64::from(f.key));
            hash.word(u64::from(f.values.ty));
            hash.word(f.values.len() as u64);
            hash.bytes(f.values.bytes);
        }
        hash.word(rec.format_fields().len() as u64);
        for f in rec.format_fields() {
            hash.word(u64::from(f.key));
            hash.word(u64::from(f.values.ty));
            hash.word(f.n as u64);
            hash.bytes(f.values.bytes);
        }
        if getters {
            let h = reader.header();
            for k in h.keys.iter().flatten() {
                hash.values(
                    rec.info_flag(h, &k.name).map(usize::from),
                    std::iter::empty(),
                );
                let r = rec.info_int32(h, &k.name, &mut ints);
                hash.values(r, ints.iter().map(|&v| v as u32));
                let r = rec.info_float(h, &k.name, &mut floats);
                hash.values(r, floats.iter().map(|v| v.to_bits()));
                let r = rec.info_string(h, &k.name, &mut strings);
                hash.values(r, strings.iter().map(|&v| u32::from(v)));
                let r = rec.format_int32(h, &k.name, &mut ints);
                hash.values(r, ints.iter().map(|&v| v as u32));
                let r = rec.format_float(h, &k.name, &mut floats);
                hash.values(r, floats.iter().map(|v| v.to_bits()));
                let r = rec.format_char(h, &k.name, &mut strings);
                hash.values(r, strings.iter().map(|&v| u32::from(v)));
            }
        }
    }
    (count, hash.0)
}

fn oracle(path: &Path, getters: bool) -> (usize, u64) {
    let path = CString::new(path.to_str().unwrap()).unwrap();
    let mut hash = Hash::default();
    let mut count = 0;
    // SAFETY: htslib owns the live handles until destruction below. Every
    // returned pointer is consumed before the next read. Zero-length spans
    // are handled without constructing Rust slices from null pointers.
    unsafe {
        htslib_sys::hts_set_log_level(htslib_sys::HTS_LOG_OFF);
        let fp = htslib_sys::hts_open(path.as_ptr(), c"r".as_ptr());
        assert!(!fp.is_null());
        let hdr = htslib_sys::bcf_hdr_read(fp);
        let rec = htslib_sys::bcf_init();
        assert!(!hdr.is_null() && !rec.is_null());
        loop {
            let r = htslib_sys::bcf_read(fp, hdr, rec);
            if r == -1 {
                break;
            }
            assert_eq!(r, 0);
            assert_eq!(htslib_sys::bcf_unpack(rec, htslib_sys::BCF_UN_ALL), 0);
            count += 1;
            hash.word(0);
            hash.word(htslib_sys::ow_bcf_pos(rec) as u64);
            hash.word(htslib_sys::ow_bcf_rlen(rec) as u64);
            hash.bytes(CStr::from_ptr(htslib_sys::ow_bcf_id(rec)).to_bytes());
            hash.word(u64::from(htslib_sys::ow_bcf_qual(rec).to_bits()));
            hash.word(htslib_sys::ow_bcf_n_sample(rec) as u64);
            let na = htslib_sys::ow_bcf_n_allele(rec);
            hash.word(na as u64);
            for i in 0..na {
                hash.bytes(CStr::from_ptr(htslib_sys::ow_bcf_allele(rec, i)).to_bytes());
            }
            let nf = htslib_sys::ow_bcf_n_flt(rec) as usize;
            hash.word(nf as u64);
            if nf != 0 {
                for &v in std::slice::from_raw_parts(htslib_sys::ow_bcf_flt(rec), nf) {
                    hash.word(v as u64);
                }
            }
            for format in [false, true] {
                let n = if format {
                    htslib_sys::ow_bcf_n_fmt(rec)
                } else {
                    htslib_sys::ow_bcf_n_info(rec)
                };
                hash.word(n as u64);
                for i in 0..n {
                    let (mut key, mut ty, mut n, mut p, mut nb) = (0, 0, 0, std::ptr::null(), 0);
                    if format {
                        htslib_sys::ow_bcf_fmt(rec, i, &mut key, &mut ty, &mut n, &mut p, &mut nb);
                    } else {
                        htslib_sys::ow_bcf_info(rec, i, &mut key, &mut ty, &mut n, &mut p, &mut nb);
                    }
                    hash.word(key as u64);
                    hash.word(ty as u64);
                    hash.word(n as u64);
                    hash.bytes(if nb == 0 {
                        &[]
                    } else {
                        std::slice::from_raw_parts(p, nb as usize)
                    });
                }
            }
            if getters {
                for i in 0..htslib_sys::ow_bcf_hdr_nkeys(hdr) {
                    let tag = htslib_sys::ow_bcf_hdr_key(hdr, i);
                    if tag.is_null() {
                        continue;
                    }
                    for (format, ty) in [
                        (false, 0),
                        (false, 1),
                        (false, 2),
                        (false, 3),
                        (true, 1),
                        (true, 2),
                        (true, 3),
                    ] {
                        let (mut p, mut cap): (*mut c_void, i32) = (std::ptr::null_mut(), 0);
                        let n = if format {
                            htslib_sys::bcf_get_format_values(hdr, rec, tag, &mut p, &mut cap, ty)
                        } else {
                            htslib_sys::bcf_get_info_values(hdr, rec, tag, &mut p, &mut cap, ty)
                        };
                        hash.word(n as u64);
                        if n > 0 && ty != 0 {
                            if ty == 3 {
                                for &v in std::slice::from_raw_parts(p as *const u8, n as usize) {
                                    hash.word(u64::from(v));
                                }
                            } else {
                                for &v in std::slice::from_raw_parts(p as *const u32, n as usize) {
                                    hash.word(u64::from(v));
                                }
                            }
                        }
                        htslib_sys::free(p);
                    }
                }
            }
        }
        htslib_sys::bcf_destroy(rec);
        htslib_sys::bcf_hdr_destroy(hdr);
        htslib_sys::hts_close(fp);
    }
    (count, hash.0)
}

fn fixture(version: &str, samples: usize, keys: usize, duplicate: bool, exponent: bool) -> String {
    let mut s = format!(
        "##fileformat=VCFv{version}\n##contig=<ID=chr1>\n##INFO=<ID=AF,Number=2,Type=Float>\n##INFO=<ID=TXT,Number=1,Type=String>\n##FORMAT=<ID=GT,Number=1,Type=String>\n##FORMAT=<ID=XF,Number=2,Type=Float>\n##FORMAT=<ID=XS,Number=1,Type=String>\n"
    );
    for i in 3..keys {
        s.push_str(&format!("##FORMAT=<ID=F{i},Number=1,Type=Integer>\n"));
    }
    s.push_str("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT");
    for i in 0..samples {
        s.push_str(&format!("\tS{i}"));
    }
    s.push('\n');
    for i in 0..256 {
        let q = if exponent { "1.25e2" } else { "50.25" };
        let af = if exponent { "1e-3,2.5e1" } else { "0.125,0.5" };
        s.push_str(&format!(
            "chr1\t{}\tid{i}\tA\tC\t{q}\tPASS\tAF={af};TXT=abc\tGT:XF:XS",
            i + 1
        ));
        for k in 3..keys {
            s.push_str(&format!(
                ":F{}",
                if duplicate { 3 + (k - 3) % 5 } else { k }
            ));
        }
        for sample in 0..samples {
            let gt = if sample % 2 == 0 { "0|1" } else { "0/0" };
            let xf = if exponent {
                "0.25e-1,2.5e+1"
            } else {
                "0.25,0.5"
            };
            s.push_str(&format!("\t{gt}:{xf}:sample"));
            for k in 3..keys {
                s.push_str(&format!(":{}", k % 127));
            }
        }
        s.push('\n');
    }
    s
}

fn measure(mut f: impl FnMut()) -> f64 {
    f();
    let mut groups = [0.0f64; 5];
    for sample in &mut groups {
        let start = Instant::now();
        let mut n = 0;
        while start.elapsed() < Duration::from_millis(100) {
            f();
            n += 1;
        }
        *sample = start.elapsed().as_secs_f64() * 1e9 / n as f64;
    }
    groups.sort_by(f64::total_cmp);
    groups[2]
}

fn main() {
    if std::env::var_os("VCF_ROUND7_BENCH").is_some() {
        round7();
        return;
    }
    let selected = std::env::args().skip(1).find(|a| !a.starts_with('-'));
    let dir = tempfile::tempdir().unwrap();
    for (name, version, samples, keys, dup, exp, binary) in [
        ("ordinary", "4.3", 4, 8, false, false, false),
        ("float-exponent", "4.4", 32, 8, false, true, false),
        ("many-samples", "4.4", 256, 8, false, false, false),
        ("keys128", "4.4", 4, 128, false, false, false),
        ("duplicates128", "4.4", 4, 128, true, false, false),
        ("keys255", "4.4", 4, 255, false, false, false),
        ("bcf43", "4.3", 32, 8, false, true, true),
        ("bcf44", "4.4", 32, 8, false, true, true),
    ] {
        if selected.as_deref().is_some_and(|s| s != name) {
            continue;
        }
        let text = fixture(version, samples, keys, dup, exp);
        let vcf = dir.path().join(format!("{name}.vcf"));
        std::fs::write(&vcf, &text).unwrap();
        let (path, bytes) = if binary {
            let path = dir.path().join(format!("{name}.bcf"));
            common::htslib_to_bcf(&vcf, &path);
            let bytes = std::fs::read(&path).unwrap();
            (path, bytes)
        } else {
            (vcf, text.into_bytes())
        };
        let raw = oracle(&path, false);
        let full = oracle(&path, true);
        assert_eq!(digest(&bytes, false), raw, "{name}: raw typed oracle");
        assert_eq!(digest(&bytes, true), full, "{name}: getters oracle");
        ALLOCS.store(0, Ordering::Relaxed);
        BYTES.store(0, Ordering::Relaxed);
        COUNT.store(true, Ordering::Relaxed);
        assert_eq!(digest(&bytes, false), raw);
        COUNT.store(false, Ordering::Relaxed);
        let allocs = ALLOCS.load(Ordering::Relaxed);
        let allocated = BYTES.load(Ordering::Relaxed);
        let ns = measure(|| {
            assert_eq!(black_box(digest(black_box(&bytes), false)), raw);
        });
        let getters_ns = measure(|| {
            assert_eq!(black_box(digest(black_box(&bytes), true)), full);
        });
        println!(
            "case={name} records={} bytes={} raw={:016x} getters={:016x} scan_ns={ns:.0} with_getters_ns={getters_ns:.0} allocation_requests={allocs} requested_bytes={allocated}",
            raw.0,
            bytes.len(),
            raw.1,
            full.1
        );
    }
}

fn round7() {
    use std::io::Read;

    let dir = tempfile::tempdir().unwrap();
    for (case, samples, width, keys, values) in [
        ("info-dictionary", 0, 1, 32, &["1", "127"][..]),
        ("scalar", 4, 1, 8, &["1", "127"][..]),
        ("samples256-i8", 256, 1, 8, &[".", "-120", "127"][..]),
        ("vectors-i8", 64, 32, 4, &[".", "-120", "127"][..]),
        ("vectors-i16", 64, 32, 4, &[".", "-32760", "32767"][..]),
        ("vectors-i32", 64, 32, 4, &[".", "-32761", "2147483647"][..]),
    ] {
        let text = round7_fixture::fixture(samples, width, keys, 64, values, "4.4");
        let vcf = dir.path().join("round7.vcf");
        std::fs::write(&vcf, &text).unwrap();
        let want = oracle(&vcf, false);
        assert_eq!(digest(text.as_bytes(), false), want);
        let parse_ns = measure(|| {
            assert_eq!(black_box(digest(black_box(text.as_bytes()), false)), want);
        });
        println!(
            "vcf.r7.{case}.read\t{parse_ns:.2}\tchecksum={:016x}",
            want.1
        );

        let mut input = text.as_bytes();
        let (header, consumed) = read_vcf_header(&mut input).unwrap();
        let encoder = BcfEncoder::new(&header);
        let mut scratch = BcfScratch::default();
        let lines: Vec<_> = input
            .split(|&b| b == b'\n')
            .filter(|l| !l.is_empty())
            .collect();
        let mut out = Vec::new();
        let mut encode = |out: &mut Vec<u8>| {
            out.clear();
            let mut checksum = 0u64;
            for (i, line) in lines.iter().enumerate() {
                let span = encoder
                    .encode_line(&mut scratch, line, consumed + 1 + i as u64, out)
                    .unwrap();
                checksum = checksum.wrapping_add(span.end as u64);
            }
            (out.len(), checksum)
        };
        let encoded = encode(&mut out);
        let bcf = dir.path().join("round7.bcf");
        common::htslib_to_bcf(&vcf, &bcf);
        let mut raw = Vec::new();
        flate2::read::MultiGzDecoder::new(std::fs::File::open(bcf).unwrap())
            .read_to_end(&mut raw)
            .unwrap();
        let header_bytes = bcf_header_bytes(&header).unwrap();
        assert_eq!(&raw[..header_bytes.len()], &header_bytes);
        assert_eq!(&raw[header_bytes.len()..], &out);
        let mut hash = Hash::default();
        hash.bytes(&out);
        let ns = measure(|| {
            assert_eq!(black_box(encode(black_box(&mut out))), encoded);
        });
        println!("vcf.r7.{case}.encode\t{ns:.2}\tchecksum={:016x}", hash.0);
    }
}
