//! ID, QUAL, FILTER and every INFO and FORMAT value of random VCFs, and of
//! the BCF htslib writes from them, as htslib and this crate decode them:
//! the stored typed encoding of each field and the results of htslib's
//! typed getters for every dictionary key must agree.

use crate::support::integer_fields;
use crate::support::rich_variants;
use crate::support::variants;

use std::collections::BTreeMap;
use std::ffi::{CStr, CString, c_int, c_void};
use std::path::Path;

use omicsweave_vcf::{GetError, Header, Record, VariantReader};
use proptest::prelude::*;
use variants::Rng;

/// A typed getter's outcome: values as 32-bit patterns (or bytes), or
/// htslib's negative return code.
type Got = Result<Vec<u32>, i32>;

#[derive(Debug, Clone, PartialEq, Eq)]
struct Snap {
    id: Vec<u8>,
    qual: u32,
    filters: Vec<i32>,
    n_sample: usize,
    /// (key, type, count, bytes)
    info: Vec<(i32, u8, usize, Vec<u8>)>,
    fmt: Vec<(i32, u8, usize, Vec<u8>)>,
    /// Per dictionary name: INFO flag, int, float, string; FORMAT int,
    /// float, string.
    got: BTreeMap<String, [Got; 7]>,
}

fn code(e: GetError) -> i32 {
    e.code()
}

fn ours_got(h: &Header, rec: &Record, name: &str) -> [Got; 7] {
    let mut ints = Vec::new();
    let mut floats = Vec::new();
    let mut bytes = Vec::new();
    let as_bytes = |b: &[u8]| b.iter().map(|&x| u32::from(x)).collect::<Vec<_>>();
    [
        rec.info_flag(h, name)
            .map(|f| vec![u32::from(f)])
            .map_err(code),
        rec.info_int32(h, name, &mut ints)
            .map(|_| ints.iter().map(|&v| v as u32).collect())
            .map_err(code),
        rec.info_float(h, name, &mut floats)
            .map(|_| floats.iter().map(|v| v.to_bits()).collect())
            .map_err(code),
        rec.info_string(h, name, &mut bytes)
            .map(|_| as_bytes(&bytes))
            .map_err(code),
        rec.format_int32(h, name, &mut ints)
            .map(|_| ints.iter().map(|&v| v as u32).collect())
            .map_err(code),
        rec.format_float(h, name, &mut floats)
            .map(|_| floats.iter().map(|v| v.to_bits()).collect())
            .map_err(code),
        rec.format_char(h, name, &mut bytes)
            .map(|_| as_bytes(&bytes))
            .map_err(code),
    ]
}

fn ours(bytes: &[u8]) -> Vec<Option<Snap>> {
    let mut reader = VariantReader::new(bytes).unwrap();
    let mut rec = Record::default();
    let mut out = Vec::new();
    loop {
        match reader.read(&mut rec) {
            Ok(false) => break,
            Err(_) => out.push(None),
            Ok(true) => {
                let h = reader.header();
                let got = h
                    .keys
                    .iter()
                    .flatten()
                    .map(|k| (k.name.clone(), ours_got(h, &rec, &k.name)))
                    .collect();
                out.push(Some(Snap {
                    id: rec.id().to_vec(),
                    qual: rec.qual_bits(),
                    filters: rec.filters().to_vec(),
                    n_sample: rec.n_sample(),
                    info: rec
                        .info_fields()
                        .map(|f| {
                            (
                                f.key as i32,
                                f.values.ty,
                                f.values.len(),
                                f.values.bytes.to_vec(),
                            )
                        })
                        .collect(),
                    fmt: rec
                        .format_fields()
                        .map(|f| (f.key as i32, f.values.ty, f.n, f.values.bytes.to_vec()))
                        .collect(),
                    got,
                }));
            }
        }
    }
    out
}

/// # Safety
/// `hdr` and `rec` are live htslib handles; `rec` is unpacked.
unsafe fn htslib_got(
    hdr: *mut htslib_sys::bcf_hdr_t,
    rec: *mut htslib_sys::bcf1_t,
    name: &str,
) -> [Got; 7] {
    let tag = CString::new(name).unwrap();
    let call = |format: bool, ty: c_int| -> Got {
        let (mut buf, mut cap): (*mut c_void, c_int) = (std::ptr::null_mut(), 0);
        // SAFETY: as the caller guarantees; `buf` is htslib-allocated.
        unsafe {
            let n = if format {
                htslib_sys::bcf_get_format_values(hdr, rec, tag.as_ptr(), &mut buf, &mut cap, ty)
            } else {
                htslib_sys::bcf_get_info_values(hdr, rec, tag.as_ptr(), &mut buf, &mut cap, ty)
            };
            let out = if n < 0 {
                Err(n)
            } else if ty == htslib_sys::BCF_HT_FLAG {
                Ok(vec![n as u32])
            } else if n == 0 {
                Ok(Vec::new())
            } else if ty == htslib_sys::BCF_HT_STR {
                Ok(std::slice::from_raw_parts(buf as *const u8, n as usize)
                    .iter()
                    .map(|&b| u32::from(b))
                    .collect())
            } else {
                Ok(std::slice::from_raw_parts(buf as *const u32, n as usize).to_vec())
            };
            htslib_sys::free(buf);
            out
        }
    };
    use htslib_sys::{BCF_HT_FLAG, BCF_HT_INT, BCF_HT_REAL, BCF_HT_STR};
    [
        call(false, BCF_HT_FLAG),
        call(false, BCF_HT_INT),
        call(false, BCF_HT_REAL),
        call(false, BCF_HT_STR),
        call(true, BCF_HT_INT),
        call(true, BCF_HT_REAL),
        call(true, BCF_HT_STR),
    ]
}

fn htslib(path: &Path) -> Vec<Option<Snap>> {
    let c = |s: &str| CString::new(s).unwrap();
    let mut out = Vec::new();
    // SAFETY: handles come from htslib and are released below; data htslib
    // returns is copied before the record is reused.
    unsafe {
        let fp = htslib_sys::hts_open(c(path.to_str().unwrap()).as_ptr(), c("r").as_ptr());
        let hdr = htslib_sys::bcf_hdr_read(fp);
        let rec = htslib_sys::bcf_init();
        loop {
            let r = htslib_sys::bcf_read(fp, hdr, rec);
            if r == -1 {
                break;
            }
            if r < -1 {
                out.push(None);
                continue;
            }
            htslib_sys::bcf_unpack(rec, htslib_sys::BCF_UN_ALL);
            let mut info = Vec::new();
            for i in 0..htslib_sys::ow_bcf_n_info(rec) {
                let (mut key, mut ty, mut len, mut p, mut nb) = (0, 0, 0, std::ptr::null(), 0);
                htslib_sys::ow_bcf_info(rec, i, &mut key, &mut ty, &mut len, &mut p, &mut nb);
                let bytes = std::slice::from_raw_parts(p, nb as usize).to_vec();
                info.push((key, ty as u8, len as usize, bytes));
            }
            let mut fmt = Vec::new();
            if htslib_sys::ow_bcf_n_sample(rec) > 0 {
                for i in 0..htslib_sys::ow_bcf_n_fmt(rec) {
                    let (mut key, mut ty, mut n, mut p, mut nb) = (0, 0, 0, std::ptr::null(), 0);
                    htslib_sys::ow_bcf_fmt(rec, i, &mut key, &mut ty, &mut n, &mut p, &mut nb);
                    let bytes = std::slice::from_raw_parts(p, nb as usize).to_vec();
                    fmt.push((key, ty as u8, n as usize, bytes));
                }
            }
            let mut got = BTreeMap::new();
            for i in 0..htslib_sys::ow_bcf_hdr_nkeys(hdr) {
                let key = htslib_sys::ow_bcf_hdr_key(hdr, i);
                if key.is_null() {
                    continue;
                }
                let name = CStr::from_ptr(key).to_str().unwrap().to_owned();
                let g = htslib_got(hdr, rec, &name);
                got.insert(name, g);
            }
            let flt = htslib_sys::ow_bcf_flt(rec);
            let n_flt = htslib_sys::ow_bcf_n_flt(rec) as usize;
            out.push(Some(Snap {
                id: CStr::from_ptr(htslib_sys::ow_bcf_id(rec))
                    .to_bytes()
                    .to_vec(),
                qual: htslib_sys::ow_bcf_qual(rec).to_bits(),
                filters: if n_flt == 0 {
                    Vec::new()
                } else {
                    std::slice::from_raw_parts(flt, n_flt).to_vec()
                },
                n_sample: htslib_sys::ow_bcf_n_sample(rec) as usize,
                info,
                fmt,
                got,
            }));
        }
        htslib_sys::bcf_destroy(rec);
        htslib_sys::bcf_hdr_destroy(hdr);
        htslib_sys::hts_close(fp);
    }
    out
}

fn same(ours: &[Option<Snap>], htslib: &[Option<Snap>], what: &str, text: &str) {
    let lines: Vec<&str> = text.lines().filter(|l| !l.starts_with('#')).collect();
    for (i, (a, b)) in ours.iter().zip(htslib).enumerate() {
        let line = lines.get(i).unwrap_or(&"?");
        match (a, b) {
            (Some(a), Some(b)) => {
                assert_eq!(
                    (&a.id, a.qual, &a.filters, a.n_sample),
                    (&b.id, b.qual, &b.filters, b.n_sample),
                    "{what} record {i} columns: {line}"
                );
                assert_eq!(a.info, b.info, "{what} record {i} INFO: {line}");
                assert_eq!(a.fmt, b.fmt, "{what} record {i} FORMAT: {line}");
                for (name, g) in &b.got {
                    assert_eq!(
                        a.got.get(name),
                        Some(g),
                        "{what} record {i} getters of {name}: {line}"
                    );
                }
                assert_eq!(
                    a.got.len(),
                    b.got.len(),
                    "{what} record {i} dictionary: {line}"
                );
            }
            (a, b) => assert_eq!(
                a.is_some(),
                b.is_some(),
                "{what} record {i} accepted: {line}"
            ),
        }
    }
    assert_eq!(ours.len(), htslib.len(), "{what}: record count");
}

fn run(seed: u64) {
    let mut rng = Rng(seed);
    let (text, _) = rich_variants::random_vcf(&mut rng);
    let dir = tempfile::tempdir().unwrap();
    let vcf = dir.path().join("t.vcf");
    std::fs::write(&vcf, &text).unwrap();
    let want = htslib(&vcf);
    same(&ours(text.as_bytes()), &want, "vcf", &text);

    let bcf = dir.path().join("t.bcf");
    let kept = rich_variants::htslib_to_bcf(&vcf, &bcf);
    let kept = rich_variants::kept_lines(&text, &kept);
    same(
        &ours(&std::fs::read(&bcf).unwrap()),
        &htslib(&bcf),
        "bcf",
        &kept,
    );
}

proptest! {
    #![proptest_config(ProptestConfig {
        cases: 400,
        failure_persistence: Some(Box::new(proptest::test_runner::FileFailurePersistence::Direct(
            concat!(env!("CARGO_MANIFEST_DIR"), "/tests/cases/fields/values.proptest-regressions")
        ))),
        ..ProptestConfig::default()
    })]

    #[test]
    fn fields_match_htslib(seed in any::<u64>()) {
        run(seed);
    }
}

#[test]
fn integer_workspace_widths_and_shared_dictionary_match_htslib() {
    let dir = tempfile::tempdir().unwrap();
    for version in ["4.3", "4.4"] {
        for samples in [1, 2, 17] {
            for values in [
                &[".", "-120", "127", "0"][..],
                &[".", "-121", "128", "-32760", "32767"][..],
                &[".", "-32761", "32768", "-2147483640", "2147483647"][..],
                &["."][..],
            ] {
                let text = integer_fields::fixture(samples, 17, 4, 20, values, version);
                let vcf = dir.path().join("integers.vcf");
                std::fs::write(&vcf, &text).unwrap();
                same(&ours(text.as_bytes()), &htslib(&vcf), "VCF integers", &text);
                let bcf = dir.path().join("integers.bcf");
                let kept = rich_variants::htslib_to_bcf(&vcf, &bcf);
                assert!(kept.iter().all(|&accepted| accepted));
                same(
                    &ours(&std::fs::read(&bcf).unwrap()),
                    &htslib(&bcf),
                    "BCF integers",
                    &text,
                );
            }
        }
    }
}

#[test]
fn dictionary_kind_completion_matches_htslib() {
    let text = "##fileformat=VCFv4.4\n##contig=<ID=chr1>\n\
        ##INFO=<ID=INFO_ONLY,Number=1,Type=Integer>\n\
        ##FORMAT=<ID=FORMAT_ONLY,Number=1,Type=Integer>\n\
        ##FILTER=<ID=FILTER_ONLY,Description=\"filter\">\n\
        #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS\n\
        chr1\t1\t.\tA\tC\t.\tFILTER_ONLY\tFORMAT_ONLY=word;FILTER_ONLY=abc;NEW=x\tINFO_ONLY:FILTER_ONLY:NEW\tfoo:bar:baz\n\
        chr1\t2\t.\tA\tC\t.\tPASS\tINFO_ONLY=127;FORMAT_ONLY=next\tFORMAT_ONLY:INFO_ONLY:NEW\t128:more:end\n";
    let dir = tempfile::tempdir().unwrap();
    let path = dir.path().join("kinds.vcf");
    std::fs::write(&path, text).unwrap();
    same(
        &ours(text.as_bytes()),
        &htslib(&path),
        "dictionary completion",
        text,
    );
}

#[test]
fn format_getters_reject_mismatched_header_sample_counts() {
    let text = b"##fileformat=VCFv4.4\n##contig=<ID=chr1>\n\
        ##FORMAT=<ID=GT,Number=1,Type=String>\n\
        ##FORMAT=<ID=DP,Number=1,Type=Integer>\n\
        ##FORMAT=<ID=S,Number=1,Type=String>\n\
        #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ts0\ts1\n\
        chr1\t1\t.\tA\tC\t.\t.\t.\tGT:DP:S\t0|1:5:a\t1:7:bc\n";
    let mut reader = VariantReader::new(&text[..]).unwrap();
    let mut rec = Record::default();
    assert!(reader.read(&mut rec).unwrap());
    let good = reader.header().clone();
    let mut ints = Vec::new();
    assert_eq!(rec.genotype_values(&good, &mut ints).unwrap(), ints.len());
    for samples in [
        vec!["s0".to_owned()],
        vec!["s0".to_owned(), "s1".to_owned(), "s2".to_owned()],
    ] {
        let mut mismatched = good.clone();
        mismatched.samples = samples;
        assert_eq!(
            rec.genotype_values(&mismatched, &mut ints),
            Err(GetError::Unsupported)
        );
        assert!(ints.is_empty());
        assert_eq!(
            rec.format_int32(&mismatched, "DP", &mut ints),
            Err(GetError::Unsupported)
        );
        assert_eq!(
            rec.format_char(&mismatched, "S", &mut Vec::new()),
            Err(GetError::Unsupported)
        );
    }
}
