//! Independent numeric and failure-state checks for the complete field reader.

use std::ffi::{CString, c_void};

use omicsweave_vcf::{FormatField, Record, Typed, VariantReader};

const HEADER: &str = "##fileformat=VCFv4.4\n##contig=<ID=chr1>\n\
    ##INFO=<ID=XI,Number=.,Type=Float,Description=\"values\">\n\
    ##FORMAT=<ID=XF,Number=.,Type=Float,Description=\"values\">\n\
    #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS\n";

#[derive(Debug, PartialEq, Eq)]
struct Numeric {
    qual: u32,
    info: Vec<u32>,
    format: Vec<u32>,
}

fn ours(text: &[u8]) -> Vec<Numeric> {
    let mut reader = VariantReader::new(text).unwrap();
    let mut record = Record::default();
    let mut out = Vec::new();
    let mut values = Vec::new();
    while reader.read(&mut record).unwrap() {
        record
            .info_float(reader.header(), "XI", &mut values)
            .unwrap();
        let info = values.iter().map(|v| v.to_bits()).collect();
        record
            .format_float(reader.header(), "XF", &mut values)
            .unwrap();
        out.push(Numeric {
            qual: record.qual_bits(),
            info,
            format: values.iter().map(|v| v.to_bits()).collect(),
        });
    }
    out
}

fn oracle(text: &[u8]) -> Vec<Numeric> {
    let directory = tempfile::tempdir().unwrap();
    let path = directory.path().join("numbers.vcf");
    std::fs::write(&path, text).unwrap();
    let path = CString::new(path.to_str().unwrap()).unwrap();
    let mut out = Vec::new();
    // SAFETY: paths and tags are terminated C strings and live through each
    // call. Handles are checked before use and destroyed after their final use.
    // Both float getters allocate an aligned buffer and report its initialized
    // element count. We copy only that count, before any reuse, and free the
    // allocation with htslib's matching free. No Rust reference aliases a write.
    unsafe {
        let input = htslib_sys::hts_open(path.as_ptr(), c"r".as_ptr());
        assert!(!input.is_null());
        let header = htslib_sys::bcf_hdr_read(input);
        assert!(!header.is_null());
        let record = htslib_sys::bcf_init();
        assert!(!record.is_null());
        let (mut buffer, mut capacity): (*mut c_void, i32) = (std::ptr::null_mut(), 0);
        loop {
            let status = htslib_sys::bcf_read(input, header, record);
            if status == -1 {
                break;
            }
            assert_eq!(status, 0);
            assert_eq!(htslib_sys::bcf_unpack(record, htslib_sys::BCF_UN_ALL), 0);
            let n = htslib_sys::bcf_get_info_values(
                header,
                record,
                c"XI".as_ptr(),
                &mut buffer,
                &mut capacity,
                htslib_sys::BCF_HT_REAL,
            );
            assert!(n > 0 && !buffer.is_null());
            let info = std::slice::from_raw_parts(buffer.cast::<u32>(), n as usize).to_vec();
            let n = htslib_sys::bcf_get_format_values(
                header,
                record,
                c"XF".as_ptr(),
                &mut buffer,
                &mut capacity,
                htslib_sys::BCF_HT_REAL,
            );
            assert!(n > 0 && !buffer.is_null());
            out.push(Numeric {
                qual: htslib_sys::ow_bcf_qual(record).to_bits(),
                info,
                format: std::slice::from_raw_parts(buffer.cast::<u32>(), n as usize).to_vec(),
            });
        }
        htslib_sys::free(buffer);
        htslib_sys::bcf_destroy(record);
        htslib_sys::bcf_hdr_destroy(header);
        assert_eq!(htslib_sys::hts_close(input), 0);
    }
    out
}

fn numeric_corpus(values: &[&str]) -> Vec<u8> {
    let mut text = HEADER.to_owned();
    for (i, value) in values.iter().enumerate() {
        text.push_str(&format!(
            "chr1\t{}\t.\tA\tC\t{value}\t.\tXI={value}\tXF\t{value}\n",
            i + 1,
        ));
    }
    text.into_bytes()
}

#[test]
fn decimal_exponents_hex_rounding_and_signed_zero_match_htslib_bits() {
    let values = [
        "1e000000000001",
        "1e-000000000001",
        "1e+0000000000000000000000002",
        "1e0000000000000000000000000",
        "0e99999999999999999999999",
        "-0",
        "-0.000e-99999999999999999999",
        "1e99999999999999999999999",
        "-1e-99999999999999999999999",
        "0x1p-149",
        "0x1p-150",
        "0x1.0000000000001p-150",
        "0x1.fffffep127",
        "0x1.ffffffp127",
        "0x1p-1074",
        "0x0.fffffffffffff8p-1022",
        "0x1.00000000000008p0",
        "0x1.0000000000000800000000000001p0",
        "1.1754943508222875e-38",
        "1.401298464324817e-45",
        "3.4028235677973366e38",
        ".5",
        "12.5",
    ];
    let mut text = numeric_corpus(&values);
    text.extend_from_slice(b"chr1\t100\t.\tA\tC\t12.5tail\t.\tXI=12.5\tXF\t12.5\n");
    assert_eq!(ours(&text), oracle(&text));
}

#[test]
fn nan_payloads_keep_the_bits_that_survive_float32_conversion() {
    let values = [
        "nan",
        "-nan",
        "nan(0x4000000000000)",
        "-nan(0x4000000000000)",
        "nan(0x2000000000000)",
        "nan(1125899906842624)",
        "nan(payload)",
        "nan()",
        "nan(08)",
        "nan(0x)",
        "nan(040000000000000000)",
        "nan(18446744073709551615)",
        "nan(18446744073709551616)",
        "nan(9999999999999999999999999999999999999999)",
    ];
    let text = numeric_corpus(&values);
    assert_eq!(ours(&text), oracle(&text));
}

#[test]
fn full_fields_are_empty_after_error_and_reusable_after_wide_record() {
    let mut text = HEADER.to_owned();
    let wide = (0..257)
        .map(|i| i.to_string())
        .collect::<Vec<_>>()
        .join(",");
    text.push_str(&format!(
        "chr1\t1\twide\tA\tC\t9\tPASS\tXI={wide}\tXF\t{wide}\n"
    ));
    text.push_str("chr1\t2\tbad\tA\tC\t8\tPASS\tXI=7\tXF\t1:q\n");
    text.push_str("chr1\t3\tshort\tA\tC\t7\t.\tXI=4\tXF\t5\n");
    let mut reader = VariantReader::new(text.as_bytes()).unwrap();
    let mut record = Record::default();
    let mut values = Vec::new();
    assert!(reader.read(&mut record).unwrap());
    assert_eq!(
        record
            .info_float(reader.header(), "XI", &mut values)
            .unwrap(),
        257
    );
    assert!(reader.read(&mut record).is_err());
    assert_eq!(record, Record::default());
    assert!(reader.read(&mut record).unwrap());
    assert_eq!(record.id(), b"short");
    assert_eq!(
        record
            .format_float(reader.header(), "XF", &mut values)
            .unwrap(),
        1
    );
    assert_eq!(values, [5.0]);
    assert_eq!(
        record
            .info_float(reader.header(), "XI", &mut values)
            .unwrap(),
        1
    );
    assert_eq!(values, [4.0]);
    assert!(!reader.read(&mut record).unwrap());
    assert_eq!(record, Record::default());
}

#[test]
fn public_typed_views_reject_overflowing_offsets_instead_of_aliasing_first_value() {
    for (ty, width) in [(2, 2), (3, 4), (4, 8), (5, 4)] {
        let bytes = [1; 8];
        let typed = Typed {
            ty,
            bytes: &bytes[..width],
        };
        let index = usize::MAX / width + 1;
        assert!(std::panic::catch_unwind(|| typed.raw(index)).is_err());
        let field = FormatField {
            key: 0,
            n: 1,
            values: typed,
        };
        assert!(std::panic::catch_unwind(|| field.sample(index)).is_err());
        let field = FormatField { n: index, ..field };
        assert!(std::panic::catch_unwind(|| field.size()).is_err());
    }
}
