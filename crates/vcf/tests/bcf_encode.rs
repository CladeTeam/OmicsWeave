//! VCF text encoded to BCF by this crate and by htslib (`vcf_hdr_read` and
//! `vcf_parse`, then `bcf_hdr_write` and `bcf_write` to `wbu`): the header
//! bytes, every record's bytes, which records are refused, and each
//! record's index span must agree.
//!
//! Records are compared against a header fixed before them. Random files
//! use each undeclared name once, so htslib's dummy declarations (added
//! while it refuses the record) never let it accept a later record this
//! crate refuses; the fixtures of the other suites, which reuse undeclared
//! names, are encoded against the header htslib completes in a first pass.

mod common;
mod rich;
mod round7_fixture;

use std::ffi::CString;
use std::io::Read;
use std::path::Path;

use common::Rng;
use omicsweave_vcf::{BcfEncoder, BcfScratch, BcfSpan, Header, bcf_header_bytes, read_vcf_header};
use proptest::prelude::*;

/// One side's output: the header bytes (`None` when the header is
/// refused) and, per data line, the record bytes and span (`None` when
/// the record is refused).
#[derive(Debug, PartialEq, Eq)]
struct Encoded {
    header: Option<Vec<u8>>,
    records: Vec<Option<(Vec<u8>, BcfSpan)>>,
}

fn data_lines(text: &[u8]) -> Vec<&[u8]> {
    let mut input = text;
    if read_vcf_header(&mut input).is_err() {
        return Vec::new();
    }
    let mut lines: Vec<&[u8]> = input.split(|&b| b == b'\n').collect();
    if lines.last() == Some(&&b""[..]) {
        lines.pop();
    }
    lines
}

fn rid(record: &[u8]) -> u32 {
    u32::from_le_bytes(record[8..12].try_into().unwrap())
}

/// htslib's encoding of the VCF at `path`. With `complete`, a first pass
/// registers undeclared names and both header and records use the header
/// it leaves.
fn htslib(path: &Path, complete: bool) -> Encoded {
    let c = |s: &str| CString::new(s).unwrap();
    let dir = tempfile::tempdir().unwrap();
    let out_path = dir.path().join("out.bcf");
    let mut accepted = Vec::new();
    // SAFETY: handles come from htslib, are checked before use and are
    // released before returning.
    unsafe {
        htslib_sys::hts_set_log_level(htslib_sys::HTS_LOG_OFF);
        let vcf = c(path.to_str().unwrap());
        let mut fp = htslib_sys::hts_open(vcf.as_ptr(), c("r").as_ptr());
        assert!(!fp.is_null());
        let hdr = htslib_sys::bcf_hdr_read(fp);
        if hdr.is_null() {
            htslib_sys::hts_close(fp);
            return Encoded {
                header: None,
                records: Vec::new(),
            };
        }
        let rec = htslib_sys::bcf_init();
        if complete {
            while htslib_sys::bcf_read(fp, hdr, rec) != -1 {}
            htslib_sys::hts_close(fp);
            fp = htslib_sys::hts_open(vcf.as_ptr(), c("r").as_ptr());
            htslib_sys::bcf_hdr_destroy(htslib_sys::bcf_hdr_read(fp));
        }
        let out = htslib_sys::hts_open(c(out_path.to_str().unwrap()).as_ptr(), c("wbu").as_ptr());
        assert!(!out.is_null());
        assert_eq!(htslib_sys::bcf_hdr_write(out, hdr), 0);
        loop {
            let r = htslib_sys::bcf_read(fp, hdr, rec);
            if r == -1 {
                break;
            }
            let span = (htslib_sys::ow_bcf_pos(rec), htslib_sys::ow_bcf_rlen(rec));
            accepted.push((r == 0 && htslib_sys::bcf_write(out, hdr, rec) == 0).then_some(span));
        }
        htslib_sys::bcf_destroy(rec);
        htslib_sys::hts_close(out);
        htslib_sys::bcf_hdr_destroy(hdr);
        htslib_sys::hts_close(fp);
    }
    // `wbu` is BGZF at level 0 for BCF.
    let mut bytes = Vec::new();
    flate2::read::MultiGzDecoder::new(std::fs::File::open(&out_path).unwrap())
        .read_to_end(&mut bytes)
        .unwrap();
    let l_text = u32::from_le_bytes(bytes[5..9].try_into().unwrap()) as usize;
    let header = bytes[..9 + l_text].to_vec();
    let mut at = 9 + l_text;
    let records = accepted
        .into_iter()
        .map(|span| {
            let (pos, rlen) = span?;
            let l_shared = u32::from_le_bytes(bytes[at..at + 4].try_into().unwrap()) as usize;
            let l_indiv = u32::from_le_bytes(bytes[at + 4..at + 8].try_into().unwrap()) as usize;
            let record = bytes[at..at + 8 + l_shared + l_indiv].to_vec();
            at += record.len();
            let span = BcfSpan {
                rid: rid(&record),
                beg: pos,
                end: pos + rlen,
            };
            Some((record, span))
        })
        .collect();
    assert_eq!(at, bytes.len(), "htslib wrote records it did not report");
    Encoded {
        header: Some(header),
        records,
    }
}

/// This crate's encoding of `text`, against `header_text` when given (as
/// VCF text) or else the header of `text` itself.
fn ours(text: &[u8], header_text: Option<&[u8]>) -> Encoded {
    let mut input = text;
    let (header, consumed) = match read_vcf_header(&mut input) {
        Ok(h) => h,
        Err(_) => {
            return Encoded {
                header: None,
                records: Vec::new(),
            };
        }
    };
    let header = match header_text {
        Some(t) => Header::from_vcf_text(t).unwrap(),
        None => header,
    };
    let encoder = BcfEncoder::new(&header);
    let mut scratch = BcfScratch::default();
    let mut out = vec![0xAA];
    let records = data_lines(text)
        .into_iter()
        .enumerate()
        .map(|(i, line)| {
            out.truncate(1);
            match encoder.encode_line(&mut scratch, line, consumed + 1 + i as u64, &mut out) {
                Ok(span) => Some((out[1..].to_vec(), span)),
                Err(e) => {
                    assert!(
                        matches!(e, omicsweave_core::Error::Format { line, .. }
                            if line == consumed + 1 + i as u64),
                        "{e:?}"
                    );
                    assert_eq!(out, [0xAA], "a refused record leaves the output as it was");
                    None
                }
            }
        })
        .collect();
    Encoded {
        header: Some(bcf_header_bytes(&header).unwrap()),
        records,
    }
}

fn hex(b: &[u8]) -> String {
    b.iter()
        .map(|x| format!("{x:02x}"))
        .collect::<Vec<_>>()
        .join(" ")
}

fn compare(ours: &Encoded, want: &Encoded, text: &str) {
    // Long fixtures are cut in messages; the line at fault is named.
    let text = text.get(..4000).unwrap_or(text);
    match (&ours.header, &want.header) {
        (Some(a), Some(b)) => assert!(
            a == b,
            "header bytes differ\nours:   {:?}\nhtslib: {:?}\ntext:\n{text}",
            String::from_utf8_lossy(a),
            String::from_utf8_lossy(b)
        ),
        (a, b) => assert_eq!(
            a.is_some(),
            b.is_some(),
            "header accepted by one side only\n{text}"
        ),
    }
    let lines = data_lines(text.as_bytes());
    for (i, (a, b)) in ours.records.iter().zip(&want.records).enumerate() {
        let line = String::from_utf8_lossy(lines.get(i).copied().unwrap_or(b"?"));
        match (a, b) {
            (Some((ra, sa)), Some((rb, sb))) => {
                assert!(
                    ra == rb,
                    "record {i} bytes differ: {line:?}\nours:   {}\nhtslib: {}\ntext:\n{text}",
                    hex(ra),
                    hex(rb)
                );
                assert_eq!(sa, sb, "record {i} span: {line:?}");
            }
            (a, b) => assert_eq!(
                a.is_some(),
                b.is_some(),
                "record {i} accepted by one side only: {line:?}\ntext:\n{text}"
            ),
        }
    }
    assert_eq!(
        ours.records.len(),
        want.records.len(),
        "record count\n{text}"
    );
}

/// htslib's BCF header bytes for header text given to `bcf_hdr_parse`;
/// `None` when it refuses the text.
fn htslib_header(text: &[u8]) -> Option<Vec<u8>> {
    let c = |s: &str| CString::new(s).unwrap();
    let dir = tempfile::tempdir().unwrap();
    let path = dir.path().join("h.bcf");
    let mut buf = text.to_vec();
    buf.push(0);
    // SAFETY: `buf` is NUL-terminated and outlives the call; the header
    // and file handles are released before returning.
    unsafe {
        htslib_sys::hts_set_log_level(htslib_sys::HTS_LOG_OFF);
        let hdr = htslib_sys::bcf_hdr_init(c("r").as_ptr());
        assert!(!hdr.is_null());
        if htslib_sys::bcf_hdr_parse(hdr, buf.as_mut_ptr().cast()) != 0 {
            htslib_sys::bcf_hdr_destroy(hdr);
            return None;
        }
        let out = htslib_sys::hts_open(c(path.to_str().unwrap()).as_ptr(), c("wbu").as_ptr());
        assert!(!out.is_null());
        assert_eq!(htslib_sys::bcf_hdr_write(out, hdr), 0);
        htslib_sys::hts_close(out);
        htslib_sys::bcf_hdr_destroy(hdr);
    }
    let mut bytes = Vec::new();
    flate2::read::MultiGzDecoder::new(std::fs::File::open(&path).unwrap())
        .read_to_end(&mut bytes)
        .unwrap();
    Some(bytes)
}

/// Header text parsed both ways must give the same BCF header, or be
/// refused by both.
fn check_header(text: &str) {
    let ours = Header::from_vcf_text(text.as_bytes())
        .ok()
        .map(|h| bcf_header_bytes(&h).unwrap());
    let want = htslib_header(text.as_bytes());
    match (&ours, &want) {
        (Some(a), Some(b)) => assert!(
            a == b,
            "header bytes differ\nours:   {:?}\nhtslib: {:?}\ntext:\n{text}",
            String::from_utf8_lossy(a),
            String::from_utf8_lossy(b)
        ),
        _ => assert_eq!(
            ours.is_some(),
            want.is_some(),
            "header accepted by one side only\n{text}"
        ),
    }
}

/// The header part of a VCF: up to the end of the first line that starts
/// with `#` but not `##`.
fn header_part(text: &str) -> &str {
    let mut end = 0;
    for line in text.split_inclusive('\n') {
        end += line.len();
        if line.starts_with('#') && !line.starts_with("##") {
            break;
        }
    }
    &text[..end]
}

/// Encodes `text` both ways against its own header.
fn check(text: &str) {
    let dir = tempfile::tempdir().unwrap();
    let path = dir.path().join("t.vcf");
    std::fs::write(&path, text).unwrap();
    compare(&ours(text.as_bytes(), None), &htslib(&path, false), text);
}

/// Encodes `text` both ways against the header htslib completes by
/// reading all records once.
fn check_completed(text: &str) {
    let dir = tempfile::tempdir().unwrap();
    let path = dir.path().join("t.vcf");
    std::fs::write(&path, text).unwrap();
    let want = htslib(&path, true);
    let Some(header) = &want.header else {
        panic!("htslib refused the header:\n{text}");
    };
    // The BCF header text, read back as VCF text (its IDX are kept).
    let completed = &header[9..header.len() - 1];
    compare(&ours(text.as_bytes(), Some(completed)), &want, text);
}

// Random files.

const INTS: &[&str] = &[
    ".",
    "0",
    "1",
    "-1",
    "-120",
    "-121",
    "127",
    "128",
    "-32760",
    "-32761",
    "32767",
    "32768",
    "2147483647",
    "-2147483640",
    "-2147483641",
    "3000000000",
    "+5",
    "12abc",
    "",
    "-",
    "007",
];
const FLOATS: &[&str] = &[
    ".", "0", "1.5", "-0", "nan", "inf", "-inf", "1e-50", "3.4e39", ".5", "abc", "1.5x", "0x1p3",
    "", "2.5E2", "-3.25",
];
const STRINGS: &[&str] = &[
    "abc",
    "a b",
    "x,y",
    ".",
    "het",
    "1|0",
    "Z",
    "longer_string_value",
];

#[derive(Clone)]
struct Decl {
    name: &'static str,
    ty: &'static str,
}

const INFO: &[(&str, &str, &str)] = &[
    ("I1", "Integer", "1"),
    ("I2", "Integer", "2"),
    ("IA", "Integer", "A"),
    ("IR", "Integer", "R"),
    ("IG", "Integer", "G"),
    ("IV", "Integer", "."),
    ("F1", "Float", "1"),
    ("FA", "Float", "A"),
    ("FV", "Float", "."),
    ("FL", "Flag", "0"),
    ("FL1", "Flag", "1"),
    ("S1", "String", "1"),
    ("SV", "String", "."),
    ("C1", "Character", "1"),
    ("END", "Integer", "1"),
    ("SVLEN", "Integer", "A"),
    ("TU", "Unknown", "1"),
    ("TQ", "\"Integer\"", "1"),
];
const FORMAT: &[(&str, &str, &str)] = &[
    ("DP", "Integer", "1"),
    ("AD", "Integer", "R"),
    ("PL", "Integer", "G"),
    ("HQ", "Integer", "2"),
    ("GL", "Float", "G"),
    ("GQ", "Float", "1"),
    ("FT", "String", "1"),
    ("XC", "Character", "1"),
    ("LEN", "Integer", "1"),
    ("PS", "Integer", "1"),
];

fn bases(rng: &mut Rng, n: u64) -> String {
    (0..n).map(|_| *rng.pick(b"ACGTN") as char).collect()
}

fn values(rng: &mut Rng, ty: &str) -> String {
    let n = *rng.pick(&[1, 1, 1, 2, 3, 6, 17]);
    (0..n)
        .map(|_| match ty {
            "Integer" => rng.pick(INTS).to_string(),
            "Float" => rng.pick(FLOATS).to_string(),
            _ => rng.pick(STRINGS).to_string(),
        })
        .collect::<Vec<_>>()
        .join(",")
}

fn genotype(rng: &mut Rng, n_alleles: u64, v44: bool) -> String {
    if rng.chance(4) {
        return rng.pick(&[".", "./.", ".|.", "10/2", "0/x"]).to_string();
    }
    let ploidy = *rng.pick(&[1, 1, 2, 2, 2, 2, 3, 4]);
    let mut s = String::new();
    if v44 && rng.chance(15) {
        s.push(*rng.pick(&['|', '/']));
    }
    for i in 0..ploidy {
        if i > 0 {
            s.push(*rng.pick(&['|', '/']));
        }
        if rng.chance(12) {
            s.push('.');
        } else {
            s += &rng.range(0, n_alleles - 1).to_string();
        }
    }
    s
}

fn format_value(rng: &mut Rng, d: &Decl) -> String {
    match d.ty {
        "Integer" => {
            if rng.chance(30) {
                rng.pick(INTS).to_string()
            } else {
                values(rng, "Integer")
            }
        }
        "Float" => values(rng, "Float"),
        // Never empty: htslib fails a first record whose FORMAT keys all
        // have zero width, depending on its history.
        _ => rng.pick(STRINGS).to_string(),
    }
}

/// A random VCF for BCF encoding: header lines of every kind htslib
/// normalises or drops, and records over all value types, sentinels and
/// widths, with refused records mixed in.
fn random_vcf(rng: &mut Rng) -> String {
    let version = *rng.pick(&["4.1", "4.2", "4.3", "4.4", "4.5"]);
    let v44 = version >= "4.4";
    let n_samples = *rng.pick(&[0, 1, 1, 2, 3, 5]);
    let mut meta: Vec<String> = Vec::new();
    for c in 1..=3 {
        meta.push(match rng.range(0, 5) {
            0 => format!("##contig=<ID=chr{c}>"),
            1 => format!("##contig=< ID = chr{c} , length = 1000000 >"),
            2 => format!("##contig=<ID=chr{c},length=\"1000000\",assembly=x>"),
            3 => format!("##contig=<ID=chr{c},IDX={}>", 20 + c),
            _ => format!("##contig=<ID=chr{c},length=1000000>"),
        });
    }
    let mut filters = vec!["PASS"];
    for f in ["q10", "s50"] {
        if rng.chance(80) {
            meta.push(format!(
                "##FILTER=<ID={f},Description=\"below \\\"{f}\\\"\">"
            ));
            filters.push(f);
        }
    }
    if rng.chance(20) {
        meta.push("##FILTER=<ID=PASS,Description=\"Mine\">".into());
    }
    let mut info = Vec::new();
    for &(name, ty, number) in INFO {
        if rng.chance(80) {
            meta.push(match rng.range(0, 6) {
                0 => format!("##INFO=<ID={name},Type={ty},Number={number}>"),
                1 => format!("##INFO=<ID={name},Number={number},Type={ty},Description=\"x, y\">  "),
                _ => format!(
                    "##INFO=<ID={name},Number={number},Type={ty},Description=\"{name} values\">"
                ),
            });
            let ty = match ty {
                "Integer" | "Float" | "Flag" => ty,
                _ => "String",
            };
            info.push(Decl { name, ty });
        }
    }
    let mut format = Vec::new();
    if n_samples > 0 || rng.chance(30) {
        if rng.chance(90) {
            meta.push("##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">".into());
            format.push(Decl {
                name: "GT",
                ty: "GT",
            });
        }
        for &(name, ty, number) in FORMAT {
            if rng.chance(70) {
                meta.push(format!(
                    "##FORMAT=<ID={name},Number={number},Type={ty},Description=\"d\">"
                ));
                let ty = if ty == "Character" { "String" } else { ty };
                format.push(Decl { name, ty });
            }
        }
    }
    for extra in [
        "##source=generator",
        "##source=generator",
        "##reference=file:///ref.fa",
        "##ALT=<ID=DEL,Description=\"Deletion\">",
        "##ALT=<ID=DEL,Description=\"Again\">",
        "##PEDIGREE=<Name_0=G0,Name_1=G1>",
        "##META=<ID=Assay,Type=String,Number=.,Values=[WholeGenome, Exome]>",
        "##broken",
        "##INFO=<ID=BAD,Number=1,Type=Integer,Description=\"unterminated>",
        "##INFO=<Number=1,Type=Integer>",
        "##INFO=<ID=DUPI,Number=1,Type=Integer>",
        "##INFO=<ID=DUPI,Number=2,Type=Float>",
        "##FORMAT=<ID=I1,Number=1,Type=Float>",
        "##FILTER=<ID=I2,Description=\"shared name\">",
    ] {
        if rng.chance(25) {
            meta.push(extra.into());
        }
    }
    // Some disorder; IDs follow declaration order.
    for _ in 0..meta.len() / 3 {
        let i = rng.range(0, meta.len() as u64 - 1) as usize;
        let j = rng.range(0, meta.len() as u64 - 1) as usize;
        meta.swap(i, j);
    }
    let fileformat = format!("##fileformat=VCFv{version}");
    match rng.range(0, 10) {
        0 => {}
        1 => {
            let i = rng.range(0, meta.len() as u64) as usize;
            meta.insert(i, fileformat);
        }
        _ => meta.insert(0, fileformat),
    }
    let mut columns = "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO".to_string();
    if n_samples > 0 {
        columns += "\tFORMAT";
        for s in 0..n_samples {
            columns += &format!("\tS{s}");
        }
    }
    meta.push(columns);

    let crlf = rng.chance(15);
    let mut v = String::new();
    let push = |v: &mut String, line: &str, rng: &mut Rng| {
        v.push_str(line);
        v.push_str(if crlf && rng.chance(50) { "\r\n" } else { "\n" });
    };
    for line in &meta {
        push(&mut v, line, rng);
    }

    for i in 0..rng.range(1, 40) {
        let chrom = if rng.chance(3) {
            format!("chrU{i}")
        } else {
            format!("chr{}", rng.range(1, 3))
        };
        let pos = match rng.range(0, 40) {
            0 => "0".to_string(),
            1 => rng
                .pick(&["2147483647", "2147483648", "+5", "x"])
                .to_string(),
            _ => rng.range(1, 1_000_000).to_string(),
        };
        let id = rng.pick(&[".", ".", "", "rs1", "rs1;rs22"]).to_string();
        let ref_len = *rng.pick(&[1, 1, 1, 2, 3, 5]);
        let reference = bases(rng, ref_len);
        let n_alt = *rng.pick(&[0, 1, 1, 1, 2, 3]);
        let alts: Vec<String> = (0..n_alt)
            .map(|_| match rng.range(0, 12) {
                0..=4 => bases(rng, 1),
                5 => format!("{reference}{}", bases(rng, 2)),
                6 => reference[..1].to_string(),
                7 => rng
                    .pick(&["<DEL>", "<DUP:TANDEM>", "<INV>", "<CNV>", "<INS>"])
                    .to_string(),
                8 => rng.pick(&["<*>", "<NON_REF>"]).to_string(),
                9 => "*".into(),
                10 => format!("{}[chr1:{}[", &reference[..1], rng.range(1, 999)),
                _ => "".into(),
            })
            .collect();
        let alt = if alts.is_empty() {
            ".".to_string()
        } else {
            alts.join(",")
        };
        let qual = rng
            .pick(&[
                ".", ".", "50", "3.25", "", "1e40", "nan", "-0", "abc", "0.1",
            ])
            .to_string();
        let filter = match rng.range(0, 10) {
            0 | 1 => ".".to_string(),
            2 => format!("UF{i}"),
            3 => format!("{};", rng.pick(&filters)),
            4 => format!("{};{}", rng.pick(&filters), rng.pick(&filters)),
            _ => rng.pick(&filters).to_string(),
        };
        let mut items: Vec<String> = Vec::new();
        for d in &info {
            if !rng.chance(40) {
                continue;
            }
            items.push(match d.ty {
                "Flag" if rng.chance(85) => d.name.to_string(),
                _ if rng.chance(5) => d.name.to_string(),
                "Flag" => format!("{}=1", d.name),
                _ if d.name == "END" && rng.chance(70) => {
                    let p: u64 = pos.parse().unwrap_or(1);
                    format!("END={}", p + rng.range(0, 3000))
                }
                ty => format!("{}={}", d.name, values(rng, ty)),
            });
        }
        if rng.chance(3) {
            items.push(format!("UI{i}=x"));
        }
        if rng.chance(4) && !items.is_empty() {
            items.push(items[0].clone());
        }
        let mut info_col = if items.is_empty() {
            ".".to_string()
        } else {
            items.join(";")
        };
        if rng.chance(5) {
            info_col += ";";
        }
        if rng.chance(3) {
            info_col = info_col.replacen(';', ";;", 1);
        }
        let mut line =
            format!("{chrom}\t{pos}\t{id}\t{reference}\t{alt}\t{qual}\t{filter}\t{info_col}");
        if n_samples > 0 && !rng.chance(3) {
            if rng.chance(3) {
                line += "\t.";
                for _ in 0..n_samples {
                    line += "\t.";
                }
            } else {
                let mut keys: Vec<Decl> =
                    format.iter().filter(|_| rng.chance(50)).cloned().collect();
                if keys.is_empty() {
                    if format.is_empty() {
                        keys.push(Decl {
                            name: "UNDECL",
                            ty: "String",
                        });
                    } else {
                        keys.push(format[0].clone());
                    }
                }
                if keys.len() > 1 && rng.chance(20) {
                    keys.swap(0, 1);
                }
                if rng.chance(3) {
                    keys.push(keys[0].clone());
                }
                let mut names: Vec<String> = keys.iter().map(|d| d.name.to_string()).collect();
                if names == ["UNDECL"] || rng.chance(3) {
                    names.push(format!("UM{i}"));
                    keys.push(Decl {
                        name: "UM",
                        ty: "String",
                    });
                }
                line += &format!("\t{}", names.join(":"));
                let columns = if rng.chance(3) {
                    n_samples + 1
                } else {
                    n_samples
                };
                for _ in 0..columns {
                    let mut fields: Vec<String> = keys
                        .iter()
                        .map(|d| match d.ty {
                            "GT" => genotype(rng, n_alt as u64 + 1, v44),
                            _ => format_value(rng, d),
                        })
                        .collect();
                    if fields.len() > 1 && rng.chance(15) {
                        fields.truncate(rng.range(1, fields.len() as u64 - 1) as usize);
                    }
                    if rng.chance(3) {
                        fields = vec![".".into()];
                    }
                    line += &format!("\t{}", fields.join(":"));
                }
            }
        }
        push(&mut v, &line, rng);
    }
    v
}

proptest! {
    #![proptest_config(ProptestConfig { cases: 300, ..ProptestConfig::default() })]

    #[test]
    fn random_files_encode_as_htslib(seed in any::<u64>()) {
        let text = random_vcf(&mut Rng(seed));
        check_header(header_part(&text));
        check(&text);
    }

    #[test]
    fn suite_fixtures_encode_as_htslib(seed in any::<u64>()) {
        let mut rng = Rng(seed);
        check_completed(&common::random_vcf(&mut rng));
        check_completed(&rich::random_vcf(&mut rng).0);
    }
}

const COLUMNS: &str = "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO";

/// Header fixtures: each is parsed both ways on its own and, when it
/// starts as htslib needs a VCF file to, encoded with a record.
fn header_fixtures() -> Vec<String> {
    let samples = format!("{COLUMNS}\tFORMAT\tS1\n");
    vec![
        // No fileformat line: the first declaration takes id 0, PASS 1.
        format!(
            "##INFO=<ID=DP,Number=1,Type=Integer>\n##FILTER=<ID=q10>\n##INFO=<ID=AF,Number=A,Type=Float>\n##FORMAT=<ID=GT,Number=1,Type=String>\n##contig=<ID=chr1>\n{samples}"
        ),
        // A user PASS replaces htslib's, keeping its place.
        format!(
            "##fileformat=VCFv4.3\n##FILTER=<ID=q10,Description=\"q\">\n##FILTER=<ID=PASS,Description=\"Mine\">\n##INFO=<ID=DP,Number=1,Type=Integer>\n##INFO=<ID=AF,Number=A,Type=Float>\n##FORMAT=<ID=GT,Number=1,Type=String>\n##FORMAT=<ID=DP,Number=1,Type=Integer>\n##contig=<ID=chr1,length=100>\n{samples}"
        ),
        // Without a first line, PASS is first and a later PASS swaps
        // with the second record when that is a FILTER.
        format!(
            "##FILTER=<ID=q10>\n##FILTER=<ID=PASS,Description=\"Mine\">\n##contig=<ID=chr1>\n{samples}"
        ),
        format!("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"),
        // A first structured line without ID is added twice.
        format!(
            "##PEDIGREE=<Name_0=G0>\n##fileformat=VCFv4.2\n##contig=<ID=chr1>\n##FILTER=<ID=q10>\n##INFO=<ID=DP,Number=1,Type=Integer>\n##INFO=<ID=AF,Number=A,Type=Float>\n##FORMAT=<ID=GT,Number=1,Type=String>\n##FORMAT=<ID=DP,Number=1,Type=Integer>\n{samples}"
        ),
        // Spacing, quoting, escapes, brackets, trailing junk, malformed lines.
        format!(
            "##fileformat=VCFv4.4\n##contig=< id = chr1 , length = \"77\" >trailing\n##FILTER=<ID=q10,Description=\"a \\\"q\\\" , b  \" >\n##INFO=<ID=DP,Number= 1 ,Type=Integer,Description=\"\">\n##INFO=<ID=AF,Number=A,Type=Float,Values=[a, b]>\n##bad line\n##INFO=<ID=X,Number=1\n##=x\n##INFO=<ID=Y,=1>\n##X=<A=[1,2>\n##Y=<A=\"\\\\\",B=<c,d>>\n##FORMAT=<ID=GT,Number=1,Type=String>\n##FORMAT=<ID=DP,Number=1,Type=Integer>\n{samples}"
        ),
        // Shared names across kinds, duplicates, IDX, odd Number and Type.
        format!(
            "##fileformat=VCFv4.2\n##contig=<ID=chr1,IDX=5>\n##contig=<ID=chr1,IDX=6>\n##FILTER=<ID=DP>\n##INFO=<ID=DP,Number=1,Type=Integer,IDX=9>\n##INFO=<ID=DP,Number=2,Type=Float>\n##FORMAT=<ID=DP,Number=1,Type=Integer>\n##INFO=<ID=AF,Number=0x,Type=Flag>\n##FILTER=<ID=q10,IDX=3>\n##FORMAT=<ID=GT,Number=1,Type=String,IDX=7>\n##INFO=<ID=Z,Number=1,Type=Integer,IDX=x>\n##contig=<ID=chr9,IDX=4294967296>\n##source=a\n##source=a\n##source=b\n##fileformat=VCFv4.4\n{samples}"
        ),
        // Empty lines and CRLF in the header.
        format!(
            "##fileformat=VCFv4.2\r\n\r\n##contig=<ID=chr1>\r\n\n##FILTER=<ID=q10>\r\n##INFO=<ID=DP,Number=1,Type=Integer>\n##INFO=<ID=AF,Number=A,Type=Float>\r\n##FORMAT=<ID=GT,Number=1,Type=String>\n##FORMAT=<ID=DP,Number=1,Type=Integer>\n{COLUMNS}\tFORMAT\tS1\r\n"
        ),
        // Refused: conflicting IDX, no sample line, bad #CHROM lines.
        format!("##FILTER=<ID=A,IDX=0>\n{COLUMNS}\n"),
        format!(
            "##fileformat=VCFv4.2\n##INFO=<ID=A,Number=1,Type=Integer,IDX=4>\n##INFO=<ID=B,Number=1,Type=Integer,IDX=4>\n{COLUMNS}\n"
        ),
        format!("##fileformat=VCFv4.2\n#comment\n{COLUMNS}\n"),
        format!("##fileformat=VCFv4.2\n{COLUMNS}\tFORMAT\n"),
        format!("##fileformat=VCFv4.2\n{COLUMNS}\tFORMAT\tA\tA\n"),
        format!("##fileformat=VCFv4.2\n{COLUMNS}\tFORMAT\tA\t\n"),
        format!("##fileformat=VCFv4.2\n#CHROM POS ID REF ALT QUAL FILTER INFO\n"),
        "##fileformat=VCFv4.2\n".to_string(),
    ]
}

#[test]
fn header_fixtures_parse_and_encode_as_htslib() {
    let record = "chr1\t10\t.\tA\tC\t.\tq10;PASS\tDP=3;AF=0.5\tGT:DP\t0/1:7\n";
    for header in header_fixtures() {
        check_header(&header);
        check(&format!("{header}{record}"));
    }
}

#[test]
fn records_without_samples_and_with_crlf() {
    check(&format!(
        "##fileformat=VCFv4.2\r\n##contig=<ID=1>\r\n##INFO=<ID=DP,Number=1,Type=Integer>\r\n\
         {COLUMNS}\r\n1\t1\t.\tA\t.\t.\t.\t.\r\n1\t2\t.\tA\tC\t5\tPASS\tDP=300\r\n\
         1\t3\t.\tA\tC\t5\tPASS\tDP=1\tGT\t0/1\r\n\r\n1\t4\t.\tA\tC\t.\t.\t.\n"
    ));
}

#[test]
fn integer_workspace_records_encode_identically_to_htslib() {
    for version in ["4.3", "4.4"] {
        for samples in [1, 2, 17] {
            for values in [
                &[".", "-120", "127", "0"][..],
                &[".", "-121", "128", "-32760", "32767"][..],
                &[".", "-32761", "32768", "-2147483640", "2147483647"][..],
                &["."][..],
            ] {
                check(&round7_fixture::fixture(
                    samples, 17, 4, 20, values, version,
                ));
            }
        }
    }
}
