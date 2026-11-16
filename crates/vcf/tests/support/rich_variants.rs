//! Random VCFs that exercise every value kind: INFO and FORMAT keys of all
//! types and `Number` kinds, missing and malformed values, undeclared keys
//! and FILTERs, multiallelic and symbolic alleles, several samples with
//! varying ploidy. Shared by field parity and malformed filter state tests.

#![allow(dead_code)]

use super::variants::Rng;

/// A declared INFO or FORMAT key.
#[derive(Debug, Clone)]
pub struct Decl {
    pub name: String,
    /// `Integer`, `Float`, `String`, `Character` or `Flag`.
    pub ty: &'static str,
    /// `1`, `2`, `A`, `R`, `G`, `.`, ...
    pub number: &'static str,
}

/// The keys a generated file uses, for building expressions over them.
#[derive(Debug, Clone, Default)]
pub struct Spec {
    pub info: Vec<Decl>,
    pub format: Vec<Decl>,
    pub filters: Vec<String>,
    pub samples: Vec<String>,
    pub undeclared_info: Vec<String>,
    pub undeclared_format: Vec<String>,
}

const INFO_NAMES: &[&str] = &[
    "DP", "AF", "MQ", "SB", "DB", "H2", "TXT", "CNT", "FR", "AA", "XI", "XF", "XS", "DP4", "AN",
    "AC", "END", "SVLEN", "Q",
];
const FORMAT_NAMES: &[&str] = &[
    "DP", "GQ", "AD", "PL", "FT", "XF", "XS", "HQ", "GL", "AF", "LEN",
];
const TYPES: &[&str] = &[
    "Integer",
    "Integer",
    "Float",
    "Float",
    "String",
    "Character",
    "Flag",
];
const NUMBERS: &[&str] = &["1", "1", "2", "A", "R", "G", ".", "3", "0"];

fn decl(rng: &mut Rng, name: &str, format: bool) -> Decl {
    // Some names get their conventional declaration more often than not.
    let fixed = match (name, format) {
        ("DP" | "AN", _) => Some(("Integer", "1")),
        ("AC", false) => Some(("Integer", "A")),
        ("AF", false) => Some(("Float", "A")),
        ("AD", true) => Some(("Integer", "R")),
        ("PL", true) => Some(("Integer", "G")),
        ("GQ", true) => Some(("Integer", "1")),
        ("DP4", false) => Some(("Integer", "4")),
        ("DB" | "H2", false) => Some(("Flag", "0")),
        ("END", false) => Some(("Integer", "1")),
        ("SVLEN", false) => Some(("Integer", "A")),
        _ => None,
    };
    let (ty, number) = match fixed {
        Some(d) if rng.chance(75) => d,
        _ => {
            let mut ty = *rng.pick(TYPES);
            if format && ty == "Flag" {
                ty = "String";
            }
            (ty, *rng.pick(NUMBERS))
        }
    };
    Decl {
        name: name.to_string(),
        ty,
        number,
    }
}

const BASES: &[u8] = b"ACGTNacgt";

fn bases(rng: &mut Rng, n: u64) -> String {
    (0..n).map(|_| *rng.pick(BASES) as char).collect()
}

fn alt(rng: &mut Rng, reference: &str) -> String {
    match rng.range(0, 12) {
        0..=3 => bases(rng, 1),
        4 => {
            let n = rng.range(1, 4);
            format!("{}{}", reference, bases(rng, n))
        }
        5 => reference[..1].to_string(),
        6 => bases(rng, reference.len() as u64),
        7 => rng
            .pick(&["<DEL>", "<DUP:TANDEM>", "<INV>", "<CNV>", "<INS>"])
            .to_string(),
        8 => rng.pick(&["<*>", "<NON_REF>", "<X>"]).to_string(),
        9 => "*".into(),
        10 => format!("{}[chr1:{}[", &reference[..1], rng.range(1, 999)),
        _ => {
            let n = rng.range(0, 3);
            format!("{}{}", &reference[..1], bases(rng, n))
        }
    }
}

pub fn int_value(rng: &mut Rng) -> String {
    match rng.range(0, 30) {
        0 => ".".into(),
        1 => format!("-{}", rng.range(1, 50)),
        2 => rng.range(100, 40_000).to_string(),
        3 => rng
            .pick(&[
                "3000000000",
                "-2147483645",
                "2147483647",
                "-120",
                "-121",
                "127",
                "128",
            ])
            .to_string(),
        4 => rng.pick(&["-", "+5", "12abc", "", "0"]).to_string(),
        _ => rng.range(0, 60).to_string(),
    }
}

pub fn float_value(rng: &mut Rng) -> String {
    match rng.range(0, 30) {
        0 => ".".into(),
        1 => format!("-{}.{}", rng.range(0, 9), rng.range(0, 99)),
        2 => rng
            .pick(&[
                "1e-3", "2.5E2", "nan", "inf", "-inf", "0x1p3", "3.4e39", ".5", "1e-50", "0.1",
            ])
            .to_string(),
        3 => rng.pick(&["", "abc", "1.5x"]).to_string(),
        4 => rng.range(0, 99).to_string(),
        _ => format!("{}.{}", rng.range(0, 99), rng.range(0, 999)),
    }
}

pub fn string_value(rng: &mut Rng) -> String {
    match rng.range(0, 14) {
        0 => ".".into(),
        1 => "".into(),
        2 => rng.pick(&["abc", "ABC", "aBc", "a.b", "x y"]).to_string(),
        3 => "abc,.,de".into(),
        4 => rng.pick(&["het", "hom", "1|0", "PASS"]).to_string(),
        _ => {
            let n = rng.range(1, 5);
            (0..n)
                .map(|_| *rng.pick(&['a', 'b', 'c', 'X', 'Y', '1', '_', '-']))
                .collect()
        }
    }
}

/// How many values a `Number` asks for at a site.
fn count(rng: &mut Rng, number: &str, n_alt: u64, ploidy: u64) -> u64 {
    let wanted = match number {
        "A" => n_alt,
        "R" => n_alt + 1,
        "G" => {
            let n = n_alt + 1;
            if ploidy == 1 { n } else { n * (n + 1) / 2 }
        }
        "." => rng.range(1, 4),
        n => n.parse().unwrap_or(1),
    };
    if rng.chance(12) {
        rng.range(1, 5)
    } else {
        wanted.max(1)
    }
}

fn values(rng: &mut Rng, d: &Decl, n_alt: u64, ploidy: u64) -> String {
    let n = count(rng, d.number, n_alt, ploidy);
    let one = |rng: &mut Rng| match d.ty {
        "Integer" => int_value(rng),
        "Float" => float_value(rng),
        _ => string_value(rng),
    };
    if rng.chance(8) {
        return ".".into();
    }
    (0..n).map(|_| one(rng)).collect::<Vec<_>>().join(",")
}

fn genotype(rng: &mut Rng, n_alleles: u64, v44: bool) -> String {
    let ploidy = *rng.pick(&[1, 2, 2, 2, 3]);
    let mut s = String::new();
    if v44 && rng.chance(15) {
        s.push(*rng.pick(&['|', '/']));
    }
    for i in 0..ploidy {
        if i > 0 {
            s.push(*rng.pick(&['|', '/', '/']));
        }
        if rng.chance(12) {
            s.push('.');
        } else {
            s.push_str(&rng.range(0, n_alleles - 1).to_string());
        }
    }
    if rng.chance(3) {
        return ".".into();
    }
    s
}

/// A random VCF and the keys it uses.
pub fn random_vcf(rng: &mut Rng) -> (String, Spec) {
    let version = *rng.pick(&["4.2", "4.3", "4.4", "4.5"]);
    let v44 = version >= "4.4";
    let mut spec = Spec::default();
    let n_samples = rng.range(0, 4);
    let mut v = format!("##fileformat=VCFv{version}\n");
    for i in 0..rng.range(0, 3) {
        v += &format!("##contig=<ID=chr{},length=1000000>\n", i + 1);
    }
    for f in ["q10", "s50", "LowQual"] {
        if rng.chance(70) {
            v += &format!("##FILTER=<ID={f},Description=\"x\">\n");
            spec.filters.push(f.to_string());
        }
    }
    spec.filters.push("PASS".into());
    let mut info_names: Vec<&str> = INFO_NAMES.to_vec();
    for _ in 0..rng.range(3, 10) {
        let i = rng.range(0, info_names.len() as u64 - 1) as usize;
        let name = info_names.remove(i);
        let d = decl(rng, name, false);
        let number = if d.ty == "Flag" && rng.chance(80) {
            "0"
        } else {
            d.number
        };
        v += &format!(
            "##INFO=<ID={},Number={},Type={},Description=\"x\">\n",
            d.name, number, d.ty
        );
        spec.info.push(d);
    }
    if n_samples > 0 {
        if rng.chance(90) {
            v += "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n";
            spec.format.push(Decl {
                name: "GT".into(),
                ty: "String",
                number: "1",
            });
        }
        let mut names: Vec<&str> = FORMAT_NAMES.to_vec();
        for _ in 0..rng.range(1, 6) {
            let i = rng.range(0, names.len() as u64 - 1) as usize;
            let name = names.remove(i);
            let d = decl(rng, name, true);
            v += &format!(
                "##FORMAT=<ID={},Number={},Type={},Description=\"x\">\n",
                d.name, d.number, d.ty
            );
            spec.format.push(d);
        }
    }
    v += "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO";
    if n_samples > 0 {
        v += "\tFORMAT";
        for s in 0..n_samples {
            let name = format!("S{s}");
            v += &format!("\t{name}");
            spec.samples.push(name);
        }
    }
    v.push('\n');

    for _ in 0..rng.range(1, 40) {
        let chrom = format!("chr{}", rng.range(1, 4));
        let pos = rng.range(1, 900_000);
        let ref_len = *rng.pick(&[1, 1, 1, 2, 3, 5]);
        let reference = bases(rng, ref_len);
        let n_alt = rng.range(0, 3);
        let alts: Vec<String> = (0..n_alt).map(|_| alt(rng, &reference)).collect();
        let alt_field = if alts.is_empty() {
            ".".to_string()
        } else {
            alts.join(",")
        };
        let id = match rng.range(0, 6) {
            0 | 1 => ".".to_string(),
            2 => format!("rs{};rs{}", rng.range(1, 9), rng.range(1, 9)),
            _ => format!("rs{}", rng.range(1, 30)),
        };
        let qual = match rng.range(0, 10) {
            0 | 1 => ".".to_string(),
            2 => rng
                .pick(&["1e3", "abc", "", "59.2", "0x10", "nan", "-3"])
                .to_string(),
            3 => format!("{}.{}", rng.range(0, 99), rng.range(0, 9)),
            _ => rng.range(0, 200).to_string(),
        };
        let filter = match rng.range(0, 8) {
            0 | 1 => ".".to_string(),
            2 | 3 => "PASS".to_string(),
            4 => "XUNDEF".to_string(),
            _ => {
                let mut fs: Vec<String> = Vec::new();
                for _ in 0..rng.range(1, 3) {
                    let f = rng.pick(&spec.filters).clone();
                    if !fs.contains(&f) {
                        fs.push(f);
                    }
                }
                fs.join(";")
            }
        };
        let mut info = Vec::new();
        for d in &spec.info {
            if !rng.chance(55) {
                continue;
            }
            if d.ty == "Flag" {
                info.push(if rng.chance(10) {
                    format!("{}=1", d.name)
                } else {
                    d.name.clone()
                });
            } else if rng.chance(4) {
                info.push(d.name.clone());
            } else {
                info.push(format!("{}={}", d.name, values(rng, d, n_alt, 2)));
            }
        }
        if rng.chance(10) {
            let name = rng.pick(&["UNDEF", "ZZ"]).to_string();
            info.push(format!("{name}={}", string_value(rng)));
            if !spec.undeclared_info.contains(&name) {
                spec.undeclared_info.push(name);
            }
        }
        if rng.chance(5) && !info.is_empty() {
            let dup = info[0].clone();
            info.push(dup);
        }
        // Shuffle a little.
        if info.len() > 1 && rng.chance(50) {
            let i = rng.range(0, info.len() as u64 - 1) as usize;
            info.swap(0, i);
        }
        let info = if info.is_empty() {
            ".".to_string()
        } else {
            info.join(";")
        };
        let mut line =
            format!("{chrom}\t{pos}\t{id}\t{reference}\t{alt_field}\t{qual}\t{filter}\t{info}");
        if n_samples > 0 && rng.chance(92) {
            let mut keys: Vec<&Decl> = spec.format.iter().filter(|_| rng.chance(70)).collect();
            if keys.is_empty() {
                keys.push(&spec.format[0]);
            }
            if keys.len() > 1 && rng.chance(20) {
                keys.swap(0, 1);
            }
            let mut names: Vec<String> = keys.iter().map(|d| d.name.clone()).collect();
            if rng.chance(5) {
                names.push("UF".into());
                if !spec.undeclared_format.contains(&"UF".to_string()) {
                    spec.undeclared_format.push("UF".into());
                }
            }
            if rng.chance(2) {
                line += "\t.";
                for _ in 0..n_samples {
                    line += "\t.";
                }
            } else {
                line += &format!("\t{}", names.join(":"));
                for _ in 0..n_samples {
                    let ploidy = rng.range(1, 2);
                    let mut fields: Vec<String> = names
                        .iter()
                        .map(|name| match spec.format.iter().find(|d| &d.name == name) {
                            Some(d) if d.name == "GT" => genotype(rng, n_alt + 1, v44),
                            Some(d) => values(rng, d, n_alt, ploidy),
                            None => string_value(rng),
                        })
                        .collect();
                    // Keep a GT field: htslib exits on a GT key no sample gives.
                    if names[0] == "GT" && fields.len() > 1 && rng.chance(12) {
                        fields.truncate(rng.range(1, fields.len() as u64 - 1) as usize);
                    }
                    line += &format!("\t{}", fields.join(":"));
                }
            }
        }
        v += &line;
        v.push('\n');
    }
    (v, spec)
}

/// Rewrites a VCF as BCF with htslib, keeping the records it both reads
/// and writes (it refuses, for one, records without sample columns when
/// the header has samples). Returns, per input record, whether it was kept.
pub fn htslib_to_bcf(vcf: &std::path::Path, bcf: &std::path::Path) -> Vec<bool> {
    use std::ffi::CString;
    let c = |s: &str| CString::new(s).unwrap();
    let mut kept = Vec::new();
    // SAFETY: handles come from htslib and are released below.
    unsafe {
        let path = c(vcf.to_str().unwrap());
        // A first pass registers undeclared contigs and keys in the header;
        // htslib refuses to write records read before they were declared.
        let fp = htslib_sys::hts_open(path.as_ptr(), c("r").as_ptr());
        let hdr = htslib_sys::bcf_hdr_read(fp);
        let rec = htslib_sys::bcf_init();
        while htslib_sys::bcf_read(fp, hdr, rec) != -1 {}
        htslib_sys::hts_close(fp);

        let fp = htslib_sys::hts_open(path.as_ptr(), c("r").as_ptr());
        htslib_sys::bcf_hdr_destroy(htslib_sys::bcf_hdr_read(fp));
        let out = htslib_sys::hts_open(c(bcf.to_str().unwrap()).as_ptr(), c("wb").as_ptr());
        assert!(!fp.is_null() && !out.is_null());
        assert_eq!(htslib_sys::bcf_hdr_write(out, hdr), 0);
        loop {
            let r = htslib_sys::bcf_read(fp, hdr, rec);
            if r == -1 {
                break;
            }
            kept.push(r == 0 && htslib_sys::bcf_write(out, hdr, rec) == 0);
        }
        htslib_sys::bcf_destroy(rec);
        htslib_sys::bcf_hdr_destroy(hdr);
        htslib_sys::hts_close(out);
        htslib_sys::hts_close(fp);
    }
    kept
}

/// The header and the data lines `kept` marks, as text.
pub fn kept_lines(text: &str, kept: &[bool]) -> String {
    text.lines()
        .filter(|l| l.starts_with('#'))
        .chain(
            text.lines()
                .filter(|l| !l.starts_with('#'))
                .zip(kept)
                .filter(|(_, k)| **k)
                .map(|(l, _)| l),
        )
        .map(|l| format!("{l}\n"))
        .collect()
}
