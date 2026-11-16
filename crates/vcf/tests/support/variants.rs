//! htslib as the reference reader, and random VCF fixtures.

#![allow(dead_code)]

use std::ffi::{CStr, CString, c_void};
use std::path::Path;

use omicsweave_vcf::{Record, VariantReader};

/// splitmix64, for deterministic fixtures from a proptest-chosen seed.
pub struct Rng(pub u64);

impl Rng {
    pub fn next(&mut self) -> u64 {
        self.0 = self.0.wrapping_add(0x9E37_79B9_7F4A_7C15);
        let mut z = self.0;
        z = (z ^ (z >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
        z = (z ^ (z >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
        z ^ (z >> 31)
    }

    pub fn range(&mut self, lo: u64, hi: u64) -> u64 {
        lo + self.next() % (hi - lo + 1)
    }

    pub fn chance(&mut self, percent: u64) -> bool {
        self.next() % 100 < percent
    }

    pub fn pick<'a, T>(&mut self, items: &'a [T]) -> &'a T {
        &items[self.range(0, items.len() as u64 - 1) as usize]
    }
}

/// What both readers report of a record.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Seen {
    pub chrom: String,
    pub pos: i64,
    pub rlen: i64,
    pub alleles: Vec<Vec<u8>>,
    pub types: i32,
    /// Genotype values, all samples, when the record has them.
    pub gt: Option<Vec<i32>>,
}

/// Records as htslib reads them; `None` for a record it rejects.
pub fn htslib(path: &Path) -> Vec<Option<Seen>> {
    let c = |s: &str| CString::new(s).unwrap();
    let mut out = Vec::new();
    // SAFETY: handles come from htslib and are released below; strings
    // htslib returns are copied before the record is reused.
    unsafe {
        let fp = htslib_sys::hts_open(c(path.to_str().unwrap()).as_ptr(), c("r").as_ptr());
        assert!(!fp.is_null());
        let hdr = htslib_sys::bcf_hdr_read(fp);
        assert!(!hdr.is_null());
        let rec = htslib_sys::bcf_init();
        let gt = c("GT");
        let (mut buf, mut cap): (*mut c_void, i32) = (std::ptr::null_mut(), 0);
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
            let n = htslib_sys::ow_bcf_n_allele(rec);
            let alleles = (0..n)
                .map(|i| {
                    CStr::from_ptr(htslib_sys::ow_bcf_allele(rec, i))
                        .to_bytes()
                        .to_vec()
                })
                .collect();
            let ngt = htslib_sys::bcf_get_format_values(
                hdr,
                rec,
                gt.as_ptr(),
                &mut buf,
                &mut cap,
                htslib_sys::BCF_HT_INT,
            );
            let gt = (ngt >= 0)
                .then(|| std::slice::from_raw_parts(buf as *const i32, ngt as usize).to_vec());
            out.push(Some(Seen {
                chrom: CStr::from_ptr(htslib_sys::ow_bcf_seqname(hdr, rec))
                    .to_string_lossy()
                    .into_owned(),
                pos: htslib_sys::ow_bcf_pos(rec),
                rlen: htslib_sys::ow_bcf_rlen(rec),
                alleles,
                types: htslib_sys::bcf_get_variant_types(rec),
                gt,
            }));
        }
        htslib_sys::free(buf);
        htslib_sys::bcf_destroy(rec);
        htslib_sys::bcf_hdr_destroy(hdr);
        htslib_sys::hts_close(fp);
    }
    out
}

/// Rewrites a VCF as BCF with htslib (records it rejects are dropped).
pub fn htslib_to_bcf(vcf: &Path, bcf: &Path) {
    let c = |s: &str| CString::new(s).unwrap();
    // SAFETY: as in `htslib`.
    unsafe {
        let path = c(vcf.to_str().unwrap());
        // A first pass registers undeclared contigs and keys in the header;
        // htslib refuses to write records read before they were declared,
        // so the second pass reads them against the completed header.
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
            if r == 0 {
                assert_eq!(htslib_sys::bcf_write(out, hdr, rec), 0);
            }
        }
        htslib_sys::bcf_destroy(rec);
        htslib_sys::bcf_hdr_destroy(hdr);
        htslib_sys::hts_close(out);
        htslib_sys::hts_close(fp);
    }
}

/// Records as this crate reads them; `None` for a record it rejects.
pub fn ours(bytes: &[u8]) -> Vec<Option<Seen>> {
    let mut reader = VariantReader::new(bytes).unwrap();
    let mut rec = Record::default();
    let mut out = Vec::new();
    loop {
        match reader.read(&mut rec) {
            Ok(false) => break,
            Ok(true) => out.push(Some(Seen {
                chrom: reader.header().contigs[rec.rid as usize].name.clone(),
                pos: rec.pos,
                rlen: rec.rlen,
                alleles: rec.alleles().map(<[u8]>::to_vec).collect(),
                types: rec.variant_types(),
                gt: rec.genotypes().map(|g| g.values.to_vec()),
            })),
            Err(_) => out.push(None),
        }
    }
    out
}

const BASES: &[u8] = b"ACGTNacgt";

fn bases(rng: &mut Rng, n: u64) -> String {
    (0..n).map(|_| *rng.pick(BASES) as char).collect()
}

fn alt(rng: &mut Rng, reference: &str) -> String {
    match rng.range(0, 14) {
        0..=3 => bases(rng, 1),
        4 => {
            let n = rng.range(1, 4);
            format!("{}{}", reference, bases(rng, n))
        }
        5 => reference[..1].to_string(),
        6 => bases(rng, reference.len() as u64),
        7 => {
            let n = rng.range(0, 3);
            format!("{}{}", &reference[..1], bases(rng, n))
        }
        8 => rng
            .pick(&["<DEL>", "<DUP:TANDEM>", "<INV>", "<CNV>", "<INS>"])
            .to_string(),
        9 => rng.pick(&["<*>", "<NON_REF>", "<X>"]).to_string(),
        10 => "*".into(),
        11 => format!("{}[chr1:{}[", &reference[..1], rng.range(1, 999)),
        12 => format!("]chr2:{}]{}", rng.range(1, 999), &reference[..1]),
        _ => format!("{}{}", reference.to_ascii_lowercase(), bases(rng, 1)),
    }
}

fn genotype(rng: &mut Rng, n_alleles: u64, version44: bool) -> String {
    let ploidy = *rng.pick(&[1, 2, 2, 2, 3]);
    let mut s = String::new();
    if version44 && rng.chance(20) {
        s.push(*rng.pick(&['|', '/']));
    }
    for i in 0..ploidy {
        if i > 0 {
            s.push(*rng.pick(&['|', '/', '|']));
        }
        if rng.chance(10) {
            s.push('.');
        } else {
            s.push_str(&rng.range(0, n_alleles - 1).to_string());
        }
    }
    s
}

/// A random VCF covering the cases htslib's record decoding depends on:
/// versions, undeclared contigs and keys, END/SVLEN/LEN of various types,
/// symbolic and breakend alleles, ploidy, phasing and missing values.
pub fn random_vcf(rng: &mut Rng) -> String {
    let version = *rng.pick(&["4.2", "4.3", "4.4", "4.5"]);
    let v44 = version >= "4.4";
    let n_samples = rng.range(0, 4);
    let mut v = format!("##fileformat=VCFv{version}\n");
    let declared_contigs = rng.range(0, 3);
    for i in 0..declared_contigs {
        v += &format!("##contig=<ID=chr{},length={}>\n", i + 1, 1_000_000);
    }
    let int_or = |rng: &mut Rng| *rng.pick(&["Integer", "Integer", "Integer", "String", "Float"]);
    if rng.chance(85) {
        v += &format!(
            "##INFO=<ID=END,Number=1,Type={},Description=\"End\">\n",
            int_or(rng)
        );
    }
    if rng.chance(70) {
        v += &format!(
            "##INFO=<ID=SVLEN,Number=A,Type={},Description=\"x\">\n",
            int_or(rng)
        );
    }
    v += "##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Depth\">\n";
    if rng.chance(90) {
        v += "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n";
    }
    if rng.chance(60) {
        v += &format!(
            "##FORMAT=<ID=LEN,Number=1,Type={},Description=\"x\">\n",
            int_or(rng)
        );
    }
    v += "##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"Depth\">\n";
    v += "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT";
    for s in 0..n_samples {
        v += &format!("\tS{s}");
    }
    v.truncate(if n_samples == 0 {
        v.len() - "\tFORMAT".len()
    } else {
        v.len()
    });
    v.push('\n');

    for _ in 0..rng.range(1, 60) {
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
        let mut info = Vec::new();
        if rng.chance(30) {
            let end = if rng.chance(80) {
                pos + rng.range(0, 500)
            } else {
                pos.saturating_sub(5)
            };
            info.push(match rng.range(0, 9) {
                0 => "END=.".to_string(),
                1 => "END".to_string(),
                2 => format!("END={end},{}", end + 1),
                _ => format!("END={end}"),
            });
        }
        if rng.chance(30) {
            let vals: Vec<String> = (0..rng.range(1, 3))
                .map(|_| match rng.range(0, 5) {
                    0 => ".".to_string(),
                    1 => format!("-{}", rng.range(1, 900)),
                    _ => rng.range(0, 900).to_string(),
                })
                .collect();
            info.push(format!("SVLEN={}", vals.join(",")));
        }
        if rng.chance(30) {
            info.push(format!("DP={}", rng.range(0, 99)));
        }
        if rng.chance(10) {
            info.push("XUNDECLARED=1".into());
        }
        let info = if info.is_empty() {
            ".".to_string()
        } else {
            info.join(";")
        };
        let mut line = format!("{chrom}\t{pos}\t.\t{reference}\t{alt_field}\t50\tPASS\t{info}");
        if n_samples > 0 {
            let keys: Vec<&str> = match rng.range(0, 5) {
                0 => vec!["GT"],
                1 => vec!["GT", "DP"],
                2 => vec!["DP", "GT"],
                3 => vec!["GT", "LEN"],
                _ => vec!["LEN", "GT", "DP"],
            };
            line += &format!("\t{}", keys.join(":"));
            for _ in 0..n_samples {
                let mut fields: Vec<String> = keys
                    .iter()
                    .map(|k| match *k {
                        "GT" => genotype(rng, n_alt + 1, v44),
                        "LEN" => {
                            if rng.chance(20) {
                                ".".into()
                            } else {
                                rng.range(0, 700).to_string()
                            }
                        }
                        _ => rng.range(0, 50).to_string(),
                    })
                    .collect();
                // Trailing fields may be omitted (keeping GT: htslib exits
                // on a GT field no sample gives).
                if keys[0] == "GT" && rng.chance(10) {
                    fields.truncate(rng.range(1, fields.len() as u64) as usize);
                }
                line += &format!("\t{}", fields.join(":"));
            }
        }
        v += &line;
        v.push('\n');
    }
    v
}
