//! VCF generators for exact BCF byte comparisons.
use crate::support::variants::Rng;

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
pub fn random_vcf(rng: &mut Rng) -> String {
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
