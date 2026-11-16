//! Type-correct inputs and typed expression families for C differential tests.

/// A typed differential corpus. String fields begin at their largest layout,
/// so filter.c never reads the unwritten holes in a later-grown allocation.
/// The malformed corpus remains separate because a C process returning a bool
/// does not establish that its memory reads were defined.
pub fn defined_vcf(
    rng: &mut crate::support::variants::Rng,
) -> (String, crate::support::rich_variants::Spec) {
    use crate::support::rich_variants::{Decl, Spec};
    let declarations = |fields: &[(&str, &'static str, &'static str)]| {
        fields
            .iter()
            .map(|&(name, ty, number)| Decl {
                name: name.into(),
                ty,
                number,
            })
            .collect::<Vec<_>>()
    };
    let spec = Spec {
        info: declarations(&[
            ("DP", "Integer", "1"),
            ("AF", "Float", "A"),
            ("DP4", "Integer", "4"),
            ("AN", "Integer", "1"),
            ("AC", "Integer", "A"),
            ("XI", "Integer", "."),
            ("XF", "Float", "."),
            ("TXT", "String", "."),
            ("XS", "String", "."),
            ("CNT", "Character", "1"),
            ("DB", "Flag", "0"),
        ]),
        format: declarations(&[
            ("DP", "Integer", "."),
            ("GQ", "Integer", "."),
            ("AD", "Integer", "R"),
            ("ADR", "Integer", "R"),
            ("PL", "Integer", "G"),
            ("LEN", "Integer", "."),
            ("XF", "Float", "."),
            ("GL", "Float", "G"),
            ("FT", "String", "."),
            ("HQ", "String", "."),
        ]),
        filters: vec!["PASS".into(), "LowQual".into(), "s50".into()],
        samples: vec!["S0".into(), "S1".into(), "S2".into()],
        ..Spec::default()
    };
    let mut text = "##fileformat=VCFv4.4\n##contig=<ID=chr1,length=1000000>\n##contig=<ID=chr2,length=1000000>\n##FILTER=<ID=LowQual,Description=\"quality\">\n##FILTER=<ID=s50,Description=\"depth\">\n".to_owned();
    for (kind, declarations) in [("INFO", &spec.info), ("FORMAT", &spec.format)] {
        for field in declarations {
            text.push_str(&format!(
                "##{kind}=<ID={},Number={},Type={}>\n",
                field.name, field.number, field.ty
            ));
        }
    }
    text.push_str("##FORMAT=<ID=GT,Number=1,Type=String>\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS0\tS1\tS2\n");
    let largest_string = "aBc,XYZ,het,ref,bnd,foo,bar";
    for row in 0..24 {
        let mut genotypes = Vec::new();
        let mut counts = [0; 3];
        for _ in 0..3 {
            if row != 0 && rng.chance(15) {
                genotypes.push("./.".to_owned());
            } else {
                let a = rng.range(0, 2) as usize;
                let b = rng.range(0, 2) as usize;
                counts[a] += 1;
                counts[b] += 1;
                genotypes.push(format!("{a}{}{b}", if rng.chance(50) { '|' } else { '/' }));
            }
        }
        let mut integer_vector = |width: usize| {
            (0..width)
                .map(|_| {
                    if rng.chance(10) {
                        ".".into()
                    } else {
                        rng.range(0, 40).to_string()
                    }
                })
                .collect::<Vec<_>>()
                .join(",")
        };
        let dp4 = integer_vector(4);
        let xi = integer_vector(4);
        let depth = rng.range(0, 80);
        let af = format!(
            "{},{}",
            rng.range(0, 10) as f64 / 10.0,
            rng.range(0, 10) as f64 / 10.0
        );
        let xf = (0..4)
            .map(|_| {
                rng.pick(&[".", "0", "-0", "0.25", "1.5", "nan", "inf", "-inf"])
                    .to_string()
            })
            .collect::<Vec<_>>()
            .join(",");
        let info_string = if row == 0 {
            largest_string
        } else {
            rng.pick(&["abc,.,de", "aBc,,het", ".,.,xyz"])
        };
        let info = format!(
            "DP={depth};AF={af};DP4={dp4};AN={};AC={},{};XI={xi};XF={xf};TXT={info_string};XS={info_string};CNT=a{}",
            counts.iter().sum::<i32>(),
            counts[1],
            counts[2],
            if rng.chance(50) { ";DB" } else { "" }
        );
        let qual = rng.pick(&[".", "0", "12.5", "30", "99", "nan", "inf"]);
        let filter = &spec.filters[rng.range(0, 2) as usize];
        text.push_str(&format!(
            "chr{}\t{}\trs{}\tA\tC,G\t{qual}\t{filter}\t{info}\tGT:DP:GQ:AD:ADR:PL:LEN:XF:GL:FT:HQ",
            1 + row / 12,
            10 + row * 101,
            row % 13
        ));
        for (sample, gt) in genotypes.iter().enumerate() {
            let mut fields = vec![gt.clone()];
            for width in [4, 4, 3, 3, 6, 4] {
                fields.push(
                    (0..width)
                        .map(|_| {
                            if rng.chance(10) {
                                ".".into()
                            } else {
                                rng.range(0, 40).to_string()
                            }
                        })
                        .collect::<Vec<_>>()
                        .join(","),
                );
            }
            for width in [4, 6] {
                fields.push(
                    (0..width)
                        .map(|_| rng.pick(&[".", "0.25", "1.5", "2.5", "10"]).to_string())
                        .collect::<Vec<_>>()
                        .join(","),
                );
            }
            // The final sample is padded inside each initialized BCF payload,
            // so C string consumers also encounter a written terminator.
            let string = if sample == 2 {
                "."
            } else if row == 0 {
                largest_string
            } else {
                rng.pick(&["abc,.,de", "aBc,,het", ".,.,xyz"])
            };
            fields.extend([string.to_owned(), string.to_owned()]);
            text.push('\t');
            text.push_str(&fields.join(":"));
        }
        text.push('\n');
    }
    (text, spec)
}

/// Each family fixes operand kinds, sample masks and vector lengths before
/// randomizing operators/constants. Broad arbitrary syntax belongs to the
/// malformed-state suite, where C memory observations are not an oracle.
pub fn defined_expressions(
    rng: &mut crate::support::variants::Rng,
    directory: &std::path::Path,
) -> Vec<String> {
    let lookup = directory.join("lookup.txt");
    std::fs::write(&lookup, "rs0\nrs3\nrs7\n").unwrap();
    let samples = directory.join("samples.txt");
    std::fs::write(&samples, "S0\nS1\n").unwrap();
    let comparison =
        |rng: &mut crate::support::variants::Rng| *rng.pick(&["==", "!=", "<", "<=", ">", ">="]);
    let mut expressions = Vec::new();
    let op = comparison(rng);
    let offset = rng.range(0, 40);
    let scale = rng.range(1, 5);
    let threshold = rng.range(0, 100);
    expressions.push(format!("(INFO/DP+{offset})*{scale}{op}{threshold}"));

    let op = comparison(rng);
    expressions.push(match rng.range(0, 3) {
        0 => format!("INFO/XI+0.25{op}INFO/XF+0.5"),
        1 => format!("INFO/AF[0]{op}INFO/AF[1]"),
        2 => format!("INFO/XI[0,1]{op}INFO/XI[2]"),
        _ => format!("MAX(INFO/XI){op}{threshold}"),
    });

    let op = comparison(rng);
    expressions.push(match rng.range(0, 4) {
        0 => format!("FMT/DP{op}FMT/GQ"),
        1 => format!("FMT/DP[1:0-1]{op}FMT/GQ[1:0-1]"),
        2 => format!("sSUM(FMT/DP){op}sSUM(FMT/GQ)"),
        3 => format!(
            "FMT/DP[@{}:1]{op}FMT/GQ[@{}:1]",
            samples.display(),
            samples.display()
        ),
        _ => format!("FMT/AD[:GT]{op}{threshold}"),
    });

    expressions.push(match rng.range(0, 6) {
        0 => format!("INFO/TXT{}\"abc\"", rng.pick(&["=", "!=", "~", "!~"])),
        1 => format!("FMT/FT[:0]{}\"abc\"", rng.pick(&["=", "!=", "~", "!~"])),
        2 => "FMT/FT=INFO/TXT".into(),
        3 => "STRLEN(INFO/TXT)>2".into(),
        4 => format!("ID=@{}", lookup.display()),
        5 => format!("FMT/FT[@{}:0]!~\"^abc$\"", samples.display()),
        _ => "STRLEN(FMT/FT)>2".into(),
    });

    expressions.push(match rng.range(0, 3) {
        0 => format!("GT=\"het\" && FMT/DP[:0]>{offset}"),
        1 => "N_PASS(GT=\"het\")>=1".into(),
        2 => format!("F_PASS(FMT/DP[:0]>{offset})>0.5"),
        _ => "GT=\"mis\" || INFO/DB=1".into(),
    });

    expressions.push(match rng.range(0, 4) {
        0 => "fisher(INFO/DP4)<0.05".into(),
        1 => "binom(FMT/AD)<0.5".into(),
        2 => "fisher(FMT/AD[:0-1],FMT/ADR[:0-1])>=0".into(),
        3 => "binom(FMT/DP[:0],FMT/GQ[:0])>=0".into(),
        _ => "sMAX(FMT/DP)>2".into(),
    });
    for expression in &mut expressions {
        if rng.chance(35) {
            // Whitespace terminates an unquoted @file before the closing parenthesis.
            *expression = format!(
                "({expression} ) {} (QUAL>=0)",
                rng.pick(&["&&", "||", "&", "|"])
            );
        }
    }
    // A known absent field gives a deterministic compile-error comparison,
    // without relying on malformed parser internals or invalid operand shapes.
    expressions.push(
        rng.pick(&["INFO/NO_SUCH_FIELD>0", "FMT/NO_SUCH_FIELD>0"])
            .to_string(),
    );
    expressions
}
