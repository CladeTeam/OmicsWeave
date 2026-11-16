//! Shared expression grammar over either valid or deliberately malformed VCFs.
use crate::support::rich_variants;
use crate::support::rich_variants::{Decl, Spec};
use crate::support::variants::Rng;
use std::path::Path;

/// Files `@name` expressions refer to.
struct Files {
    samples: String,
    ids: String,
    values: String,
}

struct Gen<'a> {
    rng: &'a mut Rng,
    spec: &'a Spec,
    files: &'a Files,
}

const NUM_CMP: &[&str] = &["=", "==", "!=", "<", "<=", ">", ">=", "=", ">"];
const STR_CMP: &[&str] = &["=", "==", "!=", "~", "!~"];
const LOGIC: &[&str] = &["&&", "||", "&", "|", "&&", "||"];
const ARITH: &[&str] = &["+", "-", "*", "/", "%"];

impl Gen<'_> {
    fn pick<'b>(&mut self, items: &'b [&'b str]) -> &'b str {
        let item: &&str = self.rng.pick(items);
        item
    }

    fn chance(&mut self, p: u64) -> bool {
        self.rng.chance(p)
    }

    fn decls(&self, format: bool, kinds: &[&str]) -> Vec<Decl> {
        let list = if format {
            &self.spec.format
        } else {
            &self.spec.info
        };
        list.iter()
            .filter(|d| kinds.contains(&d.ty) && d.name != "GT")
            .cloned()
            .collect()
    }

    /// An INFO tag reference with an optional prefix and subscript.
    fn info_ref(&mut self, kinds: &[&str]) -> String {
        let decls = self.decls(false, kinds);
        let name = if decls.is_empty() || self.chance(6) {
            let mut names: Vec<String> = self.spec.undeclared_info.clone();
            names.extend(["NOPE".into(), "DP".into(), "AC".into(), "AF".into()]);
            self.rng.pick(&names).clone()
        } else {
            self.rng.pick(&decls).name.clone()
        };
        let prefix = self.pick(&["INFO/", "INFO/", "", "info/"]);
        let sub = if self.chance(30) {
            self.pick(&[
                "[0]", "[1]", "[*]", "[0,1]", "[1-]", "[0-1]", "[2]", "[]", "[0:1]", "[GT]",
            ])
        } else {
            ""
        };
        format!("{prefix}{name}{sub}")
    }

    /// A FORMAT tag reference with an optional prefix and subscript.
    fn fmt_ref(&mut self, kinds: &[&str]) -> String {
        let decls = self.decls(true, kinds);
        let name = if decls.is_empty() || self.chance(5) {
            let mut names: Vec<String> = self.spec.undeclared_format.clone();
            names.extend(["NOPE".into(), "DP".into(), "AD".into()]);
            self.rng.pick(&names).clone()
        } else {
            self.rng.pick(&decls).name.clone()
        };
        let prefix = self.pick(&["FMT/", "FORMAT/", "FMT/", "", "fmt/"]);
        let sample_list = format!("[@{}]", self.files.samples);
        let sample_list2 = format!("[@{}:1]", self.files.samples);
        let sub = if self.chance(45) {
            match self.rng.range(0, 14) {
                0 => sample_list,
                1 => sample_list2,
                _ => self
                    .pick(&[
                        "[0]", "[1]", "[*]", "[:0]", "[0:0]", "[*:1]", "[0:*]", "[1:0-1]", "[:GT]",
                        "[0,1:0]", "[1-:1]", "[5]", "[0:1,2]", "[GT]",
                    ])
                    .to_string(),
            }
        } else {
            String::new()
        };
        format!("{prefix}{name}{sub}")
    }

    fn number(&mut self) -> String {
        self.pick(&[
            "0",
            "1",
            "2",
            "5",
            "10",
            "30",
            "-1",
            "0.5",
            "1.5",
            "59.2",
            "100",
            "1e3",
            ".5",
            "1e-3",
            "16777217",
            "3000000000",
            "2.5e1",
            "0x10",
            "1e400",
            "inf",
            "nan",
            "3",
            "0.25",
        ])
        .to_string()
    }

    fn numeric_atom(&mut self) -> String {
        match self.rng.range(0, 24) {
            0..=1 => "QUAL".into(),
            2 => "POS".into(),
            3 => self
                .pick(&["N_ALT", "N_SAMPLES", "N_MISSING", "F_MISSING"])
                .into(),
            4 => self
                .pick(&["AC", "AN", "AF", "MAC", "MAF", "AC[0]", "AF[1]"])
                .into(),
            5 => self.pick(&["ILEN", "%ILEN", "ILEN[0]"]).into(),
            6..=11 => self.info_ref(&["Integer", "Float", "Flag"]),
            12..=16 => self.fmt_ref(&["Integer", "Float"]),
            17 => "\".\"".into(),
            18 => self
                .pick(&["TYPE", "%QUAL", "q", "N_A", "INFO/END", "SVLEN"])
                .into(),
            _ => self.number(),
        }
    }

    fn function(&mut self, arg: String) -> String {
        let f = self.pick(&[
            "MAX", "MIN", "AVG", "MEAN", "MEDIAN", "STDEV", "SUM", "ABS", "COUNT", "PHRED", "sMAX",
            "sMIN", "sAVG", "sMEDIAN", "sSTDEV", "sSUM", "sCOUNT", "SMPL_MAX", "smpl_sum", "%MAX",
            "max",
        ]);
        format!("{f}({arg})")
    }

    fn numeric(&mut self, depth: u32) -> String {
        if depth > 2 {
            return self.numeric_atom();
        }
        match self.rng.range(0, 20) {
            0..=10 => self.numeric_atom(),
            11..=13 => {
                let a = self.numeric(depth + 1);
                let b = self.numeric(depth + 1);
                let op = self.pick(ARITH);
                format!("{a}{op}{b}")
            }
            14..=16 => {
                let a = self.numeric(depth + 1);
                self.function(a)
            }
            17 => {
                let a = self.numeric(depth + 1);
                format!("-{a}")
            }
            18 => match self.rng.range(0, 5) {
                0 => format!("binom({})", self.fmt_ref(&["Integer"])),
                1 => format!(
                    "binom({},{})",
                    self.info_ref(&["Integer"]),
                    self.info_ref(&["Integer"])
                ),
                2 => format!(
                    "binom({},{})",
                    self.fmt_ref(&["Integer"]),
                    self.fmt_ref(&["Integer"])
                ),
                3 => format!("fisher({})", self.info_ref(&["Integer"])),
                _ => format!(
                    "fisher({},{})",
                    self.fmt_ref(&["Integer"]),
                    self.fmt_ref(&["Integer"])
                ),
            },
            _ => {
                let s = self.string_term();
                format!("STRLEN({s})")
            }
        }
    }

    fn string_term(&mut self) -> String {
        match self.rng.range(0, 10) {
            0 => "CHROM".into(),
            1 => "REF".into(),
            2 => self
                .pick(&["ALT", "ALT[0]", "ALT[1]", "ALT[*]", "ALT[0,1]"])
                .into(),
            3..=6 => self.info_ref(&["String", "Character"]),
            _ => self.fmt_ref(&["String", "Character"]),
        }
    }

    fn regex(&mut self) -> String {
        let mut s = String::new();
        if self.chance(25) {
            s.push('^');
        }
        for _ in 0..self.rng.range(1, 4) {
            let piece = self.pick(&[
                "a",
                "b",
                "c",
                "A",
                "X",
                "1",
                ".",
                "[abc]",
                "[^a]",
                "[a-z]",
                "[[:digit:]]",
                "[[:upper:]]",
                "\\(a\\)",
                "x\\|y",
                "\\w",
                "\\b",
                "rs",
                "_",
                "-",
                ",",
                "[",
            ]);
            s.push_str(piece);
            if self.chance(25) {
                s.push_str(self.pick(&["*", "\\+", "\\?", "\\{1,2\\}", "\\{2\\}", "**"]));
            }
        }
        if self.chance(20) {
            s.push('$');
        }
        if self.chance(25) {
            s.push_str("/i");
        }
        s
    }

    fn string_value(&mut self) -> String {
        match self.rng.range(0, 10) {
            0 => ".".into(),
            1..=4 => rich_variants::string_value(self.rng),
            5 => self
                .pick(&[
                    "chr1", "chr2", "A", "C", "AT", "<DEL>", "*", "C,G", "a,.", "abc,x",
                ])
                .into(),
            _ => {
                let n = self.rng.range(1, 3);
                (0..n)
                    .map(|_| *self.rng.pick(&['a', 'b', 'c', 'A', 'X', '1', '.']))
                    .collect()
            }
        }
    }

    fn comparison(&mut self, depth: u32) -> String {
        match self.rng.range(0, 30) {
            0..=10 => {
                let a = self.numeric(depth);
                let b = self.numeric(depth + 1);
                let op = self.pick(NUM_CMP);
                format!("{a}{op}{b}")
            }
            11..=15 => {
                let a = self.string_term();
                let op = self.pick(STR_CMP);
                let v = if op.contains('~') {
                    self.regex()
                } else {
                    self.string_value()
                };
                if self.chance(15) {
                    format!("\"{v}\"{op}{a}")
                } else {
                    format!("{a}{op}\"{v}\"")
                }
            }
            16..=18 => {
                let gt = self.pick(&["GT", "FMT/GT", "GT[0]", "GT[*]", "FMT/GT[1:0]"]);
                let v = self.pick(&[
                    "het", "hom", "hap", "mis", "ref", "alt", "rr", "ra", "ar", "aa", "aA", "Aa",
                    "r", "a", "HET", "0/1", "1|0", "./.", "0", ".", "1/1",
                ]);
                let op = self.pick(&["=", "!=", "==", "~", "!~"]);
                format!("{gt}{op}\"{v}\"")
            }
            19..=20 => {
                let mut f: Vec<String> = Vec::new();
                for _ in 0..self.rng.range(1, 2) {
                    f.push(self.rng.pick(&self.spec.filters).clone());
                }
                let v = match self.rng.range(0, 8) {
                    0 => ".".to_string(),
                    1 => "XUNDEF".to_string(),
                    2 => "nope".to_string(),
                    _ => f.join(";"),
                };
                let op = self.pick(&["=", "!=", "~", "!~", "=="]);
                if self.chance(15) {
                    format!("\"{v}\"{op}FILTER")
                } else {
                    format!("FILTER{op}\"{v}\"")
                }
            }
            21..=22 => {
                let v = self.pick(&[
                    "snp", "snps", "indel", "mnp", "other", "bnd", "overlap", "ref", "SNP", "foo",
                ]);
                let op = self.pick(&["=", "!=", "~", "!~"]);
                format!("TYPE{op}\"{v}\"")
            }
            23..=25 => {
                let op = self.pick(&["=", "!=", "~", "!~"]);
                match self.rng.range(0, 4) {
                    0 => format!("ID{}@{}", self.pick(&["=", "!=", "=="]), self.files.ids),
                    1 => {
                        let r = self.regex();
                        format!("ID{op}\"{r}\"")
                    }
                    _ => {
                        let id = self.pick(&["rs1", "rs2", "rs12", ".", "rs3;rs4"]);
                        format!("ID{op}\"{id}\"")
                    }
                }
            }
            26 => {
                let s = self.string_term();
                let op = self.pick(&["=", "!="]);
                format!("{s}{op}@{}", self.files.values)
            }
            27 => {
                let inner = self.comparison(depth + 1);
                let f = self.pick(&["N_PASS", "F_PASS"]);
                let op = self.pick(NUM_CMP);
                let n = self.number();
                format!("{f}({inner}){op}{n}")
            }
            28 => self.numeric(depth),
            _ => self.info_ref(&["Flag"]),
        }
    }

    fn expr(&mut self, depth: u32) -> String {
        if depth > 2 {
            return self.comparison(depth);
        }
        match self.rng.range(0, 10) {
            0..=4 => self.comparison(depth),
            5..=7 => {
                let a = self.expr(depth + 1);
                let b = self.expr(depth + 1);
                let op = self.pick(LOGIC);
                format!("{a} {op} {b}")
            }
            8 => {
                let a = self.expr(depth + 1);
                format!("({a})")
            }
            _ => {
                let a = self.expr(depth + 1);
                let b = self.expr(depth + 1);
                let op = self.pick(LOGIC);
                format!("({a}){op}({b})")
            }
        }
    }

    /// An expression, occasionally damaged to exercise parse errors.
    fn expression(&mut self) -> String {
        if self.chance(6) {
            return self
                .pick(&[
                    "QUAL={}",
                    "{}=QUAL",
                    "QUAL={int} | CHROM={str}",
                    "PERL.f(QUAL)",
                    "perl:x.pl; QUAL>1",
                    "N_PASS(QUAL>1)",
                    "F_PASS(FMT/DP>1)>0.5",
                    "COUNT(FMT/DP>1)",
                    "sCOUNT(FMT/DP)",
                    "MAX(FMT/AD[0:*])",
                    "AVG(INFO/AF)",
                    "binom()",
                    "fisher()",
                    "GT=\"het\" && FMT/DP>5",
                    "QUAL>30 || FMT/GQ>30",
                    "-QUAL<-5",
                    "--QUAL",
                    "QUAL",
                    "1",
                    "\"a\"",
                    "()",
                    "(",
                    "QUAL<=",
                    "QUAL=<5",
                    "QUAL>=5)",
                    "ILEN>0",
                    "N_ALT>=2",
                    "%FILTER=\"PASS\"",
                    "%ID~\"rs\"",
                    "CHROM~\"chr[12]\"",
                    "ALT~\"<\"",
                    "REF=\"A\" & ALT=\"C\"",
                    "TYPE=\"snp\" && TYPE!=\"indel\"",
                    "\"snp\"=TYPE",
                    "FILTER=\"PASS\" | FILTER=\".\"",
                    "N_MISSING>0 || F_MISSING>0.5",
                    "MAC>1 && MAF<0.5",
                    "AC[0]>1",
                    "INFO/DP4[0]+INFO/DP4[1]>5",
                    "fisher(INFO/DP4)<0.05",
                    "binom(FMT/AD)<0.05",
                    "FMT/AD[*:GT]>3",
                    "FMT/AD[0:GT]>3",
                    "smpl_median(FMT/AD)>1",
                    "MEDIAN(FMT/AD)>1",
                    "STDEV(FMT/AD)>1",
                    "PHRED(INFO/AF)>3",
                    "ABS(ILEN)>1",
                    "STRLEN(REF)>1",
                    "STRLEN(INFO/TXT)>1",
                    "QUAL % 2 == 1",
                    "POS/1000 > 500",
                    "1==1",
                    "\".\"=QUAL",
                    "QUAL!=\".\"",
                    "INFO/DP=\".\"",
                    "FMT/DP=\".\"",
                    "FMT/FT=\".\"",
                    "INFO/TXT!=\".\"",
                    "GT=\"mis\"",
                    "GT!~\"0\"",
                    "ID=\"rs1\" || ID=\"rs2\"",
                    "ID!=\".\"",
                    "%TYPE=\"ref\"",
                    "N_SAMPLES>1",
                    "q>10",
                    "QUAL>10 & DP>5",
                    "MIN(FMT/GQ)>20",
                    "SUM(FMT/AD)>10",
                    "FMT/DP>10 & FMT/GQ>20",
                    "FMT/DP>10 && FMT/GQ>20",
                    "FMT/DP>10 | FMT/GQ>20",
                    "FMT/DP>10 || FMT/GQ>20",
                ])
                .to_string();
        }
        let e = self.expr(0);
        if !self.chance(10) {
            return e;
        }
        let mut b: Vec<char> = e.chars().collect();
        let at = self.rng.range(0, b.len() as u64) as usize;
        match self.rng.range(0, 3) {
            0 if at < b.len() => {
                b.remove(at);
            }
            1 => {
                let c = *self
                    .rng
                    .pick(&['(', ')', '"', '[', ']', '=', '!', '&', '-', ' ', ',', '@']);
                b.insert(at.min(b.len()), c);
            }
            _ => {
                if at < b.len() {
                    let c = b[at];
                    b.insert(at, c);
                }
            }
        }
        b.into_iter().collect()
    }
}

pub fn generate(rng: &mut Rng, spec: &Spec, directory: &Path) -> Vec<String> {
    let samples = directory.join("samples.txt");
    let mut names: Vec<String> = spec
        .samples
        .iter()
        .filter(|_| rng.chance(60))
        .cloned()
        .collect();
    if rng.chance(10) {
        names.push("NOSUCH".into());
    }
    std::fs::write(&samples, names.join("\n") + "\n").unwrap();
    let ids = directory.join("ids.txt");
    std::fs::write(&ids, "rs1\nrs3 extra\nrs12\r\n\n.\n").unwrap();
    let values = directory.join("values.txt");
    std::fs::write(&values, "abc\nX\n.\nhet\nchr1\nA\n").unwrap();
    let files = Files {
        samples: samples.to_str().unwrap().into(),
        ids: ids.to_str().unwrap().into(),
        values: values.to_str().unwrap().into(),
    };
    let mut g = Gen {
        rng,
        spec,
        files: &files,
    };
    (0..6).map(|_| g.expression()).collect()
}
