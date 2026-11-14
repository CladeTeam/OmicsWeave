//! Consensus sequences of random regions, variants and options against
//! `bcftools consensus` (run in process on the same inputs): sequences and
//! chains must be byte-identical, and an input bcftools rejects must be
//! rejected. Each region is loaded once and queried with several allele
//! choices, as a batch over samples would.

mod common;

use std::collections::HashMap;
use std::sync::Arc;

use common::{Genome, Outcome, Rng, bcftools_consensus, fasta_sequence, to_indexed_bcf};
use omicsweave_consensus::{
    Alleles, Consensus, ConsensusOptions, Mark, Mask, MaskWith, Overlap, RecordFilter, VcfStore,
};
use omicsweave_core::Region;
use proptest::prelude::*;

const BASES: &[u8] = b"ACGT";

fn random_genome(rng: &mut Rng) -> Vec<(String, Vec<u8>)> {
    (1..=rng.range(1, 2))
        .map(|i| {
            let len = rng.range(80, 600) as usize;
            let mut seq: Vec<u8> = (0..len).map(|_| *rng.pick(BASES)).collect();
            // Soft-masked runs and Ns, as reference FASTAs have.
            for _ in 0..rng.range(0, 3) {
                let start = rng.range(0, len as u64 - 1) as usize;
                let end = (start + rng.range(1, 30) as usize).min(len);
                let lower = rng.chance(70);
                for b in &mut seq[start..end] {
                    *b = if lower { b.to_ascii_lowercase() } else { b'N' };
                }
            }
            (format!("chr{i}"), seq)
        })
        .collect()
}

fn other_base(rng: &mut Rng, b: u8) -> u8 {
    loop {
        let c = *rng.pick(BASES);
        if c != b.to_ascii_uppercase() {
            return c;
        }
    }
}

fn random_bases(rng: &mut Rng, n: u64) -> String {
    (0..n).map(|_| *rng.pick(BASES) as char).collect()
}

/// One ALT allele for `reference` (the REF bases, uppercase).
fn random_alt(rng: &mut Rng, reference: &[u8]) -> String {
    let r = String::from_utf8_lossy(reference).into_owned();
    match rng.range(0, 11) {
        0..=2 if reference.len() == 1 => (other_base(rng, reference[0]) as char).to_string(),
        0..=2 => {
            let mut alt = reference.to_vec();
            let i = rng.range(0, alt.len() as u64 - 1) as usize;
            alt[i] = other_base(rng, alt[i]);
            String::from_utf8(alt).unwrap()
        }
        3 | 4 => {
            let n = rng.range(1, 5);
            format!("{}{}", &r[..1], random_bases(rng, n))
        }
        5 | 6 => r[..1].to_string(),
        7 => {
            let n = rng.range(0, 3);
            format!("{}{}", random_bases(rng, 1), random_bases(rng, n))
        }
        8 => "<DEL>".into(),
        9 => "*".into(),
        10 => "<*>".into(),
        _ => format!("{}{}", r, random_bases(rng, 1)),
    }
}

fn genotype(rng: &mut Rng, n_alleles: u64) -> String {
    let ploidy = *rng.pick(&[1, 2, 2, 2]);
    let sep = *rng.pick(&['|', '|', '/']);
    (0..ploidy)
        .map(|_| {
            if rng.chance(8) {
                ".".to_string()
            } else {
                rng.range(0, n_alleles - 1).to_string()
            }
        })
        .collect::<Vec<_>>()
        .join(&sep.to_string())
}

/// A VCF of random records over `genome` (REFs read from it, with a few
/// mismatches), sorted by position as an indexed file must be.
fn random_vcf(rng: &mut Rng, genome: &[(String, Vec<u8>)], samples: &[String]) -> String {
    let mut v = String::from("##fileformat=VCFv4.2\n");
    for (name, seq) in genome {
        v += &format!("##contig=<ID={name},length={}>\n", seq.len());
    }
    v += "##INFO=<ID=END,Number=1,Type=Integer,Description=\"End\">\n";
    v += "##INFO=<ID=SVLEN,Number=A,Type=Integer,Description=\"Length\">\n";
    v += "##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Depth\">\n";
    v += "##FILTER=<ID=q10,Description=\"Low quality\">\n";
    v += "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n";
    v += "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT";
    for s in samples {
        v += &format!("\t{s}");
    }
    v.push('\n');
    for (name, seq) in genome {
        let mut pos = rng.range(1, 10) as usize;
        while pos < seq.len() {
            let len = (*rng.pick(&[1, 1, 1, 2, 3, 4]) as usize).min(seq.len() - pos + 1);
            let mut reference: Vec<u8> = seq[pos - 1..pos - 1 + len]
                .iter()
                .map(u8::to_ascii_uppercase)
                .collect();
            if rng.chance(2) {
                reference[0] = other_base(rng, reference[0]);
            }
            let n_alt = rng.range(1, 3);
            let alts: Vec<String> = (0..n_alt).map(|_| random_alt(rng, &reference)).collect();
            let mut info = Vec::new();
            if rng.chance(70) {
                info.push(format!("DP={}", rng.range(0, 60)));
            }
            if alts.iter().any(|a| a.starts_with('<')) && rng.chance(60) {
                let end = (pos + rng.range(0, 15) as usize).min(seq.len());
                info.push(format!("END={end}"));
            }
            let info = if info.is_empty() {
                ".".into()
            } else {
                info.join(";")
            };
            let id = if rng.chance(50) {
                format!("rs{}", rng.range(1, 5))
            } else {
                ".".into()
            };
            let qual = if rng.chance(80) {
                rng.range(0, 60).to_string()
            } else {
                ".".into()
            };
            let filter = *rng.pick(&["PASS", "PASS", "q10", "."]);
            let mut line = format!(
                "{name}\t{pos}\t{id}\t{}\t{}\t{qual}\t{filter}\t{info}\tGT",
                String::from_utf8_lossy(&reference),
                alts.join(",")
            );
            for _ in samples {
                line += &format!("\t{}", genotype(rng, n_alt + 1));
            }
            v += &line;
            v.push('\n');
            // Close steps make overlapping records.
            pos += rng.range(1, 25) as usize;
        }
    }
    v
}

/// A random `-i`/`-e` expression over the fields the VCF has, now and
/// then one that does not compile.
fn random_filter(rng: &mut Rng) -> String {
    let atom = |rng: &mut Rng| -> String {
        let n = rng.range(0, 60);
        match rng.range(0, 21) {
            0 => format!("QUAL>{n}"),
            1 => format!("INFO/DP>{n}"),
            2 => format!("DP<{n}"),
            3 => rng
                .pick(&[
                    "TYPE=\"snp\"",
                    "TYPE=\"indel\"",
                    "TYPE=\"mnp\"",
                    "TYPE=\"other\"",
                ])
                .to_string(),
            4 => rng
                .pick(&[
                    "GT=\"het\"",
                    "GT=\"alt\"",
                    "GT=\"hom\"",
                    "GT=\"mis\"",
                    "GT=\"ref\"",
                ])
                .to_string(),
            5 => format!("N_ALT>{}", rng.range(0, 2)),
            6 => "ALT~\"<\"".into(),
            7 => format!("ILEN{}0", rng.pick(&["<", ">", "="])),
            8 => format!("POS<{}", rng.range(1, 300)),
            9 => rng
                .pick(&["FILTER=\"PASS\"", "FILTER=\"q10\"", "FILTER=\".\""])
                .to_string(),
            10 => format!("AC>{}", rng.range(0, 3)),
            11 => format!("AF<0.{}", rng.range(1, 9)),
            12 => format!("MAF>0.{}", rng.range(1, 4)),
            13 => format!("strlen(REF)>{}", rng.range(1, 3)),
            14 => format!("ID=\"rs{}\"", rng.range(1, 5)),
            15 => format!("REF=\"{}\"", *rng.pick(&["A", "C", "G", "T"])),
            16 => format!("INFO/END>{n}"),
            17 => "N_MISSING>0".into(),
            18 => format!("F_MISSING>0.{}", rng.range(1, 9)),
            19 => "QUAL=\".\"".into(),
            _ => rng
                .pick(&["FOO>1", "QUAL>", "(QUAL>1", "INFO/XX=1"])
                .to_string(),
        }
    };
    let mut expr = atom(rng);
    for _ in 0..rng.range(0, 2) {
        let op = *rng.pick(&["&&", "||", "&", "|"]);
        expr = if rng.chance(30) {
            format!("({expr}) {op} {}", atom(rng))
        } else {
            format!("{expr} {op} {}", atom(rng))
        };
    }
    expr
}

#[derive(Debug)]
enum Mode {
    AllAlt,
    Samples(Vec<String>),
    Haplotype(String, String),
}

fn random_mode(rng: &mut Rng, samples: &[String]) -> Mode {
    match rng.range(0, 2) {
        0 => Mode::AllAlt,
        1 => {
            let mut s: Vec<String> = samples.iter().filter(|_| rng.chance(60)).cloned().collect();
            if s.is_empty() {
                s.push(samples[0].clone());
            }
            Mode::Samples(s)
        }
        _ => Mode::Haplotype(
            rng.pick(samples).clone(),
            rng.pick(&[
                "1", "2", "R", "A", "I", "LR", "LA", "SR", "SA", "1pIu", "2pIu",
            ])
            .to_string(),
        ),
    }
}

/// One random case; `Err` describes the first difference from bcftools.
fn run(seed: u64) -> Result<(), String> {
    let mut rng = Rng(seed);
    let genome = random_genome(&mut rng);
    let samples: Vec<String> = (1..=rng.range(1, 3)).map(|i| format!("S{i}")).collect();
    let vcf_text = random_vcf(&mut rng, &genome, &samples);

    let dir = tempfile::tempdir().unwrap();
    let path = |name: &str| dir.path().join(name).to_str().unwrap().to_owned();
    std::fs::write(path("in.vcf"), &vcf_text).unwrap();
    to_indexed_bcf(
        dir.path().join("in.vcf").as_path(),
        dir.path().join("in.bcf").as_path(),
    );

    // Options shared by every query of this case.
    let mut options = ConsensusOptions::default();
    let mut args: Vec<String> = Vec::new();
    let flag = |args: &mut Vec<String>, k: &str, v: String| {
        args.push(k.to_owned());
        args.push(v);
    };
    if rng.chance(20) {
        options.missing = Some(b'?');
        flag(&mut args, "-M", "?".into());
    }
    if rng.chance(15) {
        options.absent = Some(b'.');
        flag(&mut args, "-a", ".".into());
    }
    if rng.chance(15) {
        options.mark_del = Some(b'-');
        flag(&mut args, "--mark-del", "-".into());
    }
    for (field, name) in [(0, "--mark-ins"), (1, "--mark-snv")] {
        if rng.chance(15) {
            let (mark, text) = *rng.pick(&[
                (Mark::Upper, "uc"),
                (Mark::Lower, "lc"),
                (Mark::Char(b'x'), "x"),
            ]);
            if field == 0 {
                options.mark_ins = Some(mark);
            } else {
                options.mark_snv = Some(mark);
            }
            flag(&mut args, name, text.into());
        }
    }
    if rng.chance(20) {
        let mut bed = String::new();
        for (name, seq) in &genome {
            for _ in 0..rng.range(1, 3) {
                let start = rng.range(0, seq.len() as u64 - 2);
                let end = (start + rng.range(1, 20)).min(seq.len() as u64);
                bed += &format!("{name}\t{start}\t{end}\n");
            }
        }
        std::fs::write(path("mask.bed"), &bed).unwrap();
        let (with, text) = *rng.pick(&[
            (MaskWith::Char(b'N'), "N"),
            (MaskWith::Upper, "uc"),
            (MaskWith::Lower, "lc"),
        ]);
        options.mask = Some(Arc::new(Mask::parse(bed.as_bytes(), with).unwrap()));
        flag(&mut args, "-m", path("mask.bed"));
        flag(&mut args, "--mask-with", text.into());
    }
    let (overlap, text) = *rng.pick(&[
        (Overlap::Start, "0"),
        (Overlap::Record, "1"),
        (Overlap::Variant, "2"),
    ]);
    options.regions_overlap = overlap;
    flag(&mut args, "--regions-overlap", text.into());
    options.chain = rng.chance(20);
    if rng.chance(25) {
        options.iupac_codes = true;
        args.push("-I".into());
    }
    let expr = rng.chance(30).then(|| random_filter(&mut rng));
    let include = rng.chance(50);
    let filter = expr.as_deref().map(|e| {
        if include {
            flag(&mut args, "-i", e.into());
            RecordFilter::Include(e)
        } else {
            flag(&mut args, "-e", e.into());
            RecordFilter::Exclude(e)
        }
    });

    // Both read the same BCF: htslib reads genotype phasing of BCF and VCF
    // text differently before VCF 4.4.
    let bcf = std::fs::File::open(path("in.bcf")).unwrap();
    let store = match VcfStore::read(bcf, filter) {
        Ok(store) => Arc::new(store),
        // An expression that does not compile: bcftools must reject it too.
        Err(e) => {
            let mut run_args = args.clone();
            flag(&mut run_args, "-f", path("in.vcf"));
            flag(&mut run_args, "-o", path("out.fa"));
            run_args.push(path("in.bcf"));
            return match bcftools_consensus(&run_args) {
                Outcome::Failed | Outcome::Crashed(_) | Outcome::Done(1..) => Ok(()),
                theirs => Err(format!("filter {expr:?}: ours {e}, bcftools {theirs:?}")),
            };
        }
    };
    let mut consensus = Consensus::new(store, options.clone());
    let mut reference = Genome(genome.iter().cloned().collect::<HashMap<_, _>>());

    for _ in 0..4 {
        let (name, seq) = rng.pick(&genome).clone();
        let start = rng.range(0, seq.len() as u64 - 1) as u32;
        let end = (start + rng.range(1, 200) as u32).min(seq.len() as u32);
        std::fs::write(
            path("ref.fa"),
            format!(
                ">{name}:{}-{end}\n{}\n",
                start + 1,
                String::from_utf8_lossy(&seq[start as usize..end as usize])
            ),
        )
        .unwrap();
        consensus
            .load(&mut reference, Region::new(&name, start, end))
            .unwrap();
        for _ in 0..3 {
            let mode = random_mode(&mut rng, &samples);
            let mut run_args = args.clone();
            match &mode {
                Mode::AllAlt => flag(&mut run_args, "-s", "-".into()),
                Mode::Samples(s) => flag(&mut run_args, "-s", s.join(",")),
                Mode::Haplotype(s, h) => {
                    flag(&mut run_args, "-s", s.clone());
                    flag(&mut run_args, "-H", h.clone());
                }
            }
            let _ = std::fs::remove_file(path("out.fa"));
            let _ = std::fs::remove_file(path("out.chain"));
            flag(&mut run_args, "-f", path("ref.fa"));
            flag(&mut run_args, "-o", path("out.fa"));
            if options.chain {
                flag(&mut run_args, "-c", path("out.chain"));
            }
            run_args.push(path("in.bcf"));
            let theirs = bcftools_consensus(&run_args);

            let names: Vec<&str> = match &mode {
                Mode::Samples(s) => s.iter().map(String::as_str).collect(),
                _ => Vec::new(),
            };
            let alleles = match &mode {
                Mode::AllAlt => Alleles::AllAlt,
                Mode::Samples(_) => Alleles::Samples(&names),
                Mode::Haplotype(s, h) => Alleles::Haplotype {
                    sample: s,
                    haplotype: h,
                },
            };
            let ours = consensus.query(alleles);
            let what = format!("{mode:?} {args:?} {name}:{start}-{end}\n{vcf_text}");
            match (theirs, ours) {
                (Outcome::Done(0), Ok(n)) => {
                    let mut seq = vec![0u8; n];
                    consensus.write(&mut seq).unwrap();
                    let want = fasta_sequence(&std::fs::read_to_string(path("out.fa")).unwrap());
                    if seq != want {
                        return Err(format!(
                            "sequence\n ours: {}\n bcftools: {}\n{what}",
                            String::from_utf8_lossy(&seq),
                            String::from_utf8_lossy(&want)
                        ));
                    }
                    if options.chain {
                        let want = std::fs::read_to_string(path("out.chain")).unwrap();
                        if consensus.chain().unwrap() != want {
                            return Err(format!(
                                "chain\n ours: {}\n bcftools: {want}\n{what}",
                                consensus.chain().unwrap()
                            ));
                        }
                    }
                }
                (Outcome::Failed, Err(_)) => {}
                // bcftools' own memory errors, and htslib exiting on what
                // bcftools passes it from them: nothing to compare against.
                (Outcome::Crashed(_) | Outcome::Done(1..), _) => {}
                (theirs, ours) => {
                    return Err(format!(
                        "outcome: bcftools {theirs:?}, ours {:?}\n{what}",
                        ours.map(|_| ())
                    ));
                }
            }
        }
    }
    Ok(())
}

proptest! {
    #![proptest_config(ProptestConfig { cases: 300, ..ProptestConfig::default() })]

    #[test]
    fn sequences_match_bcftools(seed in any::<u64>()) {
        if let Err(e) = run(seed) {
            panic!("{e}");
        }
    }
}
