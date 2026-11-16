use std::collections::HashMap;
use std::sync::Arc;

use crate::support::{
    Genome, Outcome, bcf_of, bcftools_consensus, consensus as open, fasta_sequence, to_indexed_bcf,
};
use omicsweave_consensus::{Alleles, Consensus, ConsensusOptions, Mark, Mask, MaskWith, Overlap};
use omicsweave_core::{Region, Strand};

const HEADER: &str = "##fileformat=VCFv4.3\n##contig=<ID=c,length=128>\n\
    ##INFO=<ID=END,Number=1,Type=Integer,Description=\"end\">\n\
    ##FORMAT=<ID=GT,Number=1,Type=String,Description=\"gt\">\n\
    #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS\n";

fn output(consensus: &Consensus, n: usize) -> Vec<u8> {
    let mut out = vec![0; n];
    consensus.write(Strand::Unstranded, &mut out).unwrap();
    out
}

fn genome() -> Genome {
    Genome(HashMap::from([("c".to_owned(), vec![b'A'; 128])]))
}

#[test]
fn empty_regions_inside_masks_are_valid_for_every_mask_mode() {
    let dir = tempfile::tempdir().unwrap();
    let bcf = bcf_of(dir.path(), "in", HEADER.as_bytes());
    for with in [MaskWith::Char(b'N'), MaskWith::Lower, MaskWith::Upper] {
        let mask = Arc::new(Mask::parse(b"c\t0\t100\n", with).unwrap());
        let mut consensus = open(
            &bcf,
            None,
            ConsensusOptions {
                mask: Some(mask),
                absent: Some(b'.'),
                chain: true,
                ..ConsensusOptions::default()
            },
        );
        for coordinate in [0, 5, 99, 100, 128] {
            consensus
                .load(&mut genome(), Region::new("c", coordinate, coordinate))
                .unwrap();
            let n = consensus.query(Alleles::AllAlt).unwrap();
            assert_eq!(n, 0, "{with:?}, {coordinate}");
            assert!(output(&consensus, n).is_empty());
            assert!(consensus.chain().is_some());
        }
    }
}

#[test]
fn targeted_valid_regions_masks_and_chains_require_successful_bcftools_oracle() {
    let dir = tempfile::tempdir().unwrap();
    let path = |name: &str| dir.path().join(name);
    let text = format!(
        "{HEADER}\
        c\t2\t.\tAAA\tA\t.\t.\t.\tGT\t1/0\n\
        c\t4\t.\tA\tATT\t.\t.\t.\tGT\t1/0\n\
        c\t10\t.\tAA\tCG\t.\t.\t.\tGT\t0/1\n\
        c\t17\t.\tA\tC\t.\t.\t.\tGT\t./.\n\
        c\t24\t.\tA\t<NON_REF>\t.\t.\tEND=31\tGT\t0/0\n\
        c\t35\t.\tAAA\tA\t.\t.\t.\tGT\t1/1\n"
    );
    std::fs::write(path("in.vcf"), text).unwrap();
    to_indexed_bcf(&path("in.vcf"), &path("in.bcf"));
    let bed = b"c\t9\t12\nc\t33\t36\n";
    std::fs::write(path("mask.bed"), bed).unwrap();
    for (start, end) in [(0, 48), (2, 40), (25, 30)] {
        for (overlap, overlap_arg) in [
            (Overlap::Start, "0"),
            (Overlap::Record, "1"),
            (Overlap::Variant, "2"),
        ] {
            for (with, mask_arg) in [(MaskWith::Lower, "lc"), (MaskWith::Char(b'N'), "N")] {
                let options = ConsensusOptions {
                    absent: Some(b'.'),
                    missing: Some(b'?'),
                    chain: true,
                    mark_ins: Some(Mark::Lower),
                    mark_del: Some(b'-'),
                    regions_overlap: overlap,
                    mask: Some(Arc::new(Mask::parse(bed, with).unwrap())),
                    ..ConsensusOptions::default()
                };
                std::fs::write(
                    path("ref.fa"),
                    format!(
                        ">c:{}-{end}\n{}\n",
                        start + 1,
                        "A".repeat((end - start) as usize)
                    ),
                )
                .unwrap();
                let args = [
                    "-f".to_owned(),
                    path("ref.fa").display().to_string(),
                    "-o".to_owned(),
                    path("out.fa").display().to_string(),
                    "-c".to_owned(),
                    path("out.chain").display().to_string(),
                    "-s".to_owned(),
                    "S".to_owned(),
                    "-H".to_owned(),
                    "1".to_owned(),
                    "-a".to_owned(),
                    ".".to_owned(),
                    "-M".to_owned(),
                    "?".to_owned(),
                    "--mark-ins".to_owned(),
                    "lc".to_owned(),
                    "--mark-del".to_owned(),
                    "-".to_owned(),
                    "--regions-overlap".to_owned(),
                    overlap_arg.to_owned(),
                    "-m".to_owned(),
                    path("mask.bed").display().to_string(),
                    "--mask-with".to_owned(),
                    mask_arg.to_owned(),
                    path("in.bcf").display().to_string(),
                ];
                assert_eq!(
                    bcftools_consensus(&args),
                    Outcome::Done(0),
                    "invalid oracle: {start}..{end}, {overlap:?}, {with:?}"
                );
                let expected = fasta_sequence(&std::fs::read_to_string(path("out.fa")).unwrap());
                let expected_chain = std::fs::read_to_string(path("out.chain")).unwrap();
                let mut consensus = open(&path("in.bcf"), None, options);
                consensus
                    .load(&mut genome(), Region::new("c", start, end))
                    .unwrap();
                let n = consensus
                    .query(Alleles::Haplotype {
                        sample: "S",
                        haplotype: "1",
                    })
                    .unwrap();
                assert_eq!(
                    output(&consensus, n),
                    expected,
                    "{start}..{end}, {overlap:?}, {with:?}"
                );
                assert_eq!(consensus.chain().unwrap(), expected_chain);
            }
        }
    }
}

#[test]
fn absent_skipped_records_and_claims_match_bcftools() {
    let dir = tempfile::tempdir().unwrap();
    let path = |name: &str| dir.path().join(name);
    let text = "##fileformat=VCFv4.3\n##contig=<ID=c,length=64>\n\
        ##INFO=<ID=END,Number=1,Type=Integer,Description=\"end\">\n\
        ##FORMAT=<ID=GT,Number=1,Type=String,Description=\"gt\">\n\
        #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS\n\
        c\t2\t.\tA\tC\t.\t.\t.\tGT\t./.\n\
        c\t3\t.\tA\tC\t.\t.\t.\tGT\t./.\n\
        c\t4\t.\tAAA\tA\t.\t.\t.\tGT\t1/1\n\
        c\t7\t.\tA\tC\t.\t.\t.\tGT\t./.\n\
        c\t8\t.\tA\tC\t.\t.\t.\tGT\t0/0\n\
        c\t9\t.\tA\tC\t.\t.\t.\tGT\t./.\n\
        c\t10\t.\tA\tATT\t.\t.\t.\tGT\t1/1\n\
        c\t12\t.\tA\tC\t.\t.\t.\tGT\t./.\n\
        c\t14\t.\tA\t<NON_REF>\t.\t.\tEND=20\tGT\t0/0\n\
        c\t21\t.\tA\tG\t.\t.\t.\tGT\t./.\n\
        c\t22\t.\tA\tC\t.\t.\t.\tGT\t./.\n\
        c\t24\t.\tA\tC\t.\t.\t.\tGT\t1/1\n";
    std::fs::write(path("input.vcf"), text).unwrap();
    to_indexed_bcf(&path("input.vcf"), &path("input.bcf"));
    for missing in [None, Some(b'?')] {
        for masked in [false, true] {
            let mut options = ConsensusOptions {
                absent: Some(b'.'),
                missing,
                chain: true,
                ..ConsensusOptions::default()
            };
            let mut args = vec![
                "-f".to_owned(),
                path("ref.fa").display().to_string(),
                "-o".to_owned(),
                path("out.fa").display().to_string(),
                "-c".to_owned(),
                path("out.chain").display().to_string(),
                "-s".to_owned(),
                "S".to_owned(),
                "-H".to_owned(),
                "1".to_owned(),
                "-a".to_owned(),
                ".".to_owned(),
            ];
            if missing.is_some() {
                args.extend(["-M".to_owned(), "?".to_owned()]);
            }
            if masked {
                let bed = b"c\t2\t3\nc\t23\t24\n";
                std::fs::write(path("mask.bed"), bed).unwrap();
                options.mask = Some(Arc::new(Mask::parse(bed, MaskWith::Char(b'N')).unwrap()));
                args.extend(["-m".to_owned(), path("mask.bed").display().to_string()]);
            }
            args.push(path("input.bcf").display().to_string());
            let bases = vec![b'A'; 64];
            std::fs::write(path("ref.fa"), format!(">c\n{}\n", "A".repeat(64))).unwrap();
            assert_eq!(bcftools_consensus(&args), Outcome::Done(0));
            let expected = fasta_sequence(&std::fs::read_to_string(path("out.fa")).unwrap());
            let mut genome = Genome(HashMap::from([("c".to_owned(), bases)]));
            let mut consensus = open(&path("input.bcf"), None, options);
            consensus
                .load(&mut genome, Region::new("c", 0, 64))
                .unwrap();
            for _ in 0..3 {
                let size = consensus
                    .query(Alleles::Haplotype {
                        sample: "S",
                        haplotype: "1",
                    })
                    .unwrap();
                let mut out = vec![0; size];
                consensus.write(Strand::Unstranded, &mut out).unwrap();
                assert_eq!(out, expected, "missing={missing:?} masked={masked}");
                assert_eq!(
                    consensus.chain().unwrap(),
                    std::fs::read_to_string(path("out.chain")).unwrap()
                );
            }
        }
    }
}

#[test]
fn repeated_masked_loads_and_gapped_chains_match_bcftools() {
    let dir = tempfile::tempdir().unwrap();
    let path = |name: &str| dir.path().join(name);
    let vcf = "##fileformat=VCFv4.3\n##contig=<ID=c,length=128>\n\
        ##FORMAT=<ID=GT,Number=1,Type=String,Description=\"gt\">\n\
        #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS\n\
        c\t5\t.\tA\tACG\t.\t.\t.\tGT\t0/1\n\
        c\t13\t.\tAAAAAAA\tA\t.\t.\t.\tGT\t0/1\n\
        c\t13\t.\tA\tC\t.\t.\t.\tGT\t1/0\n\
        c\t33\t.\tAA\tA\t.\t.\t.\tGT\t0/1\n\
        c\t49\t.\tA\tATTT\t.\t.\t.\tGT\t0/1\n\
        c\t73\t.\tA\tC\t.\t.\t.\tGT\t0/1\n";
    let bed = "c\t16\t18\nc\t72\t74\nc\t100\t110\nother\t0\t128\n";
    std::fs::write(path("in.vcf"), vcf).unwrap();
    std::fs::write(path("mask.bed"), bed).unwrap();
    to_indexed_bcf(&path("in.vcf"), &path("in.bcf"));
    let mut genome = Genome(HashMap::from([("c".to_owned(), vec![b'A'; 128])]));
    for (with, text) in [(MaskWith::Char(b'N'), "N"), (MaskWith::Lower, "lc")] {
        for absent in [false, true] {
            let mut consensus = open(
                &path("in.bcf"),
                None,
                ConsensusOptions {
                    chain: true,
                    absent: absent.then_some(b'.'),
                    mask: Some(Arc::new(Mask::parse(bed.as_bytes(), with).unwrap())),
                    ..ConsensusOptions::default()
                },
            );
            for (start, end) in [(0, 128), (24, 88), (112, 128), (0, 128)] {
                std::fs::write(
                    path("ref.fa"),
                    format!(
                        ">c:{}-{end}\n{}\n",
                        start + 1,
                        "A".repeat((end - start) as usize)
                    ),
                )
                .unwrap();
                consensus
                    .load(&mut genome, Region::new("c", start, end))
                    .unwrap();
                for haplotype in ["2", "1", "2"] {
                    let mut args: Vec<String> = ["-s", "S", "-H", haplotype]
                        .into_iter()
                        .map(str::to_owned)
                        .collect();
                    args.extend([
                        "-m".to_owned(),
                        path("mask.bed").display().to_string(),
                        "--mask-with".to_owned(),
                        text.to_owned(),
                    ]);
                    if absent {
                        args.extend(["-a".to_owned(), ".".to_owned()]);
                    }
                    for (flag, file) in [("-f", "ref.fa"), ("-o", "out.fa"), ("-c", "out.chain")] {
                        args.extend([flag.to_owned(), path(file).display().to_string()]);
                    }
                    args.push(path("in.bcf").display().to_string());
                    assert_eq!(bcftools_consensus(&args), Outcome::Done(0));
                    let len = consensus
                        .query(Alleles::Haplotype {
                            sample: "S",
                            haplotype,
                        })
                        .unwrap();
                    let mut out = vec![0; len];
                    consensus.write(Strand::Unstranded, &mut out).unwrap();
                    assert_eq!(
                        out,
                        fasta_sequence(&std::fs::read_to_string(path("out.fa")).unwrap())
                    );
                    assert_eq!(
                        consensus.chain().unwrap(),
                        std::fs::read_to_string(path("out.chain")).unwrap()
                    );
                }
            }
        }
    }
}
