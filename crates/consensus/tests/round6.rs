mod common;

use std::collections::HashMap;
use std::fmt::Write as _;
use std::sync::Arc;

use common::{Genome, Outcome, bcftools_consensus, fasta_sequence, to_indexed_bcf};
use omicsweave_consensus::{Alleles, Consensus, ConsensusOptions, Mark, VcfStore};
use omicsweave_core::Region;

#[test]
fn iupac_reused_output_matches_bcftools_across_lengths_and_symbolic_fallbacks() {
    let dir = tempfile::tempdir().unwrap();
    let path = |name: &str| dir.path().join(name);
    for width in [1, 65, 4096] {
        for reference in ["A".repeat(width), "X".repeat(width)] {
            for alts in [
                format!("{},{}", "C".repeat(width), "G".repeat(width)),
                "C,G".to_owned(),
                format!("<DEL>,{}", "G".repeat(width)),
                "<DEL>,<NON_REF>".to_owned(),
            ] {
                let mut vcf = String::from(
                    "##fileformat=VCFv4.3\n##contig=<ID=c>\n\
                 ##FORMAT=<ID=GT,Number=1,Type=String,Description=\"gt\">\n\
                 #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS\tT\n",
                );
                writeln!(vcf, "c\t1\t.\t{reference}\t{alts}\t.\t.\t.\tGT\t0/1\t2/2").unwrap();
                std::fs::write(path("in.vcf"), &vcf).unwrap();
                to_indexed_bcf(&path("in.vcf"), &path("in.bcf"));
                let store = Arc::new(
                    VcfStore::read(std::fs::File::open(path("in.bcf")).unwrap(), None).unwrap(),
                );
                for lower in [false, true] {
                    let bases = if lower {
                        reference.to_ascii_lowercase()
                    } else {
                        reference.clone()
                    };
                    std::fs::write(path("ref.fa"), format!(">c\n{bases}\n")).unwrap();
                    let mut genome = Genome(HashMap::from([("c".to_owned(), bases.into_bytes())]));
                    let mut consensus = Consensus::new(
                        store.clone(),
                        ConsensusOptions {
                            iupac_codes: true,
                            mark_snv: Some(Mark::Lower),
                            chain: true,
                            ..ConsensusOptions::default()
                        },
                    );
                    consensus
                        .load(&mut genome, Region::new("c", 0, width as u32))
                        .unwrap();
                    // Alternate growing, shrinking and symbolic choices in the
                    // same builder, with each oracle invocation starting fresh.
                    for mode in [0, 1, 2, 0, 2, 1] {
                        let mut args = vec![
                            "-I".to_owned(),
                            "--mark-snv".to_owned(),
                            "lc".to_owned(),
                            "-s".to_owned(),
                            ["-", "S,T", "T"][mode].to_owned(),
                        ];
                        if mode == 2 {
                            args.extend(["-H".to_owned(), "I".to_owned()]);
                        }
                        for (flag, name) in
                            [("-f", "ref.fa"), ("-o", "out.fa"), ("-c", "out.chain")]
                        {
                            args.extend([flag.to_owned(), path(name).display().to_string()]);
                        }
                        args.push(path("in.bcf").display().to_string());
                        assert_eq!(bcftools_consensus(&args), Outcome::Done(0), "{vcf}");
                        let chosen = match mode {
                            0 => Alleles::AllAlt,
                            1 => Alleles::Samples(&["S", "T"]),
                            _ => Alleles::Haplotype {
                                sample: "T",
                                haplotype: "I",
                            },
                        };
                        let n = consensus.query(chosen).unwrap();
                        let mut out = vec![0; n];
                        consensus.write(&mut out).unwrap();
                        assert_eq!(
                            out,
                            fasta_sequence(&std::fs::read_to_string(path("out.fa")).unwrap()),
                            "width={width}, lower={lower}, mode={mode}, alts={alts}"
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
}
