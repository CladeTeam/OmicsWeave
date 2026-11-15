mod common;

use std::collections::HashMap;

use common::{
    Genome, Outcome, bcftools_consensus, consensus as open, fasta_sequence, to_indexed_bcf,
};
use omicsweave_consensus::{Alleles, ConsensusOptions};
use omicsweave_core::{Region, Strand};

#[test]
fn case_insensitive_selectors_and_window_chains_match_bcftools() {
    let dir = tempfile::tempdir().unwrap();
    let path = |name: &str| dir.path().join(name);
    let vcf = "##fileformat=VCFv4.3\n##contig=<ID=c,length=64>\n\
        ##FORMAT=<ID=GT,Number=1,Type=String,Description=\"gt\">\n\
        #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS\n\
        c\t5\t.\tA\tC\t.\t.\t.\tGT\t0|1\n\
        c\t13\t.\tA\tATTT\t.\t.\t.\tGT\t1|0\n\
        c\t25\t.\tAAA\tA\t.\t.\t.\tGT\t0/1\n";
    std::fs::write(path("in.vcf"), vcf).unwrap();
    to_indexed_bcf(&path("in.vcf"), &path("in.bcf"));
    let mut genome = Genome(HashMap::from([("c".to_owned(), vec![b'A'; 64])]));
    let mut consensus = open(
        &path("in.bcf"),
        None,
        ConsensusOptions {
            chain: true,
            ..ConsensusOptions::default()
        },
    );
    for (start, end) in [(0, 64), (8, 48), (0, 64)] {
        std::fs::write(
            path("window.fa"),
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
        for haplotype in [
            None,
            Some("1"),
            Some("2"),
            Some("r"),
            Some("a"),
            Some("l"),
            Some("s"),
            Some("lR"),
            Some("lA"),
            Some("sR"),
            Some("sA"),
            Some("i"),
            Some("1PiU"),
            Some("2pIu"),
        ] {
            let mut args: Vec<String> = match haplotype {
                Some(h) => ["-s", "S", "-H", h]
                    .into_iter()
                    .map(str::to_owned)
                    .collect(),
                None => ["-s", "-"].into_iter().map(str::to_owned).collect(),
            };
            for (flag, file) in [("-f", "window.fa"), ("-o", "out.fa"), ("-c", "out.chain")] {
                args.extend([flag.to_owned(), path(file).display().to_string()]);
            }
            args.push(path("in.bcf").display().to_string());
            assert_eq!(bcftools_consensus(&args), Outcome::Done(0));
            let mode = match haplotype {
                Some(haplotype) => Alleles::Haplotype {
                    sample: "S",
                    haplotype,
                },
                None => Alleles::AllAlt,
            };
            let len = consensus.query(mode).unwrap();
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
