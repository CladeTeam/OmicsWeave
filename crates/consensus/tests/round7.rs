mod common;

use std::collections::HashMap;
use std::sync::Arc;

use common::{
    Genome, Outcome, bcftools_consensus, consensus as open, fasta_sequence, to_indexed_bcf,
};
use omicsweave_consensus::{Alleles, ConsensusOptions, Mask, MaskWith};
use omicsweave_core::{Region, Strand};

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
