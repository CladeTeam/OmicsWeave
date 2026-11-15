mod common;

use std::collections::HashMap;

use common::{
    Genome, Outcome, bcftools_consensus, consensus as open, fasta_sequence, to_indexed_bcf,
};
use omicsweave_consensus::{Alleles, ConsensusOptions};
use omicsweave_core::{Region, Strand};

fn byte_missing_genotype(count: u32) -> Vec<u8> {
    let header = b"##fileformat=VCFv4.3\n\
        ##FILTER=<ID=PASS,Description=\"pass\",IDX=0>\n\
        ##contig=<ID=c,length=16,IDX=0>\n\
        ##FORMAT=<ID=GT,Number=1,Type=String,Description=\"gt\",IDX=1>\n\
        #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS\n\0";
    let mut bytes = b"BCF\x02\x02".to_vec();
    bytes.extend_from_slice(&(header.len() as u32).to_le_bytes());
    bytes.extend_from_slice(header);
    let mut shared = vec![0x07];
    for i in 0..count {
        shared.extend_from_slice(&[
            0x17,
            if i == 0 {
                b'A'
            } else if i == 63 {
                b'G'
            } else {
                b'C'
            },
        ]);
    }
    shared.push(0);
    // Integer key 1, a diploid int8 vector: ordinary missing followed by the
    // BCF integer missing sentinel. Its low phase bit remains clear.
    let indiv = [0x11, 1, 0x21, 0, 0x80];
    bytes.extend_from_slice(&(24 + shared.len() as u32).to_le_bytes());
    bytes.extend_from_slice(&(indiv.len() as u32).to_le_bytes());
    for value in [0u32, 2, 1, 0x7f80_0001, count << 16, (1 << 24) | 1] {
        bytes.extend_from_slice(&value.to_le_bytes());
    }
    bytes.extend_from_slice(&shared);
    bytes.extend_from_slice(&indiv);
    bytes
}

#[test]
fn byte_missing_sentinel_keeps_its_unsigned_allele_63_behavior() {
    let dir = tempfile::tempdir().unwrap();
    let path = |name: &str| dir.path().join(name);
    std::fs::write(path("ref.fa"), ">c\nAAAAAAAAAAAAAAAA\n").unwrap();
    for count in [63, 64, 65] {
        let bytes = byte_missing_genotype(count);
        std::fs::write(path("raw.bcf"), &bytes).unwrap();
        to_indexed_bcf(&path("raw.bcf"), &path("in.bcf"));
        let mut consensus = open(&path("in.bcf"), None, ConsensusOptions::default());
        consensus
            .load(
                &mut Genome(HashMap::from([("c".to_owned(), vec![b'A'; 16])])),
                Region::new("c", 0, 16),
            )
            .unwrap();
        let args = [
            "-s".to_owned(),
            "S".to_owned(),
            "-f".to_owned(),
            path("ref.fa").display().to_string(),
            "-o".to_owned(),
            path("out.fa").display().to_string(),
            path("in.bcf").display().to_string(),
        ];
        let result = consensus.query(Alleles::Samples(&["S"]));
        if count == 63 {
            assert!(result.is_err());
            assert_eq!(bcftools_consensus(&args), Outcome::Failed);
        } else {
            let mut out = vec![0; result.unwrap()];
            consensus.write(Strand::Unstranded, &mut out).unwrap();
            let mut expected = vec![b'A'; 16];
            expected[2] = b'G';
            assert_eq!(out, expected);
            assert_eq!(bcftools_consensus(&args), Outcome::Done(0));
            assert_eq!(
                out,
                fasta_sequence(&std::fs::read_to_string(path("out.fa")).unwrap())
            );
        }
    }
}

#[test]
fn many_unselected_alleles_preserve_sample_iupac_and_full_allele_modes() {
    let dir = tempfile::tempdir().unwrap();
    let path = |name: &str| dir.path().join(name);
    std::fs::write(path("ref.fa"), ">c\nAAAAAAAAAAAAAAAA\n").unwrap();
    for count in [63, 64, 65, 4096] {
        let mut alts = vec!["C"; count - 1];
        alts[61] = "G";
        if count > 64 {
            alts[count - 2] = "AA";
        }
        for genotype in ["0/1", "0/62", "0/63"] {
            if count == 63 && genotype == "0/63" {
                continue;
            }
            let vcf = format!(
                "##fileformat=VCFv4.3\n##contig=<ID=c,length=16>\n\
                 ##FORMAT=<ID=GT,Number=1,Type=String,Description=\"gt\">\n\
                 #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS\n\
                 c\t3\t.\tA\t{}\t.\t.\t.\tGT\t{genotype}\n",
                alts.join(",")
            );
            std::fs::write(path("in.vcf"), &vcf).unwrap();
            to_indexed_bcf(&path("in.vcf"), &path("in.bcf"));
            let mut reference = Genome(HashMap::from([("c".to_owned(), vec![b'A'; 16])]));
            let mut consensus = open(&path("in.bcf"), None, ConsensusOptions::default());
            consensus
                .load(&mut reference, Region::new("c", 0, 16))
                .unwrap();
            let args = [
                "-s".to_owned(),
                "S".to_owned(),
                "-f".to_owned(),
                path("ref.fa").display().to_string(),
                "-o".to_owned(),
                path("out.fa").display().to_string(),
                path("in.bcf").display().to_string(),
            ];
            let result = consensus.query(Alleles::Samples(&["S"]));
            if genotype == "0/63" {
                assert!(result.is_err());
                assert_eq!(bcftools_consensus(&args), Outcome::Failed);
            } else {
                let mut out = vec![0; result.unwrap()];
                consensus.write(Strand::Unstranded, &mut out).unwrap();
                let mut expected = vec![b'A'; 16];
                expected[2] = if genotype == "0/1" { b'M' } else { b'R' };
                assert_eq!(out, expected);
                assert_eq!(bcftools_consensus(&args), Outcome::Done(0));
                assert_eq!(
                    out,
                    fasta_sequence(&std::fs::read_to_string(path("out.fa")).unwrap())
                );
            }

            if count > 64 {
                let mut all = open(
                    &path("in.bcf"),
                    None,
                    ConsensusOptions {
                        iupac_codes: true,
                        ..ConsensusOptions::default()
                    },
                );
                all.load(&mut reference, Region::new("c", 0, 16)).unwrap();
                let mut out = vec![0; all.query(Alleles::AllAlt).unwrap()];
                all.write(Strand::Unstranded, &mut out).unwrap();
                let mut expected = vec![b'A'; 16];
                expected[2] = b'C';
                assert_eq!(out, expected);
                let args = [
                    "-I".to_owned(),
                    "-s".to_owned(),
                    "-".to_owned(),
                    "-f".to_owned(),
                    path("ref.fa").display().to_string(),
                    "-o".to_owned(),
                    path("out.fa").display().to_string(),
                    path("in.bcf").display().to_string(),
                ];
                assert_eq!(bcftools_consensus(&args), Outcome::Done(0));
                assert_eq!(
                    out,
                    fasta_sequence(&std::fs::read_to_string(path("out.fa")).unwrap())
                );
            }
        }
    }
}
