mod common;

use std::collections::HashMap;

use common::{Genome, bcf_of, consensus as open};
use omicsweave_consensus::{Alleles, ConsensusOptions};
use omicsweave_core::orient::reverse_complement;
use omicsweave_core::{Region, Strand};

/// A soft-masked reference with SNPs, an insertion, a deletion and a
/// multi-allelic site, so forward outputs hold IUPAC codes in both cases.
const VCF: &[u8] = b"##fileformat=VCFv4.3\n##contig=<ID=c,length=32>\n\
    ##FORMAT=<ID=GT,Number=1,Type=String,Description=\"gt\">\n\
    #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS\tT\n\
    c\t3\t.\tA\tG\t.\t.\t.\tGT\t0|1\t1|1\n\
    c\t6\t.\tC\tCTTG\t.\t.\t.\tGT\t1|0\t0|1\n\
    c\t11\t.\tGAT\tG\t.\t.\t.\tGT\t1|1\t0|0\n\
    c\t20\t.\tT\tA,C\t.\t.\t.\tGT\t1|2\t0|2\n\
    c\t25\t.\tG\tT\t.\t.\t.\tGT\t0/1\t1/1\n";

fn reference() -> Genome {
    let seq = b"ACACGCTAGGGATCCATGCTtacGGTACCGTA".to_vec();
    Genome(HashMap::from([("c".to_owned(), seq)]))
}

#[test]
fn reverse_strand_is_the_reverse_complement_of_the_forward_sequence() {
    let dir = tempfile::tempdir().unwrap();
    let bcf = bcf_of(dir.path(), "in", VCF);
    for iupac_codes in [false, true] {
        let mut consensus = open(
            &bcf,
            None,
            ConsensusOptions {
                iupac_codes,
                chain: true,
                ..ConsensusOptions::default()
            },
        );
        for (start, end) in [(0, 32), (4, 22), (9, 10)] {
            consensus
                .load(&mut reference(), Region::new("c", start, end))
                .unwrap();
            for alleles in [
                Alleles::AllAlt,
                Alleles::Samples(&["S", "T"]),
                Alleles::Haplotype {
                    sample: "S",
                    haplotype: "1",
                },
                Alleles::Haplotype {
                    sample: "T",
                    haplotype: "2",
                },
                Alleles::Haplotype {
                    sample: "S",
                    haplotype: "I",
                },
            ] {
                let n = consensus.query(alleles).unwrap();
                let mut forward = vec![0; n];
                consensus.write(Strand::Forward, &mut forward).unwrap();
                for unstranded in [Strand::Unstranded, Strand::Unknown] {
                    let mut out = vec![0; n];
                    consensus.write(unstranded, &mut out).unwrap();
                    assert_eq!(out, forward);
                }
                let chain = consensus.chain().unwrap().to_owned();
                let mut reverse = vec![0; n];
                consensus.write(Strand::Reverse, &mut reverse).unwrap();
                let mut expected = forward.clone();
                reverse_complement(&mut expected);
                assert_eq!(
                    reverse,
                    expected,
                    "{alleles:?} {start}..{end} forward {}",
                    String::from_utf8_lossy(&forward)
                );
                // The chain maps forward coordinates whatever the strand.
                assert_eq!(consensus.chain().unwrap(), chain);
            }
        }
    }
}

#[test]
fn reverse_strand_complements_iupac_codes_and_keeps_case() {
    let dir = tempfile::tempdir().unwrap();
    let bcf = bcf_of(dir.path(), "in", VCF);
    let mut consensus = open(
        &bcf,
        None,
        ConsensusOptions {
            iupac_codes: true,
            ..ConsensusOptions::default()
        },
    );
    consensus
        .load(&mut reference(), Region::new("c", 18, 26))
        .unwrap();
    let n = consensus.query(Alleles::Samples(&["S"])).unwrap();
    let mut forward = vec![0; n];
    consensus.write(Strand::Forward, &mut forward).unwrap();
    // c:19-26 reads `CTtacGGT`: the multi-allelic site at 20 becomes M
    // (A/C) and the unphased SNP at 25 K (G/T).
    assert_eq!(forward, b"CMtacGKT");
    let mut reverse = vec![0; n];
    consensus.write(Strand::Reverse, &mut reverse).unwrap();
    assert_eq!(reverse, b"AMCgtaKG");
}
