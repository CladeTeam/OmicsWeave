use crate::support::{Genome, Outcome, bcftools_consensus, consensus as open, to_indexed_raw_bcf};
use omicsweave_consensus::{Alleles, ConsensusOptions};
use omicsweave_core::{Region, Strand};
use std::collections::HashMap;

fn genome() -> Genome {
    Genome(HashMap::from([("c".to_owned(), vec![b'A'; 128])]))
}

fn encoded_gt(integer_type: u8) -> Vec<u8> {
    let header = b"##fileformat=VCFv4.3\n\
        ##FILTER=<ID=PASS,Description=\"pass\",IDX=0>\n\
        ##contig=<ID=c,length=16,IDX=0>\n\
        ##FORMAT=<ID=GT,Number=1,Type=String,Description=\"gt\",IDX=1>\n\
        #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS\n\0";
    let mut bytes = b"BCF\x02\x02".to_vec();
    bytes.extend_from_slice(&(header.len() as u32).to_le_bytes());
    bytes.extend_from_slice(header);
    let shared = [0x07, 0x17, b'A', 0x17, b'C', 0x00];
    let mut indiv = vec![0x11, 1, 0x20 | integer_type];
    for value in [2_i32, 5] {
        match integer_type {
            1 => indiv.push(value as u8),
            2 => indiv.extend_from_slice(&(value as i16).to_le_bytes()),
            3 => indiv.extend_from_slice(&value.to_le_bytes()),
            _ => unreachable!(),
        }
    }
    bytes.extend_from_slice(&(24 + shared.len() as u32).to_le_bytes());
    bytes.extend_from_slice(&(indiv.len() as u32).to_le_bytes());
    for value in [0u32, 2, 1, 0x7f80_0001, 2 << 16, (1 << 24) | 1] {
        bytes.extend_from_slice(&value.to_le_bytes());
    }
    bytes.extend_from_slice(&shared);
    bytes.extend_from_slice(&indiv);
    bytes
}

#[test]
fn physical_genotype_width_matches_bcftools_even_for_small_values() {
    let dir = tempfile::tempdir().unwrap();
    let fasta = dir.path().join("ref.fa");
    std::fs::write(&fasta, ">c\nAAAAAAAAAAAAAAAA\n").unwrap();
    for (integer_type, iupac) in [1, 2, 3]
        .into_iter()
        .flat_map(|ty| [(ty, false), (ty, true)])
    {
        let bcf = dir.path().join(format!("type{integer_type}.bcf"));
        to_indexed_raw_bcf(&encoded_gt(integer_type), &bcf);
        let mut consensus = open(
            &bcf,
            None,
            ConsensusOptions {
                chain: true,
                iupac_codes: iupac,
                ..ConsensusOptions::default()
            },
        );
        let mut reference = Genome(HashMap::from([("c".to_owned(), vec![b'A'; 16])]));
        consensus
            .load(&mut reference, Region::new("c", 0, 16))
            .unwrap();
        let mut cases = vec![
            (Alleles::AllAlt, vec!["-s", "-"]),
            (Alleles::Samples(&["S"]), vec!["-s", "S"]),
        ];
        cases.extend(
            [
                "1", "2", "R", "A", "L", "S", "LR", "LA", "SR", "SA", "I", "1pIu", "2pIu",
            ]
            .map(|haplotype| {
                (
                    Alleles::Haplotype {
                        sample: "S",
                        haplotype,
                    },
                    vec!["-s", "S", "-H", haplotype],
                )
            }),
        );
        for (alleles, selector) in cases {
            let mut args: Vec<_> = selector.iter().map(|s| (*s).to_owned()).collect();
            if iupac {
                args.push("-I".to_owned());
            }
            let output = dir.path().join("out.fa");
            args.extend([
                "-f".to_owned(),
                fasta.display().to_string(),
                "-o".to_owned(),
                output.display().to_string(),
                bcf.display().to_string(),
            ]);
            let oracle = bcftools_consensus(&args);
            let actual = consensus.query(alleles);
            if integer_type == 1 || matches!(alleles, Alleles::AllAlt) {
                assert_eq!(
                    oracle,
                    Outcome::Done(0),
                    "type={integer_type}, iupac={iupac}, selector={selector:?}"
                );
                let n = actual.unwrap_or_else(|error| {
                    panic!("type={integer_type}, iupac={iupac}, selector={selector:?}: {error}")
                });
                let mut sequence = vec![0; n];
                consensus.write(Strand::Unstranded, &mut sequence).unwrap();
                assert_eq!(
                    sequence,
                    crate::support::fasta_sequence(&std::fs::read_to_string(output).unwrap()),
                    "type={integer_type}, iupac={iupac}, selector={selector:?}"
                );
            } else {
                assert_eq!(oracle, Outcome::Failed);
                assert!(actual.is_err(), "type={integer_type} selector={selector:?}");
                assert!(consensus.chain().is_none());
                let mut guard = [0xa5; 16];
                assert!(consensus.write(Strand::Unstranded, &mut guard).is_err());
                assert_eq!(guard, [0xa5; 16]);
            }
        }
        assert_eq!(consensus.query(Alleles::AllAlt).unwrap(), 16);
    }
}

/// A structurally decoded BCF with one sample and a zero-width INT8 GT.
/// No external oracle is used for malformed genotype layouts.
fn zero_width_gt() -> Vec<u8> {
    let header = b"##fileformat=VCFv4.3\n\
        ##FILTER=<ID=PASS,Description=\"pass\",IDX=0>\n\
        ##contig=<ID=c,length=128,IDX=0>\n\
        ##FORMAT=<ID=GT,Number=1,Type=String,Description=\"gt\",IDX=1>\n\
        #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS\n\0";
    let mut bytes = b"BCF\x02\x02".to_vec();
    bytes.extend_from_slice(&(header.len() as u32).to_le_bytes());
    bytes.extend_from_slice(header);
    let shared = [0x07, 0x17, b'A', 0x17, b'C', 0x00];
    let indiv = [0x11, 1, 0x01];
    bytes.extend_from_slice(&(24 + shared.len() as u32).to_le_bytes());
    bytes.extend_from_slice(&(indiv.len() as u32).to_le_bytes());
    for value in [0u32, 1, 1, 0x7f80_0001, 2 << 16, (1 << 24) | 1] {
        bytes.extend_from_slice(&value.to_le_bytes());
    }
    bytes.extend_from_slice(&shared);
    bytes.extend_from_slice(&indiv);
    bytes
}

#[test]
fn zero_width_genotypes_produce_errors_instead_of_index_panics() {
    let bytes = zero_width_gt();
    // Verify that the fixture reaches genotype decoding, independently of
    // the decision to reject that layout when a region is loaded.
    let mut reader = omicsweave_vcf::VariantReader::new(bytes.as_slice()).unwrap();
    let mut record = omicsweave_vcf::Record::default();
    assert!(reader.read(&mut record).unwrap());
    assert_eq!(record.genotypes().unwrap().width, 0);
    let dir = tempfile::tempdir().unwrap();
    let bcf = dir.path().join("in.bcf");
    to_indexed_raw_bcf(&bytes, &bcf);
    let mut consensus = open(&bcf, None, ConsensusOptions::default());
    let loaded = consensus
        .load(&mut genome(), Region::new("c", 0, 8))
        .is_ok();
    for haplotype in ["1", "2", "1pIu", "I", "R"] {
        assert!(
            consensus
                .query(Alleles::Haplotype {
                    sample: "S",
                    haplotype
                })
                .is_err(),
            "zero-width GT must be rejected for {haplotype}, loaded={loaded}"
        );
        let mut guard = [0xa5; 8];
        assert!(consensus.write(Strand::Forward, &mut guard).is_err());
        assert_eq!(guard, [0xa5; 8]);
        assert!(consensus.chain().is_none());
    }
}
