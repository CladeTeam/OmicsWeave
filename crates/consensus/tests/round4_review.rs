mod common;

use std::collections::HashMap;
use std::fmt::Write as _;
use std::sync::Arc;

use common::{
    Genome, Outcome, bcf_of, bcftools_consensus, consensus as open, fasta_sequence, to_indexed_bcf,
    to_indexed_raw_bcf,
};
use omicsweave_consensus::{
    Alleles, Consensus, ConsensusOptions, Mark, Mask, MaskWith, Overlap, RecordFilter,
};
use omicsweave_core::{Error, Locus, Region, Result, SequenceRead, Strand};

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
    if consensus
        .load(&mut genome(), Region::new("c", 0, 8))
        .is_err()
    {
        return;
    }
    for haplotype in ["1", "2", "1pIu", "I", "R"] {
        assert!(
            consensus
                .query(Alleles::Haplotype {
                    sample: "S",
                    haplotype
                })
                .is_err(),
            "zero-width GT must be rejected for {haplotype}"
        );
    }
}

#[test]
fn disjoint_edits_and_absent_claims_match_an_independent_spelling_oracle() {
    // The oracle appends whole disjoint spans; it has no index, modified
    // coordinate, freeze position, or in-place edit state.
    for phase in 0..7 {
        let mut text = HEADER.to_owned();
        let mut expected = Vec::new();
        let mut cursor = 0;
        for i in 0..16usize {
            let start = i * 7 + 2;
            let (reference, alt) = match i % 4 {
                0 => ("A", "C"),
                1 => ("A", "ATT"),
                2 => ("AAA", "A"),
                _ => ("AA", "CG"),
            };
            let gt = match (i + phase) % 3 {
                0 => "0/0",
                1 => "1/1",
                _ => "./.",
            };
            writeln!(
                text,
                "c\t{}\t.\t{reference}\t{alt}\t.\t.\t.\tGT\t{gt}",
                start + 1
            )
            .unwrap();
            expected.resize(expected.len() + start - cursor, b'.');
            match gt {
                "0/0" => expected.extend_from_slice(reference.as_bytes()),
                "1/1" => expected.extend_from_slice(alt.as_bytes()),
                _ => expected.resize(expected.len() + reference.len(), b'.'),
            }
            cursor = start + reference.len();
        }
        expected.resize(expected.len() + 128 - cursor, b'.');
        let dir = tempfile::tempdir().unwrap();
        let mut consensus = open(
            &bcf_of(dir.path(), "in", text.as_bytes()),
            None,
            ConsensusOptions {
                absent: Some(b'.'),
                ..ConsensusOptions::default()
            },
        );
        consensus
            .load(&mut genome(), Region::new("c", 0, 128))
            .unwrap();
        for _ in 0..3 {
            let n = consensus
                .query(Alleles::Haplotype {
                    sample: "S",
                    haplotype: "1",
                })
                .unwrap();
            assert_eq!(output(&consensus, n), expected, "phase={phase}");
        }
    }
}

struct FailingReference;

impl SequenceRead for FailingReference {
    fn chrom_len(&self, _: &str) -> Result<u32> {
        Ok(128)
    }

    fn sequence(&mut self, _: &Locus<'_>, _: &mut [u8]) -> Result<()> {
        Err(Error::invalid("reference", "injected read failure"))
    }
}

#[test]
fn failed_load_and_query_revoke_output_and_allow_recovery() {
    let dir = tempfile::tempdir().unwrap();
    let text = format!("{HEADER}c\t2\t.\tC\tG\t.\t.\t.\tGT\t0/1\n");
    let mut consensus = open(
        &bcf_of(dir.path(), "in", text.as_bytes()),
        None,
        ConsensusOptions {
            chain: true,
            ..ConsensusOptions::default()
        },
    );
    consensus
        .load(&mut genome(), Region::new("c", 0, 8))
        .unwrap();
    let reference_haplotype = || Alleles::Haplotype {
        sample: "S",
        haplotype: "1",
    };
    assert_eq!(consensus.query(reference_haplotype()).unwrap(), 8);
    assert!(consensus.chain().is_some());
    // The ALT path must inspect REF and reject this mismatching record.
    assert!(consensus.query(Alleles::AllAlt).is_err());
    let mut guard = [0xa5; 8];
    assert!(matches!(
        consensus.write(Strand::Unstranded, &mut guard),
        Err(Error::NotQueried)
    ));
    assert_eq!(guard, [0xa5; 8]);
    assert!(consensus.chain().is_none());
    assert_eq!(consensus.query(reference_haplotype()).unwrap(), 8);
    assert_eq!(output(&consensus, 8), b"AAAAAAAA");
    assert!(
        consensus
            .load(&mut FailingReference, Region::new("c", 0, 8))
            .is_err()
    );
    assert!(matches!(
        consensus.write(Strand::Unstranded, &mut guard),
        Err(Error::NotQueried)
    ));
    assert!(consensus.chain().is_none());
    consensus
        .load(&mut genome(), Region::new("c", 0, 3))
        .unwrap();
    assert_eq!(consensus.query(reference_haplotype()).unwrap(), 3);
    assert_eq!(output(&consensus, 3), b"AAA");
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
fn insertion_can_repopulate_a_region_emptied_by_a_symbolic_deletion() {
    let dir = tempfile::tempdir().unwrap();
    let path = |name: &str| dir.path().join(name);
    let text = "##fileformat=VCFv4.2\n##contig=<ID=c,length=512>\n\
        ##INFO=<ID=END,Number=1,Type=Integer,Description=\"end\">\n\
        ##FORMAT=<ID=GT,Number=1,Type=String,Description=\"gt\">\n\
        #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS\n\
        c\t482\t.\tGC\tG,<DEL>\t.\t.\tEND=496\tGT\t0/2\n\
        c\t490\t.\tG\tT,<DEL>,GCT\t.\t.\t.\tGT\t0/3\n";
    std::fs::write(path("in.vcf"), text).unwrap();
    to_indexed_bcf(&path("in.vcf"), &path("in.bcf"));
    // TYPE uses all ALT alleles. Each record has mixed variant types, so
    // equality to the single indel category does not match either: every
    // filter keeps both records.
    for (filter, filter_arg) in [
        (None, None),
        (
            Some(RecordFilter::Exclude("TYPE=\"indel\"")),
            Some(("-e", "TYPE=\"indel\"")),
        ),
        (
            Some(RecordFilter::Include("TYPE!=\"indel\"")),
            Some(("-i", "TYPE!=\"indel\"")),
        ),
    ] {
        for (input, spelling) in [(b"AAAG".as_slice(), b"CT".as_slice()), (b"aaaG", b"ct")] {
            std::fs::write(
                path("ref.fa"),
                format!(">c:487-490\n{}\n", String::from_utf8_lossy(input)),
            )
            .unwrap();
            let mut args = vec![
                "-f".to_owned(),
                path("ref.fa").display().to_string(),
                "-o".to_owned(),
                path("out.fa").display().to_string(),
                "-c".to_owned(),
                path("out.chain").display().to_string(),
                "--regions-overlap".to_owned(),
                "2".to_owned(),
                "-s".to_owned(),
                "S".to_owned(),
                "-H".to_owned(),
                "1".to_owned(),
                "-I".to_owned(),
            ];
            if let Some((flag, expression)) = filter_arg {
                args.extend([flag.to_owned(), expression.to_owned()]);
            }
            args.push(path("in.bcf").display().to_string());
            assert_eq!(bcftools_consensus(&args), Outcome::Done(0));
            let expected = fasta_sequence(&std::fs::read_to_string(path("out.fa")).unwrap());
            assert_eq!(
                expected, spelling,
                "the external oracle must repopulate the region"
            );
            let expected_chain = std::fs::read_to_string(path("out.chain")).unwrap();
            let mut bases = vec![b'A'; 512];
            bases[481..483].copy_from_slice(b"GC");
            bases[486..490].copy_from_slice(input);
            let mut reference = Genome(HashMap::from([("c".to_owned(), bases)]));
            let mut consensus = open(
                &path("in.bcf"),
                filter,
                ConsensusOptions {
                    iupac_codes: true,
                    regions_overlap: Overlap::Variant,
                    chain: true,
                    ..ConsensusOptions::default()
                },
            );
            consensus
                .load(&mut reference, Region::new("c", 486, 490))
                .unwrap();
            for _ in 0..3 {
                let n = consensus
                    .query(Alleles::Haplotype {
                        sample: "S",
                        haplotype: "1",
                    })
                    .unwrap();
                assert_eq!(output(&consensus, n), expected, "filter={filter:?}");
                assert_eq!(consensus.chain().unwrap(), expected_chain);
            }
        }
    }
}
