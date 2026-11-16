use std::collections::HashMap;
use std::fmt::Write as _;

use crate::support::{
    Genome, Outcome, bcf_of, bcftools_consensus, consensus as open, fasta_sequence, to_indexed_bcf,
};
use omicsweave_consensus::{Alleles, Consensus, ConsensusOptions, Overlap, RecordFilter};
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

#[test]
fn insertion_repopulates_a_symbolically_deleted_window_with_case_filters_and_chain() {
    let dir = tempfile::tempdir().unwrap();
    let path = |name: &str| dir.path().join(name);
    for (length, deletion_start) in [(20, 1), (512, 481)] {
        let insertion_start = deletion_start + 8;
        let region_start = deletion_start + 5;
        let region_end = insertion_start + 1;
        let text = format!(
            "##fileformat=VCFv4.3\n##contig=<ID=c,length={length}>\n\
             ##INFO=<ID=END,Number=1,Type=Integer,Description=\"end\">\n\
             ##FORMAT=<ID=GT,Number=1,Type=String,Description=\"gt\">\n\
             #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS\n\
             c\t{}\t.\tGC\tG,<DEL>\t.\t.\tEND={}\tGT\t0/2\n\
             c\t{}\t.\tG\tT,<DEL>,GCT\t.\t.\t.\tGT\t0/3\n",
            deletion_start + 1,
            deletion_start + 15,
            insertion_start + 1,
        );
        std::fs::write(path("in.vcf"), text).unwrap();
        to_indexed_bcf(&path("in.vcf"), &path("in.bcf"));
        // TYPE inspects all ALT alleles. Mixed types keep both records for
        // each expression, independently of which haplotype is selected.
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
            for (input, spelling, lower_reference) in [
                (b"AAAG".as_slice(), b"CT".as_slice(), false),
                (b"aaaG", b"ct", false),
                (b"aaag", b"ct", true),
            ] {
                std::fs::write(
                    path("ref.fa"),
                    format!(
                        ">c:{}-{region_end}\n{}\n",
                        region_start + 1,
                        String::from_utf8_lossy(input)
                    ),
                )
                .unwrap();
                let mut args: Vec<String> = ["--regions-overlap", "2", "-s", "S", "-H", "1", "-I"]
                    .into_iter()
                    .map(str::to_owned)
                    .collect();
                for (flag, file) in [("-f", "ref.fa"), ("-o", "out.fa"), ("-c", "out.chain")] {
                    args.extend([flag.to_owned(), path(file).display().to_string()]);
                }
                if let Some((flag, expression)) = filter_arg {
                    args.extend([flag.to_owned(), expression.to_owned()]);
                }
                args.push(path("in.bcf").display().to_string());
                assert_eq!(bcftools_consensus(&args), Outcome::Done(0));
                let expected = fasta_sequence(&std::fs::read_to_string(path("out.fa")).unwrap());
                assert_eq!(
                    expected, spelling,
                    "independent spelling: start={deletion_start}, input={input:?}"
                );
                let expected_chain = std::fs::read_to_string(path("out.chain")).unwrap();
                let mut bases = vec![b'A'; length];
                bases[deletion_start..deletion_start + 2].copy_from_slice(b"GC");
                if lower_reference {
                    bases.make_ascii_lowercase();
                }
                bases[region_start..region_end].copy_from_slice(input);
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
                    .load(
                        &mut reference,
                        Region::new("c", region_start as u32, region_end as u32),
                    )
                    .unwrap();
                for _ in 0..3 {
                    let n = consensus
                        .query(Alleles::Haplotype {
                            sample: "S",
                            haplotype: "1",
                        })
                        .unwrap();
                    assert_eq!(
                        output(&consensus, n),
                        expected,
                        "filter={filter:?}, start={deletion_start}, input={input:?}"
                    );
                    assert_eq!(consensus.chain().unwrap(), expected_chain);
                }
            }
        }
    }
}
