//! Single-record selections: one record with one ALT, every other record of
//! the region left out. A record applied alone must give the sequence a
//! BCF holding only that record (and only that ALT) gives with `AllAlt`.

#[path = "support/mod.rs"]
mod common;

use std::collections::HashMap;
use std::fmt::Write;

use common::{Genome, Rng, bcf_of, consensus as open};
use omicsweave_consensus::{Alleles, Consensus, ConsensusOptions, Mask, MaskWith};
use omicsweave_core::{Region, Strand};
use proptest::prelude::*;

const LEN: usize = 120;
const HEADER: &str = "##fileformat=VCFv4.3\n##contig=<ID=c,length=120>\n\
    ##INFO=<ID=DP,Number=1,Type=Integer,Description=\"depth\">\n\
    #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n";

struct Rec {
    pos: usize,
    reference: String,
    alts: Vec<String>,
}

fn genome(rng: &mut Rng) -> Vec<u8> {
    (0..LEN).map(|_| *rng.pick(&b"ACGT"[..])).collect()
}

fn base(rng: &mut Rng) -> char {
    *rng.pick(&['A', 'C', 'G', 'T'][..])
}

/// Overlapping SNVs, insertions, deletions and multi-allelic records.
fn records(rng: &mut Rng, seq: &[u8]) -> Vec<Rec> {
    let mut out: Vec<Rec> = (0..rng.range(1, 12))
        .map(|_| {
            let pos = rng.range(1, LEN as u64 - 8) as usize;
            let rlen = if rng.chance(30) {
                rng.range(2, 6) as usize
            } else {
                1
            };
            let reference: String = seq[pos..pos + rlen].iter().map(|&b| b as char).collect();
            let mut alts = (0..rng.range(1, 4))
                .map(|_| {
                    let mut alt = String::new();
                    alt.push(seq[pos] as char);
                    if rlen == 1 {
                        if rng.chance(30) {
                            for _ in 0..rng.range(1, 4) {
                                alt.push(base(rng));
                            }
                        } else {
                            alt.clear();
                            alt.push(base(rng));
                        }
                    }
                    alt
                })
                .filter(|a| a != &reference)
                .collect::<Vec<_>>();
            let mut seen = std::collections::HashSet::new();
            alts.retain(|a| seen.insert(a.clone()));
            Rec {
                pos,
                reference,
                alts,
            }
        })
        .filter(|r| !r.alts.is_empty())
        .collect();
    out.sort_by_key(|r| r.pos);
    out
}

fn vcf(records: &[&Rec], alt: Option<usize>) -> Vec<u8> {
    let mut text = HEADER.to_owned();
    for r in records {
        let alts = match alt {
            Some(k) => r.alts[k - 1].clone(),
            None => r.alts.join(","),
        };
        writeln!(
            text,
            "c\t{}\t.\t{}\t{alts}\t.\t.\tDP=1",
            r.pos + 1,
            r.reference
        )
        .unwrap();
    }
    text.into_bytes()
}

fn query(consensus: &mut Consensus, alleles: Alleles<'_>, strand: Strand) -> Vec<u8> {
    let n = consensus.query(alleles).unwrap();
    let mut out = vec![0; n];
    consensus.write(strand, &mut out).unwrap();
    out
}

fn reference(seq: &[u8]) -> Genome {
    Genome(HashMap::from([("c".to_owned(), seq.to_vec())]))
}

proptest! {
    #![proptest_config(ProptestConfig::with_cases(48))]

    #[test]
    fn one_record_matches_a_file_of_that_record_alone(seed in any::<u64>()) {
        let mut rng = Rng(seed);
        let seq = genome(&mut rng);
        let recs = records(&mut rng, &seq);
        let dir = tempfile::tempdir().unwrap();
        let all = bcf_of(dir.path(), "all", &vcf(&recs.iter().collect::<Vec<_>>(), None));
        let lo = rng.range(0, 40) as u32;
        let region = Region::new("c", lo, lo + rng.range(30, 80) as u32);
        let strand = if rng.chance(50) { Strand::Reverse } else { Strand::Forward };
        let mut genome = reference(&seq);
        let mut cohort = open(&all, None, ConsensusOptions::default());
        cohort.load(&mut genome, region).unwrap();
        let listed: Vec<(i64, Vec<u8>, usize)> = cohort
            .records()
            .unwrap()
            .map(|r| (r.start(), r.allele(0).to_vec(), r.n_alleles()))
            .collect();
        let taken: Vec<&Rec> = recs
            .iter()
            .filter(|r| r.pos < region.end as usize && r.pos + r.reference.len() > lo as usize)
            .collect();
        prop_assert_eq!(listed.len(), taken.len());
        for (index, rec) in taken.iter().enumerate() {
            prop_assert_eq!(listed[index].0, rec.pos as i64);
            prop_assert_eq!(&listed[index].1, rec.reference.as_bytes());
            for alt in 1..=rec.alts.len() {
                let name = format!("r{index}a{alt}");
                let alone = bcf_of(dir.path(), &name, &vcf(&[rec], Some(alt)));
                let mut single = open(&alone, None, ConsensusOptions::default());
                single.load(&mut genome, region).unwrap();
                let expected = query(&mut single, Alleles::AllAlt, strand);
                let by_index = query(&mut cohort, Alleles::Record { index, alt }, strand);
                prop_assert_eq!(&by_index, &expected, "record {} alt {}", index, alt);
                let by_content = query(
                    &mut cohort,
                    Alleles::Variant {
                        start: rec.pos as i64,
                        reference: rec.reference.as_bytes(),
                        alt: rec.alts[alt - 1].as_bytes(),
                    },
                    strand,
                );
                // The first record in file order with this REF and ALT.
                let first = taken
                    .iter()
                    .position(|r| {
                        r.pos == rec.pos
                            && r.reference == rec.reference
                            && r.alts.contains(&rec.alts[alt - 1])
                    })
                    .unwrap();
                if first == index {
                    prop_assert_eq!(&by_content, &expected);
                }
            }
        }
    }
}

#[test]
fn selections_name_what_the_region_lacks() {
    let dir = tempfile::tempdir().unwrap();
    let seq = b"ACGTACGTACGTACGTACGT";
    let text = format!("{HEADER}c\t3\t.\tG\tA,T\t.\t.\tDP=1\nc\t9\t.\tA\tC\t.\t.\tDP=1\n");
    let bcf = bcf_of(dir.path(), "v", text.as_bytes());
    let mut genome = reference(seq);
    let mut c = open(&bcf, None, ConsensusOptions::default());
    assert!(c.query(Alleles::Record { index: 0, alt: 1 }).is_err());
    c.load(&mut genome, Region::new("c", 0, 20)).unwrap();
    for alleles in [
        Alleles::Record { index: 2, alt: 1 },
        Alleles::Record { index: 0, alt: 0 },
        Alleles::Record { index: 0, alt: 3 },
        Alleles::Variant {
            start: 2,
            reference: b"G",
            alt: b"C",
        },
        Alleles::Variant {
            start: 3,
            reference: b"G",
            alt: b"A",
        },
    ] {
        assert!(c.query(alleles).is_err(), "{alleles:?}");
    }
    let lower = query(
        &mut c,
        Alleles::Variant {
            start: 2,
            reference: b"g",
            alt: b"t",
        },
        Strand::Forward,
    );
    assert_eq!(lower, b"ACTTACGTACGTACGTACGT");
    let all = query(&mut c, Alleles::AllAlt, Strand::Forward);
    assert_eq!(all, b"ACATACGTCCGTACGTACGT");
}

#[test]
fn a_masked_record_is_left_out_and_other_records_errors_do_not_count() {
    let dir = tempfile::tempdir().unwrap();
    let seq = b"ACGTACGTACGTACGTACGT";
    let text = format!("{HEADER}c\t3\t.\tG\tA\t.\t.\tDP=1\nc\t9\t.\tA\tC\t.\t.\tDP=1\n");
    let bcf = bcf_of(dir.path(), "v", text.as_bytes());
    let mut genome = reference(seq);
    let options = ConsensusOptions {
        mask: Some(std::sync::Arc::new(
            Mask::parse(b"c\t2\t3\n", MaskWith::Char(b'N')).unwrap(),
        )),
        ..ConsensusOptions::default()
    };
    let mut c = open(&bcf, None, options);
    c.load(&mut genome, Region::new("c", 0, 20)).unwrap();
    let masked: Vec<bool> = c.records().unwrap().map(|r| r.masked).collect();
    assert_eq!(masked, [true, false]);
    let first = query(
        &mut c,
        Alleles::Record { index: 0, alt: 1 },
        Strand::Forward,
    );
    assert_eq!(first, b"ACNTACGTACGTACGTACGT");
    let second = query(
        &mut c,
        Alleles::Record { index: 1, alt: 1 },
        Strand::Forward,
    );
    assert_eq!(second, b"ACNTACGTCCGTACGTACGT");
}

#[test]
fn single_record_absent_filling_and_chain_match_a_file_of_only_that_record() {
    let dir = tempfile::tempdir().unwrap();
    let seq = b"ACGTACGTACGTACGTACGT";
    let lines = [
        "c\t3\t.\tG\tGAA\t.\t.\tDP=1\n",
        "c\t9\t.\tACG\tA\t.\t.\tDP=1\n",
    ];
    let full = bcf_of(
        dir.path(),
        "all-absent",
        format!("{HEADER}{}{}", lines[0], lines[1]).as_bytes(),
    );
    let mut genome = reference(seq);
    for (index, line) in lines.iter().enumerate() {
        let alone = bcf_of(
            dir.path(),
            &format!("alone-absent-{index}"),
            format!("{HEADER}{line}").as_bytes(),
        );
        let options = ConsensusOptions {
            absent: Some(b'N'),
            chain: true,
            iupac_codes: true,
            ..ConsensusOptions::default()
        };
        let mut cohort = open(&full, None, options.clone());
        let mut single = open(&alone, None, options);
        let region = Region::new("c", 0, seq.len() as u32);
        cohort.load(&mut genome, region).unwrap();
        single.load(&mut genome, region).unwrap();
        let got = query(
            &mut cohort,
            Alleles::Record { index, alt: 1 },
            Strand::Forward,
        );
        let expected = query(&mut single, Alleles::AllAlt, Strand::Forward);
        assert_eq!(got, expected);
        assert_eq!(cohort.chain(), single.chain());
        let record = cohort.records().unwrap().next().unwrap();
        assert!(std::panic::catch_unwind(|| record.allele(record.n_alleles())).is_err());
        assert!(std::panic::catch_unwind(|| record.allele(usize::MAX)).is_err());
    }
}
