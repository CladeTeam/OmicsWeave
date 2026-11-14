mod common;

use std::collections::HashMap;
use std::sync::Arc;

use common::{Genome, Outcome, bcftools_consensus, fasta_sequence, to_indexed_bcf};
use omicsweave_consensus::{
    Alleles, Consensus, ConsensusOptions, Mask, MaskWith, Overlap, RecordFilter, VcfStore,
};
use omicsweave_core::Region;

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
    let store =
        Arc::new(VcfStore::read(std::fs::File::open(path("input.bcf")).unwrap(), None).unwrap());
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
            let mut consensus = Consensus::new(store.clone(), options);
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
                consensus.write(&mut out).unwrap();
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
fn insertion_can_restore_bases_after_the_loaded_region_is_deleted() {
    let text = b"##fileformat=VCFv4.3\n##contig=<ID=c,length=20>\n\
        ##INFO=<ID=END,Number=1,Type=Integer,Description=\"end\">\n\
        ##FORMAT=<ID=GT,Number=1,Type=String,Description=\"gt\">\n\
        #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS\n\
        c\t2\t.\tGC\tG,<DEL>\t.\t.\tEND=16\tGT\t0/2\n\
        c\t10\t.\tG\tT,<DEL>,GCT\t.\t.\t.\tGT\t0/3\n";
    for filter in [
        None,
        Some(RecordFilter::Exclude("TYPE=\"indel\"")),
        Some(RecordFilter::Include("TYPE!=\"indel\"")),
    ] {
        for lower in [false, true] {
            let store = Arc::new(VcfStore::read(&text[..], filter).unwrap());
            assert_eq!(store.n_records(), 2);
            let mut bases = vec![b'A'; 20];
            bases[1..3].copy_from_slice(b"GC");
            bases[9] = b'G';
            if lower {
                bases.make_ascii_lowercase();
            }
            let mut genome = Genome(HashMap::from([("c".to_owned(), bases)]));
            let mut consensus = Consensus::new(
                store,
                ConsensusOptions {
                    iupac_codes: true,
                    regions_overlap: Overlap::Variant,
                    ..ConsensusOptions::default()
                },
            );
            consensus
                .load(&mut genome, Region::new("c", 6, 10))
                .unwrap();
            for _ in 0..3 {
                let size = consensus
                    .query(Alleles::Haplotype {
                        sample: "S",
                        haplotype: "1",
                    })
                    .unwrap();
                let mut out = vec![0; size];
                consensus.write(&mut out).unwrap();
                assert_eq!(out, if lower { b"ct" } else { b"CT" });
            }
        }
    }
}
