mod common;

use std::sync::Arc;
use std::sync::atomic::{AtomicUsize, Ordering};

use common::{bcf_of, consensus as open};
use omicsweave_consensus::{Alleles, Consensus, ConsensusOptions, Mask, MaskWith};
use omicsweave_core::{Error, Locus, Region, Result, SequenceRead, Strand};

struct Reference(Arc<AtomicUsize>);

impl SequenceRead for Reference {
    fn chrom_len(&self, chrom: &str) -> Result<u32> {
        match chrom {
            "c" | "d" => Ok(16),
            _ => Err(Error::UnknownChrom(chrom.to_owned())),
        }
    }

    fn sequence(&mut self, locus: &Locus<'_>, out: &mut [u8]) -> Result<()> {
        locus.check_within(self.chrom_len(locus.chrom)?)?;
        assert_eq!(out.len(), locus.len() as usize);
        self.0.fetch_add(1, Ordering::Relaxed);
        out.fill(if locus.chrom == "c" { b'A' } else { b'C' });
        Ok(())
    }
}

fn query(consensus: &mut Consensus, alleles: Alleles<'_>) -> Vec<u8> {
    let n = consensus.query(alleles).unwrap();
    let mut out = vec![0; n];
    consensus.write(Strand::Unstranded, &mut out).unwrap();
    out
}

#[test]
fn cloned_contexts_share_options_and_keep_reference_output_and_chain_independent() {
    let dir = tempfile::tempdir().unwrap();
    let vcf = b"##fileformat=VCFv4.3\n\
        ##contig=<ID=c,length=16>\n##contig=<ID=d,length=16>\n\
        ##FORMAT=<ID=GT,Number=1,Type=String,Description=\"gt\">\n\
        #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS\n\
        c\t3\t.\tA\tT\t.\t.\t.\tGT\t0|1\n\
        c\t7\t.\tAAA\tA\t.\t.\t.\tGT\t1|0\n\
        d\t4\t.\tC\tCGG\t.\t.\t.\tGT\t1|1\n";
    let mask = Arc::new(Mask::parse(b"c\t10\t12\n", MaskWith::Lower).unwrap());
    let mut original = open(
        &bcf_of(dir.path(), "in", vcf),
        None,
        ConsensusOptions {
            chain: true,
            mask: Some(mask.clone()),
            ..ConsensusOptions::default()
        },
    );
    let reads = Arc::new(AtomicUsize::new(0));
    original
        .load(&mut Reference(reads.clone()), Region::new("c", 0, 16))
        .unwrap();
    let hap2 = Alleles::Haplotype {
        sample: "S",
        haplotype: "2",
    };
    let expected_original = b"AATAAAAAAAaaAAAA";
    assert_eq!(query(&mut original, hap2), expected_original);
    let original_chain = original.chain().unwrap().to_owned();

    let mut cloned = original.try_clone().unwrap();
    assert_eq!(cloned.samples(), ["S"]);
    assert!(Arc::ptr_eq(cloned.options().mask.as_ref().unwrap(), &mask));
    assert_eq!(reads.load(Ordering::Relaxed), 1);
    assert!(cloned.query(Alleles::AllAlt).is_err());
    assert!(cloned.write(Strand::Unstranded, &mut []).is_err());
    assert!(cloned.chain().is_none());

    let mut cloned = std::thread::scope(|scope| {
        let child_reads = reads.clone();
        let child = scope.spawn(move || {
            cloned
                .load(&mut Reference(child_reads), Region::new("d", 1, 12))
                .unwrap();
            assert_eq!(query(&mut cloned, hap2), b"CCCGGCCCCCCCC");
            assert!(cloned.chain().unwrap().contains(" d "));
            cloned
        });
        assert_eq!(query(&mut original, hap2), expected_original);
        assert_eq!(original.chain().unwrap(), original_chain);
        child.join().unwrap()
    });

    drop(original);
    drop(mask);
    assert_eq!(query(&mut cloned, hap2), b"CCCGGCCCCCCCC");
    cloned
        .load(&mut Reference(reads.clone()), Region::new("c", 0, 16))
        .unwrap();
    assert_eq!(query(&mut cloned, hap2), expected_original);
    assert_eq!(reads.load(Ordering::Relaxed), 3);
}

#[test]
fn sample_sets_preserve_iupac_selection_across_inline_and_large_sets() {
    let dir = tempfile::tempdir().unwrap();
    let vcf = b"##fileformat=VCFv4.3\n##contig=<ID=c,length=16>\n\
        ##FORMAT=<ID=GT,Number=1,Type=String,Description=\"gt\">\n\
        #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\tS2\tS3\tS4\tS5\n\
        c\t3\t.\tA\tC\t.\t.\t.\tGT\t0/1\t1/1\t0/0\t0|1\t1|0\n";
    let mut consensus = open(
        &bcf_of(dir.path(), "in", vcf),
        None,
        ConsensusOptions::default(),
    );
    let reads = Arc::new(AtomicUsize::new(0));
    consensus
        .load(&mut Reference(reads.clone()), Region::new("c", 0, 16))
        .unwrap();
    let names = ["S1", "S2", "S3", "S4", "S5"];
    for n in [1, 4, 5, 2, 5, 1] {
        assert_eq!(
            query(&mut consensus, Alleles::Samples(&names[..n])),
            b"AAMAAAAAAAAAAAAA"
        );
    }
    let repeated: Vec<_> = names.iter().copied().cycle().take(32_768).collect();
    assert_eq!(
        query(&mut consensus, Alleles::Samples(&repeated)),
        b"AAMAAAAAAAAAAAAA"
    );
    assert_eq!(
        query(&mut consensus, Alleles::Samples(&["S2"])),
        b"AACAAAAAAAAAAAAA"
    );
    assert!(consensus.query(Alleles::Samples(&[])).is_err());
    assert!(consensus.query(Alleles::Samples(&["unknown"])).is_err());
    assert_eq!(query(&mut consensus, Alleles::AllAlt), b"AACAAAAAAAAAAAAA");
    assert_eq!(reads.load(Ordering::Relaxed), 1);
}
