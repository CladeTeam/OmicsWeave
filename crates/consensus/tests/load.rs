mod common;

use std::cell::Cell;
use std::sync::Arc;
use std::sync::atomic::{AtomicBool, AtomicUsize, Ordering};

use common::{Outcome, bcf_of, bcftools_consensus, consensus as open, fasta_sequence};
use omicsweave_bgzf::{BinIndex, locate_index};
use omicsweave_consensus::{Alleles, Consensus, ConsensusOptions, RecordFilter};
use omicsweave_core::{Error, FileSource, Locus, ReadAt, Region, Result, SequenceRead, Strand};
use omicsweave_vcf::BcfReader;

const HEADER: &[u8] = b"##fileformat=VCFv4.3\n##contig=<ID=c,length=16>\n\
    #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n\
    c\t3\t.\tA\tC\t.\t.\t.\n";

#[derive(Default)]
struct Reference {
    lengths: Cell<usize>,
    reads: usize,
    fail: bool,
}

impl SequenceRead for Reference {
    fn chrom_len(&self, chrom: &str) -> Result<u32> {
        self.lengths.set(self.lengths.get() + 1);
        if chrom == "c" {
            Ok(16)
        } else {
            Err(Error::UnknownChrom(chrom.to_owned()))
        }
    }

    fn sequence(&mut self, _: &Locus<'_>, out: &mut [u8]) -> Result<()> {
        self.reads += 1;
        if self.fail {
            return Err(Error::invalid("reference", "injected failure"));
        }
        // Bounds belong to the public load contract even when a provider
        // can synthesize bytes for an arbitrary range.
        out.fill(b'A');
        Ok(())
    }
}

fn assert_unavailable<S: ReadAt>(consensus: &mut Consensus<S>) {
    assert!(consensus.query(Alleles::AllAlt).is_err());
    let mut out = [0xa5; 16];
    assert!(matches!(
        consensus.write(Strand::Unstranded, &mut out),
        Err(Error::NotQueried)
    ));
    assert_eq!(out, [0xa5; 16]);
    assert!(consensus.chain().is_none());
}

fn successful_batch<S: ReadAt>(consensus: &mut Consensus<S>, reference: &mut Reference) {
    let reads = reference.reads;
    consensus.load(reference, Region::new("c", 0, 16)).unwrap();
    for _ in 0..3 {
        assert_eq!(consensus.query(Alleles::AllAlt).unwrap(), 16);
        let mut out = [0; 16];
        consensus.write(Strand::Unstranded, &mut out).unwrap();
        assert_eq!(&out, b"AACAAAAAAAAAAAAA");
        assert!(consensus.chain().is_some());
    }
    assert_eq!(reference.reads, reads + 1);
}

#[test]
fn invalid_load_rejects_old_region_and_recovers_without_reference_reads() {
    let dir = tempfile::tempdir().unwrap();
    let mut consensus = open(
        &bcf_of(dir.path(), "in", HEADER),
        None,
        ConsensusOptions {
            chain: true,
            ..ConsensusOptions::default()
        },
    );
    let mut reference = Reference::default();
    for (region, length_calls) in [
        (Region::new("c", 9, 2), 0),
        (Region::new("c", 0, 17), 1),
        (Region::new("unknown", 0, 1), 1),
    ] {
        successful_batch(&mut consensus, &mut reference);
        let (lengths, reads) = (reference.lengths.get(), reference.reads);
        assert!(consensus.load(&mut reference, region).is_err());
        assert_eq!(reference.lengths.get(), lengths + length_calls);
        assert_eq!(reference.reads, reads);
        assert_unavailable(&mut consensus);
    }
    successful_batch(&mut consensus, &mut reference);
    reference.fail = true;
    let reads = reference.reads;
    assert!(
        consensus
            .load(&mut reference, Region::new("c", 0, 16))
            .is_err()
    );
    assert_eq!(reference.reads, reads + 1);
    assert_unavailable(&mut consensus);
    reference.fail = false;
    successful_batch(&mut consensus, &mut reference);
}

#[test]
fn bounds_are_checked_before_large_region_allocation_or_provider_access() {
    let dir = tempfile::tempdir().unwrap();
    let mut consensus = open(
        &bcf_of(dir.path(), "in", HEADER),
        None,
        ConsensusOptions::default(),
    );
    let mut reference = Reference::default();
    assert!(
        consensus
            .load(&mut reference, Region::new("c", 0, u32::MAX))
            .is_err()
    );
    assert_eq!(reference.lengths.get(), 1);
    assert_eq!(reference.reads, 0);
    assert_unavailable(&mut consensus);
}

struct Source {
    inner: FileSource,
    fail: Arc<AtomicBool>,
    reads: Arc<AtomicUsize>,
}

impl ReadAt for Source {
    fn len(&self) -> u64 {
        self.inner.len()
    }
    fn read_exact_at(&self, buf: &mut [u8], offset: u64) -> Result<()> {
        self.reads.fetch_add(1, Ordering::Relaxed);
        if self.fail.load(Ordering::Relaxed) {
            return Err(Error::invalid("BCF", "injected read failure"));
        }
        self.inner.read_exact_at(buf, offset)
    }
}

#[test]
fn bcf_read_failure_discards_loaded_region_and_a_retry_rereads_the_source() {
    let dir = tempfile::tempdir().unwrap();
    let bcf = bcf_of(dir.path(), "in", HEADER);
    let fail = Arc::new(AtomicBool::new(false));
    let reads = Arc::new(AtomicUsize::new(0));
    let source = Source {
        inner: FileSource::open(&bcf).unwrap(),
        fail: fail.clone(),
        reads: reads.clone(),
    };
    let reader = BcfReader::open(
        source,
        BinIndex::read_file(locate_index(&bcf).unwrap()).unwrap(),
    )
    .unwrap();
    let mut consensus = Consensus::new(
        reader,
        None,
        ConsensusOptions {
            chain: true,
            ..ConsensusOptions::default()
        },
    )
    .unwrap();
    let mut reference = Reference::default();
    successful_batch(&mut consensus, &mut reference);
    let before = reads.load(Ordering::Relaxed);
    fail.store(true, Ordering::Relaxed);
    assert!(
        consensus
            .load(&mut reference, Region::new("c", 0, 16))
            .is_err()
    );
    assert!(reads.load(Ordering::Relaxed) > before);
    assert_unavailable(&mut consensus);
    fail.store(false, Ordering::Relaxed);
    let before = reads.load(Ordering::Relaxed);
    successful_batch(&mut consensus, &mut reference);
    assert!(reads.load(Ordering::Relaxed) > before);
}

#[test]
fn filters_keep_record_state_inside_a_load_and_start_fresh_for_each_region() {
    let dir = tempfile::tempdir().unwrap();
    let text = b"##fileformat=VCFv4.3\n##contig=<ID=c,length=16>\n\
        ##INFO=<ID=AF,Number=1,Type=Float,Description=\"frequency\">\n\
        ##FORMAT=<ID=LEN,Number=1,Type=String,Description=\"text\">\n\
        #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS\n\
        c\t1\t.\tA\tC\t1\t.\tAF=5.266\tLEN\tABCDEFGH\n\
        c\t2\t.\tA\tC\t1\t.\tAF=5.266\tLEN\tY-b1\n";
    let bcf = bcf_of(dir.path(), "in", text);
    let expression = "STRLEN(FMT/LEN)>AF";
    let mut consensus = open(
        &bcf,
        Some(RecordFilter::Include(expression)),
        ConsensusOptions::default(),
    );
    let mut reference = Reference::default();
    for (start, end, expected) in [
        (0, 2, b"CC".as_slice()),
        (1, 2, b"A"),
        (0, 1, b"C"),
        (1, 2, b"A"),
    ] {
        consensus
            .load(&mut reference, Region::new("c", start, end))
            .unwrap();
        let n = consensus.query(Alleles::AllAlt).unwrap();
        let mut actual = vec![0; n];
        consensus.write(Strand::Unstranded, &mut actual).unwrap();
        assert_eq!(actual, expected);
        let fasta = dir.path().join("window.fa");
        let output = dir.path().join("out.fa");
        std::fs::write(
            &fasta,
            format!(
                ">c:{}-{end}\n{}\n",
                start + 1,
                "A".repeat((end - start) as usize)
            ),
        )
        .unwrap();
        let args = [
            "-s".to_owned(),
            "-".to_owned(),
            "-i".to_owned(),
            expression.to_owned(),
            "-f".to_owned(),
            fasta.display().to_string(),
            "-o".to_owned(),
            output.display().to_string(),
            bcf.display().to_string(),
        ];
        assert_eq!(bcftools_consensus(&args), Outcome::Done(0));
        assert_eq!(
            actual,
            fasta_sequence(&std::fs::read_to_string(&output).unwrap())
        );
    }
}

#[test]
fn clipped_insertion_after_whole_window_deletion_returns_an_error_and_recovers() {
    let dir = tempfile::tempdir().unwrap();
    let text = b"##fileformat=VCFv4.3\n##contig=<ID=c,length=64>\n\
        ##INFO=<ID=END,Number=1,Type=Integer,Description=\"end\">\n\
        #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n\
        c\t32\t.\tTACA\t*\t.\t.\t.\n\
        c\t35\t.\tA\tAA,<DEL>\t.\t.\tEND=37\n";
    let mut bases = vec![b'A'; 64];
    bases[31..35].copy_from_slice(b"TACA");
    let mut reference = common::Genome(std::collections::HashMap::from([("c".to_owned(), bases)]));
    let mut consensus = open(
        &bcf_of(dir.path(), "in", text),
        None,
        ConsensusOptions {
            chain: true,
            ..ConsensusOptions::default()
        },
    );
    consensus
        .load(&mut reference, Region::new("c", 32, 35))
        .unwrap();
    for _ in 0..2 {
        assert!(matches!(
            consensus.query(Alleles::AllAlt),
            Err(Error::Unsupported(_))
        ));
        assert!(consensus.write(Strand::Unstranded, &mut []).is_err());
        assert!(consensus.chain().is_none());
    }
    consensus
        .load(&mut reference, Region::new("c", 0, 1))
        .unwrap();
    assert_eq!(consensus.query(Alleles::AllAlt).unwrap(), 1);
    let mut output = [0];
    consensus.write(Strand::Unstranded, &mut output).unwrap();
    assert_eq!(&output, b"A");
}
