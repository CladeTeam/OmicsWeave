//! Records of random VCFs, and of the BCF htslib writes from them, as
//! htslib and this crate read them: contig, position, reference span,
//! alleles, variant types and genotype values must agree.

mod common;

use std::io::Write;

use common::{Rng, htslib, htslib_to_bcf, ours, random_vcf};
use omicsweave_vcf::{Format, VariantReader};
use proptest::prelude::*;

/// Compares record by record, naming the first difference and its line.
fn same(ours: &[Option<common::Seen>], htslib: &[Option<common::Seen>], what: &str, text: &str) {
    let lines: Vec<&str> = text.lines().filter(|l| !l.starts_with('#')).collect();
    for (i, (a, b)) in ours.iter().zip(htslib).enumerate() {
        assert_eq!(a, b, "{what} record {i}: {}", lines.get(i).unwrap_or(&"?"));
    }
    assert_eq!(ours.len(), htslib.len(), "{what}: record count");
}

fn run(seed: u64) {
    let mut rng = Rng(seed);
    let text = random_vcf(&mut rng);
    let dir = tempfile::tempdir().unwrap();
    let vcf = dir.path().join("t.vcf");
    std::fs::write(&vcf, &text).unwrap();
    let want = htslib(&vcf);
    same(&ours(text.as_bytes()), &want, "vcf", &text);

    // gzip (one member) reads the same.
    let mut z = flate2::write::GzEncoder::new(Vec::new(), flate2::Compression::fast());
    z.write_all(text.as_bytes()).unwrap();
    same(&ours(&z.finish().unwrap()), &want, "vcf.gz", &text);

    let bcf = dir.path().join("t.bcf");
    htslib_to_bcf(&vcf, &bcf);
    let bytes = std::fs::read(&bcf).unwrap();
    assert_eq!(
        VariantReader::new(&bytes[..]).unwrap().format(),
        Format::Bcf
    );
    // BCF holds only the records htslib accepted.
    let kept: String = text
        .lines()
        .filter(|l| l.starts_with('#'))
        .chain(
            text.lines()
                .filter(|l| !l.starts_with('#'))
                .zip(&want)
                .filter(|(_, w)| w.is_some())
                .map(|(l, _)| l),
        )
        .map(|l| format!("{l}\n"))
        .collect();
    same(&ours(&bytes), &htslib(&bcf), "bcf", &kept);
}

proptest! {
    #![proptest_config(ProptestConfig { cases: 400, ..ProptestConfig::default() })]

    #[test]
    fn records_match_htslib(seed in any::<u64>()) {
        run(seed);
    }
}

#[test]
fn header_and_errors() {
    let text = "##fileformat=VCFv4.3\n##contig=<ID=chr1,length=100>\n\
                ##FORMAT=<ID=GT,Number=1,Type=String,Description=\"G\">\n\
                #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tA\tB\n\
                chr1\t5\t.\tA\tC\t.\t.\t.\tGT\t0|1\t1/1\n\
                chr1\tx\t.\tA\tC\t.\t.\t.\tGT\t0|1\t1/1\n\
                chr1\t7\t.\tA\tC\t.\t.\t.\tGT\t0|1\n\
                chr2\t9\t.\tA\tC\t.\t.\t.\tGT\t0|1\t1/q\n\
                chr2\t9\t.\tAG\tA\t.\t.\t.\tGT\t0\t.\n";
    let mut reader = VariantReader::new(text.as_bytes()).unwrap();
    assert_eq!(reader.header().samples, ["A", "B"]);
    assert_eq!(reader.header().version, 4_003_000);
    assert_eq!(reader.header().contigs[0].len, Some(100));
    let mut rec = omicsweave_vcf::Record::default();
    assert!(reader.read(&mut rec).unwrap());
    assert_eq!((rec.pos, rec.rlen), (4, 1));
    // A bad position, too few samples and a bad genotype are errors;
    // reading goes on with the next line.
    for _ in 0..3 {
        assert!(matches!(
            reader.read(&mut rec),
            Err(omicsweave_core::Error::Format { .. })
        ));
    }
    assert!(reader.read(&mut rec).unwrap());
    assert_eq!(reader.header().contigs[rec.rid as usize].name, "chr2");
    assert!(!reader.read(&mut rec).unwrap());

    assert!(VariantReader::new(&b"##fileformat=VCFv4.2\n"[..]).is_err());
    assert!(VariantReader::new(&b"BCF\x02\x02\x10\x00"[..]).is_err());
}
