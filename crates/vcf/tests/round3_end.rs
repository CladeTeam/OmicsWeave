//! Explicit END is distinct from the effective span and must survive both formats.
mod common;

use omicsweave_vcf::{Record, VariantReader};

const HEADER: &str = "##fileformat=VCFv4.5\n##contig=<ID=chr1,length=100000>\n##INFO=<ID=END,Number=1,Type=Integer,Description=\"END\">\n##INFO=<ID=SVLEN,Number=A,Type=Integer,Description=\"SVLEN\">\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n";

#[test]
fn explicit_end_is_preserved_separately_from_effective_span() {
    let text = format!(
        "{HEADER}chr1\t10\t.\tA\tC\t.\tPASS\tEND=25\nchr1\t30\t.\tA\t<DEL>\t.\tPASS\tEND=31;SVLEN=1000\nchr1\t50\t.\tAT\tA\t.\tPASS\t.\n"
    );
    let dir = tempfile::tempdir().unwrap();
    let vcf = dir.path().join("end.vcf");
    let bcf = dir.path().join("end.bcf");
    std::fs::write(&vcf, &text).unwrap();
    common::htslib_to_bcf(&vcf, &bcf);
    for bytes in [text.as_bytes().to_vec(), std::fs::read(bcf).unwrap()] {
        let mut reader = VariantReader::new(&bytes[..]).unwrap();
        let mut rec = Record::default();
        // A symbolic deletion's SVLEN excludes the one-base anchor.
        for (end, span) in [(Some(25), 16), (Some(31), 1001), (None, 2)] {
            assert!(reader.read(&mut rec).unwrap());
            assert_eq!(rec.end, end);
            assert_eq!(rec.rlen, span);
        }
        assert!(!reader.read(&mut rec).unwrap());
        assert_eq!(rec, Record::default());
    }
}

#[test]
fn explicit_end_is_cleared_after_a_failed_text_record() {
    let text = format!(
        "{HEADER}chr1\t10\t.\tA\tC\t.\tPASS\tEND=25\nchr1\tbad\t.\tA\tC\t.\tPASS\tEND=50\nchr1\t50\t.\tAT\tA\t.\tPASS\t.\n"
    );
    let mut reader = VariantReader::new(text.as_bytes()).unwrap();
    let mut rec = Record::default();
    assert!(reader.read(&mut rec).unwrap());
    assert_eq!(rec.end, Some(25));
    assert!(reader.read(&mut rec).is_err());
    assert_eq!(rec, Record::default());
    assert!(reader.read(&mut rec).unwrap());
    assert_eq!(rec.end, None);
}
