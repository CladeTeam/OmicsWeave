//! FORMAT string capacity can contain bytes the C getter never initialized.

use omicsweave_core::Error;
use omicsweave_vcf::{Filter, Record, VariantReader};

#[test]
fn growing_format_string_rejects_unwritten_terminator_and_recovers_after_overwrite() {
    let text = format!(
        "##fileformat=VCFv4.4\n##contig=<ID=chr1>\n\
         ##FORMAT=<ID=GL,Number=1,Type=String>\n\
         #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS0\n\
         chr1\t1\t.\tA\tC\t.\tPASS\t.\tGL\tabc\n\
         chr1\t2\t.\tA\tC\t.\tPASS\t.\tGL\tabcde\n\
         chr1\t3\t.\tA\tC\t.\tPASS\t.\tGL\t{}\n\
         chr1\t4\t.\tA\tC\t.\tPASS\t.\tGL\tabcde\n",
        "x".repeat(64),
    );
    let mut reader = VariantReader::new(text.as_bytes()).unwrap();
    let mut filter = Filter::new(reader.header(), "FMT/GL=\"never\"").unwrap();
    let mut record = Record::default();
    assert!(reader.read(&mut record).unwrap());
    assert!(!filter.test(reader.header(), &record).unwrap());
    assert!(reader.read(&mut record).unwrap());
    assert!(matches!(
        filter.test(reader.header(), &record),
        Err(Error::Unsupported(_))
    ));
    for _ in 0..2 {
        assert!(reader.read(&mut record).unwrap());
        assert!(!filter.test(reader.header(), &record).unwrap());
    }
}
