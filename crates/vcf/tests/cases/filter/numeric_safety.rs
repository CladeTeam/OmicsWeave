//! Finite FORMAT masks must not expose C's reads beyond initialized values.

use omicsweave_core::Error;
use omicsweave_vcf::{Filter, Record, VariantReader};

#[test]
fn numeric_format_range_past_sample_width_is_rejected_and_filter_recovers() {
    for kind in ["Integer", "Float"] {
        let text = format!(
            "##fileformat=VCFv4.4\n##contig=<ID=chr1>\n\
             ##FORMAT=<ID=DP,Number=.,Type={kind}>\n\
             #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tA\tB\n\
             chr1\t1\t.\tA\tC\t.\tPASS\t.\tDP\t57\t60\n\
             chr1\t2\t.\tA\tC\t.\tPASS\t.\tDP\t1,2\t3,30\n"
        );
        let mut reader = VariantReader::new(text.as_bytes()).unwrap();
        let mut filter = Filter::new(reader.header(), "FMT/DP[1:0-1]=30").unwrap();
        let mut record = Record::default();
        assert!(reader.read(&mut record).unwrap());
        // filter.c's finite-index path reads src[1] for sample B when each
        // sample has width 1. This is past the getter's initialized buffer.
        assert!(
            matches!(
                filter.test(reader.header(), &record),
                Err(Error::Unsupported(_))
            ),
            "{kind}"
        );
        assert!(reader.read(&mut record).unwrap());
        assert!(filter.test(reader.header(), &record).unwrap(), "{kind}");
    }
}

#[test]
fn multi_sample_fisher_results_require_initialized_format_declarations() {
    let expression = "fisher(FMT/DP[:0-1],FMT/GQ[:0-1])>0.001";
    for number in [".", "R"] {
        let text = format!(
            "##fileformat=VCFv4.4\n##contig=<ID=chr1>\n\
             ##FORMAT=<ID=DP,Number={number},Type=Integer>\n\
             ##FORMAT=<ID=GQ,Number={number},Type=Integer>\n\
             #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS0\tS1\n\
             chr1\t1\t.\tA\tC,G\t30\tPASS\t.\tDP:GQ\t1,9,0:11,3,0\t1,9,0:11,3,0\n"
        );
        let mut reader = VariantReader::new(text.as_bytes()).unwrap();
        let mut filter = Filter::new(reader.header(), expression).unwrap();
        let mut record = Record::default();
        assert!(reader.read(&mut record).unwrap());
        if number == "." {
            // filter_init zeroes the first result slot. Additional sample
            // slots stay unwritten for two non-Number=R FORMAT arguments.
            assert!(matches!(
                filter.test(reader.header(), &record),
                Err(Error::Unsupported(_))
            ));
        } else {
            assert!(filter.test(reader.header(), &record).unwrap());
            let directory = tempfile::tempdir().unwrap();
            let path = directory.path().join("fisher.vcf");
            std::fs::write(&path, &text).unwrap();
            assert_eq!(
                crate::filter_oracle::c_filter(&path, expression),
                crate::filter_oracle::Run::Records(vec![crate::filter_oracle::Outcome::Pass(true)])
            );
        }
    }
}
