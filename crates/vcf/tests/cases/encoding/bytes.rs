//! Header and record encoding against independent htslib bytes.
use crate::encoding_fixture::random_vcf;
use crate::encoding_oracle::{check, check_completed, check_header, header_part};
use crate::support::{integer_fields, rich_variants, variants};
use proptest::prelude::*;
use variants::Rng;

proptest! {
    #![proptest_config(ProptestConfig {
        cases: 300,
        failure_persistence: Some(Box::new(proptest::test_runner::FileFailurePersistence::Direct(
            concat!(env!("CARGO_MANIFEST_DIR"), "/tests/cases/encoding/bytes.proptest-regressions")
        ))),
        ..ProptestConfig::default()
    })]

    #[test]
    fn random_files_encode_as_htslib(seed in any::<u64>()) {
        let text = random_vcf(&mut Rng(seed));
        check_header(header_part(&text));
        check(&text);
    }

    #[test]
    fn suite_fixtures_encode_as_htslib(seed in any::<u64>()) {
        let mut rng = Rng(seed);
        check_completed(&variants::random_vcf(&mut rng));
        check_completed(&rich_variants::random_vcf(&mut rng).0);
    }
}

const COLUMNS: &str = "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO";

/// Header fixtures: each is parsed both ways on its own and, when it
/// starts as htslib needs a VCF file to, encoded with a record.
fn header_fixtures() -> Vec<String> {
    let samples = format!("{COLUMNS}\tFORMAT\tS1\n");
    vec![
        // No fileformat line: the first declaration takes id 0, PASS 1.
        format!(
            "##INFO=<ID=DP,Number=1,Type=Integer>\n##FILTER=<ID=q10>\n##INFO=<ID=AF,Number=A,Type=Float>\n##FORMAT=<ID=GT,Number=1,Type=String>\n##contig=<ID=chr1>\n{samples}"
        ),
        // A user PASS replaces htslib's, keeping its place.
        format!(
            "##fileformat=VCFv4.3\n##FILTER=<ID=q10,Description=\"q\">\n##FILTER=<ID=PASS,Description=\"Mine\">\n##INFO=<ID=DP,Number=1,Type=Integer>\n##INFO=<ID=AF,Number=A,Type=Float>\n##FORMAT=<ID=GT,Number=1,Type=String>\n##FORMAT=<ID=DP,Number=1,Type=Integer>\n##contig=<ID=chr1,length=100>\n{samples}"
        ),
        // Without a first line, PASS is first and a later PASS swaps
        // with the second record when that is a FILTER.
        format!(
            "##FILTER=<ID=q10>\n##FILTER=<ID=PASS,Description=\"Mine\">\n##contig=<ID=chr1>\n{samples}"
        ),
        format!("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"),
        // A first structured line without ID is added twice.
        format!(
            "##PEDIGREE=<Name_0=G0>\n##fileformat=VCFv4.2\n##contig=<ID=chr1>\n##FILTER=<ID=q10>\n##INFO=<ID=DP,Number=1,Type=Integer>\n##INFO=<ID=AF,Number=A,Type=Float>\n##FORMAT=<ID=GT,Number=1,Type=String>\n##FORMAT=<ID=DP,Number=1,Type=Integer>\n{samples}"
        ),
        // Spacing, quoting, escapes, brackets, trailing junk, malformed lines.
        format!(
            "##fileformat=VCFv4.4\n##contig=< id = chr1 , length = \"77\" >trailing\n##FILTER=<ID=q10,Description=\"a \\\"q\\\" , b  \" >\n##INFO=<ID=DP,Number= 1 ,Type=Integer,Description=\"\">\n##INFO=<ID=AF,Number=A,Type=Float,Values=[a, b]>\n##bad line\n##INFO=<ID=X,Number=1\n##=x\n##INFO=<ID=Y,=1>\n##X=<A=[1,2>\n##Y=<A=\"\\\\\",B=<c,d>>\n##FORMAT=<ID=GT,Number=1,Type=String>\n##FORMAT=<ID=DP,Number=1,Type=Integer>\n{samples}"
        ),
        // Shared names across kinds, duplicates, IDX, odd Number and Type.
        format!(
            "##fileformat=VCFv4.2\n##contig=<ID=chr1,IDX=5>\n##contig=<ID=chr1,IDX=6>\n##FILTER=<ID=DP>\n##INFO=<ID=DP,Number=1,Type=Integer,IDX=9>\n##INFO=<ID=DP,Number=2,Type=Float>\n##FORMAT=<ID=DP,Number=1,Type=Integer>\n##INFO=<ID=AF,Number=0x,Type=Flag>\n##FILTER=<ID=q10,IDX=3>\n##FORMAT=<ID=GT,Number=1,Type=String,IDX=7>\n##INFO=<ID=Z,Number=1,Type=Integer,IDX=x>\n##contig=<ID=chr9,IDX=4294967296>\n##source=a\n##source=a\n##source=b\n##fileformat=VCFv4.4\n{samples}"
        ),
        // Empty lines and CRLF in the header.
        format!(
            "##fileformat=VCFv4.2\r\n\r\n##contig=<ID=chr1>\r\n\n##FILTER=<ID=q10>\r\n##INFO=<ID=DP,Number=1,Type=Integer>\n##INFO=<ID=AF,Number=A,Type=Float>\r\n##FORMAT=<ID=GT,Number=1,Type=String>\n##FORMAT=<ID=DP,Number=1,Type=Integer>\n{COLUMNS}\tFORMAT\tS1\r\n"
        ),
        // Refused: conflicting IDX, no sample line, bad #CHROM lines.
        format!("##FILTER=<ID=A,IDX=0>\n{COLUMNS}\n"),
        format!(
            "##fileformat=VCFv4.2\n##INFO=<ID=A,Number=1,Type=Integer,IDX=4>\n##INFO=<ID=B,Number=1,Type=Integer,IDX=4>\n{COLUMNS}\n"
        ),
        format!("##fileformat=VCFv4.2\n#comment\n{COLUMNS}\n"),
        format!("##fileformat=VCFv4.2\n{COLUMNS}\tFORMAT\n"),
        format!("##fileformat=VCFv4.2\n{COLUMNS}\tFORMAT\tA\tA\n"),
        format!("##fileformat=VCFv4.2\n{COLUMNS}\tFORMAT\tA\t\n"),
        format!("##fileformat=VCFv4.2\n#CHROM POS ID REF ALT QUAL FILTER INFO\n"),
        "##fileformat=VCFv4.2\n".to_string(),
    ]
}

#[test]
fn header_fixtures_parse_and_encode_as_htslib() {
    let record = "chr1\t10\t.\tA\tC\t.\tq10;PASS\tDP=3;AF=0.5\tGT:DP\t0/1:7\n";
    for header in header_fixtures() {
        check_header(&header);
        check(&format!("{header}{record}"));
    }
}

#[test]
fn records_without_samples_and_with_crlf() {
    check(&format!(
        "##fileformat=VCFv4.2\r\n##contig=<ID=1>\r\n##INFO=<ID=DP,Number=1,Type=Integer>\r\n\
         {COLUMNS}\r\n1\t1\t.\tA\t.\t.\t.\t.\r\n1\t2\t.\tA\tC\t5\tPASS\tDP=300\r\n\
         1\t3\t.\tA\tC\t5\tPASS\tDP=1\tGT\t0/1\r\n\r\n1\t4\t.\tA\tC\t.\t.\t.\n"
    ));
}

#[test]
fn integer_workspace_records_encode_identically_to_htslib() {
    for version in ["4.3", "4.4"] {
        for samples in [1, 2, 17] {
            for values in [
                &[".", "-120", "127", "0"][..],
                &[".", "-121", "128", "-32760", "32767"][..],
                &[".", "-32761", "32768", "-2147483640", "2147483647"][..],
                &["."][..],
            ] {
                check(&integer_fields::fixture(
                    samples, 17, 4, 20, values, version,
                ));
            }
        }
    }
}

#[test]
fn seeded_headers_and_records_encode_as_htslib() {
    let seed = 14_209_948_683_229_990_680;
    let text = random_vcf(&mut Rng(seed));
    check_header(header_part(&text));
    check(&text);
    let mut rng = Rng(seed);
    check_completed(&variants::random_vcf(&mut rng));
    check_completed(&rich_variants::random_vcf(&mut rng).0);
}
