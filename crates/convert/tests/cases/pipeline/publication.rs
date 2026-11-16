//! Output publication and cleanup through public conversion APIs.

use super::budget;
use crate::cool_support as coolfix;
use omicsweave_convert::{BedToBigBed, ChromSizes, Conversion, CoolToHic, convert};

#[test]
fn companion_creation_failure_keeps_previous_primary_outputs() {
    let dir = tempfile::tempdir().unwrap();
    let bed = dir.path().join("in.bed");
    std::fs::write(&bed, b"chr1\t0\t10\n").unwrap();
    let cool = dir.path().join("in.cool");
    let mut spec = coolfix::MapSpec::new(&[("chr1", 100)], 10);
    spec.pixels = vec![(0, 0, 1.0)];
    coolfix::write_cool(&cool, &spec);
    let bb = dir.path().join("out.bb");
    let hic = dir.path().join("out.hic");
    for conversion in [
        Conversion::BedToBigBed(BedToBigBed::new(&bed, &bb, ChromSizes::FromData)),
        Conversion::CoolToHic(CoolToHic::new(&cool, &hic)),
    ] {
        let dst = conversion.dst().to_owned();
        std::fs::write(&dst, b"previous primary").unwrap();
        let sidecar = omicsweave_core::sidecar::sidecar_path(&dst);
        std::fs::create_dir(&sidecar).unwrap();
        let count = std::fs::read_dir(dir.path()).unwrap().count();
        let err = convert(conversion, budget(), None).unwrap_err();
        assert!(err.to_string().contains("directory"), "{err}");
        assert_eq!(std::fs::read(&dst).unwrap(), b"previous primary");
        assert!(sidecar.is_dir());
        assert_eq!(std::fs::read_dir(dir.path()).unwrap().count(), count);
    }
}

#[test]
fn invalid_records_leave_previous_outputs_and_no_temporary_files() {
    use omicsweave_convert::{FastaToBgzipFasta, SamToBam, VcfToBcf};

    let dir = tempfile::tempdir().unwrap();
    let cases: &[(&str, &str, &[u8])] = &[
        ("input.fa", "out.fa.gz", b"bases before any FASTA header\n"),
        ("input.bed", "out.bb", b"c\t1\t5\nc\t9\t3\n"),
        ("input.vcf", "out.bcf", b"##fileformat=VCFv4.2\n##contig=<ID=c,length=100>\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\nc\tnot-a-position\t.\tA\tC\t.\t.\t.\n"),
        ("input.sam", "out.bam", b"@SQ\tSN:c\tLN:100\ninvalid record\n"),
    ];
    for &(source, output, bytes) in cases {
        let src = dir.path().join(source);
        let dst = dir.path().join(output);
        std::fs::write(&src, bytes).unwrap();
        std::fs::write(&dst, b"previous output").unwrap();
        let before = std::fs::read_dir(dir.path()).unwrap().count();
        let conversion = if source.ends_with(".fa") {
            Conversion::FastaToBgzipFasta(FastaToBgzipFasta::new(&src, &dst))
        } else if source.ends_with(".bed") {
            Conversion::BedToBigBed(BedToBigBed::new(&src, &dst, ChromSizes::FromData))
        } else if source.ends_with(".vcf") {
            Conversion::VcfToBcf(VcfToBcf::new(&src, &dst))
        } else {
            Conversion::SamToBam(SamToBam::new(&src, &dst))
        };
        assert!(convert(conversion, budget(), None).is_err(), "{source}");
        assert_eq!(std::fs::read(&dst).unwrap(), b"previous output", "{source}");
        assert_eq!(
            std::fs::read_dir(dir.path()).unwrap().count(),
            before,
            "{source}: temporary files remain"
        );
    }
}

#[test]
fn a_late_publish_failure_preserves_the_per_file_commit_contract_and_cleans_staging() {
    use omicsweave_convert::stage::Outputs;

    let dir = tempfile::tempdir().unwrap();
    let primary = dir.path().join("primary");
    let companion = dir.path().join("companion");
    std::fs::write(&primary, b"previous primary").unwrap();
    let mut outputs = Outputs::new();
    outputs.write(&primary, (), b"new primary").unwrap();
    outputs.write(&companion, (), b"new companion").unwrap();
    // Publication is atomic per file; a later rename can fail after the first succeeds.
    std::fs::create_dir(&companion).unwrap();
    assert!(outputs.commit().is_err());
    assert_eq!(std::fs::read(&primary).unwrap(), b"new primary");
    assert!(companion.is_dir());
    assert_eq!(std::fs::read_dir(dir.path()).unwrap().count(), 2);
}
