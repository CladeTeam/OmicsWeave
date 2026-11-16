//! Metadata must agree with the input before output is replaced.

use super::budget;
use crate::cool_support as coolfix;
use omicsweave_convert::{
    BedToBigBed, ChromSizes, Conversion, CoolToHic, FastaToTwoBit, convert, plan,
};
use omicsweave_fasta::{TwoBitWriter, encode_twobit};
use std::io::Cursor;

#[test]
fn fai_names_must_match_input_order_even_when_lengths_match() {
    for text in [
        b">new\nACGT\n>b\nTGCA\n".as_slice(),
        b">b\nACGT\n>a\nTGCA\n",
    ] {
        let dir = tempfile::tempdir().unwrap();
        let src = dir.path().join("in.fa");
        let dst = dir.path().join("out.2bit");
        std::fs::write(&src, text).unwrap();
        std::fs::write(
            dir.path().join("in.fa.fai"),
            b"a\t4\t3\t4\t5\nb\t4\t11\t4\t5\n",
        )
        .unwrap();
        std::fs::write(&dst, b"previous output").unwrap();
        let err = convert(
            Conversion::FastaToTwoBit(FastaToTwoBit::new(&src, &dst)),
            budget(),
            None,
        )
        .unwrap_err();
        assert!(err.to_string().contains("out of date"), "{err}");
        assert_eq!(std::fs::read(&dst).unwrap(), b"previous output");
        assert_eq!(std::fs::read_dir(dir.path()).unwrap().count(), 3);
    }
}

#[test]
fn twobit_chromosome_sizes_use_lengths_from_record_metadata() {
    let dir = tempfile::tempdir().unwrap();
    let sizes = dir.path().join("genome.2bit");
    let records = [("chrA", b"NNacgT".as_slice()), ("chrB", b"ACGTNNacg")];
    let mut writer = TwoBitWriter::create(
        Cursor::new(Vec::new()),
        records.iter().map(|(name, _)| name.to_string()).collect(),
        0,
    )
    .unwrap();
    for (name, bases) in records {
        let (record, packed) = encode_twobit(bases).unwrap();
        writer.add(name, &record, &mut packed.as_slice()).unwrap();
    }
    std::fs::write(&sizes, writer.finish().unwrap().into_inner()).unwrap();
    let src = dir.path().join("in.bed");
    let dst = dir.path().join("out.bb");
    std::fs::write(&src, b"chrB\t0\t9\nchrA\t0\t6\n").unwrap();
    convert(
        Conversion::BedToBigBed(BedToBigBed::new(
            &src,
            &dst,
            ChromSizes::File(sizes.clone()),
        )),
        budget(),
        None,
    )
    .unwrap();
    let previous = std::fs::read(&dst).unwrap();
    std::fs::write(&src, b"chrA\t0\t7\n").unwrap();
    assert!(
        convert(
            Conversion::BedToBigBed(BedToBigBed::new(&src, &dst, ChromSizes::File(sizes))),
            budget(),
            None,
        )
        .is_err()
    );
    assert_eq!(std::fs::read(&dst).unwrap(), previous);
}

#[test]
fn empty_mcool_selections_are_errors_before_execution() {
    let dir = tempfile::tempdir().unwrap();
    let src = dir.path().join("in.mcool");
    let dst = dir.path().join("out.hic");
    for specs in [vec![], vec![coolfix::MapSpec::new(&[("chr1", 100)], 10)]] {
        coolfix::write_mcool(&src, &specs);
        let mut conversion = CoolToHic::new(&src, &dst);
        if !specs.is_empty() {
            conversion.resolutions = Some(vec![]);
        }
        let err = plan(Conversion::CoolToHic(conversion), budget()).unwrap_err();
        assert!(err.to_string().contains("no resolutions"), "{err}");
        assert!(!dst.exists());
    }
}
