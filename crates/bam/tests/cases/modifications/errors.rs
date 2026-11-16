use crate::support::modifications::write_bam;
use omicsweave_bam::{BamReader, BinIndex};
use omicsweave_core::{AlignmentFields, AlignmentRead, AlignmentsOut, Error, FileSource, Region};

#[test]
fn a_malformed_record_fails_the_query_with_its_name() {
    let dir = tempfile::tempdir().unwrap();
    let lines = vec![
        "ok\t0\tchr1\t10\t60\t4M\t*\t0\t0\tACGT\tIIII\tMM:Z:C+m,0;\tML:B:C,200".to_string(),
        "bad\t0\tchr1\t20\t60\t4M\t*\t0\t0\tACGT\tIIII\tMM:Z:C+m,3;".to_string(),
    ];
    let bam = write_bam(dir.path(), "bad", &lines);
    let bai = dir.path().join("bad.bam.bai");
    crate::support::index(&bam, &bai, 0);
    let mut reader = BamReader::open(
        FileSource::open(&bam).unwrap(),
        BinIndex::read_file(&bai).unwrap(),
    )
    .unwrap();
    let err = reader
        .query(Region::new("chr1", 0, 100), AlignmentFields::MODIFICATIONS)
        .unwrap_err();
    assert!(err.to_string().contains("bad"), "{err}");
    assert!(matches!(
        reader.write(AlignmentsOut {
            entries: &mut [],
            names: None,
            seqs: None,
            quals: None,
            cigars: None,
            aux: None,
            modifications: None,
        }),
        Err(Error::NotQueried)
    ));
    // Without modifications the record reads as usual.
    assert_eq!(
        reader
            .query(Region::new("chr1", 0, 100), AlignmentFields::AUX)
            .unwrap()
            .records,
        2
    );
    let size = reader
        .query(Region::new("chr1", 0, 15), AlignmentFields::MODIFICATIONS)
        .unwrap();
    assert_eq!((size.records, size.modifications), (1, 1));
}
