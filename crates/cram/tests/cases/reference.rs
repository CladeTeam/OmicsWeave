use crate::support::*;
use omicsweave_core::{
    AlignmentFields, AlignmentRead, AlignmentsOut, CigarColumn, Error, FileSource, Region,
};
use omicsweave_cram::{CramIndex, CramOptions, CramReader, NoReference};

#[test]
fn embedded_and_absent_references() {
    let f = fixture(
        2,
        standard_refs(),
        SamSpec {
            per_ref: 1_500,
            gap: 200,
            long_cigars: false,
            unplaced: 5,
        },
    );
    let (cram, crai) = f.cram("embed", &["version=3.1", "embed_ref=1"], true);
    assert!(check_queries(&f, &cram, &crai, false, 7) > 1_000);
    let (cram, crai) = f.cram("noref", &["version=3.0", "no_ref=1"], false);
    assert!(check_queries(&f, &cram, &crai, false, 8) > 1_000);
}

#[test]
fn fixed_fields_need_no_reference() {
    let f = fixture(
        5,
        standard_refs(),
        SamSpec {
            per_ref: 1_000,
            gap: 80,
            long_cigars: true,
            unplaced: 0,
        },
    );
    let (cram, crai) = f.cram("plain", &["version=3.1"], true);
    let hts = Hts::open(&cram, &crai, Some(&f.fasta));
    let mut bare: CramReader<FileSource, NoReference> = {
        let header = omicsweave_cram::CramHeader::read(&FileSource::open(&cram).unwrap()).unwrap();
        let index = CramIndex::read_file(&crai, header.references().len()).unwrap();
        CramReader::open(
            FileSource::open(&cram).unwrap(),
            index,
            None,
            CramOptions::default(),
        )
        .unwrap()
    };
    let fields = AlignmentFields::CIGAR | AlignmentFields::QUAL;
    for (chrom, len) in &f.refs {
        let tid = hts.tid(chrom);
        let want = hts.query(tid, 0, *len as i64);
        let got = query(&mut bare, tid, Region::new(chrom, 0, *len as u32), fields);
        assert_eq!(got, only(&want, fields), "{chrom}");
    }
    let err = bare
        .query(Region::new("chr1", 0, 600_000), AlignmentFields::SEQ)
        .unwrap_err();
    assert!(matches!(err, Error::Invalid { .. }), "{err}");
}

#[test]
fn reference_checks_and_errors() {
    let f = fixture(
        6,
        standard_refs(),
        SamSpec {
            per_ref: 300,
            gap: 80,
            long_cigars: false,
            unplaced: 0,
        },
    );
    let (cram, crai) = f.cram("x", &["version=3.0"], true);
    // A different reference fails the slice MD5 check, unless told not to.
    let mut rng = Rng::new(99);
    let other = reference(&mut rng, &f.refs);
    let wrong = f.dir.path().join("wrong.fa");
    write_fasta(&wrong, &f.refs, &other);
    let mut r = open(&cram, &crai, Some(&wrong));
    let err = r
        .query(Region::new("chr1", 0, 600_000), AlignmentFields::SEQ)
        .unwrap_err();
    assert!(err.to_string().contains("MD5"), "{err}");
    assert!(
        r.query(Region::new("chr1", 0, 600_000), AlignmentFields::NONE)
            .is_ok()
    );
    assert!(matches!(
        r.query(Region::new("nope", 0, 10), AlignmentFields::NONE),
        Err(Error::UnknownChrom(_))
    ));
    assert_eq!(
        r.query(Region::new("chr1", 500, 500), AlignmentFields::ALL)
            .unwrap()
            .records,
        0
    );
    assert!(
        r.write(AlignmentsOut {
            entries: &mut [],
            names: None,
            seqs: None,
            quals: None,
            cigars: Some(CigarColumn {
                offsets: &mut [0],
                ops: &mut []
            }),
            aux: None,
            modifications: None,
        })
        .is_err()
    );
}
