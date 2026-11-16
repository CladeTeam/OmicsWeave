use crate::support::{fixture, open};
use omicsweave_core::{AlignmentFields, AlignmentRead, AlignmentsOut, Error, Region};

#[test]
fn unknown_and_empty_regions() {
    let f = fixture(8, 500, 100, 6);
    let mut r = open(&f.bam, &f.bai);
    assert!(matches!(
        r.query(Region::new("nope", 0, 10), AlignmentFields::NONE),
        Err(Error::UnknownChrom(_))
    ));
    assert_eq!(
        r.query(Region::new("chr1", 500, 500), AlignmentFields::NONE)
            .unwrap()
            .records,
        0
    );
    assert!(matches!(
        r.query(Region::new("chr1", 5, 1), AlignmentFields::NONE),
        Err(Error::InvalidRegion { .. })
    ));
    // A failed query leaves nothing to write.
    let mut entries = Vec::new();
    assert!(matches!(
        r.write(AlignmentsOut {
            entries: &mut entries,
            names: None,
            seqs: None,
            quals: None,
            cigars: None,
            aux: None,
            modifications: None,
        }),
        Err(Error::NotQueried)
    ));
    assert_eq!(r.chrom_len("chrM").unwrap(), 16_569);
}

#[test]
fn query_failures_discard_previously_staged_records_and_allow_reuse() {
    let f = fixture(83, 12, 30, 0);
    let mut reader = open(&f.bam, &f.bai);
    let valid = Region::new("chr1", 0, 3_000_000);
    for invalid in [Region::new("absent", 0, 1), Region::new("chr1", 10, 1)] {
        assert!(reader.query(valid, AlignmentFields::NONE).unwrap().records > 0);
        assert!(
            reader.query(invalid, AlignmentFields::NAME).is_err(),
            "{invalid:?}"
        );
        assert!(
            matches!(
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
            ),
            "{invalid:?}"
        );
        let size = reader.query(valid, AlignmentFields::NONE).unwrap();
        let mut entries = vec![Default::default(); size.records];
        reader
            .write(AlignmentsOut {
                entries: &mut entries,
                names: None,
                seqs: None,
                quals: None,
                cigars: None,
                aux: None,
                modifications: None,
            })
            .unwrap();
    }
}

#[test]
fn visitor_failure_stops_iteration_and_next_query_recovers() {
    let f = fixture(84, 12, 30, 0);
    let mut reader = open(&f.bam, &f.bai);
    let region = Region::new("chr1", 0, 3_000_000);
    let mut visits = 0;
    let err = reader
        .visit(region, |_| {
            visits += 1;
            Err(Error::invalid("visitor", "injected failure"))
        })
        .unwrap_err();
    assert_eq!(visits, 1);
    assert!(matches!(
        err,
        Error::Invalid {
            what: "visitor",
            ..
        }
    ));
    assert!(!crate::support::query_all(&mut reader, region).is_empty());
}
