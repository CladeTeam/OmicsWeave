use omicsweave_core::{ContactQuery, ContactRead, ContactUnit, Region};
use omicsweave_hic::{BuildOptions, Chrom, HicBuilder, HicHeader, HicReader, HicWriter};
use std::io::Cursor;

fn header() -> HicHeader {
    let mut header = HicHeader::new(
        "g",
        vec![
            Chrom {
                name: "a".into(),
                len: 100,
            },
            Chrom {
                name: "b".into(),
                len: 100,
            },
        ],
        vec![1],
    );
    header.frag_resolutions = vec![1];
    header.sites = vec![vec![10, 20], vec![10, 20]];
    header
}

#[test]
fn duplicate_resolutions_are_rejected_before_collecting_contacts() {
    for unit in [ContactUnit::Bp, ContactUnit::Frag] {
        let mut header = header();
        let resolutions = match unit {
            ContactUnit::Bp => &mut header.bp_resolutions,
            ContactUnit::Frag => &mut header.frag_resolutions,
        };
        *resolutions = vec![1, 2, 1];
        assert!(HicBuilder::new(header.clone()).is_err(), "{unit:?}");
        assert!(
            HicWriter::create(Cursor::new(Vec::new()), header, 1).is_err(),
            "{unit:?}"
        );
    }
}

#[test]
fn equal_bp_and_fragment_resolutions_store_each_contact_once() {
    let mut builder = HicBuilder::new(header()).unwrap();
    builder.add(0, 5, 1, 15, 3.0).unwrap();
    let bytes = builder
        .finish(
            Cursor::new(Vec::new()),
            &BuildOptions {
                level: 1,
                ..Default::default()
            },
        )
        .unwrap()
        .into_inner();
    let mut reader = HicReader::open(&bytes[..]).unwrap();
    // The same contact is BP cell (5,15), and fragment cell (0,1), since
    // restriction sites at 10 and 20 bracket its two original positions.
    for (unit, x, y) in [(ContactUnit::Bp, 5, 15), (ContactUnit::Frag, 0, 1)] {
        assert_eq!(reader.resolutions(unit), &[1]);
        let query = ContactQuery::new(Region::new("a", x, x + 1), Region::new("b", y, y + 1), 1)
            .with_unit(unit);
        let mut out = [0.0];
        reader.matrix(query, &mut out).unwrap();
        assert_eq!(out, [3.0], "{unit:?}");
    }
}
