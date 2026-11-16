use std::io::Cursor;

use omicsweave_core::{
    ContactEntry, ContactQuery, ContactRead, ContactUnit, Error, FileSource, Region,
};
use omicsweave_hic::{Chrom, HicHeader, HicReader, HicWriter, Record};

use crate::support::io::Counted;

fn fixture(fragments: bool) -> Vec<u8> {
    let mut header = HicHeader::new(
        "fragments",
        vec![
            Chrom {
                name: "a".into(),
                len: 100,
            },
            Chrom {
                name: "b".into(),
                len: 200,
            },
        ],
        vec![1],
    );
    if fragments {
        header.frag_resolutions = vec![1];
        header.sites = vec![vec![0, 17, 17, 100], vec![]];
    }
    let mut writer = HicWriter::create(Cursor::new(Vec::new()), header, 1).unwrap();
    writer
        .write_matrix(
            0,
            1,
            ContactUnit::Bp,
            1,
            &[Record {
                x: 2,
                y: 3,
                value: 7.0,
            }],
        )
        .unwrap();
    writer.finish().unwrap().into_inner()
}

#[test]
fn fragment_counts_use_metadata_and_sites_are_cached_per_context() {
    let source = Counted::new(fixture(true));
    let mut reader = HicReader::open(source.clone()).unwrap();
    let before = source.reads();
    assert_eq!(reader.fragment_count("a").unwrap(), 4);
    assert_eq!(reader.fragment_count("b").unwrap(), 0);
    assert!(matches!(
        reader.fragment_count("missing"),
        Err(Error::UnknownChrom(_))
    ));
    assert_eq!(source.reads(), before, "counts must not load site payloads");
    assert_eq!(reader.sites(0).unwrap(), &[0, 17, 17, 100]);
    let loaded = source.reads();
    assert!(loaded > before);
    assert_eq!(reader.sites(0).unwrap(), &[0, 17, 17, 100]);
    let mut clone = reader.clone_context();
    assert_eq!(clone.sites(0).unwrap(), &[0, 17, 17, 100]);
    assert_eq!(
        source.reads(),
        loaded,
        "loaded site vectors are shared by context clones"
    );
    assert!(reader.sites(2).is_err());
    let bytes = fixture(false);
    let mut no_fragments = HicReader::open(bytes.as_slice()).unwrap();
    assert!(no_fragments.fragment_count("a").is_err());
    assert!(no_fragments.sites(0).is_err());
}

fn assert_fresh_context<S: omicsweave_core::ReadAt>(reader: &mut HicReader<S>) {
    assert!(matches!(
        reader.write_entries(&mut []),
        Err(Error::NotQueried)
    ));
    assert_eq!(reader.block_inflations(), 0);
    assert_eq!(reader.cache_preparations(), 0);
    assert_eq!(reader.cached_records_visited(), 0);
}

#[test]
fn memory_and_file_context_clones_have_independent_staged_entries() {
    let bytes = fixture(true);
    let query = ContactQuery::new(Region::new("a", 0, 10), Region::new("b", 0, 10), 1);
    let mut reader = HicReader::open(bytes.as_slice()).unwrap();
    assert_eq!(reader.query_entries(query).unwrap(), 1);
    let mut clone = reader.clone_context();
    assert_eq!(clone.cache_namespace(), reader.cache_namespace());
    assert_fresh_context(&mut clone);
    assert_eq!(
        clone
            .query_entries(query.with_unit(ContactUnit::Frag))
            .unwrap(),
        0
    );
    clone.write_entries(&mut []).unwrap();
    let invalid = ContactQuery {
        rows: Region::new("missing", 0, 10),
        ..query
    };
    assert!(matches!(
        clone.query_entries(invalid),
        Err(Error::UnknownChrom(_))
    ));
    assert!(matches!(
        clone.write_entries(&mut []),
        Err(Error::NotQueried)
    ));
    let expected = [ContactEntry {
        row: 2,
        col: 3,
        value: 7.0,
    }];
    let mut actual = [ContactEntry::default()];
    reader.write_entries(&mut actual).unwrap();
    assert_eq!(actual, expected);
    assert_eq!(clone.query_entries(query).unwrap(), 1);
    clone.write_entries(&mut actual).unwrap();
    assert_eq!(actual, expected);

    let dir = tempfile::tempdir().unwrap();
    let path = dir.path().join("contacts.hic");
    std::fs::write(&path, bytes).unwrap();
    let mut file = HicReader::open(FileSource::open(&path).unwrap()).unwrap();
    assert_eq!(file.query_entries(query).unwrap(), 1);
    let mut copy = file.try_clone_context().unwrap();
    assert_eq!(copy.cache_namespace(), file.cache_namespace());
    assert_fresh_context(&mut copy);
    file.write_entries(&mut actual).unwrap();
    assert_eq!(actual, expected);
    assert_eq!(copy.query_entries(query).unwrap(), 1);
    copy.write_entries(&mut actual).unwrap();
    assert_eq!(actual, expected);
}
