use omicsweave_core::Strand;

use crate::support::*;

#[test]
fn projection_species_insertions_overlap_and_missing_rows() {
    let records = vec![
        (2, 6, "a;s ref.c 2 4 + 100 a-CGn;s dog.x 0 4 - 99 TtG-R;s dog.y 0 5 + 99 AAAAA;e missing.c 0 4 + 100 I;".into()),
        (4, 7, "a;s ref.c 4 3 + 100 TGA;s cat.x 0 3 + 99 CYN;".into()),
    ];
    let mut r = reader(&records, None, b'.');
    let bases = fetch(&mut r, 0, 10, Strand::Forward);
    assert_eq!(r.species(), ["ref", "dog", "cat"]);
    assert_eq!(&bases[..10], b"..aCTGA...");
    assert_eq!(&bases[10..20], b"..TG-R....");
    assert_eq!(&bases[20..], b"....CYN...");
    let mut fixed = reader(
        &records,
        Some(vec!["cat".into(), "ref".into(), "missing".into()]),
        b'.',
    );
    let bases = fetch(&mut fixed, 0, 10, Strand::Reverse);
    assert_eq!(&bases[..10], b"...NRG....");
    assert_eq!(&bases[10..20], b"...TCAGt..");
    assert_eq!(&bases[20..], b"..........");
}

#[test]
fn fill_stays_literal_while_covered_iupac_bases_complement() {
    let records = [(
        2,
        7,
        "s ref.c 2 5 + 100 AUryN;s dog.c 0 4 + 99 -kmBV;".into(),
    )];
    let mut r = reader(&records, Some(vec!["ref".into(), "absent".into()]), b'A');
    let got = fetch(&mut r, 0, 10, Strand::Reverse);
    assert_eq!(&got[..10], b"AAANryATAA");
    assert_eq!(&got[10..], b"AAAAAAAAAA");
}

#[test]
fn utf8_species_and_strand_coordinates_remain_independent() {
    let mut r = reader(
        &[(2, 5, "s ref.c 2 3 + 100 ACG;s 狼.x 5 3 - 90 TGn;".into())],
        None,
        b'.',
    );
    assert_eq!(fetch(&mut r, 2, 5, Strand::Forward), b"ACGTGn");
    assert_eq!(r.species(), ["ref", "狼"]);
    assert_eq!(fetch(&mut r, 2, 5, Strand::Reverse), b"CGTnCA");
}

#[test]
fn all_iupac_codes_reverse_with_case_and_literal_gap_preserved() {
    let mut r = reader(
        &[(0, 32, "s ref.c 0 32 + 100 ACGTURYSWKMBDHVNacgturyswkmbdhvn;s other.x 0 31 - 100 ACGTURYSWKMBDHVNacgturyswkmbdhv-;".into())],
        Some(vec!["ref".into(), "other".into()]),
        b'.',
    );
    assert_eq!(
        fetch(&mut r, 0, 32, Strand::Reverse),
        b"nbdhvkmwsryaacgtNBDHVKMWSRYAACGT-bdhvkmwsryaacgtNBDHVKMWSRYAACGT"
    );
}

#[test]
fn clipped_query_discards_reference_insertions_at_both_edges() {
    let mut r = reader(
        &[(
            5,
            9,
            "a;s ref.c 5 4 + 100 --A-CG-T--;s dog.x 7 8 - 100 ttR-yN-Agg;".into(),
        )],
        None,
        b'.',
    );
    assert_eq!(fetch(&mut r, 6, 8, Strand::Forward), b"CGyN");
    assert_eq!(r.species(), ["ref", "dog"]);
}

#[test]
fn explicit_empty_species_selects_zero_rows_without_disabling_validation() {
    let mut r = reader(
        &[(0, 2, "s ref.c 0 2 + 100 AC;".into())],
        Some(vec![]),
        b'.',
    );
    assert!(fetch(&mut r, 0, 4, Strand::Forward).is_empty());
    assert!(r.species().is_empty());
    let mut bad = reader(
        &[(0, 2, "s ref.c 0 2 + 100 AC;s ignored.x 0 2 + 10 A-;".into())],
        Some(vec![]),
        b'.',
    );
    assert!(bad.query(omicsweave_core::Region::new("c", 0, 2)).is_err());
}
