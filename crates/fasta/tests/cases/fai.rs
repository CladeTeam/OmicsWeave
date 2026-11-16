//! `.fai` building and plain FASTA reads against htslib (`samtools faidx`).

use crate::support::{Rng, build_fai, fetch_plain, htslib, random_fasta};
use omicsweave_core::{Locus, Region, SequenceRead};
use omicsweave_fasta::{Error, Fai, FaiEntry, FastaReader};
use proptest::prelude::*;

#[test]
fn built_names_round_trip_without_utf8_replacement() {
    use omicsweave_fasta::FaiBuilder;
    for text in [b">\xff\nACGT\n".as_slice(), b">\xfe\nACGT\n>\xff\nAAAA\n"] {
        let mut builder = FaiBuilder::new();
        let result = builder.push(text).and_then(|()| builder.finish());
        assert!(result.is_err());
    }
    let mut builder = FaiBuilder::new();
    builder.push(">染色体\nACGT\n".as_bytes()).unwrap();
    let fai = builder.finish().unwrap();
    assert_eq!(Fai::parse(fai.to_text().as_bytes()).unwrap(), fai);
}

#[test]
fn builder_deduplicates_names_across_arbitrary_chunks_like_htslib() {
    use omicsweave_fasta::FaiBuilder;
    let text = ">\nAC\n>染色体 desc\nACGT\n>\nTTTT\n>染色体\nCC\n>last\nG\n".as_bytes();
    let dir = tempfile::tempdir().unwrap();
    let path = dir.path().join("names.fa");
    std::fs::write(&path, text).unwrap();
    htslib::quiet();
    assert!(htslib::build_index(&path));
    let expected = std::fs::read_to_string(dir.path().join("names.fa.fai")).unwrap();
    for split in 0..=text.len() {
        let mut builder = FaiBuilder::new();
        builder.push(&text[..split]).unwrap();
        builder.push(&text[split..]).unwrap();
        let index = builder.finish().unwrap();
        assert_eq!(index.to_text(), expected, "split {split}");
        assert_eq!(index.entries().len(), 3);
        assert_eq!(index.get("").unwrap().len, 2);
        assert_eq!(index.get("染色体").unwrap().len, 4);
    }
}

proptest! {
    #![proptest_config(ProptestConfig { cases: 200, failure_persistence: Some(Box::new(proptest::test_runner::FileFailurePersistence::Direct(concat!(env!("CARGO_MANIFEST_DIR"), "/tests/cases/fai.proptest-regressions")))), ..ProptestConfig::default() })]

    /// Same `.fai` as htslib, or rejected by both.
    #[test]
    fn builds_the_same_index(seed in any::<u64>(), damage in any::<bool>()) {
        htslib::quiet();
        let mut rng = Rng::new(seed);
        let text = random_fasta(&mut rng, damage);
        let dir = tempfile::tempdir().unwrap();
        let fa = dir.path().join("t.fa");
        std::fs::write(&fa, &text).unwrap();

        let theirs = htslib::build_index(&fa).then(|| std::fs::read_to_string(dir.path().join("t.fa.fai")).unwrap());
        let ours = build_fai(&text, &mut rng);
        match (&ours, &theirs) {
            (Ok(fai), Some(expected)) => {
                prop_assert_eq!(&fai.to_text(), expected);
                prop_assert_eq!(&Fai::parse(expected.as_bytes()).unwrap(), fai);
            }
            (Err(_), None) => {}
            _ => prop_assert!(false, "ours {:?}, htslib {}", ours.as_ref().map(|f| f.to_text()), theirs.is_some()),
        }
    }

    /// Every sequence range reads the same bases as htslib.
    #[test]
    fn reads_the_same_bases(seed in any::<u64>()) {
        htslib::quiet();
        let mut rng = Rng::new(seed);
        let text = random_fasta(&mut rng, false);
        let Ok(fai) = build_fai(&text, &mut rng) else { return Ok(()) };
        let dir = tempfile::tempdir().unwrap();
        let fa = dir.path().join("t.fa");
        std::fs::write(&fa, &text).unwrap();
        prop_assert!(htslib::build_index(&fa));
        let theirs = htslib::Faidx::load(&fa).unwrap();
        let names: Vec<String> = fai.entries().iter().map(|e| e.name.clone()).collect();
        prop_assert_eq!(theirs.names(), names);

        let mut reader = FastaReader::new(&text[..], fai.clone(), None);
        for entry in fai.entries().iter().filter(|e| e.len > 0) {
            prop_assert_eq!(reader.chrom_len(&entry.name).unwrap(), entry.len);
            let mut ranges = vec![(0, entry.len)];
            for _ in 0..20 {
                let start = rng.coord(0, entry.len - 1);
                ranges.push((start, rng.coord(start + 1, entry.len)));
            }
            for (start, end) in ranges {
                let expected = theirs.fetch(&entry.name, start, end).unwrap();
                prop_assert_eq!(&fetch_plain(entry, &text, start, end), &expected, "{}:{}-{}", entry.name, start, end);
                let mut out = vec![0u8; (end - start) as usize];
                reader.sequence(&Locus::from(Region::new(&entry.name, start, end)), &mut out).unwrap();
                prop_assert_eq!(&out, &expected);
            }
        }
    }
}

#[test]
fn fai_rejects_ranges_outside_the_sequence() {
    let e = entry();
    assert!(matches!(
        e.text_range(5, 11),
        Err(Error::RegionOutOfRange { .. })
    ));
    assert!(matches!(
        e.text_range(6, 5),
        Err(Error::InvalidRegion { .. })
    ));
    assert_eq!(e.text_range(4, 9).unwrap(), 8..14);
    assert_eq!(e.text_range(7, 7).unwrap(), 3..3);
}

#[test]
fn fai_decode_checks_the_line_layout() {
    let e = entry();
    let text = b">s\nACGT\nACGT\nAC\n";
    let mut out = [0u8; 10];
    e.decode(0, &text[3..], &mut out).unwrap();
    assert_eq!(&out, b"ACGTACGTAC");

    let shifted = b">s\nACGTA\nCGT\nAC\n";
    assert!(matches!(
        e.decode(0, &shifted[3..], &mut out),
        Err(Error::Invalid { .. })
    ));
    assert!(matches!(
        e.decode(0, &text[3..8], &mut out),
        Err(Error::Truncated { .. })
    ));
}

#[test]
fn fai_parse_rejects_bad_rows() {
    assert!(Fai::parse(b"s\t10\t3\t4\n").is_err());
    assert!(Fai::parse(b"s\t10\t3\t4\t5\t9\n").is_err());
    assert!(Fai::parse(b"s\tten\t3\t4\t5\n").is_err());
    assert!(Fai::parse(b"s\t10\t3\t5\t4\n").is_err());
    assert_eq!(
        Fai::parse(b"s\t10\t3\t4\t5\n").unwrap().get("s"),
        Some(&entry())
    );
}

fn entry() -> FaiEntry {
    FaiEntry {
        name: "s".into(),
        len: 10,
        offset: 3,
        line_bases: 4,
        line_width: 5,
    }
}

#[test]
fn parsed_duplicate_names_keep_first_rows_and_still_validate_every_row() {
    let text = b"s\t12\t3\t6\t7\ns\t99\t15\t9\t10\nt\t8\t42\t4\t5\n";
    let fai = Fai::parse(text).unwrap();
    assert_eq!(fai.len(), 2);
    assert_eq!(fai.get("s").unwrap().len, 12);
    assert_eq!(fai.to_bytes(), b"s\t12\t3\t6\t7\nt\t8\t42\t4\t5\n");
    for tail in [
        b"s\tnan\t15\t9\t10\n".as_slice(),
        b"s\t99\t15\t10\t9\n",
        b"s\t99\t15\t9\t10\textra\n",
    ] {
        let mut invalid = text.to_vec();
        invalid.extend_from_slice(tail);
        assert!(Fai::parse(&invalid).is_err());
    }
}

#[test]
fn public_layout_overflow_and_short_input_leave_output_guards_intact() {
    let invalid = [
        FaiEntry {
            offset: u64::MAX,
            ..entry()
        },
        FaiEntry {
            line_bases: 0,
            ..entry()
        },
        FaiEntry {
            line_width: 3,
            ..entry()
        },
    ];
    for metadata in invalid {
        let mut out = [0xA5; 6];
        assert!(metadata.text_range(0, 4).is_err(), "{metadata:?}");
        assert!(
            metadata.decode(0, b"ACGT", &mut out[1..5]).is_err(),
            "{metadata:?}"
        );
        assert_eq!(out, [0xA5; 6]);
    }
    let mut out = [0xA5; 6];
    assert!(entry().decode(0, b"ACG", &mut out[1..5]).is_err());
    assert_eq!(out, [0xA5; 6]);
}
