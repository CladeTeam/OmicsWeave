//! `.fai` building and plain FASTA reads against htslib (`samtools faidx`).

mod common;

use common::{Rng, build_fai, fetch_plain, htslib, random_fasta};
use omicsweave_core::{Locus, Region, SequenceRead};
use omicsweave_fasta::{Fai, FastaReader};
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
    #![proptest_config(ProptestConfig { cases: 200, ..ProptestConfig::default() })]

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
