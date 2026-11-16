//! bgzip FASTA and `.gzi` against htslib, in both directions: files htslib
//! writes are read by the crate, and files the crate writes are read and
//! indexed by htslib.

use crate::support::{Bgzipped, Rng, build_fai, htslib, random_fasta};
use omicsweave_core::{Locus, Region, SequenceRead};
use omicsweave_fasta::{
    BLOCK_TEXT_LEN, BgzfWriter, BlockHeader, Error, FastaReader, Gzi, Inflater, MAX_BLOCK_LEN,
    compress_fasta,
};
use proptest::prelude::*;

/// Random valid FASTA, made larger than a few BGZF blocks.
fn big_fasta(rng: &mut Rng) -> Option<Vec<u8>> {
    let mut text = Vec::new();
    while text.len() < 3 * BLOCK_TEXT_LEN {
        let mut part = random_fasta(rng, false);
        if !part.ends_with(b"\n") {
            part.push(b'\n');
        }
        // Make names unique across the concatenated parts.
        let tag = format!(">p{}_", text.len());
        let mut renamed = Vec::with_capacity(part.len());
        for line in part.split_inclusive(|&c| c == b'\n') {
            if line.starts_with(b">") {
                renamed.extend_from_slice(tag.as_bytes());
                renamed.extend_from_slice(&line[1..]);
            } else {
                renamed.extend_from_slice(line);
            }
        }
        text.extend_from_slice(&renamed);
    }
    build_fai(&text, rng).ok().map(|_| text)
}

fn queries(rng: &mut Rng, len: u32) -> Vec<(u32, u32)> {
    let mut out = vec![(0, len)];
    for _ in 0..15 {
        let start = rng.coord(0, len - 1);
        out.push((start, rng.coord(start + 1, len)));
    }
    out
}

proptest! {
    #![proptest_config(ProptestConfig { cases: 12, failure_persistence: Some(Box::new(proptest::test_runner::FileFailurePersistence::Direct(concat!(env!("CARGO_MANIFEST_DIR"), "/tests/cases/bgzf.proptest-regressions")))), ..ProptestConfig::default() })]

    #[test]
    fn reads_bgzip_written_by_htslib(seed in any::<u64>()) {
        htslib::quiet();
        let mut rng = Rng::new(seed);
        let Some(text) = big_fasta(&mut rng) else { return Ok(()) };
        let fai = build_fai(&text, &mut rng).unwrap();
        let dir = tempfile::tempdir().unwrap();
        let gz = dir.path().join("t.fa.gz");
        htslib::bgzip(&text, &gz);
        let bytes = std::fs::read(&gz).unwrap();
        let written_gzi = Gzi::parse(&std::fs::read(dir.path().join("t.fa.gz.gzi")).unwrap()).unwrap();

        // htslib's index built by reading the file equals the one bgzip
        // wrote with it, and this crate builds the same from the file.
        prop_assert!(htslib::build_index(&gz));
        let read_gzi = Gzi::parse(&std::fs::read(dir.path().join("t.fa.gz.gzi")).unwrap()).unwrap();
        prop_assert_eq!(&read_gzi, &written_gzi);
        prop_assert_eq!(&Gzi::build(&bytes[..]).unwrap(), &written_gzi);
        prop_assert_eq!(std::fs::read_to_string(dir.path().join("t.fa.gz.fai")).unwrap(), fai.to_text());
        let theirs = htslib::Faidx::load(&gz).unwrap();

        for gzi in [written_gzi, read_gzi] {
            let mut reader = FastaReader::new(&bytes[..], fai.clone(), Some(gzi.clone()));
            let mut ours = Bgzipped::new(&bytes, gzi);
            for entry in fai.entries().iter().filter(|e| e.len > 0) {
                for (start, end) in queries(&mut rng, entry.len) {
                    let expected = theirs.fetch(&entry.name, start, end).unwrap();
                    prop_assert_eq!(&ours.fetch_span(entry, start, end), &expected, "{}:{}-{}", entry.name, start, end);
                    prop_assert_eq!(&ours.fetch_blocks(entry, start, end), &expected, "{}:{}-{}", entry.name, start, end);
                    let mut out = vec![0u8; (end - start) as usize];
                    reader.sequence(&Locus::from(Region::new(&entry.name, start, end)), &mut out).unwrap();
                    prop_assert_eq!(&out, &expected);
                }
            }
        }
    }

    #[test]
    fn writes_bgzip_htslib_reads(seed in any::<u64>()) {
        htslib::quiet();
        let mut rng = Rng::new(seed);
        let Some(text) = big_fasta(&mut rng) else { return Ok(()) };
        let level = rng.range(0, 12) as i32;
        let (bytes, fai, gzi) = compress_fasta(&text[..], Vec::new(), level).unwrap();
        prop_assert_eq!(&fai, &build_fai(&text, &mut rng).unwrap());

        let dir = tempfile::tempdir().unwrap();
        let gz = dir.path().join("t.fa.gz");
        std::fs::write(&gz, &bytes).unwrap();
        prop_assert_eq!(htslib::bgzf_read_all(&gz), text.clone());

        // Same block boundaries in the text as bgzip.
        let reference = dir.path().join("r.fa.gz");
        htslib::bgzip(&text, &reference);
        let bgzip_gzi = Gzi::parse(&std::fs::read(dir.path().join("r.fa.gz.gzi")).unwrap()).unwrap();
        let text_starts = |g: &Gzi| g.starts().iter().map(|s| s.1).collect::<Vec<_>>();
        prop_assert_eq!(text_starts(&gzi), text_starts(&bgzip_gzi));

        // htslib indexes the file to the same `.fai` and `.gzi`.
        prop_assert!(htslib::build_index(&gz));
        prop_assert_eq!(std::fs::read_to_string(dir.path().join("t.fa.gz.fai")).unwrap(), fai.to_text());
        prop_assert_eq!(std::fs::read(dir.path().join("t.fa.gz.gzi")).unwrap(), gzi.to_bytes());

        // Our `.gzi` bytes are what bgzip writes: htslib serves reads with it.
        std::fs::write(dir.path().join("t.fa.gz.gzi"), gzi.to_bytes()).unwrap();
        let theirs = htslib::Faidx::load(&gz).unwrap();
        let mut ours = Bgzipped::new(&bytes, gzi);
        for entry in fai.entries().iter().filter(|e| e.len > 0) {
            for (start, end) in queries(&mut rng, entry.len) {
                let expected = theirs.fetch(&entry.name, start, end).unwrap();
                prop_assert_eq!(&ours.fetch_span(entry, start, end), &expected, "{}:{}-{}", entry.name, start, end);
            }
        }
    }
}

#[test]
fn bgzf_detects_corruption() {
    let mut w = BgzfWriter::new(Vec::new());
    std::io::Write::write_all(&mut w, b"ACGTACGTACGT").unwrap();
    let (mut bytes, _) = w.finish().unwrap();
    let mut out = vec![0u8; MAX_BLOCK_LEN];
    let mut inflater = Inflater::new();
    assert_eq!(inflater.inflate_block(&bytes, &mut out).unwrap().1, 12);

    let header = BlockHeader::parse(&bytes).unwrap();
    bytes[header.block_len - 8] ^= 1;
    assert!(matches!(
        inflater.inflate_block(&bytes, &mut out),
        Err(Error::Inflate(_))
    ));
    assert!(matches!(
        inflater.inflate_block(&bytes[..10], &mut out),
        Err(Error::Truncated { .. })
    ));
    assert!(BlockHeader::parse(b"not a gzip block").is_err());
}

#[test]
fn gzi_rejects_unordered_and_short_input() {
    let mut b = 2u64.to_le_bytes().to_vec();
    for v in [100u64, 1000, 50, 2000] {
        b.extend_from_slice(&v.to_le_bytes());
    }
    assert!(matches!(Gzi::parse(&b), Err(Error::Invalid { .. })));
    assert!(matches!(Gzi::parse(&b[..20]), Err(Error::Truncated { .. })));
}

#[test]
fn gzi_bounds_each_block_before_reads_are_planned() {
    let encode = |points: &[(u64, u64)]| {
        let mut bytes = (points.len() as u64).to_le_bytes().to_vec();
        for &(compressed, text) in points {
            bytes.extend(compressed.to_le_bytes());
            bytes.extend(text.to_le_bytes());
        }
        bytes
    };
    let limit = MAX_BLOCK_LEN as u64;
    for points in [vec![(1, limit + 1)], vec![(28, 0), (29, limit + 1)]] {
        assert!(matches!(
            Gzi::parse(&encode(&points)),
            Err(Error::Invalid { .. })
        ));
    }
    let points = [(28, 0), (29 + limit, limit), (57 + limit, limit)];
    let index = Gzi::parse(&encode(&points)).unwrap();
    assert_eq!(&index.starts()[1..], &points);
    assert_eq!(Gzi::parse(&index.to_bytes()).unwrap(), index);
}
