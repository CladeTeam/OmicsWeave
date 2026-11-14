//! Index building from existing files (`Fai::build`, `Gzi::build`,
//! `build_index`) and while writing (`FastaWriter`), against htslib's
//! `fai_build3` (`samtools faidx`) and `bgzip -r`.

mod common;

use std::path::Path;

use common::{Rng, build_fai, htslib, random_fasta};
use omicsweave_core::{Region, SequenceRead};
use omicsweave_fasta::{
    EOF_BLOCK, Error, Fai, FastaReader, FastaWriter, Gzi, build_index, build_index_files, fai_path,
    gzi_path,
};
use proptest::prelude::*;

/// What htslib's `fai_build3` writes for the file at `path`: the `.fai`
/// text and, for compressed input, the `.gzi` bytes; `None` when it fails.
fn htslib_index(path: &Path) -> Option<(String, Option<Vec<u8>>)> {
    let _ = std::fs::remove_file(fai_path(path));
    let _ = std::fs::remove_file(gzi_path(path));
    htslib::build_index(path).then(|| {
        (
            std::fs::read_to_string(fai_path(path)).unwrap(),
            std::fs::read(gzi_path(path)).ok(),
        )
    })
}

fn starts_compressed(bytes: &[u8]) -> bool {
    bytes.len() >= 18 && bytes[..2] == [0x1f, 0x8b]
}

/// Random records for the writer: unique names (some non-ASCII), bases from
/// the printable range, lengths from one base to several BGZF blocks.
fn random_records(rng: &mut Rng) -> Vec<(String, Vec<u8>)> {
    const BASES: &[u8] = b"ACGTACGTACGTNacgtnRYKM*-.";
    (0..rng.range(1, 8))
        .map(|i| {
            let name = if rng.chance(15) {
                format!("染色体{i}|x")
            } else {
                format!("r{i}_{}", rng.next() % 1000)
            };
            let len = match rng.range(0, 9) {
                0 => rng.range(1, 3) as usize,
                1 => rng.range(60_000, 200_000) as usize,
                _ => rng.range(1, 5_000) as usize,
            };
            (name, (0..len).map(|_| *rng.pick(BASES)).collect())
        })
        .collect()
}

fn write_all<W: std::io::Write>(
    mut writer: FastaWriter<W>,
    records: &[(String, Vec<u8>)],
) -> (W, Fai, Option<Gzi>) {
    for (name, seq) in records {
        writer.add(name, seq).unwrap();
    }
    writer.finish().unwrap()
}

/// Every record reads back whole and in random slices.
fn check_reads(
    reader: &mut impl SequenceRead,
    records: &[(String, Vec<u8>)],
    rng: &mut Rng,
) -> std::result::Result<(), TestCaseError> {
    for (name, seq) in records {
        let len = seq.len() as u32;
        let mut ranges = vec![(0, len)];
        for _ in 0..5 {
            let start = rng.coord(0, len - 1);
            ranges.push((start, rng.coord(start + 1, len)));
        }
        for (start, end) in ranges {
            let mut out = vec![0u8; (end - start) as usize];
            reader
                .sequence(Region::new(name, start, end), &mut out)
                .unwrap();
            prop_assert_eq!(&out[..], &seq[start as usize..end as usize]);
        }
    }
    Ok(())
}

/// Random text split into parts bgzipped separately and concatenated, so
/// the file has short blocks, empty members between data blocks and, when
/// the first part is empty, before the first.
fn bgzip_in_parts(text: &[u8], rng: &mut Rng, dir: &Path) -> Vec<u8> {
    let mut cuts: Vec<usize> = (0..rng.range(0, 4))
        .map(|_| rng.range(0, text.len() as u64) as usize)
        .collect();
    if rng.chance(30) {
        cuts.push(0);
    }
    cuts.push(text.len());
    cuts.sort();
    let mut parts = Vec::new();
    let mut at = 0;
    for cut in cuts {
        parts.push(&text[at..cut]);
        at = cut;
    }
    htslib::bgzip_parts(&parts, dir)
}

#[test]
fn files_beside_the_fasta_match_samtools() {
    htslib::quiet();
    let dir = tempfile::tempdir().unwrap();
    let mut rng = Rng::new(7);
    let records = random_records(&mut rng);
    for (file, compressed) in [("a.fa", false), ("a.fa.gz", true)] {
        let ours = dir.path().join(file);
        let writer = if compressed {
            FastaWriter::bgzf(Vec::new(), 60, 6).unwrap()
        } else {
            FastaWriter::new(Vec::new(), 60)
        };
        let (bytes, fai, gzi) = write_all(writer, &records);
        std::fs::write(&ours, &bytes).unwrap();
        let theirs = dir.path().join(format!("htslib.{file}"));
        std::fs::write(&theirs, &bytes).unwrap();
        let (fai_text, gzi_bytes) = htslib_index(&theirs).unwrap();

        assert_eq!(
            fai_path(&ours),
            dir.path().join(format!("{file}.fai")),
            "samtools naming"
        );
        assert_eq!(gzi_path(&ours), dir.path().join(format!("{file}.gzi")));
        assert_eq!(
            build_index_files(&ours).unwrap(),
            (fai.clone(), gzi.clone())
        );
        assert_eq!(std::fs::read_to_string(fai_path(&ours)).unwrap(), fai_text);
        assert_eq!(std::fs::read(gzi_path(&ours)).ok(), gzi_bytes);
        assert_eq!(Fai::load(&ours).unwrap(), fai);

        std::fs::remove_file(fai_path(&ours)).unwrap();
        assert_eq!(Fai::build_file(&ours).unwrap(), fai);
        assert_eq!(std::fs::read_to_string(fai_path(&ours)).unwrap(), fai_text);
        if let Some(gzi) = gzi {
            std::fs::remove_file(gzi_path(&ours)).unwrap();
            assert_eq!(Gzi::build_file(&ours).unwrap(), gzi);
            assert_eq!(std::fs::read(gzi_path(&ours)).unwrap(), gzi_bytes.unwrap());
            assert_eq!(Gzi::load(&ours).unwrap(), gzi);
        } else {
            assert!(Gzi::build_file(&ours).is_err());
        }
    }
}

#[test]
fn leading_empty_members_index_and_read_like_htslib() {
    htslib::quiet();
    let dir = tempfile::tempdir().unwrap();
    let text = b">s\nACGTACGT\nAC\n>t\nNNNN\n";
    let mut bytes = Vec::new();
    for _ in 0..3 {
        bytes.extend(EOF_BLOCK);
    }
    bytes.extend(htslib::bgzip_parts(&[text, b"", text], dir.path()));
    let path = dir.path().join("lead.fa.gz");
    std::fs::write(&path, &bytes).unwrap();
    let (fai_text, gzi_bytes) = htslib_index(&path).unwrap();
    let (fai, gzi) = build_index(&bytes[..]).unwrap();
    let gzi = gzi.unwrap();
    assert_eq!(fai.to_text(), fai_text);
    assert_eq!(gzi.to_bytes(), gzi_bytes.unwrap());
    assert_eq!(Gzi::build(&bytes[..]).unwrap(), gzi);
    assert_eq!(htslib::reindex(&path).unwrap(), gzi.to_bytes());

    let reference = htslib::Faidx::load(&path).unwrap();
    let mut reader = FastaReader::new(&bytes[..], fai, Some(gzi));
    for (name, start, end) in [("s", 0, 10), ("s", 3, 9), ("t", 0, 4), ("t", 1, 2)] {
        let mut out = vec![0u8; end as usize - start as usize];
        reader
            .sequence(Region::new(name, start, end), &mut out)
            .unwrap();
        assert_eq!(out, reference.fetch(name, start, end).unwrap());
    }
}

#[test]
fn plain_gzip_is_unsupported_like_htslib() {
    htslib::quiet();
    let dir = tempfile::tempdir().unwrap();
    let text = b">s\nACGTACGT\n";
    let mut compressor = libdeflater::Compressor::new(libdeflater::CompressionLvl::default());
    let mut gzip = vec![0u8; compressor.gzip_compress_bound(text.len())];
    let n = compressor.gzip_compress(text, &mut gzip).unwrap();
    gzip.truncate(n);
    let mut bgzf_then_gzip = htslib::bgzip_parts(&[text], dir.path());
    bgzf_then_gzip.extend(&gzip);

    let path = dir.path().join("plain.fa.gz");
    std::fs::write(&path, &gzip).unwrap();
    assert!(htslib_index(&path).is_none());
    assert!(htslib::reindex(&path).is_none());

    // A gzip member after BGZF data fails `bgzip -r`. `samtools faidx`
    // stops reading there without an error and indexes only the data
    // before it; index building fails instead of describing a prefix.
    let path = dir.path().join("mixed.fa.gz");
    std::fs::write(&path, &bgzf_then_gzip).unwrap();
    assert!(htslib::reindex(&path).is_none());

    for bytes in [&gzip, &bgzf_then_gzip] {
        assert!(build_index(&bytes[..]).is_err());
        assert!(Fai::build(&bytes[..]).is_err());
        assert!(Gzi::build(&bytes[..]).is_err());
    }
    assert!(matches!(Fai::build(&gzip[..]), Err(Error::Unsupported(_))));
    assert!(matches!(Gzi::build(&gzip[..]), Err(Error::Unsupported(_))));
    assert!(Gzi::build(&text[..]).is_err(), "plain text has no .gzi");
}

#[test]
fn writer_rejects_records_that_would_not_index_as_written() {
    let mut writer = FastaWriter::new(Vec::new(), 4);
    writer.add("ok", b"ACGTA").unwrap();
    let rejected: [(&str, &[u8]); 13] = [
        ("", b"ACGT"),
        ("a b", b"ACGT"),
        ("a\tb", b"ACGT"),
        ("a\rb", b"ACGT"),
        ("a\0b", b"ACGT"),
        ("ok", b"ACGT"),
        ("e", b""),
        ("g", b"AC GT"),
        ("g", b"AC\nGT"),
        ("g", b"AC>GT"),
        ("g", b">ACGT"),
        ("g", b"AC\x7fGT"),
        ("g", "ACéGT".as_bytes()),
    ];
    for (name, seq) in rejected {
        assert!(
            matches!(writer.add(name, seq), Err(Error::Invalid { .. })),
            "{name:?} {seq:?}"
        );
    }
    writer.add("染色体", b"N").unwrap();
    let (bytes, fai, gzi) = writer.finish().unwrap();
    assert_eq!(bytes, ">ok\nACGT\nA\n>染色体\nN\n".as_bytes());
    assert_eq!(fai, Fai::build(&bytes[..]).unwrap());
    assert_eq!(gzi, None);
    assert!(FastaWriter::new(Vec::new(), 60).finish().is_err());
    assert!(FastaWriter::bgzf(Vec::new(), 60, 13).is_err());
}

proptest! {
    #![proptest_config(ProptestConfig { cases: 48, ..ProptestConfig::default() })]

    /// The writer's `.fai`/`.gzi` equal the indexes built from its output
    /// by this crate and by htslib, byte for byte, and the output reads
    /// back as written.
    #[test]
    fn writer_indexes_equal_htslib_indexes_of_its_output(seed in any::<u64>()) {
        htslib::quiet();
        let mut rng = Rng::new(seed);
        let records = random_records(&mut rng);
        let width = match rng.range(0, 5) {
            0 => 0,
            1 => records[0].1.len(),
            2 => *rng.pick(&[1, 2, 60, 61, 70, 80]),
            _ => rng.range(1, 200) as usize,
        };
        let level = rng.range(0, 12) as i32;
        let dir = tempfile::tempdir().unwrap();

        let (plain, fai, none) = write_all(FastaWriter::new(Vec::new(), width), &records);
        prop_assert!(none.is_none());
        prop_assert_eq!(&Fai::build(&plain[..]).unwrap(), &fai);
        prop_assert_eq!(&build_fai(&plain, &mut rng).unwrap(), &fai);
        prop_assert_eq!(fai.entries().len(), records.len());
        let path = dir.path().join("w.fa");
        std::fs::write(&path, &plain).unwrap();
        let (fai_text, gzi_bytes) = htslib_index(&path).unwrap();
        prop_assert_eq!(&fai.to_text(), &fai_text);
        prop_assert!(gzi_bytes.is_none());
        check_reads(&mut FastaReader::new(&plain[..], fai.clone(), None), &records, &mut rng)?;

        let (gz, gz_fai, gzi) =
            write_all(FastaWriter::bgzf(Vec::new(), width, level).unwrap(), &records);
        let gzi = gzi.unwrap();
        prop_assert_eq!(&gz_fai, &fai);
        prop_assert_eq!(&Gzi::build(&gz[..]).unwrap(), &gzi);
        prop_assert_eq!(build_index(&gz[..]).unwrap(), (fai.clone(), Some(gzi.clone())));
        let path = dir.path().join("w.fa.gz");
        std::fs::write(&path, &gz).unwrap();
        prop_assert_eq!(htslib::bgzf_read_all(&path), plain);
        prop_assert_eq!(htslib::reindex(&path).unwrap(), gzi.to_bytes());
        let (fai_text, gzi_bytes) = htslib_index(&path).unwrap();
        prop_assert_eq!(&fai.to_text(), &fai_text);
        prop_assert_eq!(gzi_bytes.unwrap(), gzi.to_bytes());
        check_reads(&mut FastaReader::new(&gz[..], fai, Some(gzi)), &records, &mut rng)?;
    }

    /// On random (possibly malformed) FASTA, plain and in bgzip files with
    /// short blocks and empty members, `build_index` succeeds exactly when
    /// htslib does, with the same `.fai` and `.gzi`; `Gzi::build` gives the
    /// `.gzi` of `bgzip -r` whatever the text.
    #[test]
    fn build_index_matches_htslib(seed in any::<u64>(), damage in any::<bool>()) {
        htslib::quiet();
        let mut rng = Rng::new(seed);
        let mut text = random_fasta(&mut rng, damage);
        if rng.chance(30) {
            let more = random_fasta(&mut rng, false);
            while text.len() < 200_000 {
                text.extend_from_slice(&more);
            }
        }
        let dir = tempfile::tempdir().unwrap();
        let plain = dir.path().join("t.fa");
        std::fs::write(&plain, &text).unwrap();
        let theirs = htslib_index(&plain);
        let ours = build_index(&text[..]);
        match (&ours, &theirs) {
            (Ok((fai, None)), Some((expected, None))) => prop_assert_eq!(&fai.to_text(), expected),
            (Err(_), None) => {}
            _ => prop_assert!(false, "plain: ours {:?}, htslib {:?}", ours, theirs),
        }
        prop_assert_eq!(Fai::build(&text[..]).ok(), ours.ok().map(|(fai, _)| fai));

        for variant in 0..3 {
            let gz = match variant {
                0 => {
                    let path = dir.path().join("whole.gz");
                    htslib::bgzip(&text, &path);
                    std::fs::read(&path).unwrap()
                }
                1 => bgzip_in_parts(&text, &mut rng, dir.path()),
                _ => {
                    let mut writer = omicsweave_fasta::BgzfWriter::new(Vec::new());
                    for chunk in text.chunks(rng.range(1, 70_000) as usize) {
                        std::io::Write::write_all(&mut writer, chunk).unwrap();
                        if rng.chance(50) {
                            std::io::Write::flush(&mut writer).unwrap();
                        }
                    }
                    writer.finish().unwrap().0
                }
            };
            let path = dir.path().join("t.fa.gz");
            std::fs::write(&path, &gz).unwrap();
            let theirs = htslib_index(&path);
            let ours = build_index(&gz[..]);
            let gzi = Gzi::build(&gz[..]).unwrap();
            prop_assert_eq!(&gzi.to_bytes(), &htslib::reindex(&path).unwrap(), "variant {}", variant);
            match (&ours, &theirs) {
                (Ok((fai, Some(built))), Some((expected, Some(expected_gzi)))) => {
                    prop_assert_eq!(&fai.to_text(), expected);
                    prop_assert_eq!(&built.to_bytes(), expected_gzi);
                    prop_assert_eq!(built, &gzi);
                    let reference = htslib::Faidx::load(&path).unwrap();
                    let mut reader = FastaReader::new(&gz[..], fai.clone(), Some(gzi));
                    for entry in fai.entries().iter().filter(|e| e.len > 0) {
                        let start = rng.coord(0, entry.len - 1);
                        for (start, end) in [(0, entry.len), (start, rng.coord(start + 1, entry.len))] {
                            let mut out = vec![0u8; (end - start) as usize];
                            reader.sequence(Region::new(&entry.name, start, end), &mut out).unwrap();
                            prop_assert_eq!(out, reference.fetch(&entry.name, start, end).unwrap());
                        }
                    }
                }
                (Err(_), None) => {}
                _ => prop_assert!(false, "variant {}: ours {:?}, htslib {:?}", variant, ours, theirs),
            }
        }
    }

    /// Truncated, extended and header-corrupted bgzip files: index building
    /// fails exactly when htslib's does.
    #[test]
    fn damaged_bgzf_fails_like_htslib(seed in any::<u64>()) {
        htslib::quiet();
        let mut rng = Rng::new(seed);
        let mut text = random_fasta(&mut rng, false);
        let more = text.clone();
        while text.len() < 150_000 && !more.is_empty() {
            text.extend_from_slice(&more);
        }
        let dir = tempfile::tempdir().unwrap();
        let mut gz = bgzip_in_parts(&text, &mut rng, dir.path());
        let headers: Vec<u64> = Gzi::build(&gz[..]).unwrap().starts().iter().map(|s| s.0).collect();
        match rng.range(0, 3) {
            0 => gz.truncate(rng.range(0, gz.len() as u64 - 1) as usize),
            1 => {
                let n = rng.range(1, 40) as usize;
                gz.extend((0..n).map(|_| rng.next() as u8));
            }
            2 => {
                let mut tail = htslib::bgzip_parts(&[b">extra\nAC\n"], dir.path());
                tail.truncate(rng.range(1, tail.len() as u64) as usize);
                gz.extend(tail);
            }
            _ => {
                // Bytes htslib checks: magic, method, XLEN, the BC subfield
                // and BSIZE.
                let block = *rng.pick(&headers) as usize;
                let at = block + *rng.pick(&[0, 1, 2, 10, 11, 12, 13, 14, 15, 16, 17]);
                gz[at] ^= 1 << rng.range(0, 7);
            }
        }
        let path = dir.path().join("d.fa.gz");
        std::fs::write(&path, &gz).unwrap();
        let theirs = htslib_index(&path);
        let ours = build_index(&gz[..]);
        match (&ours, &theirs) {
            (Ok((fai, gzi)), Some((expected, expected_gzi))) => {
                prop_assert_eq!(&fai.to_text(), expected);
                prop_assert_eq!(gzi.as_ref().map(Gzi::to_bytes), expected_gzi.clone());
            }
            (Err(_), None) => {}
            // A gzip member that is not a BGZF block after BGZF data: see
            // `plain_gzip_is_unsupported_like_htslib`. `bgzip -r` fails.
            (Err(Error::Unsupported(_)), Some(_)) if starts_compressed(&gz) => {
                prop_assert!(htslib::reindex(&path).is_none());
            }
            _ => prop_assert!(false, "ours {:?}, htslib {:?}", ours, theirs),
        }
        if starts_compressed(&gz) {
            let gzi = Gzi::build(&gz[..]);
            let reindexed = htslib::reindex(&path);
            prop_assert_eq!(gzi.as_ref().ok().map(Gzi::to_bytes), reindexed.clone(), "{:?}", gzi);
        } else {
            prop_assert!(Gzi::build(&gz[..]).is_err());
        }
    }
}
