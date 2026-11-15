//! bgzip FASTA and `.gzi` against htslib, in both directions: files htslib
//! writes are read by the crate, and files the crate writes are read and
//! indexed by htslib.

mod common;

use common::{Bgzipped, Rng, build_fai, htslib, random_fasta};
use omicsweave_core::{Locus, Region, SequenceRead};
use omicsweave_fasta::{BLOCK_TEXT_LEN, FastaReader, Gzi, compress_fasta};
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

#[test]
fn reader_carries_split_crlf_between_bgzf_blocks_and_recovers_after_error() {
    use omicsweave_fasta::{BgzfWriter, FaiBuilder};
    use std::io::Write;

    let text = b">s\r\nACGT\r\nacgt\r\nNNNN\r\n";
    let mut builder = FaiBuilder::new();
    builder.push(text).unwrap();
    let fai = builder.finish().unwrap();
    // Every position becomes a physical BGZF boundary, including both
    // halves of CRLF and the first/last base of a sequence line.
    for split in 1..text.len() {
        let mut writer = BgzfWriter::new(Vec::new());
        writer.write_all(&text[..split]).unwrap();
        writer.flush().unwrap();
        writer.write_all(&text[split..]).unwrap();
        let (bytes, gzi) = writer.finish().unwrap();
        let mut reader = FastaReader::new(bytes.as_slice(), fai.clone(), Some(gzi));
        let mut out = [0xAA; 14];
        reader
            .sequence(&Locus::from(Region::new("s", 0, 12)), &mut out[1..13])
            .unwrap();
        assert_eq!(&out, b"\xAAACGTacgtNNNN\xAA");
        assert!(
            reader
                .sequence(&Locus::from(Region::new("s", 0, 13)), &mut out[1..13])
                .is_err()
        );
        reader
            .sequence(&Locus::from(Region::new("s", 3, 5)), &mut out[1..3])
            .unwrap();
        assert_eq!(&out[1..3], b"Ta");
        reader
            .sequence(&Locus::from(Region::new("s", 0, 12)), &mut out[1..13])
            .unwrap();
        assert_eq!(&out, b"\xAAACGTacgtNNNN\xAA");
    }
}

#[test]
fn reader_bounds_reads_when_htslib_indexes_skip_empty_members() {
    use omicsweave_core::{ReadAt, Result};
    use omicsweave_fasta::{BgzfWriter, EOF_BLOCK, MAX_BLOCK_LEN};
    use std::{cell::Cell, io::Write};

    struct Counted<'a> {
        bytes: &'a [u8],
        calls: Cell<usize>,
        total: Cell<usize>,
        largest: Cell<usize>,
    }
    impl ReadAt for Counted<'_> {
        fn len(&self) -> u64 {
            self.bytes.len() as u64
        }
        fn read_exact_at(&self, out: &mut [u8], offset: u64) -> Result<()> {
            self.calls.set(self.calls.get() + 1);
            self.total.set(self.total.get() + out.len());
            self.largest.set(self.largest.get().max(out.len()));
            self.bytes.read_exact_at(out, offset)
        }
    }

    let mut writer = BgzfWriter::new(Vec::new());
    writer.write_all(b">s\nACGT\n").unwrap();
    let (mut bytes, _) = writer.finish().unwrap();
    for _ in 0..3000 {
        bytes.extend(EOF_BLOCK);
    }
    let mut writer = BgzfWriter::new(bytes);
    writer.write_all(b"acgt\nNNNN\n").unwrap();
    let (bytes, _) = writer.finish().unwrap();
    let dir = tempfile::tempdir().unwrap();
    let path = dir.path().join("empty.fa.gz");
    std::fs::write(&path, &bytes).unwrap();
    htslib::quiet();
    assert!(htslib::build_index(&path));
    let fai =
        omicsweave_fasta::Fai::parse(&std::fs::read(dir.path().join("empty.fa.gz.fai")).unwrap())
            .unwrap();
    let gzi = Gzi::parse(&std::fs::read(dir.path().join("empty.fa.gz.gzi")).unwrap()).unwrap();
    assert!(gzi.starts()[1].0 > MAX_BLOCK_LEN as u64);
    let reference = htslib::Faidx::load(&path).unwrap();
    let expected = reference.fetch("s", 0, 12).unwrap();
    assert_eq!(&expected, b"ACGTacgtNNNN");
    let source = Counted {
        bytes: &bytes,
        calls: Cell::new(0),
        total: Cell::new(0),
        largest: Cell::new(0),
    };
    let mut reader = FastaReader::new(&source, fai, Some(gzi));
    let mut out = [0; 12];
    reader
        .sequence(&Locus::from(Region::new("s", 0, 12)), &mut out)
        .unwrap();
    assert_eq!(&out, expected.as_slice());
    assert_eq!(source.calls.get(), 2);
    assert!(source.largest.get() <= MAX_BLOCK_LEN);
    assert!(source.total.get() <= 2 * MAX_BLOCK_LEN);
}

proptest! {
    #![proptest_config(ProptestConfig { cases: 12, ..ProptestConfig::default() })]

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
