//! Physical-member windows, independent metadata and true htslib results.

use std::io::Write;

use omicsweave_core::{Locus, Region, SequenceRead};
use omicsweave_fasta::{BgzfWriter, EOF_BLOCK, Fai, FaiEntry, FastaReader, Gzi, MAX_BLOCK_LEN};

use crate::support::{Source, htslib};

fn fixture(empty_members: usize, bases: &[u8]) -> (Vec<u8>, Fai, Gzi) {
    let mut bytes = Vec::new();
    for _ in 0..empty_members {
        bytes.extend_from_slice(&EOF_BLOCK);
    }
    let mut writer = BgzfWriter::new(bytes);
    writer.write_all(b">s\n").unwrap();
    writer.write_all(bases).unwrap();
    writer.write_all(b"\n").unwrap();
    let (bytes, _) = writer.finish().unwrap();
    let fai = Fai::from_entries([FaiEntry {
        name: "s".into(),
        len: bases.len() as u32,
        offset: 3,
        line_bases: bases.len() as u32,
        line_width: bases.len() as u32 + 1,
    }]);
    // htslib leaves the first non-empty member implicit, even after empties.
    let gzi = Gzi::parse(&0u64.to_le_bytes()).unwrap();
    (bytes, fai, gzi)
}

#[test]
fn leading_empty_member_windows_stay_bounded_and_match_htslib() {
    let long: Vec<u8> = (0u32..4096)
        .map(|i| {
            let v = i.wrapping_mul(0x9E37_79B9).rotate_left(13) ^ i.wrapping_mul(31);
            b"ACGT"[(v as usize >> 16) & 3]
        })
        .collect();
    let dir = tempfile::tempdir().unwrap();
    for count in [1, 2, 17, 2338, 2339, 2340, 2341, 2342, 3000, 5000] {
        for bases in [b"ACGTacgtNNNN".as_slice(), long.as_slice()] {
            let (bytes, fai, gzi) = fixture(count, bases);
            let physical = bytes.len();
            let path = dir
                .path()
                .join(format!("leading-{count}-{}.fa.gz", bases.len()));
            std::fs::write(&path, &bytes).unwrap();
            assert!(htslib::build_index(&path));
            assert_eq!(
                std::fs::read(omicsweave_fasta::fai_path(&path)).unwrap(),
                fai.to_bytes()
            );
            assert_eq!(
                std::fs::read(omicsweave_fasta::gzi_path(&path)).unwrap(),
                gzi.to_bytes()
            );
            let oracle = htslib::Faidx::load(&path).unwrap();
            let source = Source::new(bytes);
            let mut reader = FastaReader::new(&source, fai, Some(gzi));
            let len = bases.len() as u32;
            let mut windows = vec![(0, len), (1, 2)];
            if len == 4096 {
                windows.extend([(27, 100), (4000, 4096)]);
            }
            for (start, end) in windows {
                source.reset();
                let n = (end - start) as usize;
                let mut guarded = vec![0xA5; n + 2];
                reader
                    .sequence(
                        &Locus::from(Region::new("s", start, end)),
                        &mut guarded[1..n + 1],
                    )
                    .unwrap();
                assert_eq!(
                    &guarded[1..n + 1],
                    &bases[start as usize..end as usize],
                    "count={count}, range={start}..{end}"
                );
                assert_eq!(oracle.fetch("s", start, end).unwrap(), guarded[1..n + 1]);
                assert_eq!((guarded[0], guarded[n + 1]), (0xA5, 0xA5));
                assert!(source.largest.get() <= MAX_BLOCK_LEN, "count={count}");
                assert!(
                    source.calls.get() <= 3 + physical / (MAX_BLOCK_LEN - EOF_BLOCK.len()),
                    "count={count}"
                );
                assert!(
                    source.total.get() <= 2 * physical + MAX_BLOCK_LEN,
                    "count={count}, read={}",
                    source.total.get()
                );
            }
        }
    }
}

#[test]
fn failed_empty_crc_and_truncated_data_allow_repaired_source_queries() {
    let bases = b"ACGTacgtNNNN";
    let count = 2342;
    let (bytes, fai, gzi) = fixture(count, bases);
    let good = bytes.clone();
    let source = Source::new(bytes);
    let mut reader = FastaReader::new(&source, fai, Some(gzi));
    let mut out = [0xA5; 14];
    let crc = 2340 * EOF_BLOCK.len() + EOF_BLOCK.len() - 8;
    source.bytes.borrow_mut()[crc] ^= 1;
    assert!(
        reader
            .sequence(&Locus::from(Region::new("s", 0, 12)), &mut out[1..13])
            .is_err()
    );
    assert_eq!(out, [0xA5; 14]);
    source.bytes.borrow_mut()[crc] ^= 1;
    reader
        .sequence(&Locus::from(Region::new("s", 0, 12)), &mut out[1..13])
        .unwrap();
    assert_eq!(&out, b"\xA5ACGTacgtNNNN\xA5");

    source
        .bytes
        .borrow_mut()
        .truncate(count * EOF_BLOCK.len() + 20);
    assert!(
        reader
            .sequence(&Locus::from(Region::new("s", 0, 12)), &mut out[1..13])
            .is_err()
    );
    *source.bytes.borrow_mut() = good;
    reader
        .sequence(&Locus::from(Region::new("s", 0, 12)), &mut out[1..13])
        .unwrap();
    assert_eq!(&out, b"\xA5ACGTacgtNNNN\xA5");
}

#[test]
fn reader_bounds_reads_when_htslib_indexes_skip_empty_members() {
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
    let source = Source::new(bytes);
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
