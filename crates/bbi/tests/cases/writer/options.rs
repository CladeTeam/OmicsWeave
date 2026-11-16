use crate::support::libbigwig::Reader;
use crate::support::sidecar_check::check_sidecar;
use omicsweave_bbi::{BbiReader, BedSchema, BigBedWriter, BigWigWriter, Chrom, Kind, WriteOptions};
use omicsweave_core::{FileSource, Locus, Region, TrackRead};
use std::fs::File;
use std::io::{BufWriter, Cursor};

#[test]
fn writers_reject_unsorted_and_invalid_items() {
    let chroms = vec![
        Chrom {
            name: "a".into(),
            len: 1000,
        },
        Chrom {
            name: "b".into(),
            len: 1000,
        },
    ];
    let mut w = BigWigWriter::create(
        std::io::Cursor::new(Vec::new()),
        chroms.clone(),
        WriteOptions::default(),
    )
    .unwrap();
    w.add(1, 10, 20, 1.0).unwrap();
    assert!(w.add(1, 15, 30, 1.0).is_err(), "overlap");
    assert!(w.add(0, 0, 5, 1.0).is_err(), "earlier chromosome");
    assert!(w.add(1, 30, 30, 1.0).is_err(), "empty item");
    assert!(w.add(1, 30, 1001, 1.0).is_err(), "past the chromosome");
    assert!(w.add(1, 30, 40, f32::NAN).is_err(), "NaN");
    assert!(w.add(2, 0, 1, 1.0).is_err(), "unknown chromosome");

    let mut b = BigBedWriter::create(
        std::io::Cursor::new(Vec::new()),
        chroms.clone(),
        BedSchema::bed(3),
        WriteOptions::default(),
    )
    .unwrap();
    b.add(0, 10, 20, b"").unwrap();
    b.add(0, 10, 15, b"").unwrap();
    assert!(b.add(0, 5, 15, b"").is_err(), "unsorted");
    assert!(b.add(0, 30, 40, b"x\0y").is_err(), "NUL");

    let dup = vec![chroms[0].clone(), chroms[0].clone()];
    assert!(
        BigWigWriter::create(
            std::io::Cursor::new(Vec::new()),
            dup,
            WriteOptions::default()
        )
        .is_err()
    );
}

#[test]
fn empty_files_are_readable() {
    let dir = tempfile::tempdir().unwrap();
    for kind in [Kind::BigWig, Kind::BigBed] {
        let path = dir.path().join(format!("{kind:?}"));
        let out = BufWriter::new(File::create(&path).unwrap());
        let chroms = vec![Chrom {
            name: "chr1".into(),
            len: 10,
        }];
        let (out, sidecar) = match kind {
            Kind::BigWig => BigWigWriter::create(out, chroms, WriteOptions::default())
                .unwrap()
                .finish()
                .unwrap(),
            Kind::BigBed => {
                BigBedWriter::create(out, chroms, BedSchema::bed(3), WriteOptions::default())
                    .unwrap()
                    .finish()
                    .unwrap()
            }
        };
        drop(out.into_inner().unwrap());
        assert_eq!(sidecar.levels.len(), 1);
        assert!(sidecar.levels[0].blocks.is_empty());
        assert!(Reader::open(&path).is_some());
        let mut reader = BbiReader::open(FileSource::open(&path).unwrap()).unwrap();
        if kind == Kind::BigWig {
            let mut v = [0.0; 10];
            reader
                .values(&Locus::from(Region::new("chr1", 0, 10)), &mut v)
                .unwrap();
            assert!(v.iter().all(|x| x.is_nan()));
        }
        check_sidecar(&path, &sidecar);
    }
}

#[test]
fn writer_tree_fanout_fits_the_encoded_node_count() {
    for block_size in [0, 1, u16::MAX as u32 + 1, u32::MAX] {
        let mut output = Cursor::new(Vec::new());
        assert!(
            BigWigWriter::create(
                &mut output,
                vec![Chrom {
                    name: "c".into(),
                    len: 10
                }],
                WriteOptions {
                    block_size,
                    ..WriteOptions::default()
                },
            )
            .is_err()
        );
        assert!(output.into_inner().is_empty());
    }
    let count = u16::MAX as u32;
    let mut writer = BigWigWriter::create(
        Cursor::new(Vec::new()),
        vec![Chrom {
            name: "c".into(),
            len: count + 1,
        }],
        WriteOptions {
            items_per_slot: Some(1),
            block_size: count,
            compress: false,
            zoom_reductions: Some(vec![]),
            ..WriteOptions::default()
        },
    )
    .unwrap();
    for i in 0..count {
        writer.add(0, i, i + 1, i as f32).unwrap();
    }
    let (output, sidecar) = writer.finish().unwrap();
    assert_eq!(sidecar.levels[0].blocks.len(), count as usize);
    let mut reader = BbiReader::open(output.into_inner()).unwrap();
    let mut out = [0.0; 2];
    reader
        .values(&Locus::from(Region::new("c", count - 2, count)), &mut out)
        .unwrap();
    assert_eq!(out, [(count - 2) as f32, (count - 1) as f32]);
}
