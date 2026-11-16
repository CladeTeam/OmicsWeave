use omicsweave_core::{ContactQuery, ContactRead, ContactUnit, Error, Region};
use omicsweave_hic::{Chrom, HicHeader, HicReader, HicWriter, Record};
use std::io::Cursor;

#[test]
fn reader_checks_zero_resolution_and_writer_checks_duplicates() {
    let chrom = Chrom {
        name: "s".into(),
        len: 100,
    };
    let duplicate = HicHeader::new("g", vec![chrom.clone(), chrom.clone()], vec![1]);
    assert!(HicWriter::create(Cursor::new(Vec::new()), duplicate, 1).is_err());
    let mut writer = HicWriter::create(
        Cursor::new(Vec::new()),
        HicHeader::new("g", vec![chrom], vec![1]),
        1,
    )
    .unwrap();
    let overflow = [Record {
        x: 0,
        y: 0,
        value: f32::MAX,
    }; 2];
    assert!(
        writer
            .write_matrix(0, 0, ContactUnit::Bp, 1, &overflow)
            .is_err()
    );
    writer
        .write_matrix(
            0,
            0,
            ContactUnit::Bp,
            1,
            &[Record {
                x: 0,
                y: 0,
                value: 1.0,
            }],
        )
        .unwrap();
    let bytes = writer.finish().unwrap().into_inner();
    let mut reader = HicReader::open(&bytes[..]).unwrap();
    let region = Region::new("s", 0, 10);
    assert!(matches!(
        reader.matrix(ContactQuery::new(region, region, 0), &mut []),
        Err(Error::Invalid { .. })
    ));
    let mut invalid = bytes.clone();
    let pos = invalid.windows(2).position(|w| w == b"s\0").unwrap();
    invalid[pos] = 0xFF;
    assert!(HicReader::open(&invalid[..]).is_err());
    // Repeated chromosome names in a source must not silently overwrite IDs.
    let header = HicHeader::new(
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
    let writer = HicWriter::create(Cursor::new(Vec::new()), header, 1).unwrap();
    let mut duplicate = writer.finish().unwrap().into_inner();
    let pos = duplicate.windows(2).position(|w| w == b"b\0").unwrap();
    duplicate[pos] = b'a';
    assert!(HicReader::open(&duplicate[..]).is_err());
}

#[test]
fn merged_floats_follow_contact_input_order() {
    let header = HicHeader::new(
        "g",
        vec![Chrom {
            name: "s".into(),
            len: 7,
        }],
        vec![1],
    );
    let mut records = Vec::new();
    let mut want = [[0.0f32; 7]; 7];
    for i in 0..1024usize {
        let (x, y) = (i % 7, (i * 3 / 7) % 7);
        let value = [100_000_000.0, 1.0, -100_000_000.0, 3.0][i % 4];
        records.push(Record {
            x: x as u32,
            y: y as u32,
            value,
        });
        let (a, b) = (x.min(y), x.max(y));
        want[a][b] += value;
    }
    let mut writer = HicWriter::create(Cursor::new(Vec::new()), header, 1).unwrap();
    writer
        .write_matrix(0, 0, ContactUnit::Bp, 1, &records)
        .unwrap();
    let bytes = writer.finish().unwrap().into_inner();
    let mut reader = HicReader::open(&bytes[..]).unwrap();
    let region = Region::new("s", 0, 7);
    let mut out = [0.0; 49];
    reader
        .matrix(ContactQuery::new(region, region, 1), &mut out)
        .unwrap();
    for x in 0..7 {
        for y in 0..7 {
            assert_eq!(out[x * 7 + y].to_bits(), want[x.min(y)][x.max(y)].to_bits());
        }
    }
}

#[test]
fn wide_resolution_cache_keeps_unit_indexes_and_rejects_unknown_keys() {
    let header = {
        let mut h = HicHeader::new(
            "test",
            (0..3)
                .map(|i| Chrom {
                    name: format!("c{i}"),
                    len: 1000,
                })
                .collect(),
            (1..=9).rev().collect(),
        );
        h.frag_resolutions = (1..=9).rev().collect();
        h.sites = vec![(10..=990).step_by(10).collect(); 3];
        h
    };
    let mut writer = HicWriter::create(Cursor::new(Vec::new()), header, 1).unwrap();
    let records = [Record {
        x: 0,
        y: 1,
        value: 7.0,
    }];
    for unit in [ContactUnit::Bp, ContactUnit::Frag] {
        assert!(writer.write_matrix(0, 1, unit, 1000, &records).is_err());
        for res in (1..=9).rev() {
            writer.write_matrix(0, 1, unit, res, &records).unwrap();
            assert!(writer.write_matrix(0, 1, unit, res, &records).is_err());
        }
    }
    let bytes = writer.finish().unwrap().into_inner();
    let mut reader = HicReader::open(bytes.as_slice()).unwrap();
    for unit in [ContactUnit::Bp, ContactUnit::Frag] {
        for res in 1..=9 {
            let q = ContactQuery::new(
                Region::new("c0", 0, res),
                Region::new("c1", res, 2 * res),
                res,
            )
            .with_unit(unit);
            let mut out = [0.0];
            reader.matrix(q, &mut out).unwrap();
            assert_eq!(out, [7.0]);
        }
    }
}
