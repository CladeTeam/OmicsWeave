use std::collections::BTreeMap;
use std::io::Cursor;

use omicsweave_core::{ContactQuery, ContactRead, ContactUnit, Region};
use omicsweave_hic::{
    BlockLayout, BuildOptions, Chrom, HicBuilder, HicHeader, HicReader, HicSidecar, HicWriter,
    Record,
};

type Contact = (u32, u32, u32, u32, f32);

fn header(version: i32, resolutions: usize) -> HicHeader {
    let mut header = HicHeader::new(
        "test",
        (0..3)
            .map(|i| Chrom {
                name: format!("c{i}"),
                len: 1000,
            })
            .collect(),
        (1..=resolutions as u32).rev().collect(),
    );
    header.version = version;
    header.frag_resolutions = (1..=resolutions as u32).rev().collect();
    header.sites = vec![(10..=990).step_by(10).collect(); 3];
    header
}

// Bins each resolution directly from the input and hands the unmerged
// records to the public writer, in the header's unit/resolution/pair order.
fn expected(
    header: &HicHeader,
    contacts: &[Contact],
    blocks: BlockLayout,
) -> (Vec<u8>, HicSidecar) {
    let mut writer =
        HicWriter::create_with(Cursor::new(Vec::new()), header.clone(), 1, blocks).unwrap();
    for unit in [ContactUnit::Bp, ContactUnit::Frag] {
        for &res in header.resolutions(unit) {
            let mut pairs: BTreeMap<(u32, u32), Vec<Record>> = BTreeMap::new();
            for &(a, p, b, q, value) in contacts {
                let (a, p, b, q) = if a <= b { (a, p, b, q) } else { (b, q, a, p) };
                let position = |c: u32, p: u32| {
                    let p = p.min(header.chroms[c as usize].len as u32);
                    match unit {
                        ContactUnit::Bp => p,
                        ContactUnit::Frag => {
                            header.sites[c as usize].iter().filter(|&&s| s <= p).count() as u32
                        }
                    }
                };
                let (mut x, mut y) = (position(a, p) / res, position(b, q) / res);
                if a == b && x > y {
                    std::mem::swap(&mut x, &mut y);
                }
                pairs
                    .entry((a, b))
                    .or_default()
                    .push(Record { x, y, value });
            }
            for ((a, b), records) in pairs {
                writer.write_matrix(a, b, unit, res, &records).unwrap();
            }
        }
    }
    let (out, sidecar) = writer.finish_with_sidecar().unwrap();
    (out.into_inner(), sidecar)
}

#[test]
fn grouped_builder_preserves_header_order_bytes_and_float_sums() {
    let contacts = [
        (2, 19, 0, 31, 16_777_216.0),
        (0, 31, 2, 19, 1.0),
        (2, 19, 0, 31, -16_777_216.0),
        (0, 31, 2, 19, 0.5),
        (1, 9999, 1, 22, 7.0),
        (1, 22, 1, 1000, -2.0),
        (2, 3, 1, 8, 0.0),
        (0, 32, 0, 1, 3.0),
    ];
    for version in [8, 9] {
        for (count, blocks) in [1, 2, 4, 5, 8, 9, 17].into_iter().flat_map(|count| {
            [BlockLayout::default(), BlockLayout::juicer()]
                .into_iter()
                .map(move |blocks| (count, blocks))
        }) {
            let header = header(version, count);
            let (want, expected_sidecar) = expected(&header, &contacts, blocks);
            let mut builder = HicBuilder::new(header).unwrap();
            for &(a, p, b, q, v) in &contacts {
                builder.add(a, p, b, q, v).unwrap();
            }
            let (out, sidecar) = builder
                .finish_with_sidecar(Cursor::new(Vec::new()), &BuildOptions { level: 1, blocks })
                .unwrap();
            let bytes = out.into_inner();
            assert_eq!(bytes, want, "version={version}, resolutions={count}");
            assert_eq!(sidecar, expected_sidecar);
            assert_eq!(sidecar, HicSidecar::build(bytes.as_slice()).unwrap());
        }
    }
}

#[test]
fn wide_resolution_cache_keeps_unit_indexes_and_rejects_unknown_keys() {
    let header = header(9, 9);
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

#[test]
fn empty_builder_preserves_declared_resolutions_without_maps() {
    for count in [0, 1, 8, 9, 64] {
        let header = header(9, count);
        let builder = HicBuilder::new(header.clone()).unwrap();
        let (actual, sidecar) = builder
            .finish_with_sidecar(
                Cursor::new(Vec::new()),
                &BuildOptions {
                    level: 1,
                    ..Default::default()
                },
            )
            .unwrap();
        let expected = HicWriter::create(Cursor::new(Vec::new()), header, 1)
            .unwrap()
            .finish()
            .unwrap();
        assert_eq!(actual.into_inner(), expected.into_inner());
        assert!(sidecar.matrices.is_empty());
    }
}
