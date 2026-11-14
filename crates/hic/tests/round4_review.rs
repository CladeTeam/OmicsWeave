use std::cell::RefCell;
use std::collections::BTreeSet;
use std::io::Cursor;
use std::rc::Rc;

use omicsweave_core::{
    ContactEntry, ContactQuery, ContactRead, ContactUnit, Error, ReadAt, Region,
};
use omicsweave_hic::{
    BlockEntry, Chrom, Depth, HicHeader, HicReader, HicWriter, Record, ZoomData, decode_block,
};

#[derive(Clone)]
struct EditableSource(Rc<RefCell<Vec<u8>>>);

impl ReadAt for EditableSource {
    fn len(&self) -> u64 {
        self.0.borrow().len() as u64
    }

    fn read_exact_at(&self, buf: &mut [u8], offset: u64) -> omicsweave_core::Result<()> {
        self.0.borrow().as_slice().read_exact_at(buf, offset)
    }
}

#[test]
fn conflicting_matrix_identity_invalidates_entries_and_allows_retry() {
    let header = HicHeader::new(
        "g",
        ["a", "b"]
            .into_iter()
            .map(|name| Chrom {
                name: name.into(),
                len: 16,
            })
            .collect(),
        vec![1],
    );
    let mut writer = HicWriter::create(Cursor::new(Vec::new()), header, 1).unwrap();
    for (c2, x, y, value) in [(0, 1, 1, 5.0), (1, 2, 3, 7.0)] {
        writer
            .write_matrix(0, c2, ContactUnit::Bp, 1, &[Record { x, y, value }])
            .unwrap();
    }
    let original = writer.finish().unwrap().into_inner();
    let own = ContactQuery::new(Region::new("a", 0, 8), Region::new("a", 0, 8), 1);
    let inter = ContactQuery::new(Region::new("a", 0, 8), Region::new("b", 0, 8), 1);
    for bad_pair in [(0i32, 0i32), (1, 0), (0, 77)] {
        let source = EditableSource(Rc::new(RefCell::new(original.clone())));
        let mut reader = HicReader::open(source.clone()).unwrap();
        assert_eq!(reader.query_entries(own).unwrap(), 1);
        let matrix = reader.footer().matrices[&(0, 1)] as usize;
        {
            let mut bytes = source.0.borrow_mut();
            bytes[matrix..matrix + 4].copy_from_slice(&bad_pair.0.to_le_bytes());
            bytes[matrix + 4..matrix + 8].copy_from_slice(&bad_pair.1.to_le_bytes());
        }
        assert!(matches!(
            reader.query_entries(inter),
            Err(Error::Invalid { .. })
        ));
        assert!(matches!(
            reader.write_entries(&mut [ContactEntry::default()]),
            Err(Error::NotQueried)
        ));
        source.0.borrow_mut()[matrix..matrix + 8].copy_from_slice(&original[matrix..matrix + 8]);
        assert_eq!(reader.query_entries(inter).unwrap(), 1);
        let mut entries = [ContactEntry::default()];
        reader.write_entries(&mut entries).unwrap();
        assert_eq!(
            entries[0],
            ContactEntry {
                row: 2,
                col: 3,
                value: 7.0
            }
        );
    }
}

fn put_int(bytes: &mut Vec<u8>, wide: bool, value: i32) {
    if wide {
        bytes.extend_from_slice(&value.to_le_bytes());
    } else {
        bytes.extend_from_slice(&(value as i16).to_le_bytes());
    }
}

#[test]
fn sparse_coordinates_are_checked_for_every_field_layout() {
    for wide_x in [false, true] {
        for wide_y in [false, true] {
            for float in [false, true] {
                for bad_x in [false, true] {
                    let mut bytes = 1i32.to_le_bytes().to_vec();
                    let (x0, y0) = if bad_x { (i32::MAX, 0) } else { (0, i32::MIN) };
                    bytes.extend_from_slice(&x0.to_le_bytes());
                    bytes.extend_from_slice(&y0.to_le_bytes());
                    bytes.extend_from_slice(&[
                        u8::from(float),
                        u8::from(wide_x),
                        u8::from(wide_y),
                        1,
                    ]);
                    put_int(&mut bytes, wide_y, 1);
                    put_int(&mut bytes, wide_y, if bad_x { 0 } else { -1 });
                    put_int(&mut bytes, wide_x, 1);
                    put_int(&mut bytes, wide_x, i32::from(bad_x));
                    if float {
                        bytes.extend_from_slice(&3.0f32.to_le_bytes());
                    } else {
                        bytes.extend_from_slice(&3i16.to_le_bytes());
                    }
                    let mut callbacks = 0;
                    assert!(decode_block(9, &bytes, |_, _, _| callbacks += 1).is_err());
                    assert_eq!(callbacks, 0);
                }
            }
        }
    }
}

fn stored_number(
    version: i32,
    depth: Depth,
    intra: bool,
    bb: u32,
    columns: u32,
    x: u32,
    y: u32,
) -> u32 {
    let bb = u64::from(bb);
    let x = u64::from(x);
    let y = u64::from(y);
    let row = if version == 9 && intra {
        let distance = x.abs_diff(y) as f64 / std::f64::consts::SQRT_2 / bb as f64;
        match depth {
            Depth::Log(base) => (1.0 + distance).log(f64::from(base)).floor() as u64,
            Depth::Linear(step) => (distance / f64::from(step)).floor() as u64,
        }
    } else {
        y / bb
    };
    let col = if version == 9 && intra {
        (x + y) / (2 * bb)
    } else {
        x / bb
    };
    u32::try_from(row * u64::from(columns) + col).unwrap()
}

#[test]
fn candidate_blocks_cover_actual_contacts_and_their_mirrors() {
    for version in [8, 9] {
        for intra in [false, true] {
            for depth in [Depth::Log(2), Depth::Log(7), Depth::Linear(3)] {
                for (base, bb, columns) in [
                    (0, 1, 2),
                    (0, 3, 20),
                    (0, 7, 3),
                    (u32::MAX - 32, 1_000_000_000, 3),
                ] {
                    let records: Vec<_> = (base..base + 24)
                        .flat_map(|x| {
                            (base..base + 24)
                                .filter(move |&y| !intra || x <= y)
                                .map(move |y| (x, y))
                        })
                        .collect();
                    let numbers: BTreeSet<_> = records
                        .iter()
                        .map(|&(x, y)| stored_number(version, depth, intra, bb, columns, x, y))
                        .collect();
                    let zoom = ZoomData {
                        chr1: 0,
                        chr2: u32::from(!intra),
                        unit: ContactUnit::Bp,
                        bin_size: 1,
                        sum_counts: 0.0,
                        block_bin_count: bb,
                        block_column_count: columns,
                        blocks: numbers
                            .into_iter()
                            .map(|number| BlockEntry {
                                number,
                                position: 0,
                                size: 0,
                            })
                            .collect(),
                    };
                    for (x, y, width, height) in [
                        (0, 0, 1, 1),
                        (3, 9, 5, 7),
                        (9, 3, 7, 5),
                        (7, 8, 9, 2),
                        (1, 1, 21, 21),
                    ] {
                        let xr = base + x..base + x + width;
                        let yr = base + y..base + y + height;
                        let blocks = zoom.blocks_for(version, depth, xr.clone(), yr.clone());
                        assert!(blocks.windows(2).all(|w| w[0].number < w[1].number));
                        let selected: BTreeSet<_> = blocks.iter().map(|b| b.number).collect();
                        for &(rx, ry) in &records {
                            if (xr.contains(&rx) && yr.contains(&ry))
                                || (intra && xr.contains(&ry) && yr.contains(&rx))
                            {
                                assert!(
                                    selected.contains(&stored_number(
                                        version, depth, intra, bb, columns, rx, ry
                                    )),
                                    "v{version}, {depth:?}, ({rx},{ry}), {xr:?} x {yr:?}"
                                );
                            }
                        }
                    }
                }
            }
        }
    }
}
