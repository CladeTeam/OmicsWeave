use omicsweave_core::{ContactQuery, ContactRead, ContactUnit, Region};
use omicsweave_hic::{BlockBins, BlockLayout, Chrom, HicHeader, HicReader, HicWriter, Record};
use std::collections::BTreeMap;
use std::io::Cursor;

fn one_chrom() -> HicHeader {
    HicHeader::new(
        "g",
        vec![Chrom {
            name: "s".into(),
            len: 10_000,
        }],
        vec![10],
    )
}

#[test]
fn misplaced_blocks_change_nothing_and_misplaced_records_fail_the_writer() {
    let header = one_chrom();
    let mut w = HicWriter::create(Cursor::new(Vec::new()), header.clone(), 1).unwrap();
    let layout = w.matrix_layout(0, 0, ContactUnit::Bp, 10).unwrap();
    let cell = |x, y| {
        let (number, x, y) = layout.place(x, y);
        (number, Record { x, y, value: 2.0 })
    };
    let (near, a) = cell(3, 5);
    let (far, b) = cell(10, 900);
    assert!(near < far);
    assert!(w.add_block(near, vec![a]).is_err(), "no map has begun");
    w.begin_matrix(&layout).unwrap();
    assert!(w.begin_matrix(&layout).is_err(), "one map at a time");
    w.add_block(far, vec![b]).unwrap();
    assert!(
        w.add_block(near, vec![a]).is_err(),
        "block numbers increase"
    );
    w.end_matrix().unwrap();
    assert!(w.begin_matrix(&layout).is_err(), "written twice");
    let streamed = w.finish().unwrap().into_inner();
    let mut whole = HicWriter::create(Cursor::new(Vec::new()), header.clone(), 1).unwrap();
    whole.write_matrix(0, 0, ContactUnit::Bp, 10, &[b]).unwrap();
    assert_eq!(streamed, whole.finish().unwrap().into_inner());

    for bad in [
        Record {
            x: 10,
            y: 900,
            value: 1.0,
        }, // another block's cell
        Record {
            x: 900,
            y: 10,
            value: 1.0,
        }, // below the diagonal
        Record {
            x: 5,
            y: 5_000,
            value: 1.0,
        }, // outside the map
        Record {
            x: 3,
            y: 5,
            value: f32::NAN,
        },
    ] {
        let mut w = HicWriter::create(Cursor::new(Vec::new()), header.clone(), 1).unwrap();
        w.begin_matrix(&layout).unwrap();
        assert!(w.add_block(near, vec![a, bad]).is_err());
        assert!(w.end_matrix().is_err());
        assert!(w.finish().is_err());
    }
}

/// A map streamed on a pool reads back: each cell holds the sum of its
/// records, mirrored across the diagonal.
#[test]
fn streamed_map_reads_back() {
    let header = one_chrom();
    let out = Cursor::new(Vec::new());
    let pool = omicsweave_core::exec::Pool::new(4);
    let exec = omicsweave_core::exec::Exec::Pool {
        pool: &pool,
        limit: 4,
    };
    let mut w = HicWriter::create_exec(out, header, 1, BlockLayout::default(), exec).unwrap();
    let layout = w.matrix_layout(0, 0, ContactUnit::Bp, 10).unwrap();
    let mut blocks: BTreeMap<u32, Vec<Record>> = BTreeMap::new();
    let mut want = BTreeMap::new();
    for i in 0..20_000u32 {
        let (x, y) = (i * 7 % 1001, i * 13 % 1001);
        let (number, x, y) = layout.place(x, y);
        blocks
            .entry(number)
            .or_default()
            .push(Record { x, y, value: 1.0 });
        *want.entry((x, y)).or_insert(0.0f32) += 1.0;
    }
    w.begin_matrix(&layout).unwrap();
    for (number, records) in blocks {
        w.add_block(number, records).unwrap();
    }
    w.end_matrix().unwrap();
    let bytes = w.finish().unwrap().into_inner();
    let mut reader = HicReader::open(bytes.as_slice()).unwrap();
    let region = Region::new("s", 0, 10_010);
    let mut tile = vec![0.0f32; 1001 * 1001];
    reader
        .matrix(ContactQuery::new(region, region, 10), &mut tile)
        .unwrap();
    for (&(x, y), &v) in &want {
        assert_eq!(tile[x as usize * 1001 + y as usize], v);
        assert_eq!(tile[y as usize * 1001 + x as usize], v);
    }
    assert_eq!(tile.iter().filter(|&&v| v != 0.0).count(), {
        want.keys()
            .map(|&(x, y)| if x == y { 1 } else { 2 })
            .sum::<usize>()
    });
}

#[test]
fn whole_and_streamed_maps_preserve_cell_sum_order_and_wire_bytes() {
    let mut input = Vec::new();
    for y in (0..80).rev() {
        for x in (0..=y).rev() {
            for value in [16_777_216.0, 1.0, -16_777_216.0, 0.5, 0.0] {
                input.push(Record { x: y, y: x, value });
            }
        }
    }
    input.extend([
        Record {
            x: 3,
            y: 7,
            value: -0.5,
        },
        Record {
            x: 5,
            y: 5,
            value: -0.5,
        },
    ]);
    for version in [8, 9] {
        let mut header = HicHeader::new(
            "order",
            vec![Chrom {
                name: "a".into(),
                len: 100,
            }],
            vec![1],
        );
        header.version = version;
        for blocks in [
            BlockLayout::juicer(),
            BlockLayout {
                intra: BlockBins::Fixed(8),
                inter: BlockBins::Fixed(8),
            },
        ] {
            let new = || {
                HicWriter::create_with(Cursor::new(Vec::new()), header.clone(), 1, blocks).unwrap()
            };
            let mut whole = new();
            whole
                .write_matrix(0, 0, ContactUnit::Bp, 1, &input)
                .unwrap();
            let actual = whole.finish_with_sidecar().unwrap();
            let mut reader = HicReader::open(actual.0.get_ref().as_slice()).unwrap();
            let mut matrix = vec![0.0; 80 * 80];
            reader
                .matrix(
                    ContactQuery::new(Region::new("a", 0, 80), Region::new("a", 0, 80), 1),
                    &mut matrix,
                )
                .unwrap();
            for row in 0..80 {
                for col in 0..80 {
                    let canceled = matches!((row, col), (3, 7) | (7, 3) | (5, 5));
                    assert_eq!(matrix[row * 80 + col], if canceled { 0.0 } else { 0.5 });
                }
            }
            let mut stream = new();
            let layout = stream.matrix_layout(0, 0, ContactUnit::Bp, 1).unwrap();
            let mut by_block: BTreeMap<u32, Vec<Record>> = BTreeMap::new();
            for record in &input {
                let (number, x, y) = layout.place(record.x, record.y);
                by_block.entry(number).or_default().push(Record {
                    x,
                    y,
                    value: record.value,
                });
            }
            stream.begin_matrix(&layout).unwrap();
            for (number, records) in by_block {
                stream.add_block(number, records).unwrap();
            }
            stream.end_matrix().unwrap();
            let expected = stream.finish_with_sidecar().unwrap();
            assert_eq!(actual.0.into_inner(), expected.0.into_inner());
            assert_eq!(actual.1, expected.1);
        }
    }
}
