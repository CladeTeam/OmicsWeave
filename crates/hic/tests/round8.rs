use std::collections::BTreeMap;
use std::io::Cursor;

use omicsweave_core::sidecar::SourceIdentity;
use omicsweave_core::{ContactQuery, ContactRead, ContactUnit, Region};
use omicsweave_hic::{
    BlockBins, BlockLayout, BlockSize, Chrom, HicHeader, HicReader, HicSidecar, HicWriter,
    MatrixBlocks, Record,
};

#[test]
fn sidecar_checks_index_invariants_even_with_a_valid_container_checksum() {
    let identity = SourceIdentity {
        len: 100,
        fingerprint: 7,
    };
    let valid = HicSidecar {
        version: 9,
        matrices: vec![MatrixBlocks {
            chr1: 0,
            chr2: 1,
            unit: ContactUnit::Bp,
            bin_size: 1,
            block_bin_count: 16,
            block_column_count: 8,
            blocks: [1, 2]
                .into_iter()
                .map(|number| BlockSize {
                    number,
                    position: 80,
                    size: 20,
                    inflated: 64,
                    records: 2,
                })
                .collect(),
        }],
    };
    assert_eq!(
        HicSidecar::parse(&valid.to_bytes(identity), identity).unwrap(),
        valid
    );
    for case in 0..8 {
        let mut invalid = valid.clone();
        let matrix = &mut invalid.matrices[0];
        match case {
            0 => matrix.blocks.reverse(),
            1 => matrix.blocks[1].number = matrix.blocks[0].number,
            2 => matrix.blocks[0].position = 81,
            3 => matrix.blocks[0].position = u64::MAX - 10,
            4 => matrix.chr1 = 2,
            5 => matrix.bin_size = 0,
            6 => matrix.block_bin_count = 0,
            7 => matrix.block_column_count = 0,
            _ => unreachable!(),
        }
        assert!(
            HicSidecar::parse(&invalid.to_bytes(identity), identity).is_err(),
            "case {case}"
        );
    }
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
