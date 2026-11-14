use std::cell::RefCell;
use std::io::Cursor;

use omicsweave_core::{ContactUnit, ReadAt, Result};
use omicsweave_hic::{
    BlockBins, BlockEntry, BlockLayout, Chrom, Deflater, Depth, HicHeader, HicSidecar, HicWriter,
    Record, ZoomData, decode_block, encode_block,
};

#[test]
fn encoder_matches_literal_dense_and_sparse_wire_layouts() {
    for version in [8, 9] {
        for float in [false, true] {
            for dense in [false, true] {
                let records = if dense {
                    vec![(10, 20), (11, 20), (10, 21), (11, 21)]
                } else {
                    vec![(10, 20), (100, 40)]
                };
                let records: Vec<_> = records
                    .into_iter()
                    .enumerate()
                    .map(|(i, (x, y))| Record {
                        x,
                        y,
                        value: (i + 1) as f32 + if float { 0.5 } else { 0.0 },
                    })
                    .collect();
                let mut expected = (records.len() as i32).to_le_bytes().to_vec();
                expected.extend_from_slice(&10i32.to_le_bytes());
                expected.extend_from_slice(&20i32.to_le_bytes());
                expected.push(u8::from(float));
                if version == 9 {
                    expected.extend_from_slice(&[0, 0]);
                }
                if dense {
                    // Four values occupy a 2-by-2 rectangle, with no holes.
                    expected.extend_from_slice(&[2, 4, 0, 0, 0, 2, 0]);
                    for record in &records {
                        if float {
                            expected.extend_from_slice(&record.value.to_le_bytes());
                        } else {
                            expected.extend_from_slice(&(record.value as i16).to_le_bytes());
                        }
                    }
                } else {
                    // Two rows, each with one x offset and one value.
                    expected.extend_from_slice(&[1, 2, 0]);
                    for (record, offsets) in records.iter().zip([[0, 0], [90, 20]]) {
                        for field in [offsets[1], 1, offsets[0]] {
                            expected.extend_from_slice(&(field as i16).to_le_bytes());
                        }
                        if float {
                            expected.extend_from_slice(&record.value.to_le_bytes());
                        } else {
                            expected.extend_from_slice(&(record.value as i16).to_le_bytes());
                        }
                    }
                }
                let mut encoded = Vec::new();
                encode_block(version, &records, &mut encoded).unwrap();
                assert_eq!(
                    encoded, expected,
                    "v{version}, float={float}, dense={dense}"
                );
                let mut decoded = Vec::new();
                decode_block(version, &expected, |x, y, value| {
                    decoded.push(Record {
                        x: x as u32,
                        y: y as u32,
                        value,
                    });
                })
                .unwrap();
                assert_eq!(decoded, records);
            }
        }
    }
}

#[test]
fn sparse_extreme_coordinates_never_size_an_ineligible_dense_rectangle() {
    let records = [
        Record {
            x: 0,
            y: 0,
            value: 0.5,
        },
        Record {
            x: i32::MAX as u32,
            y: i32::MAX as u32,
            value: 0.5,
        },
    ];
    let mut encoded = Vec::new();
    encode_block(9, &records, &mut encoded).unwrap();
    let mut expected = Vec::new();
    for field in [2i32, 0, 0] {
        expected.extend_from_slice(&field.to_le_bytes());
    }
    expected.extend_from_slice(&[1, 1, 1, 1]);
    expected.extend_from_slice(&2i32.to_le_bytes());
    for coordinate in [0i32, i32::MAX] {
        for field in [coordinate, 1, coordinate] {
            expected.extend_from_slice(&field.to_le_bytes());
        }
        expected.extend_from_slice(&0.5f32.to_le_bytes());
    }
    assert_eq!(encoded, expected);
    assert_eq!(encoded.len(), 52);
    let mut decoded = Vec::new();
    decode_block(9, &encoded, |x, y, value| {
        decoded.push(Record {
            x: x as u32,
            y: y as u32,
            value,
        });
    })
    .unwrap();
    assert_eq!(decoded, records);
}

#[test]
fn wide_rectangles_match_the_union_of_numeric_row_intervals() {
    let blocks: Vec<_> = (0..200_000)
        .filter(|number| number % 3 != 0)
        .map(|number| BlockEntry {
            number,
            position: 0,
            size: 1,
        })
        .collect();
    for intra in [false, true] {
        let zoom = ZoomData {
            chr1: 0,
            chr2: u32::from(!intra),
            unit: ContactUnit::Bp,
            bin_size: 1,
            sum_counts: 0.0,
            block_bin_count: 1,
            block_column_count: 16,
            blocks: blocks.clone(),
        };
        for (x, y) in [(3..40, 100..200), (3..40, 100..u32::MAX)] {
            let x_end = u64::from(x.end);
            let y_end = u64::from(y.end);
            let expected: Vec<_> = blocks
                .iter()
                .copied()
                .filter(|block| {
                    let number = u64::from(block.number);
                    (u64::from(y.start) * 16 + u64::from(x.start)..=y_end * 16 + x_end)
                        .contains(&number)
                        || (intra
                            && (u64::from(x.start) * 16 + u64::from(y.start)..=x_end * 16 + y_end)
                                .contains(&number))
                })
                .collect();
            assert_eq!(zoom.blocks_for(8, Depth::Log(2), x, y), expected);
        }
    }
    let zoom = ZoomData {
        chr1: 0,
        chr2: 1,
        unit: ContactUnit::Bp,
        bin_size: 1,
        sum_counts: 0.0,
        block_bin_count: 1,
        block_column_count: u32::MAX - 3,
        blocks,
    };
    // Even the first numeric interval starts above i64::MAX; there is no
    // block number in the format's u32 range and no reason to visit rows.
    assert!(
        zoom.blocks_for(8, Depth::Log(2), 0..u32::MAX, u32::MAX - 2..u32::MAX)
            .is_empty()
    );
}

struct Counted<'a> {
    bytes: &'a [u8],
    calls: RefCell<Vec<(u64, usize)>>,
}

impl ReadAt for Counted<'_> {
    fn len(&self) -> u64 {
        self.bytes.len() as u64
    }

    fn read_exact_at(&self, out: &mut [u8], offset: u64) -> Result<()> {
        self.calls.borrow_mut().push((offset, out.len()));
        self.bytes.read_exact_at(out, offset)
    }
}

#[test]
fn sidecar_coalesces_each_adjacent_multiresolution_block_run() {
    for version in [8, 9] {
        let mut header = HicHeader::new(
            "g",
            ["a", "b"]
                .into_iter()
                .map(|name| Chrom {
                    name: name.into(),
                    len: 8192,
                })
                .collect(),
            vec![1, 2, 4],
        );
        header.version = version;
        let layout = BlockLayout {
            intra: BlockBins::Fixed(4),
            inter: BlockBins::Fixed(4),
        };
        let mut writer =
            HicWriter::create_with(Cursor::new(Vec::new()), header, 1, layout).unwrap();
        for res in [1, 2, 4] {
            let records: Vec<_> = (0..512)
                .map(|i| Record {
                    x: i,
                    y: i + 1,
                    value: (i % 17) as f32 + 0.5,
                })
                .collect();
            for chr2 in [0, 1] {
                writer
                    .write_matrix(0, chr2, ContactUnit::Bp, res, &records)
                    .unwrap();
            }
        }
        let (out, expected) = writer.finish_with_sidecar().unwrap();
        let bytes = out.into_inner();
        let source = Counted {
            bytes: &bytes,
            calls: RefCell::new(Vec::new()),
        };
        assert_eq!(HicSidecar::build(&source).unwrap(), expected);
        for matrix in &expected.matrices {
            assert!(matrix.blocks.len() > 10);
            let first = matrix.blocks.first().unwrap();
            let last = matrix.blocks.last().unwrap();
            let end = last.position + u64::from(last.size);
            let reads: Vec<_> = source
                .calls
                .borrow()
                .iter()
                .copied()
                .filter(|&(offset, _)| (first.position..end).contains(&offset))
                .collect();
            assert_eq!(reads, [(first.position, (end - first.position) as usize)]);
        }
        source.calls.borrow_mut().clear();
        let filtered =
            HicSidecar::build_filtered(&source, |_, chr2, _, res| chr2 == 1 && res == 2).unwrap();
        assert_eq!(
            filtered.matrices,
            expected
                .matrices
                .iter()
                .filter(|matrix| matrix.chr2 == 1 && matrix.bin_size == 2)
                .cloned()
                .collect::<Vec<_>>()
        );
        for matrix in &expected.matrices {
            if matrix.chr2 == 1 && matrix.bin_size == 2 {
                continue;
            }
            let start = matrix.blocks.first().unwrap().position;
            let last = matrix.blocks.last().unwrap();
            let end = last.position + u64::from(last.size);
            assert!(
                source
                    .calls
                    .borrow()
                    .iter()
                    .all(|&(offset, _)| !(start..end).contains(&offset))
            );
        }
        let mut corrupt = bytes.clone();
        let last = expected.matrices.last().unwrap().blocks.last().unwrap();
        corrupt[last.position as usize] ^= 0xff;
        assert!(HicSidecar::build(corrupt.as_slice()).is_err());
    }
}

#[test]
fn sidecar_preserves_index_order_and_bounds_physical_read_batches() {
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
    let mut writer = HicWriter::create_with(
        Cursor::new(Vec::new()),
        header,
        0,
        BlockLayout {
            intra: BlockBins::Fixed(1),
            inter: BlockBins::Fixed(1),
        },
    )
    .unwrap();
    let record = Record {
        x: 0,
        y: 0,
        value: 1.0,
    };
    writer
        .write_matrix(
            0,
            1,
            ContactUnit::Bp,
            1,
            &(0..3)
                .map(|i| Record {
                    x: i,
                    y: i,
                    ..record
                })
                .collect::<Vec<_>>(),
        )
        .unwrap();
    let (out, expected) = writer.finish_with_sidecar().unwrap();
    let original = out.into_inner();
    assert_eq!(expected.matrices[0].blocks.len(), 3);
    let mut raw = Vec::new();
    encode_block(9, &[record], &mut raw).unwrap();
    let mut packed = Vec::new();
    Deflater::new(0).deflate(&raw, &mut packed).unwrap();
    let raw_len = raw.len() as u32;
    let mut huge = raw.clone();
    // Trailing decompressed bytes are legal to this raw-block decoder.
    // Level zero makes the compressed stream itself larger than 1 MiB.
    huge.resize((1 << 20) + 32, 0);
    let mut packed_huge = Vec::new();
    Deflater::new(0).deflate(&huge, &mut packed_huge).unwrap();
    assert!(packed_huge.len() > 1 << 20);
    for case in ["gap-reverse", "overlap", "oversized", "outside"] {
        let mut bytes = original.clone();
        let base = bytes.len() as u64;
        let mut expected = expected.clone();
        let mut positions = Vec::new();
        let mut lengths = Vec::new();
        for i in 0..3 {
            let data = if case == "oversized" && i == 1 {
                &packed_huge
            } else {
                &packed
            };
            positions.push(bytes.len() as u64);
            lengths.push(data.len());
            bytes.extend_from_slice(data);
            if case == "gap-reverse" {
                bytes.extend_from_slice(&[0xaa; 17]);
            }
        }
        if case == "gap-reverse" {
            positions.reverse();
        } else if case == "overlap" {
            positions.fill(base);
        } else if case == "outside" {
            positions[1] = bytes.len() as u64 + 4096;
        }
        for (i, block) in expected.matrices[0].blocks.iter_mut().enumerate() {
            let mut key = (block.number as i32).to_le_bytes().to_vec();
            key.extend_from_slice(&block.position.to_le_bytes());
            key.extend_from_slice(&block.size.to_le_bytes());
            let field = original
                .windows(key.len())
                .position(|window| window == key)
                .unwrap();
            bytes[field + 4..field + 12].copy_from_slice(&positions[i].to_le_bytes());
            bytes[field + 12..field + 16].copy_from_slice(&(lengths[i] as u32).to_le_bytes());
            block.position = positions[i];
            block.size = lengths[i] as u32;
            block.inflated = if case == "oversized" && i == 1 {
                huge.len() as u32
            } else {
                raw_len
            };
            block.records = 1;
        }
        let source = Counted {
            bytes: &bytes,
            calls: RefCell::new(Vec::new()),
        };
        let result = HicSidecar::build(&source);
        if case == "outside" {
            assert!(result.is_err());
            assert!(
                source
                    .calls
                    .borrow()
                    .iter()
                    .all(|&(at, _)| at != positions[1])
            );
            continue;
        }
        assert_eq!(result.unwrap(), expected, "{case}");
        let actual: Vec<_> = source
            .calls
            .borrow()
            .iter()
            .copied()
            .filter(|&(at, _)| at >= base)
            .collect();
        let expected_reads: Vec<_> = positions.into_iter().zip(lengths).collect();
        assert_eq!(actual, expected_reads, "{case}");
    }
}
