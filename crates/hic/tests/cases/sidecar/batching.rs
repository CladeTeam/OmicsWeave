use crate::support::io::Counted;
use omicsweave_core::ContactUnit;
use omicsweave_hic::{
    BlockBins, BlockLayout, Chrom, Deflater, HicHeader, HicSidecar, HicWriter, Record, encode_block,
};
use std::io::Cursor;

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
        let source = Counted::new(bytes.as_slice());
        assert_eq!(HicSidecar::build(&source).unwrap(), expected);
        for matrix in &expected.matrices {
            assert!(matrix.blocks.len() > 10);
            let first = matrix.blocks.first().unwrap();
            let last = matrix.blocks.last().unwrap();
            let end = last.position + u64::from(last.size);
            let reads: Vec<_> = source
                .calls()
                .iter()
                .copied()
                .filter(|&(offset, _)| (first.position..end).contains(&offset))
                .collect();
            assert_eq!(reads, [(first.position, (end - first.position) as usize)]);
        }
        source.reset();
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
                    .calls()
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
        let source = Counted::new(bytes.as_slice());
        let result = HicSidecar::build(&source);
        if case == "outside" {
            assert!(result.is_err());
            assert!(source.calls().iter().all(|&(at, _)| at != positions[1]));
            continue;
        }
        assert_eq!(result.unwrap(), expected, "{case}");
        let actual: Vec<_> = source
            .calls()
            .iter()
            .copied()
            .filter(|&(at, _)| at >= base)
            .collect();
        let expected_reads: Vec<_> = positions.into_iter().zip(lengths).collect();
        assert_eq!(actual, expected_reads, "{case}");
    }
}
