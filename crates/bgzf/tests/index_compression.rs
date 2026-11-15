use std::io::Write;

use flate2::{Compression, GzBuilder};
use omicsweave_bgzf::{BinIndex, Chunk, EOF_BLOCK, check_level, deflate_block};

fn csi() -> Vec<u8> {
    let mut bytes = b"CSI\x01".to_vec();
    for n in [14i32, 5, 0, 1, 1, 0] {
        bytes.extend_from_slice(&n.to_le_bytes());
    }
    bytes.extend_from_slice(&0u64.to_le_bytes());
    bytes.extend_from_slice(&1i32.to_le_bytes());
    bytes.extend_from_slice(&65536u64.to_le_bytes());
    bytes.extend_from_slice(&131072u64.to_le_bytes());
    bytes
}

fn bai() -> Vec<u8> {
    let mut bytes = b"BAI\x01".to_vec();
    for n in [1i32, 1, 0, 1] {
        bytes.extend_from_slice(&n.to_le_bytes());
    }
    bytes.extend_from_slice(&65536u64.to_le_bytes());
    bytes.extend_from_slice(&131072u64.to_le_bytes());
    bytes.extend_from_slice(&0i32.to_le_bytes());
    bytes
}

fn gzip(bytes: &[u8], extra: bool) -> Vec<u8> {
    let builder = if extra {
        GzBuilder::new().extra(vec![b'X', b'Y', 1, 0, 7])
    } else {
        GzBuilder::new()
    };
    let mut writer = builder.write(Vec::new(), Compression::fast());
    writer.write_all(bytes).unwrap();
    writer.finish().unwrap()
}

#[test]
fn plain_bgzf_and_concatenated_gzip_indexes_have_identical_chunks() {
    for plain in [csi(), bai()] {
        let mut compressor = libdeflater::Compressor::new(check_level(6).unwrap());
        let mut bgzf = Vec::new();
        deflate_block(&mut compressor, &plain, &mut bgzf).unwrap();
        bgzf.extend_from_slice(&EOF_BLOCK);
        let mut concatenated = gzip(&plain[..19], false);
        concatenated.extend_from_slice(&gzip(&plain[19..], true));
        for raw in [
            plain.clone(),
            bgzf,
            gzip(&plain, false),
            gzip(&plain, true),
            concatenated,
        ] {
            let index = BinIndex::parse(&raw).unwrap();
            assert_eq!(
                index.chunks(0, 100, 200),
                vec![Chunk {
                    start: 65536,
                    end: 131072
                }]
            );
        }
    }
}

#[test]
fn compressed_index_crc_and_truncation_errors_are_not_accepted() {
    let plain = csi();
    let mut raw = gzip(&plain, false);
    let crc = raw.len() - 8;
    raw[crc] ^= 1;
    assert!(BinIndex::parse(&raw).is_err());
    let raw = gzip(&plain, true);
    assert!(BinIndex::parse(&raw[..raw.len() - 3]).is_err());
    let mut compressor = libdeflater::Compressor::new(check_level(6).unwrap());
    let mut bgzf = Vec::new();
    deflate_block(&mut compressor, &plain, &mut bgzf).unwrap();
    let mut bad_crc = bgzf.clone();
    let crc = bad_crc.len() - 8;
    bad_crc[crc] ^= 1;
    assert!(BinIndex::parse(&bad_crc).is_err());
    assert!(BinIndex::parse(&bgzf[..bgzf.len() - 1]).is_err());
    bgzf[16] ^= 1;
    assert!(BinIndex::parse(&bgzf).is_err());
}

#[test]
fn duplicate_bgzf_size_fields_cannot_fall_back_to_generic_gzip() {
    let mut compressor = libdeflater::Compressor::new(check_level(6).unwrap());
    let mut bgzf = Vec::new();
    deflate_block(&mut compressor, &csi(), &mut bgzf).unwrap();
    let mut duplicate = bgzf[..18].to_vec();
    duplicate.extend_from_slice(&bgzf[12..18]);
    duplicate.extend_from_slice(&bgzf[18..]);
    duplicate[10..12].copy_from_slice(&12u16.to_le_bytes());
    let size = (duplicate.len() - 1) as u16;
    duplicate[16..18].copy_from_slice(&size.to_le_bytes());
    duplicate[22..24].copy_from_slice(&size.to_le_bytes());
    assert!(BinIndex::parse(&duplicate).is_err());
}
