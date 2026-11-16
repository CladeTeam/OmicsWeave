//! Sidecars: built from files of every block layout (the test writer's v8 and v9
//! files and the crate's own), checked block by block against the block
//! index and an independent inflate and decode; the writers' sidecars equal
//! the ones built from what they wrote.

use crate::support as common;

use std::collections::BTreeSet;
use std::io::Cursor;

use common::{Rng, random_spec};
use omicsweave_core::sidecar::{SourceIdentity, sidecar_path};
use omicsweave_core::{ContactUnit, Error};
use omicsweave_hic::{
    BuildOptions, Chrom, HicBuilder, HicHeader, HicReader, HicSidecar, HicWriter, Inflater, Record,
    decode_block,
};
use proptest::prelude::*;

/// Checks `sidecar` against the file: one matrix per stored resolution of
/// the master index's pairs, and each block's position, size, inflated
/// length and decoded record count.
fn check_against_file(bytes: &[u8], sidecar: &HicSidecar) {
    let mut reader = HicReader::open(bytes).unwrap();
    let version = reader.header().version;
    assert_eq!(sidecar.version, version);
    let pairs: Vec<(u32, u32)> = reader.footer().matrices.keys().copied().collect();
    let mut expected = BTreeSet::new();
    for (c1, c2) in pairs {
        for unit in [ContactUnit::Bp, ContactUnit::Frag] {
            for &res in reader.header().resolutions(unit).to_vec().iter() {
                if reader.zoom(c1, c2, unit, res).unwrap().is_some() {
                    expected.insert((c1, c2, unit, res));
                }
            }
        }
    }
    let found: BTreeSet<_> = sidecar
        .matrices
        .iter()
        .map(|m| (m.chr1, m.chr2, m.unit, m.bin_size))
        .collect();
    assert_eq!(found, expected);

    let mut inflater = Inflater::new();
    let mut block = Vec::new();
    for m in &sidecar.matrices {
        let zoom = reader
            .zoom(m.chr1, m.chr2, m.unit, m.bin_size)
            .unwrap()
            .unwrap()
            .clone();
        assert_eq!(
            (m.block_bin_count, m.block_column_count),
            (zoom.block_bin_count, zoom.block_column_count)
        );
        assert_eq!(m.blocks.len(), zoom.blocks.len());
        for (b, z) in m.blocks.iter().zip(&zoom.blocks) {
            assert_eq!(
                (b.number, b.position, b.size),
                (z.number, z.position, z.size)
            );
            inflater
                .inflate(
                    &bytes[z.range().start as usize..z.range().end as usize],
                    &mut block,
                )
                .unwrap();
            assert_eq!(b.inflated as usize, block.len());
            let mut n = 0;
            decode_block(version, &block, |_, _, _| n += 1).unwrap();
            assert_eq!(b.records, n);
            assert_eq!(m.block(b.number), Some(b));
        }
        assert_eq!(sidecar.matrix(m.chr1, m.chr2, m.unit, m.bin_size), Some(m));
    }
}

fn check_encoding(bytes: &[u8], sidecar: &HicSidecar) {
    let identity = SourceIdentity::of(bytes).unwrap();
    let sealed = sidecar.to_bytes(identity);
    assert_eq!(&HicSidecar::parse(&sealed, identity).unwrap(), sidecar);
    let other = SourceIdentity {
        len: identity.len + 1,
        ..identity
    };
    assert_eq!(HicSidecar::parse(&sealed, other), Err(Error::StaleSidecar));
    let mut damaged = sealed.clone();
    damaged[60] ^= 1;
    assert!(HicSidecar::parse(&damaged, identity).is_err());
    assert!(HicSidecar::parse(&sealed[..sealed.len() - 1], identity).is_err());
}

fn from_test_writer(seed: u64, version: i32) {
    let mut rng = Rng(seed);
    let spec = random_spec(&mut rng, version);
    let bytes = common::write(&spec, &mut rng);
    let sidecar = HicSidecar::build(bytes.as_slice()).unwrap();
    check_against_file(&bytes, &sidecar);
    check_encoding(&bytes, &sidecar);

    let only_bp =
        HicSidecar::build_filtered(bytes.as_slice(), |_, _, unit, _| unit == ContactUnit::Bp)
            .unwrap();
    let want: Vec<_> = sidecar
        .matrices
        .iter()
        .filter(|m| m.unit == ContactUnit::Bp)
        .cloned()
        .collect();
    assert_eq!(only_bp.matrices, want);
}

fn random_header(rng: &mut Rng, version: i32) -> HicHeader {
    let chroms: Vec<Chrom> = (0..rng.range(1, 3))
        .map(|i| Chrom {
            name: format!("chr{}", i + 1),
            len: rng.range(20_000, 300_000),
        })
        .collect();
    let mut header = HicHeader::new("hg38", chroms, vec![10_000, 1_000, 200]);
    header.version = version;
    if version == 9 && rng.chance(50) {
        header
            .attributes
            .push(("v9-depth-base".into(), rng.range(2, 4).to_string()));
    }
    header
}

fn from_builder(seed: u64, version: i32) {
    let mut rng = Rng(seed);
    let header = random_header(&mut rng, version);
    let mut builder = HicBuilder::new(header.clone()).unwrap();
    let n = header.chroms.len() as u64;
    for _ in 0..rng.range(200, 8_000) {
        let (c1, c2) = (rng.range(0, n - 1) as u32, rng.range(0, n - 1) as u32);
        let len = |c: u32| header.chroms[c as usize].len;
        let p1 = rng.range(0, len(c1)) as u32;
        let p2 = if c1 == c2 && rng.chance(80) {
            {
                let reach = 1u64 << rng.range(4, 15);
                (p1 as u64 + rng.range(0, reach)).min(len(c2)) as u32
            }
        } else {
            rng.range(0, len(c2)) as u32
        };
        // Fractional counts force float blocks, whole ones allow short.
        let count = if rng.chance(30) {
            rng.range(1, 9) as f32 / 4.0
        } else {
            rng.range(1, 5) as f32
        };
        builder.add(c1, p1, c2, p2, count).unwrap();
    }
    let options = BuildOptions {
        level: rng.range(1, 9) as u32,
        ..Default::default()
    };
    let (out, sidecar) = builder
        .finish_with_sidecar(Cursor::new(Vec::new()), &options)
        .unwrap();
    let bytes = out.into_inner();
    assert_eq!(HicSidecar::build(bytes.as_slice()).unwrap(), sidecar);
    check_against_file(&bytes, &sidecar);
    check_encoding(&bytes, &sidecar);
}

proptest! {
    #![proptest_config(ProptestConfig { cases: 24, failure_persistence: Some(Box::new(proptest::test_runner::FileFailurePersistence::Direct(concat!(env!("CARGO_MANIFEST_DIR"), "/tests/cases/sidecar/encoding.proptest-regressions")))), ..ProptestConfig::default() })]

    #[test]
    fn built_from_files_of_every_block_layout(seed in any::<u64>(), version in 8i32..=9) {
        from_test_writer(seed, version);
    }

    #[test]
    fn writers_sidecars_match_built_ones(seed in any::<u64>(), v9 in any::<bool>()) {
        from_builder(seed, if v9 { 9 } else { 8 });
    }
}

#[test]
fn writer_sidecar_covers_dense_and_empty_maps() {
    let chroms = vec![
        Chrom {
            name: "chr1".into(),
            len: 50_000,
        },
        Chrom {
            name: "chr2".into(),
            len: 30_000,
        },
    ];
    let header = HicHeader::new("t", chroms, vec![1_000]);
    let mut w = HicWriter::create(Cursor::new(Vec::new()), header, 6).unwrap();
    // A full square: stored densely.
    let dense: Vec<Record> = (0..40u32)
        .flat_map(|x| {
            (x..40).map(move |y| Record {
                x,
                y,
                value: (x + y + 1) as f32,
            })
        })
        .collect();
    w.write_matrix(0, 0, ContactUnit::Bp, 1_000, &dense)
        .unwrap();
    w.write_matrix(0, 1, ContactUnit::Bp, 1_000, &[]).unwrap();
    w.write_matrix(
        1,
        1,
        ContactUnit::Bp,
        1_000,
        &[Record {
            x: 3,
            y: 7,
            value: 0.5,
        }],
    )
    .unwrap();
    let (out, sidecar) = w.finish_with_sidecar().unwrap();
    let bytes = out.into_inner();
    assert_eq!(HicSidecar::build(bytes.as_slice()).unwrap(), sidecar);
    check_against_file(&bytes, &sidecar);
    let empty = sidecar.matrix(0, 1, ContactUnit::Bp, 1_000).unwrap();
    assert!(empty.blocks.is_empty());
    let full = sidecar.matrix(0, 0, ContactUnit::Bp, 1_000).unwrap();
    assert_eq!(
        full.blocks
            .iter()
            .map(|b| b.records as usize)
            .sum::<usize>(),
        dense.len()
    );
}

#[test]
fn saved_beside_the_file_and_invalidated_by_changes() {
    let dir = tempfile::tempdir().unwrap();
    let path = dir.path().join("m.hic");
    let mut rng = Rng(7);
    let spec = random_spec(&mut rng, 8);
    std::fs::write(&path, common::write(&spec, &mut rng)).unwrap();
    let built = HicSidecar::build_file(&path).unwrap();
    assert!(sidecar_path(&path).exists());
    assert_eq!(HicSidecar::load(&path).unwrap(), built);
    built.save(&path).unwrap();
    assert_eq!(HicSidecar::load(&path).unwrap(), built);

    let mut bytes = std::fs::read(&path).unwrap();
    let last = bytes.len() - 1;
    bytes[last] ^= 0xff;
    std::fs::write(&path, &bytes).unwrap();
    assert_eq!(HicSidecar::load(&path), Err(Error::StaleSidecar));
}
