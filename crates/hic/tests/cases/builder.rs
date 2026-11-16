use crate::support::oracle::{Contact, binned_file};
use omicsweave_hic::{
    BlockLayout, BuildOptions, Chrom, HicBuilder, HicHeader, HicSidecar, HicWriter,
};
use std::io::Cursor;

fn resolution_header(version: i32, resolutions: usize) -> HicHeader {
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
            let header = resolution_header(version, count);
            let (want, expected_sidecar) = binned_file(&header, &contacts, blocks);
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
fn empty_builder_preserves_declared_resolutions_without_maps() {
    for count in [0, 1, 8, 9, 64] {
        let header = resolution_header(9, count);
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
fn fragment_header(bp: Vec<u32>, frag: Vec<u32>) -> HicHeader {
    let mut h = HicHeader::new(
        "oracle",
        vec![
            Chrom {
                name: "a".into(),
                len: 1000,
            },
            Chrom {
                name: "b".into(),
                len: 1000,
            },
        ],
        bp,
    );
    h.frag_resolutions = frag;
    h.sites = vec![vec![0, 17, 17, 100, 700, 1000]; 2];
    h
}

// Each resolution is independently binned from the input. Public writer
// sorting supplies the cell sums without the builder's accumulators.
fn build(h: HicHeader, contacts: &[Contact]) -> (Vec<u8>, HicSidecar) {
    let mut builder = HicBuilder::new(h).unwrap();
    // Rejected contacts must not instantiate a chromosome pair or mutate
    // one of its resolution accumulators.
    for bad in [(2, 0, 0, 0, 1.0), (0, 0, 0, 0, f32::INFINITY)] {
        assert!(builder.add(bad.0, bad.1, bad.2, bad.3, bad.4).is_err());
    }
    for &(a, p, b, q, v) in contacts {
        builder.add(a, p, b, q, v).unwrap();
    }
    let (out, sidecar) = builder
        .finish_with_sidecar(
            Cursor::new(Vec::new()),
            &BuildOptions {
                level: 1,
                blocks: BlockLayout::juicer(),
            },
        )
        .unwrap();
    (out.into_inner(), sidecar)
}

#[test]
fn pair_accumulators_preserve_float_order_across_compactions() {
    let mut contacts: Vec<_> = (0..(2 * (1 << 16) + 7))
        .map(|i| {
            let value = [16_777_216.0, 1.0, -16_777_216.0, 0.5][i % 4];
            if i % 2 == 0 {
                (1, 17, 0, 700, value)
            } else {
                (0, 700, 1, 17, value)
            }
        })
        .collect();
    contacts.extend([
        (1, u32::MAX, 1, 17, 3.5),
        (0, 1000, 0, 0, -2.0),
        (0, 0, 0, 1000, 0.0),
    ]);
    for version in [8, 9] {
        let mut h = fragment_header(vec![7, 1], vec![2, 1]);
        h.version = version;
        let expected = binned_file(&h, &contacts, BlockLayout::juicer());
        let actual = build(h, &contacts);
        assert_eq!(actual, expected);
        assert_eq!(actual.1, HicSidecar::build(actual.0.as_slice()).unwrap());
    }
}

#[test]
fn pair_accumulators_handle_independent_unit_lists_and_zero_only_pairs() {
    let contacts = [
        (1, 9, 0, 17, 0.0),
        (1, 1001, 1, 0, 1.5),
        (0, 99, 0, 100, 2.0),
    ];
    for (bp, frag) in [
        (vec![], vec![]),
        (vec![], vec![4, 1, 2]),
        (vec![9, 1, 7], vec![]),
        ((1..=9).rev().collect(), vec![3, 1]),
    ] {
        let h = fragment_header(bp, frag);
        assert_eq!(
            build(h.clone(), &contacts),
            binned_file(&h, &contacts, BlockLayout::juicer())
        );
    }
}
