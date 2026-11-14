use std::collections::BTreeMap;
use std::io::Cursor;

use omicsweave_core::ContactUnit;
use omicsweave_hic::{
    BlockLayout, BuildOptions, Chrom, HicBuilder, HicHeader, HicSidecar, HicWriter, Record,
};

type Contact = (u32, u32, u32, u32, f32);

fn header(bp: Vec<u32>, frag: Vec<u32>) -> HicHeader {
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
fn oracle(h: &HicHeader, contacts: &[Contact]) -> (Vec<u8>, HicSidecar) {
    let mut writer =
        HicWriter::create_with(Cursor::new(Vec::new()), h.clone(), 1, BlockLayout::juicer())
            .unwrap();
    for unit in [ContactUnit::Bp, ContactUnit::Frag] {
        for &res in h.resolutions(unit) {
            let mut pairs: BTreeMap<_, Vec<Record>> = BTreeMap::new();
            for &(a, p, b, q, value) in contacts {
                let (a, p, b, q) = if a <= b { (a, p, b, q) } else { (b, q, a, p) };
                let bin = |c: u32, p: u32| {
                    let p = p.min(1000);
                    let pos = if unit == ContactUnit::Bp {
                        p
                    } else {
                        h.sites[c as usize].iter().filter(|&&s| s <= p).count() as u32
                    };
                    pos / res
                };
                let (mut x, mut y) = (bin(a, p), bin(b, q));
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
        let mut h = header(vec![7, 1], vec![2, 1]);
        h.version = version;
        let expected = oracle(&h, &contacts);
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
        let h = header(bp, frag);
        assert_eq!(build(h.clone(), &contacts), oracle(&h, &contacts));
    }
}
