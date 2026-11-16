//! cool / mcool → .hic: the file and sidecar `HicWriter::write_matrix`
//! writes for the same pixels, whether the sort spilled or not.

use crate::cool_support as coolfix;
use crate::support as common;

use std::io::Cursor;

use common::rng::Rng;
use coolfix::{MapSpec, sample_pixels, write_cool, write_mcool};
use omicsweave_convert::{Budget, Conversion, CoolToHic, convert};
use omicsweave_core::ContactUnit;
use omicsweave_hic::{Chrom, HicHeader, HicSidecar, HicWriter, Record};

const CHROMS: [(&str, u32); 3] = [("chr1", 900_000), ("chr2", 610_000), ("chrX", 333_333)];

fn spec(bin: u32, n: usize, seed: u64) -> MapSpec {
    let mut s = MapSpec::new(&CHROMS, bin);
    s.pixels = sample_pixels(s.n_bins(), n, seed, true);
    s
}

/// The `.hic` of maps written whole, in the converter's order: resolutions
/// finest first, chromosome pairs ascending.
fn expected(specs: &[MapSpec]) -> (Vec<u8>, HicSidecar) {
    let mut res: Vec<u32> = specs.iter().map(|s| s.bin_size).collect();
    res.reverse();
    let header = HicHeader::new(
        "test",
        CHROMS
            .iter()
            .map(|&(n, l)| Chrom {
                name: n.into(),
                len: l as u64,
            })
            .collect(),
        res,
    );
    let mut w = HicWriter::create(Cursor::new(Vec::new()), header, 6).unwrap();
    for s in specs {
        let off = s.chrom_offsets();
        let chrom = |b: u64| off.partition_point(|&o| o <= b) as u32 - 1;
        let mut maps: std::collections::BTreeMap<(u32, u32), Vec<Record>> = Default::default();
        for &(b1, b2, c) in &s.pixels {
            let (c1, c2) = (chrom(b1), chrom(b2));
            maps.entry((c1, c2)).or_default().push(Record {
                x: (b1 - off[c1 as usize]) as u32,
                y: (b2 - off[c2 as usize]) as u32,
                value: c as f32,
            });
        }
        for ((c1, c2), recs) in maps {
            w.write_matrix(c1, c2, ContactUnit::Bp, s.bin_size, &recs)
                .unwrap();
        }
    }
    let (out, sidecar) = w.finish_with_sidecar().unwrap();
    (out.into_inner(), sidecar)
}

#[test]
fn hic_matches_the_writer() {
    let dir = tempfile::tempdir().unwrap();
    let mut rng = Rng::new(7);
    let specs = vec![
        spec(1_000, 200_000, rng.next()),
        spec(5_000, 30_000, rng.next()),
    ];
    let cases = [
        ("one.cool", specs[..1].to_vec(), None),
        ("two.mcool", specs.clone(), None),
        ("pick.mcool", specs.clone(), Some(vec![5_000])),
    ];
    for (name, maps, pick) in cases {
        let src = dir.path().join(name);
        if name.ends_with(".mcool") {
            write_mcool(&src, &maps);
        } else {
            write_cool(&src, &maps[0]);
        }
        let chosen: Vec<MapSpec> = match &pick {
            Some(p) => maps
                .iter()
                .filter(|m| p.contains(&m.bin_size))
                .cloned()
                .collect(),
            None => maps.clone(),
        };
        let (want, want_sidecar) = expected(&chosen);
        for budget in [
            Budget::new(1, 1 << 30),
            Budget::new(6, 1 << 30),
            Budget::new(3, 8 << 20),
        ] {
            let dst = dir.path().join("out.hic");
            let mut conversion = CoolToHic::new(&src, &dst);
            conversion.resolutions = pick.clone();
            let report = convert(Conversion::CoolToHic(conversion), budget.clone(), None).unwrap();
            let got = std::fs::read(&dst).unwrap();
            let first_diff = got.iter().zip(&want).position(|(a, b)| a != b);
            assert!(
                got == want,
                "{name} {budget:?}: {} vs {} bytes, first difference at {first_diff:?}",
                got.len(),
                want.len()
            );
            assert_eq!(HicSidecar::load(&dst).unwrap(), want_sidecar);
            let pixels: usize = chosen.iter().map(|m| m.pixels.len()).sum();
            assert_eq!(report.records, pixels as u64);
            if budget.memory <= 8 << 20 && chosen.iter().any(|m| m.pixels.len() > 100_000) {
                assert!(report.spilled > 0, "a tight budget spills");
            }
        }
    }
}
