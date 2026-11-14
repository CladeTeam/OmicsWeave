//! Reading: tiles of files written by an independent test writer, against
//! straw and against the records the files were written from, for v8 and
//! v9 files, base-pair and fragment units. The files also hold
//! normalization and expected vectors, which the reader skips.

mod common;

use std::collections::BTreeMap;

use common::{Rng, Spec, random_spec, straw, write};
use omicsweave_core::{ContactEntry, ContactQuery, ContactRead, ContactUnit, Error, Region};
use omicsweave_hic::HicReader;
use proptest::prelude::*;

/// A tile as bins.
#[derive(Debug, Clone, Copy)]
struct Tile {
    unit: ContactUnit,
    r: u32,
    rows: (u32, u32),
    c: u32,
    cols: (u32, u32),
    res: u32,
}

/// Cells as (row bin, col bin) -> value bits.
type Cells = BTreeMap<(u32, u32), u32>;

fn random_tile(spec: &Spec, rng: &mut Rng) -> Tile {
    let units = spec.units();
    let unit = units[rng.range(0, units.len() as u64 - 1) as usize];
    let list = spec.resolutions(unit);
    let res = list[rng.range(0, list.len() as u64 - 1) as usize];
    let n = spec.chroms.len() as u64 - 1;
    let (r, c) = (rng.range(0, n) as u32, rng.range(0, n) as u32);
    let span = |rng: &mut Rng, chrom: u32| {
        let total = spec.bins(chrom, unit, res);
        let start = rng.range(0, total as u64 - 1) as u32;
        (start, (start + rng.range(1, 120) as u32).min(total))
    };
    let rows = span(rng, r);
    let cols = if r == c && rng.chance(60) {
        // Tiles on or near the diagonal, either side.
        let total = spec.bins(c, unit, res);
        let start = rows
            .0
            .saturating_sub(rng.range(0, 40) as u32)
            .min(total - 1);
        (start, (start + rng.range(1, 120) as u32).min(total))
    } else {
        span(rng, c)
    };
    Tile {
        unit,
        r,
        rows,
        c,
        cols,
        res,
    }
}

fn query<'a>(spec: &'a Spec, t: &Tile) -> ContactQuery<'a> {
    let res = t.res;
    ContactQuery::new(
        Region::new(&spec.chroms[t.r as usize].0, t.rows.0 * res, t.rows.1 * res),
        Region::new(&spec.chroms[t.c as usize].0, t.cols.0 * res, t.cols.1 * res),
        res,
    )
    .with_unit(t.unit)
}

/// The tile cells of a stored record (x on chr1, y on chr2), mirrored
/// within a chromosome.
fn cells_of(t: &Tile, x: u32, y: u32) -> Vec<(u32, u32)> {
    let inside = |row: u32, col: u32| {
        (t.rows.0..t.rows.1).contains(&row) && (t.cols.0..t.cols.1).contains(&col)
    };
    let candidates = if t.r == t.c {
        vec![(x, y), (y, x)]
    } else if t.r < t.c {
        vec![(x, y)]
    } else {
        vec![(y, x)]
    };
    candidates
        .into_iter()
        .filter(|&(r, c)| inside(r, c))
        .collect()
}

/// Counts of the records the file was written from.
struct Oracle<'a> {
    spec: &'a Spec,
    t: Tile,
}

impl Oracle<'_> {
    fn chroms(&self) -> (u32, u32) {
        (self.t.r.min(self.t.c), self.t.r.max(self.t.c))
    }

    fn cells(&self) -> Cells {
        let (c1, c2) = self.chroms();
        let t = self.t;
        let mut cells = Cells::new();
        for map in self
            .spec
            .maps
            .iter()
            .filter(|m| m.chr1 == c1 && m.chr2 == c2 && m.unit == t.unit && m.res == t.res)
        {
            for &(x, y, v) in &map.records {
                if v != 0.0 {
                    for cell in cells_of(&t, x, y) {
                        cells.insert(cell, v.to_bits());
                    }
                }
            }
        }
        cells
    }

    fn check_matrix(&self, cells: &Cells, matrix: &[f32]) {
        let t = self.t;
        let width = (t.cols.1 - t.cols.0) as usize;
        for row in t.rows.0..t.rows.1 {
            for col in t.cols.0..t.cols.1 {
                let got = matrix[(row - t.rows.0) as usize * width + (col - t.cols.0) as usize];
                let want = cells.get(&(row, col)).map_or(0.0, |&b| f32::from_bits(b));
                assert_eq!(got.to_bits(), want.to_bits(), "cell ({row}, {col})");
            }
        }
    }

    /// straw's raw observed cells.
    fn straw(&self, path: &std::path::Path) -> Cells {
        let (spec, t) = (self.spec, self.t);
        let loc = |chrom: u32, b: (u32, u32)| {
            format!(
                "{}:{}:{}",
                spec.chroms[chrom as usize].0,
                b.0 * t.res,
                (b.1 - 1) * t.res
            )
        };
        let records = straw(path, &loc(t.r, t.rows), &loc(t.c, t.cols), t.unit, t.res);
        let mut cells = Cells::new();
        for (x, y, v) in records {
            for cell in cells_of(&t, x / t.res, y / t.res) {
                cells.insert(cell, v.to_bits());
            }
        }
        cells
    }
}

fn ours(reader: &mut HicReader<&[u8]>, q: ContactQuery<'_>) -> (Cells, Vec<f32>) {
    let n = reader.query_entries(q).unwrap();
    let mut entries = vec![ContactEntry::default(); n];
    reader.write_entries(&mut entries).unwrap();
    for pair in entries.windows(2) {
        assert!(
            (pair[0].row, pair[0].col) < (pair[1].row, pair[1].col),
            "entries sorted and unique"
        );
    }
    let cells = entries
        .iter()
        .map(|e| {
            (
                (e.row / q.resolution, e.col / q.resolution),
                e.value.to_bits(),
            )
        })
        .collect();
    let (rows, cols) = q.shape();
    let mut matrix = vec![1.5f32; rows * cols];
    reader.matrix(q, &mut matrix).unwrap();
    (cells, matrix)
}

fn run(seed: u64, version: i32) {
    let mut rng = Rng(seed);
    let spec = random_spec(&mut rng, version);
    let bytes = write(&spec, &mut rng);
    let dir = tempfile::tempdir().unwrap();
    let path = dir.path().join("t.hic");
    std::fs::write(&path, &bytes).unwrap();
    let mut reader = HicReader::open(&bytes[..]).unwrap();
    assert_eq!(
        reader.resolutions(ContactUnit::Bp),
        &spec.bp_resolutions[..]
    );
    assert_eq!(
        reader.resolutions(ContactUnit::Frag),
        &spec.frag_resolutions[..]
    );
    for (c, sites) in spec.sites.iter().enumerate() {
        assert_eq!(reader.sites(c as u32).unwrap(), &sites[..]);
    }

    for _ in 0..48 {
        let t = random_tile(&spec, &mut rng);
        let oracle = Oracle { spec: &spec, t };
        let (cells, matrix) = ours(&mut reader, query(&spec, &t));
        let what = format!("v{version} {t:?}");
        assert_eq!(cells, oracle.cells(), "{what}: written records");
        assert_eq!(cells, oracle.straw(&path), "{what}: straw");
        oracle.check_matrix(&cells, &matrix);
    }
}

proptest! {
    #![proptest_config(ProptestConfig { cases: 24, ..ProptestConfig::default() })]

    #[test]
    fn v8_matches_straw(seed in any::<u64>()) {
        run(seed, 8);
    }

    #[test]
    fn v9_matches_straw(seed in any::<u64>()) {
        run(seed, 9);
    }
}

#[test]
fn highly_compressed_blocks_decode_whole() {
    // straw inflates into 10x the compressed size; this reader has no limit.
    let mut rng = Rng(5);
    let mut spec = random_spec(&mut rng, 9);
    spec.level = 9;
    spec.short_values = true;
    for map in &mut spec.maps {
        for r in &mut map.records {
            r.2 = 7.0;
        }
    }
    let bytes = write(&spec, &mut rng);
    let mut reader = HicReader::open(&bytes[..]).unwrap();
    for _ in 0..40 {
        let t = random_tile(&spec, &mut rng);
        let (cells, _) = ours(&mut reader, query(&spec, &t));
        assert_eq!(cells, Oracle { spec: &spec, t }.cells());
    }
}

#[test]
fn rejects_bad_requests_and_input() {
    let mut rng = Rng(1);
    let spec = random_spec(&mut rng, 9);
    let bytes = write(&spec, &mut rng);
    let mut reader = HicReader::open(&bytes[..]).unwrap();
    let chrom = spec.chroms[0].0.as_str();
    let region = Region::new(chrom, 0, 10_000);
    let q = |res| ContactQuery::new(region, region, res);

    let mut none: [ContactEntry; 0] = [];
    assert_eq!(reader.write_entries(&mut none), Err(Error::NotQueried));
    assert_eq!(
        reader.chrom_len("chrUn"),
        Err(Error::UnknownChrom("chrUn".into()))
    );
    let invalid = |r: omicsweave_core::Result<usize>| matches!(r, Err(Error::Invalid { .. }));
    assert!(invalid(reader.query_entries(q(2_000))));
    if spec.frag_resolutions.is_empty() {
        assert!(invalid(
            reader.query_entries(q(1).with_unit(ContactUnit::Frag))
        ));
    }
    let mut small = vec![0f32; 3];
    assert!(matches!(
        reader.matrix(q(1_000), &mut small),
        Err(Error::OutputSize { .. })
    ));
    assert!(matches!(
        reader.query_entries(ContactQuery::new(
            Region::new("chrUn", 0, 10),
            region,
            1_000,
        )),
        Err(Error::UnknownChrom(_))
    ));
    // A failed query leaves nothing staged.
    assert_eq!(reader.write_entries(&mut none), Err(Error::NotQueried));

    assert!(HicReader::open(&bytes[..200]).is_err());
    let mut bad = bytes.clone();
    bad[0] = b'X';
    assert!(matches!(
        HicReader::open(&bad[..]),
        Err(Error::BadMagic { .. })
    ));
}

#[test]
fn versions_other_than_8_and_9_are_unsupported() {
    let mut rng = Rng(3);
    let spec = random_spec(&mut rng, 9);
    let mut bytes = write(&spec, &mut rng);
    for version in [5i32, 6, 7, 10] {
        bytes[4..8].copy_from_slice(&version.to_le_bytes());
        assert!(
            matches!(HicReader::open(&bytes[..]), Err(Error::Unsupported(_))),
            "v{version}"
        );
    }
}
