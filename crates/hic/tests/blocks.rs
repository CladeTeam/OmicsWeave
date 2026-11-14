//! Block layouts: files cut into blocks of any size read back as the
//! contacts they were built from, by `HicReader` and by straw (which takes
//! each map's block size from the file), with matching sidecars.

mod common;

use std::collections::BTreeMap;
use std::io::Cursor;

use common::{Rng, straw};
use omicsweave_core::{ContactEntry, ContactQuery, ContactRead, ContactUnit, Error, Region};
use omicsweave_hic::{
    BlockBins, BlockLayout, BuildOptions, Chrom, DEFAULT_INTRA_BLOCK_BINS, DEPTH_BASE_ATTRIBUTE,
    Grid, HicBuilder, HicHeader, HicReader, HicSidecar, HicWriter, Record,
};
use proptest::prelude::*;

/// (row bin, col bin) -> value bits.
type Cells = BTreeMap<(u32, u32), u32>;
/// (chr1, chr2, res) -> stored (x, y) -> summed count.
type Maps = BTreeMap<(u32, u32, u32), BTreeMap<(u32, u32), f32>>;

fn random_blocks(rng: &mut Rng) -> BlockBins {
    match rng.range(0, 5) {
        0 => BlockBins::Juicer,
        1 => BlockBins::Fixed(1),
        2 => BlockBins::Fixed(rng.range(2, 16) as u32),
        3 => BlockBins::Fixed(DEFAULT_INTRA_BLOCK_BINS),
        _ => BlockBins::Fixed(rng.range(17, 400) as u32),
    }
}

fn build(rng: &mut Rng, version: i32) -> (HicHeader, Maps, BlockLayout, Vec<u8>, HicSidecar) {
    let chroms: Vec<Chrom> = (0..rng.range(1, 3))
        .map(|i| Chrom {
            name: format!("chr{}", i + 1),
            len: rng.range(20_000, 120_000),
        })
        .collect();
    let mut header = HicHeader::new("t", chroms, vec![2_000, 500]);
    header.version = version;
    if version == 9 && rng.chance(40) {
        let base = ["3", "-2", "-5"][rng.range(0, 2) as usize];
        header
            .attributes
            .push((DEPTH_BASE_ATTRIBUTE.into(), base.into()));
    }
    let layout = BlockLayout {
        intra: random_blocks(rng),
        inter: random_blocks(rng),
    };
    let mut builder = HicBuilder::new(header.clone()).unwrap();
    let mut maps = Maps::new();
    let n = header.chroms.len() as u64;
    for _ in 0..rng.range(100, 4_000) {
        let (a, b) = (rng.range(0, n - 1) as u32, rng.range(0, n - 1) as u32);
        let len = |c: u32| header.chroms[c as usize].len;
        let pa = rng.range(0, len(a) - 1) as u32;
        let pb = if a == b && rng.chance(80) {
            let reach = 1u64 << rng.range(3, 15);
            (pa as u64 + rng.range(0, reach)).min(len(b) - 1) as u32
        } else {
            rng.range(0, len(b) - 1) as u32
        };
        let count = rng.range(1, 6) as f32;
        builder.add(a, pa, b, pb, count).unwrap();
        let (c1, p1, c2, p2) = if a <= b {
            (a, pa, b, pb)
        } else {
            (b, pb, a, pa)
        };
        for &res in &header.bp_resolutions {
            let (mut x, mut y) = (p1 / res, p2 / res);
            if c1 == c2 && x > y {
                std::mem::swap(&mut x, &mut y);
            }
            *maps
                .entry((c1, c2, res))
                .or_default()
                .entry((x, y))
                .or_default() += count;
        }
    }
    let options = BuildOptions {
        level: 1,
        blocks: layout,
    };
    let (out, sidecar) = builder
        .finish_with_sidecar(Cursor::new(Vec::new()), &options)
        .unwrap();
    (header, maps, layout, out.into_inner(), sidecar)
}

/// A tile: rows on chromosome `r`, columns on `c`, as bin ranges.
#[derive(Debug, Clone, Copy)]
struct Tile {
    r: u32,
    c: u32,
    rows: (u32, u32),
    cols: (u32, u32),
}

/// Puts a stored record (x on the lower chromosome) into the tile,
/// mirrored within a chromosome.
fn place(t: &Tile, x: u32, y: u32, bits: u32, out: &mut Cells) {
    let candidates = if t.r == t.c {
        vec![(x, y), (y, x)]
    } else if t.r < t.c {
        vec![(x, y)]
    } else {
        vec![(y, x)]
    };
    for (row, col) in candidates {
        if (t.rows.0..t.rows.1).contains(&row) && (t.cols.0..t.cols.1).contains(&col) {
            out.insert((row, col), bits);
        }
    }
}

fn case(seed: u64, version: i32) {
    let mut rng = Rng(seed);
    let (header, maps, layout, bytes, sidecar) = build(&mut rng, version);
    let dir = tempfile::tempdir().unwrap();
    let path = dir.path().join("b.hic");
    std::fs::write(&path, &bytes).unwrap();
    assert_eq!(HicSidecar::build(bytes.as_slice()).unwrap(), sidecar);
    let mut reader = HicReader::open(&bytes[..]).unwrap();

    // Every stored map carries the grid the layout asked for.
    for &(c1, c2, res) in maps.keys() {
        let zoom = reader.zoom(c1, c2, ContactUnit::Bp, res).unwrap().unwrap();
        let extent = header.chroms[c1 as usize]
            .len
            .max(header.chroms[c2 as usize].len);
        let want = match if c1 == c2 { layout.intra } else { layout.inter } {
            BlockBins::Juicer => Grid::juicer(ContactUnit::Bp, res, extent, c1 == c2),
            BlockBins::Fixed(b) => Grid::with_block_bins(res, extent, b).unwrap(),
        };
        assert_eq!(
            (zoom.block_bin_count, zoom.block_column_count),
            (want.block_bins, want.columns)
        );
    }

    for _ in 0..30 {
        let res = header.bp_resolutions[rng.range(0, 1) as usize];
        let n = header.chroms.len() as u64 - 1;
        let (r, c) = (rng.range(0, n) as u32, rng.range(0, n) as u32);
        let span = |rng: &mut Rng, chrom: u32| {
            let total = header.bins(chrom, ContactUnit::Bp, res) as u32;
            let start = rng.range(0, total as u64 - 1) as u32;
            (start, (start + rng.range(1, 300) as u32).min(total))
        };
        let rows = span(&mut rng, r);
        let cols = if r == c && rng.chance(60) {
            rows
        } else {
            span(&mut rng, c)
        };
        let t = Tile { r, c, rows, cols };
        let (c1, c2) = (r.min(c), r.max(c));
        let mut want = Cells::new();
        if let Some(map) = maps.get(&(c1, c2, res)) {
            for (&(x, y), &v) in map {
                place(&t, x, y, v.to_bits(), &mut want);
            }
        }
        let q = ContactQuery::new(
            Region::new(&header.chroms[r as usize].name, rows.0 * res, rows.1 * res),
            Region::new(&header.chroms[c as usize].name, cols.0 * res, cols.1 * res),
            res,
        );
        let n = reader.query_entries(q).unwrap();
        let mut entries = vec![ContactEntry::default(); n];
        reader.write_entries(&mut entries).unwrap();
        let ours: Cells = entries
            .iter()
            .map(|e| ((e.row / res, e.col / res), e.value.to_bits()))
            .collect();
        let what = format!("v{version} {layout:?} res {res} {r}x{c} {rows:?} {cols:?}");
        assert_eq!(ours, want, "{what}: reader");

        // straw computes v9 intra block numbers with a base-2 log depth
        // whatever `v9-depth-base` says, so it reads only that layout.
        if version == 9
            && header
                .attributes
                .iter()
                .any(|(k, _)| k == DEPTH_BASE_ATTRIBUTE)
        {
            continue;
        }
        let loc = |chrom: u32, b: (u32, u32)| {
            format!(
                "{}:{}:{}",
                header.chroms[chrom as usize].name,
                b.0 * res,
                (b.1 - 1) * res
            )
        };
        let mut theirs = Cells::new();
        for (x, y, v) in straw(&path, &loc(r, rows), &loc(c, cols), ContactUnit::Bp, res) {
            place(&t, x / res, y / res, v.to_bits(), &mut theirs);
        }
        assert_eq!(theirs, want, "{what}: straw");
    }
}

proptest! {
    #![proptest_config(ProptestConfig { cases: 24, ..ProptestConfig::default() })]

    #[test]
    fn any_block_size_reads_back_v9(seed in any::<u64>()) {
        case(seed, 9);
    }

    #[test]
    fn any_block_size_reads_back_v8(seed in any::<u64>()) {
        case(seed, 8);
    }
}

#[test]
fn default_layout_uses_small_blocks_within_chromosomes() {
    let chroms = vec![
        Chrom {
            name: "a".into(),
            len: 2_000_000,
        },
        Chrom {
            name: "b".into(),
            len: 1_000_000,
        },
    ];
    let header = HicHeader::new("t", chroms, vec![1_000]);
    let one = [Record {
        x: 3,
        y: 9,
        value: 2.0,
    }];
    let mut w = HicWriter::create(Cursor::new(Vec::new()), header.clone(), 1).unwrap();
    assert_eq!(w.layout(), BlockLayout::default());
    w.write_matrix(0, 0, ContactUnit::Bp, 1_000, &one).unwrap();
    w.write_matrix(0, 1, ContactUnit::Bp, 1_000, &one).unwrap();
    // One map's own block size overrides the layout.
    w.write_matrix_with_blocks(1, 1, ContactUnit::Bp, 1_000, &one, BlockBins::Fixed(32))
        .unwrap();
    let bytes = w.finish().unwrap().into_inner();
    let mut reader = HicReader::open(&bytes[..]).unwrap();
    let grid = |r: &mut HicReader<&[u8]>, a, b| {
        let z = r.zoom(a, b, ContactUnit::Bp, 1_000).unwrap().unwrap();
        (z.block_bin_count, z.block_column_count)
    };
    assert_eq!(
        grid(&mut reader, 0, 0),
        (DEFAULT_INTRA_BLOCK_BINS, 2_001 / 128 + 1)
    );
    let juicer = Grid::juicer(ContactUnit::Bp, 1_000, 2_000_000, false);
    assert_eq!(grid(&mut reader, 0, 1), (juicer.block_bins, juicer.columns));
    assert_eq!(grid(&mut reader, 1, 1), (32, 1_001 / 32 + 1));
}

#[test]
fn grids_too_fine_for_block_numbers_are_rejected() {
    let header = HicHeader::new(
        "t",
        vec![Chrom {
            name: "a".into(),
            len: 1 << 30,
        }],
        vec![10],
    );
    let mut w = HicWriter::create_with(
        Cursor::new(Vec::new()),
        header.clone(),
        1,
        BlockLayout {
            intra: BlockBins::Fixed(1),
            inter: BlockBins::Juicer,
        },
    )
    .unwrap();
    let one = [Record {
        x: 0,
        y: 0,
        value: 1.0,
    }];
    assert!(matches!(
        w.write_matrix(0, 0, ContactUnit::Bp, 10, &one),
        Err(Error::Invalid { .. })
    ));
    // The writer stays usable: the map was rejected before any output.
    w.write_matrix_with_blocks(0, 0, ContactUnit::Bp, 10, &one, BlockBins::Fixed(4096))
        .unwrap();
    assert!(w.finish().is_ok());

    let zero = BlockLayout {
        intra: BlockBins::Fixed(0),
        inter: BlockBins::Juicer,
    };
    assert!(matches!(
        HicWriter::create_with(Cursor::new(Vec::new()), header, 1, zero),
        Err(Error::Invalid { .. })
    ));
}
