//! Writing: files from `HicBuilder` and `HicWriter`, read back by straw and
//! by `HicReader`, against the contacts and records they were written from.

mod common;

use std::collections::BTreeMap;
use std::io::Cursor;

use common::{Rng, straw};
use omicsweave_core::{ContactEntry, ContactQuery, ContactRead, ContactUnit, Error, Region};
use omicsweave_hic::{
    BuildOptions, Chrom, DEPTH_BASE_ATTRIBUTE, HicBuilder, HicHeader, HicReader, HicWriter, Record,
};
use proptest::prelude::*;

/// Cells as (row bin, col bin) -> value bits.
type Cells = BTreeMap<(u32, u32), u32>;

/// A tile as bins: unit, resolution, rows chromosome and bins, columns
/// chromosome and bins.
#[derive(Debug, Clone, Copy)]
struct Tile {
    unit: ContactUnit,
    res: u32,
    r: u32,
    rows: (u32, u32),
    c: u32,
    cols: (u32, u32),
}

fn random_header(rng: &mut Rng, version: i32) -> HicHeader {
    let chroms: Vec<Chrom> = (0..rng.range(1, 3))
        .map(|i| Chrom {
            name: format!("chr{}", i + 1),
            len: rng.range(20_000, 200_000),
        })
        .collect();
    let mut header = HicHeader::new("hg38", chroms, vec![5_000, 1_000]);
    header.version = version;
    if rng.chance(50) {
        header.frag_resolutions = vec![5, 1];
        for c in &header.chroms {
            let mut s = Vec::new();
            for _ in 0..rng.range(30, 400) {
                s.push(rng.range(1, c.len - 1) as u32);
            }
            s.sort_unstable();
            s.dedup();
            header.sites.push(s);
        }
    }
    header
}

/// Contacts as (chr1, pos1, chr2, pos2, count); counts are multiples of
/// 1/2 so that sums are exact in any order.
fn random_contacts(rng: &mut Rng, header: &HicHeader) -> Vec<(u32, u32, u32, u32, f32)> {
    let n = header.chroms.len() as u64;
    (0..rng.range(500, 6_000))
        .map(|_| {
            let (c1, c2) = (rng.range(0, n - 1) as u32, rng.range(0, n - 1) as u32);
            let len = |c: u32| header.chroms[c as usize].len;
            let p1 = rng.range(0, len(c1) + 10) as u32;
            let p2 = if c1 == c2 && rng.chance(80) {
                let reach = 1u64 << rng.range(4, 16);
                (p1 as u64 + rng.range(0, reach)).min(len(c2)) as u32
            } else {
                rng.range(0, len(c2)) as u32
            };
            (c1, p1, c2, p2, rng.range(1, 6) as f32 / 2.0)
        })
        .collect()
}

/// Stored records of each map: (chr1, chr2, unit, res) -> (x, y) -> count.
type Maps = BTreeMap<(u32, u32, ContactUnit, u32), BTreeMap<(u32, u32), f32>>;

fn bin_contacts(header: &HicHeader, contacts: &[(u32, u32, u32, u32, f32)]) -> Maps {
    let mut maps = Maps::new();
    for &(a, pa, b, pb, v) in contacts {
        let (c1, p1, c2, p2) = if a <= b {
            (a, pa, b, pb)
        } else {
            (b, pb, a, pa)
        };
        let clamp = |c: u32, p: u32| p.min(header.chroms[c as usize].len as u32);
        let (p1, p2) = (clamp(c1, p1), clamp(c2, p2));
        let mut put = |unit, res: u32, p: u32, q: u32| {
            let (mut x, mut y) = (p / res, q / res);
            if c1 == c2 && x > y {
                std::mem::swap(&mut x, &mut y);
            }
            *maps
                .entry((c1, c2, unit, res))
                .or_default()
                .entry((x, y))
                .or_default() += v;
        };
        for &res in &header.bp_resolutions {
            put(ContactUnit::Bp, res, p1, p2);
        }
        for &res in &header.frag_resolutions {
            let frag =
                |c: u32, p: u32| header.sites[c as usize].partition_point(|&s| s <= p) as u32;
            put(ContactUnit::Frag, res, frag(c1, p1), frag(c2, p2));
        }
    }
    maps
}

fn units(header: &HicHeader) -> Vec<ContactUnit> {
    if header.frag_resolutions.is_empty() {
        vec![ContactUnit::Bp]
    } else {
        vec![ContactUnit::Bp, ContactUnit::Frag]
    }
}

fn random_tile(header: &HicHeader, rng: &mut Rng) -> Tile {
    let units = units(header);
    let unit = units[rng.range(0, units.len() as u64 - 1) as usize];
    let list = header.resolutions(unit);
    let res = list[rng.range(0, list.len() as u64 - 1) as usize];
    let n = header.chroms.len() as u64 - 1;
    let (r, c) = (rng.range(0, n) as u32, rng.range(0, n) as u32);
    let span = |rng: &mut Rng, chrom: u32| {
        let total = header.bins(chrom, unit, res) as u32;
        let start = rng.range(0, total as u64 - 1) as u32;
        (start, (start + rng.range(1, 150) as u32).min(total))
    };
    let rows = span(rng, r);
    let cols = if r == c && rng.chance(50) {
        rows
    } else {
        span(rng, c)
    };
    Tile {
        unit,
        res,
        r,
        rows,
        c,
        cols,
    }
}

fn query<'a>(header: &'a HicHeader, t: &Tile) -> ContactQuery<'a> {
    ContactQuery::new(
        Region::new(
            &header.chroms[t.r as usize].name,
            t.rows.0 * t.res,
            t.rows.1 * t.res,
        ),
        Region::new(
            &header.chroms[t.c as usize].name,
            t.cols.0 * t.res,
            t.cols.1 * t.res,
        ),
        t.res,
    )
    .with_unit(t.unit)
}

/// Puts a stored record (x on the lower chromosome) into the tile,
/// mirrored within a chromosome.
fn place(t: &Tile, x: u32, y: u32, bits: u32, cells: &mut Cells) {
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
    for (row, col) in candidates {
        if inside(row, col) {
            cells.insert((row, col), bits);
        }
    }
}

fn entries<S: omicsweave_core::ReadAt>(reader: &mut HicReader<S>, q: ContactQuery<'_>) -> Cells {
    let n = reader.query_entries(q).unwrap();
    let mut out = vec![ContactEntry::default(); n];
    reader.write_entries(&mut out).unwrap();
    out.iter()
        .map(|e| {
            (
                (e.row / q.resolution, e.col / q.resolution),
                e.value.to_bits(),
            )
        })
        .collect()
}

fn straw_cells(path: &std::path::Path, header: &HicHeader, t: &Tile) -> Cells {
    let loc = |chrom: u32, b: (u32, u32)| {
        format!(
            "{}:{}:{}",
            header.chroms[chrom as usize].name,
            b.0 * t.res,
            (b.1 - 1) * t.res
        )
    };
    let mut cells = Cells::new();
    for (x, y, v) in straw(path, &loc(t.r, t.rows), &loc(t.c, t.cols), t.unit, t.res) {
        place(t, x / t.res, y / t.res, v.to_bits(), &mut cells);
    }
    cells
}

fn build(seed: u64, version: i32) {
    let mut rng = Rng(seed);
    let header = random_header(&mut rng, version);
    let contacts = random_contacts(&mut rng, &header);
    let mut builder = HicBuilder::new(header.clone()).unwrap();
    for &(c1, p1, c2, p2, v) in &contacts {
        builder.add(c1, p1, c2, p2, v).unwrap();
    }
    let bytes = builder
        .finish(Cursor::new(Vec::new()), &BuildOptions::default())
        .unwrap()
        .into_inner();
    let dir = tempfile::tempdir().unwrap();
    let path = dir.path().join("b.hic");
    std::fs::write(&path, &bytes).unwrap();
    let maps = bin_contacts(&header, &contacts);
    let mut reader = HicReader::open(&bytes[..]).unwrap();

    for _ in 0..40 {
        let t = random_tile(&header, &mut rng);
        // Observed counts against the contacts.
        let (c1, c2) = (t.r.min(t.c), t.r.max(t.c));
        let mut want = Cells::new();
        if let Some(map) = maps.get(&(c1, c2, t.unit, t.res)) {
            for (&(x, y), &v) in map {
                place(&t, x, y, v.to_bits(), &mut want);
            }
        }
        let ours = entries(&mut reader, query(&header, &t));
        assert_eq!(ours, want, "v{version} {t:?}: counts");
        assert_eq!(
            straw_cells(&path, &header, &t),
            want,
            "v{version} {t:?}: straw"
        );
    }
}

proptest! {
    #![proptest_config(ProptestConfig { cases: 12, ..ProptestConfig::default() })]

    #[test]
    fn built_v9_files_read_as_straw_reads_them(seed in any::<u64>()) {
        build(seed, 9);
    }

    #[test]
    fn built_v8_files_read_as_straw_reads_them(seed in any::<u64>()) {
        build(seed, 8);
    }
}

/// Random records of one map, in any orientation, with repeats and zeros.
fn random_records(rng: &mut Rng, n1: u32, n2: u32, intra: bool, short: bool) -> Vec<Record> {
    let mut records = Vec::new();
    for _ in 0..rng.range(1, 4_000) {
        let x = rng.range(0, n1 as u64 - 1) as u32;
        let y = if intra && rng.chance(70) {
            let reach = 1 << rng.range(0, 18);
            (x as u64 + rng.range(0, reach)).min(n2 as u64 - 1) as u32
        } else {
            rng.range(0, n2 as u64 - 1) as u32
        };
        let (x, y) = if intra && rng.chance(50) {
            (y, x)
        } else {
            (x, y)
        };
        let value = match rng.range(0, 20) {
            0 => 0.0,
            _ if short => rng.range(1, 300) as f32,
            _ => rng.range(1, 1 << 20) as f32 / 64.0,
        };
        records.push(Record { x, y, value });
    }
    records
}

/// What a map stores: oriented, summed, zeros dropped.
fn stored(records: &[Record], intra: bool, flip: bool) -> BTreeMap<(u32, u32), f32> {
    let mut out = BTreeMap::new();
    for r in records {
        let (mut x, mut y) = if flip { (r.y, r.x) } else { (r.x, r.y) };
        if intra && x > y {
            std::mem::swap(&mut x, &mut y);
        }
        if r.value != 0.0 {
            *out.entry((x, y)).or_insert(0.0) += r.value;
        }
    }
    out.retain(|_, v| *v != 0.0);
    out
}

fn round_trip(seed: u64, version: i32, depth: Option<&str>, with_straw: bool) {
    let mut rng = Rng(seed);
    // Long chromosomes at a fine resolution make blocks wider than shorts.
    let chroms = vec![
        Chrom {
            name: "chrA".into(),
            len: rng.range(40_000_000, 150_000_000),
        },
        Chrom {
            name: "chrB".into(),
            len: rng.range(1_000_000, 90_000_000),
        },
    ];
    let mut header = HicHeader::new("test", chroms, vec![100_000, 2_000]);
    header.version = version;
    if let Some(d) = depth {
        header
            .attributes
            .push((DEPTH_BASE_ATTRIBUTE.into(), d.into()));
    }
    let mut writer = HicWriter::create(Cursor::new(Vec::new()), header.clone(), 1).unwrap();
    let mut written = BTreeMap::new();
    for &res in &header.bp_resolutions {
        for (a, b) in [(0u32, 0u32), (1, 0), (1, 1)] {
            let (na, nb) = (
                header.bins(a, ContactUnit::Bp, res) as u32,
                header.bins(b, ContactUnit::Bp, res) as u32,
            );
            let short = rng.chance(50);
            let records = random_records(&mut rng, na, nb, a == b, short);
            writer
                .write_matrix(a, b, ContactUnit::Bp, res, &records)
                .unwrap();
            written.insert((a.min(b), a.max(b), res), stored(&records, a == b, a > b));
        }
    }
    let bytes = writer.finish().unwrap().into_inner();
    let dir = tempfile::tempdir().unwrap();
    let path = dir.path().join("w.hic");
    std::fs::write(&path, &bytes).unwrap();
    let mut reader = HicReader::open(&bytes[..]).unwrap();

    for _ in 0..30 {
        let &(c1, c2, res) = written
            .keys()
            .nth(rng.range(0, written.len() as u64 - 1) as usize)
            .unwrap();
        let map = &written[&(c1, c2, res)];
        let (n1, n2) = (
            header.bins(c1, ContactUnit::Bp, res) as u32,
            header.bins(c2, ContactUnit::Bp, res) as u32,
        );
        let span = |rng: &mut Rng, n: u32| {
            let start = rng.range(0, n as u64 - 1) as u32;
            (start, (start + rng.range(1, 3_000) as u32).min(n))
        };
        let rows = span(&mut rng, n1);
        let cols = if c1 == c2 && rng.chance(60) {
            rows
        } else {
            span(&mut rng, n2)
        };
        let t = Tile {
            unit: ContactUnit::Bp,
            res,
            r: c1,
            rows,
            c: c2,
            cols,
        };
        let mut want = Cells::new();
        for (&(x, y), &v) in map {
            place(&t, x, y, v.to_bits(), &mut want);
        }
        let got = entries(&mut reader, query(&header, &t));
        if got != want {
            let differences: Vec<_> = got
                .keys()
                .chain(want.keys())
                .filter(|k| got.get(k) != want.get(k))
                .take(10)
                .map(|k| (*k, got.get(k), want.get(k)))
                .collect();
            panic!(
                "v{version} depth {depth:?} {t:?}: got {} cells, want {}; first differences {differences:?}",
                got.len(),
                want.len()
            );
        }
        if with_straw {
            assert_eq!(
                straw_cells(&path, &header, &t),
                want,
                "v{version} {t:?}: straw"
            );
        }
    }
}

proptest! {
    #![proptest_config(ProptestConfig { cases: 8, ..ProptestConfig::default() })]

    #[test]
    fn written_maps_read_back(seed in any::<u64>()) {
        round_trip(seed, 9, None, true);
        round_trip(seed, 8, None, true);
    }

    #[test]
    fn other_v9_depth_layouts_read_back(seed in any::<u64>()) {
        round_trip(seed, 9, Some("3"), false);
        round_trip(seed, 9, Some("-4"), false);
    }
}

#[test]
fn writer_rejects_bad_input() {
    let chroms = vec![Chrom {
        name: "chr1".into(),
        len: 10_000,
    }];
    let header = HicHeader::new("g", chroms, vec![1_000]);
    let create = |h: HicHeader| HicWriter::create(Cursor::new(Vec::new()), h, 6);
    let invalid = |r: omicsweave_core::Result<()>| matches!(r, Err(Error::Invalid { .. }));

    for version in [7, 10] {
        let mut other = header.clone();
        other.version = version;
        assert!(matches!(create(other), Err(Error::Unsupported(_))));
    }
    let mut frag = header.clone();
    frag.frag_resolutions = vec![1];
    assert!(matches!(create(frag), Err(Error::Invalid { .. })));
    let mut long = header.clone();
    long.version = 8;
    long.chroms[0].len = 1 << 31;
    assert!(matches!(create(long), Err(Error::Invalid { .. })));

    let mut w = create(header).unwrap();
    let one = [Record {
        x: 1,
        y: 2,
        value: 1.0,
    }];
    assert!(invalid(w.write_matrix(0, 1, ContactUnit::Bp, 1_000, &one)));
    assert!(invalid(w.write_matrix(0, 0, ContactUnit::Bp, 500, &one)));
    assert!(invalid(w.write_matrix(
        0,
        0,
        ContactUnit::Frag,
        1_000,
        &one
    )));
    let outside = [Record {
        x: 11,
        y: 0,
        value: 1.0,
    }];
    assert!(invalid(w.write_matrix(
        0,
        0,
        ContactUnit::Bp,
        1_000,
        &outside
    )));
    let nan = [Record {
        x: 0,
        y: 0,
        value: f32::NAN,
    }];
    assert!(invalid(w.write_matrix(0, 0, ContactUnit::Bp, 1_000, &nan)));
    w.write_matrix(0, 0, ContactUnit::Bp, 1_000, &one).unwrap();
    assert!(invalid(w.write_matrix(0, 0, ContactUnit::Bp, 1_000, &one)));
    let bytes = w.finish().unwrap().into_inner();
    let reader = HicReader::open(&bytes[..]).unwrap();
    assert_eq!(reader.header().version, 9);
}
