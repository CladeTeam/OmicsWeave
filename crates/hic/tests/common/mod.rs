//! A minimal `.hic` writer for reader test inputs, independent of the crate's
//! own writer; straw as the reference reader; random fixtures. Fixtures
//! carry normalization and expected vectors, which readers must skip.

#![allow(dead_code)]

use std::collections::BTreeMap;
use std::ffi::CString;
use std::io::Write;
use std::path::Path;

use omicsweave_core::ContactUnit;

/// splitmix64, for deterministic fixtures from a proptest-chosen seed.
pub struct Rng(pub u64);

impl Rng {
    pub fn next(&mut self) -> u64 {
        self.0 = self.0.wrapping_add(0x9E37_79B9_7F4A_7C15);
        let mut z = self.0;
        z = (z ^ (z >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
        z = (z ^ (z >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
        z ^ (z >> 31)
    }

    pub fn range(&mut self, lo: u64, hi: u64) -> u64 {
        lo + self.next() % (hi - lo + 1)
    }

    pub fn chance(&mut self, percent: u64) -> bool {
        self.next() % 100 < percent
    }
}

/// One chromosome pair's map at one resolution. `chr1 <= chr2`; within a
/// chromosome, records keep `x <= y` (the stored half).
#[derive(Debug, Clone)]
pub struct Map {
    pub chr1: u32,
    pub chr2: u32,
    pub unit: ContactUnit,
    pub res: u32,
    pub records: Vec<(u32, u32, f32)>,
}

#[derive(Debug, Clone)]
pub struct Norm {
    pub kind: String,
    pub chrom: u32,
    pub unit: ContactUnit,
    pub res: u32,
    pub values: Vec<f64>,
}

/// An expected vector; `kind` "NONE" for raw counts.
#[derive(Debug, Clone)]
pub struct Expected {
    pub kind: String,
    pub unit: ContactUnit,
    pub res: u32,
    pub values: Vec<f64>,
    pub factors: Vec<(u32, f64)>,
}

#[derive(Debug, Clone)]
pub struct Spec {
    pub version: i32,
    pub chroms: Vec<(String, u64)>,
    pub bp_resolutions: Vec<u32>,
    pub frag_resolutions: Vec<u32>,
    pub sites: Vec<Vec<u32>>,
    pub maps: Vec<Map>,
    pub norms: Vec<Norm>,
    pub expected: Vec<Expected>,
    pub block_bins: u32,
    /// Integral values stored as short; otherwise float.
    pub short_values: bool,
    pub level: u32,
}

pub const NORMS: [&str; 2] = ["KR", "VC"];

impl Spec {
    pub fn resolutions(&self, unit: ContactUnit) -> &[u32] {
        match unit {
            ContactUnit::Bp => &self.bp_resolutions,
            ContactUnit::Frag => &self.frag_resolutions,
        }
    }

    /// A chromosome's span in `unit`.
    pub fn extent(&self, chrom: u32, unit: ContactUnit) -> u64 {
        match unit {
            ContactUnit::Bp => self.chroms[chrom as usize].1,
            ContactUnit::Frag => self.sites[chrom as usize].len() as u64,
        }
    }

    /// Bins of a chromosome (the last one partial).
    pub fn bins(&self, chrom: u32, unit: ContactUnit, res: u32) -> u32 {
        (self.extent(chrom, unit) / res as u64) as u32 + 1
    }

    pub fn units(&self) -> Vec<ContactUnit> {
        if self.frag_resolutions.is_empty() {
            vec![ContactUnit::Bp]
        } else {
            vec![ContactUnit::Bp, ContactUnit::Frag]
        }
    }

    pub fn norm(&self, kind: &str, chrom: u32, unit: ContactUnit, res: u32) -> Option<&[f64]> {
        self.norms
            .iter()
            .find(|n| n.kind == kind && n.chrom == chrom && n.unit == unit && n.res == res)
            .map(|n| n.values.as_slice())
    }

    pub fn expected(&self, kind: &str, unit: ContactUnit, res: u32) -> Option<&Expected> {
        self.expected
            .iter()
            .find(|e| e.kind == kind && e.unit == unit && e.res == res)
    }

    /// A vector value as stored (float from v9).
    pub fn stored(&self, v: f64) -> f64 {
        if self.version > 8 { v as f32 as f64 } else { v }
    }
}

pub fn random_spec(rng: &mut Rng, version: i32) -> Spec {
    let chroms: Vec<(String, u64)> = (0..rng.range(1, 3))
        .map(|i| (format!("chr{}", i + 1), rng.range(20_000, 200_000)))
        .collect();
    let (frag_resolutions, sites) = if rng.chance(50) {
        let mut sites = Vec::new();
        for (_, len) in &chroms {
            let mut s = Vec::new();
            for _ in 0..rng.range(30, 600) {
                s.push(rng.range(1, *len - 1) as u32);
            }
            s.sort_unstable();
            s.dedup();
            sites.push(s);
        }
        (vec![5, 1], sites)
    } else {
        (Vec::new(), Vec::new())
    };
    let mut spec = Spec {
        version,
        chroms,
        bp_resolutions: vec![5_000, 1_000],
        frag_resolutions,
        sites,
        maps: Vec::new(),
        norms: Vec::new(),
        expected: Vec::new(),
        block_bins: rng.range(4, 60) as u32,
        short_values: rng.chance(50),
        level: 0,
    };
    // As the Juicer tools write them: a pair present in the file has a
    // record at every resolution (possibly without blocks).
    let n = spec.chroms.len() as u32;
    let mut pairs = Vec::new();
    for c1 in 0..n {
        for c2 in c1..n {
            if c1 == c2 || rng.chance(70) {
                pairs.push((c1, c2));
            }
        }
    }
    for unit in spec.units() {
        for res in spec.resolutions(unit).to_vec() {
            for &(c1, c2) in &pairs {
                let (n1, n2) = (spec.bins(c1, unit, res) - 1, spec.bins(c2, unit, res) - 1);
                let mut cells = BTreeMap::new();
                for _ in 0..rng.range(0, 3_000) {
                    let (mut x, mut y) = (
                        rng.range(0, n1 as u64) as u32,
                        rng.range(0, n2 as u64) as u32,
                    );
                    if c1 == c2 {
                        // Contacts concentrate near the diagonal.
                        let reach = 1 << rng.range(0, 7);
                        y = (x + rng.range(0, reach) as u32).min(n2);
                        if y < x {
                            std::mem::swap(&mut x, &mut y);
                        }
                    }
                    let v = if spec.short_values {
                        rng.range(1, 500) as f32
                    } else {
                        rng.range(1, 100_000) as f32 / 37.0
                    };
                    cells.insert((x, y), v);
                }
                spec.maps.push(Map {
                    chr1: c1,
                    chr2: c2,
                    unit,
                    res,
                    records: cells.into_iter().map(|((x, y), v)| (x, y, v)).collect(),
                });
            }
            for kind in NORMS {
                for chrom in 0..n {
                    let values = (0..spec.bins(chrom, unit, res))
                        .map(|_| match rng.range(0, 30) {
                            0 => f64::NAN,
                            1 => 0.0,
                            _ => rng.range(1, 4_000) as f64 / 1000.0,
                        })
                        .collect();
                    spec.norms.push(Norm {
                        kind: kind.into(),
                        chrom,
                        unit,
                        res,
                        values,
                    });
                }
            }
            let longest = (0..n).map(|c| spec.bins(c, unit, res)).max().unwrap();
            for kind in ["NONE", NORMS[0], NORMS[1]] {
                // Shorter vectors than a chromosome too: distances past the
                // end read the last value.
                let len = rng.range(1, longest as u64 + 3);
                let values = (0..len)
                    .map(|_| rng.range(1, 90_000) as f64 / 7.0)
                    .collect();
                let mut factors = Vec::new();
                for c in 0..n {
                    if rng.chance(70) {
                        factors.push((c, rng.range(1, 3_000) as f64 / 1000.0));
                    }
                }
                spec.expected.push(Expected {
                    kind: kind.into(),
                    unit,
                    res,
                    values,
                    factors,
                });
            }
        }
    }
    spec
}

struct Out(Vec<u8>);

impl Out {
    fn i16(&mut self, v: i16) {
        self.0.extend_from_slice(&v.to_le_bytes());
    }
    fn i32(&mut self, v: i32) {
        self.0.extend_from_slice(&v.to_le_bytes());
    }
    fn i64(&mut self, v: i64) {
        self.0.extend_from_slice(&v.to_le_bytes());
    }
    fn f32(&mut self, v: f32) {
        self.0.extend_from_slice(&v.to_le_bytes());
    }
    fn f64(&mut self, v: f64) {
        self.0.extend_from_slice(&v.to_le_bytes());
    }
    fn str(&mut self, s: &str) {
        self.0.extend_from_slice(s.as_bytes());
        self.0.push(0);
    }
    /// `int` before v9, `long` from v9.
    fn count(&mut self, version: i32, v: i64) {
        if version > 8 {
            self.i64(v)
        } else {
            self.i32(v as i32)
        }
    }
    /// float from v9, double before.
    fn value(&mut self, version: i32, v: f64) {
        if version > 8 {
            self.f32(v as f32)
        } else {
            self.f64(v)
        }
    }
    fn pos(&self) -> i64 {
        self.0.len() as i64
    }
}

/// The block a record belongs to, as the Juicer tools number them.
fn block_number(spec: &Spec, map: &Map, columns: i64, x: u32, y: u32) -> i64 {
    let bb = spec.block_bins as i64;
    let (x, y) = (x as i64, y as i64);
    if spec.version > 8 && map.chr1 == map.chr2 {
        let pad = (x + y) / 2 / bb;
        let depth =
            (1.0 + (x - y).abs() as f64 / std::f64::consts::SQRT_2 / bb as f64).log2() as i64;
        depth * columns + pad
    } else {
        (y / bb) * columns + x / bb
    }
}

/// Encodes one block's records, choosing a representation and field widths
/// at random within what the version allows.
fn encode_block(spec: &Spec, records: &[(u32, u32, f32)], rng: &mut Rng) -> Vec<u8> {
    let x0 = records.iter().map(|r| r.0).min().unwrap();
    let y0 = records.iter().map(|r| r.1).min().unwrap();
    let w = records.iter().map(|r| r.0 - x0).max().unwrap() + 1;
    let h = records.iter().map(|r| r.1 - y0).max().unwrap() + 1;
    let mut b = Out(Vec::new());
    b.i32(records.len() as i32);
    b.i32(x0 as i32);
    b.i32(y0 as i32);
    b.0.push(!spec.short_values as u8);
    let wide = spec.version > 8 && rng.chance(50);
    let (wide_x, wide_y) = (wide && rng.chance(50), wide && rng.chance(50));
    if spec.version > 8 {
        b.0.push(wide_x as u8);
        b.0.push(wide_y as u8);
    }
    let value = |b: &mut Out, v: f32| {
        if spec.short_values {
            b.i16(v as i16)
        } else {
            b.f32(v)
        }
    };
    if (w as u64) * (h as u64) <= 2_000 && w < 32_000 && rng.chance(40) {
        b.0.push(2);
        b.i32((w * h) as i32);
        b.i16(w as i16);
        let mut grid = vec![None; (w * h) as usize];
        for &(x, y, v) in records {
            grid[((y - y0) * w + (x - x0)) as usize] = Some(v);
        }
        for cell in grid {
            match cell {
                Some(v) => value(&mut b, v),
                None if spec.short_values => b.i16(i16::MIN),
                None => b.f32(f32::NAN),
            }
        }
    } else {
        b.0.push(1);
        let mut rows: BTreeMap<u32, Vec<(u32, f32)>> = BTreeMap::new();
        for &(x, y, v) in records {
            rows.entry(y - y0).or_default().push((x - x0, v));
        }
        let int = |b: &mut Out, wide: bool, v: u32| {
            if wide {
                b.i32(v as i32)
            } else {
                b.i16(v as i16)
            }
        };
        int(&mut b, wide_y, rows.len() as u32);
        for (row, cols) in rows {
            int(&mut b, wide_y, row);
            int(&mut b, wide_x, cols.len() as u32);
            for (col, v) in cols {
                int(&mut b, wide_x, col);
                value(&mut b, v);
            }
        }
    }
    let mut z = flate2::write::ZlibEncoder::new(Vec::new(), flate2::Compression::new(spec.level));
    z.write_all(&b.0).unwrap();
    z.finish().unwrap()
}

/// A map whose blocks are written: (map index, block columns, block index
/// entries as (number, position, size)).
type WrittenMap = (usize, i64, Vec<(i64, i64, i32)>);

fn expected_list(o: &mut Out, spec: &Spec, normalized: bool) {
    let list: Vec<&Expected> = spec
        .expected
        .iter()
        .filter(|e| (e.kind != "NONE") == normalized)
        .collect();
    let v = spec.version;
    o.i32(list.len() as i32);
    for e in list {
        if normalized {
            o.str(&e.kind);
        }
        o.str(e.unit.as_str());
        o.i32(e.res as i32);
        o.count(v, e.values.len() as i64);
        for &x in &e.values {
            o.value(v, x);
        }
        o.i32(e.factors.len() as i32);
        for &(c, f) in &e.factors {
            o.i32(c as i32);
            o.value(v, f);
        }
    }
}

pub fn write(spec: &Spec, rng: &mut Rng) -> Vec<u8> {
    let v = spec.version;
    let mut o = Out(Vec::new());
    o.0.extend_from_slice(b"HIC\0");
    o.i32(v);
    let footer_at = o.0.len();
    o.i64(0);
    o.str("hg38");
    let nvi_at = o.0.len();
    if v > 8 {
        o.i64(0);
        o.i64(0);
    }
    o.i32(1);
    o.str("software");
    o.str("omicsweave test writer");
    o.i32(spec.chroms.len() as i32);
    for (name, len) in &spec.chroms {
        o.str(name);
        o.count(v, *len as i64);
    }
    for list in [&spec.bp_resolutions, &spec.frag_resolutions] {
        o.i32(list.len() as i32);
        for &r in list.iter() {
            o.i32(r as i32);
        }
    }
    if !spec.frag_resolutions.is_empty() {
        for sites in &spec.sites {
            o.i32(sites.len() as i32);
            for &s in sites {
                o.i32(s as i32);
            }
        }
    }

    // Blocks, then each pair's matrix record.
    let mut pairs: BTreeMap<(u32, u32), Vec<WrittenMap>> = BTreeMap::new();
    for (mi, map) in spec.maps.iter().enumerate() {
        let n = spec
            .bins(map.chr1, map.unit, map.res)
            .max(spec.bins(map.chr2, map.unit, map.res)) as i64;
        let columns = n / spec.block_bins as i64 + 2;
        let mut blocks: BTreeMap<i64, Vec<(u32, u32, f32)>> = BTreeMap::new();
        for &r in &map.records {
            blocks
                .entry(block_number(spec, map, columns, r.0, r.1))
                .or_default()
                .push(r);
        }
        let mut index = Vec::new();
        for (number, records) in blocks {
            let bytes = encode_block(spec, &records, rng);
            index.push((number, o.pos(), bytes.len() as i32));
            o.0.extend_from_slice(&bytes);
        }
        pairs
            .entry((map.chr1, map.chr2))
            .or_default()
            .push((mi, columns, index));
    }
    let mut matrices = Vec::new();
    for ((c1, c2), resolutions) in &pairs {
        let start = o.pos();
        o.i32(*c1 as i32);
        o.i32(*c2 as i32);
        o.i32(resolutions.len() as i32);
        for (mi, columns, index) in resolutions {
            let map = &spec.maps[*mi];
            o.str(map.unit.as_str());
            let list = spec.resolutions(map.unit);
            o.i32(list.iter().position(|&r| r == map.res).unwrap() as i32);
            o.f32(map.records.iter().map(|r| r.2).sum());
            o.f32(0.0);
            o.f32(0.0);
            o.f32(0.0);
            o.i32(map.res as i32);
            o.i32(spec.block_bins as i32);
            o.i32(*columns as i32);
            o.i32(index.len() as i32);
            for &(number, position, size) in index {
                o.i32(number as i32);
                o.i64(position);
                o.i32(size);
            }
        }
        matrices.push((*c1, *c2, start, (o.pos() - start) as i32));
    }
    let mut norms = Vec::new();
    for n in &spec.norms {
        let start = o.pos();
        o.count(v, n.values.len() as i64);
        for &x in &n.values {
            o.value(v, x);
        }
        norms.push((n, start, o.pos() - start));
    }

    // Footer: master index and raw expected vectors ("nBytesV5"), then the
    // normalized expected vectors and the norm vector index.
    let footer = o.pos();
    let mut f = Out(Vec::new());
    f.i32(matrices.len() as i32);
    for (c1, c2, position, size) in &matrices {
        f.str(&format!("{c1}_{c2}"));
        f.i64(*position);
        f.i32(*size);
    }
    expected_list(&mut f, spec, false);
    o.count(v, f.0.len() as i64);
    o.0.extend_from_slice(&f.0);
    expected_list(&mut o, spec, true);
    let nvi = o.pos();
    o.i32(norms.len() as i32);
    for (n, position, nbytes) in &norms {
        o.str(&n.kind);
        o.i32(n.chrom as i32);
        o.str(n.unit.as_str());
        o.i32(n.res as i32);
        o.i64(*position);
        o.count(v, *nbytes);
    }
    let nvi_len = o.pos() - nvi;
    o.0[footer_at..footer_at + 8].copy_from_slice(&footer.to_le_bytes());
    if v > 8 {
        o.0[nvi_at..nvi_at + 8].copy_from_slice(&nvi.to_le_bytes());
        o.0[nvi_at + 8..nvi_at + 16].copy_from_slice(&nvi_len.to_le_bytes());
    }
    o.0
}

/// straw's raw observed records (`observed`, `NONE`) for
/// `chrom1:start1:end1` x `chrom2:start2:end2` (inclusive bin-start bounds
/// in `unit`): (x, y, count), x on the lower-indexed chromosome.
pub fn straw(
    path: &Path,
    loc1: &str,
    loc2: &str,
    unit: ContactUnit,
    res: u32,
) -> Vec<(u32, u32, f32)> {
    let c = |s: &str| CString::new(s).unwrap();
    let (t, n, f, l1, l2, u) = (
        c("observed"),
        c("NONE"),
        c(path.to_str().unwrap()),
        c(loc1),
        c(loc2),
        c(unit.as_str()),
    );
    // SAFETY: valid C strings; the result is released after copying.
    unsafe {
        let r = straw_sys::straw_query(
            t.as_ptr(),
            n.as_ptr(),
            f.as_ptr(),
            l1.as_ptr(),
            l2.as_ptr(),
            u.as_ptr(),
            res as i32,
        );
        assert!(!r.is_null(), "straw failed");
        let n = (*r).n as usize;
        let out = (0..n)
            .map(|i| {
                (
                    *(*r).bin_x.add(i) as u32,
                    *(*r).bin_y.add(i) as u32,
                    *(*r).counts.add(i),
                )
            })
            .collect();
        straw_sys::straw_free(r);
        out
    }
}
