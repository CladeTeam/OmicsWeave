//! Independent .hic byte writer for reader fixtures.
use super::fixtures::{Expected, Map, Rng, Spec};
use std::collections::BTreeMap;
use std::io::Write;

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
