//! Contacts → Juicer `.hic` (v9, raw counts, with an optional sidecar)
//! through the [`omicsweave_hic::HicWriter`]: every contact binned at each
//! resolution as the Juicer tools' `pre` bins a pairs file (the two ends
//! by reference index, positions clamped to the reference length, `bin =
//! position / resolution`, within a reference the upper triangle), cells
//! sorted out of core per resolution and folded, maps written coarsest
//! first, each in block order.

use std::fs::File;
use std::io::{BufWriter, Read};
use std::path::Path;

use omicsweave_convert::stage::{Sorter, Spill, read_or_end};
use omicsweave_core::exec::Exec;
use omicsweave_core::{ContactUnit, Error, Result};
use omicsweave_hic::{
    BlockLayout, Chrom, HicHeader, HicWriter, MatrixLayout, MatrixLayouts, Record,
};

use crate::accum::Delivery;

/// zlib level of blocks, as the Juicer tools default to.
const LEVEL: u32 = 6;

/// A cell of one map, placed in its block (stored orientation).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
struct Pixel {
    c1: u32,
    c2: u32,
    block: u32,
    x: u32,
    y: u32,
    count: u32,
}

impl Spill for Pixel {
    /// File order within one resolution.
    type Key = (u32, u32, u32, u32, u32);
    fn key(&self) -> Self::Key {
        (self.c1, self.c2, self.block, self.y, self.x)
    }
    fn mem_size(&self) -> usize {
        0
    }
    fn encode(&self, out: &mut Vec<u8>) {
        for v in [self.c1, self.c2, self.block, self.x, self.y, self.count] {
            out.extend_from_slice(&v.to_le_bytes());
        }
    }
    fn decode(input: &mut impl Read) -> Result<Option<Pixel>> {
        let mut b = [0u8; 24];
        if !read_or_end(input, &mut b)? {
            return Ok(None);
        }
        let u = |i: usize| u32::from_le_bytes(b[4 * i..4 * i + 4].try_into().unwrap());
        Ok(Some(Pixel {
            c1: u(0),
            c2: u(1),
            block: u(2),
            x: u(3),
            y: u(4),
            count: u(5),
        }))
    }
    const ABSORBS: bool = true;
    fn absorb(&mut self, next: Pixel) {
        self.count = self.count.saturating_add(next.count);
    }
}

/// Cells gathered per resolution before they enter its sorter.
const STAGE: usize = 1 << 16;

fn fold(v: &mut Vec<Pixel>) {
    v.sort_unstable_by_key(Pixel::key);
    let mut w = 0;
    for i in 0..v.len() {
        if w > 0 && v[w - 1].key() == v[i].key() {
            v[w - 1].count = v[w - 1].count.saturating_add(v[i].count);
        } else {
            v[w] = v[i];
            w += 1;
        }
    }
    v.truncate(w);
}

/// Contacts binned at every resolution of a `.hic` header.
pub struct HicPixels<'p> {
    header: HicHeader,
    layouts: MatrixLayouts,
    lens: Vec<u32>,
    /// Per resolution (header order): its sorter, staged cells and the
    /// layout of the last map placed in.
    maps: Vec<(Sorter<'p, Pixel>, Vec<Pixel>, Option<MatrixLayout>)>,
}

impl<'p> HicPixels<'p> {
    /// A file over `chroms` at `resolutions` (written coarsest first), its
    /// sorters sharing the delivery's buffer.
    pub fn new(
        chroms: &[(String, u32)],
        genome: &str,
        resolutions: &[u32],
        d: &Delivery<'p>,
    ) -> Result<HicPixels<'p>> {
        let mut res = resolutions.to_vec();
        res.sort_unstable_by(|a, b| b.cmp(a));
        res.dedup();
        if res.contains(&0) {
            return Err(Error::invalid("hic", "a resolution must be at least 1 bp"));
        }
        let header = HicHeader::new(
            genome,
            chroms
                .iter()
                .map(|(n, l)| Chrom {
                    name: n.clone(),
                    len: *l as u64,
                })
                .collect(),
            res.clone(),
        );
        let layouts = MatrixLayouts::new(&header, BlockLayout::default())?;
        let share = d.buffer / res.len().max(1) as u64;
        let maps = res
            .iter()
            .map(|_| (Sorter::new(d.pool, share, d.tmp), Vec::new(), None))
            .collect();
        Ok(HicPixels {
            header,
            layouts,
            lens: chroms.iter().map(|c| c.1).collect(),
            maps,
        })
    }

    /// Adds a contact between position `p1` of `c1` and `p2` of `c2`
    /// (`c1 <= c2`).
    pub fn push(&mut self, c1: u32, p1: u32, c2: u32, p2: u32) -> Result<()> {
        let p1 = p1.min(self.lens[c1 as usize]);
        let p2 = p2.min(self.lens[c2 as usize]);
        for (i, &res) in self.header.bp_resolutions.iter().enumerate() {
            let (sorter, staged, last) = &mut self.maps[i];
            let layout = match last {
                Some(l) if l.chroms() == (c1, c2) => *l,
                _ => {
                    let l = self.layouts.layout(c1, c2, ContactUnit::Bp, res)?;
                    *last = Some(l);
                    l
                }
            };
            let (block, x, y) = layout.place(p1 / res, p2 / res);
            staged.push(Pixel {
                c1,
                c2,
                block,
                x,
                y,
                count: 1,
            });
            if staged.len() >= STAGE {
                fold(staged);
                sorter.extend(std::mem::take(staged))?;
            }
        }
        Ok(())
    }

    /// Writes the file at `path` (and its sidecar); returns the number of
    /// cells written over all resolutions.
    pub fn finish(self, path: &Path, sidecar: bool, d: &mut Delivery<'_>) -> Result<u64> {
        let HicPixels {
            header,
            layouts,
            maps,
            ..
        } = self;
        let resolutions = header.bp_resolutions.clone();
        let file = d.outputs.create(path, "hic")?;
        let mut w: HicWriter<BufWriter<File>> = HicWriter::create_exec(
            BufWriter::with_capacity(1 << 20, file),
            header,
            LEVEL,
            BlockLayout::default(),
            Exec::Pool {
                pool: d.pool,
                limit: 4 * d.pool.threads(),
            },
        )?;
        let mut cells = 0u64;
        for ((mut sorter, mut staged, _), res) in maps.into_iter().zip(resolutions) {
            fold(&mut staged);
            sorter.extend(staged)?;
            let mut merged = sorter.finish()?;
            let mut map: Option<(u32, u32)> = None;
            let mut block: Option<u32> = None;
            let mut records: Vec<Record> = Vec::new();
            let mut cur: Option<Pixel> = None;
            let flush = |p: Pixel,
                         w: &mut HicWriter<BufWriter<File>>,
                         map: &mut Option<(u32, u32)>,
                         block: &mut Option<u32>,
                         records: &mut Vec<Record>|
             -> Result<()> {
                if *map != Some((p.c1, p.c2)) {
                    if let Some(n) = block.take() {
                        w.add_block(n, std::mem::take(records))?;
                    }
                    if map.is_some() {
                        w.end_matrix()?;
                    }
                    w.begin_matrix(&layouts.layout(p.c1, p.c2, ContactUnit::Bp, res)?)?;
                    *map = Some((p.c1, p.c2));
                }
                if *block != Some(p.block) {
                    if let Some(n) = block.take() {
                        w.add_block(n, std::mem::take(records))?;
                    }
                    *block = Some(p.block);
                }
                records.push(Record {
                    x: p.x,
                    y: p.y,
                    value: p.count as f32,
                });
                Ok(())
            };
            while let Some(p) = merged.next_record()? {
                match &mut cur {
                    Some(c) if c.key() == p.key() => c.count = c.count.saturating_add(p.count),
                    _ => {
                        if let Some(c) = cur.replace(p) {
                            cells += 1;
                            flush(c, &mut w, &mut map, &mut block, &mut records)?;
                        }
                    }
                }
            }
            if let Some(c) = cur.take() {
                cells += 1;
                flush(c, &mut w, &mut map, &mut block, &mut records)?;
            }
            if let Some(n) = block.take() {
                w.add_block(n, std::mem::take(&mut records))?;
            }
            if map.is_some() {
                w.end_matrix()?;
            }
        }
        let (out, car) = w.finish_with_sidecar()?;
        out.into_inner().map_err(|e| e.into_error())?;
        if sidecar {
            let source = d.outputs.source(path)?;
            let file = d
                .outputs
                .create(&omicsweave_core::sidecar::sidecar_path(path), "sidecar")?;
            car.write_to(&source, file)?;
        }
        Ok(cells)
    }
}
