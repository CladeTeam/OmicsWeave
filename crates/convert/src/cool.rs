//! cool / mcool → .hic (+ sidecar).

use std::collections::HashMap;
use std::io::{BufWriter, Read};
use std::path::PathBuf;
use std::sync::Arc;

use omicsweave_cool::{Cool, PixelBatch, PixelScratch, StorageMode, is_mcool, mcool_resolutions};
use omicsweave_core::exec::{Exec, Ordered};
use omicsweave_core::{ContactUnit, Error, Result};
use omicsweave_hic::{
    BlockLayout, Chrom, HicHeader, HicWriter, MatrixLayout, MatrixLayouts, Record,
};

use crate::budget::{Budget, MemoryPlan, WRITE_BUF, sizing};
use crate::conversion::{Context, Report};
use crate::format::Format;
use crate::pipeline::sort::{Sorter, Spill, read_or_end};

/// cooler `.cool` / `.mcool` to Juicer `.hic` (v9) with the raw counts, one
/// `.hic` resolution per cool resolution.
///
/// Pixels are read in parallel, placed into the `.hic` block that holds
/// them, sorted by (chromosome pair, block) in memory or out of core, and
/// written block by block, blocks compressed in parallel. Balancing weights
/// are not read (normalization is analysis, not format). Only the upper
/// triangle of a `square` cooler is taken: `.hic` stores each cell once.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct CoolToHic {
    pub src: PathBuf,
    pub dst: PathBuf,
    /// Resolutions of an `.mcool` to convert; `None` takes all.
    pub resolutions: Option<Vec<u32>>,
    /// zlib level of blocks.
    pub level: u32,
    /// Also write `<dst>.sidecar`.
    pub sidecar: bool,
}

/// Bytes per pixel in flight: the decoded columns and the placed record.
const PIXEL_BYTES: u64 = 24 + 24;

impl CoolToHic {
    pub fn new(src: impl Into<PathBuf>, dst: impl Into<PathBuf>) -> CoolToHic {
        CoolToHic {
            src: src.into(),
            dst: dst.into(),
            resolutions: None,
            level: 6,
            sidecar: true,
        }
    }

    /// The cooler groups to convert and their bin sizes, finest first.
    fn groups(&self) -> Result<Vec<(String, u32)>> {
        if !is_mcool(&self.src)? {
            if self.resolutions.is_some() {
                return Err(Error::invalid(
                    "cool",
                    "resolutions apply to .mcool files only",
                ));
            }
            let cool = Cool::open(&self.src, "/")?;
            return Ok(vec![("/".into(), cool.bin_size())]);
        }
        let all = mcool_resolutions(&self.src)?;
        let chosen = match &self.resolutions {
            None => all,
            Some(want) => {
                for r in want {
                    if !all.contains(r) {
                        return Err(Error::invalid(
                            "cool",
                            format!("resolution {r} is not in the file (it has {all:?})"),
                        ));
                    }
                }
                let mut w = want.clone();
                w.sort_unstable();
                w.dedup();
                w
            }
        };
        if chosen.is_empty() {
            return Err(Error::invalid("cool", "no resolutions selected"));
        }
        Ok(chosen
            .into_iter()
            .map(|r| (format!("/resolutions/{r}"), r))
            .collect())
    }

    fn rows(budget: &Budget) -> u64 {
        (budget.memory
            / (budget.threads.max(1) as u64)
                .saturating_mul(16)
                .saturating_mul(PIXEL_BYTES))
        .clamp(1 << 12, 1 << 20)
    }

    fn span_sizing(budget: &Budget) -> crate::budget::Sizing {
        sizing(
            budget.threads,
            Self::rows(budget) * PIXEL_BYTES,
            budget.memory / 4,
        )
    }

    fn block_sizing(budget: &Budget) -> crate::budget::Sizing {
        sizing(budget.threads, 1 << 20, budget.memory / 16)
    }

    pub(crate) fn plan(&self, budget: &Budget) -> Result<MemoryPlan> {
        let mut fixed = 0u64;
        for (group, _) in self.groups()? {
            let cool = Cool::open(&self.src, &group)?;
            // The bins table is read and checked when a group opens.
            fixed = fixed.max(cool.n_bins().saturating_mul(24));
        }
        let pipeline = Self::span_sizing(budget).per_queue
            + Self::block_sizing(budget).per_queue
            + WRITE_BUF as u64;
        Ok(MemoryPlan {
            fixed,
            pipeline,
            buffer: budget.memory.saturating_sub(fixed.saturating_add(pipeline)),
        })
    }

    pub(crate) fn run(&self, ctx: &Context<'_>) -> Result<Report> {
        let groups = self.groups()?;
        let first = Cool::open(&self.src, &groups[0].0)?;
        let chroms: Vec<(String, u32)> = first
            .chroms()
            .iter()
            .map(|c| (c.name.clone(), c.length))
            .collect();
        let genome = first.assembly().unwrap_or("unknown").to_owned();
        drop(first);
        // Coarsest first, as the Juicer tools list resolutions.
        let mut resolutions: Vec<u32> = groups.iter().map(|g| g.1).collect();
        resolutions.reverse();
        let header = HicHeader::new(
            &genome,
            chroms
                .iter()
                .map(|(n, l)| Chrom {
                    name: n.clone(),
                    len: *l as u64,
                })
                .collect(),
            resolutions,
        );
        let mut outputs = crate::pipeline::output::Outputs::new();
        let file = outputs.create(&self.dst, Format::Hic)?;
        let b = Self::block_sizing(ctx.budget);
        let mut writer = HicWriter::create_exec(
            BufWriter::with_capacity(WRITE_BUF, file),
            header,
            self.level,
            BlockLayout::default(),
            Exec::Pool {
                pool: &ctx.pool,
                limit: b.limit,
            },
        )?;
        let layouts = Arc::new(writer.layouts());
        let tmp = ctx.tmp_parent(&self.dst);
        let (mut records, mut spilled) = (0u64, 0u64);
        let total: u64 = groups.len() as u64;
        for (done, (group, res)) in groups.iter().enumerate() {
            let cool = Arc::new(Cool::open(&self.src, group)?);
            let same = cool.chroms().len() == chroms.len()
                && cool
                    .chroms()
                    .iter()
                    .zip(&chroms)
                    .all(|(c, (n, l))| &c.name == n && c.length == *l);
            if !same {
                return Err(Error::invalid(
                    "cool",
                    format!("{group}: chromosomes differ from the other resolutions"),
                ));
            }
            let mut sorter = Sorter::new(&ctx.pool, ctx.memory.buffer, &tmp);
            let s = Self::span_sizing(ctx.budget);
            let mut exec = ctx.pool.ordered::<Vec<Placed>>(s.limit);
            let spans = cool.pixel_spans(Self::rows(ctx.budget));
            for span in spans {
                while exec.full() {
                    sorter.extend(exec.next(true)?.expect("pending"))?;
                }
                let (cool, layouts, res) = (Arc::clone(&cool), Arc::clone(&layouts), *res);
                exec.submit(Box::new(move || place(&cool, &layouts, res, span)))?;
                while let Some(ps) = exec.next(false)? {
                    sorter.extend(ps)?;
                }
                ctx.progress("read", done as u64, total);
            }
            while let Some(ps) = exec.next(true)? {
                sorter.extend(ps)?;
            }
            spilled += sorter.spilled();
            ctx.progress("sort", done as u64, total);
            let mut merged = sorter.finish()?;
            ctx.progress("write", done as u64, total);
            let mut current: Option<(u32, u32)> = None;
            let mut block: Option<u32> = None;
            let mut cells: Vec<Record> = Vec::new();
            while let Some(p) = merged.next_record()? {
                records += 1;
                if current != Some((p.c1, p.c2)) {
                    if let Some(n) = block.take() {
                        writer.add_block(n, std::mem::take(&mut cells))?;
                    }
                    if current.is_some() {
                        writer.end_matrix()?;
                    }
                    let layout: MatrixLayout = layouts.layout(p.c1, p.c2, ContactUnit::Bp, *res)?;
                    writer.begin_matrix(&layout)?;
                    current = Some((p.c1, p.c2));
                }
                if block != Some(p.block) {
                    if let Some(n) = block {
                        writer.add_block(n, std::mem::take(&mut cells))?;
                    }
                    block = Some(p.block);
                }
                cells.push(Record {
                    x: p.x,
                    y: p.y,
                    value: p.value,
                });
            }
            if let Some(n) = block {
                writer.add_block(n, cells)?;
            }
            if current.is_some() {
                writer.end_matrix()?;
            }
        }
        let (out, sidecar) = writer.finish_with_sidecar()?;
        out.into_inner().map_err(|e| e.into_error())?;
        if self.sidecar {
            let source = outputs.source(&self.dst)?;
            let file = outputs.create(
                &omicsweave_core::sidecar::sidecar_path(&self.dst),
                Format::Sidecar,
            )?;
            sidecar.write_to(&source, file)?;
        }
        let written = outputs.commit()?;
        Ok(Report {
            files: written,
            records,
            bytes_read: std::fs::metadata(&self.src).map_or(0, |m| m.len()),
            spilled,
            elapsed: Default::default(),
        })
    }
}

/// A pixel placed in its `.hic` map and block (stored orientation).
struct Placed {
    c1: u32,
    c2: u32,
    block: u32,
    x: u32,
    y: u32,
    value: f32,
}

impl Spill for Placed {
    type Key = (u32, u32, u32, u32, u32);
    fn key(&self) -> Self::Key {
        (self.c1, self.c2, self.block, self.y, self.x)
    }
    fn mem_size(&self) -> usize {
        0
    }
    fn encode(&self, out: &mut Vec<u8>) {
        for v in [
            self.c1,
            self.c2,
            self.block,
            self.x,
            self.y,
            self.value.to_bits(),
        ] {
            out.extend_from_slice(&v.to_le_bytes());
        }
    }
    fn decode(input: &mut impl Read) -> Result<Option<Placed>> {
        let mut b = [0u8; 24];
        if !read_or_end(input, &mut b)? {
            return Ok(None);
        }
        let u = |i: usize| u32::from_le_bytes(b[4 * i..4 * i + 4].try_into().unwrap());
        Ok(Some(Placed {
            c1: u(0),
            c2: u(1),
            block: u(2),
            x: u(3),
            y: u(4),
            value: f32::from_bits(u(5)),
        }))
    }
}

/// Reads one span of pixels and places each in its map and block.
fn place(
    cool: &Cool,
    layouts: &MatrixLayouts,
    res: u32,
    rows: std::ops::Range<u64>,
) -> Result<Vec<Placed>> {
    let mut batch = PixelBatch::new();
    cool.read_pixels(rows, &mut PixelScratch::new(), &mut batch)?;
    let square = cool.storage_mode() == StorageMode::Square;
    let mut cache: HashMap<(u32, u32), MatrixLayout> = HashMap::new();
    let mut out = Vec::with_capacity(batch.len());
    for i in 0..batch.len() {
        let (b1, b2, count) = (batch.bin1[i], batch.bin2[i], batch.count[i]);
        if square && b1 > b2 {
            continue;
        }
        if count == 0.0 {
            continue;
        }
        if !count.is_finite() {
            return Err(Error::invalid(
                "cool",
                format!("pixel ({b1}, {b2}) has a non-finite count"),
            ));
        }
        let ((c1, x), (c2, y)) = (cool.locate(b1), cool.locate(b2));
        let layout = match cache.get(&(c1, c2)) {
            Some(l) => *l,
            None => {
                let l = layouts.layout(c1, c2, ContactUnit::Bp, res)?;
                cache.insert((c1, c2), l);
                l
            }
        };
        let (block, x, y) = layout.place(x, y);
        let (c1, c2) = layout.chroms();
        out.push(Placed {
            c1,
            c2,
            block,
            x,
            y,
            value: count as f32,
        });
    }
    Ok(out)
}
