//! The last stage of every conversion into `.hic`: placed pixels, in file
//! order, cut into maps and blocks and written with the file's sidecar.

use std::fs::File;
use std::io::{BufWriter, Read};
use std::path::{Path, PathBuf};

use omicsweave_core::exec::Exec;
use omicsweave_core::{ContactUnit, Result};
use omicsweave_hic::{BlockLayout, HicHeader, HicWriter, MatrixLayouts, Record};

use super::output::Outputs;
use super::sort::{Spill, read_or_end};
use crate::budget::{Budget, Sizing, WRITE_BUF, sizing};
use crate::conversion::Context;
use crate::format::Format;

/// A pixel placed in its map and block, in stored orientation (`c1 <= c2`,
/// and `x <= y` within a chromosome).
#[derive(Debug, Clone, Copy, PartialEq)]
pub(crate) struct Pixel {
    pub c1: u32,
    pub c2: u32,
    pub block: u32,
    pub x: u32,
    pub y: u32,
    pub value: f32,
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
            value: f32::from_bits(u(5)),
        }))
    }
}

/// Blocks in flight between the sink and the file.
pub(crate) fn block_sizing(budget: &Budget) -> Sizing {
    sizing(budget.threads, 1 << 20, budget.memory / 16)
}

/// Writes a `.hic` file (and its `.sidecar`) from pixels given map by map
/// (resolution, chromosome pair) and, within a map, in [`Pixel`] key
/// order. Blocks are encoded and compressed on the pool and written in the
/// order given, so the file depends only on the pixels and their order.
pub(crate) struct HicSink {
    writer: HicWriter<BufWriter<File>>,
    layouts: MatrixLayouts,
    outputs: Outputs,
    dst: PathBuf,
    sidecar: bool,
    map: Option<(u32, u32, u32)>,
    block: Option<u32>,
    cells: Vec<Record>,
    /// Pixels taken so far.
    pub pixels: u64,
}

impl HicSink {
    /// Creates `dst` (published by [`HicSink::finish`]) and writes its header;
    /// maps are cut into blocks as [`BlockLayout::default`] says.
    pub fn create(
        ctx: &Context<'_>,
        dst: &Path,
        header: HicHeader,
        level: u32,
        sidecar: bool,
    ) -> Result<HicSink> {
        let mut outputs = Outputs::new();
        let file = outputs.create(dst, Format::Hic)?;
        let writer = HicWriter::create_exec(
            BufWriter::with_capacity(WRITE_BUF, file),
            header,
            level,
            BlockLayout::default(),
            Exec::Pool {
                pool: &ctx.pool,
                limit: block_sizing(ctx.budget).limit,
            },
        )?;
        Ok(HicSink {
            layouts: writer.layouts(),
            writer,
            outputs,
            dst: dst.to_owned(),
            sidecar,
            map: None,
            block: None,
            cells: Vec::new(),
            pixels: 0,
        })
    }

    /// Where the file's maps put their cells, to place pixels on workers.
    pub fn layouts(&self) -> &MatrixLayouts {
        &self.layouts
    }

    /// Adds a pixel of the map of `res` base pairs.
    pub fn push(&mut self, res: u32, p: Pixel) -> Result<()> {
        if self.map != Some((res, p.c1, p.c2)) {
            self.end_map()?;
            let layout = self.layouts.layout(p.c1, p.c2, ContactUnit::Bp, res)?;
            self.writer.begin_matrix(&layout)?;
            self.map = Some((res, p.c1, p.c2));
        }
        if self.block != Some(p.block) {
            if let Some(n) = self.block {
                self.writer.add_block(n, std::mem::take(&mut self.cells))?;
            }
            self.block = Some(p.block);
        }
        self.cells.push(Record {
            x: p.x,
            y: p.y,
            value: p.value,
        });
        self.pixels += 1;
        Ok(())
    }

    fn end_map(&mut self) -> Result<()> {
        if let Some(n) = self.block.take() {
            self.writer.add_block(n, std::mem::take(&mut self.cells))?;
        }
        if self.map.take().is_some() {
            self.writer.end_matrix()?;
        }
        Ok(())
    }

    /// Writes the indexes and the sidecar, and publishes the files.
    pub fn finish(mut self) -> Result<Vec<(Format, PathBuf)>> {
        self.end_map()?;
        let HicSink {
            writer,
            mut outputs,
            dst,
            sidecar,
            ..
        } = self;
        let (out, built) = writer.finish_with_sidecar()?;
        out.into_inner().map_err(|e| e.into_error())?;
        if sidecar {
            let source = outputs.source(&dst)?;
            let file = outputs.create(
                &omicsweave_core::sidecar::sidecar_path(&dst),
                Format::Sidecar,
            )?;
            built.write_to(&source, file)?;
        }
        outputs.commit()
    }
}
