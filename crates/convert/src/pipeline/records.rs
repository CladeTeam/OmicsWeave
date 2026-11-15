//! The pipeline shared by text-to-BGZF record conversions (VCF → BCF,
//! SAM → BAM): line batches encoded in parallel, records framed into BGZF
//! blocks in order exactly as htslib frames them, blocks compressed in
//! parallel, and the binning index built from the records' positions as
//! they are written.

use std::io::{BufRead, BufWriter};
use std::path::{Path, PathBuf};
use std::sync::Arc;
use std::sync::atomic::{AtomicU64, Ordering};

use omicsweave_bgzf::{BinningIndex, IndexFormat, Writer};
use omicsweave_core::Result;
use omicsweave_core::exec::Ordered;

use super::input::{LineBatch, LineBatches, READ_CHUNK, lines};
use super::output::Outputs;
use crate::budget::{Budget, MemoryPlan, Sizing, WRITE_BUF, sizing, text_batch};
use crate::conversion::{Context, Report};
use crate::format::Format;

/// Where a record goes in the index (`hts_idx_push`'s arguments).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) struct Entry {
    pub tid: i32,
    pub beg: i64,
    pub end: i64,
    pub mapped: bool,
}

/// Encodes one line of text into one binary record.
pub(crate) trait LineEncoder: Send + Sync + 'static {
    type Scratch: Send;
    fn scratch(&self) -> Self::Scratch;
    /// Appends the record of `line` (no `\n`) to `out`; `None` skips the
    /// line (blank).
    fn encode(
        &self,
        scratch: &mut Self::Scratch,
        line: &[u8],
        line_no: u64,
        out: &mut Vec<u8>,
    ) -> Result<Option<Entry>>;
}

/// How records are framed.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) enum Framing {
    /// `bgzf_write` only (BCF): records may straddle blocks.
    Stream,
    /// `bgzf_flush_try` before each record (BAM).
    Record,
}

struct Encoded {
    bytes: Vec<u8>,
    /// Length and index entry of each record.
    records: Vec<(u32, Entry)>,
}

fn encode_batch<E: LineEncoder>(enc: &E, batch: &LineBatch) -> Result<Encoded> {
    let mut scratch = enc.scratch();
    let mut out = Encoded {
        bytes: Vec::with_capacity(batch.text.len()),
        records: Vec::new(),
    };
    for (no, line) in lines(batch) {
        let at = out.bytes.len();
        if let Some(entry) = enc.encode(&mut scratch, line, no, &mut out.bytes)? {
            out.records.push(((out.bytes.len() - at) as u32, entry));
        }
    }
    Ok(out)
}

pub(crate) fn batch_sizing(budget: &Budget) -> Sizing {
    // Text and its encoding.
    sizing(
        budget.threads,
        2 * text_batch(budget) as u64,
        budget.memory / 4,
    )
}

pub(crate) fn block_sizing(budget: &Budget) -> Sizing {
    sizing(
        budget.threads,
        3 * omicsweave_bgzf::MAX_BLOCK_LEN as u64,
        budget.memory / 8,
    )
}

pub(crate) fn plan(budget: &Budget, fixed: u64) -> MemoryPlan {
    MemoryPlan {
        fixed,
        pipeline: batch_sizing(budget).per_queue
            + 2 * block_sizing(budget).per_queue
            + 2 * READ_CHUNK as u64
            + WRITE_BUF as u64,
        buffer: 0,
    }
}

/// What a record conversion writes.
pub(crate) struct Spec<'a, E> {
    pub dst: &'a Path,
    /// The format `dst` holds.
    pub format: Format,
    pub level: i32,
    pub header: Vec<u8>,
    pub framing: Framing,
    pub encoder: Arc<E>,
    /// The index to build, its reference count and its path.
    pub index: Option<(IndexFormat, usize, PathBuf)>,
    /// Text after the header, its first line number, and progress counters.
    pub text: Box<dyn BufRead + Send>,
    pub first_line: u64,
    pub read: Arc<AtomicU64>,
    pub file_len: u64,
}

pub(crate) fn run<E: LineEncoder>(ctx: &Context<'_>, spec: Spec<'_, E>) -> Result<Report> {
    let mut outputs = Outputs::new();
    let file = outputs.create(spec.dst, spec.format)?;
    let blocks = block_sizing(ctx.budget);
    let mut writer = Writer::new(
        BufWriter::with_capacity(WRITE_BUF, file),
        spec.level,
        Box::new(ctx.pool.ordered(blocks.limit)),
    )?;
    writer.write(&spec.header)?;
    writer.flush_block()?;
    let mut index = spec
        .index
        .as_ref()
        .map(|(format, n_refs, _)| BinningIndex::new(*format, *n_refs, writer.tell()));

    let s = batch_sizing(ctx.budget);
    let mut batches = LineBatches::new(spec.text, text_batch(ctx.budget), spec.first_line);
    let mut exec = ctx.pool.ordered::<Encoded>(s.limit);
    let mut records = 0u64;
    let mut take =
        |e: Encoded, writer: &mut Writer<_>, index: &mut Option<BinningIndex>| -> Result<()> {
            let mut at = 0usize;
            for &(len, entry) in &e.records {
                let rec = &e.bytes[at..at + len as usize];
                at += len as usize;
                if spec.framing == Framing::Record {
                    writer.flush_try(rec.len())?;
                    if let Some(ix) = index.as_mut() {
                        ix.amend_last(writer.tell());
                    }
                }
                writer.write(rec)?;
                if let Some(ix) = index.as_mut() {
                    ix.push(entry.tid, entry.beg, entry.end, writer.tell(), entry.mapped)?;
                }
            }
            records += e.records.len() as u64;
            Ok(())
        };
    while let Some(batch) = batches.next_batch()? {
        while exec.full() {
            take(exec.next(true)?.expect("pending"), &mut writer, &mut index)?;
        }
        let enc = Arc::clone(&spec.encoder);
        exec.submit(Box::new(move || encode_batch(&*enc, &batch)))?;
        while let Some(e) = exec.next(false)? {
            take(e, &mut writer, &mut index)?;
        }
        ctx.progress("encode", spec.read.load(Ordering::Relaxed), spec.file_len);
    }
    while let Some(e) = exec.next(true)? {
        take(e, &mut writer, &mut index)?;
    }
    writer.flush_block()?;
    let end = writer.tell();
    let done = writer.finish()?;
    done.out.into_inner().map_err(|e| e.into_error())?;
    if let (Some(ix), Some((kind, _, path))) = (index, &spec.index) {
        let bytes = ix.finish(end, &done.blocks)?;
        let format = match kind {
            IndexFormat::Bai => Format::Bai,
            IndexFormat::Csi { .. } => Format::Csi,
        };
        outputs.write(path, format, &bytes)?;
    }
    Ok(Report {
        files: outputs.commit()?,
        records,
        bytes_read: spec.read.load(Ordering::Relaxed),
        spilled: 0,
        elapsed: Default::default(),
        warnings: Vec::new(),
    })
}
