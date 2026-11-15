//! bedGraph and wig → bigWig (+ sidecar).
//!
//! Both read text the way UCSC `bedGraphToBigWig` and `wigToBigWig` do and
//! write the items with [`BigWigWriter`]; values are stored as given.
//!
//! Wig is stateful: a `fixedStep` or `variableStep` declaration governs the
//! lines after it, and a batch of lines may begin inside a section declared
//! in an earlier batch. Batches are still parsed in parallel: lines before a
//! batch's first declaration (or bedGraph line) are kept as a *lead* whose
//! chromosome and positions are not yet known, and the serial gather
//! resolves it with the section the previous batch ended in. Every line is
//! decoded once.

use std::collections::HashMap;
use std::io::{BufWriter, Read};
use std::path::{Path, PathBuf};
use std::sync::Arc;
use std::sync::atomic::Ordering;

use omicsweave_bbi::{BigWigWriter, Chrom, WriteOptions};
use omicsweave_core::exec::{Exec, Ordered};
use omicsweave_core::{Error, Result};

use crate::bed::{ChromSizes, read_chrom_sizes};
use crate::budget::{Budget, MemoryPlan, Sizing, WRITE_BUF, sizing, text_batch};
use crate::conversion::{Context, Report};
use crate::format::Format;
use crate::pipeline::input::{LineBatch, LineBatches, READ_CHUNK, lines, open_text};
use crate::pipeline::output::Outputs;
use crate::pipeline::sort::{Sorter, Spill, read_or_end};

/// bedGraph (plain, gzip or BGZF) to bigWig, with its sidecar.
///
/// Every data line is `chrom start end value`, whitespace-separated, with
/// 0-based half-open coordinates; `track`, `browser`, `#` and blank lines
/// are skipped. Lines need not be sorted: items are sorted by chromosome
/// (in order of first appearance) and start, in memory or out of core
/// within the budget. Overlapping items, empty items and values that are
/// not finite 32-bit floats are errors.
#[derive(Debug, Clone, PartialEq)]
pub struct BedGraphToBigWig {
    pub src: PathBuf,
    pub dst: PathBuf,
    pub chrom_sizes: ChromSizes,
    pub options: WriteOptions,
    /// Also write `<dst>.sidecar`.
    pub sidecar: bool,
}

/// Wig (plain, gzip or BGZF) to bigWig, with its sidecar.
///
/// Accepts what `wigToBigWig` does: `fixedStep chrom= start= [step=]
/// [span=]` sections of one value per line, `variableStep chrom= [span=]`
/// sections of `position value` lines (positions 1-based), and bedGraph
/// lines, which end a step section. `step` and `span` default to 1. Items
/// are sorted and checked as in [`BedGraphToBigWig`]; a `span` larger
/// than the `step` makes items overlap and is an error.
#[derive(Debug, Clone, PartialEq)]
pub struct WigToBigWig {
    pub src: PathBuf,
    pub dst: PathBuf,
    pub chrom_sizes: ChromSizes,
    pub options: WriteOptions,
    /// Also write `<dst>.sidecar`.
    pub sidecar: bool,
}

impl BedGraphToBigWig {
    pub fn new(
        src: impl Into<PathBuf>,
        dst: impl Into<PathBuf>,
        chrom_sizes: ChromSizes,
    ) -> BedGraphToBigWig {
        BedGraphToBigWig {
            src: src.into(),
            dst: dst.into(),
            chrom_sizes,
            options: WriteOptions::default(),
            sidecar: true,
        }
    }

    pub(crate) fn plan(&self, budget: &Budget) -> Result<MemoryPlan> {
        plan(&self.chrom_sizes, Dialect::BedGraph, budget)
    }

    pub(crate) fn run(&self, ctx: &Context<'_>) -> Result<Report> {
        run(
            Job {
                src: &self.src,
                dst: &self.dst,
                chrom_sizes: &self.chrom_sizes,
                options: &self.options,
                sidecar: self.sidecar,
                dialect: Dialect::BedGraph,
            },
            ctx,
        )
    }
}

impl WigToBigWig {
    pub fn new(
        src: impl Into<PathBuf>,
        dst: impl Into<PathBuf>,
        chrom_sizes: ChromSizes,
    ) -> WigToBigWig {
        WigToBigWig {
            src: src.into(),
            dst: dst.into(),
            chrom_sizes,
            options: WriteOptions::default(),
            sidecar: true,
        }
    }

    pub(crate) fn plan(&self, budget: &Budget) -> Result<MemoryPlan> {
        plan(&self.chrom_sizes, Dialect::Wig, budget)
    }

    pub(crate) fn run(&self, ctx: &Context<'_>) -> Result<Report> {
        run(
            Job {
                src: &self.src,
                dst: &self.dst,
                chrom_sizes: &self.chrom_sizes,
                options: &self.options,
                sidecar: self.sidecar,
                dialect: Dialect::Wig,
            },
            ctx,
        )
    }
}

/// Which lines a source may hold.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum Dialect {
    /// bedGraph lines only.
    BedGraph,
    /// Step declarations and their lines, and bedGraph lines.
    Wig,
}

struct Job<'a> {
    src: &'a PathBuf,
    dst: &'a PathBuf,
    chrom_sizes: &'a ChromSizes,
    options: &'a WriteOptions,
    sidecar: bool,
    dialect: Dialect,
}

fn parse_sizing(budget: &Budget, dialect: Dialect) -> Sizing {
    // A batch in flight holds its text and what it parses into: at most
    // 2.5 times the text for bedGraph lines (a 20-byte item per line of at
    // least 8 bytes), 5 times for wig (`variableStep` lines of 4 bytes).
    let factor = match dialect {
        Dialect::BedGraph => 4,
        Dialect::Wig => 6,
    };
    sizing(
        budget.threads,
        factor * text_batch(budget) as u64,
        budget.memory / 4,
    )
}

fn block_sizing(budget: &Budget) -> Sizing {
    sizing(budget.threads, 256 << 10, budget.memory / 16)
}

/// Parse batches in flight for a text of signal items read as bedGraph lines.
pub(crate) fn bedgraph_sizing(budget: &Budget) -> Sizing {
    parse_sizing(budget, Dialect::BedGraph)
}

/// Memory of a text of signal items as bedGraph lines plan it: parse
/// batches in flight, blocks in flight, input and output buffers, and the
/// sort buffer.
pub(crate) fn plan_bedgraph(chrom_sizes: &ChromSizes, budget: &Budget) -> Result<MemoryPlan> {
    plan(chrom_sizes, Dialect::BedGraph, budget)
}

fn plan(chrom_sizes: &ChromSizes, dialect: Dialect, budget: &Budget) -> Result<MemoryPlan> {
    if let ChromSizes::File(p) = chrom_sizes
        && !p.exists()
    {
        return Err(Error::Io {
            kind: std::io::ErrorKind::NotFound,
            detail: format!("chromosome sizes {} do not exist", p.display()),
        });
    }
    let pipeline = parse_sizing(budget, dialect).per_queue
        + block_sizing(budget).per_queue
        + 2 * READ_CHUNK as u64
        + WRITE_BUF as u64;
    Ok(MemoryPlan {
        fixed: 0,
        pipeline,
        buffer: budget.memory.saturating_sub(pipeline),
    })
}

fn run(job: Job<'_>, ctx: &Context<'_>) -> Result<Report> {
    let sizes = match job.chrom_sizes {
        ChromSizes::File(p) => Some(read_chrom_sizes(p)?),
        ChromSizes::FromData => None,
    };
    let s = parse_sizing(ctx.budget, job.dialect);
    let input = open_text(job.src, &ctx.pool, s.limit)?;
    let (read, file_len) = (Arc::clone(&input.read), input.file_len);
    let mut batches = LineBatches::new(input.text, text_batch(ctx.budget), 1);
    let mut exec = ctx.pool.ordered::<Parsed>(s.limit);
    let tmp = ctx.tmp_parent(job.dst);
    let mut gather = Gather {
        sorter: Sorter::new(&ctx.pool, ctx.memory.buffer, &tmp),
        chroms: Registry::new(sizes.as_ref()),
        section: Section::None,
        items: 0,
    };
    let dialect = job.dialect;
    while let Some(batch) = batches.next_batch()? {
        while exec.full() {
            gather.take(exec.next(true)?.expect("pending"))?;
        }
        exec.submit(Box::new(move || parse(&batch, dialect)))?;
        while let Some(p) = exec.next(false)? {
            gather.take(p)?;
        }
        ctx.progress("read", read.load(Ordering::Relaxed), file_len);
    }
    while let Some(p) = exec.next(true)? {
        gather.take(p)?;
    }
    let records = gather.items;
    let spilled = gather.sorter.spilled();
    let (names, chroms) = gather.chroms.into_chroms("bigwig")?;
    ctx.progress("sort", read.load(Ordering::Relaxed), file_len);
    let mut merged = gather.sorter.finish()?;
    let mut last: Option<WigRec> = None;
    let written = write_bigwig(
        ctx,
        job.dst,
        chroms,
        job.options,
        job.sidecar,
        records,
        move || {
            let Some(r) = merged.next_record()? else {
                return Ok(None);
            };
            if let Some(p) = &last
                && p.chrom == r.chrom
                && r.start < p.end
            {
                return Err(Error::Format {
                    line: r.line,
                    detail: format!(
                        "{}:{}-{} overlaps {}-{} of line {}",
                        names[r.chrom as usize], r.start, r.end, p.start, p.end, p.line
                    ),
                });
            }
            last = Some(r);
            Ok(Some((r.chrom, r.start, r.end, r.value)))
        },
    )?;
    Ok(Report {
        files: written,
        records,
        bytes_read: read.load(Ordering::Relaxed),
        spilled,
        elapsed: Default::default(),
        warnings: Vec::new(),
    })
}

/// Writes the items `next` yields, sorted by (chromosome, start) and
/// checked, as `dst` (and its `.sidecar`) with blocks compressed on the
/// pool; returns the files once published. `total` items are expected, for
/// progress.
pub(crate) fn write_bigwig(
    ctx: &Context<'_>,
    dst: &Path,
    chroms: Vec<Chrom>,
    options: &WriteOptions,
    sidecar: bool,
    total: u64,
    mut next: impl FnMut() -> Result<Option<(u32, u32, u32, f32)>>,
) -> Result<Vec<(Format, PathBuf)>> {
    let mut outputs = Outputs::new();
    let file = outputs.create(dst, Format::BigWig)?;
    let b = block_sizing(ctx.budget);
    let mut writer = BigWigWriter::create_with(
        BufWriter::with_capacity(WRITE_BUF, file),
        chroms,
        options.clone(),
        Exec::Pool {
            pool: &ctx.pool,
            limit: b.limit,
        },
    )?;
    let mut n = 0u64;
    while let Some((chrom, start, end, value)) = next()? {
        writer.add(chrom, start, end, value)?;
        n += 1;
        if n.is_multiple_of(1 << 16) {
            ctx.progress("write", n, total);
        }
    }
    // The sorted stream's buffers go before the zoom levels are finished.
    drop(next);
    let (out, built) = writer.finish()?;
    out.into_inner().map_err(|e| e.into_error())?;
    if sidecar {
        let source = outputs.source(dst)?;
        let file = outputs.create(
            &omicsweave_core::sidecar::sidecar_path(dst),
            Format::Sidecar,
        )?;
        built.write_to(&source, file)?;
    }
    outputs.commit()
}

/// The section a line belongs to, with its chromosome named.
#[derive(Debug, Clone, PartialEq, Eq)]
enum Section {
    /// Outside any step section: only bedGraph lines.
    None,
    /// `fixedStep`: the next value covers `[next, next + span)`.
    Fixed {
        chrom: String,
        next: u64,
        step: u32,
        span: u32,
    },
    /// `variableStep`: a line at position `p` covers `[p - 1, p - 1 + span)`.
    Variable { chrom: String, span: u32 },
}

/// An item of a batch: chromosome local to the batch, line relative to
/// the batch's first.
#[derive(Debug, Clone, Copy)]
struct Item {
    chrom: u32,
    start: u32,
    end: u32,
    value: f32,
    line: u32,
}

/// Consecutive `fixedStep` values: the `i`-th of `values` covers
/// `[start + i * step, start + i * step + span)`.
struct Run {
    chrom: u32,
    start: u64,
    step: u32,
    span: u32,
    values: std::ops::Range<usize>,
}

/// What one batch parses into.
struct Parsed {
    first_line: u64,
    /// Names of the local chromosome ids, in order of first item.
    chroms: Vec<String>,
    /// Lines before the batch's first declaration or bedGraph line, which
    /// continue the section the previous batch ended in: `fixedStep` values
    /// (the first `lead_fixed` of `values`) or `variableStep` items (start
    /// set, chromosome and end not yet known). At most one kind is present.
    lead_fixed: usize,
    lead_variable: Vec<Item>,
    runs: Vec<Run>,
    /// (value, relative line) of every `fixedStep` line.
    values: Vec<(f32, u32)>,
    items: Vec<Item>,
    /// The section the batch ends in; `None` when the whole batch is lead.
    end: Option<Section>,
}

/// The section state while parsing a batch; chromosomes as slices of the
/// batch, numbered once they have data.
enum Local<'a> {
    None,
    Fixed {
        name: &'a [u8],
        id: Option<u32>,
        next: u64,
        step: u32,
        span: u32,
        run: Option<usize>,
    },
    Variable {
        name: &'a [u8],
        id: Option<u32>,
        span: u32,
    },
}

/// Local chromosome ids of a batch, in order of first item.
#[derive(Default)]
struct LocalChroms<'a> {
    ids: HashMap<&'a [u8], u32>,
    names: Vec<String>,
}

impl<'a> LocalChroms<'a> {
    fn id(&mut self, name: &'a [u8]) -> u32 {
        let next = self.names.len() as u32;
        *self.ids.entry(name).or_insert_with(|| {
            self.names.push(String::from_utf8_lossy(name).into_owned());
            next
        })
    }

    fn cached(&mut self, name: &'a [u8], id: &mut Option<u32>) -> u32 {
        *id.get_or_insert_with(|| self.id(name))
    }
}

fn bad(line: u64, detail: String) -> Error {
    Error::Format { line, detail }
}

/// The line without its `\r` and surrounding whitespace.
pub(crate) fn trim(line: &[u8]) -> &[u8] {
    let line = line.strip_suffix(b"\r").unwrap_or(line);
    let start = line
        .iter()
        .position(|c| !c.is_ascii_whitespace())
        .unwrap_or(line.len());
    let end = line
        .iter()
        .rposition(|c| !c.is_ascii_whitespace())
        .map_or(start, |i| i + 1);
    &line[start..end]
}

/// Whether `line` starts with the word `word`.
pub(crate) fn starts_with_word(line: &[u8], word: &[u8]) -> bool {
    line.starts_with(word) && line.get(word.len()).is_none_or(|c| c.is_ascii_whitespace())
}

/// Whether a trimmed line carries no data: blank, a comment, a `track` or
/// `browser` line.
pub(crate) fn is_meta(line: &[u8]) -> bool {
    line.is_empty()
        || line[0] == b'#'
        || starts_with_word(line, b"track")
        || starts_with_word(line, b"browser")
}

fn words(line: &[u8]) -> impl Iterator<Item = &[u8]> {
    line.split(|c| c.is_ascii_whitespace())
        .filter(|w| !w.is_empty())
}

fn text(field: &[u8]) -> std::borrow::Cow<'_, str> {
    String::from_utf8_lossy(field)
}

fn coordinate(field: &[u8], line: u64, what: &str) -> Result<u32> {
    std::str::from_utf8(field)
        .ok()
        .and_then(|s| s.parse::<u32>().ok())
        .ok_or_else(|| {
            bad(
                line,
                format!("{what} {:?} is not a coordinate", text(field)),
            )
        })
}

/// A value as `bedGraphToBigWig` reads it: a double, stored as a float.
pub(crate) fn value(field: &[u8], line: u64) -> Result<f32> {
    let v = std::str::from_utf8(field)
        .ok()
        .and_then(|s| s.parse::<f64>().ok())
        .ok_or_else(|| bad(line, format!("value {:?} is not a number", text(field))))?;
    let f = v as f32;
    if !f.is_finite() {
        return Err(bad(
            line,
            format!("value {:?} is not a finite 32-bit float", text(field)),
        ));
    }
    Ok(f)
}

fn setting(field: &[u8], line: u64, what: &str) -> Result<u32> {
    match std::str::from_utf8(field)
        .ok()
        .and_then(|s| s.parse::<u32>().ok())
    {
        Some(v) if v > 0 => Ok(v),
        _ => Err(bad(
            line,
            format!("{what} {:?} is not a positive integer", text(field)),
        )),
    }
}

/// Parses a `fixedStep` / `variableStep` declaration.
fn declaration(line: &[u8], no: u64) -> Result<Local<'_>> {
    let mut w = words(line);
    let fixed = w.next() == Some(b"fixedStep".as_slice());
    let kind = if fixed { "fixedStep" } else { "variableStep" };
    let (mut name, mut start, mut step, mut span) = (None, None, 1u32, 1u32);
    for word in w {
        let Some(eq) = memchr::memchr(b'=', word) else {
            return Err(bad(
                no,
                format!("{kind}: {:?} is not a setting", text(word)),
            ));
        };
        let (key, val) = (&word[..eq], &word[eq + 1..]);
        match key {
            b"chrom" if !val.is_empty() => name = Some(val),
            b"span" => span = setting(val, no, "span")?,
            b"start" if fixed => start = Some(setting(val, no, "start")?),
            b"step" if fixed => step = setting(val, no, "step")?,
            _ => return Err(bad(no, format!("{kind}: unknown setting {:?}", text(word)))),
        }
    }
    let Some(name) = name else {
        return Err(bad(no, format!("{kind} without chrom=")));
    };
    if !fixed {
        return Ok(Local::Variable {
            name,
            id: None,
            span,
        });
    }
    let Some(start) = start else {
        return Err(bad(no, "fixedStep without start=".into()));
    };
    Ok(Local::Fixed {
        name,
        id: None,
        next: start as u64 - 1,
        step,
        span,
        run: None,
    })
}

/// `[start, start + span)`, checked to fit coordinates.
fn span_end(start: u64, span: u32, line: u64) -> Result<(u32, u32)> {
    let end = start + span as u64;
    if end > u32::MAX as u64 {
        return Err(bad(
            line,
            format!("item {start}-{end} is past the largest coordinate"),
        ));
    }
    Ok((start as u32, end as u32))
}

fn bedgraph_fields(n: usize) -> String {
    format!(
        "{n} field{}; bedGraph lines have 4 (chrom start end value)",
        if n == 1 { "" } else { "s" }
    )
}

fn parse(batch: &LineBatch, dialect: Dialect) -> Result<Parsed> {
    let mut out = Parsed {
        first_line: batch.first_line,
        chroms: Vec::new(),
        lead_fixed: 0,
        lead_variable: Vec::new(),
        runs: Vec::new(),
        values: Vec::new(),
        items: Vec::new(),
        end: None,
    };
    let mut local = LocalChroms::default();
    // `None` while the batch's lines continue an earlier section.
    let mut state: Option<Local<'_>> = match dialect {
        Dialect::BedGraph => Some(Local::None),
        Dialect::Wig => None,
    };
    let mut fields: [&[u8]; 4] = [&[]; 4];
    for (no, raw) in lines(batch) {
        let line = trim(raw);
        if is_meta(line) {
            continue;
        }
        let rel = (no - batch.first_line) as u32;
        if dialect == Dialect::Wig
            && (starts_with_word(line, b"fixedStep") || starts_with_word(line, b"variableStep"))
        {
            state = Some(declaration(line, no)?);
            continue;
        }
        let mut n = 0;
        for w in words(line) {
            if n < 4 {
                fields[n] = w;
            }
            n += 1;
            if n > 4 {
                break;
            }
        }
        match (&mut state, n) {
            (
                Some(Local::Fixed {
                    name,
                    id,
                    next,
                    step,
                    span,
                    run,
                }),
                1,
            ) => {
                let v = value(fields[0], no)?;
                span_end(*next, *span, no)?;
                let at = out.values.len();
                match run {
                    Some(r) => out.runs[*r].values.end = at + 1,
                    None => {
                        *run = Some(out.runs.len());
                        out.runs.push(Run {
                            chrom: local.cached(name, id),
                            start: *next,
                            step: *step,
                            span: *span,
                            values: at..at + 1,
                        });
                    }
                }
                out.values.push((v, rel));
                *next += *step as u64;
            }
            (Some(Local::Variable { name, id, span }), 2) => {
                let pos = setting(fields[0], no, "position")?;
                let v = value(fields[1], no)?;
                let (start, end) = span_end(pos as u64 - 1, *span, no)?;
                out.items.push(Item {
                    chrom: local.cached(name, id),
                    start,
                    end,
                    value: v,
                    line: rel,
                });
            }
            (Some(Local::Variable { .. }), 1) => {
                return Err(bad(
                    no,
                    "1 field; variableStep lines have 2 (position value)".into(),
                ));
            }
            (_, 4) => {
                let start = coordinate(fields[1], no, "start")?;
                let end = coordinate(fields[2], no, "end")?;
                if start >= end {
                    return Err(bad(
                        no,
                        if start == end {
                            format!("empty item {start}-{end}")
                        } else {
                            format!("start {start} after end {end}")
                        },
                    ));
                }
                let v = value(fields[3], no)?;
                out.items.push(Item {
                    chrom: local.id(fields[0]),
                    start,
                    end,
                    value: v,
                    line: rel,
                });
                state = Some(Local::None);
            }
            (None, 1) if out.lead_variable.is_empty() => {
                out.values.push((value(fields[0], no)?, rel));
                out.lead_fixed += 1;
            }
            (None, 2) if out.lead_fixed == 0 => {
                let pos = setting(fields[0], no, "position")?;
                out.lead_variable.push(Item {
                    chrom: 0,
                    start: pos - 1,
                    end: 0,
                    value: value(fields[1], no)?,
                    line: rel,
                });
            }
            (None, n @ (1 | 2)) => {
                return Err(bad(
                    no,
                    format!("{n} fields among lines of {} in one section", 3 - n),
                ));
            }
            (_, n) => return Err(bad(no, bedgraph_fields(n))),
        }
    }
    out.end = state.map(|s| match s {
        Local::None => Section::None,
        Local::Fixed {
            name,
            next,
            step,
            span,
            ..
        } => Section::Fixed {
            chrom: text(name).into_owned(),
            next,
            step,
            span,
        },
        Local::Variable { name, span, .. } => Section::Variable {
            chrom: text(name).into_owned(),
            span,
        },
    });
    out.chroms = local.names;
    Ok(out)
}

/// An item in the sort buffer; `line` breaks ties in input order.
#[derive(Debug, Clone, Copy)]
pub(crate) struct WigRec {
    pub chrom: u32,
    pub start: u32,
    pub end: u32,
    pub value: f32,
    pub line: u64,
}

impl Spill for WigRec {
    type Key = (u32, u32, u64);
    fn key(&self) -> Self::Key {
        (self.chrom, self.start, self.line)
    }
    fn mem_size(&self) -> usize {
        0
    }
    fn encode(&self, out: &mut Vec<u8>) {
        out.extend_from_slice(&self.chrom.to_le_bytes());
        out.extend_from_slice(&self.start.to_le_bytes());
        out.extend_from_slice(&self.end.to_le_bytes());
        out.extend_from_slice(&self.value.to_le_bytes());
        out.extend_from_slice(&self.line.to_le_bytes());
    }
    fn decode(input: &mut impl Read) -> Result<Option<WigRec>> {
        let mut h = [0u8; 24];
        if !read_or_end(input, &mut h)? {
            return Ok(None);
        }
        let u = |i: usize| u32::from_le_bytes(h[i..i + 4].try_into().unwrap());
        Ok(Some(WigRec {
            chrom: u(0),
            start: u(4),
            end: u(8),
            value: f32::from_bits(u(12)),
            line: u64::from_le_bytes(h[16..24].try_into().unwrap()),
        }))
    }
}

/// The output's chromosomes, numbered in order of first item.
pub(crate) struct Registry<'s> {
    ids: HashMap<String, u32>,
    names: Vec<String>,
    max_end: Vec<u32>,
    sizes: Option<&'s HashMap<String, u32>>,
}

impl<'s> Registry<'s> {
    /// No chromosomes yet; with `sizes`, only those are known.
    pub(crate) fn new(sizes: Option<&'s HashMap<String, u32>>) -> Registry<'s> {
        Registry {
            ids: HashMap::new(),
            names: Vec::new(),
            max_end: Vec::new(),
            sizes,
        }
    }

    /// The names, and the chromosomes to write: each with its size, or its
    /// largest item end; none is an error of `what`.
    pub(crate) fn into_chroms(self, what: &'static str) -> Result<(Vec<String>, Vec<Chrom>)> {
        let chroms: Vec<Chrom> = self
            .names
            .iter()
            .zip(&self.max_end)
            .map(|(name, &end)| Chrom {
                name: name.clone(),
                len: self.sizes.map_or(end, |s| s[name]),
            })
            .collect();
        if chroms.is_empty() {
            return Err(Error::invalid(
                what,
                "no data; a bigWig needs at least one chromosome",
            ));
        }
        Ok((self.names, chroms))
    }

    pub(crate) fn id(&mut self, name: &str) -> Result<u32> {
        if let Some(&id) = self.ids.get(name) {
            return Ok(id);
        }
        if let Some(s) = self.sizes
            && !s.contains_key(name)
        {
            return Err(Error::invalid(
                "bigwig",
                format!("chromosome {name:?} is not in the chromosome sizes"),
            ));
        }
        let id = self.names.len() as u32;
        self.ids.insert(name.to_owned(), id);
        self.names.push(name.to_owned());
        self.max_end.push(0);
        Ok(id)
    }

    /// Checks that an item ending at `end` fits its chromosome.
    pub(crate) fn cover(&mut self, id: u32, end: u32, line: u64) -> Result<()> {
        if let Some(s) = self.sizes {
            let name = &self.names[id as usize];
            let len = s[name];
            if end > len {
                return Err(bad(
                    line,
                    format!("end {end} past the end of {name} ({len})"),
                ));
            }
        }
        let m = &mut self.max_end[id as usize];
        *m = (*m).max(end);
        Ok(())
    }
}

struct Gather<'p, 's> {
    sorter: Sorter<'p, WigRec>,
    chroms: Registry<'s>,
    /// The section the batches taken so far end in.
    section: Section,
    items: u64,
}

impl Gather<'_, '_> {
    fn take(&mut self, p: Parsed) -> Result<()> {
        let mut out = Vec::with_capacity(p.values.len() + p.items.len() + p.lead_variable.len());
        let line = |rel: u32| p.first_line + rel as u64;
        let mut push = |chroms: &mut Registry<'_>, r: WigRec| -> Result<()> {
            chroms.cover(r.chrom, r.end, r.line)?;
            out.push(r);
            Ok(())
        };
        if p.lead_fixed > 0 {
            let first = line(p.values[0].1);
            match &mut self.section {
                Section::Fixed {
                    chrom,
                    next,
                    step,
                    span,
                } => {
                    let id = self.chroms.id(chrom)?;
                    for &(value, rel) in &p.values[..p.lead_fixed] {
                        let (start, end) = span_end(*next, *span, line(rel))?;
                        push(
                            &mut self.chroms,
                            WigRec {
                                chrom: id,
                                start,
                                end,
                                value,
                                line: line(rel),
                            },
                        )?;
                        *next += *step as u64;
                    }
                }
                Section::Variable { .. } => {
                    return Err(bad(
                        first,
                        "1 field; variableStep lines have 2 (position value)".into(),
                    ));
                }
                Section::None => return Err(bad(first, bedgraph_fields(1))),
            }
        }
        if let Some(it) = p.lead_variable.first() {
            let Section::Variable { chrom, span } = &self.section else {
                return Err(bad(line(it.line), bedgraph_fields(2)));
            };
            let id = self.chroms.id(chrom)?;
            for it in &p.lead_variable {
                let (start, end) = span_end(it.start as u64, *span, line(it.line))?;
                push(
                    &mut self.chroms,
                    WigRec {
                        chrom: id,
                        start,
                        end,
                        value: it.value,
                        line: line(it.line),
                    },
                )?;
            }
        }
        let mut map = Vec::with_capacity(p.chroms.len());
        for name in &p.chroms {
            map.push(self.chroms.id(name)?);
        }
        for run in &p.runs {
            let id = map[run.chrom as usize];
            for (i, &(value, rel)) in p.values[run.values.clone()].iter().enumerate() {
                let (start, end) =
                    span_end(run.start + i as u64 * run.step as u64, run.span, line(rel))?;
                push(
                    &mut self.chroms,
                    WigRec {
                        chrom: id,
                        start,
                        end,
                        value,
                        line: line(rel),
                    },
                )?;
            }
        }
        for it in &p.items {
            push(
                &mut self.chroms,
                WigRec {
                    chrom: map[it.chrom as usize],
                    start: it.start,
                    end: it.end,
                    value: it.value,
                    line: line(it.line),
                },
            )?;
        }
        if let Some(end) = p.end {
            self.section = end;
        }
        self.items += out.len() as u64;
        self.sorter.extend(out)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn batch(text: &str, first_line: u64) -> LineBatch {
        LineBatch {
            text: text.as_bytes().to_vec(),
            first_line,
        }
    }

    #[test]
    fn a_batch_inside_a_section_keeps_its_lead() {
        let p = parse(&batch("1.5\n# c\n2\nchr2 0 5 1\n", 10), Dialect::Wig).unwrap();
        assert_eq!(p.lead_fixed, 2);
        assert_eq!(p.values, [(1.5, 0), (2.0, 2)]);
        assert_eq!(p.items.len(), 1);
        assert_eq!(p.end, Some(Section::None));

        let p = parse(&batch("5 1\n9 2\n", 1), Dialect::Wig).unwrap();
        assert_eq!(p.lead_variable.len(), 2);
        assert_eq!(p.end, None);
    }

    #[test]
    fn declarations_set_the_section_the_batch_ends_in() {
        let text = "fixedStep chrom=chr1 start=11 step=5 span=2\n1\n2\nvariableStep chrom=chr3\n";
        let p = parse(&batch(text, 1), Dialect::Wig).unwrap();
        assert_eq!(p.runs.len(), 1);
        assert_eq!(
            (p.runs[0].start, p.runs[0].step, p.runs[0].span),
            (10, 5, 2)
        );
        assert_eq!(p.chroms, ["chr1"]);
        assert_eq!(
            p.end,
            Some(Section::Variable {
                chrom: "chr3".into(),
                span: 1
            })
        );
    }

    #[test]
    fn bad_lines_are_reported_with_their_numbers() {
        for (dialect, text, line) in [
            (Dialect::BedGraph, "chr1 0 10 1\nchr1 5 2\n", 2),
            (Dialect::BedGraph, "fixedStep chrom=chr1 start=1\n", 1),
            (Dialect::BedGraph, "chr1 5 5 1\n", 1),
            (Dialect::BedGraph, "chr1 0 5 nan\n", 1),
            (Dialect::Wig, "fixedStep chrom=chr1\n", 1),
            (Dialect::Wig, "variableStep chrom=chr1 step=2\n1 1\n", 1),
            (Dialect::Wig, "variableStep chrom=chr1\n1 1\n7\n", 3),
            (Dialect::Wig, "fixedStep chrom=c start=1\n1\n2 3\n", 3),
            (Dialect::Wig, "1\n2 3\n", 2),
        ] {
            let Err(Error::Format { line: got, .. }) = parse(&batch(text, 1), dialect) else {
                panic!("{text:?} parses");
            };
            assert_eq!(got, line, "{text:?}");
        }
    }
}
