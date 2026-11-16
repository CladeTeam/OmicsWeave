//! 4DN pairs, Juicer short and HiC-Pro text → `.hic` (+ sidecar).
//!
//! Each row is one contact, binned as [`omicsweave_hic::HicBuilder::add`]
//! (and the Juicer tools' `pre`) bins it: the two ends ordered by
//! chromosome index, positions taken as written and clamped to the
//! chromosome's length, `bin = position / resolution`, and within a
//! chromosome the cell moved to the upper triangle.
//!
//! The rows are never held whole:
//!
//! 1. serial: the header, then batches of lines (BGZF inflated on the pool);
//! 2. parallel: every row parsed, binned at the finest resolution and placed
//!    in its block;
//! 3. an external sort by (chromosome pair, block, y, x) that folds the
//!    contacts of one cell into a pixel;
//! 4. the finest pixels, streamed once, are kept in a spill file and, in
//!    chunks of a fixed size, summed on the pool into every coarser
//!    resolution, each sorted on its own;
//! 5. maps written coarsest first through the stage shared with cool → hic.
//!
//! Counts are summed exactly (integers in `f64`), so no fold depends on
//! where a sort spilled. Scores of the short format may be fractional:
//! their contacts are then kept apart by row (finest) and chunk (coarser)
//! until the last merge, so every sum is taken in one order, whatever the
//! budget.

use std::collections::{BTreeMap, HashMap};
use std::fs::File;
use std::io::{BufRead, BufReader, BufWriter, Read, Write};
use std::path::PathBuf;
use std::sync::Arc;
use std::sync::atomic::Ordering;

use omicsweave_core::exec::Ordered;
use omicsweave_core::{ContactUnit, Error, Result};
use omicsweave_hic::{Chrom, HicHeader, MatrixLayout, MatrixLayouts};

use crate::bed::read_chrom_list;
use crate::budget::{Budget, MemoryPlan, Sizing, WRITE_BUF, sizing, text_batch};
use crate::conversion::{Context, Report};
use crate::pipeline::hic::{HicSink, Pixel, block_sizing};
use crate::pipeline::input::{LineBatch, LineBatches, READ_CHUNK, lines, open_text, with_path};
use crate::pipeline::output::TempDir;
use crate::pipeline::sort::{Sorter, Spill, read_or_end};

/// The base-pair resolutions the Juicer tools' `pre` bins at by default.
pub const JUICER_RESOLUTIONS: [u32; 9] = [
    2_500_000, 1_000_000, 500_000, 250_000, 100_000, 50_000, 25_000, 10_000, 5_000,
];

/// Which contact text a source holds.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum PairsDialect {
    /// 4DN pairs: `#` header lines (`#chromsize:`, `#columns:`), then one
    /// row per pair with `chr1 pos1 chr2 pos2` among its columns.
    Pairs,
    /// Juicer short: `str1 chr1 pos1 frag1 str2 chr2 pos2 frag2`, with an
    /// optional ninth column, the score, which is then the row's count.
    Short,
    /// HiC-Pro validPairs/allValidPairs: read name, two chromosome/position/
    /// strand triples, fragment information and mapping qualities.
    HicPro,
}

/// 4DN pairs, Juicer short or HiC-Pro text (plain, gzip or BGZF; a pairix index is
/// not read) to Juicer `.hic` (v9) at base-pair resolutions, with its
/// sidecar.
///
/// Every row counts once (in short-with-score text, its score); nothing is
/// filtered (no MAPQ, pair type or same-fragment rule). Rows naming a
/// chromosome outside the chromosome list are skipped and reported as a
/// warning with their count per chromosome. The chromosomes, in order, come
/// from `chrom_sizes` or else the pairs header's `#chromsize:` lines; short
/// and HiC-Pro text have no header and need `chrom_sizes`. Coarser resolutions are sums
/// of the finest, so each must be a multiple of the finest. Fragment
/// resolutions (which need restriction-site files) are not written.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct PairsToHic {
    pub src: PathBuf,
    pub dst: PathBuf,
    pub dialect: PairsDialect,
    /// `chrom.sizes`, `.fai` or `.2bit`: the chromosomes and their order.
    pub chrom_sizes: Option<PathBuf>,
    /// Base-pair resolutions; `None` takes [`JUICER_RESOLUTIONS`].
    pub resolutions: Option<Vec<u32>>,
    /// zlib level of blocks.
    pub level: u32,
    /// Also write `<dst>.sidecar`.
    pub sidecar: bool,
}

/// Finest pixels per coarsening task. Fixed, so that sums of fractional
/// scores are grouped the same way whatever the budget.
const CHUNK: usize = 1 << 12;

/// The 4DN pairs columns when the header has no `#columns:` line.
const PAIRS_COLUMNS: [&str; 7] = [
    "readID", "chr1", "pos1", "chr2", "pos2", "strand1", "strand2",
];

impl PairsToHic {
    pub fn new(
        src: impl Into<PathBuf>,
        dst: impl Into<PathBuf>,
        dialect: PairsDialect,
    ) -> PairsToHic {
        PairsToHic {
            src: src.into(),
            dst: dst.into(),
            dialect,
            chrom_sizes: None,
            resolutions: None,
            level: 6,
            sidecar: true,
        }
    }

    /// The resolutions, coarsest first as the Juicer tools list them.
    fn resolutions(&self) -> Result<Vec<u32>> {
        let mut r = self
            .resolutions
            .clone()
            .unwrap_or_else(|| JUICER_RESOLUTIONS.to_vec());
        if r.is_empty() {
            return Err(Error::invalid("pairs", "no resolutions given"));
        }
        if r.contains(&0) {
            return Err(Error::invalid("pairs", "a resolution of 0 base pairs"));
        }
        r.sort_unstable_by(|a, b| b.cmp(a));
        r.dedup();
        let finest = *r.last().expect("not empty");
        if let Some(bad) = r.iter().find(|&&x| x % finest != 0) {
            return Err(Error::invalid(
                "pairs",
                format!(
                    "resolution {bad} is not a multiple of the finest, {finest}: coarser maps \
                     are summed from the finest one"
                ),
            ));
        }
        Ok(r)
    }

    fn parse_sizing(budget: &Budget) -> Sizing {
        // A batch's text and its contacts (40 bytes per row of at least
        // about 16 bytes, often far more).
        sizing(
            budget.threads,
            3 * text_batch(budget) as u64,
            budget.memory / 4,
        )
    }

    fn chunk_sizing(budget: &Budget, coarse: usize) -> Sizing {
        let unit = CHUNK * (std::mem::size_of::<Fine>() + coarse * std::mem::size_of::<Contact>());
        sizing(budget.threads, unit as u64, budget.memory / 8)
    }

    pub(crate) fn plan(&self, budget: &Budget) -> Result<MemoryPlan> {
        if let Some(p) = &self.chrom_sizes
            && !p.exists()
        {
            return Err(Error::Io {
                kind: std::io::ErrorKind::NotFound,
                detail: format!("chromosome sizes {} do not exist", p.display()),
            });
        }
        let coarse = self.resolutions()?.len() - 1;
        let pipeline = Self::parse_sizing(budget).per_queue
            + Self::chunk_sizing(budget, coarse).per_queue
            + block_sizing(budget).per_queue
            + 2 * READ_CHUNK as u64
            + 2 * WRITE_BUF as u64;
        Ok(MemoryPlan {
            fixed: 0,
            pipeline,
            buffer: budget.memory.saturating_sub(pipeline),
        })
    }

    pub(crate) fn run(&self, ctx: &Context<'_>) -> Result<Report> {
        let resolutions = self.resolutions()?;
        let finest = *resolutions.last().expect("not empty");
        let s = Self::parse_sizing(ctx.budget);
        let input = open_text(&self.src, &ctx.pool, s.limit)?;
        let (read, file_len) = (Arc::clone(&input.read), input.file_len);
        let mut text = input.text;
        let head_limit = (ctx.budget.memory / 8).clamp(64 << 10, 8 << 20) as usize;
        let mut head = read_head(&mut text, self.dialect, head_limit)?;
        let chroms: Vec<(String, u64)> = match &self.chrom_sizes {
            Some(p) => read_chrom_list(p)?
                .into_iter()
                .map(|(n, l)| (n, l as u64))
                .collect(),
            None => head.chroms.clone(),
        };
        if chroms.is_empty() {
            return Err(Error::invalid(
                "pairs",
                match self.dialect {
                    PairsDialect::Pairs => format!(
                        "{} has no #chromsize lines: give the chromosome sizes",
                        self.src.display()
                    ),
                    PairsDialect::Short => {
                        "Juicer short text names no chromosome lengths: give the chromosome \
                         sizes"
                            .to_owned()
                    }
                    PairsDialect::HicPro => {
                        "HiC-Pro text names no chromosome lengths: give the chromosome sizes"
                            .to_owned()
                    }
                },
            ));
        }
        let header = HicHeader::new(
            head.genome.as_deref().unwrap_or("unknown"),
            chroms
                .iter()
                .map(|(name, len)| Chrom {
                    name: name.clone(),
                    len: *len,
                })
                .collect(),
            resolutions.clone(),
        );
        let mut sink = HicSink::create(ctx, &self.dst, header, self.level, self.sidecar)?;
        let layouts = sink.layouts().clone();
        let columns = Columns::of(self.dialect, &head)?;
        let mut ids = HashMap::with_capacity(chroms.len());
        for (i, (name, _)) in chroms.iter().enumerate() {
            ids.insert(name.as_bytes().to_vec(), i as u32);
        }
        let rows = Arc::new(Rows {
            columns,
            ids,
            lens: chroms
                .iter()
                .map(|(_, l)| (*l).min(u32::MAX as u64) as u32)
                .collect(),
            layouts: layouts.clone(),
            finest,
        });
        let scored = rows.columns.score.is_some();

        // Read, parse and bin: contacts into the finest sorter.
        let tmp = ctx.tmp_parent(&self.dst);
        let coarse: Vec<(u32, u32)> = resolutions[..resolutions.len() - 1]
            .iter()
            .map(|&r| (r, r / finest))
            .collect();
        let buffer = ctx.memory.buffer;
        let fine_limit = if coarse.is_empty() {
            buffer
        } else {
            buffer / 2
        };
        let mut fine_sorter: Sorter<'_, Contact> = Sorter::new(&ctx.pool, fine_limit, &tmp);
        let mut batches = LineBatches::with_prefix(
            text,
            text_batch(ctx.budget),
            head.lines + 1,
            std::mem::take(&mut head.first_row),
        );
        let mut exec = ctx.pool.ordered::<Parsed>(s.limit);
        let mut taken = Taken::default();
        while let Some(batch) = batches.next_batch()? {
            while exec.full() {
                taken.take(exec.next(true)?.expect("pending"), &mut fine_sorter)?;
            }
            let rows = Arc::clone(&rows);
            exec.submit(Box::new(move || parse(&rows, &batch)))?;
            while let Some(p) = exec.next(false)? {
                taken.take(p, &mut fine_sorter)?;
            }
            ctx.progress("read", read.load(Ordering::Relaxed), file_len);
        }
        while let Some(p) = exec.next(true)? {
            taken.take(p, &mut fine_sorter)?;
        }
        if let Some(w) = taken.warning(&chroms) {
            ctx.warn(w);
        }
        let mut spilled = fine_sorter.spilled();
        ctx.progress("sort", read.load(Ordering::Relaxed), file_len);

        // Finest pixels: kept for last, and summed into the coarser maps.
        let mut fine = fine_sorter.finish()?;
        let mut coarse_sorters: Vec<Sorter<'_, Contact>> = coarse
            .iter()
            .map(|_| Sorter::new(&ctx.pool, buffer / 2 / coarse.len().max(1) as u64, &tmp))
            .collect();
        let mut kept: Option<(TempDir, PathBuf, BufWriter<File>)> = None;
        if !coarse.is_empty() {
            let dir = TempDir::new(&tmp)?;
            let path = dir.path().join("finest");
            let file = File::create(&path).map_err(|e| with_path(e, &path))?;
            kept = Some((dir, path, BufWriter::with_capacity(WRITE_BUF, file)));
        }
        let c = Self::chunk_sizing(ctx.budget, coarse.len());
        let mut chunks = ctx.pool.ordered::<Vec<Vec<Contact>>>(c.limit);
        let shared = Arc::new((coarse.clone(), layouts));
        let mut chunk: Vec<Fine> = Vec::new();
        let mut chunk_no = 0u64;
        let mut buf = Vec::with_capacity(24);
        let mut current: Option<Fine> = None;
        let mut pixels = 0u64;
        loop {
            let next = fine.next_record()?;
            if let (Some(f), Some(c)) = (&mut current, &next)
                && f.cell() == c.cell()
            {
                f.value += c.value;
                continue;
            }
            if let Some(f) = current.take() {
                pixels += 1;
                match &mut kept {
                    None => sink.push(finest, f.pixel())?,
                    Some((_, _, out)) => {
                        buf.clear();
                        f.pixel().encode(&mut buf);
                        out.write_all(&buf)?;
                        spilled += buf.len() as u64;
                        chunk.push(f);
                    }
                }
                if chunk.len() == CHUNK || (next.is_none() && !chunk.is_empty()) {
                    while chunks.full() {
                        feed(&mut coarse_sorters, chunks.next(true)?.expect("pending"))?;
                    }
                    let job = std::mem::take(&mut chunk);
                    let shared = Arc::clone(&shared);
                    let tie = if scored { chunk_no } else { 0 };
                    chunk_no += 1;
                    chunks.submit(Box::new(move || coarsen(&job, &shared.0, &shared.1, tie)))?;
                    while let Some(done) = chunks.next(false)? {
                        feed(&mut coarse_sorters, done)?;
                    }
                }
                if pixels.is_multiple_of(1 << 16) {
                    ctx.progress("bin", pixels, 0);
                }
            }
            match next {
                Some(c) => current = Some(Fine::of(&c)),
                None => break,
            }
        }
        drop(fine);
        while let Some(done) = chunks.next(true)? {
            feed(&mut coarse_sorters, done)?;
        }

        // Write coarsest first, the finest last.
        ctx.progress("write", 0, resolutions.len() as u64);
        for (i, (sorter, &(res, _))) in coarse_sorters.into_iter().zip(&coarse).enumerate() {
            spilled += sorter.spilled();
            let mut merged = sorter.finish()?;
            let mut current: Option<Fine> = None;
            while let Some(c) = merged.next_record()? {
                match &mut current {
                    Some(f) if f.cell() == c.cell() => f.value += c.value,
                    _ => {
                        if let Some(f) = current.replace(Fine::of(&c)) {
                            sink.push(res, f.pixel())?;
                        }
                    }
                }
            }
            if let Some(f) = current {
                sink.push(res, f.pixel())?;
            }
            ctx.progress("write", i as u64 + 1, resolutions.len() as u64);
        }
        if let Some((dir, path, out)) = kept {
            out.into_inner().map_err(|e| e.into_error())?;
            let file = File::open(&path).map_err(|e| with_path(e, &path))?;
            let mut input = BufReader::with_capacity(WRITE_BUF, file);
            while let Some(p) = Pixel::decode(&mut input)? {
                sink.push(finest, p)?;
            }
            drop(dir);
        }
        let written = sink.finish()?;
        Ok(Report {
            files: written,
            records: taken.rows,
            bytes_read: read.load(Ordering::Relaxed),
            spilled,
            elapsed: Default::default(),
            warnings: Vec::new(),
        })
    }
}

/// What the header says, and how many lines it takes.
#[derive(Debug, Default)]
struct Head {
    chroms: Vec<(String, u64)>,
    genome: Option<String>,
    columns: Option<Vec<String>>,
    lines: u64,
    /// Short text: the first row has nine fields, the ninth a score.
    scored: bool,
    /// The first data row, read while telling the short dialect.
    first_row: Vec<u8>,
}

/// Reads the `#` and blank lines that open the text, and looks at the first
/// row.
fn read_head(text: &mut dyn BufRead, dialect: PairsDialect, limit: usize) -> Result<Head> {
    let mut head = Head::default();
    let mut line = Vec::new();
    loop {
        line.clear();
        if text.take(limit as u64 + 1).read_until(b'\n', &mut line)? == 0 {
            break;
        }
        if line.len() > limit {
            return Err(Error::invalid(
                "pairs",
                format!(
                    "header or first data line {} exceeds the {limit}-byte probe limit; \
                     shorten the line or increase memory (the probe is capped at 8 MiB)",
                    head.lines + 1
                ),
            ));
        }
        let header = line.first() == Some(&b'#');
        if !header && !line.iter().all(u8::is_ascii_whitespace) {
            if dialect == PairsDialect::Short {
                head.scored = words(&line).count() == 9;
            }
            head.first_row = line;
            break;
        }
        head.lines += 1;
        if !header || dialect != PairsDialect::Pairs {
            continue;
        }
        let l = std::str::from_utf8(&line)
            .map_err(|_| bad(head.lines, "header line is not UTF-8".into()))?
            .trim_end_matches(['\n', '\r']);
        if let Some(rest) = l.strip_prefix("#chromsize:") {
            let mut w = rest.split_ascii_whitespace();
            let (Some(name), Some(len), None) = (w.next(), w.next(), w.next()) else {
                return Err(bad(
                    head.lines,
                    format!("{l:?}: expected a name and a length"),
                ));
            };
            let len = len
                .parse::<u64>()
                .map_err(|_| bad(head.lines, format!("{l:?}: bad length")))?;
            head.chroms.push((name.to_owned(), len));
        } else if let Some(rest) = l.strip_prefix("#columns:") {
            head.columns = Some(rest.split_ascii_whitespace().map(str::to_owned).collect());
        } else if let Some(rest) = l.strip_prefix("#genome_assembly:") {
            let g = rest.trim();
            if !g.is_empty() {
                head.genome = Some(g.to_owned());
            }
        }
    }
    Ok(head)
}

fn bad(line: u64, detail: String) -> Error {
    Error::Format { line, detail }
}

/// Where a row's fields are.
#[derive(Debug, Clone, Copy)]
struct Columns {
    dialect: PairsDialect,
    chr1: usize,
    pos1: usize,
    chr2: usize,
    pos2: usize,
    score: Option<usize>,
    /// Fields a row has at least (pairs) or exactly (short).
    fields: usize,
}

impl Columns {
    fn of(dialect: PairsDialect, head: &Head) -> Result<Columns> {
        if dialect == PairsDialect::HicPro {
            return Ok(Columns {
                dialect,
                chr1: 1,
                pos1: 2,
                chr2: 4,
                pos2: 5,
                score: None,
                fields: 12,
            });
        }
        if dialect == PairsDialect::Short {
            return Ok(Columns {
                dialect,
                chr1: 1,
                pos1: 2,
                chr2: 5,
                pos2: 6,
                score: head.scored.then_some(8),
                fields: if head.scored { 9 } else { 8 },
            });
        }
        let names: Vec<&str> = match &head.columns {
            Some(c) => c.iter().map(String::as_str).collect(),
            None => PAIRS_COLUMNS.to_vec(),
        };
        let find = |alts: &[&str]| {
            names.iter().position(|n| alts.contains(n)).ok_or_else(|| {
                Error::invalid(
                    "pairs",
                    format!("the #columns line names no {} column", alts[0]),
                )
            })
        };
        let (chr1, pos1) = (find(&["chr1"])?, find(&["pos1", "position1"])?);
        let (chr2, pos2) = (find(&["chr2"])?, find(&["pos2", "position2"])?);
        Ok(Columns {
            dialect,
            chr1,
            pos1,
            chr2,
            pos2,
            score: None,
            fields: chr1.max(pos1).max(chr2).max(pos2) + 1,
        })
    }
}

/// What the parse tasks share.
struct Rows {
    columns: Columns,
    ids: HashMap<Vec<u8>, u32>,
    /// Chromosome lengths, saturated to `u32` (positions clamp to them).
    lens: Vec<u32>,
    layouts: MatrixLayouts,
    finest: u32,
}

/// What one batch parses into.
struct Parsed {
    contacts: Vec<Contact>,
    rows: u64,
    /// Rows skipped, and per unknown chromosome the rows naming it.
    skipped: u64,
    unknown: Vec<(Vec<u8>, u64)>,
}

/// Whitespace-separated words.
fn words(line: &[u8]) -> impl Iterator<Item = &[u8]> {
    line.split(|c| c.is_ascii_whitespace())
        .filter(|w| !w.is_empty())
}

/// A non-negative integer of ASCII digits.
fn digits(field: &[u8]) -> Option<u64> {
    if field.is_empty() {
        return None;
    }
    if field.len() > 19 {
        return field.iter().try_fold(0u64, |v, &c| {
            if !c.is_ascii_digit() {
                return None;
            }
            v.checked_mul(10)?.checked_add(u64::from(c - b'0'))
        });
    }
    // Every integer of at most 19 digits fits in u64, so ordinary coordinates need
    // no overflow check per digit. Longer spellings use checked arithmetic.
    let mut v = 0u64;
    for &c in field {
        if !c.is_ascii_digit() {
            return None;
        }
        v = v * 10 + (c - b'0') as u64;
    }
    Some(v)
}

/// A Juicer short row: `str1 chr1 pos1 frag1 str2 chr2 pos2 frag2
/// [score]`, the strands and fragments integers, the chromosomes not
/// BED-style strands.
pub(crate) fn is_short_line(line: &[u8]) -> bool {
    let mut f: [&[u8]; 10] = [&[]; 10];
    let mut n = 0;
    for w in words(line) {
        if n == 10 {
            return false;
        }
        f[n] = w;
        n += 1;
    }
    let int = |w: &[u8]| {
        let w = w.strip_prefix(b"-").unwrap_or(w);
        digits(w).is_some()
    };
    let chrom = |w: &[u8]| !matches!(w, b"+" | b"-" | b".");
    (n == 8 || n == 9)
        && [0, 2, 3, 4, 6, 7].iter().all(|&i| int(f[i]))
        && chrom(f[1])
        && chrom(f[5])
        && (n == 8
            || std::str::from_utf8(f[8])
                .ok()
                .is_some_and(|s| s.parse::<f64>().is_ok()))
}

/// HiC-Pro allows an optional allele column and a trailing tab.
pub(crate) fn is_hicpro_line(line: &[u8]) -> bool {
    let mut fields: [&[u8]; 13] = [&[]; 13];
    let mut n = 0;
    for field in words(line) {
        if n == fields.len() {
            return false;
        }
        fields[n] = field;
        n += 1;
    }
    hicpro_positions(&fields[..n]).is_some()
}

fn hicpro_positions(fields: &[&[u8]]) -> Option<(u64, u64)> {
    if !matches!(fields.len(), 12 | 13)
        || !matches!(fields[3], b"+" | b"-")
        || !matches!(fields[6], b"+" | b"-")
    {
        return None;
    }
    Some((digits(fields[2])?, digits(fields[5])?))
}

/// The id of chromosome `name`, remembering the last one asked for.
fn chrom_id<'a>(rows: &Rows, name: &'a [u8], last: &mut (&'a [u8], Option<u32>)) -> Option<u32> {
    if last.0 != name {
        *last = (name, rows.ids.get(name).copied());
    }
    last.1
}

fn parse(rows: &Rows, batch: &LineBatch) -> Result<Parsed> {
    let cols = rows.columns;
    let mut out = Parsed {
        contacts: Vec::with_capacity(batch.text.len() / 24),
        rows: 0,
        skipped: 0,
        unknown: Vec::new(),
    };
    let mut unknown: HashMap<&[u8], u64> = HashMap::new();
    let field_limit = cols.fields
        + if cols.dialect == PairsDialect::HicPro {
            2
        } else {
            1
        };
    let mut fields: Vec<&[u8]> = Vec::with_capacity(field_limit);
    let mut last = [(&b""[..], None), (&b""[..], None)];
    let mut layout: Option<((u32, u32), MatrixLayout)> = None;
    for (no, line) in lines(batch) {
        let line = line.strip_suffix(b"\r").unwrap_or(line);
        if line.first().is_none_or(|&c| c == b'#') || line.iter().all(u8::is_ascii_whitespace) {
            continue;
        }
        fields.clear();
        let tabs = cols.dialect == PairsDialect::Pairs && memchr::memchr(b'\t', line).is_some();
        if tabs {
            fields.extend(line.split(|&c| c == b'\t').take(cols.fields));
        } else {
            fields.extend(words(line).take(field_limit));
        }
        let hicpro = if cols.dialect == PairsDialect::HicPro {
            hicpro_positions(&fields)
        } else {
            None
        };
        let ok = match cols.dialect {
            PairsDialect::Pairs => fields.len() >= cols.fields,
            PairsDialect::Short => fields.len() == cols.fields,
            PairsDialect::HicPro => hicpro.is_some(),
        };
        if !ok {
            let n = match cols.dialect {
                PairsDialect::Short | PairsDialect::HicPro => words(line).count(),
                PairsDialect::Pairs => fields.len(),
            };
            return Err(bad(
                no,
                match cols.dialect {
                    PairsDialect::Pairs => {
                        format!("{n} fields; the columns name {} or more", cols.fields)
                    }
                    PairsDialect::Short => format!(
                        "{n} fields; Juicer short rows have {} (str1 chr1 pos1 frag1 str2 chr2 \
                         pos2 frag2{}) as the first row has; the medium and long formats are \
                         not read",
                        cols.fields,
                        if cols.fields == 9 { " score" } else { "" }
                    ),
                    PairsDialect::HicPro => format!(
                        "{n} fields; HiC-Pro rows require 12 or 13 fields, integer positions \
                         in columns 3 and 6 and +/- strands in columns 4 and 7"
                    ),
                },
            ));
        }
        let position = |i: usize| {
            digits(fields[i]).ok_or_else(|| {
                bad(
                    no,
                    format!(
                        "position {:?} is not a coordinate",
                        String::from_utf8_lossy(fields[i])
                    ),
                )
            })
        };
        let (pa, pb) = match hicpro {
            Some(positions) => positions,
            None => (position(cols.pos1)?, position(cols.pos2)?),
        };
        let count = match cols.score {
            Some(i) => crate::wig::value(fields[i], no)? as f64,
            None => 1.0,
        };
        let (na, nb) = (fields[cols.chr1], fields[cols.chr2]);
        let [la, lb] = &mut last;
        let (a, b) = (chrom_id(rows, na, la), chrom_id(rows, nb, lb));
        let (Some(a), Some(b)) = (a, b) else {
            out.skipped += 1;
            if a.is_none() {
                *unknown.entry(na).or_default() += 1;
            }
            if b.is_none() && (a.is_some() || na != nb) {
                *unknown.entry(nb).or_default() += 1;
            }
            continue;
        };
        let (c1, p1, c2, p2) = if a <= b {
            (a, pa, b, pb)
        } else {
            (b, pb, a, pa)
        };
        let clamp = |c: u32, p: u64| {
            // SAFETY: c comes only from rows.ids; run builds ids and lens
            // from the same chromosome list. Rows is immutable during parsing.
            let len = unsafe { *rows.lens.get_unchecked(c as usize) };
            p.min(u64::from(len)) as u32
        };
        let (mut x, mut y) = (clamp(c1, p1) / rows.finest, clamp(c2, p2) / rows.finest);
        if c1 == c2 && x > y {
            std::mem::swap(&mut x, &mut y);
        }
        let l = match layout {
            Some((pair, l)) if pair == (c1, c2) => l,
            _ => {
                let l = rows.layouts.layout(c1, c2, ContactUnit::Bp, rows.finest)?;
                layout = Some(((c1, c2), l));
                l
            }
        };
        let (block, x, y) = l.place(x, y);
        out.contacts.push(Contact {
            c1,
            c2,
            block,
            y,
            x,
            tie: if cols.score.is_some() { no } else { 0 },
            value: count,
        });
        out.rows += 1;
    }
    if cols.score.is_none() {
        // Counts of one cell fold here already: integer sums are exact.
        out.contacts.sort_unstable_by_key(Contact::key);
        out.contacts.dedup_by(|b, a| {
            let same = a.key() == b.key();
            if same {
                a.value += b.value;
            }
            same
        });
    }
    out.unknown = unknown.into_iter().map(|(k, v)| (k.to_vec(), v)).collect();
    Ok(out)
}

/// What the serial stage gathered from the batches.
#[derive(Default)]
struct Taken {
    rows: u64,
    skipped: u64,
    unknown: BTreeMap<Vec<u8>, u64>,
}

impl Taken {
    fn take(&mut self, p: Parsed, sorter: &mut Sorter<'_, Contact>) -> Result<()> {
        self.rows += p.rows;
        self.skipped += p.skipped;
        for (name, n) in p.unknown {
            *self.unknown.entry(name).or_default() += n;
        }
        sorter.extend(p.contacts)
    }

    /// The warning for the rows skipped, if any: per chromosome, most rows
    /// first.
    fn warning(&self, chroms: &[(String, u64)]) -> Option<String> {
        if self.skipped == 0 {
            return None;
        }
        let mut by: Vec<(&Vec<u8>, u64)> = self.unknown.iter().map(|(k, &v)| (k, v)).collect();
        by.sort_by(|a, b| b.1.cmp(&a.1).then_with(|| a.0.cmp(b.0)));
        const SHOWN: usize = 10;
        let mut list: Vec<String> = by
            .iter()
            .take(SHOWN)
            .map(|(name, n)| format!("{} ({n})", String::from_utf8_lossy(name)))
            .collect();
        if by.len() > SHOWN {
            list.push(format!("{} more", by.len() - SHOWN));
        }
        Some(format!(
            "skipped {} of {} rows on chromosomes not among the {} listed (rows per \
             chromosome): {}",
            self.skipped,
            self.skipped + self.rows,
            chroms.len(),
            list.join(", ")
        ))
    }
}

/// Hands each coarser resolution its contacts of a chunk.
fn feed(sorters: &mut [Sorter<'_, Contact>], done: Vec<Vec<Contact>>) -> Result<()> {
    for (sorter, contacts) in sorters.iter_mut().zip(done) {
        sorter.extend(contacts)?;
    }
    Ok(())
}

/// Sums a chunk of finest pixels into the cells of each coarser resolution
/// (`(resolution, multiple of the finest)`), folded in chunk order.
fn coarsen(
    chunk: &[Fine],
    coarse: &[(u32, u32)],
    layouts: &MatrixLayouts,
    tie: u64,
) -> Result<Vec<Vec<Contact>>> {
    coarse
        .iter()
        .map(|&(res, factor)| {
            let mut out = Vec::with_capacity(chunk.len());
            let mut layout: Option<((u32, u32), MatrixLayout)> = None;
            for p in chunk {
                let l = match layout {
                    Some((pair, l)) if pair == (p.c1, p.c2) => l,
                    _ => {
                        let l = layouts.layout(p.c1, p.c2, ContactUnit::Bp, res)?;
                        layout = Some(((p.c1, p.c2), l));
                        l
                    }
                };
                let (block, x, y) = l.place(p.x / factor, p.y / factor);
                out.push(Contact {
                    c1: p.c1,
                    c2: p.c2,
                    block,
                    y,
                    x,
                    tie,
                    value: p.value,
                });
            }
            // Stable: the pixels of one cell are summed in chunk order.
            out.sort_by_key(Contact::key);
            out.dedup_by(|b, a| {
                let same = a.key() == b.key();
                if same {
                    a.value += b.value;
                }
                same
            });
            Ok(out)
        })
        .collect()
}

/// A finest pixel and its sum so far.
#[derive(Debug, Clone, Copy)]
struct Fine {
    c1: u32,
    c2: u32,
    block: u32,
    x: u32,
    y: u32,
    value: f64,
}

impl Fine {
    fn of(c: &Contact) -> Fine {
        Fine {
            c1: c.c1,
            c2: c.c2,
            block: c.block,
            x: c.x,
            y: c.y,
            value: c.value,
        }
    }

    fn cell(&self) -> (u32, u32, u32, u32, u32) {
        (self.c1, self.c2, self.block, self.y, self.x)
    }

    fn pixel(&self) -> Pixel {
        Pixel {
            c1: self.c1,
            c2: self.c2,
            block: self.block,
            x: self.x,
            y: self.y,
            value: self.value as f32,
        }
    }
}

/// A contact (or a sum of contacts) placed in its cell at one resolution.
/// `tie` keeps fractional scores apart until the last merge (the row, or
/// the chunk); counts have none, and fold wherever they meet.
#[derive(Debug, Clone, Copy)]
struct Contact {
    c1: u32,
    c2: u32,
    block: u32,
    y: u32,
    x: u32,
    tie: u64,
    value: f64,
}

impl Contact {
    fn key(&self) -> (u32, u32, u32, u32, u32, u64) {
        (self.c1, self.c2, self.block, self.y, self.x, self.tie)
    }

    fn cell(&self) -> (u32, u32, u32, u32, u32) {
        (self.c1, self.c2, self.block, self.y, self.x)
    }
}

fn put_varint(out: &mut Vec<u8>, mut v: u64) {
    while v >= 0x80 {
        out.push(v as u8 | 0x80);
        v >>= 7;
    }
    out.push(v as u8);
}

fn get_varint(input: &mut impl Read, first: u8) -> Result<u64> {
    let (mut v, mut shift, mut b) = (0u64, 0u32, first);
    loop {
        v |= ((b & 0x7f) as u64) << shift;
        if b & 0x80 == 0 {
            return Ok(v);
        }
        shift += 7;
        if shift > 63 {
            return Err(Error::invalid("spill", "bad varint in a run file"));
        }
        let mut one = [0u8];
        if !read_or_end(input, &mut one)? {
            return Err(Error::invalid("spill", "truncated run file"));
        }
        b = one[0];
    }
}

fn next_varint(input: &mut impl Read) -> Result<u64> {
    let mut one = [0u8];
    if !read_or_end(input, &mut one)? {
        return Err(Error::invalid("spill", "truncated run file"));
    }
    get_varint(input, one[0])
}

/// Largest count stored as an integer in a run file.
const INTEGER_VALUES: f64 = (1u64 << 52) as f64;

impl Spill for Contact {
    type Key = (u32, u32, u32, u32, u32, u64);
    const ABSORBS: bool = true;

    fn key(&self) -> Self::Key {
        Contact::key(self)
    }
    fn mem_size(&self) -> usize {
        0
    }
    /// Varints: counts and small numbers take a byte or two.
    fn encode(&self, out: &mut Vec<u8>) {
        for v in [self.c1, self.c2, self.block, self.y, self.x] {
            put_varint(out, v as u64);
        }
        put_varint(out, self.tie);
        let v = self.value;
        if (0.0..INTEGER_VALUES).contains(&v) && (v as u64 as f64).to_bits() == v.to_bits() {
            put_varint(out, (v as u64) << 1);
        } else {
            out.push(1);
            out.extend_from_slice(&v.to_bits().to_le_bytes());
        }
    }
    fn decode(input: &mut impl Read) -> Result<Option<Contact>> {
        let mut one = [0u8];
        if !read_or_end(input, &mut one)? {
            return Ok(None);
        }
        let c1 = get_varint(input, one[0])? as u32;
        let c2 = next_varint(input)? as u32;
        let block = next_varint(input)? as u32;
        let y = next_varint(input)? as u32;
        let x = next_varint(input)? as u32;
        let tie = next_varint(input)?;
        let tag = next_varint(input)?;
        let value = if tag & 1 == 0 {
            (tag >> 1) as f64
        } else {
            let mut b = [0u8; 8];
            if !read_or_end(input, &mut b)? {
                return Err(Error::invalid("spill", "truncated run file"));
            }
            f64::from_bits(u64::from_le_bytes(b))
        };
        Ok(Some(Contact {
            c1,
            c2,
            block,
            y,
            x,
            tie,
            value,
        }))
    }
    fn absorb(&mut self, next: Contact) {
        self.value += next.value;
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Cursor;

    #[test]
    fn coordinate_digits_cover_u64_without_wrapping() {
        assert_eq!(digits(b"0"), Some(0));
        assert_eq!(
            digits(b"9999999999999999999"),
            Some(9_999_999_999_999_999_999)
        );
        assert_eq!(digits(b"18446744073709551615"), Some(u64::MAX));
        assert_eq!(digits(b"00018446744073709551615"), Some(u64::MAX));
        for text in [
            &b""[..],
            b"18446744073709551616",
            b"99999999999999999999",
            b"-1",
            b"1a",
        ] {
            assert_eq!(digits(text), None, "{text:?}");
        }
    }

    #[test]
    fn hicpro_field_validation_reuses_coordinates_and_rejects_extra_fields() {
        for ending in ["", "\t", "\t1", "\t1\t", "\t1\t\r\n"] {
            let line =
                format!("r\tchr1\t18446744073709551615\t+\tchr2\t0\t-\t0\tf1\tf2\t0\t0{ending}");
            assert!(is_hicpro_line(line.as_bytes()));
            let fields: Vec<_> = words(line.as_bytes()).collect();
            assert_eq!(hicpro_positions(&fields), Some((u64::MAX, 0)));
        }
        for line in [
            "r chr1 1 + chr2 2 - 0 f1 f2 0",
            "r chr1 1 + chr2 2 - 0 f1 f2 0 0 1 extra",
            "r chr1 18446744073709551616 + chr2 2 - 0 f1 f2 0 0",
            "r chr1 -1 + chr2 2 - 0 f1 f2 0 0",
            "r chr1 1 . chr2 2 - 0 f1 f2 0 0",
        ] {
            assert!(!is_hicpro_line(line.as_bytes()), "{line}");
        }
    }

    #[test]
    fn header_and_first_row_are_independent_of_reader_boundaries() {
        let source = b" \t\r\n#chromsize: chr1 1000\r\n#columns: readID chr1 pos1 chr2 pos2\nr\tchr1\t2\tchr1\t5\n";
        for capacity in [1, 2, 7, 64, 4096] {
            let mut input = BufReader::with_capacity(capacity, Cursor::new(source));
            let mut head = read_head(&mut input, PairsDialect::Pairs, 4096).unwrap();
            assert_eq!(head.chroms, [("chr1".into(), 1000)]);
            assert_eq!(head.lines, 3);
            let mut restored = LineBatches::with_prefix(
                input,
                16,
                head.lines + 1,
                std::mem::take(&mut head.first_row),
            );
            let mut rest = Vec::new();
            while let Some(batch) = restored.next_batch().unwrap() {
                rest.extend_from_slice(&batch.text);
            }
            assert_eq!(rest, b"r\tchr1\t2\tchr1\t5\n");
        }
        for ending in ["", "\n", "\r\n"] {
            let source = format!(
                " \t\r\n#ignored\n0 {}chr1 2 0 16 chr1 5 1 2.5{ending}",
                " ".repeat(8192)
            );
            for capacity in [1, 7, 1024] {
                let mut input = BufReader::with_capacity(capacity, Cursor::new(&source));
                let mut head = read_head(&mut input, PairsDialect::Short, 16384).unwrap();
                assert!(head.scored);
                assert_eq!(head.lines, 2);
                let mut restored = LineBatches::with_prefix(
                    input,
                    16,
                    head.lines + 1,
                    std::mem::take(&mut head.first_row),
                );
                let mut text = Vec::new();
                while let Some(batch) = restored.next_batch().unwrap() {
                    text.extend_from_slice(&batch.text);
                }
                assert_eq!(text, source.splitn(3, '\n').nth(2).unwrap().as_bytes());
            }
        }
    }

    #[test]
    fn first_line_probe_has_a_bounded_extent() {
        let source = vec![b'#'; 1024];
        let mut input = Cursor::new(&source);
        assert!(read_head(&mut input, PairsDialect::Pairs, 63).is_err());
        assert_eq!(input.position(), 64);
    }

    #[test]
    fn contacts_round_trip_through_run_files() {
        for value in [0.0, 1.0, 7.0, 0.5, -3.0, -0.0, 1e300, (1u64 << 53) as f64] {
            let c = Contact {
                c1: 3,
                c2: 70_000,
                block: 1 << 30,
                y: 5,
                x: 0,
                tie: u64::MAX,
                value,
            };
            let mut buf = Vec::new();
            c.encode(&mut buf);
            let back = Contact::decode(&mut buf.as_slice()).unwrap().unwrap();
            assert_eq!(back.key(), c.key());
            assert_eq!(back.value.to_bits(), value.to_bits());
        }
    }

    #[test]
    fn short_rows_are_told_from_bed() {
        assert!(is_short_line(b"0 chr1 100 0 16 chr2 500 1"));
        assert!(is_short_line(b"0\t1\t100\t0\t0\t1\t500\t1\t2.5"));
        assert!(!is_short_line(b"chr1\t100\t200\tname\t0\t+\t100\t200"));
        assert!(!is_short_line(b"0 chr1 100 0 16 chr2 500"));
    }
}
