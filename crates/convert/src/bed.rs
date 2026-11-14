//! BED → bigBed (+ sidecar).

use std::collections::HashMap;
use std::io::{BufRead, BufReader, BufWriter, Read};
use std::path::{Path, PathBuf};
use std::sync::Arc;
use std::sync::atomic::Ordering;

use omicsweave_bbi::{BedSchema, BigBedWriter, Chrom, WriteOptions};
use omicsweave_core::exec::{Exec, Ordered};
use omicsweave_core::{Error, FileSource, Result};

use crate::budget::{Budget, MemoryPlan, WRITE_BUF, sizing, text_batch};
use crate::conversion::{Context, Report};
use crate::format::Format;
use crate::pipeline::input::{LineBatch, LineBatches, READ_CHUNK, open_text, with_path};
use crate::pipeline::output::Outputs;
use crate::pipeline::sort::{Sorter, Spill, read_or_end};

/// Where chromosome lengths come from.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum ChromSizes {
    /// A two-column `name<TAB>length` file (`chrom.sizes`, or a `.fai`,
    /// whose first two columns are the same), or a `.2bit` file.
    File(PathBuf),
    /// The largest end of each chromosome's intervals. The file is then
    /// self-consistent but its chromosome lengths are not the genome's.
    FromData,
}

/// BED (plain, gzip or BGZF) to bigBed, with its sidecar.
///
/// Records need not be sorted: they are sorted by chromosome (in order of
/// first appearance) and start, ties kept in input order, in memory or out
/// of core within the budget. Fields may be tab- or whitespace-separated;
/// fields after the third are stored tab-joined. `track`, `browser`, `#`
/// and blank lines are skipped. Every record must have the same number of
/// fields.
#[derive(Debug, Clone, PartialEq)]
pub struct BedToBigBed {
    pub src: PathBuf,
    pub dst: PathBuf,
    pub chrom_sizes: ChromSizes,
    pub options: WriteOptions,
    /// Also write `<dst>.sidecar`.
    pub sidecar: bool,
}

impl BedToBigBed {
    pub fn new(
        src: impl Into<PathBuf>,
        dst: impl Into<PathBuf>,
        chrom_sizes: ChromSizes,
    ) -> BedToBigBed {
        BedToBigBed {
            src: src.into(),
            dst: dst.into(),
            chrom_sizes,
            options: WriteOptions::default(),
            sidecar: true,
        }
    }

    fn parse_sizing(budget: &Budget) -> crate::budget::Sizing {
        // Text plus parsed records, about twice the text.
        sizing(
            budget.threads,
            3 * text_batch(budget) as u64,
            budget.memory / 4,
        )
    }

    fn block_sizing(budget: &Budget) -> crate::budget::Sizing {
        sizing(budget.threads, 256 << 10, budget.memory / 16)
    }

    pub(crate) fn plan(&self, budget: &Budget) -> Result<MemoryPlan> {
        if let ChromSizes::File(p) = &self.chrom_sizes
            && !p.exists()
        {
            return Err(Error::Io {
                kind: std::io::ErrorKind::NotFound,
                detail: format!("chromosome sizes {} do not exist", p.display()),
            });
        }
        let pipeline = Self::parse_sizing(budget).per_queue
            + Self::block_sizing(budget).per_queue
            + 2 * READ_CHUNK as u64
            + WRITE_BUF as u64;
        Ok(MemoryPlan {
            fixed: 0,
            pipeline,
            buffer: budget.memory.saturating_sub(pipeline),
        })
    }

    pub(crate) fn run(&self, ctx: &Context<'_>) -> Result<Report> {
        let sizes = match &self.chrom_sizes {
            ChromSizes::File(p) => Some(read_chrom_sizes(p)?),
            ChromSizes::FromData => None,
        };
        let s = Self::parse_sizing(ctx.budget);
        let input = open_text(&self.src, &ctx.pool, s.limit)?;
        let (read, file_len) = (Arc::clone(&input.read), input.file_len);
        let mut batches = LineBatches::new(input.text, text_batch(ctx.budget), 1);
        let mut exec = ctx.pool.ordered::<Parsed>(s.limit);
        let tmp = ctx.tmp_parent(&self.dst);
        let mut gather = Gather {
            sorter: Sorter::new(&ctx.pool, ctx.memory.buffer, &tmp),
            ids: HashMap::new(),
            names: Vec::new(),
            max_end: Vec::new(),
            fields: None,
            seq: 0,
            sizes: sizes.as_ref(),
        };
        while let Some(batch) = batches.next_batch()? {
            while exec.full() {
                gather.take(exec.next(true)?.expect("pending"))?;
            }
            exec.submit(Box::new(move || parse(&batch)))?;
            while let Some(p) = exec.next(false)? {
                gather.take(p)?;
            }
            ctx.progress("read", read.load(Ordering::Relaxed), file_len);
        }
        while let Some(p) = exec.next(true)? {
            gather.take(p)?;
        }
        let records = gather.seq;
        let spilled = gather.sorter.spilled();
        let chroms: Vec<Chrom> = gather
            .names
            .iter()
            .zip(&gather.max_end)
            .map(|(name, &end)| Chrom {
                name: name.clone(),
                len: sizes.as_ref().map_or(end, |s| s[name]),
            })
            .collect();
        if chroms.is_empty() {
            return Err(Error::invalid(
                "bed",
                "no intervals; a bigBed needs at least one chromosome",
            ));
        }
        let fields = gather.fields.unwrap_or(3);
        let schema = BedSchema {
            field_count: fields,
            defined_field_count: fields.min(12),
            auto_sql: None,
        };
        ctx.progress("sort", read.load(Ordering::Relaxed), file_len);
        let mut merged = gather.sorter.finish()?;
        let mut outputs = Outputs::new();
        let file = outputs.create(&self.dst, Format::BigBed)?;
        let b = Self::block_sizing(ctx.budget);
        let mut writer = BigBedWriter::create_with(
            BufWriter::with_capacity(WRITE_BUF, file),
            chroms,
            schema,
            self.options.clone(),
            Exec::Pool {
                pool: &ctx.pool,
                limit: b.limit,
            },
        )?;
        let mut n = 0u64;
        while let Some(r) = merged.next_record()? {
            writer.add(r.chrom, r.start, r.end, &r.rest)?;
            n += 1;
            if n % (1 << 16) == 0 {
                ctx.progress("write", n, records);
            }
        }
        drop(merged);
        let (out, sidecar) = writer.finish()?;
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
            bytes_read: read.load(Ordering::Relaxed),
            spilled,
            elapsed: Default::default(),
        })
    }
}

/// Reads `name<TAB>length` lines, or the directory of a `.2bit` file.
pub fn read_chrom_sizes(path: &Path) -> Result<HashMap<String, u32>> {
    let mut input = BufReader::new(std::fs::File::open(path).map_err(|e| with_path(e, path))?);
    let head = input.fill_buf()?;
    if head.len() >= 4
        && u32::from_le_bytes(head[..4].try_into().unwrap()) == omicsweave_fasta::TWOBIT_MAGIC
    {
        let reader = omicsweave_fasta::TwoBitReader::open(FileSource::open(path)?, false)?;
        let mut out = HashMap::new();
        for name in reader.names() {
            let len = omicsweave_core::SequenceRead::chrom_len(&reader, name)?;
            out.insert(name.to_owned(), len);
        }
        return Ok(out);
    }
    let mut out = HashMap::new();
    for (i, line) in input.lines().enumerate() {
        let line = line?;
        let line = line.trim_end_matches('\r');
        if line.is_empty() || line.starts_with('#') {
            continue;
        }
        let mut cols = line.split(['\t', ' ']).filter(|c| !c.is_empty());
        let (Some(name), Some(len)) = (cols.next(), cols.next()) else {
            return Err(Error::Format {
                line: i as u64 + 1,
                detail: format!("{}: expected a name and a length", path.display()),
            });
        };
        let len: u32 = len.parse().map_err(|_| Error::Format {
            line: i as u64 + 1,
            detail: format!("{}: bad length {len:?}", path.display()),
        })?;
        if out.insert(name.to_owned(), len).is_some() {
            return Err(Error::Format {
                line: i as u64 + 1,
                detail: format!("{}: duplicate chromosome {name}", path.display()),
            });
        }
    }
    Ok(out)
}

/// One parsed interval of a batch, its chromosome local to the batch.
struct Item {
    chrom: u32,
    start: u32,
    end: u32,
    fields: u16,
    rest: (u32, u32),
    line: u64,
}

struct Parsed {
    chroms: Vec<String>,
    items: Vec<Item>,
    rest: Vec<u8>,
}

fn parse_u32(field: &[u8], line: u64, what: &str) -> Result<u32> {
    std::str::from_utf8(field)
        .ok()
        .and_then(|s| s.parse::<u32>().ok())
        .ok_or_else(|| Error::Format {
            line,
            detail: format!(
                "{what} {:?} is not a coordinate",
                String::from_utf8_lossy(field)
            ),
        })
}

fn parse(batch: &LineBatch) -> Result<Parsed> {
    let mut out = Parsed {
        chroms: Vec::new(),
        items: Vec::new(),
        rest: Vec::with_capacity(batch.text.len() / 2),
    };
    let mut local: HashMap<&[u8], u32> = HashMap::new();
    let mut fields: Vec<&[u8]> = Vec::with_capacity(16);
    for (no, line) in crate::pipeline::input::lines(batch) {
        let line = line.strip_suffix(b"\r").unwrap_or(line);
        if line.is_empty()
            || line[0] == b'#'
            || line.starts_with(b"track")
            || line.starts_with(b"browser")
            || line.iter().all(|c| c.is_ascii_whitespace())
        {
            continue;
        }
        fields.clear();
        if memchr::memchr(b'\t', line).is_some() {
            fields.extend(line.split(|&c| c == b'\t'));
        } else {
            fields.extend(
                line.split(|c| c.is_ascii_whitespace())
                    .filter(|f| !f.is_empty()),
            );
        }
        if fields.len() < 3 {
            return Err(Error::Format {
                line: no,
                detail: format!("{} fields; BED needs at least 3", fields.len()),
            });
        }
        let start = parse_u32(fields[1], no, "start")?;
        let end = parse_u32(fields[2], no, "end")?;
        if start > end {
            return Err(Error::Format {
                line: no,
                detail: format!("start {start} after end {end}"),
            });
        }
        let next = local.len() as u32;
        let chrom = *local.entry(fields[0]).or_insert_with(|| {
            out.chroms
                .push(String::from_utf8_lossy(fields[0]).into_owned());
            next
        });
        let at = out.rest.len();
        for (i, f) in fields[3..].iter().enumerate() {
            if i > 0 {
                out.rest.push(b'\t');
            }
            out.rest.extend_from_slice(f);
        }
        out.items.push(Item {
            chrom,
            start,
            end,
            fields: fields.len().min(u16::MAX as usize) as u16,
            rest: (at as u32, (out.rest.len() - at) as u32),
            line: no,
        });
    }
    Ok(out)
}

/// A record in the sort buffer.
struct BedRec {
    chrom: u32,
    start: u32,
    end: u32,
    seq: u64,
    rest: Vec<u8>,
}

impl Spill for BedRec {
    type Key = (u32, u32, u64);
    fn key(&self) -> Self::Key {
        (self.chrom, self.start, self.seq)
    }
    fn mem_size(&self) -> usize {
        self.rest.capacity()
    }
    fn encode(&self, out: &mut Vec<u8>) {
        out.extend_from_slice(&self.chrom.to_le_bytes());
        out.extend_from_slice(&self.start.to_le_bytes());
        out.extend_from_slice(&self.end.to_le_bytes());
        out.extend_from_slice(&self.seq.to_le_bytes());
        out.extend_from_slice(&(self.rest.len() as u32).to_le_bytes());
        out.extend_from_slice(&self.rest);
    }
    fn decode(input: &mut impl Read) -> Result<Option<BedRec>> {
        let mut h = [0u8; 24];
        if !read_or_end(input, &mut h)? {
            return Ok(None);
        }
        let u = |i: usize| u32::from_le_bytes(h[i..i + 4].try_into().unwrap());
        let n = u(20) as usize;
        let mut rest = vec![0u8; n];
        input.read_exact(&mut rest)?;
        Ok(Some(BedRec {
            chrom: u(0),
            start: u(4),
            end: u(8),
            seq: u64::from_le_bytes(h[12..20].try_into().unwrap()),
            rest,
        }))
    }
}

struct Gather<'p, 's> {
    sorter: Sorter<'p, BedRec>,
    ids: HashMap<String, u32>,
    names: Vec<String>,
    max_end: Vec<u32>,
    fields: Option<u16>,
    seq: u64,
    sizes: Option<&'s HashMap<String, u32>>,
}

impl Gather<'_, '_> {
    fn take(&mut self, p: Parsed) -> Result<()> {
        let mut map = Vec::with_capacity(p.chroms.len());
        for name in p.chroms {
            let id = match self.ids.get(&name) {
                Some(&id) => id,
                None => {
                    if let Some(s) = self.sizes
                        && !s.contains_key(&name)
                    {
                        return Err(Error::invalid(
                            "bed",
                            format!("chromosome {name:?} is not in the chromosome sizes"),
                        ));
                    }
                    let id = self.names.len() as u32;
                    self.ids.insert(name.clone(), id);
                    self.names.push(name);
                    self.max_end.push(0);
                    id
                }
            };
            map.push(id);
        }
        let mut out = Vec::with_capacity(p.items.len());
        for it in p.items {
            let chrom = map[it.chrom as usize];
            match self.fields {
                None => self.fields = Some(it.fields),
                Some(f) if f != it.fields => {
                    return Err(Error::Format {
                        line: it.line,
                        detail: format!("{} fields; earlier records have {f}", it.fields),
                    });
                }
                Some(_) => {}
            }
            if let Some(s) = self.sizes {
                let len = s[&self.names[chrom as usize]];
                if it.end > len {
                    return Err(Error::Format {
                        line: it.line,
                        detail: format!(
                            "end {} past the end of {} ({len})",
                            it.end, self.names[chrom as usize]
                        ),
                    });
                }
            }
            let m = &mut self.max_end[chrom as usize];
            *m = (*m).max(it.end);
            let (at, n) = (it.rest.0 as usize, it.rest.1 as usize);
            out.push(BedRec {
                chrom,
                start: it.start,
                end: it.end,
                seq: self.seq,
                rest: p.rest[at..at + n].to_vec(),
            });
            self.seq += 1;
        }
        self.sorter.extend(out)
    }
}
