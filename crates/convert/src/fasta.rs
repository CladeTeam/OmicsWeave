//! FASTA → bgzip FASTA, FASTA → 2bit.

use std::collections::VecDeque;
use std::fs::File;
use std::io::{BufWriter, Cursor, Read, Seek, SeekFrom, Write};
use std::path::PathBuf;

use omicsweave_core::exec::Ordered;
use omicsweave_core::{Error, Result};
use omicsweave_fasta::{
    BgzfWriter, Fai, FaiBuilder, TwoBitPiece, TwoBitRecord, TwoBitWriter, fai_path, gzi_path,
    twobit_base_size,
};

use crate::budget::{Budget, MemoryPlan, WRITE_BUF, sizing, text_batch};
use crate::conversion::{Context, Report};
use crate::format::Format;
use crate::pipeline::input::{Compression, READ_CHUNK, detect, open_text, with_path};
use crate::pipeline::output::{Outputs, TempDir};

/// BGZF block, compressed and uncompressed, with slack.
const BLOCK_UNIT: u64 = 3 * omicsweave_bgzf::MAX_BLOCK_LEN as u64;

/// FASTA (plain, gzip or BGZF) to bgzip FASTA with its `.fai` and `.gzi`,
/// byte-identical to `bgzip` at the same level (and to
/// [`omicsweave_fasta::compress_fasta`]) whatever the thread count.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct FastaToBgzipFasta {
    pub src: PathBuf,
    pub dst: PathBuf,
    /// libdeflate level, 0 to 12.
    pub level: i32,
}

impl FastaToBgzipFasta {
    pub fn new(src: impl Into<PathBuf>, dst: impl Into<PathBuf>) -> FastaToBgzipFasta {
        FastaToBgzipFasta {
            src: src.into(),
            dst: dst.into(),
            level: omicsweave_bgzf::DEFAULT_LEVEL,
        }
    }

    pub(crate) fn plan(&self, budget: &Budget) -> Result<MemoryPlan> {
        omicsweave_bgzf::check_level(self.level)?;
        let s = sizing(budget.threads, BLOCK_UNIT, budget.memory / 4);
        Ok(MemoryPlan {
            fixed: 0,
            pipeline: 2 * s.per_queue + 2 * READ_CHUNK as u64 + WRITE_BUF as u64,
            buffer: 0,
        })
    }

    pub(crate) fn run(&self, ctx: &Context<'_>) -> Result<Report> {
        let s = sizing(ctx.budget.threads, BLOCK_UNIT, ctx.budget.memory / 4);
        let mut input = open_text(&self.src, &ctx.pool, s.limit)?;
        let mut outputs = Outputs::new();
        let file = outputs.create(&self.dst, Format::BgzipFasta)?;
        let mut writer = BgzfWriter::with_exec(
            BufWriter::with_capacity(WRITE_BUF, file),
            self.level,
            Box::new(ctx.pool.ordered(s.limit)),
        )?;
        let mut fai = FaiBuilder::new();
        let mut chunk = vec![0u8; READ_CHUNK];
        loop {
            let n = input.text.read(&mut chunk)?;
            if n == 0 {
                break;
            }
            fai.push(&chunk[..n])?;
            writer.write_all(&chunk[..n])?;
            ctx.progress(
                "compress",
                input.read.load(std::sync::atomic::Ordering::Relaxed),
                input.file_len,
            );
        }
        let (out, gzi) = writer.finish()?;
        out.into_inner().map_err(|e| e.into_error())?;
        let fai = fai.finish()?;
        outputs.write(&fai_path(&self.dst), Format::Fai, &fai.to_bytes())?;
        outputs.write(&gzi_path(&self.dst), Format::Gzi, &gzi.to_bytes())?;
        let records = fai.len() as u64;
        Ok(Report {
            files: outputs.commit()?,
            records,
            bytes_read: input.read.load(std::sync::atomic::Ordering::Relaxed),
            spilled: 0,
            elapsed: Default::default(),
            warnings: Vec::new(),
        })
    }
}

/// FASTA (plain, gzip or BGZF) to UCSC 2bit, as `faToTwoBit` writes it:
/// bases other than ACGT become N, lowercase is kept as soft-mask blocks.
///
/// The directory precedes the records, so names come first: from the
/// source's `.fai` when it has one (plain or BGZF sources), otherwise from
/// a first pass over the input. Each record's packed bases are held until
/// its N and mask blocks are known; a sequence too large for the budget is
/// held in a temporary file instead.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct FastaToTwoBit {
    pub src: PathBuf,
    pub dst: PathBuf,
}

/// Bases per encoding task (a multiple of 4).
fn piece_len(budget: &Budget) -> usize {
    text_batch(budget) & !3
}

impl FastaToTwoBit {
    pub fn new(src: impl Into<PathBuf>, dst: impl Into<PathBuf>) -> FastaToTwoBit {
        FastaToTwoBit {
            src: src.into(),
            dst: dst.into(),
        }
    }

    fn piece_sizing(budget: &Budget) -> crate::budget::Sizing {
        // Bases and packed bases of a piece.
        let piece = piece_len(budget) as u64;
        sizing(budget.threads, piece + piece / 4, budget.memory / 4)
    }

    pub(crate) fn plan(&self, budget: &Budget) -> Result<MemoryPlan> {
        let s = Self::piece_sizing(budget);
        let pipeline = s.per_queue + 2 * READ_CHUNK as u64 + WRITE_BUF as u64;
        Ok(MemoryPlan {
            fixed: 0,
            pipeline,
            buffer: budget.memory.saturating_sub(pipeline) / 2,
        })
    }

    /// Names and lengths of the sequences, from a usable `.fai` or a pass
    /// over the input.
    fn names(&self, ctx: &Context<'_>) -> Result<(Vec<(String, u64)>, bool)> {
        let fai = fai_path(&self.src);
        if fai.exists() && detect(&self.src)? != Compression::Gzip {
            let fai = Fai::load(&self.src)?;
            let names = fai
                .entries()
                .iter()
                .map(|e| (e.name.clone(), e.len as u64))
                .collect();
            return Ok((names, true));
        }
        let s = Self::piece_sizing(ctx.budget);
        let mut input = open_text(&self.src, &ctx.pool, s.limit)?;
        let mut names = Vec::new();
        let mut parser = Parser::new();
        let mut chunk = vec![0u8; READ_CHUNK];
        loop {
            let n = input.text.read(&mut chunk)?;
            if n == 0 {
                break;
            }
            parser.feed(&chunk[..n], &mut |ev| {
                match ev {
                    Event::Start(name) => names.push((name, 0)),
                    Event::Bases(b) => {
                        names.last_mut().expect("bases follow a name").1 += b.len() as u64
                    }
                    Event::End => {}
                }
                Ok(())
            })?;
            ctx.progress(
                "scan",
                input.read.load(std::sync::atomic::Ordering::Relaxed),
                input.file_len,
            );
        }
        parser.finish(&mut |_| Ok(()))?;
        Ok((names, false))
    }

    pub(crate) fn run(&self, ctx: &Context<'_>) -> Result<Report> {
        let (names, from_fai) = self.names(ctx)?;
        let version =
            if twobit_base_size(names.iter().map(|(n, l)| (n.as_str(), *l))) > u32::MAX as u64 {
                1
            } else {
                0
            };
        let mut outputs = Outputs::new();
        let file = outputs.create(&self.dst, Format::TwoBit)?;
        let writer = TwoBitWriter::create(
            BufWriter::with_capacity(WRITE_BUF, file),
            names.iter().map(|(n, _)| n.clone()).collect(),
            version,
        )?;
        let s = Self::piece_sizing(ctx.budget);
        let mut input = open_text(&self.src, &ctx.pool, s.limit)?;
        let mut sink = Sink {
            writer,
            names: &names,
            from_fai,
            record: 0,
            current: TwoBitRecord::new(),
            packed: Packed::new(ctx.memory.buffer, ctx.tmp_parent(&self.dst)),
            pending: VecDeque::new(),
        };
        let mut exec = ctx.pool.ordered::<Option<TwoBitPiece>>(s.limit);
        let mut parser = Parser::new();
        let piece_size = piece_len(ctx.budget);
        let mut piece: Vec<u8> = Vec::with_capacity(piece_size);
        let mut piece_start = 0u64;
        let mut started = 0usize;
        let mut chunk = vec![0u8; READ_CHUNK];
        let submit = |exec: &mut dyn Ordered<Option<TwoBitPiece>>,
                      sink: &mut Sink<'_>,
                      task: Option<(u64, Vec<u8>)>|
         -> Result<()> {
            while exec.full() {
                let done = exec.next(true)?.expect("pending result");
                sink.take(done)?;
            }
            sink.pending.push_back(task.is_none());
            match task {
                Some((start, bases)) => exec.submit(Box::new(move || {
                    TwoBitPiece::encode(start, &bases).map(Some)
                }))?,
                None => exec.submit(Box::new(|| Ok(None)))?,
            }
            while let Some(done) = exec.next(false)? {
                sink.take(done)?;
            }
            Ok(())
        };
        loop {
            let n = input.text.read(&mut chunk)?;
            let mut handle = |ev: Event<'_>| -> Result<()> {
                match ev {
                    Event::Start(name) => {
                        let Some((expected, _)) = names.get(started) else {
                            return Err(stale(from_fai, "more sequences than expected".into()));
                        };
                        if &name != expected {
                            return Err(stale(
                                from_fai,
                                format!("sequence {name:?}, expected {expected:?}"),
                            ));
                        }
                        started += 1;
                        piece_start = 0;
                    }
                    Event::Bases(mut b) => {
                        while !b.is_empty() {
                            let take = (piece_size - piece.len()).min(b.len());
                            piece.extend_from_slice(&b[..take]);
                            b = &b[take..];
                            if piece.len() == piece_size {
                                let full =
                                    std::mem::replace(&mut piece, Vec::with_capacity(piece_size));
                                submit(&mut exec, &mut sink, Some((piece_start, full)))?;
                                piece_start += piece_size as u64;
                            }
                        }
                    }
                    Event::End => {
                        if !piece.is_empty() {
                            let last =
                                std::mem::replace(&mut piece, Vec::with_capacity(piece_size));
                            submit(&mut exec, &mut sink, Some((piece_start, last)))?;
                        }
                        submit(&mut exec, &mut sink, None)?;
                    }
                }
                Ok(())
            };
            if n == 0 {
                parser.finish(&mut handle)?;
                break;
            }
            parser.feed(&chunk[..n], &mut handle)?;
            ctx.progress(
                "encode",
                input.read.load(std::sync::atomic::Ordering::Relaxed),
                input.file_len,
            );
        }
        while let Some(done) = exec.next(true)? {
            sink.take(done)?;
        }
        if sink.record != names.len() {
            return Err(stale(
                from_fai,
                format!(
                    "{} sequences in the input, {} expected",
                    sink.record,
                    names.len()
                ),
            ));
        }
        let spilled = sink.packed.spilled;
        let out = sink.writer.finish()?;
        out.into_inner().map_err(|e| e.into_error())?;
        Ok(Report {
            files: outputs.commit()?,
            records: names.len() as u64,
            bytes_read: input.read.load(std::sync::atomic::Ordering::Relaxed),
            spilled,
            elapsed: Default::default(),
            warnings: Vec::new(),
        })
    }
}

fn stale(from_fai: bool, detail: String) -> Error {
    if from_fai {
        Error::invalid(
            "fasta",
            format!("{detail}; the .fai is out of date, rebuild it"),
        )
    } else {
        Error::invalid("fasta", detail)
    }
}

struct Sink<'n> {
    writer: TwoBitWriter<BufWriter<File>>,
    names: &'n [(String, u64)],
    from_fai: bool,
    record: usize,
    current: TwoBitRecord,
    packed: Packed,
    /// Whether each pending result ends a record.
    pending: VecDeque<bool>,
}

impl Sink<'_> {
    fn take(&mut self, done: Option<TwoBitPiece>) -> Result<()> {
        let ends = self.pending.pop_front().expect("pending entry per result");
        match done {
            Some(piece) => {
                self.current.push(&piece)?;
                self.packed.write(&piece.packed)?;
            }
            None => {
                debug_assert!(ends);
                let Some((name, len)) = self.names.get(self.record) else {
                    return Err(stale(self.from_fai, "more sequences than expected".into()));
                };
                if self.current.len() != *len {
                    return Err(stale(
                        self.from_fai,
                        format!("{name}: {} bases, {len} expected", self.current.len()),
                    ));
                }
                let record = std::mem::take(&mut self.current);
                let mut reader = self.packed.reader()?;
                self.writer.add(name, &record, &mut reader)?;
                drop(reader);
                self.packed.reset()?;
                self.record += 1;
            }
        }
        Ok(())
    }
}

/// Packed bases of the record being assembled: in memory up to a limit,
/// then in a temporary file.
struct Packed {
    mem: Vec<u8>,
    limit: u64,
    parent: PathBuf,
    file: Option<(TempDir, File)>,
    spilled: u64,
}

impl Packed {
    fn new(limit: u64, parent: PathBuf) -> Packed {
        Packed {
            mem: Vec::new(),
            limit,
            parent,
            file: None,
            spilled: 0,
        }
    }

    fn write(&mut self, bytes: &[u8]) -> Result<()> {
        if self.file.is_none() && (self.mem.len() + bytes.len()) as u64 > self.limit {
            let dir = TempDir::new(&self.parent)?;
            let path = dir.path().join("packed");
            let mut f = File::options()
                .read(true)
                .write(true)
                .create_new(true)
                .open(&path)
                .map_err(|e| with_path(e, &path))?;
            f.write_all(&self.mem)?;
            self.spilled += self.mem.len() as u64;
            self.mem = Vec::new();
            self.file = Some((dir, f));
        }
        match &mut self.file {
            Some((_, f)) => {
                f.write_all(bytes)?;
                self.spilled += bytes.len() as u64;
            }
            None => self.mem.extend_from_slice(bytes),
        }
        Ok(())
    }

    fn reader(&mut self) -> Result<Box<dyn Read + '_>> {
        Ok(match &mut self.file {
            Some((_, f)) => {
                f.seek(SeekFrom::Start(0))?;
                Box::new(std::io::BufReader::with_capacity(READ_CHUNK, f))
            }
            None => Box::new(Cursor::new(&self.mem[..])),
        })
    }

    fn reset(&mut self) -> Result<()> {
        self.mem.clear();
        if let Some((_, f)) = &mut self.file {
            f.set_len(0)?;
            f.seek(SeekFrom::Start(0))?;
        }
        Ok(())
    }
}

/// What the FASTA parser reports.
pub(crate) enum Event<'a> {
    /// A record starts; its name is the header's first word.
    Start(String),
    /// Bases of the current record, whitespace removed.
    Bases(&'a [u8]),
    /// The current record ends.
    End,
}

/// A streaming FASTA parser over chunks of any size.
pub(crate) struct Parser {
    at_line_start: bool,
    in_header: bool,
    header: Vec<u8>,
    in_record: bool,
    line: u64,
    clean: Vec<u8>,
}

impl Parser {
    pub(crate) fn new() -> Parser {
        Parser {
            at_line_start: true,
            in_header: false,
            header: Vec::new(),
            in_record: false,
            line: 1,
            clean: Vec::new(),
        }
    }

    fn start(&mut self, emit: &mut dyn FnMut(Event<'_>) -> Result<()>) -> Result<()> {
        if self.in_record {
            emit(Event::End)?;
        }
        let text = std::str::from_utf8(&self.header).map_err(|_| Error::Format {
            line: self.line,
            detail: "sequence name is not UTF-8".into(),
        })?;
        let name = text
            .split(|c: char| c.is_ascii_whitespace())
            .next()
            .unwrap_or("");
        if name.is_empty() {
            return Err(Error::Format {
                line: self.line,
                detail: "empty sequence name".into(),
            });
        }
        emit(Event::Start(name.to_owned()))?;
        self.header.clear();
        self.in_record = true;
        Ok(())
    }

    pub(crate) fn feed(
        &mut self,
        mut chunk: &[u8],
        emit: &mut dyn FnMut(Event<'_>) -> Result<()>,
    ) -> Result<()> {
        while !chunk.is_empty() {
            if self.in_header {
                match memchr::memchr(b'\n', chunk) {
                    Some(i) => {
                        self.header.extend_from_slice(&chunk[..i]);
                        chunk = &chunk[i + 1..];
                        self.start(emit)?;
                        self.in_header = false;
                        self.at_line_start = true;
                        self.line += 1;
                    }
                    None => {
                        self.header.extend_from_slice(chunk);
                        return Ok(());
                    }
                }
                continue;
            }
            if self.at_line_start && chunk[0] == b'>' {
                self.in_header = true;
                chunk = &chunk[1..];
                continue;
            }
            // Sequence text up to the end of the line (or chunk).
            let (body, rest, ended) = match memchr::memchr(b'\n', chunk) {
                Some(i) => (&chunk[..i], &chunk[i + 1..], true),
                None => (chunk, &chunk[chunk.len()..], false),
            };
            self.clean.clear();
            self.clean
                .extend(body.iter().copied().filter(|c| !c.is_ascii_whitespace()));
            if !self.clean.is_empty() {
                if !self.in_record {
                    return Err(Error::Format {
                        line: self.line,
                        detail: "sequence before the first '>' header".into(),
                    });
                }
                if let Some(c) = self.clean.iter().find(|c| !c.is_ascii_graphic()) {
                    return Err(Error::Format {
                        line: self.line,
                        detail: format!("byte {c:#04x} in a sequence"),
                    });
                }
                let clean = std::mem::take(&mut self.clean);
                emit(Event::Bases(&clean))?;
                self.clean = clean;
            }
            self.at_line_start = ended;
            if ended {
                self.line += 1;
            }
            chunk = rest;
        }
        Ok(())
    }

    pub(crate) fn finish(&mut self, emit: &mut dyn FnMut(Event<'_>) -> Result<()>) -> Result<()> {
        if self.in_header {
            self.start(emit)?;
            self.in_header = false;
        }
        if self.in_record {
            emit(Event::End)?;
            self.in_record = false;
        }
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn parse(text: &[u8], chunk: usize) -> Vec<(String, Vec<u8>)> {
        let mut p = Parser::new();
        let mut out: Vec<(String, Vec<u8>)> = Vec::new();
        let mut emit = |ev: Event<'_>| {
            match ev {
                Event::Start(n) => out.push((n, Vec::new())),
                Event::Bases(b) => out.last_mut().unwrap().1.extend_from_slice(b),
                Event::End => {}
            }
            Ok(())
        };
        for c in text.chunks(chunk) {
            p.feed(c, &mut emit).unwrap();
        }
        p.finish(&mut emit).unwrap();
        out
    }

    #[test]
    fn parses_across_any_chunking() {
        let text = b">a desc\nAC GT\r\nnn\n>b\n>c x y\nACGT\nA";
        let want = vec![
            ("a".to_string(), b"ACGTnn".to_vec()),
            ("b".to_string(), Vec::new()),
            ("c".to_string(), b"ACGTA".to_vec()),
        ];
        for chunk in [1, 2, 3, 7, 100] {
            assert_eq!(parse(text, chunk), want, "chunk {chunk}");
        }
    }
}
