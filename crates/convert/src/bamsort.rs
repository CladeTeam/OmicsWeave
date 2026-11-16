//! BAM → coordinate-sorted BAM (+ BAI / CSI), as `samtools sort` sorts.

use std::collections::HashMap;
use std::fs::File;
use std::io::{BufReader, BufWriter, Read};
use std::ops::Range;
use std::path::{Path, PathBuf};
use std::sync::Arc;

use omicsweave_bam::{BAM_MAGIC, Record, Reference, Scanner, flag};
use omicsweave_bgzf::{BinningIndex, IndexFormat, Writer, try_csi_levels};
use omicsweave_core::exec::Ordered;
use omicsweave_core::{Error, Result};

use crate::budget::{self, Budget, MemoryPlan, Sizing, WRITE_BUF, sizing, text_batch};
use crate::conversion::{Context, Report};
use crate::format::Format;
use crate::pipeline::input::{READ_CHUNK, with_path};
use crate::pipeline::output::Outputs;
use crate::pipeline::records::block_sizing;
use crate::pipeline::sort::{Sorter, Spill, read_or_end};
use crate::sam::BamIndex;

/// What the record scanner holds: up to 4 MiB of consumed text before it
/// drops it, in a buffer that may have doubled past that.
const SCAN_TEXT: u64 = 8 << 20;

/// The smallest sort buffer (`samtools sort` refuses less than 1 MiB per
/// thread): below it nearly every batch would become a run file.
const MIN_BUFFER: u64 = 1 << 20;

/// `SAM_FORMAT_VERSION`, the `VN` of an `@HD` line htslib adds.
const SAM_VERSION: &[u8] = b"1.6";

/// A reference of the written header: its name and length.
type OutRef = (Vec<u8>, u32);

/// BAM to BAM sorted by coordinate, optionally indexed.
///
/// Records are ordered as `samtools sort` orders them: by reference, those
/// without one (`tid` -1) last, then by `pos + 1` (as an unsigned number
/// of 63 bits), then forward before reverse strand, ties kept in input
/// order. Unmapped reads placed at their mate's position sort there. The
/// header is the one `samtools sort --no-PG` writes: the text as htslib
/// parses and rebuilds it, `@HD` first with `SO:coordinate` (added with
/// `VN:1.6` when missing) and without `GO` or `SS`, an `@SQ` line for every
/// reference. Records are copied as they are stored, with the `bin` htslib
/// recomputes for a record with a CIGAR; they are cut into BGZF blocks where
/// htslib cuts them, so at level 0 the file is byte for byte the one
/// `samtools sort -l 0 --no-PG` writes, and an index has the bins and
/// chunks `samtools index` builds. The buffer spills to sorted runs when
/// full; the output does not depend on it, nor on the threads.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct BamSort {
    pub src: PathBuf,
    pub dst: PathBuf,
    /// libdeflate level, 0 to 12.
    pub level: i32,
    pub index: BamIndex,
}

/// Records in flight between the reader and the sort, as raw bytes (each
/// with its `block_size`) of about one text batch.
fn batch_sizing(budget: &Budget) -> Sizing {
    // The batch, the records copied out of it and their sort entries.
    sizing(
        budget.threads,
        3 * text_batch(budget) as u64,
        budget.memory / 4,
    )
}

fn pipeline(budget: &Budget) -> u64 {
    // Blocks inflating ahead of the scanner and compressing behind the
    // writer.
    batch_sizing(budget).per_queue
        + 2 * block_sizing(budget).per_queue
        + SCAN_TEXT
        + READ_CHUNK as u64
        + WRITE_BUF as u64
}

impl BamSort {
    pub fn new(src: impl Into<PathBuf>, dst: impl Into<PathBuf>) -> BamSort {
        BamSort {
            src: src.into(),
            dst: dst.into(),
            level: omicsweave_bgzf::DEFAULT_LEVEL,
            index: BamIndex::None,
        }
    }

    pub(crate) fn plan(&self, budget: &Budget) -> Result<MemoryPlan> {
        omicsweave_bgzf::check_level(self.level)?;
        let pipeline = pipeline(budget);
        budget::require(budget, pipeline.saturating_add(MIN_BUFFER), "sorting a BAM")?;
        Ok(MemoryPlan {
            fixed: 0,
            pipeline,
            buffer: budget.memory - pipeline,
        })
    }

    pub(crate) fn run(&self, ctx: &Context<'_>) -> Result<Report> {
        let file = File::open(&self.src).map_err(|e| with_path(e, &self.src))?;
        let file_len = file.metadata()?.len();
        let blocks = block_sizing(ctx.budget);
        let mut scan = Scanner::new_streaming(
            BufReader::with_capacity(READ_CHUNK, file),
            Box::new(ctx.pool.ordered(blocks.limit)),
        );
        let header = scan.header()?;
        let refs = header.references();
        if let Some(r) = refs.iter().find(|r| r.len > i32::MAX as u32) {
            return Err(Error::invalid(
                "bam sort",
                format!(
                    "reference {} of {} bases is longer than BAM positions reach",
                    r.name, r.len
                ),
            ));
        }
        let (text, out_refs) = sorted_header(header.text(), refs)?;
        let index = index_format(self.index, &out_refs)?;
        let n_refs = refs.len();
        let mut outputs = Outputs::new();
        let file = outputs.create(&self.dst, Format::Bam)?;

        let tmp = ctx.tmp_parent(&self.dst);
        let mut sorter = Sorter::<Sorted>::new(&ctx.pool, ctx.memory.buffer, &tmp);
        let s = batch_sizing(ctx.budget);
        let target = text_batch(ctx.budget);
        let mut exec = ctx.pool.ordered::<Vec<Sorted>>(s.limit);
        let mut records = 0u64;
        let mut batch = Vec::with_capacity(target);
        let mut first = 0u64;
        loop {
            let range = scan.next_record()?;
            let done = range.is_none();
            if let Some(range) = range {
                let data = scan.slice(range);
                batch.extend_from_slice(&(data.len() as u32).to_le_bytes());
                batch.extend_from_slice(data);
                records += 1;
            }
            if batch.len() >= target || (done && !batch.is_empty()) {
                while exec.full() {
                    sorter.extend(exec.next(true)?.expect("pending"))?;
                }
                let bytes = std::mem::replace(&mut batch, Vec::with_capacity(target));
                let seq = first;
                first = records;
                exec.submit(Box::new(move || parse_batch(bytes, seq, n_refs)))?;
                while let Some(parsed) = exec.next(false)? {
                    sorter.extend(parsed)?;
                }
                ctx.progress("read", scan.compressed_read(), file_len);
            }
            if done {
                break;
            }
        }
        while let Some(parsed) = exec.next(true)? {
            sorter.extend(parsed)?;
        }
        drop(exec);
        let bytes_read = scan.compressed_read();
        drop(scan);
        let spilled = sorter.spilled();
        ctx.progress("sort", bytes_read, file_len);
        let mut merged = sorter.finish()?;

        let out = BufWriter::with_capacity(WRITE_BUF, file);
        let exec: Box<dyn Ordered<Vec<u8>>> = Box::new(ctx.pool.ordered(blocks.limit));
        let mut writer = if index.is_some() {
            Writer::new(out, self.level, exec)?
        } else {
            Writer::new_streaming(out, self.level, exec)?
        };
        writer.write(&header_bytes(&text, &out_refs))?;
        writer.flush_block()?;
        let mut ix = index
            .as_ref()
            .map(|(format, _)| BinningIndex::new(*format, out_refs.len(), writer.tell()));
        let mut n = 0u64;
        while let Some(r) = merged.next_record()? {
            writer.flush_try(r.data.bytes().len())?;
            if let Some(ix) = ix.as_mut() {
                ix.amend_last(writer.tell());
            }
            writer.write(r.data.bytes())?;
            if let Some(ix) = ix.as_mut() {
                let d = &r.data.bytes()[4..];
                let flags = (u32_at(d, 12) >> 16) as u16;
                ix.push(
                    i32_at(d, 0),
                    i32_at(d, 4) as i64,
                    r.end,
                    writer.tell(),
                    flags & flag::UNMAPPED == 0,
                )?;
            }
            n += 1;
            if n.is_multiple_of(1 << 16) {
                ctx.progress("write", n, records);
            }
        }
        drop(merged);
        writer.flush_block()?;
        let end = writer.tell();
        if let (Some(ix), Some((kind, suffix))) = (ix, index) {
            let done = writer.finish()?;
            done.out.into_inner().map_err(|e| e.into_error())?;
            let bytes = ix.finish(end, &done.blocks)?;
            let format = match kind {
                IndexFormat::Bai => Format::Bai,
                IndexFormat::Csi { .. } => Format::Csi,
            };
            outputs.write(&suffixed(&self.dst, suffix), format, &bytes)?;
        } else {
            writer
                .finish_streaming()?
                .into_inner()
                .map_err(|e| e.into_error())?;
        }
        Ok(Report {
            files: outputs.commit()?,
            records,
            bytes_read,
            spilled,
            elapsed: Default::default(),
            warnings: Vec::new(),
        })
    }
}

fn suffixed(path: &Path, suffix: &str) -> PathBuf {
    let mut s = path.as_os_str().to_owned();
    s.push(suffix);
    PathBuf::from(s)
}

/// The index to write for references `refs`, and its suffix.
fn index_format(index: BamIndex, refs: &[OutRef]) -> Result<Option<(IndexFormat, &'static str)>> {
    let max_len = refs.iter().map(|r| r.1 as u64).max().unwrap_or(0);
    let fits_bai = max_len <= 1 << 29;
    Ok(match index {
        BamIndex::None => None,
        BamIndex::Bai if !fits_bai => {
            return Err(Error::invalid(
                "bam index",
                format!("a reference of {max_len} bases does not fit a BAI index; use CSI"),
            ));
        }
        BamIndex::Bai => Some((IndexFormat::Bai, ".bai")),
        BamIndex::Auto if fits_bai => Some((IndexFormat::Bai, ".bai")),
        BamIndex::Auto => {
            let (min_shift, levels) = try_csi_levels(max_len, 14)?;
            Some((IndexFormat::Csi { min_shift, levels }, ".csi"))
        }
        BamIndex::Csi { min_shift } => {
            let (min_shift, levels) = try_csi_levels(max_len, min_shift)?;
            Some((IndexFormat::Csi { min_shift, levels }, ".csi"))
        }
    })
}

#[inline]
fn u32_at(b: &[u8], at: usize) -> u32 {
    u32::from_le_bytes(b[at..at + 4].try_into().unwrap())
}

#[inline]
fn i32_at(b: &[u8], at: usize) -> i32 {
    u32_at(b, at) as i32
}

/// A record in the sort buffer: its bytes (from its `block_size` on) and
/// where it ends on the reference, for the index.
struct Sorted {
    /// `tid` as unsigned: no reference (-1) after every reference.
    tid: u32,
    /// `(pos + 1) mod 2^63`, shifted left, with the reverse-strand flag.
    pos: u64,
    /// Input order.
    seq: u64,
    end: i64,
    data: RecordData,
}

/// Input batches enter the sorter whole. Sharing their immutable bytes
/// avoids an allocation and payload copy for every buffered record.
enum RecordData {
    Owned(Vec<u8>),
    Shared(Arc<Vec<u8>>, Range<u32>),
}

impl RecordData {
    fn bytes(&self) -> &[u8] {
        match self {
            Self::Owned(bytes) => bytes,
            Self::Shared(bytes, range) => &bytes[range.start as usize..range.end as usize],
        }
    }

    fn heap_bytes(&self) -> usize {
        match self {
            Self::Owned(bytes) => bytes.capacity(),
            Self::Shared(bytes, range) => {
                // All ranges of a batch enter the buffer together. Charge
                // excess capacity and the Arc allocation once, to its first
                // row; the full-buffer reservation lasts through spilling.
                (range.end - range.start) as usize
                    + if range.start == 0 {
                        bytes.capacity() - bytes.len()
                            + size_of::<Vec<u8>>()
                            + 2 * size_of::<usize>()
                    } else {
                        0
                    }
            }
        }
    }
}

impl Sorted {
    /// The sort fields of `data`, a record with its `block_size`.
    fn new(seq: u64, end: i64, data: RecordData) -> Sorted {
        let d = &data.bytes()[4..];
        let flags = (u32_at(d, 12) >> 16) as u16;
        let pos = ((i32_at(d, 4) as i64 + 1) as u64) & (u64::MAX >> 1);
        Sorted {
            tid: u32_at(d, 0),
            pos: pos << 1 | u64::from(flags & flag::REVERSE != 0),
            seq,
            end,
            data,
        }
    }
}

impl Spill for Sorted {
    type Key = (u32, u64, u64);
    fn key(&self) -> Self::Key {
        (self.tid, self.pos, self.seq)
    }
    fn mem_size(&self) -> usize {
        self.data.heap_bytes()
    }
    fn encode(&self, out: &mut Vec<u8>) {
        out.extend_from_slice(&self.seq.to_le_bytes());
        out.extend_from_slice(&self.end.to_le_bytes());
        out.extend_from_slice(self.data.bytes());
    }
    fn decode(input: &mut impl Read) -> Result<Option<Sorted>> {
        let mut h = [0u8; 20];
        if !read_or_end(input, &mut h)? {
            return Ok(None);
        }
        let seq = u64::from_le_bytes(h[..8].try_into().unwrap());
        let end = i64::from_le_bytes(h[8..16].try_into().unwrap());
        let size = u32_at(&h, 16) as usize;
        let mut data = vec![0u8; 4 + size];
        data[..4].copy_from_slice(&h[16..]);
        input.read_exact(&mut data[4..])?;
        Ok(Some(Sorted::new(seq, end, RecordData::Owned(data))))
    }
}

/// `hts_reg2bin(beg, end, 14, 5)`, with htslib's arithmetic shifts of
/// negative positions.
fn reg2bin(beg: i64, end: i64) -> u16 {
    let end = end - 1;
    let mut s = 14;
    let mut t = ((1i64 << 15) - 1) / 7;
    let mut l = 5;
    while l > 0 {
        if beg >> s == end >> s {
            return (t + (beg >> s)) as u16;
        }
        l -= 1;
        s += 3;
        t -= 1 << (3 * l);
    }
    0
}

/// Checks the records of `bytes` (each after its `block_size`) as
/// `sam_read1` checks them, sets the `bin` it recomputes, and gives them
/// input numbers from `first` on.
fn parse_batch(mut bytes: Vec<u8>, first: u64, n_refs: usize) -> Result<Vec<Sorted>> {
    if bytes.len() > u32::MAX as usize {
        return Err(Error::Unsupported("BAM sort batch exceeds 4 GiB"));
    }
    if bytes.is_empty() {
        return Ok(Vec::new());
    }
    if 4 + u32_at(&bytes, 0) as usize == bytes.len() {
        let (_, end) = prepare_record(&mut bytes, 0, n_refs)?;
        bytes.shrink_to_fit();
        return Ok(vec![Sorted::new(first, end, RecordData::Owned(bytes))]);
    }
    let mut ends = Vec::new();
    let mut at = 0;
    while at < bytes.len() {
        let (total, end) = prepare_record(&mut bytes, at, n_refs)?;
        ends.push(end);
        at += total;
    }
    bytes.shrink_to_fit();
    let bytes = Arc::new(bytes);
    let mut out = Vec::with_capacity(ends.len());
    at = 0;
    for end in ends {
        let total = 4 + u32_at(&bytes, at) as usize;
        out.push(Sorted::new(
            first + out.len() as u64,
            end,
            RecordData::Shared(Arc::clone(&bytes), at as u32..(at + total) as u32),
        ));
        at += total;
    }
    Ok(out)
}

fn prepare_record(bytes: &mut [u8], at: usize, n_refs: usize) -> Result<(usize, i64)> {
    let total = 4 + u32_at(bytes, at) as usize;
    let (end, bin) = {
        let rec = Record::parse(&bytes[at + 4..at + total], n_refs)?;
        let bin =
            (u32_at(bytes, at + 16) & 0xffff != 0).then(|| reg2bin(rec.pos as i64, rec.end()));
        (rec.end(), bin)
    };
    if let Some(bin) = bin {
        bytes[at + 14..at + 16].copy_from_slice(&bin.to_le_bytes());
    }
    Ok((total, end))
}

fn header_error(detail: impl Into<String>) -> Error {
    Error::invalid("BAM header", detail)
}

/// One header line: its type and its fields (`@CO`: the comment).
struct Line {
    kind: [u8; 2],
    tags: Vec<Vec<u8>>,
}

impl Line {
    /// The first field of key `key`.
    fn find(&self, key: &[u8; 2]) -> Option<usize> {
        self.tags.iter().position(|t| t.starts_with(key))
    }

    /// The value of the last field of key `key`.
    fn last(&self, key: &[u8; 2]) -> Option<&[u8]> {
        self.tags
            .iter()
            .rev()
            .find(|t| t.starts_with(key))
            .map(|t| &t[3..])
    }

    /// Sets the first field of key `key`, or appends one.
    fn set(&mut self, key: &[u8; 2], value: &[u8]) {
        let field = [&key[..], b":", value].concat();
        match self.find(key) {
            Some(i) => self.tags[i] = field,
            None => self.tags.push(field),
        }
    }

    fn remove(&mut self, key: &[u8; 2]) {
        if let Some(i) = self.find(key) {
            self.tags.remove(i);
        }
    }
}

/// `strtoll(s, NULL, 10)`: leading blanks, a sign, digits; 0 for none,
/// saturated when out of range.
fn strtoll(s: &[u8]) -> i64 {
    let mut i = s
        .iter()
        .position(|c| !matches!(c, b' ' | b'\t' | b'\n' | b'\x0b' | b'\x0c' | b'\r'))
        .unwrap_or(s.len());
    let negative = s.get(i) == Some(&b'-');
    if matches!(s.get(i), Some(b'-' | b'+')) {
        i += 1;
    }
    let mut v: i128 = 0;
    while let Some(d) = s.get(i).filter(|c| c.is_ascii_digit()) {
        v = (v * 10 + i128::from(d - b'0')).min(1 << 64);
        i += 1;
    }
    let v = if negative { -v } else { v };
    v.clamp(i64::MIN as i128, i64::MAX as i128) as i64
}

/// `sam_hdr_sanitise`, which `sam_hdr_read` applies to a BAM's text: up to
/// the first NUL every line must start with `@`, and the last one gets the
/// newline it lacks.
fn sanitise(text: &[u8]) -> Result<std::borrow::Cow<'_, [u8]>> {
    let end = memchr::memchr(0, text).unwrap_or(text.len());
    let mut last = b'\n';
    let mut lno = 0;
    for &c in &text[..end] {
        if last == b'\n' {
            lno += 1;
            if c != b'@' {
                return Err(header_error(format!("malformed header text at line {lno}")));
            }
        }
        last = c;
    }
    if text.is_empty() || last == b'\n' {
        return Ok(text.into());
    }
    let mut fixed = text.to_vec();
    match fixed.get_mut(end) {
        Some(c) => *c = b'\n',
        None => fixed.push(b'\n'),
    }
    Ok(fixed.into())
}

/// The lines of a header text as `sam_hrecs_parse_lines` reads them.
fn parse_lines(text: &[u8]) -> Result<Vec<Line>> {
    let len = text.len();
    if len < 3 {
        if len == 0 || text[0] == 0 {
            return Ok(Vec::new());
        }
        return Err(header_error("header line too short"));
    }
    let mut lines = Vec::new();
    let mut i = 0;
    let mut lno = 1;
    while i < len - 3 && text[i] != 0 {
        let h = &text[i..];
        if h[0] != b'@' {
            return Err(header_error(format!(
                "header line {lno} does not start with '@'"
            )));
        }
        if !h[1].is_ascii_alphabetic() || !h[2].is_ascii_alphabetic() {
            return Err(header_error(format!(
                "header line {lno} does not have a two character key"
            )));
        }
        if h[3] == b'\n' {
            return Err(header_error(format!("header line {lno} has no fields")));
        }
        let kind = [h[1], h[2]];
        let body = &h[3..];
        let mut tags = Vec::new();
        let used = if kind == *b"CO" {
            if body[0] != b'\t' {
                return Err(header_error(format!("header line {lno}: missing tab")));
            }
            let j = 1 + body[1..]
                .iter()
                .position(|&c| c == 0 || c == b'\n')
                .unwrap_or(body.len() - 1);
            tags.push(body[1..j].to_vec());
            j
        } else {
            let mut k = 0;
            loop {
                if k == body.len() || body[k] != b'\t' {
                    return Err(header_error(format!("header line {lno}: missing tab")));
                }
                k += 1;
                let j = k + body[k..]
                    .iter()
                    .position(|&c| c == 0 || c == b'\n' || c == b'\t')
                    .unwrap_or(body.len() - k);
                if j - k < 3 || body[k + 2] != b':' {
                    return Err(header_error(format!(
                        "header line {lno}: malformed key:value pair"
                    )));
                }
                tags.push(body[k..j].to_vec());
                k = j;
                if !(k < body.len() && body[k] != 0 && body[k] != b'\n') {
                    break k;
                }
            }
        };
        lines.push(Line { kind, tags });
        i += 3 + used + 1;
        lno += 1;
    }
    Ok(lines)
}

/// The header text and references `samtools sort --no-PG` writes for a BAM
/// whose header is `text` and `refs`: htslib's header read and parse
/// (`sam_hdr_sanitise`, `sam_hdr_fill_hrecs`), the sort order set
/// (`set_sort_order`), and the text rebuilt (`sam_hrecs_rebuild_text`).
fn sorted_header(text: &[u8], refs: &[Reference]) -> Result<(Vec<u8>, Vec<OutRef>)> {
    struct Ref {
        name: Vec<u8>,
        len: i64,
        /// The `@SQ` line, by position in parse order.
        line: Option<usize>,
    }
    let mut names: HashMap<Vec<u8>, usize> = HashMap::new();
    let mut all: Vec<Ref> = Vec::with_capacity(refs.len());
    for r in refs {
        let name = r.name.as_bytes().to_vec();
        if names.insert(name.clone(), all.len()).is_some() {
            return Err(header_error(format!(
                "duplicate reference {:?} in the reference list",
                r.name
            )));
        }
        all.push(Ref {
            name,
            len: r.len as i64,
            line: None,
        });
    }
    let mut lines = parse_lines(&sanitise(text)?)?;
    /// `AN` names of reference `idx`, unless the names are taken.
    fn add_altnames(names: &mut HashMap<Vec<u8>, usize>, idx: usize, list: Option<&[u8]>) {
        for alt in list.into_iter().flat_map(|l| l.split(|&c| c == b',')) {
            if !alt.is_empty() {
                names.entry(alt.to_vec()).or_insert(idx);
            }
        }
    }
    for (no, line) in lines.iter_mut().enumerate() {
        match &line.kind {
            b"SQ" => {
                let mut len = -1i64;
                let mut lens_differ = false;
                for t in line.tags.iter().filter(|t| t.starts_with(b"LN")) {
                    let v = strtoll(&t[3..]);
                    lens_differ |= len != -1 && len != v;
                    len = v;
                }
                let Some(name) = line.last(b"SN").map(<[u8]>::to_vec) else {
                    return Err(header_error("an @SQ line has no SN field"));
                };
                let shown = String::from_utf8_lossy(&name).into_owned();
                if len <= 0 {
                    return Err(header_error(format!(
                        "@SQ line {shown:?} has an invalid or no LN field"
                    )));
                }
                if lens_differ {
                    return Err(header_error(format!(
                        "@SQ line {shown:?} has LN fields of different values"
                    )));
                }
                let altnames = line.last(b"AN").map(<[u8]>::to_vec);
                if let Some(&idx) = names.get(&name) {
                    if all[idx].line.is_none() {
                        all[idx].line = Some(no);
                        if len != all[idx].len {
                            line.set(b"LN", all[idx].len.to_string().as_bytes());
                        }
                        add_altnames(&mut names, idx, altnames.as_deref());
                        continue;
                    }
                    if all[idx].name == name {
                        return Err(header_error(format!(
                            "duplicate @SQ line {shown:?} in the header"
                        )));
                    }
                }
                let idx = all.len();
                names.insert(name.clone(), idx);
                all.push(Ref {
                    name,
                    len,
                    line: Some(no),
                });
                add_altnames(&mut names, idx, altnames.as_deref());
            }
            b"RG" if line.find(b"ID").is_none() => {
                return Err(header_error("an @RG line has no ID field"));
            }
            b"PG" if line.find(b"ID").is_none() => {
                return Err(header_error("a @PG line has no ID field"));
            }
            _ => {}
        }
    }
    // The global line order: the first @HD goes to the top.
    let mut order: Vec<Line> = Vec::with_capacity(lines.len() + 1);
    for line in lines {
        if line.kind == *b"HD" && order.first().is_some_and(|f| f.kind != *b"HD") {
            order.insert(0, line);
        } else {
            order.push(line);
        }
    }
    // Lines for references without one, after the last @SQ line.
    let at = order
        .iter()
        .rposition(|l| l.kind == *b"SQ")
        .map_or(order.len(), |i| i + 1);
    let mut missing = Vec::new();
    for r in all.iter().filter(|r| r.line.is_none()) {
        if r.name.is_empty() {
            return Err(header_error("a reference has an empty name"));
        }
        missing.push(Line {
            kind: *b"SQ",
            tags: vec![
                [&b"SN:"[..], &r.name].concat(),
                format!("LN:{}", r.len).into_bytes(),
            ],
        });
    }
    if !missing.is_empty() {
        // Move the trailing non-SQ lines once, independently of the number
        // of synthesized references.
        order.splice(at..at, missing);
    }
    match order.first_mut().filter(|l| l.kind == *b"HD") {
        Some(hd) => {
            hd.set(b"SO", b"coordinate");
            hd.remove(b"GO");
            hd.remove(b"SS");
        }
        None => order.insert(
            0,
            Line {
                kind: *b"HD",
                tags: vec![
                    [&b"VN:"[..], SAM_VERSION].concat(),
                    b"SO:coordinate".to_vec(),
                ],
            },
        ),
    }
    let mut out = Vec::new();
    for line in &order {
        out.push(b'@');
        out.extend_from_slice(&line.kind);
        for t in &line.tags {
            out.push(b'\t');
            out.extend_from_slice(t);
        }
        out.push(b'\n');
    }
    let refs = all
        .into_iter()
        .map(|r| (r.name, r.len.min(u32::MAX as i64) as u32))
        .collect();
    Ok((out, refs))
}

/// The BAM header block: magic, text, references.
fn header_bytes(text: &[u8], refs: &[OutRef]) -> Vec<u8> {
    let mut out = Vec::with_capacity(12 + text.len() + refs.len() * 16);
    out.extend_from_slice(&BAM_MAGIC);
    out.extend_from_slice(&(text.len() as u32).to_le_bytes());
    out.extend_from_slice(text);
    out.extend_from_slice(&(refs.len() as u32).to_le_bytes());
    for (name, len) in refs {
        out.extend_from_slice(&(name.len() as u32 + 1).to_le_bytes());
        out.extend_from_slice(name);
        out.push(0);
        out.extend_from_slice(&len.to_le_bytes());
    }
    out
}

#[cfg(test)]
mod tests {
    use super::*;

    fn refs(list: &[(&str, u32)]) -> Vec<Reference> {
        list.iter()
            .map(|&(name, len)| Reference {
                name: name.into(),
                len,
            })
            .collect()
    }

    fn text_of(text: &[u8], list: &[(&str, u32)]) -> String {
        String::from_utf8(sorted_header(text, &refs(list)).unwrap().0).unwrap()
    }

    #[test]
    fn hd_moves_first_and_takes_the_sort_order() {
        assert_eq!(
            text_of(
                b"@SQ\tSN:a\tLN:5\n@HD\tVN:1.4\tGO:query\tSO:unsorted\tSS:x\n",
                &[("a", 5)]
            ),
            "@HD\tVN:1.4\tSO:coordinate\n@SQ\tSN:a\tLN:5\n"
        );
        assert_eq!(
            text_of(b"@CO\tfree text\there", &[]),
            "@HD\tVN:1.6\tSO:coordinate\n@CO\tfree text\there\n"
        );
    }

    #[test]
    fn references_without_lines_get_them_after_the_last_sq() {
        assert_eq!(
            text_of(
                b"@SQ\tSN:b\tLN:9\n@RG\tID:r\n\0\0\0",
                &[("a", 5), ("b", 7), ("c", 3)]
            ),
            "@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:b\tLN:7\n@SQ\tSN:a\tLN:5\n\
             @SQ\tSN:c\tLN:3\n@RG\tID:r\n"
        );
        assert_eq!(
            text_of(b"", &[("a", 5)]),
            "@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:a\tLN:5\n"
        );
    }

    #[test]
    fn many_missing_references_keep_the_trailing_header_order() {
        let references: Vec<_> = (0..4096)
            .map(|i| Reference {
                name: format!("r{i}"),
                len: 5,
            })
            .collect();
        let mut source = b"@HD\tVN:1.6\tSO:unsorted\n@SQ\tSN:r0\tLN:5\n".to_vec();
        let mut tail = Vec::new();
        for i in 0..8192 {
            tail.extend_from_slice(format!("@RG\tID:g{i}\n").as_bytes());
        }
        source.extend_from_slice(&tail);
        let (actual, out_refs) = sorted_header(&source, &references).unwrap();
        let mut expected = b"@HD\tVN:1.6\tSO:coordinate\n".to_vec();
        for reference in &references {
            expected.extend_from_slice(format!("@SQ\tSN:{}\tLN:5\n", reference.name).as_bytes());
        }
        expected.extend_from_slice(&tail);
        assert_eq!(actual, expected);
        assert_eq!(
            out_refs,
            references
                .iter()
                .map(|r| (r.name.as_bytes().to_vec(), r.len))
                .collect::<Vec<_>>()
        );
    }

    #[test]
    fn malformed_headers_fail_as_in_htslib() {
        for bad in [
            &b"@SQ\tSN:a\n"[..],
            b"@SQ\tLN:5\n",
            b"@SQ\tSN:a\tLN:5\n@SQ\tSN:a\tLN:5\n",
            b"@RG\tSM:x\n",
            b"@HD\tVN:1.6\t\n",
            b"@HD VN:1.6\n",
            b"SQ\n",
            b"@HD\tVN:1.6\n\n",
            b"@HD\tVN:1.6\nab\n",
        ] {
            assert!(sorted_header(bad, &refs(&[("a", 5)])).is_err(), "{bad:?}");
        }
        assert!(sorted_header(b"", &refs(&[("a", 5), ("a", 6)])).is_err());
    }

    #[test]
    fn bins_match_htslib() {
        assert_eq!(reg2bin(-1, 0), 4680);
        assert_eq!(reg2bin(0, 1), 4681);
        assert_eq!(reg2bin(0, 1 << 14), 4681);
        assert_eq!(reg2bin(0, (1 << 14) + 1), 585);
        assert_eq!(reg2bin(0, 1 << 29), 0);
    }

    #[test]
    fn batches_share_checked_payload_and_charge_its_full_capacity() {
        let mut bytes = Vec::with_capacity(16384);
        for i in 0..32 {
            let mut record = vec![0; 32];
            record[4..8].copy_from_slice(&(i % 11i32).to_le_bytes());
            record[8..12].copy_from_slice(&(0xbeef_1e02u32).to_le_bytes());
            record[12..16].copy_from_slice(&(1u32 | (((i & 1) * 16) as u32) << 16).to_le_bytes());
            record[16..20].copy_from_slice(&1u32.to_le_bytes());
            record[20..28].fill(255);
            record.extend_from_slice(b"r\0");
            record.extend_from_slice(&16u32.to_le_bytes());
            record.extend_from_slice(&[0x10, 30]);
            bytes.extend_from_slice(&(record.len() as u32).to_le_bytes());
            bytes.extend_from_slice(&record);
        }
        let records = parse_batch(bytes, 7, 1).unwrap();
        assert_eq!(records.len(), 32);
        let RecordData::Shared(first, _) = &records[0].data else {
            panic!("a parsed batch shares its input");
        };
        let expected_heap = first.capacity() + size_of::<Vec<u8>>() + 2 * size_of::<usize>();
        assert_eq!(
            records.iter().map(Spill::mem_size).sum::<usize>(),
            expected_heap
        );
        let one = parse_batch(records[0].data.bytes().to_vec(), 100, 1).unwrap();
        assert_eq!(one.len(), 1);
        assert_eq!(one[0].seq, 100);
        assert_eq!(one[0].end, records[0].end);
        assert_eq!(one[0].data.bytes(), records[0].data.bytes());
        assert!(parse_batch(Vec::new(), 0, 0).unwrap().is_empty());
        let mut expected = Vec::new();
        for (i, row) in records.iter().enumerate() {
            let RecordData::Shared(owner, _) = &row.data else {
                panic!("a parsed batch shares its input");
            };
            assert!(Arc::ptr_eq(first, owner));
            let d = row.data.bytes();
            let rec = Record::parse(&d[4..], 1).unwrap();
            assert_eq!(row.seq, 7 + i as u64);
            assert_eq!(row.end, i as i64 % 11 + 1);
            assert_eq!(
                u16::from_le_bytes(d[14..16].try_into().unwrap()),
                reg2bin(rec.pos as i64, rec.end())
            );
            let mut encoded = Vec::new();
            row.encode(&mut encoded);
            let restored = Sorted::decode(&mut encoded.as_slice()).unwrap().unwrap();
            assert_eq!(restored.key(), row.key());
            assert_eq!(restored.end, row.end);
            assert_eq!(restored.data.bytes(), d);
            expected.push((row.key(), row.end, d.to_vec()));
        }
        expected.sort_unstable_by_key(|row| row.0);
        let dir = tempfile::tempdir().unwrap();
        let pool = omicsweave_core::exec::Pool::new(1);
        let mut sorter = Sorter::new(&pool, 256, dir.path());
        sorter.extend(records).unwrap();
        assert!(sorter.spilled() > 0);
        let mut merged = sorter.finish().unwrap();
        for (key, end, data) in expected {
            let row = merged.next_record().unwrap().unwrap();
            assert_eq!(row.key(), key);
            assert_eq!(row.end, end);
            assert_eq!(row.data.bytes(), data);
        }
        assert!(merged.next_record().unwrap().is_none());
        drop(merged);
        assert_eq!(std::fs::read_dir(dir.path()).unwrap().count(), 0);
    }
}
