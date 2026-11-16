//! bedMethyl → bigBed (every column) and → bigWig (one column).
//!
//! Two dialects, told apart by their field count:
//!
//! - **modkit** `pileup` (18 fields): BED9 (`name` is the modification
//!   code, with the motif when there are several; `score` is `Nvalid_cov`;
//!   strand `+`, `-` or `.` when strands are combined), then `Nvalid_cov`,
//!   `percent_modified` (`Nmod / Nvalid_cov * 100`), `Nmod`, `Ncanonical`,
//!   `Nother_mod`, `Ndelete`, `Nfail`, `Ndiff`, `Nnocall`. Older modkit
//!   versions separate the last nine with spaces; both read the same.
//! - **ENCODE** WGBS (11 fields): BED9 (`score` is the coverage capped at
//!   1000), then the coverage and the percentage of methylated reads.
//!
//! Values are stored as given: a percentage stays 0-100, and nothing is
//! derived from the counts.

use std::collections::{BTreeSet, HashMap, HashSet};
use std::io::Read;
use std::path::PathBuf;
use std::sync::Arc;
use std::sync::atomic::Ordering;

use omicsweave_bbi::{BedSchema, WriteOptions};
use omicsweave_core::exec::Ordered;
use omicsweave_core::{Error, Result};

use crate::bed::{self, BedJob, ChromSizes, Item, Parsed, read_chrom_sizes};
use crate::budget::{Budget, MemoryPlan, text_batch};
use crate::conversion::{Context, Report};
use crate::pipeline::input::{LineBatch, LineBatches, lines, open_text};
use crate::pipeline::sort::{Sorter, Spill, read_or_end};
use crate::wig::{Registry, bedgraph_sizing, plan_bedgraph, write_bigwig};

/// The two bedMethyl layouts.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum Dialect {
    Modkit,
    Encode,
}

impl Dialect {
    fn name(self) -> &'static str {
        match self {
            Dialect::Modkit => "modkit",
            Dialect::Encode => "ENCODE",
        }
    }

    /// Names of the fields, as messages call them.
    fn columns(self) -> &'static [&'static str] {
        match self {
            Dialect::Modkit => &MODKIT_COLUMNS,
            Dialect::Encode => &ENCODE_COLUMNS,
        }
    }
}

const MODKIT_COLUMNS: [&str; 18] = [
    "chrom",
    "start",
    "end",
    "code",
    "score",
    "strand",
    "thickStart",
    "thickEnd",
    "color",
    "Nvalid_cov",
    "percent_modified",
    "Nmod",
    "Ncanonical",
    "Nother_mod",
    "Ndelete",
    "Nfail",
    "Ndiff",
    "Nnocall",
];

const ENCODE_COLUMNS: [&str; 11] = [
    "chrom",
    "start",
    "end",
    "name",
    "score",
    "strand",
    "thickStart",
    "thickEnd",
    "itemRgb",
    "coverage",
    "percent",
];

/// The autoSql of modkit's bedMethyl, its column names modkit's own.
const MODKIT_AUTOSQL: &str = r#"table bedMethyl
"modkit pileup: modified base calls per reference position"
    (
    string chrom;             "Reference sequence name"
    uint   chromStart;        "Start position (0-based)"
    uint   chromEnd;          "End position (exclusive)"
    string name;              "Modified base code, with the motif when there are several"
    uint   score;             "Valid coverage (Nvalid_cov)"
    char[1] strand;           "+ or -, or . when strands are combined"
    uint   thickStart;        "Start position"
    uint   thickEnd;          "End position"
    uint   reserved;          "Color (255,0,0)"
    uint   Nvalid_cov;        "Valid coverage: Nmod + Nother_mod + Ncanonical"
    float  percent_modified;  "Nmod / Nvalid_cov * 100"
    uint   Nmod;              "Calls of this modification"
    uint   Ncanonical;        "Calls of the canonical base"
    uint   Nother_mod;        "Calls of another modification of the same base"
    uint   Ndelete;           "Reads with a deletion at this position"
    uint   Nfail;             "Calls below the probability threshold"
    uint   Ndiff;             "Reads with another base at this position"
    uint   Nnocall;           "Reads with the canonical base and no call"
    )
"#;

/// The autoSql of ENCODE's bedMethyl, its column names ENCODE's own.
const ENCODE_AUTOSQL: &str = r#"table bedMethyl
"ENCODE bedMethyl: DNA methylation per position"
    (
    string chrom;       "Reference chromosome or scaffold"
    uint   chromStart;  "Start position in chromosome"
    uint   chromEnd;    "End position in chromosome"
    string name;        "Name of item"
    uint   score;       "Score from 0-1000. Capped number of reads"
    char[1] strand;     "+ or - or . for unknown"
    uint   thickStart;  "Start of where display should be thick"
    uint   thickEnd;    "End of where display should be thick"
    uint   reserved;    "Color value R,G,B"
    uint   readCount;   "Number of reads or coverage"
    float  percentMeth; "Percentage of reads that show methylation at this position in the genome"
    )
"#;

/// The column of a bedMethyl item a bigWig stores.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Value {
    /// Column 5: modkit's valid coverage; ENCODE's coverage capped at 1000.
    Score,
    /// Column 10: modkit `Nvalid_cov`; ENCODE coverage.
    Coverage,
    /// Column 11: modkit `percent_modified`; ENCODE percentage, 0-100.
    Percent,
    /// modkit `Nmod`.
    Nmod,
    /// modkit `Ncanonical`.
    Ncanonical,
    /// modkit `Nother_mod`.
    NotherMod,
    /// modkit `Ndelete`.
    Ndelete,
    /// modkit `Nfail`.
    Nfail,
    /// modkit `Ndiff`.
    Ndiff,
    /// modkit `Nnocall`.
    Nnocall,
}

impl Value {
    pub const ALL: [Value; 10] = [
        Value::Score,
        Value::Coverage,
        Value::Percent,
        Value::Nmod,
        Value::Ncanonical,
        Value::NotherMod,
        Value::Ndelete,
        Value::Nfail,
        Value::Ndiff,
        Value::Nnocall,
    ];

    /// `score`, `coverage`, `percent`, `nmod`, `ncanonical`, `nother_mod`,
    /// `ndelete`, `nfail`, `ndiff` or `nnocall`.
    pub fn name(self) -> &'static str {
        match self {
            Value::Score => "score",
            Value::Coverage => "coverage",
            Value::Percent => "percent",
            Value::Nmod => "nmod",
            Value::Ncanonical => "ncanonical",
            Value::NotherMod => "nother_mod",
            Value::Ndelete => "ndelete",
            Value::Nfail => "nfail",
            Value::Ndiff => "ndiff",
            Value::Nnocall => "nnocall",
        }
    }

    /// The value called `name` ([`Value::name`], any case), or by its
    /// column's name: `valid_coverage` / `nvalid_cov` for coverage,
    /// `percent_modified` for percent.
    pub fn from_name(name: &str) -> Result<Value> {
        let n = name.to_ascii_lowercase();
        let alias = match n.as_str() {
            "valid_coverage" | "nvalid_cov" => Some(Value::Coverage),
            "percent_modified" => Some(Value::Percent),
            _ => None,
        };
        alias
            .or_else(|| Value::ALL.into_iter().find(|v| v.name() == n))
            .ok_or_else(|| {
                let names: Vec<&str> = Value::ALL.iter().map(|v| v.name()).collect();
                Error::invalid(
                    "bedmethyl",
                    format!("unknown value {name:?}: use one of {}", names.join(", ")),
                )
            })
    }

    /// The 0-based field the value is in.
    fn field(self) -> usize {
        match self {
            Value::Score => 4,
            Value::Coverage => 9,
            Value::Percent => 10,
            Value::Nmod => 11,
            Value::Ncanonical => 12,
            Value::NotherMod => 13,
            Value::Ndelete => 14,
            Value::Nfail => 15,
            Value::Ndiff => 16,
            Value::Nnocall => 17,
        }
    }
}

/// The strand of the bedMethyl items a bigWig keeps.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Strand {
    /// `+`
    Forward,
    /// `-`
    Reverse,
    /// `.`: strands combined (modkit `--combine-strands`).
    Combined,
}

impl Strand {
    /// `+`, `-` or `.`.
    pub fn from_name(name: &str) -> Result<Strand> {
        match name {
            "+" => Ok(Strand::Forward),
            "-" => Ok(Strand::Reverse),
            "." => Ok(Strand::Combined),
            _ => Err(Error::invalid(
                "bedmethyl",
                format!("unknown strand {name:?}: use \"+\", \"-\" or \".\""),
            )),
        }
    }

    fn byte(self) -> u8 {
        match self {
            Strand::Forward => b'+',
            Strand::Reverse => b'-',
            Strand::Combined => b'.',
        }
    }
}

/// bedMethyl (plain, gzip or BGZF) to bigBed, with its sidecar: every
/// field kept as written (BED9 + 9 for modkit, BED9 + 2 for ENCODE, the
/// space-separated counts of older modkit tab-separated), with the
/// dialect's autoSql. Items are sorted as [`crate::BedToBigBed`] sorts them.
#[derive(Debug, Clone, PartialEq)]
pub struct BedMethylToBigBed {
    pub src: PathBuf,
    pub dst: PathBuf,
    pub chrom_sizes: ChromSizes,
    pub options: WriteOptions,
    /// Also write `<dst>.sidecar`.
    pub sidecar: bool,
}

impl BedMethylToBigBed {
    pub fn new(
        src: impl Into<PathBuf>,
        dst: impl Into<PathBuf>,
        chrom_sizes: ChromSizes,
    ) -> BedMethylToBigBed {
        BedMethylToBigBed {
            src: src.into(),
            dst: dst.into(),
            chrom_sizes,
            options: WriteOptions::default(),
            sidecar: true,
        }
    }

    pub(crate) fn plan(&self, budget: &Budget) -> Result<MemoryPlan> {
        BedJob::plan(&self.chrom_sizes, budget)
    }

    pub(crate) fn run(&self, ctx: &Context<'_>) -> Result<Report> {
        bed::run(
            BedJob {
                src: &self.src,
                dst: &self.dst,
                chrom_sizes: &self.chrom_sizes,
                options: &self.options,
                sidecar: self.sidecar,
                parse: Arc::new(parse_bed),
                paragraphs: false,
                disjoint: false,
                what: "bed",
                schema,
            },
            ctx,
        )
    }
}

/// bedMethyl (plain, gzip or BGZF) to bigWig, with its sidecar: one
/// column ([`Value`], `percent` by default) of the items of one
/// modification code (`mod_code`, column 4; required when the file has
/// several) and, given `strand`, one strand. Values are stored as written,
/// as bedGraph values are ([`crate::BedGraphToBigWig`]); selected items
/// that overlap are an error that names the selection to add.
#[derive(Debug, Clone, PartialEq)]
pub struct BedMethylToBigWig {
    pub src: PathBuf,
    pub dst: PathBuf,
    pub chrom_sizes: ChromSizes,
    pub options: WriteOptions,
    /// Also write `<dst>.sidecar`.
    pub sidecar: bool,
    pub value: Value,
    /// Column 4 of the items kept (`m`, `h`, `a`, a ChEBI code, ...).
    pub mod_code: Option<String>,
    /// The strand of the items kept; `None` keeps all.
    pub strand: Option<Strand>,
}

impl BedMethylToBigWig {
    pub fn new(
        src: impl Into<PathBuf>,
        dst: impl Into<PathBuf>,
        chrom_sizes: ChromSizes,
    ) -> BedMethylToBigWig {
        BedMethylToBigWig {
            src: src.into(),
            dst: dst.into(),
            chrom_sizes,
            options: WriteOptions::default(),
            sidecar: true,
            value: Value::Percent,
            mod_code: None,
            strand: None,
        }
    }

    pub(crate) fn plan(&self, budget: &Budget) -> Result<MemoryPlan> {
        plan_bedgraph(&self.chrom_sizes, budget)
    }

    pub(crate) fn run(&self, ctx: &Context<'_>) -> Result<Report> {
        let sizes = match &self.chrom_sizes {
            ChromSizes::File(p) => Some(read_chrom_sizes(p)?),
            ChromSizes::FromData => None,
        };
        let s = bedgraph_sizing(ctx.budget);
        let input = open_text(&self.src, &ctx.pool, s.limit)?;
        let (read, file_len) = (Arc::clone(&input.read), input.file_len);
        let mut batches = LineBatches::new(input.text, text_batch(ctx.budget), 1);
        let mut exec = ctx.pool.ordered::<Selected>(s.limit);
        let tmp = ctx.tmp_parent(&self.dst);
        let select = Arc::new(Select {
            value: self.value,
            mod_code: self.mod_code.as_ref().map(|m| m.as_bytes().to_vec()),
            strand: self.strand.map(Strand::byte),
        });
        let mut gather = Gather {
            sorter: Sorter::new(&ctx.pool, ctx.memory.buffer, &tmp),
            chroms: Registry::new(sizes.as_ref()),
            codes: BTreeSet::new(),
            dialect: None,
            items: 0,
        };
        while let Some(batch) = batches.next_batch()? {
            while exec.full() {
                gather.take(exec.next(true)?.expect("pending"))?;
            }
            let select = Arc::clone(&select);
            exec.submit(Box::new(move || parse_selected(&batch, &select)))?;
            while let Some(p) = exec.next(false)? {
                gather.take(p)?;
            }
            ctx.progress("read", read.load(Ordering::Relaxed), file_len);
        }
        while let Some(p) = exec.next(true)? {
            gather.take(p)?;
        }
        let codes: Vec<String> = gather
            .codes
            .iter()
            .map(|c| String::from_utf8_lossy(c).into_owned())
            .collect();
        match &self.mod_code {
            None if codes.len() > 1 => {
                return Err(Error::invalid(
                    "bedmethyl",
                    format!(
                        "the items have {} modification codes ({}): choose one with mod",
                        codes.len(),
                        codes.join(", ")
                    ),
                ));
            }
            Some(m) if !codes.contains(m) => {
                return Err(Error::invalid(
                    "bedmethyl",
                    format!(
                        "no item has modification code {m:?}; the codes are {}",
                        codes.join(", ")
                    ),
                ));
            }
            _ => {}
        }
        if gather.items == 0 && !codes.is_empty() {
            let strand = self.strand.map_or(b'?', Strand::byte) as char;
            return Err(Error::invalid(
                "bedmethyl",
                format!("no item is on strand {strand}"),
            ));
        }
        let records = gather.items;
        let spilled = gather.sorter.spilled();
        let (names, chroms) = gather.chroms.into_chroms("bedmethyl")?;
        ctx.progress("sort", read.load(Ordering::Relaxed), file_len);
        let mut merged = gather.sorter.finish()?;
        let mut last: Option<MethRec> = None;
        let any_strand = self.strand.is_none();
        let written = write_bigwig(
            ctx,
            &self.dst,
            chroms,
            &self.options,
            self.sidecar,
            records,
            move || {
                let Some(r) = merged.next_record()? else {
                    return Ok(None);
                };
                if let Some(p) = &last
                    && p.chrom == r.chrom
                    && r.start < p.end
                {
                    let place = format!(
                        "{}:{}-{} (line {}) overlaps {}-{} (line {})",
                        names[r.chrom as usize], r.start, r.end, r.line, p.start, p.end, p.line
                    );
                    return Err(Error::Format {
                        line: r.line,
                        detail: if any_strand && p.strand != r.strand {
                            format!(
                                "{place} on the other strand: choose a strand with strand \
                                 (+ or -)"
                            )
                        } else {
                            format!("{place}: a bigWig holds one value per base")
                        },
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
}

/// Whitespace-separated words.
fn words(line: &[u8]) -> impl Iterator<Item = &[u8]> {
    line.split(|c| c.is_ascii_whitespace())
        .filter(|w| !w.is_empty())
}

/// The fields of a line: split at tabs, older modkit's space-separated
/// counts (after the ninth field) split at spaces too; split at whitespace
/// when the line has no tab.
fn split<'a>(line: &'a [u8], out: &mut Vec<&'a [u8]>) {
    out.clear();
    if memchr::memchr(b'\t', line).is_none() {
        out.extend(words(line));
        return;
    }
    for (i, f) in line.split(|&c| c == b'\t').enumerate() {
        if i >= 9 && memchr::memchr(b' ', f).is_some() {
            out.extend(f.split(|&c| c == b' ').filter(|w| !w.is_empty()));
        } else {
            out.push(f);
        }
    }
}

fn is_count(field: &[u8]) -> bool {
    !field.is_empty() && field.len() <= 19 && field.iter().all(u8::is_ascii_digit)
}

fn is_number(field: &[u8]) -> bool {
    std::str::from_utf8(field)
        .ok()
        .and_then(|s| s.parse::<f64>().ok())
        .is_some_and(f64::is_finite)
}

fn is_strand(field: &[u8]) -> bool {
    matches!(field, b"+" | b"-" | b".")
}

/// The dialect whose shape `f` has, if any.
fn shape(f: &[&[u8]]) -> Option<Dialect> {
    let d = match f.len() {
        18 => Dialect::Modkit,
        11 => Dialect::Encode,
        _ => return None,
    };
    let ok = [1, 2, 4, 6, 7, 9].iter().all(|&i| is_count(f[i]))
        && !f[3].is_empty()
        && is_strand(f[5])
        && is_number(f[10])
        && f[11..].iter().all(|w| is_count(w));
    ok.then_some(d)
}

/// Whether a data line has the shape of a bedMethyl item.
pub(crate) fn is_bedmethyl_line(line: &[u8]) -> bool {
    let mut f = Vec::with_capacity(18);
    split(line, &mut f);
    shape(&f).is_some()
}

/// Whether a line carries no item, as BED reads them.
fn is_skipped(line: &[u8]) -> bool {
    line.is_empty()
        || line[0] == b'#'
        || line.starts_with(b"track")
        || line.starts_with(b"browser")
        || line.iter().all(u8::is_ascii_whitespace)
}

/// Checks the fields of line `no`: their count tells the dialect, and every
/// field holds what the dialect says. Returns the dialect, start and end.
fn check(f: &[&[u8]], no: u64) -> Result<(Dialect, u32, u32)> {
    let bad = |detail: String| Error::Format { line: no, detail };
    let d = match f.len() {
        18 => Dialect::Modkit,
        11 => Dialect::Encode,
        n => {
            return Err(bad(format!(
                "{n} fields; bedMethyl lines have 18 (modkit) or 11 (ENCODE)"
            )));
        }
    };
    let names = d.columns();
    let text = |i: usize| String::from_utf8_lossy(f[i]).into_owned();
    let coordinate = |i: usize| {
        std::str::from_utf8(f[i])
            .ok()
            .and_then(|s| s.parse::<u32>().ok())
            .ok_or_else(|| bad(format!("{} {:?} is not a coordinate", names[i], text(i))))
    };
    let (start, end) = (coordinate(1)?, coordinate(2)?);
    if start >= end {
        return Err(bad(if start == end {
            format!("empty item {start}-{end}")
        } else {
            format!("start {start} after end {end}")
        }));
    }
    if f[3].is_empty() {
        return Err(bad(format!("empty {}", names[3])));
    }
    if !is_strand(f[5]) {
        return Err(bad(format!("strand {:?} is not +, - or .", text(5))));
    }
    for i in [4, 6, 7, 9].into_iter().chain(11..f.len()) {
        if !is_count(f[i]) {
            return Err(bad(format!("{} {:?} is not a count", names[i], text(i))));
        }
    }
    if !is_number(f[10]) {
        return Err(bad(format!("{} {:?} is not a number", names[10], text(10))));
    }
    Ok((d, start, end))
}

/// The bigBed schema of items with `fields` fields.
fn schema(fields: u16) -> Result<BedSchema> {
    let auto_sql = match fields {
        18 => MODKIT_AUTOSQL,
        11 => ENCODE_AUTOSQL,
        _ => return Err(Error::invalid("bedmethyl", "no bedMethyl items")),
    };
    Ok(BedSchema {
        field_count: fields,
        defined_field_count: 9,
        auto_sql: Some(auto_sql.to_owned()),
    })
}

/// A batch of lines as BED items for bigBed: every field after the third,
/// tab-joined.
fn parse_bed(batch: &LineBatch) -> Result<Parsed> {
    let mut out = Parsed {
        chroms: Vec::new(),
        lengths: Vec::new(),
        items: Vec::new(),
        rest: Vec::with_capacity(batch.text.len()),
    };
    let mut local: HashMap<&[u8], u32> = HashMap::new();
    let mut fields: Vec<&[u8]> = Vec::with_capacity(18);
    for (no, line) in lines(batch) {
        let line = line.strip_suffix(b"\r").unwrap_or(line);
        if is_skipped(line) {
            continue;
        }
        split(line, &mut fields);
        let (_, start, end) = check(&fields, no)?;
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
            fields: fields.len() as u16,
            rest: (at as u32, (out.rest.len() - at) as u32),
            line: no,
        });
    }
    Ok(out)
}

/// What the bigWig keeps of the items.
struct Select {
    value: Value,
    mod_code: Option<Vec<u8>>,
    strand: Option<u8>,
}

/// An item kept, its chromosome local to the batch.
struct Kept {
    chrom: u32,
    start: u32,
    end: u32,
    value: f32,
    line: u64,
    strand: u8,
}

/// What a batch of lines gives the bigWig.
struct Selected {
    /// The dialect and the line of the batch's first item.
    dialect: Option<(Dialect, u64)>,
    chroms: Vec<String>,
    items: Vec<Kept>,
    /// Column 4 of every item, kept or not.
    codes: Vec<Vec<u8>>,
}

/// Common files have only a few modification codes; larger dictionaries
/// use hashing rather than comparing every row with every known code.
#[derive(Default)]
struct Codes<'a> {
    values: Vec<&'a [u8]>,
    seen: Option<HashSet<&'a [u8]>>,
    #[cfg(test)]
    comparisons: usize,
}

impl<'a> Codes<'a> {
    fn insert(&mut self, code: &'a [u8]) {
        if let Some(seen) = &mut self.seen {
            if seen.insert(code) {
                self.values.push(code);
            }
            return;
        }
        for &value in &self.values {
            #[cfg(test)]
            {
                self.comparisons += 1;
            }
            if value == code {
                return;
            }
        }
        if self.values.len() == 8 {
            let mut seen: HashSet<_> = self.values.iter().copied().collect();
            seen.insert(code);
            self.seen = Some(seen);
        }
        self.values.push(code);
    }
}

fn parse_selected(batch: &LineBatch, select: &Select) -> Result<Selected> {
    let mut out = Selected {
        dialect: None,
        chroms: Vec::new(),
        items: Vec::new(),
        codes: Vec::new(),
    };
    let mut local: HashMap<&[u8], u32> = HashMap::new();
    let mut codes = Codes::default();
    let mut fields: Vec<&[u8]> = Vec::with_capacity(18);
    for (no, line) in lines(batch) {
        let line = line.strip_suffix(b"\r").unwrap_or(line);
        if is_skipped(line) {
            continue;
        }
        split(line, &mut fields);
        let (d, start, end) = check(&fields, no)?;
        match out.dialect {
            None => {
                if d == Dialect::Encode && select.value.field() > 10 {
                    return Err(Error::invalid(
                        "bedmethyl",
                        format!(
                            "value {} is a modkit column; ENCODE bedMethyl has score, \
                             coverage and percent",
                            select.value.name()
                        ),
                    ));
                }
                out.dialect = Some((d, no));
            }
            Some((first, at)) if first != d => {
                return Err(Error::Format {
                    line: no,
                    detail: format!(
                        "a {} item after the {} item of line {at}",
                        d.name(),
                        first.name()
                    ),
                });
            }
            Some(_) => {}
        }
        let code = fields[3];
        codes.insert(code);
        if select.mod_code.as_deref().is_some_and(|m| m != code)
            || select.strand.is_some_and(|s| fields[5] != [s])
        {
            continue;
        }
        let value = crate::wig::value(fields[select.value.field()], no)?;
        let next = local.len() as u32;
        let chrom = *local.entry(fields[0]).or_insert_with(|| {
            out.chroms
                .push(String::from_utf8_lossy(fields[0]).into_owned());
            next
        });
        out.items.push(Kept {
            chrom,
            start,
            end,
            value,
            line: no,
            strand: fields[5][0],
        });
    }
    out.codes = codes.values.into_iter().map(<[u8]>::to_vec).collect();
    Ok(out)
}

/// An item in the sort buffer; `line` breaks ties in input order.
#[derive(Debug, Clone, Copy)]
struct MethRec {
    chrom: u32,
    start: u32,
    end: u32,
    value: f32,
    line: u64,
    strand: u8,
}

impl Spill for MethRec {
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
        out.push(self.strand);
    }
    fn decode(input: &mut impl Read) -> Result<Option<MethRec>> {
        let mut h = [0u8; 25];
        if !read_or_end(input, &mut h)? {
            return Ok(None);
        }
        let u = |i: usize| u32::from_le_bytes(h[i..i + 4].try_into().unwrap());
        Ok(Some(MethRec {
            chrom: u(0),
            start: u(4),
            end: u(8),
            value: f32::from_bits(u(12)),
            line: u64::from_le_bytes(h[16..24].try_into().unwrap()),
            strand: h[24],
        }))
    }
}

struct Gather<'p, 's> {
    sorter: Sorter<'p, MethRec>,
    chroms: Registry<'s>,
    codes: BTreeSet<Vec<u8>>,
    dialect: Option<(Dialect, u64)>,
    items: u64,
}

impl Gather<'_, '_> {
    fn take(&mut self, p: Selected) -> Result<()> {
        match (self.dialect, p.dialect) {
            (None, d) => self.dialect = d,
            (Some((first, at)), Some((d, no))) if first != d => {
                return Err(Error::Format {
                    line: no,
                    detail: format!(
                        "a {} item after the {} item of line {at}",
                        d.name(),
                        first.name()
                    ),
                });
            }
            _ => {}
        }
        self.codes.extend(p.codes);
        let mut map = Vec::with_capacity(p.chroms.len());
        for name in &p.chroms {
            map.push(self.chroms.id(name)?);
        }
        let mut out = Vec::with_capacity(p.items.len());
        for it in p.items {
            let chrom = map[it.chrom as usize];
            self.chroms.cover(chrom, it.end, it.line)?;
            out.push(MethRec {
                chrom,
                start: it.start,
                end: it.end,
                value: it.value,
                line: it.line,
                strand: it.strand,
            });
        }
        self.items += out.len() as u64;
        self.sorter.extend(out)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    const MODKIT: &[u8] = b"chr1\t10\t11\tm\t5\t+\t10\t11\t255,0,0\t5\t60.00\t3\t2\t0\t0\t1\t0\t0";
    const OLD_MODKIT: &[u8] = b"chr1\t10\t11\tm\t5\t+\t10\t11\t255,0,0\t5 60.00 3 2 0 0 1 0 0";
    const ENCODE: &[u8] = b"chr1\t10\t11\t.\t7\t-\t10\t11\t0,255,0\t7\t43";

    #[test]
    fn dialects_are_told_by_their_fields() {
        let mut f = Vec::new();
        for (line, d) in [
            (MODKIT, Dialect::Modkit),
            (OLD_MODKIT, Dialect::Modkit),
            (ENCODE, Dialect::Encode),
        ] {
            split(line, &mut f);
            assert_eq!(shape(&f), Some(d));
            assert_eq!(check(&f, 1).unwrap(), (d, 10, 11));
        }
        split(b"chr1\t10\t11\tname\t0\t+", &mut f);
        assert_eq!(shape(&f), None);
        assert!(check(&f, 3).is_err());
    }

    #[test]
    fn values_have_names_and_aliases() {
        assert_eq!(
            Value::from_name("percent_modified").unwrap(),
            Value::Percent
        );
        assert_eq!(Value::from_name("Nvalid_cov").unwrap(), Value::Coverage);
        assert_eq!(Value::from_name("nother_mod").unwrap(), Value::NotherMod);
        assert!(Value::from_name("fraction").is_err());
        for v in Value::ALL {
            assert_eq!(Value::from_name(v.name()).unwrap(), v);
        }
    }

    #[test]
    fn code_dictionary_preserves_first_occurrence_and_bounds_linear_work() {
        let names: Vec<_> = (0..16384).map(|i| format!("m{i}")).collect();
        let mut codes = Codes::default();
        for name in &names {
            codes.insert(name.as_bytes());
        }
        for name in names.iter().rev() {
            codes.insert(name.as_bytes());
        }
        let expected: Vec<_> = names.iter().map(|n| n.as_bytes()).collect();
        assert_eq!(codes.values, expected);
        assert_eq!(codes.comparisons, 36);
        let mut small = Codes::default();
        for _ in 0..1000 {
            for name in &names[..8] {
                small.insert(name.as_bytes());
            }
        }
        assert_eq!(small.values, &expected[..8]);
        assert!(small.seen.is_none());
    }
}
