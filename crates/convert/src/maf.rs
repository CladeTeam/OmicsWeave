//! MAF → bigMaf: UCSC's `mafToBigMaf <ref> in.maf | sort -k1,1 -k2,2n`
//! then `bedToBigBed -type=bed3+1 -as=bigMaf.as -tab`, in one pass.
//!
//! Each alignment block becomes one bigBed record: the reference
//! component's chromosome, start and end, and the block rewritten as
//! `mafToBigMaf` writes it (kent's `mafWriteDelimiter`): `a score=` with
//! six decimals, columns padded to the block's widest, `q` and `i` lines
//! after their `s` line, every line ended by `;` instead of a newline.
//! Comments, the header and paragraphs other than `a` blocks are dropped,
//! as kent's MAF reader drops them.
//!
//! The reference is the species (the `src` before its first dot, or its
//! second when the species name has a dot) whose component opens every
//! block, on `+`; by default the species of the first block's first
//! component. Its `srcSize` gives each chromosome's length. Blocks need
//! not be sorted, but may not overlap on the reference.

use std::collections::HashMap;
use std::fs::File;
use std::io::{BufRead, BufReader, Read, Write as _};
use std::path::PathBuf;
use std::sync::Arc;

use flate2::read::MultiGzDecoder;
use omicsweave_bbi::{BedSchema, WriteOptions};
use omicsweave_core::{Error, Result};

use crate::bed::{self, BedJob, ChromSizes, Item, Parsed};
use crate::budget::{Budget, MemoryPlan};
use crate::conversion::{Context, Report};
use crate::pipeline::input::{LineBatch, lines, with_path};

/// UCSC's `bigMaf.as`, as `bedToBigBed -as=bigMaf.as` stores it.
pub const BIGMAF_AUTOSQL: &str = r#"table bedMaf
"Bed3 with MAF block"
    (
    string chrom;      "Reference sequence chromosome or scaffold"
    uint   chromStart; "Start position in chromosome"
    uint   chromEnd;   "End position in chromosome"
    lstring mafBlock;   "MAF block"
    )
"#;

/// MAF (plain, gzip or BGZF) to bigMaf, with its sidecar.
#[derive(Debug, Clone, PartialEq)]
pub struct MafToBigMaf {
    pub src: PathBuf,
    pub dst: PathBuf,
    /// The reference species (kent's `referenceDb`); `None` takes the
    /// species of the first block's first component.
    pub ref_species: Option<String>,
    /// Checked against the reference's `srcSize` when a file; the
    /// `srcSize`s themselves otherwise.
    pub chrom_sizes: ChromSizes,
    pub options: WriteOptions,
    /// Also write `<dst>.sidecar`.
    pub sidecar: bool,
}

impl MafToBigMaf {
    pub fn new(src: impl Into<PathBuf>, dst: impl Into<PathBuf>) -> MafToBigMaf {
        MafToBigMaf {
            src: src.into(),
            dst: dst.into(),
            ref_species: None,
            chrom_sizes: ChromSizes::FromData,
            options: WriteOptions::default(),
            sidecar: true,
        }
    }

    pub(crate) fn plan(&self, budget: &Budget) -> Result<MemoryPlan> {
        BedJob::plan(&self.chrom_sizes, budget)
    }

    pub(crate) fn run(&self, ctx: &Context<'_>) -> Result<Report> {
        let first = read_head(&self.src)?;
        let reference = match &self.ref_species {
            Some(r) if r.is_empty() => {
                return Err(Error::invalid("maf", "the reference species is empty"));
            }
            Some(r) => r.clone().into_bytes(),
            None => {
                let Some(src) = first else {
                    return Err(Error::invalid("maf", "no alignment blocks"));
                };
                match memchr::memchr(b'.', &src) {
                    Some(dot) => src[..dot].to_vec(),
                    None => {
                        return Err(Error::invalid(
                            "maf",
                            format!(
                                "the first component {:?} names no species (species.chromosome); \
                                 name the reference species",
                                String::from_utf8_lossy(&src)
                            ),
                        ));
                    }
                }
            }
        };
        let reference: Arc<[u8]> = reference.into();
        bed::run(
            BedJob {
                src: &self.src,
                dst: &self.dst,
                chrom_sizes: &self.chrom_sizes,
                options: &self.options,
                sidecar: self.sidecar,
                parse: Arc::new(move |batch| parse(batch, &reference)),
                paragraphs: true,
                disjoint: true,
                what: "maf",
                schema: |_| {
                    Ok(BedSchema {
                        field_count: 4,
                        defined_field_count: 3,
                        auto_sql: Some(BIGMAF_AUTOSQL.to_owned()),
                    })
                },
            },
            ctx,
        )
    }
}

/// Whether text starts as a MAF file does.
pub(crate) fn is_maf(head: &[u8]) -> bool {
    head.starts_with(b"##maf")
}

fn words(line: &[u8]) -> impl Iterator<Item = &[u8]> {
    line.split(|c| c.is_ascii_whitespace())
        .filter(|w| !w.is_empty())
}

/// Checks the `##maf` line as kent does (`name=value` words, a version)
/// and returns the `src` of the first block's first component, if any.
fn read_head(path: &std::path::Path) -> Result<Option<Vec<u8>>> {
    let mut file = File::open(path).map_err(|e| with_path(e, path))?;
    let mut magic = [0u8; 2];
    let n = file.read(&mut magic)?;
    let file = File::open(path).map_err(|e| with_path(e, path))?;
    let input: Box<dyn Read> = if n == 2 && magic == [0x1f, 0x8b] {
        Box::new(MultiGzDecoder::new(file))
    } else {
        Box::new(file)
    };
    let mut input = BufReader::with_capacity(64 << 10, input);
    let mut line = Vec::new();
    input.read_until(b'\n', &mut line)?;
    let Some(header) = line.strip_prefix(b"##maf") else {
        return Err(Error::Format {
            line: 1,
            detail: "not MAF: the first line does not start with ##maf".into(),
        });
    };
    let mut version = false;
    for w in words(header) {
        let Some(eq) = memchr::memchr(b'=', w) else {
            return Err(Error::Format {
                line: 1,
                detail: format!("missing = after {}", String::from_utf8_lossy(w)),
            });
        };
        if &w[..eq] == b"version" {
            version = atoi(&w[eq + 1..]) != 0;
        }
    }
    if !version {
        return Err(Error::Format {
            line: 1,
            detail: "no version on the ##maf line".into(),
        });
    }
    let mut in_block = false;
    let mut between = true;
    loop {
        line.clear();
        if input.read_until(b'\n', &mut line)? == 0 {
            return Ok(None);
        }
        if line.first() == Some(&b'#') {
            continue;
        }
        let mut w = words(&line);
        match (w.next(), in_block) {
            (None, _) => {
                in_block = false;
                between = true;
            }
            (Some(b"a"), false) if between => {
                in_block = true;
                between = false;
            }
            (Some(b"s" | b"e"), true) => return Ok(w.next().map(<[u8]>::to_vec)),
            _ => between = false,
        }
    }
}

/// C `atoi`: an optional sign and the digits that follow, 0 for none.
fn atoi(s: &[u8]) -> i64 {
    let (neg, digits) = match s.first() {
        Some(b'-') => (true, &s[1..]),
        Some(b'+') => (false, &s[1..]),
        _ => (false, s),
    };
    let mut v: i64 = 0;
    for &c in digits.iter().take_while(|c| c.is_ascii_digit()) {
        v = v.saturating_mul(10).saturating_add(i64::from(c - b'0'));
    }
    if neg { -v } else { v }
}

/// C `atof`: the longest leading decimal number, 0 for none.
fn atof(s: &[u8]) -> f64 {
    let mut i = 0;
    if matches!(s.first(), Some(b'-' | b'+')) {
        i = 1;
    }
    let int = s[i..].iter().take_while(|c| c.is_ascii_digit()).count();
    i += int;
    let mut frac = 0;
    if s.get(i) == Some(&b'.') {
        frac = s[i + 1..].iter().take_while(|c| c.is_ascii_digit()).count();
        if int + frac > 0 {
            i += 1 + frac;
        }
    }
    if int + frac == 0 {
        return 0.0;
    }
    if matches!(s.get(i), Some(b'e' | b'E')) {
        let mut j = i + 1;
        if matches!(s.get(j), Some(b'-' | b'+')) {
            j += 1;
        }
        let exp = s[j..].iter().take_while(|c| c.is_ascii_digit()).count();
        if exp > 0 {
            i = j + exp;
        }
    }
    std::str::from_utf8(&s[..i])
        .ok()
        .and_then(|t| t.parse().ok())
        .unwrap_or(0.0)
}

/// kent's `digitsBaseTen`.
fn digits(x: i64) -> usize {
    let mut n = if x < 0 { 2 } else { 1 };
    let mut x = x.unsigned_abs();
    while x >= 10 {
        n += 1;
        x /= 10;
    }
    n
}

/// One component of a block (`s` or `e`, with its `q` and `i`).
struct Comp<'a> {
    src: &'a [u8],
    start: i64,
    size: i64,
    strand: u8,
    src_size: i64,
    text: Option<&'a [u8]>,
    quality: Option<&'a [u8]>,
    left_status: u8,
    left_len: i64,
    right_status: u8,
    right_len: i64,
}

/// A block being read: its `a` line's number and score, `r` line and
/// components.
struct Block<'a> {
    line: u64,
    score: f64,
    region: Option<(&'a [u8], i64, &'a [u8])>,
    comps: Vec<Comp<'a>>,
    text_size: usize,
}

fn bad(line: u64, detail: String) -> Error {
    Error::Format { line, detail }
}

/// kent's `lineFileNeedNum`: a whole integer.
fn need_num(w: &[u8], line: u64, what: &str) -> Result<i64> {
    std::str::from_utf8(w)
        .ok()
        .and_then(|t| t.parse::<i64>().ok())
        .filter(|v| i32::try_from(*v).is_ok())
        .ok_or_else(|| {
            bad(
                line,
                format!("{what} {:?} is not a number", String::from_utf8_lossy(w)),
            )
        })
}

impl<'a> Block<'a> {
    fn line(&mut self, first: &[u8], rest: &[&'a [u8]], no: u64) -> Result<()> {
        let count = |n: usize, kind: &str| -> Result<()> {
            if rest.len() + 1 == n {
                Ok(())
            } else {
                Err(bad(
                    no,
                    format!("{} words on an {kind} line; it has {n}", rest.len() + 1),
                ))
            }
        };
        let last_src = |comps: &[Comp<'a>], kind: &str| -> Result<()> {
            match comps.last() {
                Some(c) if c.src == rest[0] => Ok(()),
                Some(c) => Err(bad(
                    no,
                    format!(
                        "{kind} line of {} follows the line of {}",
                        String::from_utf8_lossy(rest[0]),
                        String::from_utf8_lossy(c.src)
                    ),
                )),
                None => Err(bad(no, format!("{kind} line before any s line"))),
            }
        };
        match first {
            b"s" | b"e" => {
                count(7, if first == b"s" { "s" } else { "e" })?;
                let src_size = need_num(rest[4], no, "srcSize")?;
                let start = need_num(rest[1], no, "start")?;
                let mut c = Comp {
                    src: rest[0],
                    start,
                    size: 0,
                    strand: rest[3][0],
                    src_size,
                    text: None,
                    quality: None,
                    left_status: 0,
                    left_len: 0,
                    right_status: 0,
                    right_len: 0,
                };
                if first == b"e" {
                    let len = need_num(rest[2], no, "size")?;
                    (c.left_len, c.right_len) = (len, len);
                    (c.left_status, c.right_status) = (rest[5][0], rest[5][0]);
                } else {
                    c.size = need_num(rest[2], no, "size")?;
                    let text = rest[5];
                    if self.text_size == 0 {
                        self.text_size = text.len();
                    } else if self.text_size != text.len() {
                        return Err(bad(
                            no,
                            format!(
                                "alignment text of {} columns; the block's earlier rows have {}",
                                text.len(),
                                self.text_size
                            ),
                        ));
                    }
                    c.text = Some(text);
                }
                if c.src_size < 0 || c.size < 0 {
                    return Err(bad(no, "negative size".into()));
                }
                if c.start < 0 || c.start + c.size > c.src_size {
                    return Err(bad(
                        no,
                        format!(
                            "{}-{} is outside the source's {} bases",
                            c.start,
                            c.start + c.size,
                            c.src_size
                        ),
                    ));
                }
                self.comps.push(c);
            }
            b"i" => {
                count(6, "i")?;
                last_src(&self.comps, "i")?;
                let c = self.comps.last_mut().expect("checked");
                c.left_status = rest[1][0];
                c.left_len = atoi(rest[2]);
                c.right_status = rest[3][0];
                c.right_len = atoi(rest[4]);
            }
            b"q" => {
                count(3, "q")?;
                last_src(&self.comps, "q")?;
                self.comps.last_mut().expect("checked").quality = Some(rest[1]);
            }
            b"r" => {
                if self.region.is_some() {
                    return Err(bad(no, "a second r line in the block".into()));
                }
                count(4, "r")?;
                let size = need_num(rest[1], no, "r size")?;
                self.region = Some((rest[0], size, rest[2]));
            }
            _ => {}
        }
        Ok(())
    }

    /// The reference's chromosome, start, end and length, and the block
    /// written as `mafToBigMaf` writes it, into `out`.
    fn finish(&self, reference: &[u8], out: &mut Vec<u8>) -> Result<(&'a [u8], u32, u32, u32)> {
        let Some(master) = self.comps.first() else {
            return Err(bad(
                self.line,
                "an alignment block without components".into(),
            ));
        };
        let src = master.src;
        let split = |dot: Option<usize>| dot.filter(|&d| &src[..d] == reference);
        let first = memchr::memchr(b'.', src);
        let second = first.and_then(|d| memchr::memchr(b'.', &src[d + 1..]).map(|e| d + 1 + e));
        let Some(dot) = split(first).or_else(|| split(second)) else {
            return Err(bad(
                self.line,
                format!(
                    "the block begins with {}, not with the reference {}: the reference must \
                     be the first component of every block",
                    String::from_utf8_lossy(src),
                    String::from_utf8_lossy(reference)
                ),
            ));
        };
        if master.strand != b'+' {
            return Err(bad(
                self.line,
                "the reference is not on the + strand".into(),
            ));
        }
        let end = master.start + master.size;
        let fit = |v: i64| u32::try_from(v).expect("checked to fit in i32 and be >= 0");
        self.write(out);
        Ok((
            &src[dot + 1..],
            fit(master.start),
            fit(end),
            fit(master.src_size),
        ))
    }

    /// kent's `mafWriteDelimiter(f, ali, ';')`.
    fn write(&self, out: &mut Vec<u8>) {
        let name =
            |c: &Comp<'a>| -> &'a [u8] { if c.src == b"." { b"defaultName" } else { c.src } };
        let (mut src_w, mut start_w, mut size_w, mut src_size_w) = (0, 0, 0, 0);
        for c in &self.comps {
            src_w = src_w.max(name(c).len());
            start_w = start_w.max(digits(c.start));
            size_w = size_w.max(digits(c.size));
            src_size_w = src_size_w.max(digits(c.src_size));
        }
        let pad = |out: &mut Vec<u8>, s: &[u8], w: usize| {
            out.extend_from_slice(s);
            out.resize(out.len() + w.saturating_sub(s.len()), b' ');
        };
        write!(out, "a score={:.6};", self.score).unwrap();
        if let Some((kind, size, id)) = self.region {
            out.extend_from_slice(b"r ");
            out.extend_from_slice(kind);
            write!(out, " {size} ").unwrap();
            out.extend_from_slice(id);
            out.push(b';');
        }
        for c in &self.comps {
            let src = name(c);
            if c.size == 0 && c.left_status != 0 {
                out.extend_from_slice(b"e ");
                pad(out, src, src_w);
                write!(
                    out,
                    " {:>start_w$} {:>size_w$} {} {:>src_size_w$} {};",
                    c.start, c.left_len, c.strand as char, c.src_size, c.left_status as char
                )
                .unwrap();
                continue;
            }
            out.extend_from_slice(b"s ");
            pad(out, src, src_w);
            write!(
                out,
                " {:>start_w$} {:>size_w$} {} {:>src_size_w$} ",
                c.start, c.size, c.strand as char, c.src_size
            )
            .unwrap();
            out.extend_from_slice(c.text.unwrap_or_default());
            out.push(b';');
            if let Some(q) = c.quality {
                out.extend_from_slice(b"q ");
                pad(out, src, src_w + start_w + size_w + src_size_w + 5);
                out.push(b' ');
                out.extend_from_slice(q);
                out.push(b';');
            }
            if c.left_status != 0 {
                out.extend_from_slice(b"i ");
                pad(out, src, src_w);
                write!(
                    out,
                    " {} {} {} {};",
                    c.left_status as char, c.left_len, c.right_status as char, c.right_len
                )
                .unwrap();
            }
        }
    }
}

/// Where a parse is: between paragraphs, in an `a` block, or passing over
/// a paragraph of another kind.
enum State<'a> {
    Between,
    Block(Block<'a>),
    Skip,
}

/// Adds `block` to `out` as a record, its chromosome a local id.
fn emit<'a>(
    block: &Block<'a>,
    reference: &[u8],
    local: &mut HashMap<&'a [u8], u32>,
    out: &mut Parsed,
) -> Result<()> {
    let at = out.rest.len();
    let (chrom, start, end, len) = block.finish(reference, &mut out.rest)?;
    let next = local.len() as u32;
    let id = *local.entry(chrom).or_insert_with(|| {
        out.chroms.push(String::from_utf8_lossy(chrom).into_owned());
        out.lengths.push(len);
        next
    });
    if out.lengths[id as usize] != len {
        return Err(bad(
            block.line,
            format!(
                "{} is {len} long here but {} in an earlier block",
                String::from_utf8_lossy(chrom),
                out.lengths[id as usize]
            ),
        ));
    }
    out.items.push(Item {
        chrom: id,
        start,
        end,
        fields: 4,
        rest: (at as u32, (out.rest.len() - at) as u32),
        line: block.line,
    });
    Ok(())
}

/// A batch of whole paragraphs as bigMaf records.
fn parse(batch: &LineBatch, reference: &[u8]) -> Result<Parsed> {
    let mut out = Parsed {
        chroms: Vec::new(),
        lengths: Vec::new(),
        items: Vec::new(),
        rest: Vec::with_capacity(batch.text.len()),
    };
    let mut local: HashMap<&[u8], u32> = HashMap::new();
    let mut state = State::Between;
    let mut comps = Vec::new();
    let mut rest: Vec<&[u8]> = Vec::with_capacity(8);
    for (no, line) in lines(batch) {
        if line.first() == Some(&b'#') {
            continue;
        }
        let mut w = words(line);
        let first = w.next();
        state = match (state, first) {
            (State::Between, None) => State::Between,
            (State::Between, Some(b"a")) => {
                let mut score = 0.0;
                for word in w {
                    let Some(eq) = memchr::memchr(b'=', word) else {
                        return Err(bad(
                            no,
                            format!("missing = after {}", String::from_utf8_lossy(word)),
                        ));
                    };
                    if &word[..eq] == b"score" {
                        score = atof(&word[eq + 1..]);
                    }
                }
                State::Block(Block {
                    line: no,
                    score,
                    region: None,
                    comps: std::mem::take(&mut comps),
                    text_size: 0,
                })
            }
            (State::Between | State::Skip, Some(_)) => State::Skip,
            (State::Skip, None) => State::Between,
            (State::Block(block), None) => {
                emit(&block, reference, &mut local, &mut out)?;
                comps = block.comps;
                comps.clear();
                State::Between
            }
            (State::Block(mut block), Some(first)) => {
                rest.clear();
                rest.extend(w);
                block.line(first, &rest, no)?;
                State::Block(block)
            }
        };
    }
    if let State::Block(block) = state {
        emit(&block, reference, &mut local, &mut out)?;
    }
    Ok(out)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn batch(text: &str) -> LineBatch {
        LineBatch {
            text: text.as_bytes().to_vec(),
            first_line: 1,
        }
    }

    fn records(text: &str, reference: &str) -> Result<Vec<String>> {
        let p = parse(&batch(text), reference.as_bytes())?;
        Ok(p.items
            .iter()
            .map(|i| {
                format!(
                    "{}\t{}\t{}\t{}",
                    p.chroms[i.chrom as usize],
                    i.start,
                    i.end,
                    String::from_utf8_lossy(&p.rest[i.rest.0 as usize..][..i.rest.1 as usize])
                )
            })
            .collect())
    }

    #[test]
    fn blocks_are_written_as_kent_writes_them() {
        let text = "##maf version=1\n# a comment\n\
                    a score=23262.0 pass=2\n\
                    s hg38.chr1 10 5 + 1000 AC-GTA\n\
                    s panTro4.chr1   7 6 - 99999 ACTGTA\n\
                    q panTro4.chr1 999999\n\
                    i panTro4.chr1 C 0 I 12\n\
                    e mm10.chr2 100 30 + 5000 I\n\
                    \n\
                    x not a block\n\
                    s ignored\n\
                    \n\
                    a score=-1.5\n\
                    s hg38.chr1.alt 0 2 + 9 AC\n";
        let got = records(text, "hg38").unwrap();
        assert_eq!(
            got,
            [
                "chr1\t10\t15\ta score=23262.000000;\
                 s hg38.chr1     10 5 +  1000 AC-GTA;\
                 s panTro4.chr1   7 6 - 99999 ACTGTA;\
                 q panTro4.chr1               999999;\
                 i panTro4.chr1 C 0 I 12;\
                 e mm10.chr2    100 30 +  5000 I;",
                "chr1.alt\t0\t2\ta score=-1.500000;s hg38.chr1.alt 0 2 + 9 AC;",
            ]
        );
    }

    #[test]
    fn species_with_a_dot_and_bad_blocks() {
        let got = records("a\ns hg19.1.chrX 1 1 + 5 A\n", "hg19.1").unwrap();
        assert!(got[0].starts_with("chrX\t1\t2\t"), "{got:?}");
        for (text, says) in [
            ("a\ns mm10.chr1 1 1 + 5 A\n", "not with the reference"),
            ("a\ns hg38.chr1 1 1 - 5 A\n", "+ strand"),
            ("a\ns hg38.chr1 1 1 + 5 A\ns mm10.c 1 1 + 5 AA\n", "columns"),
            ("a\ns hg38.chr1 4 2 + 5 AA\n", "outside"),
            ("a\ns hg38.chr1 1 1 + 5\n", "words"),
            ("a\ni hg38.chr1 C 0 C 0\n", "before any s line"),
            ("a x\ns hg38.chr1 1 1 + 5 A\n", "missing ="),
            (
                "a\ns hg38.chr1 1 1 + 5 A\n\na\ns hg38.chr1 2 1 + 6 A\n",
                "6 long here but 5",
            ),
        ] {
            let err = records(text, "hg38").unwrap_err().to_string();
            assert!(err.contains(says), "{text:?}: {err}");
        }
    }

    #[test]
    fn c_number_parsing() {
        assert_eq!(atof(b"12.5x"), 12.5);
        assert_eq!(atof(b"-3e2"), -300.0);
        assert_eq!(atof(b"1e"), 1.0);
        assert_eq!(atof(b"x"), 0.0);
        assert_eq!(atoi(b"-12abc"), -12);
        assert_eq!(digits(0), 1);
        assert_eq!(digits(-5), 2);
        assert_eq!(digits(1000), 4);
    }
}
