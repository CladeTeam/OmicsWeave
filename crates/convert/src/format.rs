//! Every file format OmicsWeave reads or writes, told by name, by
//! extension or by content. Conversions are chosen by format.

use std::fmt;
use std::fs::File;
use std::io::{ErrorKind, Read};
use std::path::Path;

use flate2::read::MultiGzDecoder;
use omicsweave_bbi::{BIGBED_MAGIC, BIGWIG_MAGIC};
use omicsweave_core::sidecar::SIDECAR_MAGIC;
use omicsweave_core::{Error, Result};
use omicsweave_fasta::TWOBIT_MAGIC;

use crate::pipeline::input::with_path;
use crate::wig::{is_meta, starts_with_word, trim};

/// A file format. Formats stored more than one way (plain text, gzip,
/// bgzip) have one variant per storage, and indexes and companion files
/// are formats in their own right.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum Format {
    Fasta,
    GzipFasta,
    BgzipFasta,
    TwoBit,
    Fai,
    Gzi,
    Vcf,
    GzipVcf,
    Bcf,
    Sam,
    GzipSam,
    Bam,
    Cram,
    Csi,
    Bai,
    Crai,
    Bed,
    GzipBed,
    BedGraph,
    GzipBedGraph,
    Wig,
    GzipWig,
    BigWig,
    BigBed,
    Cool,
    Mcool,
    Hic,
    Sidecar,
    Pairs,
    GzipPairs,
    JuicerShort,
    GzipJuicerShort,
    BedMethyl,
    GzipBedMethyl,
}

/// Leading bytes enough to tell every format apart (a BGZF header is the
/// longest signature), except the interval and signal texts.
const SNIFF: usize = 18;

/// Leading bytes of text looked at to tell BED, bedGraph and wig apart.
const TEXT_SNIFF: usize = 64 << 10;

/// The HDF5 superblock signature, at the start of a cooler file.
const HDF5_MAGIC: &[u8; 8] = b"\x89HDF\r\n\x1a\n";

/// The first line of a 4DN pairs file begins so (`## pairs format v1.0`).
const PAIRS_MAGIC: &[u8] = b"## pairs format";

const SAM_HEADER_LINES: [&[u8]; 5] = [b"@HD\t", b"@SQ\t", b"@RG\t", b"@PG\t", b"@CO\t"];

impl Format {
    pub const ALL: [Format; 34] = [
        Format::Fasta,
        Format::GzipFasta,
        Format::BgzipFasta,
        Format::TwoBit,
        Format::Fai,
        Format::Gzi,
        Format::Vcf,
        Format::GzipVcf,
        Format::Bcf,
        Format::Sam,
        Format::GzipSam,
        Format::Bam,
        Format::Cram,
        Format::Csi,
        Format::Bai,
        Format::Crai,
        Format::Bed,
        Format::GzipBed,
        Format::BedGraph,
        Format::GzipBedGraph,
        Format::Wig,
        Format::GzipWig,
        Format::BigWig,
        Format::BigBed,
        Format::Cool,
        Format::Mcool,
        Format::Hic,
        Format::Sidecar,
        Format::Pairs,
        Format::GzipPairs,
        Format::JuicerShort,
        Format::GzipJuicerShort,
        Format::BedMethyl,
        Format::GzipBedMethyl,
    ];

    /// The format's name, as the Python layer names it.
    pub fn name(self) -> &'static str {
        match self {
            Format::Fasta => "fasta",
            Format::GzipFasta => "gzip-fasta",
            Format::BgzipFasta => "bgzip-fasta",
            Format::TwoBit => "2bit",
            Format::Fai => "fai",
            Format::Gzi => "gzi",
            Format::Vcf => "vcf",
            Format::GzipVcf => "gzip-vcf",
            Format::Bcf => "bcf",
            Format::Sam => "sam",
            Format::GzipSam => "gzip-sam",
            Format::Bam => "bam",
            Format::Cram => "cram",
            Format::Crai => "crai",
            Format::Csi => "csi",
            Format::Bai => "bai",
            Format::Bed => "bed",
            Format::GzipBed => "gzip-bed",
            Format::BedGraph => "bedgraph",
            Format::GzipBedGraph => "gzip-bedgraph",
            Format::Wig => "wig",
            Format::GzipWig => "gzip-wig",
            Format::BigWig => "bigwig",
            Format::BigBed => "bigbed",
            Format::Cool => "cool",
            Format::Mcool => "mcool",
            Format::Hic => "hic",
            Format::Sidecar => "sidecar",
            Format::Pairs => "pairs",
            Format::GzipPairs => "gzip-pairs",
            Format::JuicerShort => "juicer-short",
            Format::GzipJuicerShort => "gzip-juicer-short",
            Format::BedMethyl => "bedmethyl",
            Format::GzipBedMethyl => "gzip-bedmethyl",
        }
    }

    /// The format called `name`.
    pub fn from_name(name: &str) -> Result<Format> {
        Format::ALL
            .into_iter()
            .find(|f| f.name() == name)
            .ok_or_else(|| Error::invalid("format", format!("unknown format {name:?}")))
    }

    /// The format a file named `path` conventionally holds: `.fa`, `.fasta`,
    /// `.fna`, `.fas`, `.2bit`, `.fai`, `.gzi`, `.vcf`, `.bcf`, `.sam`,
    /// `.bam`, `.cram`, `.crai`, `.csi`, `.bai`, `.bed`, `.bedgraph`, `.bdg`, `.bg`, `.wig`,
    /// `.bw`, `.bigwig`, `.bb`, `.bigbed`, `.cool`, `.mcool`, `.hic`,
    /// `.sidecar`, `.pairs`, `.bedmethyl`; text formats may end in `.gz` or
    /// `.bgz` (compressed FASTA is taken to be bgzip FASTA, the storage that
    /// is written). Juicer short text has no extension of its own.
    pub fn from_extension(path: &Path) -> Result<Format> {
        let name = path
            .file_name()
            .map(|n| n.to_string_lossy().to_ascii_lowercase())
            .unwrap_or_default();
        let (base, compressed) = match name
            .strip_suffix(".gz")
            .or_else(|| name.strip_suffix(".bgz"))
        {
            Some(base) => (base, true),
            None => (name.as_str(), false),
        };
        let ext = base.rsplit_once('.').map_or("", |(_, e)| e);
        let format = match (ext, compressed) {
            ("fa" | "fasta" | "fna" | "fas", false) => Format::Fasta,
            ("fa" | "fasta" | "fna" | "fas", true) => Format::BgzipFasta,
            ("vcf", false) => Format::Vcf,
            ("vcf", true) => Format::GzipVcf,
            ("sam", false) => Format::Sam,
            ("sam", true) => Format::GzipSam,
            ("bed", false) => Format::Bed,
            ("bed", true) => Format::GzipBed,
            ("bedgraph" | "bdg" | "bg", false) => Format::BedGraph,
            ("bedgraph" | "bdg" | "bg", true) => Format::GzipBedGraph,
            ("wig", false) => Format::Wig,
            ("wig", true) => Format::GzipWig,
            ("pairs", false) => Format::Pairs,
            ("pairs", true) => Format::GzipPairs,
            ("bedmethyl", false) => Format::BedMethyl,
            ("bedmethyl", true) => Format::GzipBedMethyl,
            ("2bit", false) => Format::TwoBit,
            ("fai", false) => Format::Fai,
            ("gzi", false) => Format::Gzi,
            ("bcf", false) => Format::Bcf,
            ("bam", false) => Format::Bam,
            ("cram", false) => Format::Cram,
            ("crai", false) => Format::Crai,
            ("csi", false) => Format::Csi,
            ("bai", false) => Format::Bai,
            ("bw" | "bigwig", false) => Format::BigWig,
            ("bb" | "bigbed", false) => Format::BigBed,
            ("cool", false) => Format::Cool,
            ("mcool", false) => Format::Mcool,
            ("hic", false) => Format::Hic,
            ("sidecar", false) => Format::Sidecar,
            _ => {
                return Err(Error::invalid(
                    "format",
                    format!(
                        "cannot tell the format of {} from its name; name it",
                        path.display()
                    ),
                ));
            }
        };
        Ok(format)
    }

    /// The format the file at `path` holds, from its content.
    ///
    /// Binary formats are told by their magic numbers (BGZF ones once
    /// inflated). Text, inflated first when gzip compressed, is VCF when it
    /// starts with `##fileformat=VCF`, SAM with a header line (`@HD`,
    /// `@SQ`, `@RG`, `@PG`, `@CO`), FASTA with `>` (bgzip FASTA when the
    /// gzip stream is BGZF), 4DN pairs with `## pairs format`. Other text is
    /// wig, bedGraph, Juicer short, bedMethyl or BED, told from its first
    /// 64 KiB: wig when a `track` line says `type=wiggle_0` or a line
    /// declares a `fixedStep` / `variableStep` section; bedGraph when a
    /// `track` line says `type=bedGraph`, or when every data line has four
    /// fields, the last a number; Juicer short when every data line is
    /// `str1 chr1 pos1 frag1 str2 chr2 pos2 frag2 [score]`; bedMethyl when
    /// every data line is BED9 followed by modkit's nine counts or ENCODE's
    /// coverage and percentage; BED otherwise. `.fai` and `.gzi` files have
    /// no signature and are not told apart.
    pub fn detect(path: &Path) -> Result<Format> {
        let mut file = File::open(path).map_err(|e| with_path(e, path))?;
        let head = read_head(&mut file)?;
        if head.len() >= 4 {
            let magic = u32::from_le_bytes(head[..4].try_into().unwrap());
            match magic {
                BIGWIG_MAGIC => return Ok(Format::BigWig),
                BIGBED_MAGIC => return Ok(Format::BigBed),
                _ if magic == TWOBIT_MAGIC || magic.swap_bytes() == TWOBIT_MAGIC => {
                    return Ok(Format::TwoBit);
                }
                _ => {}
            }
            match &head[..4] {
                b"HIC\0" => return Ok(Format::Hic),
                b"BAI\x01" => return Ok(Format::Bai),
                b"CRAM" => return Ok(Format::Cram),
                _ => {}
            }
        }
        if head.starts_with(&SIDECAR_MAGIC) {
            return Ok(Format::Sidecar);
        }
        if head.starts_with(HDF5_MAGIC) {
            return Ok(if omicsweave_cool::is_mcool(path)? {
                Format::Mcool
            } else {
                Format::Cool
            });
        }
        if !head.starts_with(&[0x1f, 0x8b]) {
            return Ok(match Content::of(&head) {
                Content::Vcf => Format::Vcf,
                Content::Sam => Format::Sam,
                Content::Fasta => Format::Fasta,
                Content::Pairs => Format::Pairs,
                Content::Bcf | Content::Bam | Content::Csi | Content::Bed => {
                    match Text::of(head, &mut file) {
                        Text::Bed => Format::Bed,
                        Text::BedGraph => Format::BedGraph,
                        Text::Wig => Format::Wig,
                        Text::JuicerShort => Format::JuicerShort,
                        Text::BedMethyl => Format::BedMethyl,
                    }
                }
            });
        }
        let mut inflated = MultiGzDecoder::new(File::open(path).map_err(|e| with_path(e, path))?);
        let text = read_head(&mut inflated)?;
        Ok(match Content::of(&text) {
            Content::Bcf => Format::Bcf,
            Content::Bam => Format::Bam,
            Content::Csi => Format::Csi,
            Content::Vcf => Format::GzipVcf,
            Content::Sam => Format::GzipSam,
            Content::Fasta if omicsweave_bgzf::is_bgzf(&head) => Format::BgzipFasta,
            Content::Fasta => Format::GzipFasta,
            Content::Pairs => Format::GzipPairs,
            Content::Bed => {
                let text = first_line(text, &mut inflated);
                if is_crai(&text) {
                    Format::Crai
                } else {
                    match Text::of(text, &mut inflated) {
                        Text::Bed => Format::GzipBed,
                        Text::BedGraph => Format::GzipBedGraph,
                        Text::Wig => Format::GzipWig,
                        Text::JuicerShort => Format::GzipJuicerShort,
                        Text::BedMethyl => Format::GzipBedMethyl,
                    }
                }
            }
        })
    }
}

impl fmt::Display for Format {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.write_str(self.name())
    }
}

/// What the (inflated) start of a file shows.
enum Content {
    Bcf,
    Bam,
    Csi,
    Vcf,
    Sam,
    Fasta,
    Pairs,
    Bed,
}

impl Content {
    fn of(head: &[u8]) -> Content {
        if head.starts_with(b"BCF\x02") {
            Content::Bcf
        } else if head.starts_with(b"BAM\x01") {
            Content::Bam
        } else if head.starts_with(b"CSI\x01") {
            Content::Csi
        } else if head.starts_with(b"##fileformat=VCF") {
            Content::Vcf
        } else if SAM_HEADER_LINES.iter().any(|tag| head.starts_with(tag)) {
            Content::Sam
        } else if head.starts_with(b">") {
            Content::Fasta
        } else if head.starts_with(PAIRS_MAGIC) {
            Content::Pairs
        } else {
            Content::Bed
        }
    }
}

/// Which of the interval, signal and contact texts a file holds.
enum Text {
    Bed,
    BedGraph,
    Wig,
    JuicerShort,
    BedMethyl,
}

impl Text {
    /// Tells `head` apart, read further from `input` up to [`TEXT_SNIFF`]
    /// bytes; a read error ends the look early.
    fn of(mut head: Vec<u8>, input: &mut impl Read) -> Text {
        let mut whole = false;
        let mut filled = head.len();
        head.resize(TEXT_SNIFF, 0);
        while filled < TEXT_SNIFF {
            match input.read(&mut head[filled..]) {
                Ok(0) => {
                    whole = true;
                    break;
                }
                Ok(n) => filled += n,
                Err(e) if e.kind() == ErrorKind::Interrupted => {}
                Err(_) => break,
            }
        }
        head.truncate(filled);
        if !whole {
            // Only whole lines.
            let end = memchr::memrchr(b'\n', &head).map_or(0, |i| i + 1);
            head.truncate(end);
        }
        let (mut data, mut bedgraph, mut short, mut methyl) = (0usize, true, true, true);
        for line in head.split(|&c| c == b'\n') {
            let line = trim(line);
            if starts_with_word(line, b"track") {
                match track_type(line) {
                    Some(b"wiggle_0") => return Text::Wig,
                    Some(t) if t.eq_ignore_ascii_case(b"bedGraph") => return Text::BedGraph,
                    _ => continue,
                }
            }
            if is_meta(line) {
                continue;
            }
            if starts_with_word(line, b"fixedStep") || starts_with_word(line, b"variableStep") {
                return Text::Wig;
            }
            data += 1;
            bedgraph &= is_bedgraph_line(line);
            short &= crate::pairs::is_short_line(line);
            methyl &= crate::bedmethyl::is_bedmethyl_line(line);
        }
        match (data > 0, bedgraph, short, methyl) {
            (false, ..) => Text::Bed,
            (true, true, _, _) => Text::BedGraph,
            (true, _, true, _) => Text::JuicerShort,
            (true, _, _, true) => Text::BedMethyl,
            _ => Text::Bed,
        }
    }
}

/// The `type=` of a `track` line.
fn track_type(line: &[u8]) -> Option<&[u8]> {
    line.split(|c| c.is_ascii_whitespace())
        .find_map(|w| w.strip_prefix(b"type="))
        .map(|t| t.strip_prefix(b"\"").unwrap_or(t))
        .map(|t| t.strip_suffix(b"\"").unwrap_or(t))
}

/// `chrom start end value`.
fn is_bedgraph_line(line: &[u8]) -> bool {
    let mut w = line
        .split(|c| c.is_ascii_whitespace())
        .filter(|w| !w.is_empty());
    let parses = |f: Option<&[u8]>, number: fn(&str) -> bool| {
        f.and_then(|f| std::str::from_utf8(f).ok())
            .is_some_and(number)
    };
    w.next().is_some()
        && parses(w.next(), |s| s.parse::<u32>().is_ok())
        && parses(w.next(), |s| s.parse::<u32>().is_ok())
        && parses(w.next(), |s| s.parse::<f64>().is_ok())
        && w.next().is_none()
}

/// `head` read further from `input` until it holds its first line (at
/// most 4 KiB); a read error ends it early.
fn first_line(mut head: Vec<u8>, input: &mut impl Read) -> Vec<u8> {
    let mut byte = [0u8];
    while !head.contains(&b'\n') && head.len() < 4096 {
        match input.read(&mut byte) {
            Ok(1) => head.push(byte[0]),
            Ok(_) => break,
            Err(e) if e.kind() == ErrorKind::Interrupted => {}
            Err(_) => break,
        }
    }
    head
}

/// A `.crai` line: six tab-separated integers (the first may be -1).
fn is_crai(text: &[u8]) -> bool {
    let Some(line) = text.split(|&c| c == b'\n').next() else {
        return false;
    };
    let fields: Vec<&[u8]> = line.split(|&c| c == b'\t').collect();
    fields.len() == 6
        && fields.iter().enumerate().all(|(i, f)| {
            let digits = if i == 0 {
                f.strip_prefix(b"-").unwrap_or(f)
            } else {
                f
            };
            !digits.is_empty() && digits.iter().all(u8::is_ascii_digit)
        })
}

/// The first [`SNIFF`] bytes of `input` (fewer when it is shorter).
fn read_head(input: &mut impl Read) -> Result<Vec<u8>> {
    let mut buf = vec![0; SNIFF];
    let mut filled = 0;
    while filled < SNIFF {
        match input.read(&mut buf[filled..]) {
            Ok(0) => break,
            Ok(n) => filled += n,
            Err(e) if e.kind() == ErrorKind::Interrupted => {}
            Err(e) if matches!(e.kind(), ErrorKind::InvalidData | ErrorKind::InvalidInput) => {
                return Err(Error::Inflate("gzip data does not inflate"));
            }
            Err(e) => return Err(e.into()),
        }
    }
    buf.truncate(filled);
    Ok(buf)
}
