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
    Csi,
    Bai,
    Bed,
    GzipBed,
    BigWig,
    BigBed,
    Cool,
    Mcool,
    Hic,
    Sidecar,
}

/// Leading bytes enough to tell every format apart (a BGZF header is the
/// longest signature).
const SNIFF: usize = 18;

/// The HDF5 superblock signature, at the start of a cooler file.
const HDF5_MAGIC: &[u8; 8] = b"\x89HDF\r\n\x1a\n";

const SAM_HEADER_LINES: [&[u8]; 5] = [b"@HD\t", b"@SQ\t", b"@RG\t", b"@PG\t", b"@CO\t"];

impl Format {
    pub const ALL: [Format; 22] = [
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
        Format::Csi,
        Format::Bai,
        Format::Bed,
        Format::GzipBed,
        Format::BigWig,
        Format::BigBed,
        Format::Cool,
        Format::Mcool,
        Format::Hic,
        Format::Sidecar,
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
            Format::Csi => "csi",
            Format::Bai => "bai",
            Format::Bed => "bed",
            Format::GzipBed => "gzip-bed",
            Format::BigWig => "bigwig",
            Format::BigBed => "bigbed",
            Format::Cool => "cool",
            Format::Mcool => "mcool",
            Format::Hic => "hic",
            Format::Sidecar => "sidecar",
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
    /// `.bam`, `.csi`, `.bai`, `.bed`, `.bw`, `.bigwig`, `.bb`, `.bigbed`,
    /// `.cool`, `.mcool`, `.hic`, `.sidecar`; text formats may end in `.gz`
    /// or `.bgz` (compressed FASTA is taken to be bgzip FASTA, the storage
    /// that is written).
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
            ("2bit", false) => Format::TwoBit,
            ("fai", false) => Format::Fai,
            ("gzi", false) => Format::Gzi,
            ("bcf", false) => Format::Bcf,
            ("bam", false) => Format::Bam,
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
    /// gzip stream is BGZF), and BED otherwise. `.fai` and `.gzi` files
    /// have no signature and are not told apart.
    pub fn detect(path: &Path) -> Result<Format> {
        let head = read_head(&mut File::open(path).map_err(|e| with_path(e, path))?)?;
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
                Content::Bcf | Content::Bam | Content::Csi | Content::Bed => Format::Bed,
            });
        }
        let text = read_head(&mut MultiGzDecoder::new(
            File::open(path).map_err(|e| with_path(e, path))?,
        ))?;
        Ok(match Content::of(&text) {
            Content::Bcf => Format::Bcf,
            Content::Bam => Format::Bam,
            Content::Csi => Format::Csi,
            Content::Vcf => Format::GzipVcf,
            Content::Sam => Format::GzipSam,
            Content::Fasta if omicsweave_bgzf::is_bgzf(&head) => Format::BgzipFasta,
            Content::Fasta => Format::GzipFasta,
            Content::Bed => Format::GzipBed,
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
        } else {
            Content::Bed
        }
    }
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
