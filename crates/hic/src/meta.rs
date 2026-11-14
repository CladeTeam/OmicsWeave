//! File-level metadata: the header (chromosomes, resolutions, fragment
//! sites) and the footer's master index (where each chromosome pair's
//! matrix lies). Expected and normalization vectors are not read.

use std::collections::HashMap;

use omicsweave_core::{ContactUnit, Error, ReadAt, Result};

use crate::cursor::Cursor;

/// Versions this crate reads and writes.
pub const MIN_VERSION: i32 = 8;
pub const MAX_VERSION: i32 = 9;

/// Attribute selecting the v9 depth layout of intra-chromosomal blocks.
pub const DEPTH_BASE_ATTRIBUTE: &str = "v9-depth-base";

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Chrom {
    pub name: String,
    pub len: u64,
}

/// Where a chromosome's restriction site positions lie in the header.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct SiteList {
    pub position: u64,
    pub count: u32,
}

/// How far from the diagonal a v9 intra-chromosomal block row lies, as a
/// function of a record's distance from the diagonal `d = |x - y| / √2 /
/// blockBinCount`.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Depth {
    /// `floor(log_base(1 + d))`; base 2 is the default layout.
    Log(u32),
    /// `floor(d / step)`.
    Linear(u32),
}

impl Depth {
    /// The layout a `v9-depth-base` attribute value selects, as the Juicer
    /// tools read it: above 1 a log base, below 0 a linear step, otherwise
    /// (and when absent or unparsable) base 2.
    pub fn from_attribute(value: Option<&str>) -> Depth {
        match value.and_then(|v| v.trim().parse::<i64>().ok()) {
            Some(b) if b > 1 => Depth::Log(b.min(u32::MAX as i64) as u32),
            Some(b) if b < 0 => Depth::Linear(b.unsigned_abs().min(u32::MAX as u64) as u32),
            _ => Depth::Log(2),
        }
    }

    /// The attribute value selecting this layout; `None` for the default.
    pub fn attribute(self) -> Option<String> {
        match self {
            Depth::Log(2) => None,
            Depth::Log(b) => Some(b.to_string()),
            Depth::Linear(s) => Some(format!("-{s}")),
        }
    }

    /// Depth of a record `distance` bins from the diagonal.
    pub fn of(self, distance: i64, block_bins: i64) -> i64 {
        let d = distance.unsigned_abs() as f64 / std::f64::consts::SQRT_2 / block_bins as f64;
        match self {
            Depth::Log(2) => (1.0 + d).log2() as i64,
            Depth::Log(b) => ((1.0 + d).ln() / (b as f64).ln()) as i64,
            Depth::Linear(s) => (d / s as f64) as i64,
        }
    }
}

#[derive(Debug, Clone, PartialEq)]
pub struct Header {
    pub version: i32,
    pub footer_position: u64,
    pub genome: String,
    pub attributes: Vec<(String, String)>,
    pub chroms: Vec<Chrom>,
    /// Base-pair bin sizes.
    pub bp_resolutions: Vec<u32>,
    /// Fragment bin sizes (fragments per bin).
    pub frag_resolutions: Vec<u32>,
    /// One per chromosome; empty without fragment resolutions.
    pub sites: Vec<SiteList>,
}

impl Header {
    pub(crate) fn parse<S: ReadAt + ?Sized>(c: &mut Cursor<'_, S>) -> Result<Header> {
        let magic = c.bytes(4)?;
        if magic != b"HIC\0" {
            let found = u32::from_le_bytes(magic.try_into().unwrap());
            return Err(Error::BadMagic {
                what: "hic header",
                found,
            });
        }
        let version = c.i32()?;
        if !(MIN_VERSION..=MAX_VERSION).contains(&version) {
            return Err(Error::Unsupported(".hic versions other than 8 and 9"));
        }
        let footer_position = c.offset("footer position")?;
        let genome = c.string()?;
        if version > 8 {
            // The normalization vector index: position and length.
            c.skip(16)?;
        }
        let attributes = (0..c.count("attributes")?)
            .map(|_| Ok((c.string()?, c.string()?)))
            .collect::<Result<Vec<_>>>()?;
        let chroms = (0..c.count("chromosomes")?)
            .map(|_| {
                let name = c.string()?;
                let len = if version > 8 {
                    c.i64()?
                } else {
                    c.i32()? as i64
                };
                let len = u64::try_from(len)
                    .map_err(|_| Error::invalid("chromosome", "negative length"))?;
                Ok(Chrom { name, len })
            })
            .collect::<Result<Vec<_>>>()?;
        let mut sizes = |what| -> Result<Vec<u32>> {
            (0..c.count(what)?)
                .map(|_| bin_size(c.i32()?, what))
                .collect()
        };
        let bp_resolutions = sizes("bp resolutions")?;
        let frag_resolutions = sizes("fragment resolutions")?;
        let mut sites = Vec::new();
        if !frag_resolutions.is_empty() {
            for _ in &chroms {
                let count = c.count("fragment sites")?;
                sites.push(SiteList {
                    position: c.position(),
                    count: count as u32,
                });
                c.skip(count as u64 * 4)?;
            }
        }
        Ok(Header {
            version,
            footer_position,
            genome,
            attributes,
            chroms,
            bp_resolutions,
            frag_resolutions,
            sites,
        })
    }

    pub fn attribute(&self, key: &str) -> Option<&str> {
        self.attributes
            .iter()
            .find(|(k, _)| k == key)
            .map(|(_, v)| v.as_str())
    }

    /// The v9 intra-chromosomal block layout.
    pub fn depth(&self) -> Depth {
        Depth::from_attribute(self.attribute(DEPTH_BASE_ATTRIBUTE))
    }

    pub fn resolutions(&self, unit: ContactUnit) -> &[u32] {
        match unit {
            ContactUnit::Bp => &self.bp_resolutions,
            ContactUnit::Frag => &self.frag_resolutions,
        }
    }
}

/// The footer's master index.
#[derive(Debug, Clone, Default, PartialEq)]
pub struct Footer {
    /// (chr1 index, chr2 index) -> file position of the matrix record.
    pub matrices: HashMap<(u32, u32), u64>,
}

fn bin_size(v: i32, what: &'static str) -> Result<u32> {
    u32::try_from(v)
        .ok()
        .filter(|&v| v > 0)
        .ok_or_else(|| Error::invalid(what, format!("bad bin size {v}")))
}

fn parse_key(key: &str) -> Option<(u32, u32)> {
    let (a, b) = key.split_once('_')?;
    Some((a.parse().ok()?, b.parse().ok()?))
}

impl Footer {
    /// Reads the master index at the footer position, after the declared
    /// `nBytesV5` (`int` in v8, `long` in v9), which must fit the file.
    pub(crate) fn parse<S: ReadAt + ?Sized>(
        c: &mut Cursor<'_, S>,
        header: &Header,
        source_len: u64,
    ) -> Result<Footer> {
        c.seek(header.footer_position);
        let v5_bytes = c.count_v9(header.version, "footer length")?;
        let width = if header.version > 8 { 8 } else { 4 };
        let v5_end = header
            .footer_position
            .checked_add(width)
            .and_then(|p| p.checked_add(v5_bytes))
            .ok_or_else(|| Error::invalid("hic footer", "length overflows u64"))?;
        if v5_end > source_len {
            return Err(Error::Truncated {
                what: "hic footer",
                need: v5_end,
                have: source_len,
            });
        }
        let mut matrices = HashMap::new();
        for _ in 0..c.count("master index")? {
            let key = c.string()?;
            let position = c.offset("matrix position")?;
            c.i32()?;
            let pair = parse_key(&key)
                .ok_or_else(|| Error::invalid("master index", format!("bad key {key:?}")))?;
            matrices.insert(pair, position);
        }
        Ok(Footer { matrices })
    }
}
