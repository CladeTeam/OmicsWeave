//! [`Cool`]: one cooler contact map (a `.cool` file, or one resolution of
//! a `.mcool`), its bin layout and its raw pixels.

use std::collections::HashMap;
use std::ops::Range;
use std::path::Path;

use hdf5_metno as hdf5;
use omicsweave_core::{Error, FileSource, ReadAt, Result};

use crate::column::{Column, PixelScratch, ReadPath};
use crate::h5;

/// The group of an `.mcool` holding one group per resolution.
const RESOLUTIONS: &str = "resolutions";

/// cooler layout versions this reader handles.
const FORMAT_VERSIONS: std::ops::RangeInclusive<i64> = 2..=3;

/// A chromosome of the map's genome.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct CoolChrom {
    pub name: String,
    pub length: u32,
}

/// Which pixels a map stores (cooler's `storage-mode`).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum StorageMode {
    /// Only pixels with `bin1_id <= bin2_id`; the matrix is symmetric.
    SymmetricUpper,
    /// Every nonzero pixel of a (possibly asymmetric) matrix.
    Square,
}

/// Pixel rows decoded by [`Cool::read_pixels`], one column per field.
#[derive(Debug, Clone, Default, PartialEq)]
pub struct PixelBatch {
    pub bin1: Vec<u64>,
    pub bin2: Vec<u64>,
    pub count: Vec<f64>,
}

impl PixelBatch {
    pub fn new() -> PixelBatch {
        PixelBatch::default()
    }

    pub fn len(&self) -> usize {
        self.count.len()
    }

    pub fn is_empty(&self) -> bool {
        self.count.is_empty()
    }

    pub fn clear(&mut self) {
        self.bin1.clear();
        self.bin2.clear();
        self.count.clear();
    }
}

fn open_file(path: &Path) -> Result<(FileSource, hdf5::File)> {
    // Opening the byte source first reports a missing or unreadable file
    // as I/O rather than as an HDF5 failure.
    let source = FileSource::open(path)?;
    let file = hdf5::File::open(path).map_err(h5::err("hdf5 open"))?;
    Ok((source, file))
}

/// Whether a file is a multi-resolution `.mcool` (it has a `/resolutions`
/// group).
pub fn is_mcool(path: &Path) -> Result<bool> {
    let (_, file) = open_file(path)?;
    Ok(file.link_exists(RESOLUTIONS) && file.group(RESOLUTIONS).is_ok())
}

/// The resolutions (bin sizes) of an `.mcool`, ascending; each is the map at
/// group `/resolutions/<res>`. Fails for a file that is not an `.mcool`.
pub fn mcool_resolutions(path: &Path) -> Result<Vec<u32>> {
    let (_, file) = open_file(path)?;
    if !file.link_exists(RESOLUTIONS) {
        return Err(Error::invalid(
            "mcool",
            "no /resolutions group: not an .mcool",
        ));
    }
    let group = file.group(RESOLUTIONS).map_err(h5::err("mcool"))?;
    let mut out = Vec::new();
    for name in group.member_names().map_err(h5::err("mcool"))? {
        let res = name.parse::<u32>().ok().filter(|&r| r > 0).ok_or_else(|| {
            Error::invalid(
                "mcool",
                format!("resolution group {name:?} is not a bin size"),
            )
        })?;
        out.push(res);
    }
    out.sort_unstable();
    Ok(out)
}

/// A cooler contact map: a `.cool` file's root group, or one resolution
/// group of an `.mcool`.
///
/// Opening reads the metadata and the bin table through libhdf5 and checks
/// that bins are fixed-width, then records where every pixel chunk lies in
/// the file. Balancing weights are not read. Pixels are read on demand by
/// row range ([`Cool::read_pixels`]); `&Cool` may be shared by any number
/// of threads reading different ranges at once.
pub struct Cool {
    source: FileSource,
    chroms: Vec<CoolChrom>,
    chrom_offsets: Vec<u64>,
    bin_size: u32,
    storage_mode: StorageMode,
    assembly: Option<String>,
    n_bins: u64,
    n_pixels: u64,
    bin1: Column,
    bin2: Column,
    count: Column,
    // Kept open for the library fallback; the columns' datasets hold the
    // file open as well.
    _file: hdf5::File,
}

impl Cool {
    /// Opens the map at `group`: `"/"` for a `.cool`, `"/resolutions/<res>"`
    /// for one resolution of an `.mcool`.
    pub fn open(path: &Path, group: &str) -> Result<Cool> {
        let (source, file) = open_file(path)?;
        let name = group.trim_start_matches('/');
        let sub;
        let grp: &hdf5::Group = if name.is_empty() {
            &file
        } else {
            if !file.link_exists(name) {
                return Err(Error::invalid("cool group", format!("no group {group:?}")));
            }
            sub = file.group(name).map_err(h5::err("cool group"))?;
            &sub
        };
        check_format(grp)?;

        let bin_size = bin_size(grp)?;
        let storage_mode = match h5::attr_string(grp, "storage-mode")?.as_deref() {
            None | Some("symmetric-upper") => StorageMode::SymmetricUpper,
            Some("square") => StorageMode::Square,
            Some(other) => {
                return Err(Error::invalid(
                    "cool storage-mode",
                    format!("unknown storage mode {other:?}"),
                ));
            }
        };
        let assembly = h5::attr_string(grp, "genome-assembly")?;

        let chroms = read_chroms(grp)?;
        let (chrom_offsets, n_bins) = read_bins(grp, &chroms, bin_size)?;

        let file_len = source.len();
        let bin1 = Column::open(grp, "pixels/bin1_id", file_len, true)?;
        let bin2 = Column::open(grp, "pixels/bin2_id", file_len, true)?;
        let count = Column::open(grp, "pixels/count", file_len, false)?;
        let n_pixels = bin1.len();
        if bin2.len() != n_pixels || count.len() != n_pixels {
            return Err(Error::invalid(
                "cool pixels",
                format!(
                    "column lengths differ: bin1_id {}, bin2_id {}, count {}",
                    n_pixels,
                    bin2.len(),
                    count.len()
                ),
            ));
        }

        Ok(Cool {
            source,
            chroms,
            chrom_offsets,
            bin_size,
            storage_mode,
            assembly,
            n_bins,
            n_pixels,
            bin1,
            bin2,
            count,
            _file: file,
        })
    }

    pub fn chroms(&self) -> &[CoolChrom] {
        &self.chroms
    }

    /// The bin width in base pairs (the last bin of a chromosome is cut at
    /// its end).
    pub fn bin_size(&self) -> u32 {
        self.bin_size
    }

    pub fn storage_mode(&self) -> StorageMode {
        self.storage_mode
    }

    /// The `genome-assembly` attribute as stored (cooler writes `"unknown"`
    /// when none was given).
    pub fn assembly(&self) -> Option<&str> {
        self.assembly.as_deref()
    }

    pub fn n_bins(&self) -> u64 {
        self.n_bins
    }

    pub fn n_pixels(&self) -> u64 {
        self.n_pixels
    }

    /// The first bin id of each chromosome, plus [`Cool::n_bins`] at the end
    /// (`n_chroms + 1` entries).
    pub fn chrom_offsets(&self) -> &[u64] {
        &self.chrom_offsets
    }

    /// A bin id's chromosome index and bin index within that chromosome.
    ///
    /// # Panics
    ///
    /// If `bin_id >= n_bins()`.
    pub fn locate(&self, bin_id: u64) -> (u32, u32) {
        assert!(
            bin_id < self.n_bins,
            "bin id {bin_id} out of range (n_bins {})",
            self.n_bins
        );
        // The last chromosome whose first bin is at or before `bin_id`;
        // empty chromosomes share their offset with the next one and are
        // skipped this way.
        let c = self.chrom_offsets.partition_point(|&o| o <= bin_id) - 1;
        (c as u32, (bin_id - self.chrom_offsets[c]) as u32)
    }

    /// How each pixel column (`bin1_id`, `bin2_id`, `count`) is read.
    pub fn read_paths(&self) -> [ReadPath; 3] {
        [self.bin1.path(), self.bin2.path(), self.count.path()]
    }

    /// Cuts `[0, n_pixels)` into consecutive row ranges of about
    /// `target_rows` rows each, for spreading [`Cool::read_pixels`] over
    /// threads.
    ///
    /// Boundaries fall on chunk boundaries of every chunked column when
    /// those are close together (their least common multiple is at most
    /// `target_rows` or the largest chunk), so no chunk is decoded by two
    /// ranges; otherwise on the largest chunk's boundaries, so only smaller
    /// chunks straddle two ranges.
    pub fn pixel_spans(&self, target_rows: u64) -> Vec<Range<u64>> {
        let n = self.n_pixels;
        if n == 0 {
            return Vec::new();
        }
        let lens: Vec<u64> = [&self.bin1, &self.bin2, &self.count]
            .iter()
            .filter_map(|c| c.chunk_len())
            .collect();
        let largest = lens.iter().copied().max().unwrap_or(1);
        let common = lens
            .iter()
            .try_fold(1u64, |acc, &l| acc.checked_mul(l / gcd(acc, l)));
        let unit = match common {
            Some(l) if l <= target_rows.max(largest) => l,
            _ => largest,
        };
        let step = target_rows.max(1).div_ceil(unit).saturating_mul(unit);
        (0..n.div_ceil(step))
            .map(|k| {
                let start = k * step;
                start..start.saturating_add(step).min(n)
            })
            .collect()
    }

    /// Decodes pixel rows `rows` into `out` (cleared first): bin ids and
    /// raw counts, in file order. Each chunk the range touches is read and
    /// decoded once; `scratch` holds the decoding buffers and is meant to
    /// be kept by the calling thread across calls.
    pub fn read_pixels(
        &self,
        rows: Range<u64>,
        scratch: &mut PixelScratch,
        out: &mut PixelBatch,
    ) -> Result<()> {
        self.check_rows(&rows)?;
        out.clear();
        let n = (rows.end - rows.start) as usize;
        out.bin1.reserve(n);
        out.bin2.reserve(n);
        out.count.reserve(n);
        self.bin1
            .read(&self.source, rows.clone(), scratch, &mut out.bin1)?;
        self.bin2
            .read(&self.source, rows.clone(), scratch, &mut out.bin2)?;
        self.count
            .read(&self.source, rows, scratch, &mut out.count)?;
        self.check_batch(out)
    }

    /// [`Cool::read_pixels`] through libhdf5's own reader for every column
    /// (serialized across threads): the reference the parallel path must
    /// match, and a way around it should a file defeat it.
    pub fn read_pixels_library(&self, rows: Range<u64>, out: &mut PixelBatch) -> Result<()> {
        self.check_rows(&rows)?;
        out.clear();
        self.bin1.read_library(rows.clone(), &mut out.bin1)?;
        self.bin2.read_library(rows.clone(), &mut out.bin2)?;
        self.count.read_library(rows, &mut out.count)?;
        self.check_batch(out)
    }

    fn check_rows(&self, rows: &Range<u64>) -> Result<()> {
        if rows.start > rows.end || rows.end > self.n_pixels {
            return Err(Error::OutOfRange {
                what: "pixel rows",
                offset: if rows.start > rows.end {
                    rows.start
                } else {
                    rows.end
                },
                lo: 0,
                hi: self.n_pixels,
            });
        }
        Ok(())
    }

    fn check_batch(&self, out: &PixelBatch) -> Result<()> {
        let n = out.count.len();
        if out.bin1.len() != n || out.bin2.len() != n {
            return Err(Error::invalid(
                "cool pixels",
                "columns decoded to different lengths",
            ));
        }
        if let Some(&b) = out
            .bin1
            .iter()
            .chain(&out.bin2)
            .find(|&&b| b >= self.n_bins)
        {
            return Err(Error::invalid(
                "cool pixels",
                format!("bin id {b} out of range (n_bins {})", self.n_bins),
            ));
        }
        Ok(())
    }
}

fn gcd(mut a: u64, mut b: u64) -> u64 {
    while b != 0 {
        (a, b) = (b, a % b);
    }
    a
}

fn check_format(grp: &hdf5::Group) -> Result<()> {
    match h5::attr_string(grp, "format")?.as_deref() {
        None | Some("HDF5::Cooler") => {}
        Some("HDF5::MCOOL") => {
            return Err(Error::invalid(
                "cool format",
                "this is an .mcool; open one of its /resolutions/<res> groups",
            ));
        }
        Some(other) => {
            return Err(Error::invalid(
                "cool format",
                format!("format attribute {other:?} is not HDF5::Cooler"),
            ));
        }
    }
    let version = h5::attr_int(grp, "format-version")?
        .ok_or_else(|| Error::invalid("cool format", "missing format-version attribute"))?;
    if !FORMAT_VERSIONS.contains(&version) {
        return Err(Error::Unsupported(
            "cooler format version other than 2 or 3",
        ));
    }
    Ok(())
}

fn bin_size(grp: &hdf5::Group) -> Result<u32> {
    if h5::attr_string(grp, "bin-type")?.is_some_and(|t| t != "fixed") {
        return Err(Error::Unsupported("cooler bins of variable width"));
    }
    let size = h5::attr_int(grp, "bin-size")?
        .ok_or(Error::Unsupported("cooler bins of variable width"))?;
    u32::try_from(size)
        .ok()
        .filter(|&s| s > 0)
        .ok_or_else(|| Error::invalid("cool bin-size", format!("bin size {size} out of range")))
}

fn read_chroms(grp: &hdf5::Group) -> Result<Vec<CoolChrom>> {
    let names = h5::read_strings(&h5::dataset(grp, "chroms/name")?, "chroms/name")?;
    let lengths = h5::read_ints(&h5::dataset(grp, "chroms/length")?, "chroms/length")?;
    if names.len() != lengths.len() {
        return Err(Error::invalid(
            "cool chroms",
            format!("{} names, {} lengths", names.len(), lengths.len()),
        ));
    }
    names
        .into_iter()
        .zip(lengths)
        .map(|(name, len)| {
            let length = u32::try_from(len).map_err(|_| {
                Error::invalid(
                    "chroms/length",
                    format!("{name}: length {len} out of range"),
                )
            })?;
            Ok(CoolChrom { name, length })
        })
        .collect()
}

/// Reads the bin table and checks it is the fixed-width tiling of every
/// chromosome that `indexes/chrom_offset` describes. Returns the offsets
/// and the bin count.
fn read_bins(grp: &hdf5::Group, chroms: &[CoolChrom], bin_size: u32) -> Result<(Vec<u64>, u64)> {
    let chrom_ds = h5::dataset(grp, "bins/chrom")?;
    let chrom_col = h5::read_ints(&chrom_ds, "bins/chrom")?;
    let start = h5::read_ints(&h5::dataset(grp, "bins/start")?, "bins/start")?;
    let end = h5::read_ints(&h5::dataset(grp, "bins/end")?, "bins/end")?;
    let n_bins = chrom_col.len();
    if start.len() != n_bins || end.len() != n_bins {
        return Err(Error::invalid(
            "cool bins",
            format!(
                "column lengths differ: chrom {n_bins}, start {}, end {}",
                start.len(),
                end.len()
            ),
        ));
    }

    // Stored chromosome values -> chromosome index. An enum column maps
    // through its member names; a plain integer column holds the index.
    let by_value: Option<HashMap<i64, usize>> = match h5::enum_members(&chrom_ds, "bins/chrom")? {
        None => None,
        Some(members) => {
            let by_name: HashMap<&str, usize> = chroms
                .iter()
                .enumerate()
                .map(|(i, c)| (c.name.as_str(), i))
                .collect();
            let mut map = HashMap::new();
            for (name, value) in members {
                if let Some(&i) = by_name.get(name.as_str()) {
                    map.insert(value, i);
                }
            }
            Some(map)
        }
    };

    let offsets = h5::read_ints(
        &h5::dataset(grp, "indexes/chrom_offset")?,
        "indexes/chrom_offset",
    )?;
    if offsets.len() != chroms.len() + 1 {
        return Err(Error::invalid(
            "indexes/chrom_offset",
            format!("{} entries for {} chromosomes", offsets.len(), chroms.len()),
        ));
    }
    let offsets: Vec<u64> = offsets
        .iter()
        .map(|&o| u64::try_from(o))
        .collect::<std::result::Result<_, _>>()
        .map_err(|_| Error::invalid("indexes/chrom_offset", "negative offset"))?;
    if offsets[0] != 0 || offsets[chroms.len()] != n_bins as u64 {
        return Err(Error::invalid(
            "indexes/chrom_offset",
            format!(
                "offsets span {}..{}, bins table has {n_bins} rows",
                offsets[0],
                offsets[chroms.len()]
            ),
        ));
    }
    // Check the whole index before any row access: with the endpoints
    // above, monotonic offsets keep every chromosome inside the bin table.
    if offsets.windows(2).any(|pair| pair[0] > pair[1]) {
        return Err(Error::invalid(
            "indexes/chrom_offset",
            "chromosome offsets are not monotonic",
        ));
    }

    let size = u64::from(bin_size);
    for (c, chrom) in chroms.iter().enumerate() {
        let (lo, hi) = (offsets[c], offsets[c + 1]);
        let len = u64::from(chrom.length);
        if hi - lo != len.div_ceil(size) {
            return Err(Error::invalid(
                "cool bins",
                format!(
                    "{}: {} bins in chrom_offset, a {len} bp chromosome has {} bins of {size} bp",
                    chrom.name,
                    hi - lo,
                    len.div_ceil(size)
                ),
            ));
        }
        for i in lo..hi {
            let local = i - lo;
            let s = local * size;
            let e = (s + size).min(len);
            let row = i as usize;
            let stored = chrom_col[row];
            let chrom_index = match &by_value {
                Some(map) => map.get(&stored).copied(),
                None => usize::try_from(stored).ok(),
            };
            if chrom_index != Some(c) || start[row] != s as i64 || end[row] != e as i64 {
                return Err(Error::invalid(
                    "cool bins",
                    format!(
                        "bin {i} is (chrom {stored}, {}, {}), expected ({}, {s}, {e}) \
                         for fixed {size} bp bins",
                        start[row], end[row], chrom.name
                    ),
                ));
            }
        }
    }
    Ok((offsets, n_bins as u64))
}
