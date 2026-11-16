//! A minimal cooler writer: HDF5 files with cooler's layout and the storage
//! options cooler (through h5py) and other writers use.

use std::path::Path;

use hdf5_metno as hdf5;
use hdf5_metno_sys as sys;

use hdf5::types::{EnumMember, EnumType, FixedAscii, IntSize, TypeDescriptor, VarLenUnicode};
use hdf5::{Group, H5Type};

/// How one dataset is stored.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum Storage {
    /// Resizable and chunked, as h5py writes cooler's columns.
    Chunked {
        chunk: usize,
        shuffle: bool,
        gzip: Option<u8>,
        fletcher: bool,
    },
    /// Fixed-size, one contiguous extent.
    Contiguous,
    /// Chunked with the lossless integer scale-offset filter, which the
    /// parallel path does not decode.
    ScaleOffset { chunk: usize },
}

impl Storage {
    /// h5py's `compression="gzip", shuffle=True`.
    pub const fn gzip_shuffle(chunk: usize) -> Storage {
        Storage::Chunked {
            chunk,
            shuffle: true,
            gzip: Some(6),
            fletcher: false,
        }
    }

    pub const fn plain(chunk: usize) -> Storage {
        Storage::Chunked {
            chunk,
            shuffle: false,
            gzip: None,
            fletcher: false,
        }
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum CountType {
    I32,
    F64,
}

/// One map's content and how to store it.
#[derive(Debug, Clone)]
pub struct MapSpec {
    pub chroms: Vec<(String, u32)>,
    pub bin_size: u32,
    /// Sorted by `(bin1, bin2)`.
    pub pixels: Vec<(u64, u64, f64)>,
    pub count_type: CountType,
    /// `bin1_id`, `bin2_id`, `count`.
    pub storage: [Storage; 3],
    /// Store `bins/chrom` as an enum of chromosome names (cooler's layout)
    /// rather than plain indices.
    pub chrom_enum: bool,
    /// Fixed-length (`|S64`) chromosome names rather than variable-length.
    pub fixed_names: bool,
    pub storage_mode: Option<&'static str>,
    pub format_version: i64,
}

impl MapSpec {
    pub fn new(chroms: &[(&str, u32)], bin_size: u32) -> MapSpec {
        MapSpec {
            chroms: chroms.iter().map(|&(n, l)| (n.to_owned(), l)).collect(),
            bin_size,
            pixels: Vec::new(),
            count_type: CountType::I32,
            storage: [Storage::gzip_shuffle(64); 3],
            chrom_enum: true,
            fixed_names: true,
            storage_mode: None,
            format_version: 3,
        }
    }

    /// Bins per chromosome and in total.
    pub fn n_bins(&self) -> u64 {
        self.chroms
            .iter()
            .map(|&(_, l)| u64::from(l).div_ceil(u64::from(self.bin_size)))
            .sum()
    }

    pub fn chrom_offsets(&self) -> Vec<u64> {
        let mut out = vec![0];
        for &(_, l) in &self.chroms {
            let n = u64::from(l).div_ceil(u64::from(self.bin_size));
            out.push(out.last().unwrap() + n);
        }
        out
    }
}

fn str_attr(loc: &Group, name: &str, value: &str) -> hdf5::Result<()> {
    let v: VarLenUnicode = value.parse().unwrap();
    loc.new_attr::<VarLenUnicode>()
        .create(name)?
        .write_scalar(&v)
}

fn int_attr(loc: &Group, name: &str, value: i64) -> hdf5::Result<()> {
    loc.new_attr::<i64>().create(name)?.write_scalar(&value)
}

/// Writes a 1-D column with the given storage.
fn write_col<T: H5Type>(
    g: &Group,
    name: &str,
    data: &[T],
    storage: Storage,
) -> hdf5::Result<hdf5::Dataset> {
    let n = data.len();
    let b = g.new_dataset_builder().empty::<T>();
    let ds = match storage {
        Storage::Contiguous => b.shape(n).create(name)?,
        Storage::Chunked {
            chunk,
            shuffle,
            gzip,
            fletcher,
        } => {
            // h5py's filter order: shuffle, deflate, fletcher32.
            let mut b = b.shape(n..).chunk(chunk);
            if shuffle {
                b = b.shuffle();
            }
            if let Some(level) = gzip {
                b = b.deflate(level);
            }
            if fletcher {
                b = b.fletcher32();
            }
            b.create(name)?
        }
        Storage::ScaleOffset { chunk } => b
            .shape(n..)
            .chunk(chunk)
            .scale_offset(hdf5::filters::ScaleOffset::Integer(0))
            .create(name)?,
    };
    if n > 0 {
        ds.write_raw(data)?;
    }
    Ok(ds)
}

/// Writes an enum-typed `bins/chrom` column holding chromosome indices.
fn write_enum_chrom(g: &Group, names: &[String], data: &[i32]) -> hdf5::Result<()> {
    let desc = TypeDescriptor::Enum(EnumType {
        size: IntSize::U4,
        signed: true,
        members: names
            .iter()
            .enumerate()
            .map(|(i, n)| EnumMember {
                name: n.clone(),
                value: i as u64,
            })
            .collect(),
    });
    let ds = g
        .new_dataset_builder()
        .empty_as(&desc)
        .shape(data.len()..)
        .chunk(data.len().clamp(1, 1000))
        .deflate(6)
        .create("chrom")?;
    if !data.is_empty() {
        let dtype = ds.dtype()?;
        let _lock = sys::LOCK.lock();
        // SAFETY: the buffer holds `data.len()` 4-byte values, the dataset's
        // own enum type (whose base is a native 32-bit integer).
        let status = unsafe {
            sys::h5d::H5Dwrite(
                ds.id(),
                dtype.id(),
                sys::h5s::H5S_ALL,
                sys::h5s::H5S_ALL,
                sys::h5p::H5P_DEFAULT,
                data.as_ptr().cast(),
            )
        };
        assert!(status >= 0, "H5Dwrite of enum column failed");
    }
    Ok(())
}

/// Writes one cooler map into `g`.
fn write_map(g: &Group, spec: &MapSpec) -> hdf5::Result<()> {
    str_attr(g, "format", "HDF5::Cooler")?;
    int_attr(g, "format-version", spec.format_version)?;
    str_attr(g, "bin-type", "fixed")?;
    int_attr(g, "bin-size", i64::from(spec.bin_size))?;
    str_attr(g, "genome-assembly", "test")?;
    if let Some(mode) = spec.storage_mode {
        str_attr(g, "storage-mode", mode)?;
    }

    let names: Vec<String> = spec.chroms.iter().map(|(n, _)| n.clone()).collect();
    let chroms = g.create_group("chroms")?;
    if spec.fixed_names {
        let fixed: Vec<FixedAscii<64>> = names
            .iter()
            .map(|n| FixedAscii::from_ascii(n.as_bytes()).unwrap())
            .collect();
        write_col(&chroms, "name", &fixed, Storage::Contiguous)?;
    } else {
        let var: Vec<VarLenUnicode> = names.iter().map(|n| n.parse().unwrap()).collect();
        write_col(&chroms, "name", &var, Storage::Contiguous)?;
    }
    let lengths: Vec<i32> = spec.chroms.iter().map(|&(_, l)| l as i32).collect();
    write_col(&chroms, "length", &lengths, Storage::Contiguous)?;

    let (mut chrom, mut start, mut end) = (Vec::new(), Vec::new(), Vec::new());
    for (c, &(_, len)) in spec.chroms.iter().enumerate() {
        let mut s = 0u32;
        while s < len {
            chrom.push(c as i32);
            start.push(s as i32);
            end.push(s.saturating_add(spec.bin_size).min(len) as i32);
            s = s.saturating_add(spec.bin_size);
        }
    }
    let bins = g.create_group("bins")?;
    let bin_storage = Storage::gzip_shuffle(chrom.len().clamp(1, 1000));
    if spec.chrom_enum {
        write_enum_chrom(&bins, &names, &chrom)?;
    } else {
        write_col(&bins, "chrom", &chrom, bin_storage)?;
    }
    write_col(&bins, "start", &start, bin_storage)?;
    write_col(&bins, "end", &end, bin_storage)?;

    let pixels = g.create_group("pixels")?;
    let bin1: Vec<i64> = spec.pixels.iter().map(|p| p.0 as i64).collect();
    let bin2: Vec<i64> = spec.pixels.iter().map(|p| p.1 as i64).collect();
    write_col(&pixels, "bin1_id", &bin1, spec.storage[0])?;
    write_col(&pixels, "bin2_id", &bin2, spec.storage[1])?;
    match spec.count_type {
        CountType::I32 => {
            let c: Vec<i32> = spec.pixels.iter().map(|p| p.2 as i32).collect();
            write_col(&pixels, "count", &c, spec.storage[2])?;
        }
        CountType::F64 => {
            let c: Vec<f64> = spec.pixels.iter().map(|p| p.2).collect();
            write_col(&pixels, "count", &c, spec.storage[2])?;
        }
    }

    let indexes = g.create_group("indexes")?;
    let offsets: Vec<i64> = spec.chrom_offsets().iter().map(|&o| o as i64).collect();
    write_col(&indexes, "chrom_offset", &offsets, Storage::Contiguous)?;
    let n_bins = spec.n_bins() as usize;
    let mut bin1_offset = vec![0i64; n_bins + 1];
    for p in &spec.pixels {
        bin1_offset[p.0 as usize + 1] += 1;
    }
    for i in 0..n_bins {
        bin1_offset[i + 1] += bin1_offset[i];
    }
    write_col(
        &indexes,
        "bin1_offset",
        &bin1_offset,
        Storage::gzip_shuffle(1000.min(n_bins + 1)),
    )?;
    Ok(())
}

/// Writes a single-resolution `.cool`.
pub fn write_cool(path: &Path, spec: &MapSpec) {
    let file = hdf5::File::create(path).unwrap();
    write_map(&file, spec).unwrap();
}

/// Writes a `.mcool` holding each map at `/resolutions/<bin size>`.
pub fn write_mcool(path: &Path, specs: &[MapSpec]) {
    let file = hdf5::File::create(path).unwrap();
    str_attr(&file, "format", "HDF5::MCOOL").unwrap();
    int_attr(&file, "format-version", 2).unwrap();
    let res = file.create_group("resolutions").unwrap();
    for spec in specs {
        let g = res.create_group(&spec.bin_size.to_string()).unwrap();
        write_map(&g, spec).unwrap();
    }
}
