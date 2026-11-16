//! Small cooler inputs for contact conversion, written with the HDF5 library.
#![allow(dead_code)]

use std::path::Path;

use hdf5_metno::types::{EnumMember, EnumType, FixedAscii, IntSize, TypeDescriptor, VarLenUnicode};
use hdf5_metno::{Group, H5Type};
use hdf5_metno_sys as sys;

#[derive(Clone)]
pub struct MapSpec {
    pub chroms: Vec<(String, u32)>,
    pub bin_size: u32,
    pub pixels: Vec<(u64, u64, f64)>,
}

impl MapSpec {
    pub fn new(chroms: &[(&str, u32)], bin_size: u32) -> Self {
        Self {
            chroms: chroms
                .iter()
                .map(|&(name, len)| (name.to_owned(), len))
                .collect(),
            bin_size,
            pixels: Vec::new(),
        }
    }

    pub fn chrom_offsets(&self) -> Vec<u64> {
        let mut offsets = vec![0];
        for &(_, len) in &self.chroms {
            offsets.push(offsets.last().unwrap() + u64::from(len).div_ceil(self.bin_size.into()));
        }
        offsets
    }

    pub fn n_bins(&self) -> u64 {
        *self.chrom_offsets().last().unwrap()
    }
}

fn text_attr(group: &Group, name: &str, value: &str) {
    group
        .new_attr::<VarLenUnicode>()
        .create(name)
        .unwrap()
        .write_scalar(&value.parse::<VarLenUnicode>().unwrap())
        .unwrap();
}

fn int_attr(group: &Group, name: &str, value: i64) {
    group
        .new_attr::<i64>()
        .create(name)
        .unwrap()
        .write_scalar(&value)
        .unwrap();
}

fn column<T: H5Type>(group: &Group, name: &str, values: &[T], compressed: bool) {
    let builder = group.new_dataset_builder().empty::<T>();
    let dataset = if compressed {
        builder
            .shape(values.len()..)
            .chunk(values.len().clamp(1, 1000))
            .shuffle()
            .deflate(6)
            .create(name)
            .unwrap()
    } else {
        builder.shape(values.len()).create(name).unwrap()
    };
    if !values.is_empty() {
        dataset.write_raw(values).unwrap();
    }
}

/// Writes an enum-typed `bins/chrom` column holding chromosome indices.
fn write_enum_chrom(g: &Group, names: &[String], data: &[i32]) -> hdf5_metno::Result<()> {
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

fn write_map(group: &Group, spec: &MapSpec) {
    text_attr(group, "format", "HDF5::Cooler");
    text_attr(group, "bin-type", "fixed");
    text_attr(group, "genome-assembly", "test");
    int_attr(group, "format-version", 3);
    int_attr(group, "bin-size", spec.bin_size.into());

    let chroms = group.create_group("chroms").unwrap();
    let names: Vec<FixedAscii<64>> = spec
        .chroms
        .iter()
        .map(|(name, _)| FixedAscii::from_ascii(name.as_bytes()).unwrap())
        .collect();
    column(&chroms, "name", &names, false);
    column(
        &chroms,
        "length",
        &spec
            .chroms
            .iter()
            .map(|&(_, len)| len as i32)
            .collect::<Vec<_>>(),
        false,
    );

    let (mut ids, mut starts, mut ends) = (Vec::new(), Vec::new(), Vec::new());
    for (id, &(_, len)) in spec.chroms.iter().enumerate() {
        for start in (0..len).step_by(spec.bin_size as usize) {
            ids.push(id as i32);
            starts.push(start as i32);
            ends.push(start.saturating_add(spec.bin_size).min(len) as i32);
        }
    }
    let bins = group.create_group("bins").unwrap();
    let chrom_names: Vec<String> = spec.chroms.iter().map(|(name, _)| name.clone()).collect();
    write_enum_chrom(&bins, &chrom_names, &ids).unwrap();
    column(&bins, "start", &starts, true);
    column(&bins, "end", &ends, true);

    let pixels = group.create_group("pixels").unwrap();
    column(
        &pixels,
        "bin1_id",
        &spec.pixels.iter().map(|p| p.0 as i64).collect::<Vec<_>>(),
        true,
    );
    column(
        &pixels,
        "bin2_id",
        &spec.pixels.iter().map(|p| p.1 as i64).collect::<Vec<_>>(),
        true,
    );
    column(
        &pixels,
        "count",
        &spec.pixels.iter().map(|p| p.2 as i32).collect::<Vec<_>>(),
        true,
    );

    let indexes = group.create_group("indexes").unwrap();
    column(
        &indexes,
        "chrom_offset",
        &spec
            .chrom_offsets()
            .iter()
            .map(|&n| n as i64)
            .collect::<Vec<_>>(),
        false,
    );
    let mut offsets = vec![0i64; spec.n_bins() as usize + 1];
    for &(bin1, _, _) in &spec.pixels {
        offsets[bin1 as usize + 1] += 1;
    }
    for i in 1..offsets.len() {
        offsets[i] += offsets[i - 1];
    }
    column(&indexes, "bin1_offset", &offsets, true);
}

pub fn write_cool(path: &Path, spec: &MapSpec) {
    let file = hdf5_metno::File::create(path).unwrap();
    write_map(&file, spec);
}

pub fn write_mcool(path: &Path, specs: &[MapSpec]) {
    let file = hdf5_metno::File::create(path).unwrap();
    text_attr(&file, "format", "HDF5::MCOOL");
    int_attr(&file, "format-version", 2);
    let resolutions = file.create_group("resolutions").unwrap();
    for spec in specs {
        write_map(
            &resolutions
                .create_group(&spec.bin_size.to_string())
                .unwrap(),
            spec,
        );
    }
}

/// Sorted, unique pixels from `(bin1, bin2, count)` triples; with
/// `upper`, each pair is ordered so `bin1 <= bin2`.
pub fn normalize_pixels(raw: &[(u64, u64, f64)], upper: bool) -> Vec<(u64, u64, f64)> {
    let mut v: Vec<(u64, u64, f64)> = raw
        .iter()
        .map(|&(a, b, c)| if upper && a > b { (b, a, c) } else { (a, b, c) })
        .collect();
    v.sort_by_key(|x| (x.0, x.1));
    v.dedup_by(|x, y| x.0 == y.0 && x.1 == y.1);
    v
}

/// Dense-ish pseudo-random pixels over `n_bins` bins.
pub fn sample_pixels(n_bins: u64, n: usize, seed: u64, integer: bool) -> Vec<(u64, u64, f64)> {
    let mut state = seed
        .wrapping_mul(6364136223846793005)
        .wrapping_add(1442695040888963407);
    let mut next = || {
        state = state
            .wrapping_mul(6364136223846793005)
            .wrapping_add(1442695040888963407);
        state >> 33
    };
    let raw: Vec<(u64, u64, f64)> = (0..n)
        .map(|_| {
            let a = next() % n_bins;
            let b = next() % n_bins;
            let c = if integer {
                (next() % 1000 + 1) as f64
            } else {
                (next() % 100_000) as f64 / 7.0
            };
            (a, b, c)
        })
        .collect();
    normalize_pixels(&raw, true)
}
