use crate::support::{
    CHROMS, CountType, MapSpec, Storage, check_all, expected, read, sample_pixels, spec,
    write_and_open, write_cool,
};
use omicsweave_cool::{Cool, PixelBatch, PixelScratch, ReadPath};

#[test]
fn fletcher32_chunks() {
    let dir = tempfile::tempdir().unwrap();
    let mut s = spec(2000, 3);
    s.storage = [Storage::Chunked {
        chunk: 100,
        shuffle: true,
        gzip: Some(4),
        fletcher: true,
    }; 3];
    s.storage[2] = Storage::Chunked {
        chunk: 77,
        shuffle: false,
        gzip: None,
        fletcher: true,
    };
    let cool = write_and_open(dir.path(), &s);
    assert_eq!(cool.read_paths(), [ReadPath::Chunked; 3]);
    check_all(&cool, &s, 250);
}

#[test]
fn checksum_inside_deflate_stream_matches_library() {
    let dir = tempfile::tempdir().unwrap();
    let path = dir.path().join("inner-checksum.cool");
    let s = spec(1000, 95);
    write_cool(&path, &s);
    {
        let file = hdf5_metno::File::open_rw(&path).unwrap();
        file.unlink("pixels/count").unwrap();
        let ds = file
            .new_dataset::<i32>()
            .shape(s.pixels.len()..)
            .chunk(64)
            .fletcher32()
            .deflate(6)
            .create("pixels/count")
            .unwrap();
        let counts: Vec<i32> = s.pixels.iter().map(|p| p.2 as i32).collect();
        ds.write_raw(&counts).unwrap();
    }
    let cool = Cool::open(&path, "/").unwrap();
    assert_eq!(cool.read_paths(), [ReadPath::Chunked; 3]);
    check_all(&cool, &s, 91);
}

#[test]
fn shuffle_only_partial_chunks_and_checksum_before_shuffle_match_library() {
    let dir = tempfile::tempdir().unwrap();
    let path = dir.path().join("shuffle-planes.cool");
    let mut s = spec(1000, 195);
    s.storage = [Storage::Chunked {
        chunk: 73,
        shuffle: true,
        gzip: None,
        fletcher: false,
    }; 3];
    write_cool(&path, &s);
    let cool = Cool::open(&path, "/").unwrap();
    assert_eq!(cool.read_paths(), [ReadPath::Chunked; 3]);
    let n = s.pixels.len() as u64;
    for rows in [0..1, 17..19, 72..75, n - 2..n, 0..n] {
        let fast = read(&cool, rows.clone());
        let mut oracle = PixelBatch::new();
        cool.read_pixels_library(rows.clone(), &mut oracle).unwrap();
        assert_eq!(fast, oracle);
        assert_eq!(fast, expected(&s, rows));
    }
    drop(cool);
    {
        let file = hdf5_metno::File::open_rw(&path).unwrap();
        file.unlink("pixels/count").unwrap();
        let ds = file
            .new_dataset::<i32>()
            .shape(s.pixels.len()..)
            .chunk(73)
            .fletcher32()
            .shuffle()
            .deflate(6)
            .create("pixels/count")
            .unwrap();
        ds.write_raw(&s.pixels.iter().map(|p| p.2 as i32).collect::<Vec<_>>())
            .unwrap();
    }
    let cool = Cool::open(&path, "/").unwrap();
    check_all(&cool, &s, 91);
}

#[test]
fn failed_chunk_decode_allows_scratch_reuse() {
    use std::io::{Seek, SeekFrom, Write};

    let dir = tempfile::tempdir().unwrap();
    let path = dir.path().join("corrupt.cool");
    let s = spec(1000, 96);
    write_cool(&path, &s);
    let addr = {
        let file = hdf5_metno::File::open(&path).unwrap();
        file.dataset("pixels/count")
            .unwrap()
            .chunk_info(0)
            .unwrap()
            .addr
    };
    {
        let mut file = std::fs::OpenOptions::new().write(true).open(&path).unwrap();
        file.seek(SeekFrom::Start(addr)).unwrap();
        // A zero zlib CMF byte is invalid, independent of payload contents.
        file.write_all(&[0]).unwrap();
    }
    let cool = Cool::open(&path, "/").unwrap();
    let mut scratch = PixelScratch::new();
    let mut batch = PixelBatch::new();
    assert!(cool.read_pixels(0..1, &mut scratch, &mut batch).is_err());
    cool.read_pixels(64..80, &mut scratch, &mut batch).unwrap();
    assert_eq!(batch, expected(&s, 64..80));
    let mut oracle = PixelBatch::new();
    cool.read_pixels_library(64..80, &mut oracle).unwrap();
    assert_eq!(batch, oracle);
}

#[test]
fn skipped_filters_partial_reads_match_library() {
    use hdf5_metno_sys as sys;

    let dir = tempfile::tempdir().unwrap();
    let path = dir.path().join("skipped.cool");
    let mut s = spec(1000, 92);
    s.storage = [Storage::gzip_shuffle(64); 3];
    write_cool(&path, &s);
    {
        let file = hdf5_metno::File::open_rw(&path).unwrap();
        let ds = file.dataset("pixels/bin1_id").unwrap();
        let raw: Vec<u8> = s.pixels[..64]
            .iter()
            .flat_map(|p| (p.0 as i64).to_ne_bytes())
            .collect();
        let offset = [0u64];
        let _lock = sys::LOCK.lock();
        // SAFETY: ds is a live 1-D i64 dataset with 64-row chunks;
        // offset is one aligned coordinate and raw holds exactly 64
        // native-endian i64 values. Mask 3 skips shuffle and deflate.
        let status = unsafe {
            sys::h5d::H5Dwrite_chunk(
                ds.id(),
                sys::h5p::H5P_DEFAULT,
                3,
                offset.as_ptr(),
                raw.len(),
                raw.as_ptr().cast(),
            )
        };
        assert!(status >= 0);
    }
    let cool = Cool::open(&path, "/").unwrap();
    assert_eq!(cool.read_paths(), [ReadPath::Chunked; 3]);
    for rows in [0..1, 17..19, 63..65, 1..80] {
        let fast = read(&cool, rows.clone());
        let mut oracle = PixelBatch::new();
        cool.read_pixels_library(rows.clone(), &mut oracle).unwrap();
        assert_eq!(fast, oracle);
        assert_eq!(fast, expected(&s, rows));
    }
}

#[test]
fn unwritten_chunks_use_library_fill_values() {
    let dir = tempfile::tempdir().unwrap();
    let path = dir.path().join("fill.cool");
    let s = spec(1000, 93);
    write_cool(&path, &s);
    {
        let file = hdf5_metno::File::open_rw(&path).unwrap();
        file.unlink("pixels/count").unwrap();
        file.new_dataset::<i32>()
            .shape(s.pixels.len()..)
            .chunk(64)
            .fill_value(7)
            .create("pixels/count")
            .unwrap();
    }
    let cool = Cool::open(&path, "/").unwrap();
    assert_eq!(cool.read_paths()[2], ReadPath::Library);
    let fast = read(&cool, 17..85);
    let mut oracle = PixelBatch::new();
    cool.read_pixels_library(17..85, &mut oracle).unwrap();
    assert_eq!(fast, oracle);
    assert!(fast.count.iter().all(|&c| c == 7.0));
}

#[test]
fn large_unwritten_pixel_extent_uses_fill_without_a_logical_chunk_directory() {
    let dir = tempfile::tempdir().unwrap();
    let path = dir.path().join("sparse-extent.cool");
    write_cool(&path, &MapSpec::new(&[("chr1", 10)], 1));
    const N: usize = 1 << 26;
    {
        let file = hdf5_metno::File::open_rw(&path).unwrap();
        for name in ["pixels/bin1_id", "pixels/bin2_id", "pixels/count"] {
            file.unlink(name).unwrap();
            file.new_dataset::<i32>()
                .shape(N..)
                .chunk(16)
                .fill_value(0)
                .create(name)
                .unwrap();
        }
    }
    assert!(std::fs::metadata(&path).unwrap().len() < 1 << 20);
    let cool = Cool::open(&path, "/").unwrap();
    assert_eq!(cool.n_pixels(), N as u64);
    assert_eq!(cool.read_paths(), [ReadPath::Library; 3]);
    let fast = read(&cool, N as u64 - 2..N as u64);
    assert_eq!(fast.bin1, [0, 0]);
    assert_eq!(fast.bin2, [0, 0]);
    assert_eq!(fast.count, [0.0, 0.0]);
}

#[test]
fn unsupported_filter_falls_back_to_libhdf5() {
    let dir = tempfile::tempdir().unwrap();
    let mut s = spec(2000, 7);
    s.storage[1] = Storage::ScaleOffset { chunk: 64 };
    s.storage[2] = Storage::ScaleOffset { chunk: 64 };
    let cool = write_and_open(dir.path(), &s);
    assert_eq!(
        cool.read_paths(),
        [ReadPath::Chunked, ReadPath::Library, ReadPath::Library]
    );
    check_all(&cool, &s, 200);
}

#[test]
fn float_counts() {
    let dir = tempfile::tempdir().unwrap();
    let mut s = MapSpec::new(CHROMS, 100);
    s.pixels = sample_pixels(s.n_bins(), 2000, 8, false);
    s.count_type = CountType::F64;
    let cool = write_and_open(dir.path(), &s);
    assert!(s.pixels.iter().any(|p| p.2.fract() != 0.0));
    check_all(&cool, &s, 256);
}

#[test]
fn float_layouts_preserve_ieee_fast_paths_and_convert_custom_bias() {
    use hdf5_metno as hdf5;
    use hdf5_metno_sys as sys;

    for (size, big_endian, custom_bias) in [
        (4, false, false),
        (4, true, false),
        (8, false, false),
        (8, true, false),
        (4, false, true),
        (4, true, true),
    ] {
        let dir = tempfile::tempdir().unwrap();
        let path = dir.path().join("float-layout.cool");
        let mut s = MapSpec::new(&[("chr1", 10)], 1);
        s.pixels = vec![(0, 1, 1.5), (2, 3, 2.5)];
        write_cool(&path, &s);
        let offset = {
            let file = hdf5::File::open_rw(&path).unwrap();
            file.unlink("pixels/count").unwrap();
            let dtype = {
                let _lock = sys::LOCK.lock();
                let ieee = match (size, big_endian) {
                    (4, false) => *sys::h5t::H5T_IEEE_F32LE,
                    (4, true) => *sys::h5t::H5T_IEEE_F32BE,
                    (8, false) => *sys::h5t::H5T_IEEE_F64LE,
                    (8, true) => *sys::h5t::H5T_IEEE_F64BE,
                    _ => unreachable!(),
                };
                // SAFETY: HDF5 is initialized by the open file, ieee is a
                // predefined datatype, and sys::LOCK serializes the copy.
                let id = unsafe { sys::h5t::H5Tcopy(ieee) };
                assert!(id >= 0);
                // SAFETY: id is a newly owned, valid datatype id; ownership
                // is transferred once to its RAII wrapper under sys::LOCK.
                let dtype: hdf5::Datatype = unsafe { hdf5::from_id(id) }.unwrap();
                if custom_bias {
                    // SAFETY: this mutable transient copy is owned by dtype,
                    // has not been used by a dataset, and sys::LOCK is held.
                    assert!(unsafe { sys::h5t::H5Tset_ebias(dtype.id(), 126) } >= 0);
                }
                dtype
            };
            let ds = file
                .new_dataset_builder()
                .empty_as(&dtype)
                .shape(2)
                .create("pixels/count")
                .unwrap();
            ds.write_raw(&[1.5f32, 2.5]).unwrap();
            file.flush().unwrap();
            ds.offset().unwrap() as usize
        };
        if custom_bias {
            let bytes = std::fs::read(&path).unwrap();
            let first: [u8; 4] = bytes[offset..offset + 4].try_into().unwrap();
            let raw = if big_endian {
                f32::from_be_bytes(first)
            } else {
                f32::from_le_bytes(first)
            };
            assert_ne!(raw, 1.5);
        }
        let cool = Cool::open(&path, "/").unwrap();
        assert_eq!(
            cool.read_paths()[2],
            if custom_bias {
                ReadPath::Library
            } else {
                ReadPath::Contiguous
            }
        );
        check_all(&cool, &s, 1);
    }
}
