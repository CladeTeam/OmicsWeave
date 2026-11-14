mod common;

use std::ops::Range;
use std::path::Path;

use common::{CountType, MapSpec, Storage, sample_pixels, write_cool, write_mcool};
use omicsweave_cool::{
    Cool, PixelBatch, PixelScratch, ReadPath, StorageMode, is_mcool, mcool_resolutions,
};
use omicsweave_core::Error;
use proptest::prelude::*;

const CHROMS: &[(&str, u32)] = &[("chr1", 10_000), ("chr2", 7_350), ("chrM", 16)];

fn spec(pixels: usize, seed: u64) -> MapSpec {
    let mut s = MapSpec::new(CHROMS, 100);
    s.pixels = sample_pixels(s.n_bins(), pixels, seed, true);
    s
}

fn expected(spec: &MapSpec, rows: Range<u64>) -> PixelBatch {
    let mut out = PixelBatch::new();
    for &(a, b, c) in &spec.pixels[rows.start as usize..rows.end as usize] {
        out.bin1.push(a);
        out.bin2.push(b);
        out.count.push(c);
    }
    out
}

fn read(cool: &Cool, rows: Range<u64>) -> PixelBatch {
    let mut out = PixelBatch::new();
    cool.read_pixels(rows, &mut PixelScratch::new(), &mut out)
        .unwrap();
    out
}

/// Reads every span with one scratch and checks the spans tile the
/// pixels, each span matching the written data and libhdf5's reader.
fn check_all(cool: &Cool, spec: &MapSpec, target: u64) {
    let spans = cool.pixel_spans(target);
    let mut scratch = PixelScratch::new();
    let mut batch = PixelBatch::new();
    let mut library = PixelBatch::new();
    let mut next = 0;
    for span in spans {
        assert_eq!(span.start, next);
        assert!(!span.is_empty());
        next = span.end;
        cool.read_pixels(span.clone(), &mut scratch, &mut batch)
            .unwrap();
        assert_eq!(batch, expected(spec, span.clone()));
        cool.read_pixels_library(span, &mut library).unwrap();
        assert_eq!(library, batch);
    }
    assert_eq!(next, spec.pixels.len() as u64);
}

fn write_and_open(dir: &Path, spec: &MapSpec) -> Cool {
    let path = dir.join("m.cool");
    write_cool(&path, spec);
    Cool::open(&path, "/").unwrap()
}

#[test]
fn metadata() {
    let dir = tempfile::tempdir().unwrap();
    let s = spec(500, 1);
    let cool = write_and_open(dir.path(), &s);
    let chroms: Vec<(&str, u32)> = cool
        .chroms()
        .iter()
        .map(|c| (c.name.as_str(), c.length))
        .collect();
    assert_eq!(chroms, CHROMS);
    assert_eq!(cool.bin_size(), 100);
    assert_eq!(cool.n_bins(), 100 + 74 + 1);
    assert_eq!(cool.chrom_offsets(), &[0, 100, 174, 175]);
    assert_eq!(cool.n_pixels(), s.pixels.len() as u64);
    assert_eq!(cool.storage_mode(), StorageMode::SymmetricUpper);
    assert_eq!(cool.assembly(), Some("test"));
    assert_eq!(cool.locate(0), (0, 0));
    assert_eq!(cool.locate(99), (0, 99));
    assert_eq!(cool.locate(100), (1, 0));
    assert_eq!(cool.locate(173), (1, 73));
    assert_eq!(cool.locate(174), (2, 0));
}

#[test]
fn locate_skips_empty_chromosomes() {
    let dir = tempfile::tempdir().unwrap();
    let mut s = MapSpec::new(&[("a", 250), ("empty", 0), ("b", 100)], 100);
    s.pixels = sample_pixels(s.n_bins(), 20, 3, true);
    let cool = write_and_open(dir.path(), &s);
    assert_eq!(cool.chrom_offsets(), &[0, 3, 3, 4]);
    assert_eq!(cool.locate(2), (0, 2));
    assert_eq!(cool.locate(3), (2, 0));
}

#[test]
fn gzip_shuffle_chunks() {
    let dir = tempfile::tempdir().unwrap();
    let s = spec(3000, 2);
    let cool = write_and_open(dir.path(), &s);
    assert_eq!(cool.read_paths(), [ReadPath::Chunked; 3]);
    for target in [1, 64, 100, 1000, 1 << 20] {
        check_all(&cool, &s, target);
    }
    // Ranges that start and end inside chunks.
    for rows in [0..1, 5..70, 63..65, 100..1999, 0..s.pixels.len() as u64] {
        assert_eq!(read(&cool, rows.clone()), expected(&s, rows));
    }
}

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
fn uncompressed_chunks() {
    let dir = tempfile::tempdir().unwrap();
    let mut s = spec(2000, 4);
    s.storage = [Storage::plain(128); 3];
    let cool = write_and_open(dir.path(), &s);
    assert_eq!(cool.read_paths(), [ReadPath::Chunked; 3]);
    check_all(&cool, &s, 300);
    assert_eq!(read(&cool, 127..129), expected(&s, 127..129));
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
fn largest_span_target_covers_all_pixels() {
    let dir = tempfile::tempdir().unwrap();
    let s = spec(500, 94);
    let cool = write_and_open(dir.path(), &s);
    assert_eq!(cool.pixel_spans(u64::MAX), vec![0..cool.n_pixels()]);
}

#[test]
fn contiguous_columns() {
    let dir = tempfile::tempdir().unwrap();
    let mut s = spec(2000, 5);
    s.storage = [Storage::Contiguous; 3];
    let cool = write_and_open(dir.path(), &s);
    assert_eq!(cool.read_paths(), [ReadPath::Contiguous; 3]);
    // Without chunks there are no boundaries to respect.
    let spans = cool.pixel_spans(300);
    assert_eq!(spans[0], 0..300);
    check_all(&cool, &s, 300);
}

#[test]
fn chunk_sizes_differ_between_columns() {
    let dir = tempfile::tempdir().unwrap();
    let mut s = spec(3000, 6);
    // As h5py chunks cooler's columns by bytes: the 4-byte counts get
    // twice the rows per chunk of the 8-byte ids.
    s.storage = [
        Storage::gzip_shuffle(50),
        Storage::gzip_shuffle(50),
        Storage::gzip_shuffle(100),
    ];
    let cool = write_and_open(dir.path(), &s);
    assert_eq!(cool.read_paths(), [ReadPath::Chunked; 3]);
    for span in cool.pixel_spans(120) {
        assert!(span.start % 100 == 0, "{span:?} splits a count chunk");
    }
    check_all(&cool, &s, 120);

    // Unrelated sizes: spans follow the largest chunk.
    s.storage = [
        Storage::gzip_shuffle(37),
        Storage::plain(53),
        Storage::gzip_shuffle(91),
    ];
    let path = dir.path().join("unrelated.cool");
    write_cool(&path, &s);
    let cool = Cool::open(&path, "/").unwrap();
    for span in cool.pixel_spans(200) {
        assert!(span.start % 91 == 0, "{span:?} splits a count chunk");
    }
    check_all(&cool, &s, 200);
    check_all(&cool, &s, 1);
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
fn plain_integer_bins_and_variable_length_names() {
    let dir = tempfile::tempdir().unwrap();
    let mut s = spec(500, 9);
    s.chrom_enum = false;
    s.fixed_names = false;
    s.format_version = 2;
    let cool = write_and_open(dir.path(), &s);
    assert_eq!(cool.chroms()[1].name, "chr2");
    check_all(&cool, &s, 100);
}

#[test]
fn square_storage_mode() {
    let dir = tempfile::tempdir().unwrap();
    let mut s = MapSpec::new(CHROMS, 100);
    let mut raw = sample_pixels(s.n_bins(), 800, 10, true);
    // Lower-triangle pixels exist in square maps.
    raw.extend(raw.clone().iter().map(|&(a, b, c)| (b, a, c + 1.0)));
    s.pixels = common::normalize_pixels(&raw, false);
    s.storage_mode = Some("square");
    let cool = write_and_open(dir.path(), &s);
    assert_eq!(cool.storage_mode(), StorageMode::Square);
    assert!(s.pixels.iter().any(|p| p.0 > p.1));
    check_all(&cool, &s, 100);
}

#[test]
fn mcool_resolution_groups() {
    let dir = tempfile::tempdir().unwrap();
    let path = dir.path().join("m.mcool");
    let mut fine = MapSpec::new(CHROMS, 100);
    fine.pixels = sample_pixels(fine.n_bins(), 1500, 11, true);
    let mut coarse = MapSpec::new(CHROMS, 1000);
    coarse.pixels = sample_pixels(coarse.n_bins(), 100, 12, true);
    // Written coarse first: resolutions come back sorted regardless.
    write_mcool(&path, &[coarse.clone(), fine.clone()]);

    assert!(is_mcool(&path).unwrap());
    assert_eq!(mcool_resolutions(&path).unwrap(), vec![100, 1000]);
    for s in [&fine, &coarse] {
        let cool = Cool::open(&path, &format!("/resolutions/{}", s.bin_size)).unwrap();
        assert_eq!(cool.bin_size(), s.bin_size);
        assert_eq!(cool.n_bins(), s.n_bins());
        check_all(&cool, s, 128);
    }
    // The root of an .mcool is not a map.
    assert!(matches!(Cool::open(&path, "/"), Err(Error::Invalid { .. })));
    assert!(matches!(
        Cool::open(&path, "/resolutions/5"),
        Err(Error::Invalid { .. })
    ));

    let single = dir.path().join("s.cool");
    write_cool(&single, &fine);
    assert!(!is_mcool(&single).unwrap());
    assert!(mcool_resolutions(&single).is_err());
}

#[test]
fn empty_pixels() {
    let dir = tempfile::tempdir().unwrap();
    let s = MapSpec::new(CHROMS, 100);
    let cool = write_and_open(dir.path(), &s);
    assert_eq!(cool.n_pixels(), 0);
    assert!(cool.pixel_spans(1000).is_empty());
    let mut out = PixelBatch {
        bin1: vec![1],
        bin2: vec![2],
        count: vec![3.0],
    };
    cool.read_pixels(0..0, &mut PixelScratch::new(), &mut out)
        .unwrap();
    assert!(out.is_empty() && out.bin1.is_empty() && out.bin2.is_empty());
}

#[test]
fn rows_out_of_range() {
    let dir = tempfile::tempdir().unwrap();
    let s = spec(100, 13);
    let cool = write_and_open(dir.path(), &s);
    let n = cool.n_pixels();
    let mut out = PixelBatch::new();
    let mut scratch = PixelScratch::new();
    assert!(matches!(
        cool.read_pixels(n - 1..n + 1, &mut scratch, &mut out),
        Err(Error::OutOfRange { .. })
    ));
    #[allow(clippy::reversed_empty_ranges)]
    let backwards = 5..2;
    assert!(matches!(
        cool.read_pixels(backwards, &mut scratch, &mut out),
        Err(Error::OutOfRange { .. })
    ));
}

#[test]
fn concurrent_reads_match_serial() {
    let dir = tempfile::tempdir().unwrap();
    let s = spec(20_000, 14);
    let cool = write_and_open(dir.path(), &s);
    let spans = cool.pixel_spans(700);
    assert!(spans.len() > 8);
    let next = std::sync::atomic::AtomicUsize::new(0);
    let results = std::sync::Mutex::new(vec![PixelBatch::new(); spans.len()]);
    std::thread::scope(|scope| {
        for _ in 0..8 {
            scope.spawn(|| {
                let mut scratch = PixelScratch::new();
                loop {
                    let i = next.fetch_add(1, std::sync::atomic::Ordering::Relaxed);
                    let Some(span) = spans.get(i) else { break };
                    let mut batch = PixelBatch::new();
                    cool.read_pixels(span.clone(), &mut scratch, &mut batch)
                        .unwrap();
                    results.lock().unwrap()[i] = batch;
                }
            });
        }
    });
    for (span, batch) in spans.iter().zip(results.into_inner().unwrap()) {
        assert_eq!(batch, expected(&s, span.clone()));
    }
}

#[test]
fn cool_is_send_and_sync() {
    fn check<T: Send + Sync>() {}
    check::<Cool>();
    fn send<T: Send>() {}
    send::<PixelScratch>();
}

#[test]
fn variable_bins_are_unsupported() {
    let dir = tempfile::tempdir().unwrap();
    let path = dir.path().join("v.cool");
    write_cool(&path, &spec(10, 15));
    {
        let file = hdf5_metno::File::open_rw(&path).unwrap();
        file.attr("bin-type")
            .unwrap()
            .write_scalar(
                &"variable"
                    .parse::<hdf5_metno::types::VarLenUnicode>()
                    .unwrap(),
            )
            .unwrap();
    }
    assert!(matches!(Cool::open(&path, "/"), Err(Error::Unsupported(_))));
}

#[test]
fn bins_off_the_fixed_grid_are_rejected() {
    let dir = tempfile::tempdir().unwrap();
    let path = dir.path().join("b.cool");
    write_cool(&path, &spec(10, 16));
    {
        let file = hdf5_metno::File::open_rw(&path).unwrap();
        let start = file.dataset("bins/start").unwrap();
        start.write_slice(&[150i32], 1..2).unwrap();
    }
    let err = Cool::open(&path, "/").err().unwrap();
    assert!(
        matches!(
            err,
            Error::Invalid {
                what: "cool bins",
                ..
            }
        ),
        "{err}"
    );
}

#[test]
fn chromosome_offsets_are_checked_before_bin_table_access() {
    let dir = tempfile::tempdir().unwrap();
    let path = dir.path().join("invalid-offset.cool");
    let s = MapSpec::new(&[("first", 1), ("empty", 0)], 1);
    write_cool(&path, &s);
    {
        let file = hdf5_metno::File::open_rw(&path).unwrap();
        // The declared first chromosome spans two bins, but the table
        // holds one row and the index returns to that final row count.
        file.dataset("chroms/length")
            .unwrap()
            .write_raw(&[2i32, 0])
            .unwrap();
        file.dataset("indexes/chrom_offset")
            .unwrap()
            .write_raw(&[0i64, 2, 1])
            .unwrap();
    }
    let err = Cool::open(&path, "/").err().unwrap();
    assert!(
        matches!(
            err,
            Error::Invalid {
                what: "indexes/chrom_offset",
                ..
            }
        ),
        "{err}"
    );
}

fn storage_strategy() -> impl Strategy<Value = Storage> {
    prop_oneof![
        (
            1usize..200,
            any::<bool>(),
            proptest::option::of(1u8..9),
            any::<bool>()
        )
            .prop_map(|(chunk, shuffle, gzip, fletcher)| Storage::Chunked {
                chunk,
                shuffle,
                gzip,
                fletcher,
            }),
        Just(Storage::Contiguous),
        (1usize..100).prop_map(|chunk| Storage::ScaleOffset { chunk }),
    ]
}

proptest! {
    #![proptest_config(ProptestConfig::with_cases(48))]

    #[test]
    fn random_pixels_round_trip(
        raw in proptest::collection::vec((0u64..175, 0u64..175, 1u32..100_000), 0..1500),
        storage in [storage_strategy(), storage_strategy(), storage_strategy()],
        float in any::<bool>(),
        target in 1u64..600,
        ranges in proptest::collection::vec((0.0f64..1.0, 0.0f64..1.0), 4),
    ) {
        let dir = tempfile::tempdir().unwrap();
        let mut s = MapSpec::new(CHROMS, 100);
        let raw: Vec<(u64, u64, f64)> = raw
            .iter()
            .map(|&(a, b, c)| (a, b, if float { f64::from(c) / 8.0 } else { f64::from(c) }))
            .collect();
        s.pixels = common::normalize_pixels(&raw, true);
        s.count_type = if float { CountType::F64 } else { CountType::I32 };
        // Scale-offset is an integer filter.
        if float && matches!(storage[2], Storage::ScaleOffset { .. }) {
            s.storage = [storage[0], storage[1], Storage::Contiguous];
        } else {
            s.storage = storage;
        }
        let cool = write_and_open(dir.path(), &s);
        check_all(&cool, &s, target);
        let n = s.pixels.len() as u64;
        let mut scratch = PixelScratch::new();
        let mut out = PixelBatch::new();
        for (a, b) in ranges {
            let (a, b) = ((a * n as f64) as u64, (b * n as f64) as u64);
            let rows = a.min(b)..a.max(b);
            cool.read_pixels(rows.clone(), &mut scratch, &mut out).unwrap();
            prop_assert_eq!(&out, &expected(&s, rows));
        }
    }
}
