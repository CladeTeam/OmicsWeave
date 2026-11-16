use crate::support::{
    CHROMS, MapSpec, Storage, check_all, expected, read, spec, write_and_open, write_cool,
};
use omicsweave_cool::{Cool, PixelBatch, PixelScratch, ReadPath};
use omicsweave_core::Error;

#[test]
fn storage_layouts_tile_pixels_and_read_partial_ranges() {
    for (storage, path) in [
        (Storage::gzip_shuffle(64), ReadPath::Chunked),
        (Storage::plain(128), ReadPath::Chunked),
        (Storage::Contiguous, ReadPath::Contiguous),
    ] {
        let dir = tempfile::tempdir().unwrap();
        let mut s = spec(3000, 2);
        s.storage = [storage; 3];
        let cool = write_and_open(dir.path(), &s);
        assert_eq!(cool.read_paths(), [path; 3], "{storage:?}");
        if path == ReadPath::Contiguous {
            assert_eq!(cool.pixel_spans(300)[0], 0..300);
        }
        for target in [0, 1, 64, 100, 300, 1000, 1 << 20, u64::MAX] {
            check_all(&cool, &s, target);
        }
        assert_eq!(cool.pixel_spans(u64::MAX), vec![0..cool.n_pixels()]);
        for rows in [0..1, 5..70, 63..65, 127..129, 100..1999, 0..cool.n_pixels()] {
            let fast = read(&cool, rows.clone());
            let mut library = PixelBatch::new();
            cool.read_pixels_library(rows.clone(), &mut library)
                .unwrap();
            assert_eq!(fast, library, "storage {storage:?}, rows {rows:?}");
            assert_eq!(
                fast,
                expected(&s, rows.clone()),
                "storage {storage:?}, rows {rows:?}"
            );
        }
    }
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
fn empty_pixels() {
    let dir = tempfile::tempdir().unwrap();
    let s = MapSpec::new(CHROMS, 100);
    let cool = write_and_open(dir.path(), &s);
    assert_eq!(cool.n_pixels(), 0);
    for target in [0, 1, 1000, u64::MAX] {
        assert!(cool.pixel_spans(target).is_empty());
    }
    for library in [false, true] {
        let mut out = PixelBatch {
            bin1: vec![1],
            bin2: vec![2],
            count: vec![3.0],
        };
        if library {
            cool.read_pixels_library(0..0, &mut out).unwrap();
        } else {
            cool.read_pixels(0..0, &mut PixelScratch::new(), &mut out)
                .unwrap();
        }
        assert_eq!(out, PixelBatch::new());
    }
}

#[test]
fn row_boundaries_and_rejected_ranges_allow_buffer_reuse() {
    let dir = tempfile::tempdir().unwrap();
    let s = spec(100, 13);
    let cool = write_and_open(dir.path(), &s);
    let n = cool.n_pixels();
    let mut out = PixelBatch::new();
    let mut scratch = PixelScratch::new();
    let backwards = std::ops::Range { start: 5, end: 2 };
    for library in [false, true] {
        for rows in [n - 1..n + 1, backwards.clone(), n + 1..n + 1, 0..u64::MAX] {
            let result = if library {
                cool.read_pixels_library(rows.clone(), &mut out)
            } else {
                cool.read_pixels(rows.clone(), &mut scratch, &mut out)
            };
            assert!(
                matches!(result, Err(Error::OutOfRange { what: "pixel rows", hi, .. }) if hi == n),
                "library {library}, rows {rows:?}: {result:?}"
            );
            for valid in [0..1, n - 1..n, 0..0, n / 2..n / 2, n..n] {
                if library {
                    cool.read_pixels_library(valid.clone(), &mut out).unwrap();
                } else {
                    cool.read_pixels(valid.clone(), &mut scratch, &mut out)
                        .unwrap();
                }
                assert_eq!(
                    out,
                    expected(&s, valid.clone()),
                    "library {library}, rows {valid:?}"
                );
            }
        }
    }
}

#[test]
fn invalid_pixel_ids_fail_on_read_and_allow_scratch_reuse() {
    for storage in [
        Storage::Contiguous,
        Storage::plain(2),
        Storage::gzip_shuffle(2),
        Storage::ScaleOffset { chunk: 2 },
    ] {
        for name in ["pixels/bin1_id", "pixels/bin2_id"] {
            for bad in [-1i64, 10] {
                let dir = tempfile::tempdir().unwrap();
                let path = dir.path().join("invalid-pixel.cool");
                let mut s = MapSpec::new(&[("chr1", 10)], 1);
                s.storage = [storage; 3];
                s.pixels = vec![(0, 1, 2.0), (2, 3, 4.0)];
                write_cool(&path, &s);
                {
                    let file = hdf5_metno::File::open_rw(&path).unwrap();
                    file.dataset(name)
                        .unwrap()
                        .write_slice(&[bad], 0..1)
                        .unwrap();
                }
                let cool = Cool::open(&path, "/").unwrap();
                let mut scratch = PixelScratch::new();
                let mut out = PixelBatch::new();
                for library in [false, true] {
                    let result = if library {
                        cool.read_pixels_library(0..1, &mut out)
                    } else {
                        cool.read_pixels(0..1, &mut scratch, &mut out)
                    };
                    let what = if bad < 0 { name } else { "cool pixels" };
                    assert!(
                        matches!(&result, Err(Error::Invalid { what: actual, .. }) if *actual == what),
                        "storage {storage:?}, dataset {name}, id {bad}, library {library}: {result:?}"
                    );
                    if library {
                        cool.read_pixels_library(1..2, &mut out).unwrap();
                    } else {
                        cool.read_pixels(1..2, &mut scratch, &mut out).unwrap();
                    }
                    assert_eq!(
                        out,
                        expected(&s, 1..2),
                        "storage {storage:?}, dataset {name}, id {bad}, library {library}"
                    );
                }
            }
        }
    }
}
