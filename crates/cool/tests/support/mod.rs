mod pixels;
mod writer;

pub use pixels::{normalize_pixels, sample_pixels};
pub use writer::{CountType, MapSpec, Storage, write_cool, write_mcool};

use omicsweave_cool::{Cool, PixelBatch, PixelScratch};
use std::ops::Range;
use std::path::Path;

pub const CHROMS: &[(&str, u32)] = &[("chr1", 10_000), ("chr2", 7_350), ("chrM", 16)];

pub fn spec(pixels: usize, seed: u64) -> MapSpec {
    let mut s = MapSpec::new(CHROMS, 100);
    s.pixels = sample_pixels(s.n_bins(), pixels, seed, true);
    s
}

pub fn expected(spec: &MapSpec, rows: Range<u64>) -> PixelBatch {
    let mut out = PixelBatch::new();
    for &(a, b, c) in &spec.pixels[rows.start as usize..rows.end as usize] {
        out.bin1.push(a);
        out.bin2.push(b);
        out.count.push(c);
    }
    out
}

pub fn read(cool: &Cool, rows: Range<u64>) -> PixelBatch {
    let mut out = PixelBatch::new();
    cool.read_pixels(rows, &mut PixelScratch::new(), &mut out)
        .unwrap();
    out
}

/// Reads every span with one scratch and checks the spans tile the
/// pixels, each span matching the written data and libhdf5's reader.
pub fn check_all(cool: &Cool, spec: &MapSpec, target: u64) {
    let spans = cool.pixel_spans(target);
    let mut scratch = PixelScratch::new();
    let mut batch = PixelBatch::new();
    let mut library = PixelBatch::new();
    let mut next = 0;
    for span in spans {
        assert_eq!(span.start, next, "target {target}, span {span:?}");
        assert!(!span.is_empty(), "target {target}, span {span:?}");
        next = span.end;
        cool.read_pixels(span.clone(), &mut scratch, &mut batch)
            .unwrap_or_else(|err| {
                panic!(
                    "storage {:?}, target {target}, span {span:?}: {err}",
                    spec.storage
                )
            });
        assert_eq!(
            batch,
            expected(spec, span.clone()),
            "storage {:?}, target {target}, span {span:?}",
            spec.storage
        );
        cool.read_pixels_library(span.clone(), &mut library)
            .unwrap_or_else(|err| panic!("libhdf5, target {target}, span {span:?}: {err}"));
        assert_eq!(library, batch, "libhdf5, target {target}, span {span:?}");
    }
    assert_eq!(next, spec.pixels.len() as u64);
}

pub fn write_and_open(dir: &Path, spec: &MapSpec) -> Cool {
    let path = dir.join("m.cool");
    write_cool(&path, spec);
    Cool::open(&path, "/").unwrap()
}
