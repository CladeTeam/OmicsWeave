//! Same public conversion API and fixtures for both source snapshots.
use std::fs::File;
use std::hint::black_box;
use std::io::Write;
use std::path::Path;
use std::time::{Duration, Instant};

use omicsweave_convert::{BedToBigBed, Budget, ChromSizes, Conversion, convert};

fn sparse_twobit(path: &Path, bases: u32) {
    let mut file = File::create(path).unwrap();
    let offset = 16 + 1 + 4 + 4;
    for word in [0x1A41_2743u32, 0, 1, 0] {
        file.write_all(&word.to_le_bytes()).unwrap();
    }
    file.write_all(b"\x04chr1").unwrap();
    file.write_all(&(offset as u32).to_le_bytes()).unwrap();
    for word in [bases, 0, 0, 0] {
        file.write_all(&word.to_le_bytes()).unwrap();
    }
    file.set_len(offset + 16 + (bases as u64).div_ceil(4))
        .unwrap();
}

fn run(src: &Path, dst: &Path, sizes: &Path) -> u64 {
    let mut conversion = BedToBigBed::new(src, dst, ChromSizes::File(sizes.to_owned()));
    conversion.sidecar = false;
    let report = convert(
        Conversion::BedToBigBed(conversion),
        Budget::new(1, 32 << 20),
        None,
    )
    .unwrap();
    assert_eq!(report.records, 1);
    black_box(report.records)
}

fn main() {
    let root = std::env::args()
        .nth(1)
        .expect("pass an isolated fixture directory");
    let dir = Path::new(&root).join(format!("convert-r8-{}", std::process::id()));
    std::fs::create_dir_all(&dir).unwrap();
    let src = dir.join("in.bed");
    let dst = dir.join("out.bb");
    std::fs::write(&src, b"chr1\t1\t5\n").unwrap();
    let small = dir.join("small.2bit");
    let large = dir.join("large.2bit");
    let text = dir.join("chrom.sizes");
    sparse_twobit(&small, 1000);
    sparse_twobit(&large, 1 << 28);
    std::fs::write(&text, b"chr1\t1000\n").unwrap();
    let mut small_checksum = None;
    for (name, sizes) in [
        ("2bit_small", small),
        ("2bit_64mib_packed", large),
        ("text_sizes", text),
    ] {
        run(&src, &dst, &sizes);
        let output = std::fs::read(&dst).unwrap();
        let checksum = omicsweave_core::sidecar::checksum(&output);
        if name == "2bit_small" {
            small_checksum = Some(checksum);
        } else if name == "text_sizes" {
            assert_eq!(Some(checksum), small_checksum);
        }
        println!(
            "fixture {name} input_bytes={} output_bytes={} checksum={checksum:016x}",
            std::fs::metadata(&sizes).unwrap().len(),
            output.len()
        );
        for sample in 0..5 {
            let start = Instant::now();
            let mut iterations = 0u64;
            while start.elapsed() < Duration::from_millis(100) {
                black_box(run(&src, &dst, &sizes));
                iterations += 1;
            }
            println!(
                "{name} sample={sample} ns/iter={:.0} iterations={iterations}",
                start.elapsed().as_nanos() as f64 / iterations as f64
            );
        }
        assert_eq!(std::fs::read(&dst).unwrap(), output);
    }
    std::fs::remove_dir_all(dir).unwrap();
}
