//! Sampled identity of memory and warm files; fixtures are outside timing.
use omicsweave_core::sidecar::SourceIdentity;
use omicsweave_core::{FileSource, ReadAt};
use std::hint::black_box;
use std::io::Write;
use std::time::{Duration, Instant};

fn measure(name: &str, source: &(impl ReadAt + ?Sized)) {
    let expected = SourceIdentity::of(source).unwrap();
    println!("{name}.fingerprint\t{:016x}", expected.fingerprint);
    let mut samples = [0.0; 5];
    for sample in &mut samples {
        let start = Instant::now();
        let mut count = 0;
        while start.elapsed() < Duration::from_millis(75) {
            for _ in 0..8 {
                black_box(SourceIdentity::of(black_box(source)).unwrap());
            }
            count += 8;
        }
        *sample = start.elapsed().as_nanos() as f64 / count as f64;
    }
    samples.sort_by(f64::total_cmp);
    assert_eq!(SourceIdentity::of(source).unwrap(), expected);
    println!("{name}\t{:.2}", samples[2]);
}

fn main() {
    for len in [4096, 65_536, 131_072, 131_073, 1_048_576] {
        let bytes: Vec<u8> = (0..len).map(|i| (i % 251) as u8).collect();
        measure(&format!("identity.memory.{len}"), bytes.as_slice());
        let path = std::env::temp_dir().join(format!(
            "omicsweave-core-identity-{}-{len}",
            std::process::id()
        ));
        let mut file = std::fs::OpenOptions::new()
            .write(true)
            .create_new(true)
            .open(&path)
            .unwrap();
        file.write_all(&bytes).unwrap();
        drop(file);
        let source = FileSource::open(&path).unwrap();
        measure(&format!("identity.file.{len}"), &source);
        drop(source);
        std::fs::remove_file(path).unwrap();
    }
}
