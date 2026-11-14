//! Identical-source baseline/candidate FAI parsing benchmark; no fixture I/O is timed.
use omicsweave_fasta::Fai;
use std::{
    hint::black_box,
    time::{Duration, Instant},
};

fn main() {
    for duplicate in [false, true] {
        let mut text = String::new();
        let count = 32_768;
        for i in 0..count {
            let name = if duplicate { i % 16 } else { i };
            text.push_str(&format!("chr{name}\t1000000\t123456789\t60\t61\n"));
        }
        let expected = if duplicate { 16 } else { count };
        assert_eq!(Fai::parse(text.as_bytes()).unwrap().len(), expected);
        let checksum = text
            .as_bytes()
            .iter()
            .fold(0xcbf29ce484222325_u64, |h, &b| {
                (h ^ b as u64).wrapping_mul(0x100000001b3)
            });
        let mut samples = [0.0_f64; 7];
        for sample in &mut samples {
            let start = Instant::now();
            let mut n = 0;
            while start.elapsed() < Duration::from_millis(200) {
                black_box(Fai::parse(black_box(text.as_bytes())).unwrap());
                n += 1;
            }
            *sample = start.elapsed().as_nanos() as f64 / n as f64;
        }
        samples.sort_by(f64::total_cmp);
        println!(
            "fai.parse.duplicate{duplicate}\trows={count}\tbytes={}\tfnv1a64={checksum:016x}\tmedian_ns={:.2}",
            text.len(),
            samples[3]
        );
    }
}
