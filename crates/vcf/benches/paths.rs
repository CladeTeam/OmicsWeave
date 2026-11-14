//! In-memory whole-stream and open-only workloads. All fixture generation,
//! htslib conversion and output validation happen before the timed groups.

#[path = "../tests/common/mod.rs"]
mod common;

use std::alloc::{GlobalAlloc, Layout, System};
use std::hint::black_box;
use std::sync::atomic::{AtomicBool, AtomicU64, Ordering};
use std::time::{Duration, Instant};

use omicsweave_vcf::{Record, VariantReader};

struct CountingAllocator;
static COUNT: AtomicBool = AtomicBool::new(false);
static ALLOCS: AtomicU64 = AtomicU64::new(0);
static BYTES: AtomicU64 = AtomicU64::new(0);

fn allocated(size: usize) {
    if COUNT.load(Ordering::Relaxed) {
        ALLOCS.fetch_add(1, Ordering::Relaxed);
        BYTES.fetch_add(size as u64, Ordering::Relaxed);
    }
}

// SAFETY: all operations delegate to System with the exact caller-provided
// pointers and layouts; atomics only count requests and never access storage.
unsafe impl GlobalAlloc for CountingAllocator {
    unsafe fn alloc(&self, layout: Layout) -> *mut u8 {
        allocated(layout.size());
        // SAFETY: GlobalAlloc's caller supplies a valid allocation layout.
        unsafe { System.alloc(layout) }
    }
    unsafe fn alloc_zeroed(&self, layout: Layout) -> *mut u8 {
        allocated(layout.size());
        // SAFETY: GlobalAlloc's caller supplies a valid allocation layout.
        unsafe { System.alloc_zeroed(layout) }
    }
    unsafe fn dealloc(&self, ptr: *mut u8, layout: Layout) {
        // SAFETY: the caller supplies the System allocation and its layout.
        unsafe { System.dealloc(ptr, layout) }
    }
    unsafe fn realloc(&self, ptr: *mut u8, layout: Layout, size: usize) -> *mut u8 {
        allocated(size);
        // SAFETY: the caller supplies a live System allocation, its layout,
        // and a valid replacement size under GlobalAlloc's contract.
        unsafe { System.realloc(ptr, layout, size) }
    }
}

#[global_allocator]
static ALLOCATOR: CountingAllocator = CountingAllocator;

fn fixture(version: &str, samples: usize, wide: bool, symbolic: bool) -> String {
    let mut text = format!(
        "##fileformat=VCFv{version}\n##contig=<ID=chr1,length=1000000>\n\
         ##INFO=<ID=END,Number=1,Type=Integer,Description=\"END\">\n\
         ##INFO=<ID=SVLEN,Number=A,Type=Integer,Description=\"SVLEN\">\n\
         ##FORMAT=<ID=GT,Number=1,Type=String,Description=\"GT\">\n\
         ##FORMAT=<ID=LEN,Number=1,Type=Integer,Description=\"LEN\">\n"
    );
    for key in 0..6 {
        text.push_str(&format!(
            "##FORMAT=<ID=F{key},Number=1,Type=Integer,Description=\"F\">\n"
        ));
    }
    text.push_str("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT");
    for s in 0..samples {
        text.push_str(&format!("\tS{s}"));
    }
    text.push('\n');
    for record in 0..1024 {
        let (alleles, info) = if symbolic {
            (
                "C,C,C,C,C,C,<DEL>,<DEL>,<NON_REF>",
                "END=50;SVLEN=1,2,3,4,5,6,70,900,10",
            )
        } else {
            ("C", ".")
        };
        let format = if wide {
            "F0:F1:GT:F2:F3:F4:F5:LEN"
        } else {
            "GT"
        };
        text.push_str(&format!(
            "chr1\t{}\t.\tA\t{alleles}\t50\tPASS\t{info}\t{format}",
            record + 1
        ));
        for s in 0..samples {
            let gt = match s % 4 {
                0 => "0|1",
                1 => "0/0",
                2 => "1",
                _ => "./.",
            };
            if wide {
                text.push_str(&format!("\t1:2:{gt}:3:4:5:6:90"));
            } else {
                text.push('\t');
                text.push_str(gt);
            }
        }
        text.push('\n');
    }
    text
}

fn digest(bytes: &[u8]) -> (usize, u64) {
    let mut reader = VariantReader::new(bytes).unwrap();
    let mut rec = Record::default();
    let mut count = 0;
    let mut hash = 0u64;
    let mut push = |v: u64| hash = hash.wrapping_mul(0x100000001b3).wrapping_add(v);
    while reader.read(&mut rec).unwrap() {
        count += 1;
        push(rec.rid as u64);
        push(rec.pos as u64);
        push(rec.rlen as u64);
        push(rec.variant_types() as u64);
        push(rec.n_alleles() as u64);
        for allele in rec.alleles() {
            push(allele.len() as u64);
            for &b in allele {
                push(b as u64);
            }
        }
        if let Some(gt) = rec.genotypes() {
            push(gt.width as u64);
            for &v in gt.values {
                push(v as u64);
            }
        }
    }
    (count, hash)
}

fn measure(mut operation: impl FnMut()) -> f64 {
    operation();
    let mut groups = Vec::new();
    for _ in 0..5 {
        let start = Instant::now();
        let mut repeats = 0u64;
        while start.elapsed() < Duration::from_millis(100) {
            operation();
            repeats += 1;
        }
        groups.push(start.elapsed().as_secs_f64() * 1e9 / repeats as f64);
    }
    groups.sort_by(f64::total_cmp);
    groups[2]
}

fn main() {
    let selected = std::env::args().skip(1).find(|arg| !arg.starts_with('-'));
    let dir = tempfile::tempdir().unwrap();
    for (name, version, samples, wide, symbolic, binary) in [
        ("normal", "4.3", 4, false, false, false),
        ("wide-format", "4.3", 256, true, false, false),
        ("svlen", "4.5", 16, true, true, false),
        ("bcf43", "4.3", 256, true, false, true),
        ("bcf44", "4.4", 256, true, false, true),
    ] {
        if selected.as_deref().is_some_and(|s| s != name) {
            continue;
        }
        let text = fixture(version, samples, wide, symbolic);
        let path = dir.path().join(format!("{name}.vcf"));
        std::fs::write(&path, &text).unwrap();
        let expected = common::htslib(&path);
        assert_eq!(
            common::ours(text.as_bytes()),
            expected,
            "{name}: VCF oracle"
        );
        let expected_digest = digest(text.as_bytes());
        let bytes = if binary {
            let bcf = dir.path().join(format!("{name}.bcf"));
            common::htslib_to_bcf(&path, &bcf);
            let bytes = std::fs::read(bcf).unwrap();
            assert_eq!(common::ours(&bytes), expected, "{name}: BCF oracle");
            bytes
        } else {
            text.into_bytes()
        };
        assert_eq!(digest(&bytes), expected_digest);
        ALLOCS.store(0, Ordering::Relaxed);
        BYTES.store(0, Ordering::Relaxed);
        COUNT.store(true, Ordering::Relaxed);
        assert_eq!(black_box(digest(black_box(&bytes))), expected_digest);
        COUNT.store(false, Ordering::Relaxed);
        let allocations = ALLOCS.load(Ordering::Relaxed);
        let allocated_bytes = BYTES.load(Ordering::Relaxed);
        let ns = measure(|| {
            assert_eq!(black_box(digest(black_box(&bytes))), expected_digest);
        });
        let open_ns = measure(|| {
            black_box(VariantReader::new(black_box(&bytes[..])).unwrap());
        });
        println!(
            "case={name} records={} input_bytes={} digest={:016x} ns_per_scan={ns:.0} ns_per_record={:.2} open_ns={open_ns:.0} allocation_requests={allocations} requested_bytes={allocated_bytes}",
            expected_digest.0,
            bytes.len(),
            expected_digest.1,
            ns / expected_digest.0 as f64
        );
    }
}
