//! Reused string-token evaluation; fixture/oracle checks precede timing.

#[path = "../tests/cfilter/mod.rs"]
mod cfilter;

use std::alloc::{GlobalAlloc, Layout, System};
use std::hint::black_box;
use std::sync::atomic::{AtomicBool, AtomicU64, Ordering};
use std::time::{Duration, Instant};

use omicsweave_vcf::{Filter, Header, Record, VariantReader};

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

// SAFETY: every operation forwards the same allocation, pointer and layout
// to System. The atomics count requests without accessing allocated memory.
unsafe impl GlobalAlloc for CountingAllocator {
    unsafe fn alloc(&self, layout: Layout) -> *mut u8 {
        allocated(layout.size());
        // SAFETY: GlobalAlloc's caller supplies a valid layout.
        unsafe { System.alloc(layout) }
    }

    unsafe fn alloc_zeroed(&self, layout: Layout) -> *mut u8 {
        allocated(layout.size());
        // SAFETY: GlobalAlloc's caller supplies a valid layout.
        unsafe { System.alloc_zeroed(layout) }
    }

    unsafe fn dealloc(&self, ptr: *mut u8, layout: Layout) {
        // SAFETY: GlobalAlloc's caller supplies the live allocation/layout.
        unsafe { System.dealloc(ptr, layout) }
    }

    unsafe fn realloc(&self, ptr: *mut u8, layout: Layout, size: usize) -> *mut u8 {
        allocated(size);
        // SAFETY: the caller supplies a live allocation and valid new size.
        unsafe { System.realloc(ptr, layout, size) }
    }
}

#[global_allocator]
static ALLOCATOR: CountingAllocator = CountingAllocator;

fn checksum(bytes: &[u8]) -> u64 {
    bytes.iter().fold(0u64, |hash, &byte| {
        hash.wrapping_mul(0x100000001b3)
            .wrapping_add(u64::from(byte))
    })
}

fn scan(filter: &mut Filter, header: &Header, records: &[Record]) -> u64 {
    records.iter().fold(0u64, |hash, record| {
        hash.wrapping_mul(0x100000001b3)
            .wrapping_add(u64::from(filter.test(header, black_box(record)).unwrap()))
    })
}

fn main() {
    let dir = tempfile::tempdir().unwrap();
    for (name, width, numeric) in [
        ("short-literal", 5, false),
        ("long-literal", 256, false),
        ("numeric-control", 5, true),
    ] {
        let literal = "a".repeat(width);
        let expr = if numeric {
            "QUAL>20".to_owned()
        } else {
            format!("INFO/TXT=\"{literal}\"")
        };
        let mut text = "##fileformat=VCFv4.4\n##contig=<ID=chr1>\n\
            ##INFO=<ID=TXT,Number=1,Type=String>\n\
            #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
            .to_owned();
        for i in 0..128 {
            let value = if i % 2 == 0 { &literal } else { "beta" };
            let qual = if i % 2 == 0 { 30 } else { 10 };
            text.push_str(&format!(
                "chr1\t{}\t.\tA\tC\t{qual}\tPASS\tTXT={value}\n",
                i + 1
            ));
        }
        let path = dir.path().join(format!("{name}.vcf"));
        std::fs::write(&path, &text).unwrap();
        let expected = cfilter::c_filter(&path, &expr);
        assert!(matches!(&expected, cfilter::Run::Records(outcomes)
            if outcomes.len() == 128
                && outcomes.iter().enumerate().all(|(i, outcome)|
                    *outcome == cfilter::Outcome::Pass(i % 2 == 0))));
        assert_eq!(cfilter::ours(text.as_bytes(), &expr), expected);

        let mut reader = VariantReader::new(text.as_bytes()).unwrap();
        let mut records = Vec::new();
        let mut record = Record::default();
        while reader.read(&mut record).unwrap() {
            records.push(record.clone());
        }
        let header = reader.header();
        let mut filter = Filter::new(header, &expr).unwrap();
        let digest = scan(&mut filter, header, &records);
        ALLOCS.store(0, Ordering::Relaxed);
        BYTES.store(0, Ordering::Relaxed);
        COUNT.store(true, Ordering::Relaxed);
        let counted_digest = scan(&mut filter, header, &records);
        COUNT.store(false, Ordering::Relaxed);
        assert_eq!(counted_digest, digest);
        let allocations = ALLOCS.load(Ordering::Relaxed);
        let bytes = BYTES.load(Ordering::Relaxed);
        let mut samples = [0.0f64; 5];
        for sample in &mut samples {
            let start = Instant::now();
            let mut repeats = 0u64;
            while start.elapsed() < Duration::from_millis(100) {
                assert_eq!(black_box(scan(&mut filter, header, &records)), digest);
                repeats += 1;
            }
            *sample = start.elapsed().as_secs_f64() * 1e9 / repeats as f64;
        }
        samples.sort_by(f64::total_cmp);
        println!(
            "case={name} records={} input_bytes={} input_checksum={:016x} expression_checksum={:016x} output_checksum={digest:016x} ns_per_scan={:.2} allocations={allocations} allocated_bytes={bytes}",
            records.len(),
            text.len(),
            checksum(text.as_bytes()),
            checksum(expr.as_bytes()),
            samples[2]
        );
    }
}
