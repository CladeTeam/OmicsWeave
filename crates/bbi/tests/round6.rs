//! Independent resource accounting and byte equality for writer/sidecar boundaries.
use std::alloc::{GlobalAlloc, Layout, System};
use std::cell::Cell;
use std::hint::black_box;
use std::io::{self, Cursor, Seek, SeekFrom, Write};
use std::time::{Duration, Instant};

use omicsweave_bbi::{
    BBI_SIDECAR_VERSION, BIGWIG_MAGIC, BbiSidecar, BedSchema, BigBedWriter, BigWigWriter,
    BlockInfo, Chrom, Error, Kind, LevelBlocks, WriteOptions,
};
use omicsweave_core::sidecar::{self, SidecarKind, SourceIdentity, checksum};
use omicsweave_core::{ReadAt, Result};

struct CountingAllocator;
thread_local! {
    static TRACK: Cell<bool> = const { Cell::new(false) };
    static LARGEST: Cell<usize> = const { Cell::new(0) };
    static ALLOCATED: Cell<usize> = const { Cell::new(0) };
}
fn record(size: usize) {
    let _ = TRACK.try_with(|track| {
        if track.get() {
            let _ = LARGEST.try_with(|n| n.set(n.get().max(size)));
            let _ = ALLOCATED.try_with(|n| n.set(n.get().saturating_add(size)));
        }
    });
}
// SAFETY: Every allocator operation forwards the original pointer and layout
// to System. Counters use allocation-free thread-local Cells and observe sizes.
unsafe impl GlobalAlloc for CountingAllocator {
    unsafe fn alloc(&self, layout: Layout) -> *mut u8 {
        record(layout.size());
        // SAFETY: Forward the allocator caller's valid layout unchanged.
        unsafe { System.alloc(layout) }
    }
    unsafe fn alloc_zeroed(&self, layout: Layout) -> *mut u8 {
        record(layout.size());
        // SAFETY: Forward the allocator caller's valid layout unchanged.
        unsafe { System.alloc_zeroed(layout) }
    }
    unsafe fn dealloc(&self, ptr: *mut u8, layout: Layout) {
        // SAFETY: Forward the matching original allocation and layout unchanged.
        unsafe { System.dealloc(ptr, layout) }
    }
    unsafe fn realloc(&self, ptr: *mut u8, layout: Layout, new_size: usize) -> *mut u8 {
        record(new_size);
        // SAFETY: Forward the caller's allocation, layout and new-size contract.
        unsafe { System.realloc(ptr, layout, new_size) }
    }
}
#[global_allocator]
static ALLOCATOR: CountingAllocator = CountingAllocator;

fn counted<T>(f: impl FnOnce() -> T) -> (T, usize, usize) {
    struct Stop;
    impl Drop for Stop {
        fn drop(&mut self) {
            TRACK.with(|track| track.set(false));
        }
    }
    LARGEST.with(|n| n.set(0));
    ALLOCATED.with(|n| n.set(0));
    TRACK.with(|track| track.set(true));
    let stop = Stop;
    let out = f();
    drop(stop);
    (out, LARGEST.with(Cell::get), ALLOCATED.with(Cell::get))
}

fn chroms() -> Vec<Chrom> {
    vec![Chrom {
        name: "c".into(),
        len: 200_100,
    }]
}
fn options(explicit: bool, zooms: usize) -> WriteOptions {
    WriteOptions {
        compress: false,
        block_size: 4,
        zoom_levels: zooms,
        zoom_reductions: explicit.then(Vec::new),
        ..WriteOptions::default()
    }
}
fn wig(n: u32, explicit: bool, zooms: usize) -> (Vec<u8>, BbiSidecar) {
    let mut writer =
        BigWigWriter::create(Cursor::new(Vec::new()), chroms(), options(explicit, zooms)).unwrap();
    for i in 0..n {
        writer
            .add(0, i * 8, i * 8 + 4, (i % 97) as f32 - 48.0)
            .unwrap();
    }
    let (out, sidecar) = writer.finish().unwrap();
    (out.into_inner(), sidecar)
}
fn bed(n: u32, explicit: bool) -> (Vec<u8>, BbiSidecar) {
    let mut writer = BigBedWriter::create(
        Cursor::new(Vec::new()),
        chroms(),
        BedSchema::bed(4),
        options(explicit, 0),
    )
    .unwrap();
    for i in 0..n {
        writer.add(0, i * 3, i * 3 + 50, b"x").unwrap();
    }
    let (out, sidecar) = writer.finish().unwrap();
    (out.into_inner(), sidecar)
}

#[test]
fn disabled_derived_zooms_match_explicit_empty_zooms_byte_for_byte() {
    for n in [0, 1, 9_999, 10_000, 20_000] {
        for (bytes, expected) in [
            (wig(n, false, 0), wig(n, true, 0)),
            (bed(n, false), bed(n, true)),
        ] {
            assert_eq!(bytes, expected);
            assert_eq!(BbiSidecar::build(&bytes.0).unwrap(), bytes.1);
        }
    }
}

#[derive(Default)]
struct Discard(u64);
impl Write for Discard {
    fn write(&mut self, bytes: &[u8]) -> io::Result<usize> {
        self.0 += bytes.len() as u64;
        Ok(bytes.len())
    }
    fn flush(&mut self) -> io::Result<()> {
        Ok(())
    }
}
impl Seek for Discard {
    fn seek(&mut self, from: SeekFrom) -> io::Result<u64> {
        match from {
            SeekFrom::Start(at) => self.0 = at,
            _ => return Err(io::Error::other("only absolute seeks are used")),
        }
        Ok(self.0)
    }
}

#[test]
fn disabled_zooms_do_not_retain_the_sampling_runs() {
    let mut wig = BigWigWriter::create(Discard::default(), chroms(), options(false, 0)).unwrap();
    let (_, largest, _) = counted(|| {
        for i in 0..20_000 {
            wig.add(0, i * 8, i * 8 + 4, i as f32).unwrap();
        }
    });
    assert!(largest < 128 * 1024, "largest allocation: {largest}");
    let mut bed = BigBedWriter::create(
        Discard::default(),
        chroms(),
        BedSchema::bed(4),
        options(false, 0),
    )
    .unwrap();
    let (_, largest, _) = counted(|| {
        for i in 0..20_000 {
            bed.add(0, i * 3, i * 3 + 50, b"x").unwrap();
        }
    });
    assert!(largest < 128 * 1024, "largest allocation: {largest}");
}

struct CountReads<'a> {
    bytes: &'a [u8],
    calls: Cell<usize>,
}
impl ReadAt for CountReads<'_> {
    fn len(&self) -> u64 {
        self.bytes.len() as u64
    }
    fn read_exact_at(&self, out: &mut [u8], at: u64) -> Result<()> {
        self.calls.set(self.calls.get() + 1);
        self.bytes.read_exact_at(out, at)
    }
}

#[test]
fn zoom_table_must_fit_before_allocating_or_reading_it() {
    let mut bytes = [0u8; 64];
    bytes[..4].copy_from_slice(&BIGWIG_MAGIC.to_le_bytes());
    bytes[6..8].copy_from_slice(&u16::MAX.to_le_bytes());
    let source = CountReads {
        bytes: &bytes,
        calls: Cell::new(0),
    };
    let (result, largest, _) = counted(|| BbiSidecar::build(&source));
    assert!(matches!(result, Err(Error::OutOfRange { .. })));
    assert_eq!(source.calls.get(), 1);
    assert!(largest <= 128, "largest allocation: {largest}");
}

#[test]
fn sidecar_headers_and_block_ranges_are_validated_against_available_bytes() {
    let source = SourceIdentity {
        len: 100,
        fingerprint: 1,
    };
    let mut body = vec![0, 0];
    body.extend(u16::MAX.to_le_bytes());
    body.extend(0u32.to_le_bytes());
    let sealed = sidecar::seal(SidecarKind::Bbi, BBI_SIDECAR_VERSION, source, &body);
    let (result, largest, _) = counted(|| BbiSidecar::parse(&sealed, source));
    assert!(result.is_err());
    assert!(largest < 4096, "largest allocation: {largest}");
    for (offset, size, valid) in [
        (90, 10, true),
        (100, 0, true),
        (99, 2, false),
        (101, 0, false),
        (u64::MAX, 1, false),
    ] {
        let sidecar = BbiSidecar {
            kind: Kind::BigWig,
            uncompress_buf_size: 0,
            levels: vec![LevelBlocks {
                reduction: 0,
                blocks: vec![BlockInfo {
                    start: (0, 0),
                    end: (0, 1),
                    offset,
                    size,
                    inflated: size,
                    items: 1,
                }],
            }],
        };
        let result = BbiSidecar::parse(&sidecar.to_bytes(source), source);
        if valid {
            assert_eq!(result.unwrap(), sidecar);
            assert_eq!(
                sidecar.levels[0].blocks[0].range(),
                offset..offset + size as u64
            );
        } else {
            assert!(matches!(result, Err(Error::OutOfRange { .. })));
        }
    }
}

#[test]
fn sidecar_v1_encodes_the_maximum_level_count_and_rejects_overflow() {
    let mut sidecar = BbiSidecar {
        kind: Kind::BigWig,
        uncompress_buf_size: 0,
        levels: vec![
            LevelBlocks {
                reduction: 0,
                blocks: Vec::new()
            };
            u16::MAX as usize
        ],
    };
    let body = sidecar.body();
    assert_eq!(u16::from_le_bytes(body[2..4].try_into().unwrap()), u16::MAX);
    assert_eq!(body.len(), 8 + u16::MAX as usize * 16);
    let source = SourceIdentity {
        len: 0,
        fingerprint: 1,
    };
    let sealed = sidecar::seal(SidecarKind::Bbi, BBI_SIDECAR_VERSION, source, &body);
    assert_eq!(BbiSidecar::parse(&sealed, source).unwrap(), sidecar);
    sidecar.levels.push(LevelBlocks {
        reduction: 1,
        blocks: Vec::new(),
    });
    assert!(std::panic::catch_unwind(|| sidecar.body()).is_err());
    assert!(std::panic::catch_unwind(|| sidecar.to_bytes(source)).is_err());
}

#[test]
fn caller_constructed_block_range_rejects_overflow_in_every_profile() {
    let block = BlockInfo {
        start: (0, 0),
        end: (0, 1),
        offset: u64::MAX,
        size: 1,
        inflated: 1,
        items: 1,
    };
    assert!(std::panic::catch_unwind(|| block.range()).is_err());
    let block = BlockInfo { size: 0, ..block };
    assert_eq!(block.range(), u64::MAX..u64::MAX);
}

fn measure<T>(name: &str, mut f: impl FnMut() -> T) {
    for _ in 0..8 {
        black_box(f());
    }
    let mut samples = [0.0; 5];
    for sample in &mut samples {
        let start = Instant::now();
        let mut n = 0;
        while start.elapsed() < Duration::from_millis(100) {
            black_box(f());
            n += 1;
        }
        *sample = start.elapsed().as_nanos() as f64 / n as f64;
    }
    samples.sort_by(f64::total_cmp);
    let (_, largest, allocated) = counted(|| black_box(f()));
    println!(
        "{name}\tns={:.2}\tlargest={largest}\tallocated={allocated}",
        samples[2]
    );
}

#[test]
#[ignore = "paired performance measurement: root runs baseline and candidate serially"]
fn paired_writer_benchmark() {
    for zooms in [0, 2] {
        measure(&format!("writer.create.zooms{zooms}"), || {
            BigWigWriter::create(Discard::default(), chroms(), options(false, zooms)).unwrap()
        });
        for n in [1, 9_999, 10_000, 20_000] {
            let fixture = wig(n, false, zooms);
            println!(
                "writer.wig.n{n}.zooms{zooms}.digest\t{:016x}",
                checksum(&fixture.0)
            );
            measure(&format!("writer.wig.n{n}.zooms{zooms}"), || {
                wig(black_box(n), false, zooms)
            });
        }
    }
    for n in [1, 9_999, 10_000, 20_000] {
        let fixture = bed(n, false);
        println!(
            "writer.bed.n{n}.zooms0.digest\t{:016x}",
            checksum(&fixture.0)
        );
        measure(&format!("writer.bed.n{n}.zooms0"), || {
            bed(black_box(n), false)
        });
    }
}

#[test]
#[ignore = "root compares the same allocator probe across writer snapshots"]
fn paired_warm_flush_allocations() {
    for compress in [false, true] {
        let mut write_options = options(false, 0);
        write_options.compress = compress;
        write_options.items_per_slot = Some(16);
        let mut writer = BigWigWriter::create(Discard::default(), chroms(), write_options).unwrap();
        let warm = 2_050 * 16;
        for i in 0..warm {
            writer.add(0, i * 4, i * 4 + 2, (i % 31) as f32).unwrap();
        }
        // Retained block metadata has ample spare capacity after warm-up;
        // the counted phase isolates repeated block cutting and compression.
        let (_, largest, allocated) = counted(|| {
            for i in warm..warm + 512 * 16 {
                writer.add(0, i * 4, i * 4 + 2, (i % 31) as f32).unwrap();
            }
        });
        println!(
            "writer.warm512flushes.compressed{compress}\tlargest={largest}\tallocated={allocated}"
        );
        black_box(writer.finish().unwrap());
    }
}
