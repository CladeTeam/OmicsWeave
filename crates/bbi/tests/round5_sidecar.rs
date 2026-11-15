use std::alloc::{GlobalAlloc, Layout, System};
use std::cell::Cell;
use std::io::Cursor;

use omicsweave_bbi::{BbiReader, BbiSidecar, BigWigWriter, Chrom, Level, WriteOptions};
use omicsweave_core::{Locus, Region, TrackRead};

struct Counted;
thread_local! {
    static TRACK: Cell<bool> = const { Cell::new(false) };
    static LARGEST: Cell<usize> = const { Cell::new(0) };
}

fn record(size: usize) {
    let _ = TRACK.try_with(|track| {
        if track.get() {
            let _ = LARGEST.try_with(|largest| largest.set(largest.get().max(size)));
        }
    });
}

// SAFETY: Every allocation/deallocation is forwarded to System unchanged;
// thread-local counters only observe the requested size and own no allocation.
unsafe impl GlobalAlloc for Counted {
    unsafe fn alloc(&self, layout: Layout) -> *mut u8 {
        record(layout.size());
        // SAFETY: The caller supplies the GlobalAlloc layout contract.
        unsafe { System.alloc(layout) }
    }
    unsafe fn alloc_zeroed(&self, layout: Layout) -> *mut u8 {
        record(layout.size());
        // SAFETY: Forward the unchanged layout and zeroed allocation contract.
        unsafe { System.alloc_zeroed(layout) }
    }
    unsafe fn dealloc(&self, ptr: *mut u8, layout: Layout) {
        // SAFETY: Allocation provenance and layout are forwarded unchanged.
        unsafe { System.dealloc(ptr, layout) }
    }
    unsafe fn realloc(&self, ptr: *mut u8, layout: Layout, new_size: usize) -> *mut u8 {
        record(new_size);
        // SAFETY: The caller supplies System's pointer, layout and size contract.
        unsafe { System.realloc(ptr, layout, new_size) }
    }
}

#[global_allocator]
static ALLOCATOR: Counted = Counted;

fn measured_build(bytes: &[u8]) -> (omicsweave_bbi::Result<BbiSidecar>, usize) {
    struct Stop;
    impl Drop for Stop {
        fn drop(&mut self) {
            TRACK.with(|track| track.set(false));
        }
    }
    LARGEST.with(|largest| largest.set(0));
    TRACK.with(|track| track.set(true));
    let stop = Stop;
    let result = BbiSidecar::build(bytes);
    drop(stop);
    (result, LARGEST.with(Cell::get))
}

fn fixture(n: u32, items: u32, zooms: Vec<u32>) -> (Vec<u8>, BbiSidecar) {
    let mut writer = BigWigWriter::create(
        Cursor::new(Vec::new()),
        vec![Chrom {
            name: "c".into(),
            len: n * 8 + 8,
        }],
        WriteOptions {
            items_per_slot: Some(items),
            block_size: 4,
            zoom_reductions: Some(zooms),
            ..WriteOptions::default()
        },
    )
    .unwrap();
    for i in 0..n {
        writer.add(0, i * 8, i * 8 + 4, (i % 97) as f32).unwrap();
    }
    let (out, sidecar) = writer.finish().unwrap();
    (out.into_inner(), sidecar)
}

#[test]
fn sidecar_small_blocks_ignore_a_large_allocation_hint_and_keep_zoom_items() {
    let (mut bytes, mut expected) = fixture(64, 4, vec![16, 128]);
    assert_eq!(BbiSidecar::build(&bytes).unwrap(), expected);
    for ceiling in [u32::MAX, 256 * 1024, 256 * 1024 + 1] {
        bytes[52..56].copy_from_slice(&ceiling.to_le_bytes());
        expected.uncompress_buf_size = ceiling;
        let (got, largest) = measured_build(&bytes);
        assert_eq!(got.unwrap(), expected);
        assert!(largest <= 256 * 1024, "requested {largest} bytes");
    }
    assert!(expected.levels.iter().all(|level| !level.blocks.is_empty()));
}

#[test]
fn sidecar_grows_to_actual_output_and_keeps_the_declared_ceiling() {
    let (bytes, expected) = fixture(12_000, 12_000, vec![]);
    assert!(expected.levels[0].blocks[0].inflated > 128 * 1024);
    let (got, largest) = measured_build(&bytes);
    assert_eq!(got.unwrap(), expected);
    assert!(largest <= 512 * 1024, "requested {largest} bytes");

    let mut too_small = bytes.clone();
    too_small[52..56].copy_from_slice(&32_768u32.to_le_bytes());
    assert!(BbiSidecar::build(&too_small).is_err());
    let block = &expected.levels[0].blocks[0];
    let mut corrupt = bytes.clone();
    corrupt[block.offset as usize] ^= 0xff;
    assert!(BbiSidecar::build(&corrupt).is_err());
    assert_eq!(BbiSidecar::build(&bytes).unwrap(), expected);
}

#[test]
fn empty_sidecar_does_not_allocate_an_inflate_buffer() {
    let (mut bytes, mut expected) = fixture(0, 4, vec![]);
    bytes[52..56].copy_from_slice(&u32::MAX.to_le_bytes());
    expected.uncompress_buf_size = u32::MAX;
    let (got, largest) = measured_build(&bytes);
    assert_eq!(got.unwrap(), expected);
    assert!(largest < 64 * 1024, "requested {largest} bytes");
    assert!(expected.level(Level::Zoom(usize::MAX)).is_none());
    assert!(expected.level(Level::Zoom(0)).is_none());
}

#[test]
fn writer_tree_fanout_fits_the_encoded_node_count() {
    for block_size in [0, 1, u16::MAX as u32 + 1, u32::MAX] {
        let mut output = Cursor::new(Vec::new());
        assert!(
            BigWigWriter::create(
                &mut output,
                vec![Chrom {
                    name: "c".into(),
                    len: 10
                }],
                WriteOptions {
                    block_size,
                    ..WriteOptions::default()
                },
            )
            .is_err()
        );
        assert!(output.into_inner().is_empty());
    }
    let count = u16::MAX as u32;
    let mut writer = BigWigWriter::create(
        Cursor::new(Vec::new()),
        vec![Chrom {
            name: "c".into(),
            len: count + 1,
        }],
        WriteOptions {
            items_per_slot: Some(1),
            block_size: count,
            compress: false,
            zoom_reductions: Some(vec![]),
            ..WriteOptions::default()
        },
    )
    .unwrap();
    for i in 0..count {
        writer.add(0, i, i + 1, i as f32).unwrap();
    }
    let (output, sidecar) = writer.finish().unwrap();
    assert_eq!(sidecar.levels[0].blocks.len(), count as usize);
    let mut reader = BbiReader::open(output.into_inner()).unwrap();
    let mut out = [0.0; 2];
    reader
        .values(&Locus::from(Region::new("c", count - 2, count)), &mut out)
        .unwrap();
    assert_eq!(out, [(count - 2) as f32, (count - 1) as f32]);
}

#[test]
fn packed_block_write_failure_keeps_the_sink_closed_to_later_writes() {
    use std::io::{self, Seek, SeekFrom, Write};
    use std::sync::Arc;
    use std::sync::atomic::{AtomicBool, AtomicUsize, Ordering};

    struct ToggleSink {
        bytes: Cursor<Vec<u8>>,
        fail: Arc<AtomicBool>,
        calls: Arc<AtomicUsize>,
    }
    impl Write for ToggleSink {
        fn write(&mut self, bytes: &[u8]) -> io::Result<usize> {
            self.calls.fetch_add(1, Ordering::Relaxed);
            if self.fail.load(Ordering::Relaxed) {
                return Err(io::Error::other("sink unavailable"));
            }
            self.bytes.write(bytes)
        }
        fn flush(&mut self) -> io::Result<()> {
            Ok(())
        }
    }
    impl Seek for ToggleSink {
        fn seek(&mut self, from: SeekFrom) -> io::Result<u64> {
            self.bytes.seek(from)
        }
    }
    for compress in [false, true] {
        let fail = Arc::new(AtomicBool::new(false));
        let calls = Arc::new(AtomicUsize::new(0));
        let mut writer = BigWigWriter::create(
            ToggleSink {
                bytes: Cursor::new(Vec::new()),
                fail: fail.clone(),
                calls: calls.clone(),
            },
            vec![Chrom {
                name: "c".into(),
                len: 10,
            }],
            WriteOptions {
                items_per_slot: Some(1),
                compress,
                zoom_reductions: Some(vec![]),
                ..WriteOptions::default()
            },
        )
        .unwrap();
        writer.add(0, 0, 1, 1.0).unwrap();
        fail.store(true, Ordering::Relaxed);
        assert!(writer.add(0, 1, 2, 2.0).is_err());
        let failed_calls = calls.load(Ordering::Relaxed);
        fail.store(false, Ordering::Relaxed);
        assert!(writer.finish().is_err());
        assert_eq!(calls.load(Ordering::Relaxed), failed_calls);
    }
}
