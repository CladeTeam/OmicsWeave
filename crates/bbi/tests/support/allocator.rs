use std::alloc::{GlobalAlloc, Layout, System};
use std::cell::Cell;

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

pub(crate) fn counted<T>(f: impl FnOnce() -> T) -> (T, usize, usize) {
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
