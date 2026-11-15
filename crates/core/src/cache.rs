//! Caller-owned, explicitly budgeted storage for immutable decoded blocks.
//!
//! Readers opt in with a shared cache and a source session token; arbitrary
//! [`crate::ReadAt`] sources have no immutability guarantee. The caller must
//! keep source contents and indexes unchanged throughout a cached session.
//! Loading runs outside the cache lock. Payload buffers and prepared indexes
//! are owned by this external object, rather than a reader's private cache.
//!
//! Budgets cover admitted payload capacities, their inline value and entry
//! header, including evicted entries with live leases. Loader scratch and
//! allocations made before admission, map bookkeeping and allocator overhead
//! are not a process-RSS budget. FIFO admission is a convenience policy;
//! planners can explicitly evict keys, clear sources and change limits.

use std::any::{Any, TypeId};
use std::collections::{BTreeMap, HashMap};
use std::marker::PhantomData;
use std::ops::Deref;
use std::sync::atomic::{AtomicU64, AtomicUsize, Ordering};
use std::sync::{Arc, Condvar, Mutex, MutexGuard};
use std::thread::ThreadId;

use crate::{Error, Result};

/// A process-local source session. Independent opens get independent tokens;
/// clones of one immutable source session explicitly copy the same token.
#[derive(Clone, Copy, Debug, Eq, Hash, PartialEq)]
pub struct SourceToken(u64);

impl SourceToken {
    pub fn new() -> Self {
        static NEXT: AtomicU64 = AtomicU64::new(1);
        let id = NEXT
            .fetch_update(Ordering::Relaxed, Ordering::Relaxed, |n| n.checked_add(1))
            .expect("source session identifiers exhausted");
        Self(id)
    }

    pub fn get(self) -> u64 {
        self.0
    }
}

impl Default for SourceToken {
    fn default() -> Self {
        Self::new()
    }
}

/// Format-independent identity. Formats assign their own meaning to stream,
/// block and variant; all decoding options affecting payload must be included.
/// Kinds 1/2/3/4/5/6 are FASTA BGZF, 2bit, HIC, BBI, indexed BGZF query blocks
/// and CRAM container/slice data. Each format defines its remaining fields.
#[derive(Clone, Copy, Debug, Eq, Hash, PartialEq)]
pub struct BlockKey {
    pub source: SourceToken,
    pub kind: u32,
    pub stream: [u64; 2],
    pub block: u64,
    pub variant: u64,
}

impl BlockKey {
    pub const fn new(
        source: SourceToken,
        kind: u32,
        stream: [u64; 2],
        block: u64,
        variant: u64,
    ) -> Self {
        Self {
            source,
            kind,
            stream,
            block,
            variant,
        }
    }
}

#[derive(Clone, Copy, Debug, Default, Eq, PartialEq)]
pub struct CacheStats {
    pub hits: u64,
    pub misses: u64,
    pub loads: u64,
    pub waits: u64,
    pub admission_rejections: u64,
    pub evictions: u64,
    pub resident_bytes: usize,
    pub resident_entries: usize,
    /// All charged entries, including those evicted but held by a lease.
    pub live_bytes: usize,
    pub live_entries: usize,
    pub in_flight: usize,
    pub max_bytes: usize,
    pub max_entries: usize,
}

#[derive(Default)]
struct Budget {
    bytes: AtomicUsize,
    entries: AtomicUsize,
    loading: AtomicUsize,
}

struct Entry {
    key: BlockKey,
    payload: Box<dyn Any + Send + Sync>,
    charge: usize,
    // Fields drop in declaration order, including during unwinding: the
    // payload's allocation is released before its accounting permit.
    permit: Option<Permit>,
}

struct Permit {
    budget: Arc<Budget>,
    charge: usize,
}

impl Drop for Permit {
    fn drop(&mut self) {
        self.budget.bytes.fetch_sub(self.charge, Ordering::AcqRel);
        self.budget.entries.fetch_sub(1, Ordering::AcqRel);
    }
}

/// Immutable payload access. Eviction cannot invalidate a lease or release its
/// charged bytes early. A rejected admission still returns an uncached lease.
pub struct CacheLease<T> {
    entry: Arc<Entry>,
    marker: PhantomData<T>,
}

impl<T> Clone for CacheLease<T> {
    fn clone(&self) -> Self {
        Self {
            entry: self.entry.clone(),
            marker: PhantomData,
        }
    }
}

impl<T: Any + Send + Sync> CacheLease<T> {
    fn checked(entry: Arc<Entry>) -> Result<Self> {
        if !entry.payload.is::<T>() {
            return Err(Error::invalid(
                "block cache",
                "key used with a different payload type",
            ));
        }
        Ok(Self {
            entry,
            marker: PhantomData,
        })
    }

    pub fn key(&self) -> BlockKey {
        self.entry.key
    }

    /// Whether this payload was admitted when created; it may since have been
    /// evicted. This is not a live membership test.
    pub fn is_cached(&self) -> bool {
        self.entry.permit.is_some()
    }
}

impl<T: Any + Send + Sync> Deref for CacheLease<T> {
    type Target = T;

    fn deref(&self) -> &T {
        self.entry
            .payload
            .downcast_ref()
            .expect("lease type checked at construction")
    }
}

struct Slot {
    entry: Arc<Entry>,
    ticket: u64,
}

struct Flight {
    payload_type: TypeId,
    owner: ThreadId,
    result: Mutex<Option<Result<Arc<Entry>>>>,
    ready: Condvar,
}

struct State {
    entries: HashMap<BlockKey, Slot>,
    order: BTreeMap<u64, BlockKey>,
    pending: HashMap<BlockKey, Arc<Flight>>,
    ticket: u64,
    stats: CacheStats,
}

impl Drop for State {
    fn drop(&mut self) {
        // Release table storage before calling payload destructors. Vec's
        // drop cleanup also releases the remaining slots if one panics.
        let retired: Vec<_> = std::mem::take(&mut self.entries).into_values().collect();
        drop(retired);
    }
}

/// An external cache shared explicitly by readers. Hits are O(1) average;
/// admission and FIFO eviction are O(log N), with bounded order metadata.
pub struct BlockCache {
    state: Mutex<State>,
    budget: Arc<Budget>,
}

fn lock<T>(mutex: &Mutex<T>) -> MutexGuard<'_, T> {
    mutex.lock().unwrap_or_else(|poison| poison.into_inner())
}

impl BlockCache {
    pub fn new(max_bytes: usize, max_entries: usize) -> Self {
        Self {
            state: Mutex::new(State {
                entries: HashMap::new(),
                order: BTreeMap::new(),
                pending: HashMap::new(),
                ticket: 0,
                stats: CacheStats {
                    max_bytes,
                    max_entries,
                    ..CacheStats::default()
                },
            }),
            budget: Arc::new(Budget::default()),
        }
    }

    pub fn get<T: Any + Send + Sync>(&self, key: BlockKey) -> Result<Option<CacheLease<T>>> {
        let mut state = lock(&self.state);
        let entry = state.entries.get(&key).map(|slot| slot.entry.clone());
        match entry {
            Some(entry) => {
                state.stats.hits = state.stats.hits.saturating_add(1);
                CacheLease::checked(entry).map(Some)
            }
            None => {
                state.stats.misses = state.stats.misses.saturating_add(1);
                Ok(None)
            }
        }
    }

    /// Gets a payload, or runs one loader for concurrent requests of this key.
    /// The loader reports heap capacities of its buffers/prepared metadata;
    /// inline T and entry header are charged additionally by the cache.
    /// Different types for one key return an error, never a type-punned view.
    /// Recursive loading of the same key on one thread is rejected.
    pub fn get_or_load<T, F>(&self, key: BlockKey, load: F) -> Result<CacheLease<T>>
    where
        T: Any + Send + Sync,
        F: FnOnce() -> Result<(T, usize)>,
    {
        let flight = {
            let mut state = lock(&self.state);
            if let Some(slot) = state.entries.get(&key) {
                let entry = slot.entry.clone();
                state.stats.hits = state.stats.hits.saturating_add(1);
                return CacheLease::checked(entry);
            }
            state.stats.misses = state.stats.misses.saturating_add(1);
            if let Some(flight) = state.pending.get(&key).cloned() {
                if flight.payload_type != TypeId::of::<T>() {
                    return Err(Error::invalid(
                        "block cache",
                        "in-flight payload type mismatch",
                    ));
                }
                if flight.owner == std::thread::current().id() {
                    return Err(Error::invalid(
                        "block cache",
                        "recursive load of the same key",
                    ));
                }
                state.stats.waits = state.stats.waits.saturating_add(1);
                drop(state);
                let mut result = lock(&flight.result);
                while result.is_none() {
                    result = flight.ready.wait(result).unwrap_or_else(|p| p.into_inner());
                }
                return CacheLease::checked(result.as_ref().expect("completed flight").clone()?);
            }
            let flight = Arc::new(Flight {
                payload_type: TypeId::of::<T>(),
                owner: std::thread::current().id(),
                result: Mutex::new(None),
                ready: Condvar::new(),
            });
            state.pending.insert(key, flight.clone());
            self.budget.loading.fetch_add(1, Ordering::AcqRel);
            state.stats.loads = state.stats.loads.saturating_add(1);
            flight
        };
        let mut guard = Loading {
            cache: self,
            key,
            flight: flight.clone(),
            complete: false,
        };
        let loaded = load().map(|(value, heap_bytes)| self.admit(key, &flight, value, heap_bytes));
        guard.finish(loaded.clone());
        CacheLease::checked(loaded?)
    }

    fn admit<T: Any + Send + Sync>(
        &self,
        key: BlockKey,
        flight: &Arc<Flight>,
        value: T,
        heap_bytes: usize,
    ) -> Arc<Entry> {
        let charge = heap_bytes
            .checked_add(size_of::<T>())
            .and_then(|n| n.checked_add(size_of::<Entry>()));
        loop {
            let mut state = lock(&self.state);
            let current = state
                .pending
                .get(&key)
                .is_some_and(|f| Arc::ptr_eq(f, flight));
            let fits_limits = current
                && charge.is_some_and(|n| n <= state.stats.max_bytes)
                && state.stats.max_entries > 0
                && state.ticket < u64::MAX;
            let admitted = fits_limits
                && self.budget.bytes.load(Ordering::Acquire)
                    <= state.stats.max_bytes - charge.expect("charge fits")
                && self.budget.entries.load(Ordering::Acquire) < state.stats.max_entries;
            if fits_limits
                && !admitted
                && let Some((_, victim)) = state.order.pop_first()
            {
                let slot = state.entries.remove(&victim).expect("FIFO key exists");
                state.stats.resident_bytes -= slot.entry.charge;
                state.stats.evictions = state.stats.evictions.saturating_add(1);
                // A payload destructor is caller code and may reenter the
                // cache. Release the lock before dropping any payload.
                drop(state);
                drop(slot);
                continue;
            }
            let entry = Arc::new(Entry {
                key,
                payload: Box::new(value),
                charge: charge.unwrap_or(0),
                permit: admitted.then(|| Permit {
                    budget: self.budget.clone(),
                    charge: charge.expect("admitted charge fits"),
                }),
            });
            if admitted {
                self.budget.bytes.fetch_add(entry.charge, Ordering::AcqRel);
                self.budget.entries.fetch_add(1, Ordering::AcqRel);
                let ticket = state.ticket;
                state.ticket += 1;
                state.order.insert(ticket, key);
                state.stats.resident_bytes += entry.charge;
                state.entries.insert(
                    key,
                    Slot {
                        entry: entry.clone(),
                        ticket,
                    },
                );
            } else {
                state.stats.admission_rejections =
                    state.stats.admission_rejections.saturating_add(1);
            }
            return entry;
        }
    }

    /// Configured admission limits, excluding scratch and in-flight decoding.
    pub fn limits(&self) -> (usize, usize) {
        let state = lock(&self.state);
        (state.stats.max_bytes, state.stats.max_entries)
    }

    pub fn stats(&self) -> CacheStats {
        let state = lock(&self.state);
        CacheStats {
            resident_entries: state.entries.len(),
            live_bytes: self.budget.bytes.load(Ordering::Acquire),
            live_entries: self.budget.entries.load(Ordering::Acquire),
            in_flight: self.budget.loading.load(Ordering::Acquire),
            ..state.stats
        }
    }

    pub fn keys(&self) -> Vec<BlockKey> {
        lock(&self.state).entries.keys().copied().collect()
    }

    pub fn evict(&self, key: BlockKey) -> bool {
        let mut state = lock(&self.state);
        let pending = state.pending.remove(&key).is_some();
        if let Some(slot) = state.entries.remove(&key) {
            state.order.remove(&slot.ticket);
            state.stats.resident_bytes -= slot.entry.charge;
            state.stats.evictions = state.stats.evictions.saturating_add(1);
            drop(state);
            drop(slot);
            true
        } else {
            pending
        }
    }

    pub fn clear(&self) {
        let mut state = lock(&self.state);
        state.pending.clear();
        state.order.clear();
        state.stats.resident_bytes = 0;
        let entries = std::mem::take(&mut state.entries);
        drop(state);
        // Consume the table completely before any caller payload is dropped,
        // keeping bucket reclamation independent of payload unwinding.
        let retired: Vec<_> = entries.into_values().collect();
        drop(retired);
    }

    pub fn clear_source(&self, source: SourceToken) {
        let mut state = lock(&self.state);
        state.pending.retain(|key, _| key.source != source);
        let keys: Vec<_> = state
            .entries
            .keys()
            .filter(|k| k.source == source)
            .copied()
            .collect();
        let mut removed = Vec::with_capacity(keys.len());
        for key in keys {
            let slot = state.entries.remove(&key).expect("collected key exists");
            state.order.remove(&slot.ticket);
            state.stats.resident_bytes -= slot.entry.charge;
            removed.push(slot);
        }
        drop(state);
        drop(removed);
    }

    /// Changes admission limits and evicts FIFO entries until possible. Live
    /// leases can temporarily exceed smaller limits; new admissions then fail.
    pub fn resize(&self, max_bytes: usize, max_entries: usize) {
        {
            let mut state = lock(&self.state);
            state.stats.max_bytes = max_bytes;
            state.stats.max_entries = max_entries;
        }
        loop {
            let mut state = lock(&self.state);
            if self.budget.bytes.load(Ordering::Acquire) <= state.stats.max_bytes
                && self.budget.entries.load(Ordering::Acquire) <= state.stats.max_entries
            {
                break;
            }
            let Some((_, key)) = state.order.pop_first() else {
                break;
            };
            let slot = state.entries.remove(&key).expect("FIFO key exists");
            state.stats.resident_bytes -= slot.entry.charge;
            state.stats.evictions = state.stats.evictions.saturating_add(1);
            drop(state);
            drop(slot);
        }
    }
}

struct Loading<'a> {
    cache: &'a BlockCache,
    key: BlockKey,
    flight: Arc<Flight>,
    complete: bool,
}

impl Loading<'_> {
    fn finish(&mut self, result: Result<Arc<Entry>>) {
        let mut state = lock(&self.cache.state);
        if state
            .pending
            .get(&self.key)
            .is_some_and(|f| Arc::ptr_eq(f, &self.flight))
        {
            state.pending.remove(&self.key);
        }
        // Publish completion before new requests can become loaders. A clear
        // can detach this flight, but new requests never wait on that flight.
        *lock(&self.flight.result) = Some(result);
        self.complete = true;
        self.cache.budget.loading.fetch_sub(1, Ordering::AcqRel);
        drop(state);
        self.flight.ready.notify_all();
    }
}

impl Drop for Loading<'_> {
    fn drop(&mut self) {
        if !self.complete {
            self.finish(Err(Error::invalid(
                "block cache",
                "loader unwound before completion",
            )));
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::sync::mpsc;

    fn key(source: SourceToken, block: u64) -> BlockKey {
        BlockKey::new(source, 1, [0, 0], block, 0)
    }

    #[test]
    fn session_and_type_isolation() {
        let cache = BlockCache::new(4096, 8);
        let a = key(SourceToken::new(), 0);
        let b = key(SourceToken::new(), 0);
        let lease = cache.get_or_load(a, || Ok((vec![1u8, 2, 3], 3))).unwrap();
        assert_eq!(&**lease, &[1, 2, 3]);
        assert!(cache.get::<Vec<u8>>(b).unwrap().is_none());
        assert!(cache.get::<String>(a).is_err());
        assert!(cache.get_or_load(a, || Ok((String::new(), 0))).is_err());
        assert_eq!(&**cache.get::<Vec<u8>>(a).unwrap().unwrap(), &[1, 2, 3]);
    }

    #[test]
    fn evicted_leases_keep_their_budget_until_dropped() {
        let cache = BlockCache::new(4096, 1);
        let source = SourceToken::new();
        let a = cache
            .get_or_load(key(source, 0), || Ok((vec![7u8; 16], 16)))
            .unwrap();
        let charge = cache.stats().live_bytes;
        assert!(a.is_cached());
        cache.evict(a.key());
        assert_eq!(cache.stats().resident_entries, 0);
        assert_eq!(cache.stats().live_bytes, charge);
        let b = cache
            .get_or_load(key(source, 1), || Ok((vec![8u8; 16], 16)))
            .unwrap();
        assert!(!b.is_cached());
        assert_eq!(a[0], 7);
        drop(a);
        assert_eq!(cache.stats().live_bytes, 0);
        let c = cache
            .get_or_load(key(source, 2), || Ok((vec![9u8; 16], 16)))
            .unwrap();
        assert!(c.is_cached());
        assert_eq!(b[0], 8);
    }

    #[test]
    fn source_clear_detaches_old_load_without_blocking_new_one() {
        let cache = Arc::new(BlockCache::new(4096, 8));
        let block = key(SourceToken::new(), 0);
        let (started_tx, started_rx) = mpsc::channel();
        let (resume_tx, resume_rx) = mpsc::channel();
        let old_cache = cache.clone();
        let old = std::thread::spawn(move || {
            old_cache
                .get_or_load(block, || {
                    started_tx.send(()).unwrap();
                    resume_rx.recv().unwrap();
                    Ok((vec![1u8], 1))
                })
                .unwrap()
        });
        started_rx.recv().unwrap();
        cache.clear_source(block.source);
        let fresh = cache.get_or_load(block, || Ok((vec![2u8], 1))).unwrap();
        resume_tx.send(()).unwrap();
        let stale = old.join().unwrap();
        assert_eq!(stale[0], 1);
        assert!(!stale.is_cached());
        assert_eq!(fresh[0], 2);
        assert_eq!(cache.get::<Vec<u8>>(block).unwrap().unwrap()[0], 2);
    }

    #[test]
    fn concurrent_missing_key_is_loaded_once() {
        let cache = Arc::new(BlockCache::new(4096, 8));
        let block = key(SourceToken::new(), 0);
        let calls = AtomicUsize::new(0);
        std::thread::scope(|scope| {
            let threads: Vec<_> = (0..12)
                .map(|_| {
                    let cache = cache.clone();
                    let calls = &calls;
                    scope.spawn(move || {
                        cache
                            .get_or_load(block, || {
                                calls.fetch_add(1, Ordering::Relaxed);
                                Ok((vec![42u8], 1))
                            })
                            .unwrap()
                    })
                })
                .collect();
            for thread in threads {
                assert_eq!(thread.join().unwrap()[0], 42);
            }
        });
        assert_eq!(calls.load(Ordering::Relaxed), 1);
        assert_eq!(cache.stats().loads, 1);
    }

    #[test]
    fn waiting_request_observes_loader_failure_then_can_retry() {
        let cache = Arc::new(BlockCache::new(4096, 8));
        let block = key(SourceToken::new(), 0);
        let (started_tx, started_rx) = mpsc::channel();
        let (resume_tx, resume_rx) = mpsc::channel();
        let loading_cache = cache.clone();
        let loader = std::thread::spawn(move || {
            loading_cache
                .get_or_load::<Vec<u8>, _>(block, || {
                    started_tx.send(()).unwrap();
                    resume_rx.recv().unwrap();
                    Err(Error::Inflate("damaged member"))
                })
                .is_err()
        });
        started_rx.recv().unwrap();
        let waiting_cache = cache.clone();
        let waiter = std::thread::spawn(move || {
            waiting_cache
                .get_or_load(block, || panic!("waiter must not load"))
                .map(|_: CacheLease<Vec<u8>>| ())
                .is_err()
        });
        let deadline = std::time::Instant::now() + std::time::Duration::from_secs(2);
        while cache.stats().waits == 0 {
            assert!(
                std::time::Instant::now() < deadline,
                "waiter did not reach flight"
            );
            std::thread::yield_now();
        }
        resume_tx.send(()).unwrap();
        assert!(loader.join().unwrap());
        assert!(waiter.join().unwrap());
        assert_eq!(cache.stats().in_flight, 0);
        let value = cache.get_or_load(block, || Ok((vec![9u8], 1))).unwrap();
        assert_eq!(value[0], 9);
    }

    #[test]
    fn loader_error_panic_and_recursive_load_allow_retry() {
        let cache = BlockCache::new(4096, 8);
        let block = key(SourceToken::new(), 0);
        assert!(
            cache
                .get_or_load::<Vec<u8>, _>(block, || Err(Error::Inflate("bad block")))
                .is_err()
        );
        let panicked = std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| {
            let _ = cache.get_or_load::<Vec<u8>, _>(block, || panic!("loader"));
        }));
        assert!(panicked.is_err());
        assert_eq!(cache.stats().in_flight, 0);
        let value = cache
            .get_or_load(block, || {
                assert!(cache.get_or_load(block, || Ok((vec![0u8], 1))).is_err());
                Ok((vec![6u8], 1))
            })
            .unwrap();
        assert_eq!(value[0], 6);
    }

    #[test]
    fn capacity_overflow_zero_limit_and_resize_never_overadmit() {
        let cache = BlockCache::new(4096, 2);
        let source = SourceToken::new();
        let huge = cache
            .get_or_load(key(source, 0), || Ok((1u8, usize::MAX)))
            .unwrap();
        assert!(!huge.is_cached());
        let first = cache.get_or_load(key(source, 1), || Ok((2u8, 0))).unwrap();
        cache.resize(0, 0);
        assert_eq!(cache.stats().resident_entries, 0);
        assert_eq!(cache.stats().live_entries, 1);
        let second = cache.get_or_load(key(source, 2), || Ok((3u8, 0))).unwrap();
        assert!(!second.is_cached());
        assert_eq!(*first, 2);
        drop(first);
        assert_eq!(cache.stats().live_entries, 0);
    }

    #[test]
    fn destructors_run_without_the_cache_lock() {
        struct Reenter(Arc<BlockCache>);
        impl Drop for Reenter {
            fn drop(&mut self) {
                assert_eq!(self.0.stats().live_entries, 1);
            }
        }
        let cache = Arc::new(BlockCache::new(4096, 1));
        let source = SourceToken::new();
        drop(
            cache
                .get_or_load(key(source, 0), || Ok((Reenter(cache.clone()), 0)))
                .unwrap(),
        );
        drop(cache.get_or_load(key(source, 1), || Ok((0u8, 0))).unwrap());
        assert_eq!(cache.stats().resident_entries, 1);
        cache.clear();
        assert_eq!(cache.stats().live_entries, 0);
    }

    #[test]
    fn panicking_payload_still_releases_its_permit() {
        struct Panicking;
        impl Drop for Panicking {
            fn drop(&mut self) {
                panic!("payload destructor");
            }
        }
        let cache = BlockCache::new(4096, 1);
        drop(
            cache
                .get_or_load(key(SourceToken::new(), 0), || Ok((Panicking, 0)))
                .unwrap(),
        );
        assert_eq!(cache.stats().live_entries, 1);
        let result = std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| cache.clear()));
        assert!(result.is_err());
        assert_eq!(cache.stats().live_entries, 0);
        assert_eq!(cache.stats().live_bytes, 0);
    }

    #[test]
    fn cache_drop_releases_all_payloads_when_one_destructor_panics() {
        struct Tracked {
            panic: bool,
            dropped: Arc<AtomicUsize>,
        }
        impl Drop for Tracked {
            fn drop(&mut self) {
                self.dropped.fetch_add(1, Ordering::Relaxed);
                assert!(!self.panic, "payload destructor");
            }
        }
        let cache = BlockCache::new(4096, 4);
        let budget = cache.budget.clone();
        let dropped = Arc::new(AtomicUsize::new(0));
        let source = SourceToken::new();
        for block in 0..3 {
            drop(
                cache
                    .get_or_load(key(source, block), || {
                        Ok((
                            Tracked {
                                panic: block == 1,
                                dropped: dropped.clone(),
                            },
                            0,
                        ))
                    })
                    .unwrap(),
            );
        }
        assert_eq!(budget.entries.load(Ordering::Acquire), 3);
        let result = std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| drop(cache)));
        assert!(result.is_err());
        assert_eq!(dropped.load(Ordering::Relaxed), 3);
        assert_eq!(budget.entries.load(Ordering::Acquire), 0);
        assert_eq!(budget.bytes.load(Ordering::Acquire), 0);
    }
}
