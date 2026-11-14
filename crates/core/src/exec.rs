//! Ordered execution of independent jobs: the seam between a serial writer
//! and the threads that do its heavy, pure work (compressing blocks,
//! encoding records).
//!
//! A writer cuts its output into units on one thread, submits each unit's
//! work as a [`Job`], and takes the results back **in submission order**, so
//! what it writes does not depend on how many threads ran the jobs.
//! [`Inline`] runs every job at submission; [`Pool`] runs them on worker
//! threads shared by any number of [`PoolQueue`]s.
//!
//! The usual loop of a writer:
//!
//! ```
//! use omicsweave_core::exec::{Inline, Ordered};
//! # fn main() -> omicsweave_core::Result<()> {
//! let mut exec: Box<dyn Ordered<Vec<u8>>> = Box::new(Inline::new());
//! let mut written = Vec::new();
//! for unit in [b"ab".to_vec(), b"cd".to_vec()] {
//!     while exec.full() {
//!         written.extend(exec.next(true)?.unwrap());
//!     }
//!     exec.submit(Box::new(move || Ok(unit)))?;
//!     while let Some(done) = exec.next(false)? {
//!         written.extend(done);
//!     }
//! }
//! while let Some(done) = exec.next(true)? {
//!     written.extend(done);
//! }
//! assert_eq!(written, b"abcd");
//! # Ok(())
//! # }
//! ```

use std::cell::Cell;
use std::collections::VecDeque;
use std::panic::{AssertUnwindSafe, catch_unwind};
use std::sync::{Arc, Condvar, Mutex, MutexGuard};
use std::thread::JoinHandle;

use crate::error::{Error, Result};

/// A unit of work: pure, so it may run on any thread.
pub type Job<T> = Box<dyn FnOnce() -> Result<T> + Send + 'static>;

/// Runs jobs and hands their results back in submission order.
pub trait Ordered<T: Send + 'static>: Send {
    /// Accepts `job` without waiting for capacity; an inline queue may
    /// finish it before returning. Callers keep [`Ordered::pending`] under
    /// the limit by taking results while [`Ordered::full`]. An error means
    /// the job was not accepted and pending is unchanged.
    fn submit(&mut self, job: Job<T>) -> Result<()>;

    /// The result of the oldest job not yet returned. With `wait`, blocks
    /// until it is done; without, returns `None` while it is still running.
    /// `None` also when nothing is pending. A job's error is returned in
    /// its turn.
    fn next(&mut self, wait: bool) -> Result<Option<T>>;

    /// Jobs submitted and not yet returned by [`Ordered::next`].
    fn pending(&self) -> usize;

    /// Whether the caller should take a result before submitting again.
    fn full(&self) -> bool;
}

/// Runs each job when it is submitted, on the caller's thread.
pub struct Inline<T> {
    done: VecDeque<Result<T>>,
}

impl<T> Inline<T> {
    pub fn new() -> Inline<T> {
        Inline {
            done: VecDeque::new(),
        }
    }
}

impl<T> Default for Inline<T> {
    fn default() -> Inline<T> {
        Inline::new()
    }
}

impl<T: Send + 'static> Ordered<T> for Inline<T> {
    fn submit(&mut self, job: Job<T>) -> Result<()> {
        self.done.push_back(job());
        Ok(())
    }

    fn next(&mut self, _wait: bool) -> Result<Option<T>> {
        self.done.pop_front().transpose()
    }

    fn pending(&self) -> usize {
        self.done.len()
    }

    fn full(&self) -> bool {
        !self.done.is_empty()
    }
}

type Task = Box<dyn FnOnce() + Send + 'static>;

struct Tasks {
    queue: VecDeque<Task>,
    closed: bool,
}

struct PoolInner {
    tasks: Mutex<Tasks>,
    ready: Condvar,
}

thread_local! {
    /// The pool whose worker this thread is, compared by address only; the
    /// worker's own `Arc` keeps that address alive while it is set.
    static WORKER_POOL: Cell<*const PoolInner> = const { Cell::new(std::ptr::null()) };
}

impl PoolInner {
    fn on_worker(&self) -> bool {
        WORKER_POOL.with(|pool| pool.get() == std::ptr::from_ref(self))
    }

    fn close(&self) {
        lock(&self.tasks).closed = true;
        self.ready.notify_all();
    }
}

/// Restores the thread's previous worker identity when a worker exits.
struct WorkerContext(*const PoolInner);

impl Drop for WorkerContext {
    fn drop(&mut self) {
        WORKER_POOL.with(|pool| pool.set(self.0));
    }
}

fn lock<T>(m: &Mutex<T>) -> MutexGuard<'_, T> {
    // Jobs run under catch_unwind and no lock is held while one runs, so a
    // poisoned lock only means another thread panicked between two plain
    // field updates; the data stays consistent.
    m.lock().unwrap_or_else(|e| e.into_inner())
}

/// Worker threads shared by any number of [`PoolQueue`]s, taking jobs in
/// the order they were submitted across queues; each queue returns its own
/// results in its submission order. A job that submits to its own pool runs
/// that work inline, so waiting for it never needs an idle worker.
///
/// Dropping the pool closes it to new jobs. Dropped from outside its
/// workers, it finishes the jobs already accepted and joins the threads;
/// dropped by one of its own jobs, it leaves the threads to finish the
/// accepted jobs on their own, since joining could wait on that very job.
pub struct Pool {
    inner: Arc<PoolInner>,
    threads: Vec<JoinHandle<()>>,
}

impl Pool {
    /// A pool of `threads` workers (at least one).
    pub fn new(threads: usize) -> Pool {
        let inner = Arc::new(PoolInner {
            tasks: Mutex::new(Tasks {
                queue: VecDeque::new(),
                closed: false,
            }),
            ready: Condvar::new(),
        });
        let mut workers = Vec::with_capacity(threads.max(1));
        for i in 0..threads.max(1) {
            let worker_inner = Arc::clone(&inner);
            match std::thread::Builder::new()
                .name(format!("omicsweave-worker-{i}"))
                .spawn(move || worker(&worker_inner))
            {
                Ok(w) => workers.push(w),
                Err(e) => {
                    inner.close();
                    for w in workers {
                        let _ = w.join();
                    }
                    panic!("spawn worker thread: {e}");
                }
            }
        }
        Pool {
            inner,
            threads: workers,
        }
    }

    /// Number of worker threads.
    pub fn threads(&self) -> usize {
        self.threads.len()
    }

    /// A queue whose results come back in its own submission order, full
    /// once `limit` (at least one) results are pending.
    pub fn ordered<T: Send + 'static>(&self, limit: usize) -> PoolQueue<T> {
        PoolQueue {
            pool: Arc::clone(&self.inner),
            state: Arc::new(QueueState {
                slots: Mutex::new(Slots {
                    first: 0,
                    entries: VecDeque::new(),
                }),
                done: Condvar::new(),
            }),
            next_submit: 0,
            next_return: 0,
            limit: limit.max(1),
        }
    }

    /// Runs `f` on a worker and waits for it; for one-off work outside a
    /// queue. From a worker of this pool, `f` runs inline.
    pub fn run<T: Send + 'static>(
        &self,
        f: impl FnOnce() -> Result<T> + Send + 'static,
    ) -> Result<T> {
        let mut q = self.ordered::<T>(1);
        q.submit(Box::new(f))?;
        q.next(true)?
            .ok_or_else(|| Error::invalid("worker pool", "job lost"))
    }
}

/// Queues `task`, or runs it at once when called from one of the pool's
/// workers; fails once the pool is closed, as no worker would run it.
fn push_task(inner: &PoolInner, task: Task) -> Result<()> {
    let mut tasks = lock(&inner.tasks);
    if tasks.closed {
        return Err(Error::invalid("worker pool", "pool has shut down"));
    }
    if inner.on_worker() {
        // No pool lock may be held while nested work runs (it may drop the
        // pool or submit again).
        drop(tasks);
        task();
    } else {
        tasks.queue.push_back(task);
        drop(tasks);
        inner.ready.notify_one();
    }
    Ok(())
}

impl Drop for Pool {
    fn drop(&mut self) {
        self.inner.close();
        if self.inner.on_worker() {
            // The handles detach; the workers hold the shared state until
            // they have run every accepted job.
            return;
        }
        for t in self.threads.drain(..) {
            let _ = t.join();
        }
    }
}

fn worker(inner: &PoolInner) {
    let _context = WorkerContext(WORKER_POOL.with(|pool| pool.replace(std::ptr::from_ref(inner))));
    loop {
        let task = {
            let mut tasks = lock(&inner.tasks);
            loop {
                if let Some(t) = tasks.queue.pop_front() {
                    break t;
                }
                if tasks.closed {
                    return;
                }
                tasks = inner.ready.wait(tasks).unwrap_or_else(|e| e.into_inner());
            }
        };
        task();
    }
}

struct QueueState<T> {
    slots: Mutex<Slots<T>>,
    done: Condvar,
}

struct Slots<T> {
    /// Sequence of the front slot; all accepted jobs have a reserved slot.
    first: u64,
    entries: VecDeque<Option<Result<T>>>,
}

/// An [`Ordered`] queue running on a [`Pool`].
///
/// Dropping the queue neither waits for nor cancels its accepted jobs: they
/// own what they captured and run to completion, their results discarded.
/// Take results with [`Ordered::next`] to wait for them and see their
/// errors. A queue may outlive its pool; it then refuses new jobs.
/// Result slots retain capacity up to the queue's peak pending count;
/// callers should observe [`Ordered::full`] to bound this allocation.
pub struct PoolQueue<T> {
    pool: Arc<PoolInner>,
    state: Arc<QueueState<T>>,
    next_submit: u64,
    next_return: u64,
    limit: usize,
}

impl<T: Send + 'static> Ordered<T> for PoolQueue<T> {
    fn submit(&mut self, job: Job<T>) -> Result<()> {
        let seq = self.next_submit;
        let next_submit = seq
            .checked_add(1)
            .ok_or_else(|| Error::invalid("worker queue", "submission sequence exhausted"))?;
        usize::try_from(next_submit - self.next_return)
            .map_err(|_| Error::invalid("worker queue", "pending count exceeds usize"))?;
        let state = Arc::clone(&self.state);
        // Reserve before dispatch: a nested task can finish on this thread
        // before push_task returns. VecDeque reuses slots as results drain.
        lock(&state.slots).entries.push_back(None);
        let accepted = push_task(
            &self.pool,
            Box::new(move || {
                let result = catch_unwind(AssertUnwindSafe(job))
                    .unwrap_or_else(|_| Err(Error::invalid("worker", "job panicked")));
                let mut slots = lock(&state.slots);
                // The reserved slot cannot be consumed before completion;
                // submit checked that the whole pending range fits usize.
                let index = (seq - slots.first) as usize;
                slots.entries[index] = Some(result);
                // Only the oldest pending result can wake next(true).
                // Later completions become visible when it is consumed.
                let is_first = index == 0;
                drop(slots);
                if is_first {
                    state.done.notify_one();
                }
            }),
        );
        if let Err(error) = accepted {
            // A rejected task never runs, so only this last reservation is
            // removed. Earlier jobs may already have completed meanwhile.
            let _ = lock(&self.state.slots).entries.pop_back();
            return Err(error);
        }
        self.next_submit = next_submit;
        Ok(())
    }

    fn next(&mut self, wait: bool) -> Result<Option<T>> {
        if self.next_return == self.next_submit {
            return Ok(None);
        }
        let mut slots = lock(&self.state.slots);
        loop {
            if slots.entries.front().is_some_and(Option::is_some) {
                let r = slots.entries.pop_front().unwrap().unwrap();
                slots.first += 1;
                self.next_return += 1;
                return r.map(Some);
            }
            if !wait {
                return Ok(None);
            }
            slots = self
                .state
                .done
                .wait(slots)
                .unwrap_or_else(|e| e.into_inner());
        }
    }

    fn pending(&self) -> usize {
        (self.next_submit - self.next_return) as usize
    }

    fn full(&self) -> bool {
        self.pending() >= self.limit
    }
}

/// Where a writer runs its jobs: inline, or on a shared pool with a limit
/// on pending results, observed through [`Ordered::full`].
#[derive(Clone, Copy)]
pub enum Exec<'p> {
    Inline,
    Pool { pool: &'p Pool, limit: usize },
}

impl Exec<'_> {
    /// An [`Ordered`] queue for this execution mode.
    pub fn queue<T: Send + 'static>(self) -> Box<dyn Ordered<T>> {
        match self {
            Exec::Inline => Box::new(Inline::new()),
            Exec::Pool { pool, limit } => Box::new(pool.ordered(limit)),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn exhausted_sequence_preserves_the_last_accepted_slot() {
        let pool = Pool::new(1);
        let mut queue = pool.ordered::<usize>(1);
        queue.next_submit = u64::MAX - 1;
        queue.next_return = u64::MAX - 1;
        lock(&queue.state.slots).first = u64::MAX - 1;
        queue.submit(Box::new(|| Ok(7))).unwrap();
        assert!(queue.submit(Box::new(|| Ok(9))).is_err());
        assert_eq!(queue.pending(), 1);
        assert_eq!(queue.next(true).unwrap(), Some(7));
        assert_eq!(queue.pending(), 0);
        assert!(queue.submit(Box::new(|| Ok(11))).is_err());
        assert_eq!(queue.next(true).unwrap(), None);
    }

    #[test]
    fn results_come_back_in_submission_order() {
        let pool = Pool::new(4);
        let mut q = pool.ordered::<usize>(8);
        let mut out = Vec::new();
        for i in 0..200usize {
            while q.full() {
                out.push(q.next(true).unwrap().unwrap());
            }
            q.submit(Box::new(move || {
                std::thread::sleep(std::time::Duration::from_micros(
                    ((i * 7919) % 13) as u64 * 50,
                ));
                Ok(i)
            }))
            .unwrap();
        }
        while let Some(v) = q.next(true).unwrap() {
            out.push(v);
        }
        assert_eq!(out, (0..200).collect::<Vec<_>>());
    }

    #[test]
    fn errors_and_panics_surface_in_turn() {
        let pool = Pool::new(2);
        let mut q = pool.ordered::<u32>(4);
        q.submit(Box::new(|| Ok(1))).unwrap();
        q.submit(Box::new(|| Err(Error::invalid("job", "bad"))))
            .unwrap();
        q.submit(Box::new(|| panic!("boom"))).unwrap();
        assert_eq!(q.next(true).unwrap(), Some(1));
        assert!(q.next(true).is_err());
        assert!(q.next(true).is_err());
        assert_eq!(q.next(true).unwrap(), None);
    }

    #[test]
    fn queues_share_a_pool() {
        let pool = Pool::new(3);
        let mut a = pool.ordered::<u32>(2);
        let mut b = pool.ordered::<String>(2);
        a.submit(Box::new(|| Ok(7))).unwrap();
        b.submit(Box::new(|| Ok("x".into()))).unwrap();
        assert_eq!(b.next(true).unwrap().as_deref(), Some("x"));
        assert_eq!(a.next(true).unwrap(), Some(7));
        assert_eq!(pool.run(|| Ok(5u8)).unwrap(), 5);
    }
}
