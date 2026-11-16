mod support;

use std::process::{Command, Stdio};
use std::sync::{Arc, mpsc};
use std::time::{Duration, Instant};

use omicsweave_core::exec::{Ordered, Pool};

// Lifecycle failures can block a destructor. A separate process lets the
// parent report that failure without retaining a blocked test thread.
fn scenario(name: &str) {
    let directory = support::TempDirectory::new();
    let completed = directory.path().join("completed");
    let mut child = Command::new(std::env::current_exe().unwrap())
        .args([
            "--ignored",
            "--exact",
            "exec_lifecycle_child",
            "--nocapture",
        ])
        .env("OMICSWEAVE_CORE_EXEC_SCENARIO", name)
        .env("OMICSWEAVE_CORE_EXEC_COMPLETED", &completed)
        .stdout(Stdio::inherit())
        .stderr(Stdio::inherit())
        .spawn()
        .unwrap();
    let deadline = Instant::now() + Duration::from_secs(5);
    loop {
        if let Some(status) = child.try_wait().unwrap() {
            assert!(status.success(), "exec lifecycle {name}: {status}");
            assert_eq!(
                std::fs::read_to_string(&completed).unwrap_or_default(),
                name,
                "exec lifecycle {name}: child did not complete its selected scenario"
            );
            return;
        }
        if Instant::now() >= deadline {
            child.kill().unwrap();
            let _ = child.wait();
            panic!("exec lifecycle {name} did not finish within five seconds");
        }
        std::thread::sleep(Duration::from_millis(10));
    }
}

#[test]
fn surviving_queue_rejects_submission_after_workers_exit() {
    scenario("closed_pool");
}

#[test]
fn final_pool_owner_can_be_released_by_its_job() {
    scenario("worker_drops_pool");
}

#[test]
fn same_pool_worker_can_wait_for_nested_jobs() {
    scenario("nested_jobs");
}

#[test]
fn dropping_queue_does_not_wait_for_channel_blocked_jobs() {
    scenario("nonblocking_queue_drop");
}

#[test]
fn worker_can_drop_queue_with_work_behind_it() {
    scenario("worker_drops_queue");
}

#[test]
fn worker_releasing_pool_does_not_join_dependent_worker() {
    scenario("worker_drops_pool_with_dependency");
}

#[test]
#[ignore = "invoked by a bounded parent with an explicit scenario"]
fn exec_lifecycle_child() {
    let name = std::env::var("OMICSWEAVE_CORE_EXEC_SCENARIO")
        .expect("parent supplies the executor scenario");
    match name.as_str() {
        "closed_pool" => {
            let pool = Pool::new(1);
            let mut queue = pool.ordered::<u32>(1);
            queue.submit(Box::new(|| Ok(3))).unwrap();
            drop(pool);
            assert!(queue.submit(Box::new(|| Ok(5))).is_err());
            assert_eq!(queue.pending(), 1);
            // Closing the pool drains accepted jobs, even if their results
            // are still held by a queue that outlives the pool owner.
            assert_eq!(queue.next(true).unwrap(), Some(3));
            assert!(queue.submit(Box::new(|| Ok(7))).is_err());
            assert_eq!(queue.pending(), 0);
            assert_eq!(queue.next(true).unwrap(), None);
        }
        "worker_drops_pool" => {
            let pool = Arc::new(Pool::new(1));
            let mut queue = pool.ordered::<u32>(1);
            let owner = Arc::clone(&pool);
            let (started_tx, started_rx) = mpsc::channel();
            let (continue_tx, continue_rx) = mpsc::channel();
            queue
                .submit(Box::new(move || {
                    started_tx.send(()).unwrap();
                    continue_rx.recv().unwrap();
                    drop(owner);
                    Ok(11)
                }))
                .unwrap();
            started_rx.recv().unwrap();
            drop(pool);
            continue_tx.send(()).unwrap();
            assert_eq!(queue.next(true).unwrap(), Some(11));
            assert_eq!(queue.pending(), 0);
        }
        "nested_jobs" => {
            let pool = Arc::new(Pool::new(1));
            let nested = Arc::clone(&pool);
            let other = Arc::new(Pool::new(1));
            let other_worker = other.run(|| Ok(std::thread::current().id())).unwrap();
            let actual = pool
                .run(move || {
                    assert_eq!(other.run(|| Ok(std::thread::current().id()))?, other_worker);
                    assert_eq!(nested.run(|| Ok(7))?, 7);
                    let mut queue = nested.ordered::<u32>(4);
                    queue.submit(Box::new(|| Ok(11)))?;
                    queue.submit(Box::new(|| {
                        Err(omicsweave_core::Error::invalid("nested", "expected"))
                    }))?;
                    queue.submit(Box::new(|| panic!("nested job panic")))?;
                    queue.submit(Box::new(|| Ok(13)))?;
                    assert_eq!(queue.pending(), 4);
                    assert_eq!(queue.next(true)?, Some(11));
                    assert!(queue.next(true).is_err());
                    assert!(queue.next(true).is_err());
                    assert_eq!(queue.next(true)?, Some(13));
                    assert_eq!(queue.pending(), 0);
                    Ok(17)
                })
                .unwrap();
            assert_eq!(actual, 17);
        }
        "nonblocking_queue_drop" => {
            let pool = Pool::new(1);
            let mut queue = pool.ordered::<()>(1);
            let (started_tx, started_rx) = mpsc::channel();
            let (continue_tx, continue_rx) = mpsc::channel();
            let (done_tx, done_rx) = mpsc::channel();
            queue
                .submit(Box::new(move || {
                    started_tx.send(()).unwrap();
                    continue_rx.recv().unwrap();
                    done_tx.send(()).unwrap();
                    Ok(())
                }))
                .unwrap();
            started_rx.recv().unwrap();
            drop(queue);
            continue_tx.send(()).unwrap();
            done_rx.recv().unwrap();
            drop(pool);
        }
        "worker_drops_queue" => {
            let pool = Pool::new(1);
            let mut outer = pool.ordered::<u32>(1);
            let mut behind = pool.ordered::<u32>(1);
            let (queue_tx, queue_rx) = mpsc::channel();
            let (started_tx, started_rx) = mpsc::channel();
            let (done_tx, done_rx) = mpsc::channel();
            outer
                .submit(Box::new(move || {
                    started_tx.send(()).unwrap();
                    let queue: omicsweave_core::exec::PoolQueue<u32> = queue_rx.recv().unwrap();
                    drop(queue);
                    Ok(5)
                }))
                .unwrap();
            started_rx.recv().unwrap();
            behind
                .submit(Box::new(move || {
                    done_tx.send(()).unwrap();
                    Ok(7)
                }))
                .unwrap();
            queue_tx.send(behind).unwrap();
            assert_eq!(outer.next(true).unwrap(), Some(5));
            done_rx.recv().unwrap();
        }
        "worker_drops_pool_with_dependency" => {
            let pool = Arc::new(Pool::new(2));
            let mut waiting = pool.ordered::<()>(1);
            let mut closing = pool.ordered::<u32>(1);
            let (started_tx, started_rx) = mpsc::channel();
            let (release_tx, release_rx) = mpsc::channel();
            waiting
                .submit(Box::new(move || {
                    started_tx.send(()).unwrap();
                    release_rx.recv().unwrap();
                    Ok(())
                }))
                .unwrap();
            started_rx.recv().unwrap();
            let owner = Arc::clone(&pool);
            let (close_tx, close_rx) = mpsc::channel();
            closing
                .submit(Box::new(move || {
                    close_rx.recv().unwrap();
                    drop(owner);
                    release_tx.send(()).unwrap();
                    Ok(19)
                }))
                .unwrap();
            drop(pool);
            close_tx.send(()).unwrap();
            assert_eq!(closing.next(true).unwrap(), Some(19));
            assert_eq!(waiting.next(true).unwrap(), Some(()));
        }
        "ordered_slots_reuse" => ordered_slots_reuse_case(),
        "ordered_stream" => ordered_stream_case(),
        other => panic!("unknown exec lifecycle scenario: {other}"),
    }
    let completed = std::env::var_os("OMICSWEAVE_CORE_EXEC_COMPLETED")
        .expect("parent supplies the completion path");
    std::fs::write(completed, &name).unwrap();
}

#[test]
fn completed_jobs_wait_for_the_head_and_slots_survive_reuse() {
    scenario("ordered_slots_reuse");
}

fn ordered_slots_reuse_case() {
    let pool = Pool::new(4);
    let mut queue = pool.ordered::<usize>(64);
    for round in 0..8 {
        let (release_tx, release_rx) = mpsc::channel();
        queue
            .submit(Box::new(move || {
                release_rx.recv().unwrap();
                Ok(round * 64)
            }))
            .unwrap();
        for i in 1..64 {
            queue.submit(Box::new(move || Ok(round * 64 + i))).unwrap();
        }
        assert!(queue.full());
        assert_eq!(queue.pending(), 64);
        assert_eq!(queue.next(false).unwrap(), None);
        release_tx.send(()).unwrap();
        for i in 0..64 {
            assert_eq!(queue.next(true).unwrap(), Some(round * 64 + i));
        }
        assert_eq!(queue.pending(), 0);
        assert_eq!(queue.next(false).unwrap(), None);
    }
}

#[test]
fn a_stream_reuses_slots_while_jobs_are_still_pending() {
    scenario("ordered_stream");
}

fn ordered_stream_case() {
    let pool = Pool::new(2);
    let mut queue = pool.ordered::<usize>(7);
    for i in 0..7 {
        queue.submit(Box::new(move || Ok(i))).unwrap();
    }
    for i in 0..1024 {
        assert_eq!(queue.next(true).unwrap(), Some(i));
        queue.submit(Box::new(move || Ok(i + 7))).unwrap();
        assert_eq!(queue.pending(), 7);
    }
    for i in 1024..1031 {
        assert_eq!(queue.next(true).unwrap(), Some(i));
    }
    assert_eq!(queue.pending(), 0);
}
