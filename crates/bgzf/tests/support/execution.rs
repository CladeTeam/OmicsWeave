//! Queues with explicit readiness and submission behavior.
use omicsweave_core::exec::{Job, Ordered};
use omicsweave_core::{Error, Result};
use std::collections::VecDeque;

pub struct Ready {
    done: VecDeque<Result<Vec<u8>>>,
    limit: usize,
    accepted: usize,
    reject_after: Option<usize>,
}

impl Ready {
    pub fn boxed(limit: usize, reject_after: Option<usize>) -> Box<dyn Ordered<Vec<u8>>> {
        Box::new(Self {
            done: VecDeque::new(),
            limit,
            accepted: 0,
            reject_after,
        })
    }
}

impl Ordered<Vec<u8>> for Ready {
    fn submit(&mut self, job: Job<Vec<u8>>) -> Result<()> {
        if self
            .reject_after
            .is_some_and(|limit| self.accepted >= limit)
        {
            return Err(Error::invalid("fixture queue", "submission rejected"));
        }
        self.accepted += 1;
        self.done.push_back(job());
        Ok(())
    }
    fn next(&mut self, _: bool) -> Result<Option<Vec<u8>>> {
        self.done.pop_front().transpose()
    }
    fn pending(&self) -> usize {
        self.done.len()
    }
    fn full(&self) -> bool {
        self.pending() >= self.limit
    }
}

pub struct Deferred(pub VecDeque<Job<Vec<u8>>>);
impl Ordered<Vec<u8>> for Deferred {
    fn submit(&mut self, job: Job<Vec<u8>>) -> Result<()> {
        self.0.push_back(job);
        Ok(())
    }
    fn next(&mut self, wait: bool) -> Result<Option<Vec<u8>>> {
        if wait {
            self.0.pop_front().map(|job| job()).transpose()
        } else {
            Ok(None)
        }
    }
    fn pending(&self) -> usize {
        self.0.len()
    }
    fn full(&self) -> bool {
        self.pending() >= 3
    }
}
