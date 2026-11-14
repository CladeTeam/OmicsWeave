//! Bounds over immutable records, certified once before binary searching.

use std::ops::Range;

pub(crate) type Point = (u32, u32);

#[derive(Debug, Clone, Copy)]
pub(crate) struct Order {
    pub starts: bool,
    pub ends: bool,
    pub overlapping: bool,
    max_span: Option<u32>,
}

impl Default for Order {
    fn default() -> Self {
        Self::UNKNOWN
    }
}

impl Order {
    pub const UNKNOWN: Self = Self {
        starts: false,
        ends: false,
        overlapping: false,
        max_span: None,
    };

    pub fn scan(len: usize, bounds: impl Fn(usize) -> (Point, Point)) -> Self {
        let mut scan = OrderScan::new();
        for i in 0..len {
            let (start, end) = bounds(i);
            scan.push(start, end);
        }
        scan.finish()
    }

    #[inline]
    pub fn candidates(
        self,
        len: usize,
        lo: Point,
        hi: Point,
        bounds: impl Fn(usize) -> (Point, Point),
    ) -> Range<usize> {
        if lo >= hi {
            return 0..0;
        }
        let last = if self.starts {
            partition(len, |i| bounds(i).0 < hi)
        } else {
            len
        };
        let first = if self.ends {
            partition(last, |i| bounds(i).1 <= lo)
        } else if self.starts {
            if let Some(span) = self.max_span {
                let cutoff = lo.1.checked_sub(span);
                partition(last, |i| {
                    let start = bounds(i).0;
                    start.0 < lo.0 || (start.0 == lo.0 && cutoff.is_some_and(|c| start.1 <= c))
                })
            } else {
                0
            }
        } else {
            0
        };
        first..last
    }
}

/// Certifies an [`Order`] one record at a time, so a pass that visits the
/// records anyway (validation, offset building) needs no second pass.
#[derive(Debug, Clone, Copy)]
pub(crate) struct OrderScan {
    order: Order,
    previous: Option<(Point, Point)>,
}

impl OrderScan {
    pub fn new() -> Self {
        Self {
            order: Order {
                starts: true,
                ends: true,
                overlapping: false,
                max_span: Some(0),
            },
            previous: None,
        }
    }

    #[inline]
    pub fn push(&mut self, start: Point, end: Point) {
        let order = &mut self.order;
        if let Some((s, e)) = self.previous {
            order.starts &= s <= start;
            order.ends &= e <= end;
            order.overlapping |= e > start;
        }
        order.max_span = order
            .max_span
            .and_then(|span| (start.0 == end.0).then(|| span.max(end.1.saturating_sub(start.1))));
        self.previous = Some((start, end));
    }

    pub fn finish(self) -> Order {
        self.order
    }
}

#[inline]
fn partition(len: usize, mut predicate: impl FnMut(usize) -> bool) -> usize {
    let (mut lo, mut hi) = (0, len);
    while lo < hi {
        let mid = lo + (hi - lo) / 2;
        if predicate(mid) {
            lo = mid + 1;
        } else {
            hi = mid;
        }
    }
    lo
}
