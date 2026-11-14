//! The budget of a conversion and how it is divided.
//!
//! Memory is planned, not metered: before anything is read, a conversion states
//! what it keeps for the whole run (indexes, headers, per-file tables), and
//! the rest is split between the pipeline (a bounded number of batches in
//! flight per queue, each of a known size) and, for conversions that must reorder
//! their input, a sort buffer that spills to disk when full. The accounting
//! covers these buffers; allocator overhead, thread stacks and small
//! metadata are not counted and may take the process a little past the
//! budget (see [`RSS_SLACK`]).

use std::path::PathBuf;

use omicsweave_core::{Error, Result};

/// Memory a run may use beyond its budget for what is not accounted:
/// allocator slack, stacks, small per-record metadata.
pub const RSS_SLACK: &str = "10% of the budget plus 64 MiB";

/// Bytes of a text batch handed to one encoding task: large enough to
/// amortize a task, small enough that a few per thread fit the budget.
pub fn text_batch(budget: &Budget) -> usize {
    (budget.memory / (budget.threads.max(1) as u64).saturating_mul(24)).clamp(64 << 10, 4 << 20)
        as usize
}

/// Output buffer of a file being written.
pub const WRITE_BUF: usize = 1 << 20;

/// Threads and memory a conversion may use.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Budget {
    /// Worker threads for the parallel stages (at least one). The reading
    /// and writing stages run on the calling thread, plus one thread to
    /// inflate plain gzip input.
    pub threads: usize,
    /// Bytes of memory for buffers.
    pub memory: u64,
    /// Where spill files go; by default beside the output (never `/tmp`,
    /// which may be small or memory-backed).
    pub tmp_dir: Option<PathBuf>,
}

impl Budget {
    pub fn new(threads: usize, memory: u64) -> Budget {
        Budget {
            threads: threads.max(1),
            memory,
            tmp_dir: None,
        }
    }

    pub fn with_tmp_dir(mut self, dir: impl Into<PathBuf>) -> Budget {
        self.tmp_dir = Some(dir.into());
        self
    }
}

impl Default for Budget {
    /// Every available core and 4 GiB.
    fn default() -> Budget {
        Budget::new(
            std::thread::available_parallelism().map_or(1, |n| n.get()),
            4 << 30,
        )
    }
}

/// How a conversion uses its memory budget.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub struct MemoryPlan {
    /// Kept for the whole run: indexes, headers, tables.
    pub fixed: u64,
    /// Buffers in flight between the stages.
    pub pipeline: u64,
    /// Sort / reorder buffer, spilled to disk when full (0 when the conversion
    /// streams in order).
    pub buffer: u64,
}

impl MemoryPlan {
    pub fn total(&self) -> u64 {
        self.fixed
            .saturating_add(self.pipeline)
            .saturating_add(self.buffer)
    }
}

/// Queue limits derived from a budget.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Sizing {
    /// Batches in flight per queue.
    pub limit: usize,
    /// Bytes held in flight per queue at that limit.
    pub per_queue: u64,
}

/// In-flight limit for a queue of units of `unit` bytes (all buffers a unit
/// holds while in flight), given `share` bytes: between 2 and 4 per thread,
/// never fewer than 2.
pub fn sizing(threads: usize, unit: u64, share: u64) -> Sizing {
    let fits = usize::try_from(share / unit.max(1)).unwrap_or(usize::MAX);
    let limit = fits.clamp(2, threads.saturating_mul(4).max(2));
    Sizing {
        limit,
        per_queue: (limit as u64).saturating_mul(unit),
    }
}

/// Fails when `need` exceeds the budget.
pub fn require(budget: &Budget, need: u64, what: &str) -> Result<()> {
    if need > budget.memory {
        return Err(Error::invalid(
            "budget",
            format!(
                "{what} needs at least {} of memory; the budget is {}",
                human(need),
                human(budget.memory)
            ),
        ));
    }
    Ok(())
}

/// `1.5 GiB`-style sizes for messages.
pub fn human(bytes: u64) -> String {
    const UNITS: [&str; 5] = ["B", "KiB", "MiB", "GiB", "TiB"];
    let mut v = bytes as f64;
    let mut u = 0;
    while v >= 1024.0 && u + 1 < UNITS.len() {
        v /= 1024.0;
        u += 1;
    }
    if u == 0 {
        format!("{bytes} B")
    } else {
        format!("{v:.1} {}", UNITS[u])
    }
}

/// Parses sizes such as `512M`, `8G`, `1.5GiB` or a plain byte count.
pub fn parse_size(text: &str) -> Result<u64> {
    let t = text.trim();
    let split = t
        .find(|c: char| !(c.is_ascii_digit() || c == '.'))
        .unwrap_or(t.len());
    let (num, unit) = t.split_at(split);
    let value: f64 = num
        .parse()
        .map_err(|_| Error::invalid("size", format!("{text:?} is not a size")))?;
    let mult: u64 = match unit.trim().to_ascii_lowercase().as_str() {
        "" | "b" => 1,
        "k" | "kb" | "kib" => 1 << 10,
        "m" | "mb" | "mib" => 1 << 20,
        "g" | "gb" | "gib" => 1 << 30,
        "t" | "tb" | "tib" => 1 << 40,
        _ => return Err(Error::invalid("size", format!("unknown unit in {text:?}"))),
    };
    if !num.contains('.') {
        return num
            .parse::<u64>()
            .ok()
            .and_then(|n| n.checked_mul(mult))
            .ok_or_else(|| Error::invalid("size", format!("{text:?} is out of range")));
    }
    let bytes = value * mult as f64;
    if !bytes.is_finite() || bytes < 0.0 || bytes >= u64::MAX as f64 {
        return Err(Error::invalid("size", format!("{text:?} is out of range")));
    }
    Ok(bytes as u64)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn sizes_parse() {
        assert_eq!(parse_size("8G").unwrap(), 8 << 30);
        assert_eq!(parse_size("512MiB").unwrap(), 512 << 20);
        assert_eq!(parse_size("1.5g").unwrap(), 3 << 29);
        assert_eq!(parse_size("100").unwrap(), 100);
        assert!(parse_size("x").is_err());
    }

    #[test]
    fn limits_stay_within_bounds() {
        assert_eq!(sizing(8, 1 << 20, 1 << 40).limit, 32);
        assert_eq!(sizing(8, 1 << 20, 0).limit, 2);
        assert_eq!(sizing(8, 1 << 20, 10 << 20).limit, 10);
    }

    #[test]
    fn extreme_sizes_and_threads_do_not_wrap() {
        assert_eq!(parse_size("18446744073709551615").unwrap(), u64::MAX);
        for invalid in ["18446744073709551616", "18446744073709551615G", "1e999G"] {
            assert!(parse_size(invalid).is_err());
        }
        let budget = Budget::new(usize::MAX, u64::MAX);
        assert_eq!(text_batch(&budget), 64 << 10);
        assert_eq!(sizing(usize::MAX, u64::MAX, u64::MAX).per_queue, u64::MAX);
        assert_eq!(sizing(0, 1, 0).limit, 2);
        assert_eq!(
            MemoryPlan {
                fixed: u64::MAX,
                pipeline: 1,
                buffer: 1
            }
            .total(),
            u64::MAX
        );
    }
}
