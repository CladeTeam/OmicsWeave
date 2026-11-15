//! Signal summarized in fixed-width bins along a [`Locus`].
//!
//! Bins tile the locus in its own orientation: with resolution `r` and
//! phase `p`, local position `x` (0 at the locus' 5' end) falls in bin
//! `(x + p) / r`. With `p = 0` and a length that is a multiple of `r` every
//! bin is full; otherwise the first and last bins are partial, so any
//! window, offset or grid can be expressed (to align bins to the genome
//! grid of a forward one-span locus, pass `p = start % r`).
//!
//! Statistics follow libBigWig's `bwStats` over full-resolution data
//! (`exact`), bit for bit on full bins: every stat is computed in `f64` the
//! way libBigWig computes it and then rounded to `f32`. An empty bin reads
//! NaN, except coverage, which reads 0.

use crate::error::{Error, Result};
use crate::locus::Locus;

/// A summary statistic of a bin.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, PartialOrd, Ord)]
#[repr(u8)]
pub enum Stat {
    /// Mean over the bases with data.
    Mean = 0,
    /// Smallest value of the data overlapping the bin.
    Min = 1,
    /// Largest value of the data overlapping the bin.
    Max = 2,
    /// Sum over the bases with data of their value.
    Sum = 3,
    /// Sample standard deviation over the bases with data (two-pass); 0 for
    /// a single base.
    Std = 4,
    /// Fraction of the bin's bases with data.
    Coverage = 5,
}

impl Stat {
    pub const ALL: [Stat; 6] = [
        Stat::Mean,
        Stat::Min,
        Stat::Max,
        Stat::Sum,
        Stat::Std,
        Stat::Coverage,
    ];

    pub fn as_str(self) -> &'static str {
        match self {
            Stat::Mean => "mean",
            Stat::Min => "min",
            Stat::Max => "max",
            Stat::Sum => "sum",
            Stat::Std => "std",
            Stat::Coverage => "coverage",
        }
    }

    pub fn parse(s: &str) -> Option<Stat> {
        Stat::ALL.into_iter().find(|stat| stat.as_str() == s)
    }
}

/// Which data a bin query summarizes.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Default)]
pub enum BinSource {
    /// The full-resolution data: exact.
    #[default]
    Exact,
    /// The coarsest stored zoom level finer than half a bin, as libBigWig
    /// picks it, with each zoom record weighted by the fraction of it
    /// inside the bin; the full-resolution data when no level qualifies.
    /// An estimate: far less to read for wide bins.
    Zoom,
}

/// A bin query: `resolution` positions per bin, shifted by `phase`, the
/// statistics to write (one output row each, in this order) and the data
/// they summarize.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct BinSpec<'a> {
    pub resolution: u32,
    pub phase: u32,
    pub stats: &'a [Stat],
    pub source: BinSource,
}

impl<'a> BinSpec<'a> {
    pub fn new(resolution: u32, stats: &'a [Stat]) -> BinSpec<'a> {
        BinSpec {
            resolution,
            phase: 0,
            stats,
            source: BinSource::Exact,
        }
    }

    /// Fails unless `0 < resolution`, `phase < resolution`, and the stats
    /// are non-empty and distinct.
    pub fn check(&self) -> Result<()> {
        if self.resolution == 0 {
            return Err(Error::invalid("bins", "resolution must be positive"));
        }
        if self.phase >= self.resolution {
            return Err(Error::invalid(
                "bins",
                format!(
                    "phase {} must be smaller than the resolution {}",
                    self.phase, self.resolution
                ),
            ));
        }
        if self.stats.is_empty() {
            return Err(Error::invalid("bins", "no statistic asked for"));
        }
        for (i, stat) in self.stats.iter().enumerate() {
            if self.stats[..i].contains(stat) {
                return Err(Error::invalid(
                    "bins",
                    format!("{} asked for twice", stat.as_str()),
                ));
            }
        }
        Ok(())
    }

    /// Bins over a locus of `len` positions: `ceil((len + phase) /
    /// resolution)`, and none for an empty locus.
    pub fn count(&self, len: u32) -> usize {
        if len == 0 {
            return 0;
        }
        usize::try_from((len as u64 + self.phase as u64).div_ceil(self.resolution as u64))
            .expect("bin count exceeds usize")
    }

    /// Output elements: one row of [`BinSpec::count`] per statistic.
    pub fn output_len(&self, len: u32) -> usize {
        self.count(len)
            .checked_mul(self.stats.len())
            .expect("bin output length exceeds usize")
    }
}

/// The bins of a locus in forward coordinates (0 at the forward start of
/// the concatenated spans). For a reverse locus the phase is mirrored, so
/// that reversing the forward bins gives the locus' own bins.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct BinLayout {
    pub len: u32,
    pub resolution: u32,
    /// Phase in forward coordinates.
    pub phase: u32,
    pub count: usize,
    pub reverse: bool,
}

impl BinLayout {
    pub fn new(locus: &Locus<'_>, spec: &BinSpec<'_>) -> Result<BinLayout> {
        spec.check()?;
        locus.check()?;
        let len = locus.len();
        let count = spec.count(len);
        let reverse = locus.is_reverse();
        let res = spec.resolution as u64;
        // Bin k of the reverse locus covers forward [len + p - (k+1)r,
        // len + p - kr); in forward order these are bins of phase
        // n*r - len - p, which lies in [0, r).
        let phase = if reverse && count > 0 {
            (count as u64 * res - len as u64 - spec.phase as u64) as u32
        } else {
            spec.phase
        };
        Ok(BinLayout {
            len,
            resolution: spec.resolution,
            phase,
            count,
            reverse,
        })
    }

    /// Forward bin of forward position `x`.
    #[inline]
    pub fn bin(&self, x: u32) -> usize {
        ((x as u64 + self.phase as u64) / self.resolution as u64) as usize
    }

    /// Forward positions `[start, end)` of forward bin `k`, clipped to the
    /// locus.
    #[inline]
    pub fn bounds(&self, k: usize) -> (u32, u32) {
        let res = self.resolution as u64;
        let lo = (k as u64 * res).saturating_sub(self.phase as u64);
        let hi = ((k as u64 + 1) * res - self.phase as u64).min(self.len as u64);
        (lo as u32, hi as u32)
    }

    /// Bases of forward bin `k` inside the locus.
    #[inline]
    pub fn width(&self, k: usize) -> u32 {
        let (lo, hi) = self.bounds(k);
        hi - lo
    }
}

#[derive(Debug, Clone, Copy)]
struct Acc {
    bases: u64,
    /// `Σ len · value` in `f64`, as libBigWig's mean.
    weighted: f64,
    /// `Σ f32(len) · value`, rounded to `f32` before summing, as
    /// libBigWig's `intSum` multiplies a `uint32_t` by a `float`.
    sum: f64,
    min: f32,
    max: f32,
    any: bool,
}

const EMPTY: Acc = Acc {
    bases: 0,
    weighted: 0.0,
    sum: 0.0,
    min: f32::NAN,
    max: f32::NAN,
    any: false,
};

#[derive(Debug, Clone, Copy)]
struct ZoomAcc {
    /// `Σ validCount · f`
    bases: f64,
    /// `Σ sum · f`
    sum: f64,
    /// `Σ sumSquares · f`
    squares: f64,
    min: f32,
    max: f32,
    any: bool,
}

const ZOOM_EMPTY: ZoomAcc = ZoomAcc {
    bases: 0.0,
    sum: 0.0,
    squares: 0.0,
    min: f32::NAN,
    max: f32::NAN,
    any: false,
};

/// One zoom record: its span in forward locus coordinates (clipped to a
/// span of the locus) and the record's full span, which its statistics
/// describe.
#[derive(Debug, Clone, Copy)]
pub struct ZoomPiece {
    pub start: u32,
    pub end: u32,
    /// Bases of the record itself (`end - start` in the file).
    pub record_len: u32,
    pub bases: u32,
    pub min: f32,
    pub max: f32,
    pub sum: f32,
    pub squares: f32,
}

/// Accumulates data into bins, then writes the asked-for statistics.
/// Reusable work state: the buffers grow to the largest query and are kept.
#[derive(Debug, Default)]
pub struct BinAccumulator {
    layout: Option<BinLayout>,
    acc: Vec<Acc>,
    zoom: Vec<ZoomAcc>,
    zoomed: bool,
    /// `(bin, len, value)` pieces, kept for the standard deviation's second
    /// pass.
    pieces: Vec<(u32, u32, f32)>,
    /// The second pass's squared deviations, reused across queries.
    deviations: Vec<f64>,
    keep_pieces: bool,
}

impl BinAccumulator {
    pub fn new() -> BinAccumulator {
        BinAccumulator::default()
    }

    /// Starts a query of `layout` that will write `stats`; `zoom` selects
    /// [`BinAccumulator::add_zoom`] input instead of
    /// [`BinAccumulator::add`].
    pub fn reset(&mut self, layout: BinLayout, stats: &[Stat], zoom: bool) {
        self.layout = Some(layout);
        self.zoomed = zoom;
        self.keep_pieces = !zoom && stats.contains(&Stat::Std);
        self.pieces.clear();
        if zoom {
            self.zoom.clear();
            self.zoom.resize(layout.count, ZOOM_EMPTY);
        } else {
            self.acc.clear();
            self.acc.resize(layout.count, EMPTY);
        }
    }

    pub fn layout(&self) -> Option<&BinLayout> {
        self.layout.as_ref()
    }

    /// Adds `value` at forward positions `[start, end)` (within the locus,
    /// non-empty). Data must arrive in file order for results to match
    /// libBigWig bit for bit.
    #[inline]
    pub fn add(&mut self, start: u32, end: u32, value: f32) {
        let layout = self.layout.expect("reset before add");
        debug_assert!(!self.zoomed && start < end && end <= layout.len);
        let mut at = start;
        let mut k = layout.bin(start);
        while at < end {
            let (_, hi) = layout.bounds(k);
            let to = hi.min(end);
            let len = to - at;
            let a = &mut self.acc[k];
            a.bases += len as u64;
            a.weighted += len as f64 * value as f64;
            a.sum += (len as f32 * value) as f64;
            if !a.any {
                a.min = value;
                a.max = value;
                a.any = true;
            } else {
                if value < a.min {
                    a.min = value;
                }
                if value > a.max {
                    a.max = value;
                }
            }
            if self.keep_pieces {
                self.pieces.push((k as u32, len, value));
            }
            at = to;
            k += 1;
        }
    }

    /// Adds a zoom record clipped to `[piece.start, piece.end)` of the
    /// locus; each bin takes the fraction of the record inside it.
    pub fn add_zoom(&mut self, piece: ZoomPiece) {
        let layout = self.layout.expect("reset before add_zoom");
        debug_assert!(self.zoomed && piece.start < piece.end && piece.end <= layout.len);
        let mut at = piece.start;
        let mut k = layout.bin(piece.start);
        while at < piece.end {
            let (_, hi) = layout.bounds(k);
            let to = hi.min(piece.end);
            let f = (to - at) as f64 / piece.record_len.max(1) as f64;
            let z = &mut self.zoom[k];
            z.bases += piece.bases as f64 * f;
            z.sum += piece.sum as f64 * f;
            z.squares += piece.squares as f64 * f;
            if !z.any {
                z.min = piece.min;
                z.max = piece.max;
                z.any = true;
            } else {
                if piece.min < z.min {
                    z.min = piece.min;
                }
                if piece.max > z.max {
                    z.max = piece.max;
                }
            }
            at = to;
            k += 1;
        }
    }

    /// Writes one row per stat (`out.len() == stats.len() * count`), in the
    /// locus' orientation.
    ///
    /// # Panics
    ///
    /// Panics if the output shape differs, or exact standard deviation
    /// is requested without selecting it at [`BinAccumulator::reset`].
    pub fn finish(&mut self, stats: &[Stat], out: &mut [f32]) {
        let layout = self.layout.expect("reset before finish");
        let n = layout.count;
        assert_eq!(
            out.len(),
            stats
                .len()
                .checked_mul(n)
                .expect("bin output length exceeds usize"),
            "bin output size"
        );
        assert!(
            self.zoomed || self.keep_pieces || !stats.contains(&Stat::Std),
            "standard deviation must be selected at reset"
        );
        if self.zoomed {
            self.finish_zoom(&layout, stats, out);
        } else {
            self.finish_exact(&layout, stats, out);
        }
        if layout.reverse {
            for row in out.chunks_exact_mut(n.max(1)) {
                row.reverse();
            }
        }
    }

    fn finish_exact(&mut self, layout: &BinLayout, stats: &[Stat], out: &mut [f32]) {
        let n = layout.count;
        if stats.contains(&Stat::Std) {
            // Second pass of libBigWig's intDev: squared deviations from
            // the bin mean, weighted by length.
            self.deviations.clear();
            self.deviations.resize(n, 0.0);
            for &(k, len, value) in &self.pieces {
                let a = &self.acc[k as usize];
                let mean = a.weighted / a.bases as f64;
                let d = value as f64 - mean;
                self.deviations[k as usize] += len as f64 * (d * d);
            }
        }
        for (row, stat) in out.chunks_exact_mut(n.max(1)).zip(stats) {
            for (k, slot) in row.iter_mut().enumerate().take(n) {
                let a = &self.acc[k];
                let v = match (stat, a.any) {
                    (Stat::Coverage, _) => a.bases as f64 / layout.width(k) as f64,
                    (_, false) => f64::NAN,
                    (Stat::Mean, true) => a.weighted / a.bases as f64,
                    (Stat::Min, true) => a.min as f64,
                    (Stat::Max, true) => a.max as f64,
                    (Stat::Sum, true) => a.sum,
                    (Stat::Std, true) => match a.bases {
                        0 => f64::NAN,
                        1 => self.deviations[k].sqrt(),
                        _ => (self.deviations[k] / (a.bases - 1) as f64).sqrt(),
                    },
                };
                *slot = v as f32;
            }
        }
    }

    fn finish_zoom(&self, layout: &BinLayout, stats: &[Stat], out: &mut [f32]) {
        let n = layout.count;
        for (row, stat) in out.chunks_exact_mut(n.max(1)).zip(stats) {
            for (k, slot) in row.iter_mut().enumerate().take(n) {
                let z = &self.zoom[k];
                let v = match stat {
                    Stat::Coverage => z.bases / layout.width(k) as f64,
                    _ if !z.any => f64::NAN,
                    Stat::Min => z.min as f64,
                    Stat::Max => z.max as f64,
                    _ if z.bases == 0.0 => f64::NAN,
                    Stat::Mean => z.sum / z.bases,
                    Stat::Sum => z.sum,
                    Stat::Std => {
                        // libBigWig's blockDev.
                        if z.bases <= 1.0 {
                            f64::NAN
                        } else {
                            let diff = (z.squares - z.sum * z.sum / z.bases) / (z.bases - 1.0);
                            if diff.abs() > 1e-8 { diff.sqrt() } else { 0.0 }
                        }
                    }
                };
                *slot = v as f32;
            }
        }
    }
}
