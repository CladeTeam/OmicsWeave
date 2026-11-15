//! Bin statistics: bit-for-bit against libBigWig's `bwStats` on full bins,
//! and against per-base values on any locus (spans, strand, phase, partial
//! bins).

mod common;

use std::path::Path;

use common::libbigwig::{Reader, Writer};
use common::{Rng, same};
use omicsweave_bbi::BbiReader;
use omicsweave_core::{BinSource, BinSpec, Error, Locus, Region, Span, Stat, Strand, TrackRead};
use proptest::prelude::*;

/// A random bigWig of bedGraph, varStep and fixedStep sections, with zoom
/// levels; returns its chromosomes.
fn write_random(path: &Path, rng: &mut Rng) -> Vec<(String, u32)> {
    let chroms: Vec<(String, u32)> = (0..rng.range(1, 3))
        .map(|i| (format!("chr{}", i + 1), rng.range(5_000, 200_000) as u32))
        .collect();
    let refs: Vec<(&str, u32)> = chroms.iter().map(|(n, l)| (n.as_str(), *l)).collect();
    let mut w = Writer::create(path, rng.range(2, 10) as i32, &refs).expect("create");
    for (name, len) in &chroms {
        let mut pos = rng.range(0, 500) as u32;
        while pos < *len {
            let count = rng.range(1, 200);
            match rng.range(0, 2) {
                0 => {
                    let mut items = Vec::new();
                    for _ in 0..count {
                        if pos >= *len {
                            break;
                        }
                        let end = (pos + rng.range(1, 80) as u32).min(*len);
                        items.push((pos, end, rng.value()));
                        pos = end + rng.range(0, 40) as u32;
                    }
                    assert!(w.add_bedgraph(name, &items));
                }
                1 => {
                    let span = rng.range(1, 20) as u32;
                    let mut items = Vec::new();
                    for _ in 0..count {
                        if pos + span > *len {
                            pos = *len;
                            break;
                        }
                        items.push((pos, rng.value()));
                        pos += span + rng.range(0, 30) as u32;
                    }
                    if !items.is_empty() {
                        assert!(w.add_varstep(name, span, &items));
                    }
                }
                _ => {
                    let span = rng.range(1, 20) as u32;
                    let step = span + rng.range(0, 10) as u32;
                    if pos + span > *len {
                        break;
                    }
                    let k = (count as u32).min((*len - pos - span) / step + 1);
                    let values: Vec<f32> = (0..k).map(|_| rng.value()).collect();
                    assert!(w.add_fixedstep(name, pos, span, step, &values));
                    pos += k * step + rng.range(0, 30) as u32;
                }
            }
            pos += rng.range(0, 3_000) as u32;
        }
    }
    drop(w);
    chroms
}

/// libBigWig's `enum bwStatsType` value of a statistic.
fn stats_type(stat: Stat) -> i32 {
    match stat {
        Stat::Mean => 0,
        Stat::Std => 1,
        Stat::Max => 2,
        Stat::Min => 3,
        Stat::Coverage => 4,
        Stat::Sum => 5,
    }
}

fn bins(reader: &mut BbiReader<Vec<u8>>, locus: &Locus<'_>, spec: &BinSpec<'_>) -> Vec<f32> {
    let mut out = vec![f32::INFINITY; spec.output_len(locus.len())];
    reader.bins(locus, spec, &mut out).expect("bins");
    out
}

/// libBigWig reads an empty bin's coverage as NaN; ours is 0.
fn matches_libbigwig(stat: Stat, ours: f32, theirs: f64) -> bool {
    if stat == Stat::Coverage && theirs.is_nan() {
        return ours == 0.0;
    }
    same(ours, theirs as f32) && (ours.is_nan() || ours.to_bits() == (theirs as f32).to_bits())
}

fn compare_with_libbigwig(seed: u64, source: BinSource) -> Result<(), TestCaseError> {
    let mut rng = Rng::new(seed);
    let dir = tempfile::tempdir().unwrap();
    let path = dir.path().join("t.bw");
    let chroms = write_random(&path, &mut rng);
    let mut theirs = Reader::open(&path).expect("libBigWig open");
    let mut ours = BbiReader::open(std::fs::read(&path).unwrap()).unwrap();
    for _ in 0..12 {
        let (name, len) = &chroms[rng.range(0, chroms.len() as u64 - 1) as usize];
        let resolution = [1, 3, 64, 128, 1000, 2048, 10_000][rng.range(0, 6) as usize];
        let count = rng.range(1, 40) as u32;
        let width = (resolution * count).min(*len / resolution * resolution);
        if width == 0 {
            continue;
        }
        let count = width / resolution;
        let start = rng.range(0, (*len - width) as u64) as u32;
        let locus = Locus::from(Region::new(name, start, start + width));
        let spec = BinSpec {
            source,
            ..BinSpec::new(resolution, &Stat::ALL)
        };
        let got = bins(&mut ours, &locus, &spec);
        let zoomed = source == BinSource::Zoom && ours.zoom_level_for(resolution).is_some();
        for (row, &stat) in Stat::ALL.iter().enumerate() {
            // libBigWig's zoom sum truncates each record's weight to an
            // integer; ours weights by the fraction inside the bin.
            if zoomed && stat == Stat::Sum {
                continue;
            }
            let want = theirs
                .stats(
                    name,
                    start,
                    start + width,
                    count,
                    stats_type(stat),
                    source == BinSource::Exact,
                )
                .expect("libBigWig stats");
            let row = &got[row * count as usize..(row + 1) * count as usize];
            for (k, (&a, &b)) in row.iter().zip(&want).enumerate() {
                prop_assert!(
                    matches_libbigwig(stat, a, b),
                    "{} bin {k} of {name}:{start}+{count}x{resolution} ({source:?}): ours {a}, libBigWig {b}",
                    stat.as_str()
                );
            }
        }
    }
    Ok(())
}

/// Statistics of positions `[lo, hi)` of per-base values (NaN = no data).
fn brute(values: &[f32], stat: Stat) -> f64 {
    let data: Vec<f64> = values
        .iter()
        .filter(|v| !v.is_nan())
        .map(|&v| v as f64)
        .collect();
    if stat == Stat::Coverage {
        return data.len() as f64 / values.len() as f64;
    }
    if data.is_empty() {
        return f64::NAN;
    }
    let n = data.len() as f64;
    let mean = data.iter().sum::<f64>() / n;
    match stat {
        Stat::Mean => mean,
        Stat::Sum => data.iter().sum(),
        Stat::Min => data.iter().copied().fold(f64::INFINITY, f64::min),
        Stat::Max => data.iter().copied().fold(f64::NEG_INFINITY, f64::max),
        Stat::Std => {
            let dev: f64 = data.iter().map(|v| (v - mean) * (v - mean)).sum();
            if data.len() == 1 {
                0.0
            } else {
                (dev / (n - 1.0)).sqrt()
            }
        }
        Stat::Coverage => unreachable!(),
    }
}

/// Ours against the per-base reference. `slack` adds to the relative
/// tolerance: for the sum, libBigWig rounds each run's `len * value` to
/// `f32` before adding, an error of up to 2^-24 of every product, which
/// cancellation can leave larger than the result itself.
fn close(a: f32, b: f64, slack: f64) -> bool {
    if b.is_nan() {
        return a.is_nan();
    }
    (a as f64 - b).abs() <= 1e-4 * b.abs().max(1.0) + slack
}

/// The rounding allowance of a bin's sum: 2^-23 of the sum of magnitudes.
fn sum_slack(values: &[f32], stat: Stat) -> f64 {
    if stat != Stat::Sum {
        return 0.0;
    }
    values
        .iter()
        .filter(|v| !v.is_nan())
        .map(|v| v.abs() as f64)
        .sum::<f64>()
        * 2f64.powi(-23)
}

/// A random locus of 1–4 spans of `chrom`, possibly running past its end.
fn random_locus(rng: &mut Rng, len: u32) -> (Vec<Span>, Strand) {
    let mut spans = Vec::new();
    let mut at = rng.range(0, len as u64) as u32;
    for _ in 0..rng.range(1, 4) {
        let width = rng.range(0, 3_000) as u32;
        spans.push(Span::new(at, at + width));
        at += width + rng.range(0, 5_000) as u32;
    }
    let strand = [Strand::Forward, Strand::Reverse, Strand::Unstranded][rng.range(0, 2) as usize];
    (spans, strand)
}

proptest! {
    #![proptest_config(ProptestConfig::with_cases(24))]

    #[test]
    fn exact_bins_match_libbigwig(seed in any::<u64>()) {
        compare_with_libbigwig(seed, BinSource::Exact)?;
    }

    #[test]
    fn zoom_bins_match_libbigwig(seed in any::<u64>()) {
        compare_with_libbigwig(seed, BinSource::Zoom)?;
    }

    /// Bins of any locus summarize its oriented per-base values: bin `k`
    /// covers local positions `[k*r - p, (k+1)*r - p)` clipped to the locus.
    #[test]
    fn bins_summarize_oriented_values(seed in any::<u64>()) {
        let mut rng = Rng::new(seed);
        let dir = tempfile::tempdir().unwrap();
        let path = dir.path().join("t.bw");
        let chroms = write_random(&path, &mut rng);
        let mut ours = BbiReader::open(std::fs::read(&path).unwrap()).unwrap();
        for _ in 0..8 {
            let (name, len) = &chroms[rng.range(0, chroms.len() as u64 - 1) as usize];
            let (spans, strand) = random_locus(&mut rng, *len);
            let locus = Locus::new(name, &spans, strand).unwrap();
            let resolution = rng.range(1, 600) as u32;
            let phase = rng.range(0, resolution as u64 - 1) as u32;
            let spec = BinSpec { phase, ..BinSpec::new(resolution, &Stat::ALL) };
            let got = bins(&mut ours, &locus, &spec);
            let mut values = vec![0f32; locus.len() as usize];
            ours.values(&locus, &mut values).unwrap();
            let count = spec.count(locus.len());
            prop_assert_eq!(got.len(), count * Stat::ALL.len());
            for k in 0..count {
                let lo = (k * resolution as usize).saturating_sub(phase as usize);
                let hi = ((k + 1) * resolution as usize - phase as usize).min(values.len());
                for (row, &stat) in Stat::ALL.iter().enumerate() {
                    let a = got[row * count + k];
                    let b = brute(&values[lo..hi], stat);
                    let slack = sum_slack(&values[lo..hi], stat);
                    prop_assert!(close(a, b, slack), "{} bin {k}: ours {a}, per-base {b}", stat.as_str());
                }
            }
        }
    }

    /// Values of a multi-span locus are the spans' values concatenated, then
    /// reversed on the reverse strand.
    #[test]
    fn values_concatenate_and_orient(seed in any::<u64>()) {
        let mut rng = Rng::new(seed);
        let dir = tempfile::tempdir().unwrap();
        let path = dir.path().join("t.bw");
        let chroms = write_random(&path, &mut rng);
        let mut ours = BbiReader::open(std::fs::read(&path).unwrap()).unwrap();
        let mut theirs = Reader::open(&path).expect("libBigWig open");
        for _ in 0..8 {
            let (name, len) = &chroms[rng.range(0, chroms.len() as u64 - 1) as usize];
            let (spans, strand) = random_locus(&mut rng, *len);
            let locus = Locus::new(name, &spans, strand).unwrap();
            let mut got = vec![0f32; locus.len() as usize];
            ours.values(&locus, &mut got).unwrap();
            let mut want = Vec::new();
            for span in &spans {
                let end = span.end.min(*len);
                if span.start < end {
                    want.extend(theirs.values(name, span.start, end).unwrap());
                }
                want.extend(std::iter::repeat_n(f32::NAN, (span.end - end.max(span.start)) as usize));
            }
            if strand == Strand::Reverse {
                want.reverse();
            }
            prop_assert_eq!(got.len(), want.len());
            prop_assert!(got.iter().zip(&want).all(|(&a, &b)| same(a, b)));
        }
    }
}

#[test]
fn bin_specs_are_checked() {
    let dir = tempfile::tempdir().unwrap();
    let path = dir.path().join("t.bw");
    let mut rng = Rng::new(7);
    let chroms = write_random(&path, &mut rng);
    let mut reader = BbiReader::open(std::fs::read(&path).unwrap()).unwrap();
    let locus = Locus::from(Region::new(&chroms[0].0, 0, 1000));
    let mut out = vec![0f32; 10];
    for spec in [
        BinSpec::new(0, &[Stat::Mean]),
        BinSpec {
            phase: 100,
            ..BinSpec::new(100, &[Stat::Mean])
        },
        BinSpec::new(100, &[]),
        BinSpec::new(100, &[Stat::Mean, Stat::Mean]),
    ] {
        assert!(matches!(
            reader.bins(&locus, &spec, &mut out),
            Err(Error::Invalid { .. })
        ));
    }
    assert!(matches!(
        reader.bins(
            &locus,
            &BinSpec::new(100, &[Stat::Mean, Stat::Max]),
            &mut out
        ),
        Err(Error::OutputSize { need: 20, .. })
    ));
    // Partial bins at both ends: 1000 bases with phase 30 at 100 per bin.
    let spec = BinSpec {
        phase: 30,
        ..BinSpec::new(100, &[Stat::Coverage])
    };
    assert_eq!(spec.count(1000), 11);
    let mut out = vec![0f32; 11];
    reader.bins(&locus, &spec, &mut out).unwrap();
    assert!(out.iter().all(|c| (0.0..=1.0).contains(c)));
}
