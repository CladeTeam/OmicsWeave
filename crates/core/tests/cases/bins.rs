//! Binned statistics and contact-grid geometry.

use std::panic::{AssertUnwindSafe, catch_unwind};

use proptest::prelude::*;

use omicsweave_core::{
    BinAccumulator, BinLayout, BinSpec, ContactQuery, Error, Locus, Region, Stat, Strand,
};

#[test]
fn bins_reject_reversed_converted_loci() {
    let locus = Locus::from(Region::new("c", 9, 3));
    assert!(matches!(
        BinLayout::new(&locus, &BinSpec::new(4, &[Stat::Mean])),
        Err(Error::InvalidRegion { .. })
    ));
}

#[test]
fn output_shape_and_std_input_selection_are_checked_before_writing() {
    let locus = Locus::from(Region::new("c", 0, 8));
    let spec = BinSpec::new(4, &[Stat::Mean]);
    let layout = BinLayout::new(&locus, &spec).unwrap();
    let mut acc = BinAccumulator::new();
    acc.reset(layout, spec.stats, false);
    acc.add(0, 8, 2.0);
    let mut bad = [123.0];
    assert!(catch_unwind(AssertUnwindSafe(|| acc.finish(spec.stats, &mut bad))).is_err());
    assert_eq!(bad, [123.0]);
    let mut out = [456.0; 2];
    assert!(catch_unwind(AssertUnwindSafe(|| acc.finish(&[Stat::Std], &mut out))).is_err());
    assert_eq!(out, [456.0; 2]);
    acc.finish(spec.stats, &mut out);
    assert_eq!(out, [2.0; 2]);
}

#[test]
fn repeated_std_queries_match_per_position_statistics() {
    let mut acc = BinAccumulator::new();
    for len in [33u32, 2, 0, 127, 9, 1] {
        for resolution in [1, 3, 16] {
            let locus = Locus::from(Region::new("c", 0, len));
            let spec = BinSpec::new(resolution, &[Stat::Mean, Stat::Std, Stat::Coverage]);
            let layout = BinLayout::new(&locus, &spec).unwrap();
            acc.reset(layout, spec.stats, false);
            let values: Vec<f32> = (0..len)
                .map(|i| (((i / 5) * 7) % 19) as f32 - 9.0)
                .collect();
            // Contiguous runs exercise splitting across bins without
            // changing the independent per-position reference below.
            for start in (0..len).step_by(5) {
                let end = (start + 5).min(len);
                acc.add(start, end, values[start as usize]);
            }
            let n = layout.count;
            let mut guarded = vec![1234.0f32; spec.output_len(len) + 2];
            let last = guarded.len() - 1;
            acc.finish(spec.stats, &mut guarded[1..last]);
            assert_eq!((guarded[0], guarded[last]), (1234.0, 1234.0));
            let out = &guarded[1..last];
            for k in 0..n {
                let (start, end) = layout.bounds(k);
                let data = &values[start as usize..end as usize];
                let mean = data.iter().map(|&v| v as f64).sum::<f64>() / data.len() as f64;
                let dev = data.iter().map(|&v| (v as f64 - mean).powi(2)).sum::<f64>();
                let std = if data.len() == 1 {
                    0.0
                } else {
                    (dev / (data.len() - 1) as f64).sqrt()
                };
                assert_eq!(out[k], mean as f32);
                assert_eq!(out[n + k], std as f32);
                assert_eq!(out[2 * n + k], 1.0);
            }
        }
    }
}

#[test]
fn bins_match_bins_of_individual_bases() {
    for resolution in 1..=9 {
        for start in 0..=20 {
            for end in 0..=20 {
                let region = Region::new("s", start, end);
                let mut expected: Vec<_> = (start..end).map(|base| base / resolution).collect();
                expected.dedup();
                assert_eq!(
                    ContactQuery::bins(&region, resolution).collect::<Vec<_>>(),
                    expected
                );
            }
        }
    }
}

#[test]
// The oracle uses widened arithmetic independently of production div_ceil.
#[allow(clippy::manual_div_ceil)]
fn maximum_coordinates_do_not_overflow_rounding() {
    for start in [0, 1, u32::MAX - 2, u32::MAX] {
        let region = Region::new("s", start, u32::MAX);
        assert!(region.check_within(u32::MAX).is_ok());
        for resolution in [1, 2, 3, 1 << 31, u32::MAX] {
            let actual = ContactQuery::bins(&region, resolution);
            let lo = u64::from(start) / u64::from(resolution);
            let hi = if start == u32::MAX {
                lo
            } else {
                (u64::from(u32::MAX) + u64::from(resolution) - 1) / u64::from(resolution)
            };
            assert_eq!(u64::from(actual.start), lo);
            assert_eq!(u64::from(actual.end), hi);
            let query = ContactQuery::new(region, region, resolution);
            let bins = usize::try_from(hi - lo).unwrap();
            assert_eq!(query.shape(), (bins, bins));
        }
    }
    assert!(matches!(
        Region::new("s", u32::MAX, 0).check_within(0),
        Err(Error::InvalidRegion { .. })
    ));
    assert!(matches!(
        Region::new("s", u32::MAX, u32::MAX).check_within(u32::MAX - 1),
        Err(Error::RegionOutOfRange { .. })
    ));
}

#[test]
fn zero_resolution_is_rejected_even_for_empty_regions() {
    for region in [Region::new("s", 1, 1), Region::new("s", 1, 2)] {
        let query = ContactQuery::new(region, region, 0);
        assert!(query.check().is_err());
        assert!(std::panic::catch_unwind(|| query.shape()).is_err());
    }
}

proptest! {
    /// A reverse locus' bins are the forward bins of the mirrored phase,
    /// reversed: summing position indices shows which positions each bin
    /// holds, in the locus' own coordinates.
    #[test]
    fn reverse_bins_cover_local_positions(
        len in 0u32..500,
        resolution in 1u32..60,
        phase_seed in any::<u32>(),
        reverse in any::<bool>(),
    ) {
        let phase = phase_seed % resolution;
        let strand = if reverse { Strand::Reverse } else { Strand::Forward };
        let locus = Locus::region(Region::new("c", 1000, 1000 + len), strand).unwrap();
        let spec = BinSpec { phase, ..BinSpec::new(resolution, &[Stat::Sum, Stat::Coverage]) };
        let layout = BinLayout::new(&locus, &spec).unwrap();
        let n = spec.count(len);
        prop_assert_eq!(layout.count, n);
        let mut acc = BinAccumulator::new();
        acc.reset(layout, spec.stats, false);
        // Forward position f carries the value of its local position.
        for f in 0..len {
            let local = if reverse { len - 1 - f } else { f };
            acc.add(f, f + 1, local as f32);
        }
        let mut out = vec![0f32; spec.output_len(len)];
        acc.finish(spec.stats, &mut out);
        for k in 0..n {
            let lo = (k as u32 * resolution).saturating_sub(phase);
            let hi = ((k as u32 + 1) * resolution - phase).min(len);
            let want: u32 = (lo..hi).sum();
            prop_assert_eq!(out[k], want as f32);
            prop_assert_eq!(out[n + k], 1.0);
        }
    }
}
