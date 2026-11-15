use std::panic::{AssertUnwindSafe, catch_unwind};

use omicsweave_core::{BinAccumulator, BinLayout, BinSpec, Error, Locus, Region, Stat};

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
