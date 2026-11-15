//! Loci, orientation helpers and the bin layout.

use omicsweave_core::orient::{flip_tile, reverse_complement};
use omicsweave_core::{
    BinAccumulator, BinLayout, BinSpec, Error, Locus, Region, Span, Stat, Strand,
};
use proptest::prelude::*;

#[test]
fn loci_validate_and_describe_their_spans() {
    let spans = [
        Span::new(10, 20),
        Span::new(20, 25),
        Span::new(40, 40),
        Span::new(50, 60),
    ];
    let locus = Locus::new("c", &spans, Strand::Reverse).unwrap();
    assert_eq!(locus.len(), 25);
    assert!(locus.is_reverse());
    assert_eq!(locus.hull(), Region::new("c", 10, 60));
    let pieces: Vec<_> = locus.pieces().collect();
    assert_eq!(
        pieces,
        [
            (0, Region::new("c", 10, 20)),
            (10, Region::new("c", 20, 25)),
            (15, Region::new("c", 50, 60))
        ]
    );
    assert!(matches!(
        locus.check_within(59),
        Err(Error::RegionOutOfRange { .. })
    ));
    assert!(locus.check_within(60).is_ok());
    assert!(matches!(
        Locus::new(
            "c",
            &[Span::new(10, 20), Span::new(19, 30)],
            Strand::Forward
        ),
        Err(Error::Invalid { .. })
    ));
    let empty = Locus::new("c", &[], Strand::Forward).unwrap();
    assert!(empty.is_empty() && empty.hull().is_empty());
    assert!(matches!(
        Locus::from(Region::new("c", 5, 3)).check(),
        Err(Error::InvalidRegion { .. })
    ));
    assert!(matches!(
        Locus::from(Region::new("c", 5, 3)).check_within(10),
        Err(Error::InvalidRegion { .. })
    ));
}

#[test]
fn within_checks_converted_keys_and_extreme_spans() {
    for region in [Region::new("c", 7, 3), Region::new("c", u32::MAX, 2)] {
        assert!(matches!(
            Locus::from(region).check_within(0),
            Err(Error::InvalidRegion { .. })
        ));
    }
    let empty = Locus::new("c", &[], Strand::Reverse).unwrap();
    assert!(empty.check_within(0).is_ok());
    let at_end = Locus::from(Region::new("c", u32::MAX, u32::MAX));
    assert!(at_end.check_within(u32::MAX).is_ok());
    assert!(matches!(
        at_end.check_within(u32::MAX - 1),
        Err(Error::RegionOutOfRange { .. })
    ));
    let spans = [
        Span::new(0, 5),
        Span::new(10, 20),
        Span::new(u32::MAX, u32::MAX),
    ];
    let spliced = Locus::new("c", &spans, Strand::Forward).unwrap();
    assert!(spliced.check().is_ok());
    assert!(spliced.check_within(u32::MAX).is_ok());
    assert!(matches!(
        spliced.check_within(20),
        Err(Error::RegionOutOfRange { .. })
    ));
    let whole_span = [Span::new(0, u32::MAX)];
    let whole = Locus::new("c", &whole_span, Strand::Forward).unwrap();
    assert_eq!(whole.len(), u32::MAX);
    assert!(whole.check_within(u32::MAX).is_ok());
}

#[test]
fn reverse_complement_handles_odd_lengths() {
    let mut odd = *b"ACGTN";
    reverse_complement(&mut odd);
    assert_eq!(&odd, b"NACGT");
    let mut even = *b"aaCG";
    reverse_complement(&mut even);
    assert_eq!(&even, b"CGtt");
}

#[test]
fn tiles_flip_by_axis() {
    let tile: Vec<u32> = (0..12).collect();
    let mut t = tile.clone();
    flip_tile(&mut t, 3, 4, true, false);
    assert_eq!(t, [8, 9, 10, 11, 4, 5, 6, 7, 0, 1, 2, 3]);
    let mut t = tile.clone();
    flip_tile(&mut t, 3, 4, false, true);
    assert_eq!(t, [3, 2, 1, 0, 7, 6, 5, 4, 11, 10, 9, 8]);
    let mut t = tile.clone();
    flip_tile(&mut t, 3, 4, true, true);
    assert_eq!(t, tile.iter().rev().copied().collect::<Vec<_>>());
}

#[test]
fn tile_axes_match_coordinate_mapping_for_empty_odd_and_even_shapes() {
    for rows in 0..12 {
        for cols in 0..12 {
            for flip_rows in [false, true] {
                for flip_cols in [false, true] {
                    let mut tile: Vec<usize> = (0..rows * cols).collect();
                    flip_tile(&mut tile, rows, cols, flip_rows, flip_cols);
                    for r in 0..rows {
                        for c in 0..cols {
                            let from_row = if flip_rows { rows - 1 - r } else { r };
                            let from_col = if flip_cols { cols - 1 - c } else { c };
                            assert_eq!(tile[r * cols + c], from_row * cols + from_col);
                        }
                    }
                }
            }
        }
    }
}

#[test]
fn tile_shape_multiplication_cannot_wrap_into_a_valid_buffer() {
    assert!(
        std::panic::catch_unwind(|| {
            flip_tile(&mut [0u8], usize::MAX, usize::MAX, false, false);
        })
        .is_err()
    );
    flip_tile::<u8>(&mut [], 0, usize::MAX, true, true);
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
