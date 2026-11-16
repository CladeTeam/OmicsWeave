//! Sequence and tile orientation.

use omicsweave_core::orient::{flip_tile, reverse_complement};

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
