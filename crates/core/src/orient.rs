//! Turning forward-decoded dense output to a key's orientation, in place:
//! readers decode 5' to 3' of the forward strand into the caller's buffer
//! and call these when the key is on the reverse strand.

use crate::locus::Strand;

const fn complement_table() -> [u8; 256] {
    let mut t = [0u8; 256];
    let mut i = 0;
    while i < 256 {
        t[i] = i as u8;
        i += 1;
    }
    // IUPAC: each code maps to the code of the complementary base set.
    // S, W and N are their own complements; U complements like T.
    let pairs: [(u8, u8); 7] = [
        (b'A', b'T'),
        (b'C', b'G'),
        (b'R', b'Y'),
        (b'K', b'M'),
        (b'B', b'V'),
        (b'D', b'H'),
        (b'U', b'A'),
    ];
    let mut p = 0;
    while p < pairs.len() {
        let (a, b) = pairs[p];
        t[a as usize] = b;
        t[(a | 0x20) as usize] = b | 0x20;
        if a != b'U' {
            t[b as usize] = a;
            t[(b | 0x20) as usize] = a | 0x20;
        }
        p += 1;
    }
    t
}

static COMPLEMENT: [u8; 256] = complement_table();

/// The complement of an ASCII base: IUPAC codes map to the code of the
/// complementary set, case is kept (soft-masking survives), `U` reads as
/// `T`, and every other byte is unchanged.
#[inline]
pub fn complement(base: u8) -> u8 {
    COMPLEMENT[base as usize]
}

/// Reverses and complements `bases` in place.
pub fn reverse_complement(bases: &mut [u8]) {
    let n = bases.len();
    for i in 0..n / 2 {
        let (a, b) = (bases[i], bases[n - 1 - i]);
        bases[i] = complement(b);
        bases[n - 1 - i] = complement(a);
    }
    if n % 2 == 1 {
        bases[n / 2] = complement(bases[n / 2]);
    }
}

/// Puts forward-decoded bases in `strand`'s orientation.
pub fn orient_bases(strand: Strand, bases: &mut [u8]) {
    if strand.is_reverse() {
        reverse_complement(bases);
    }
}

/// Puts forward-decoded per-position values in `strand`'s orientation.
pub fn orient_values<T>(strand: Strand, values: &mut [T]) {
    if strand.is_reverse() {
        values.reverse();
    }
}

/// Flips a row-major `rows × cols` tile in place: row order when
/// `flip_rows`, column order within each row when `flip_cols`.
///
/// # Panics
///
/// Panics if `tile.len() != rows * cols`, or the shape exceeds `usize`.
pub fn flip_tile<T>(tile: &mut [T], rows: usize, cols: usize, flip_rows: bool, flip_cols: bool) {
    assert_eq!(
        tile.len(),
        rows.checked_mul(cols).expect("tile shape exceeds usize"),
        "tile shape"
    );
    if cols == 0 {
        return;
    }
    if flip_rows && flip_cols {
        // Flipping both axes is the flat tile's reversal, so each pair
        // is swapped only once.
        tile.reverse();
        return;
    }
    if flip_rows {
        for r in 0..rows / 2 {
            let (top, bottom) = tile.split_at_mut((rows - 1 - r) * cols);
            top[r * cols..(r + 1) * cols].swap_with_slice(&mut bottom[..cols]);
        }
    }
    if flip_cols {
        for row in tile.chunks_exact_mut(cols) {
            row.reverse();
        }
    }
}
