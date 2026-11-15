//! IUPAC ambiguity codes of several alleles (bcftools.h `iupac2bitmask`,
//! `bitmask2iupac`; consensus.c `iupac_set_allele`).

const A: u8 = 1;
const C: u8 = 2;
const G: u8 = 4;
const T: u8 = 8;

/// The bases an IUPAC code stands for; `None` for any other character.
fn bases(code: u8) -> Option<u8> {
    Some(match code.to_ascii_uppercase() {
        b'A' => A,
        b'C' => C,
        b'G' => G,
        b'T' => T,
        b'M' => A | C,
        b'R' => A | G,
        b'W' => A | T,
        b'S' => C | G,
        b'Y' => C | T,
        b'K' => G | T,
        b'V' => A | C | G,
        b'H' => A | C | T,
        b'D' => A | G | T,
        b'B' => C | G | T,
        b'N' => A | C | G | T,
        _ => return None,
    })
}

fn code(bases: u8) -> u8 {
    b".ACMGRSVTWYHKDBN"[bases as usize]
}

/// Whether every character of `allele` is an IUPAC code (not symbolic, a
/// breakpoint or `*`).
fn plain(allele: &[u8]) -> bool {
    allele.iter().all(|&b| bases(b).is_some())
}

/// Combines the plain alleles among `selected` into IUPAC codes, position
/// by position, written over the longest selected ALT (the first of equal
/// length) into `out`. Returns that ALT, or else the first selected ALT, or
/// else REF, unchanged in `out`; `None` when nothing is selected.
pub(crate) fn mix<'a>(
    n_allele: usize,
    allele: impl Fn(usize) -> &'a [u8],
    selected: impl Fn(usize) -> bool,
    out: &mut Vec<u8>,
) -> Option<usize> {
    // Alleles are borrowed from the loaded records, so the reusable output
    // can hold their bitmasks until the chosen ALT length is known.
    out.clear();
    let mut longest = None;
    let mut fallback = None;
    for i in (0..n_allele).filter(|&i| selected(i)) {
        if fallback.is_none_or(|f| f == 0) {
            fallback = Some(i);
        }
        let a = allele(i);
        if !plain(a) {
            continue;
        }
        if a.len() > out.len() {
            out.resize(a.len(), 0);
        }
        if i > 0 && longest.is_none_or(|l: usize| a.len() > allele(l).len()) {
            longest = Some(i);
        }
        for (m, &b) in out.iter_mut().zip(a) {
            *m |= bases(b).unwrap_or(0);
        }
    }
    if let Some(l) = longest {
        out.truncate(allele(l).len());
        for m in out.iter_mut() {
            *m = code(*m);
        }
        return Some(l);
    }
    out.clear();
    let f = fallback?;
    out.extend_from_slice(allele(f));
    Some(f)
}

/// `-I` without samples: the codes of every plain allele, written over the
/// longest plain ALT. bcftools then applies ALT 1 (consensus.c shadows the
/// chosen index), so the codes show only when the longest is ALT 1; `out`
/// holds ALT 1 as applied.
pub(crate) fn mix_ref_alt<'a>(
    n_allele: usize,
    allele: impl Fn(usize) -> &'a [u8],
    out: &mut Vec<u8>,
) {
    let chosen = mix(n_allele, &allele, |_| true, out);
    if chosen != Some(1) {
        out.clear();
        out.extend_from_slice(allele(1));
    }
}
