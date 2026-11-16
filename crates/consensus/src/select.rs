//! Which allele of a record a sequence takes, as consensus.c chooses it
//! from `-s`, `-H`, `-I` and the genotypes (`apply_variant`, lines
//! 602-758).

use omicsweave_core::{Error, Result};
use smallvec::SmallVec;

use crate::iupac;
use crate::records::{GT_END, Gt, Records, Site};

const PICK_REF: u8 = 1;
const PICK_ALT: u8 = 2;
const PICK_LONG: u8 = 4;
const PICK_SHORT: u8 = 8;
const PICK_IUPAC: u8 = 16;
// Readable genotypes contain values 0..=127 plus the int8 missing/vector-end
// sentinels. The unsigned missing sentinel can address allele 63; larger
// allele numbers require a wide genotype and are rejected before selection.
const BYTE_GT_ALLELES: usize = 64;

/// A parsed `-H`: pick rules and a 1-based haplotype (0 for none).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) struct Haplotype {
    pick: u8,
    index: usize,
}

impl Haplotype {
    /// `1`, `2`, ... ; `R`, `A`, `L`/`LR`, `S`/`SR`, `LA`, `SA`, `I`;
    /// `<N>pIu`. Case-insensitive.
    pub fn parse(s: &str) -> Option<Haplotype> {
        let pick = |pick| Some(Haplotype { pick, index: 0 });
        match s.as_bytes() {
            [b'R' | b'r'] => return pick(PICK_REF),
            [b'A' | b'a'] => return pick(PICK_ALT),
            [b'L' | b'l'] | [b'L' | b'l', b'R' | b'r'] => return pick(PICK_LONG | PICK_REF),
            [b'S' | b's'] | [b'S' | b's', b'R' | b'r'] => return pick(PICK_SHORT | PICK_REF),
            [b'L' | b'l', b'A' | b'a'] => return pick(PICK_LONG | PICK_ALT),
            [b'S' | b's', b'A' | b'a'] => return pick(PICK_SHORT | PICK_ALT),
            [b'I' | b'i'] => return pick(PICK_IUPAC),
            _ => {}
        }
        let digits = s.bytes().take_while(u8::is_ascii_digit).count();
        let index: usize = s[..digits].parse().ok().filter(|&n| n > 0)?;
        let pick = match &s[digits..] {
            "" => 0,
            rest if rest.eq_ignore_ascii_case("piu") => PICK_IUPAC,
            _ => return None,
        };
        Some(Haplotype { pick, index })
    }
}

/// How a query chooses alleles, resolved against the file's samples.
#[derive(Debug, Clone)]
pub(crate) enum Mode {
    /// Every first ALT (`-s -`).
    AllAlt,
    /// IUPAC codes of all plain alleles (`-s - -I`).
    IupacRefAlt,
    /// IUPAC codes of these samples' genotypes (`-s a,b`).
    IupacSamples(SmallVec<[usize; 4]>),
    /// One sample with `-H` (`iupac`: `-I` given as well).
    Sample {
        sample: usize,
        haplotype: Haplotype,
        iupac: bool,
    },
}

/// The choice for one record.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) enum Choice {
    /// Leave the record out.
    Skip,
    /// The sample's own REF: nothing to apply, but with `-a` its span
    /// counts as covered.
    Ref,
    /// A missing genotype, written with `-M`.
    Missing,
    /// Apply this allele, whose bytes are in the scratch buffer (they
    /// differ from the record's when IUPAC codes were mixed in).
    Allele(usize),
}

/// Scratch reused across records.
#[derive(Default)]
pub(crate) struct Scratch {
    /// The chosen allele's bytes.
    pub allele: Vec<u8>,
    selected: Vec<bool>,
}

fn is_missing(v: u8) -> bool {
    v >> 1 == 0
}

fn is_phased(v: u8) -> bool {
    v & 1 != 0
}

/// The allele of a value read unsigned (consensus.c's `uint8_t *ptr`).
fn allele_of(v: u8) -> i64 {
    (v >> 1) as i64 - 1
}

fn invalid(site: &Site, what: &str) -> Error {
    Error::invalid("consensus", format!("{what} at position {}", site.pos + 1))
}

/// Chooses the allele of `site` under `mode`; `missing` is `-M` given.
pub(crate) fn choose(
    records: &Records,
    site: &Site,
    mode: &Mode,
    missing: bool,
    scratch: &mut Scratch,
) -> Result<Choice> {
    let n_allele = site.n_allele as usize;
    let allele = |i: usize| records.allele(site, i);
    let ialt = match mode {
        Mode::AllAlt => 1,
        Mode::IupacRefAlt if n_allele > 1 => {
            iupac::mix_ref_alt(n_allele, allele, &mut scratch.allele);
            return Ok(Choice::Allele(1));
        }
        Mode::IupacRefAlt => 1,
        Mode::IupacSamples(samples) => {
            if !readable_gt(site)? {
                return Ok(Choice::Skip);
            }
            scratch.selected.clear();
            scratch
                .selected
                .resize(n_allele.min(BYTE_GT_ALLELES), false);
            let mut any = false;
            for &s in samples {
                let gt = records.genotype(site, s).unwrap_or_default();
                any |= add_alleles(site, gt, &mut scratch.selected)?;
            }
            if !any {
                return Ok(if missing {
                    Choice::Missing
                } else {
                    Choice::Skip
                });
            }
            let selected = &scratch.selected;
            return Ok(
                match iupac::mix(selected.len(), allele, |i| selected[i], &mut scratch.allele) {
                    Some(i) => {
                        if n_allele == 1 {
                            Choice::Ref
                        } else {
                            Choice::Allele(i)
                        }
                    }
                    None => Choice::Missing,
                },
            );
        }
        Mode::Sample {
            sample,
            haplotype,
            iupac,
        } => {
            if !readable_gt(site)? {
                return Ok(Choice::Skip);
            }
            let gt = records.genotype(site, *sample).unwrap_or_default();
            let width = gt.len();
            let use_iupac = if haplotype.pick == PICK_IUPAC {
                haplotype.index == 0 || (!is_phased(gt[0]) && !is_phased(gt[width - 1]))
            } else {
                *iupac
            };
            let ialt = if use_iupac {
                scratch.selected.clear();
                scratch
                    .selected
                    .resize(n_allele.min(BYTE_GT_ALLELES), false);
                if !add_alleles(site, gt, &mut scratch.selected)? {
                    return Ok(if missing {
                        Choice::Missing
                    } else {
                        Choice::Skip
                    });
                }
                let selected = &scratch.selected;
                match iupac::mix(selected.len(), allele, |i| selected[i], &mut scratch.allele) {
                    Some(i) => i as i64,
                    None => -1,
                }
            } else if haplotype.index == 0 {
                match pick_one(site, records, gt, haplotype.pick, missing)? {
                    Some(i) => i,
                    None => return Ok(Choice::Skip),
                }
            } else {
                match use_hap(gt, haplotype.index, missing) {
                    Some(i) => i,
                    None => return Ok(Choice::Skip),
                }
            };
            match ialt {
                0 => return Ok(Choice::Ref),
                -1 => return Ok(Choice::Missing),
                i if i < 0 || i >= n_allele as i64 => {
                    return Err(invalid(site, "genotype allele beyond the ALT alleles"));
                }
                _ => {}
            }
            if use_iupac {
                return Ok(Choice::Allele(ialt as usize));
            }
            ialt
        }
    };
    if n_allele == 1 {
        return Ok(Choice::Ref);
    }
    scratch.allele.clear();
    scratch.allele.extend_from_slice(allele(ialt as usize));
    Ok(Choice::Allele(ialt as usize))
}

/// Whether the record has genotypes bcftools reads; an error if the storage
/// is not `int8` or its values do not fit it.
fn readable_gt(site: &Site) -> Result<bool> {
    match site.gt {
        Gt::Absent => Ok(false),
        Gt::Wide => Err(invalid(site, "genotypes not represented with int8")),
        Gt::Bytes => Ok(true),
    }
}

/// Marks the alleles of one genotype (`iupac_add_gt`); whether any was.
fn add_alleles(site: &Site, gt: &[u8], selected: &mut [bool]) -> Result<bool> {
    let mut any = false;
    for &v in gt {
        if is_missing(v) {
            continue;
        }
        if v == GT_END {
            break;
        }
        let a = allele_of(v) as usize;
        if a >= selected.len() {
            return Err(invalid(site, "too few ALT alleles"));
        }
        selected[a] = true;
        any = true;
    }
    Ok(any)
}

/// `-H N` (consensus.c 640-669): the N-th value; `None` to skip, -1 for a
/// missing allele written with `-M`. Values are read signed here.
fn use_hap(gt: &[u8], index: usize, missing: bool) -> Option<i64> {
    let missing_or_skip = missing.then_some(-1);
    if index > gt.len() {
        if is_missing(gt[gt.len() - 1]) || is_missing(gt[0]) {
            return missing_or_skip;
        }
        return None;
    }
    let v = gt[index - 1];
    if is_missing(v) || v == GT_END {
        return missing_or_skip;
    }
    Some(((v as i8) >> 1) as i64 - 1)
}

/// `-H R/A/LR/LA/SR/SA` (consensus.c 679-721): the allele of a homozygous
/// genotype; for a heterozygous one, the pick rules over its alleles. The
/// homozygosity scan leaves the choice at the first allele that differs,
/// and the rules replace it only where they apply.
fn pick_one(
    site: &Site,
    records: &Records,
    gt: &[u8],
    pick: u8,
    missing: bool,
) -> Result<Option<i64>> {
    let mut ialt: i64 = 1;
    let mut hom = true;
    for (i, &v) in gt.iter().enumerate() {
        if is_missing(v) {
            if !missing {
                return Ok(None);
            }
            ialt = -1;
            break;
        }
        if v == GT_END {
            break;
        }
        ialt = allele_of(v);
        if i > 0 && ialt != allele_of(gt[i - 1]) {
            hom = false;
            break;
        }
    }
    if hom {
        return Ok(Some(ialt));
    }
    let n_allele = site.n_allele as i64;
    let mut prev_len = 0;
    for (i, &v) in gt.iter().enumerate() {
        if v == GT_END {
            break;
        }
        let jalt = allele_of(v);
        if jalt >= n_allele {
            return Err(invalid(site, "too few ALT alleles"));
        }
        if pick & (PICK_LONG | PICK_SHORT) != 0 {
            let len = match jalt {
                0 => site.rlen,
                j if j > 0 => records.allele(site, j as usize).len() as i64,
                _ => 0,
            };
            let take = if i == 0 {
                true
            } else if len == prev_len {
                (pick & PICK_REF != 0 && jalt == 0) || (pick & PICK_ALT != 0 && ialt == 0)
            } else {
                (pick & PICK_LONG != 0 && len > prev_len)
                    || (pick & PICK_SHORT != 0 && len < prev_len)
            };
            if take {
                ialt = jalt;
                prev_len = len;
            }
        } else if (pick & PICK_REF != 0 && jalt == 0) || (pick & PICK_ALT != 0 && ialt == 0) {
            ialt = jalt;
        }
    }
    Ok(Some(ialt))
}
