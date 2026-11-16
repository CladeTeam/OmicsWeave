//! [`Records`]: the records a loaded region takes, decoded once from the
//! BCF and kept with what consensus building reads of each: position,
//! reference span, alleles and genotypes.

use omicsweave_core::{Error, Result};
use omicsweave_vcf::{BT_INT8, GT_MISSING_VALUE, GT_VECTOR_END, Record, VCF_REF};

/// Genotype values are kept as bcftools reads them: BCF `int8` bytes,
/// `(allele + 1) << 1 | phased`, with these two specials.
pub(crate) const GT_END: u8 = 0x81;
const GT_MISSING_SAMPLE: u8 = 0x80;

/// How a record's genotypes are held.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) enum Gt {
    /// No FORMAT/GT.
    Absent,
    /// One byte per value.
    Bytes,
    /// Non-`int8` storage or values, which bcftools does not read.
    Wide,
}

/// One record.
#[derive(Debug, Clone)]
pub(crate) struct Site {
    /// 0-based position of the first REF base.
    pub pos: i64,
    /// Reference span (htslib's `rlen`).
    pub rlen: i64,
    /// INFO/END, 1-based inclusive, when given.
    pub info_end: Option<i64>,
    /// A multi-base block whose alleles all are reference (gVCF).
    pub ref_block: bool,
    pub n_allele: u32,
    pub gt: Gt,
    /// Genotype values per sample.
    pub gt_width: u16,
    /// Start of the allele bytes in [`Records::allele_bytes`].
    bytes: usize,
    /// Start of the allele ends in [`Records::allele_ends`].
    ends: usize,
    /// Start of the genotype bytes in [`Records::genotypes`].
    genotypes: usize,
    /// The filter's error on this record, in [`Records::errors`]: bcftools
    /// reports it when it reaches the record, and so do queries.
    error: Option<u32>,
}

impl Site {
    /// The last reference base, 0-based inclusive.
    pub fn last(&self) -> i64 {
        self.pos + self.rlen - 1
    }
}

/// Which records a region takes (`--regions-overlap`).
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum Overlap {
    /// The record starts in the region (`0`, `pos`).
    Start,
    /// The record's reference span overlaps the region (`1`, `record`;
    /// bcftools' default).
    #[default]
    Record,
    /// The variant part overlaps: the span after the bases every allele
    /// shares with REF (`2`, `variant`).
    Variant,
}

impl Overlap {
    /// Whether a record the index returns for the 0-based inclusive range
    /// `[start, end]` is taken, as htslib's synced reader decides it from
    /// the decoded span.
    pub(crate) fn takes(self, rec: &Record, start: i64, end: i64) -> bool {
        match self {
            Overlap::Start => rec.pos >= start && rec.pos <= end,
            Overlap::Record => rec.pos <= end && rec.pos + rec.rlen > start,
            Overlap::Variant => {
                rec.pos <= end && rec.pos + rec.rlen > start && rec.pos + shared_prefix(rec) <= end
            }
        }
    }
}

/// `-i`/`-e`: keep only the records a bcftools filter expression includes,
/// or drop those it matches.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum RecordFilter<'a> {
    Include(&'a str),
    Exclude(&'a str),
}

/// How many leading reference bases every ALT shares with REF, at most the
/// reference span; the span itself for a record without ALT (htslib's
/// synced reader for `--regions-overlap 2`).
fn shared_prefix(rec: &Record) -> i64 {
    let reference = rec.allele(0);
    if rec.n_alleles() == 1 {
        return rec.rlen;
    }
    let mut shared = rec.rlen.min(reference.len() as i64);
    for allele in rec.alleles().skip(1) {
        if shared == 0 {
            break;
        }
        shared = reference
            .iter()
            .zip(allele)
            .take(shared as usize)
            .take_while(|(r, a)| r == a)
            .count() as i64;
    }
    shared
}

/// The records of one loaded region, in file order, reused across loads.
#[derive(Default)]
pub(crate) struct Records {
    pub sites: Vec<Site>,
    allele_bytes: Vec<u8>,
    /// Per record: each allele's end relative to its first allele byte.
    allele_ends: Vec<u32>,
    genotypes: Vec<u8>,
    errors: Vec<Error>,
}

impl Records {
    pub fn clear(&mut self) {
        self.sites.clear();
        self.allele_bytes.clear();
        self.allele_ends.clear();
        self.genotypes.clear();
        self.errors.clear();
    }

    /// Appends `rec` of a file with `n_samples` samples; `error` is the
    /// filter's error on it.
    pub fn push(&mut self, rec: &Record, n_samples: usize, error: Option<Error>) -> Result<()> {
        // Editing and mask overlaps use the exclusive reference end,
        // including for gVCF claims that extend beyond the query.
        if rec.pos < -1 || rec.rlen < 0 || rec.pos.checked_add(rec.rlen).is_none() {
            return Err(Error::invalid(
                "vcf record",
                "reference span is outside signed 64-bit coordinates",
            ));
        }
        if rec.n_alleles() == 0 {
            return Err(Error::invalid("vcf record", "no reference allele"));
        }
        let bytes = self.allele_bytes.len();
        let ends = self.allele_ends.len();
        for allele in rec.alleles() {
            self.allele_bytes.extend_from_slice(allele);
            let end = u32::try_from(self.allele_bytes.len() - bytes)
                .map_err(|_| Error::invalid("vcf record", "alleles longer than 4 GiB"))?;
            self.allele_ends.push(end);
        }
        let genotypes = self.genotypes.len();
        let (gt, gt_width) = match rec.genotypes() {
            None => (Gt::Absent, 0),
            Some(g) => {
                if g.width == 0 || g.width.checked_mul(n_samples) != Some(g.values.len()) {
                    return Err(Error::invalid("vcf record", "invalid genotype layout"));
                }
                let width = u16::try_from(g.width)
                    .map_err(|_| Error::invalid("vcf record", "ploidy above 65535"))?;
                let gt = if rec.genotype_type().is_some_and(|ty| ty != BT_INT8) {
                    Gt::Wide
                } else {
                    self.push_genotypes(g.values)
                };
                (gt, width)
            }
        };
        let error = match error {
            None => None,
            Some(e) => {
                let i = u32::try_from(self.errors.len())
                    .map_err(|_| Error::invalid("vcf filter", "saved error index exceeds u32"))?;
                self.errors.push(e);
                Some(i)
            }
        };
        self.sites.push(Site {
            pos: rec.pos,
            rlen: rec.rlen,
            info_end: rec.end,
            ref_block: rec.rlen > 1 && rec.variant_types() == VCF_REF,
            n_allele: rec.n_alleles() as u32,
            gt,
            gt_width,
            bytes,
            ends,
            genotypes,
            error,
        });
        Ok(())
    }

    fn push_genotypes(&mut self, values: &[i32]) -> Gt {
        let start = self.genotypes.len();
        for &v in values {
            let byte = match v {
                GT_VECTOR_END => GT_END,
                GT_MISSING_VALUE => GT_MISSING_SAMPLE,
                0..=127 => v as u8,
                _ => {
                    self.genotypes.truncate(start);
                    return Gt::Wide;
                }
            };
            self.genotypes.push(byte);
        }
        Gt::Bytes
    }

    /// The filter's error on this record.
    pub fn error(&self, site: &Site) -> Option<&Error> {
        site.error.map(|i| &self.errors[i as usize])
    }

    /// Allele `i` (0 is REF).
    pub fn allele(&self, site: &Site, i: usize) -> &[u8] {
        let begin = if i == 0 {
            0
        } else {
            self.allele_ends[site.ends + i - 1] as usize
        };
        let end = self.allele_ends[site.ends + i] as usize;
        &self.allele_bytes[site.bytes + begin..site.bytes + end]
    }

    /// The genotype values of `sample` (`gt_width` of them); `None` without
    /// genotypes in `int8`.
    pub fn genotype(&self, site: &Site, sample: usize) -> Option<&[u8]> {
        match site.gt {
            Gt::Absent | Gt::Wide => None,
            Gt::Bytes => {
                let width = site.gt_width as usize;
                let first = site.genotypes + sample * width;
                Some(&self.genotypes[first..first + width])
            }
        }
    }
}
