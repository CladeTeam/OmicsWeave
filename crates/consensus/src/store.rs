//! [`VcfStore`]: the records of a variant file, read once and held resident,
//! sorted by position within each contig, with what consensus building
//! reads of each: position, reference span, alleles and genotypes.

use std::collections::HashMap;
use std::io::Read;

use omicsweave_core::{Error, Result};
use omicsweave_vcf::{Filter, GT_MISSING_VALUE, GT_VECTOR_END, Record, VCF_REF, VariantReader};

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
    /// Two bits per value (value - 2): every value is allele 0 or 1, with
    /// or without phase, which covers most records in a small fraction of
    /// the bytes.
    Packed,
    /// Values beyond `int8`, which bcftools does not read.
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
    /// Start of the allele bytes in [`VcfStore::allele_bytes`].
    bytes: usize,
    /// Start of the allele ends in [`VcfStore::allele_ends`].
    ends: usize,
    /// Start of the genotype bytes in [`VcfStore::genotypes`].
    genotypes: usize,
    /// The filter's error on this record, in [`VcfStore::errors`]: bcftools
    /// reports it when it reaches the record, and so do queries.
    error: Option<u32>,
}

impl Site {
    /// The last reference base, 0-based inclusive.
    pub fn last(&self) -> i64 {
        self.pos + self.rlen - 1
    }
}

/// The records of one contig, sorted by position.
#[derive(Default)]
struct Contig {
    sites: Vec<Site>,
    /// Empty when record ends are nondecreasing. Otherwise the sites form
    /// an implicit balanced search tree: each slice's middle is its root,
    /// and `subtree_last[middle]` is the furthest last base of that slice.
    subtree_last: Vec<i64>,
}

impl Contig {
    fn index(&mut self) {
        // Stable, so records at one position keep their file order.
        self.sites.sort_by_key(|s| s.pos);
        if self.sites.windows(2).all(|s| s[0].last() <= s[1].last()) {
            return;
        }
        self.subtree_last.resize(self.sites.len(), i64::MIN);
        fn build(sites: &[Site], last: &mut [i64]) -> i64 {
            if sites.is_empty() {
                return i64::MIN;
            }
            let mid = sites.len() / 2;
            let (left, root_and_right) = last.split_at_mut(mid);
            let (root, right) = root_and_right.split_first_mut().unwrap();
            *root = sites[mid]
                .last()
                .max(build(&sites[..mid], left))
                .max(build(&sites[mid + 1..], right));
            *root
        }
        build(&self.sites, &mut self.subtree_last);
    }

    /// Enumerates record-span overlaps in stable position order. A long
    /// record cannot force scans of disjoint records in other subtrees.
    fn overlapping(&self, start: i64, end: i64, out: &mut Vec<usize>) {
        if self.subtree_last.is_empty() {
            let first = self.sites.partition_point(|s| s.last() < start);
            let after = self.sites.partition_point(|s| s.pos <= end);
            out.extend(first..after);
            return;
        }
        fn visit(
            sites: &[Site],
            last: &[i64],
            offset: usize,
            start: i64,
            end: i64,
            out: &mut Vec<usize>,
        ) {
            if sites.is_empty() {
                return;
            }
            let mid = sites.len() / 2;
            if last[mid] < start || sites[0].pos > end {
                return;
            }
            visit(&sites[..mid], &last[..mid], offset, start, end, out);
            let site = &sites[mid];
            if site.pos > end {
                return;
            }
            if site.last() >= start {
                out.push(offset + mid);
            }
            visit(
                &sites[mid + 1..],
                &last[mid + 1..],
                offset + mid + 1,
                start,
                end,
                out,
            );
        }
        visit(&self.sites, &self.subtree_last, 0, start, end, out);
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

/// `-i`/`-e`: keep only the records a bcftools filter expression includes,
/// or drop those it matches.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum RecordFilter<'a> {
    Include(&'a str),
    Exclude(&'a str),
}

/// The records of a VCF or BCF file, resident and shared: one store serves
/// every [`crate::Consensus`] built from it, on any thread.
pub struct VcfStore {
    samples: Vec<String>,
    sample_ids: HashMap<String, usize>,
    contig_ids: HashMap<String, usize>,
    contigs: Vec<Contig>,
    allele_bytes: Vec<u8>,
    /// Per record: its shared variant prefix, then each allele's end
    /// relative to its first allele byte. `Site::ends` points past the
    /// prefix, so prefix caching does not enlarge every `Site`.
    allele_ends: Vec<u32>,
    genotypes: Vec<u8>,
    errors: Vec<Error>,
    skipped: usize,
}

impl VcfStore {
    /// Reads a VCF or BCF file (plain, gzip or BGZF) whole, keeping the
    /// records `filter` passes. Malformed VCF lines are skipped, as
    /// htslib-based readers skip them, and counted. An expression that does
    /// not compile is an error here.
    pub fn read(source: impl Read, filter: Option<RecordFilter<'_>>) -> Result<VcfStore> {
        let mut reader = VariantReader::new(source)?;
        let mut filter = match filter {
            None => None,
            Some(RecordFilter::Include(e)) => Some((Filter::new(reader.header(), e)?, true)),
            Some(RecordFilter::Exclude(e)) => Some((Filter::new(reader.header(), e)?, false)),
        };
        let samples = reader.header().samples.clone();
        let mut store = VcfStore {
            sample_ids: samples
                .iter()
                .enumerate()
                .map(|(i, s)| (s.clone(), i))
                .collect(),
            samples,
            contig_ids: HashMap::new(),
            contigs: Vec::new(),
            allele_bytes: Vec::new(),
            allele_ends: Vec::new(),
            genotypes: Vec::new(),
            errors: Vec::new(),
            skipped: 0,
        };
        let mut rec = Record::default();
        loop {
            match reader.read(&mut rec) {
                Ok(false) => break,
                Ok(true) => {}
                Err(Error::Format { .. }) => {
                    store.skipped += 1;
                    continue;
                }
                Err(e) => return Err(e),
            }
            let mut error = None;
            if let Some((filter, include)) = &mut filter {
                match filter.test(reader.header(), &rec) {
                    Ok(pass) if pass != *include => continue,
                    Ok(_) => {}
                    Err(e) => {
                        error = Some(u32::try_from(store.errors.len()).map_err(|_| {
                            Error::invalid("vcf filter", "saved error index exceeds u32")
                        })?);
                        store.errors.push(e);
                    }
                }
            }
            let name = &reader.header().contigs[rec.rid as usize].name;
            let contig = match store.contig_ids.get(name) {
                Some(&id) => id,
                None => {
                    store.contig_ids.insert(name.clone(), store.contigs.len());
                    store.contigs.push(Contig::default());
                    store.contigs.len() - 1
                }
            };
            let mut site = store.push(&rec)?;
            site.error = error;
            store.contigs[contig].sites.push(site);
        }
        for contig in &mut store.contigs {
            contig.index();
        }
        Ok(store)
    }

    fn push(&mut self, rec: &Record) -> Result<Site> {
        // Editing and overlap filtering both use the exclusive reference
        // end, including for gVCF claims that extend beyond the query.
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
        self.allele_ends.push(0);
        let ends = self.allele_ends.len();
        let reference = rec.allele(0);
        let mut shared = rec.rlen.min(reference.len() as i64);
        for (i, allele) in rec.alleles().enumerate() {
            self.allele_bytes.extend_from_slice(allele);
            let end = u32::try_from(self.allele_bytes.len() - bytes)
                .map_err(|_| Error::invalid("vcf record", "alleles longer than 4 GiB"))?;
            self.allele_ends.push(end);
            if i != 0 && shared != 0 {
                shared = reference
                    .iter()
                    .zip(allele)
                    .take(shared as usize)
                    .take_while(|(r, a)| r == a)
                    .count() as i64;
            }
        }
        if rec.n_alleles() > 1 {
            // An ALT bounds the shared prefix by the checked total allele
            // byte length. A record without ALT instead uses its rlen.
            self.allele_ends[ends - 1] = shared as u32;
        }
        let genotypes = self.genotypes.len();
        let (gt, gt_width) = match rec.genotypes() {
            None => (Gt::Absent, 0),
            Some(g) => {
                if g.width == 0 || g.width.checked_mul(self.samples.len()) != Some(g.values.len()) {
                    return Err(Error::invalid("vcf record", "invalid genotype layout"));
                }
                let width = u16::try_from(g.width)
                    .map_err(|_| Error::invalid("vcf record", "ploidy above 65535"))?;
                (self.push_genotypes(g.values), width)
            }
        };
        Ok(Site {
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
            error: None,
        })
    }

    fn push_genotypes(&mut self, values: &[i32]) -> Gt {
        let byte = |v: i32| match v {
            GT_VECTOR_END => Some(GT_END),
            GT_MISSING_VALUE => Some(GT_MISSING_SAMPLE),
            0..=127 => Some(v as u8),
            _ => None,
        };
        if values.iter().all(|&v| (2..=5).contains(&v)) {
            let start = self.genotypes.len();
            self.genotypes.resize(start + values.len().div_ceil(4), 0);
            for (i, &v) in values.iter().enumerate() {
                self.genotypes[start + i / 4] |= ((v - 2) as u8) << (2 * (i % 4));
            }
            return Gt::Packed;
        }
        let start = self.genotypes.len();
        for &v in values {
            match byte(v) {
                Some(b) => self.genotypes.push(b),
                None => {
                    self.genotypes.truncate(start);
                    return Gt::Wide;
                }
            }
        }
        Gt::Bytes
    }

    /// The samples, in file order.
    pub fn samples(&self) -> &[String] {
        &self.samples
    }

    pub fn sample_index(&self, name: &str) -> Option<usize> {
        self.sample_ids.get(name).copied()
    }

    /// Records read.
    pub fn n_records(&self) -> usize {
        self.contigs.iter().map(|c| c.sites.len()).sum()
    }

    /// Malformed VCF lines skipped while reading.
    pub fn skipped_records(&self) -> usize {
        self.skipped
    }

    pub(crate) fn contig_id(&self, name: &str) -> Option<usize> {
        self.contig_ids.get(name).copied()
    }

    pub(crate) fn site(&self, contig: usize, i: usize) -> &Site {
        &self.contigs[contig].sites[i]
    }

    /// Indices of the records of `contig` that the 0-based inclusive range
    /// `[start, end]` takes, in order.
    pub(crate) fn overlapping(
        &self,
        contig: usize,
        start: i64,
        end: i64,
        overlap: Overlap,
        out: &mut Vec<usize>,
    ) {
        let c = &self.contigs[contig];
        if end < start {
            return;
        }
        if overlap == Overlap::Start {
            let first = c.sites.partition_point(|s| s.pos < start);
            let after = c.sites.partition_point(|s| s.pos <= end);
            out.extend(first..after);
            return;
        }
        let first = out.len();
        c.overlapping(start, end, out);
        if overlap == Overlap::Variant {
            self.retain_variants(c, end, first, out);
        }
    }

    // Keep variant filtering outside the Start/Record dispatch body.
    #[inline(never)]
    fn retain_variants(&self, c: &Contig, end: i64, first: usize, out: &mut Vec<usize>) {
        // Only newly appended records belong to this query; callers can
        // append several intervals to the same output buffer.
        let mut to = first;
        for from in first..out.len() {
            let i = out[from];
            if self.variant_begin(&c.sites[i]) <= end {
                out[to] = i;
                to += 1;
            }
        }
        out.truncate(to);
    }

    /// Where a record's variant part begins: after the leading bases REF
    /// shares with every ALT, at most its reference span (htslib's synced
    /// reader for `--regions-overlap 2`).
    fn variant_begin(&self, site: &Site) -> i64 {
        let shared = if site.n_allele == 1 {
            site.rlen
        } else {
            self.allele_ends[site.ends - 1] as i64
        };
        site.pos + shared
    }

    /// The filter's error on this record.
    pub(crate) fn error(&self, site: &Site) -> Option<&Error> {
        site.error.map(|i| &self.errors[i as usize])
    }

    /// Allele `i` (0 is REF).
    pub(crate) fn allele(&self, site: &Site, i: usize) -> &[u8] {
        let begin = if i == 0 {
            0
        } else {
            self.allele_ends[site.ends + i - 1] as usize
        };
        let end = self.allele_ends[site.ends + i] as usize;
        &self.allele_bytes[site.bytes + begin..site.bytes + end]
    }

    /// The genotype values of `sample` (`gt_width` of them) into `out`;
    /// false without genotypes in `int8`.
    pub(crate) fn genotype(&self, site: &Site, sample: usize, out: &mut Vec<u8>) -> bool {
        out.clear();
        let width = site.gt_width as usize;
        let first = sample * width;
        match site.gt {
            Gt::Absent | Gt::Wide => return false,
            Gt::Bytes => out.extend_from_slice(
                &self.genotypes[site.genotypes + first..site.genotypes + first + width],
            ),
            Gt::Packed => out.extend(
                (first..first + width)
                    .map(|i| 2 + ((self.genotypes[site.genotypes + i / 4] >> (2 * (i % 4))) & 3)),
            ),
        }
        true
    }
}

#[cfg(test)]
#[path = "store_tests.rs"]
mod tests;
