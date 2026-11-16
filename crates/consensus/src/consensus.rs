//! [`Consensus`]: personal sequences of a region, as `bcftools consensus`
//! builds them, from the records a BCF index gives for the region and a
//! reference read through [`SequenceRead`].

use std::collections::{HashMap, HashSet};
use std::sync::Arc;

use omicsweave_core::orient::complement;
use omicsweave_core::{
    Error, FileSource, Locus, ReadAt, Region, Result, SequenceRead, Strand, check_output,
};
use omicsweave_vcf::{BcfReader, Filter, Header, Record};

use crate::chain::Chain;
use crate::edit::{Edit, Marks};
use crate::mask::Mask;
use crate::records::{Overlap, RecordFilter, Records};
use crate::select::{Haplotype, Mode, Scratch, choose};

/// IUPAC selection is a union: inspect each sample only once per record,
/// retaining first-occurrence order so malformed genotypes fail in order.
fn select_samples(
    names: &[&str],
    mut sample: impl FnMut(&str) -> Result<usize>,
) -> Result<smallvec::SmallVec<[usize; 4]>> {
    let mut selected = smallvec::SmallVec::<[usize; 4]>::new();
    let mut seen = None::<HashSet<usize>>;
    let mut increasing = true;
    for name in names {
        let id = sample(name)?;
        let fresh = match &mut seen {
            Some(seen) => seen.insert(id),
            None if selected.len() <= 4 && selected.contains(&id) => false,
            None if selected.len() < 4 => true,
            // Lists in header order are already unique without a hash table.
            None if increasing && id >= selected[selected.len() - 1] => {
                id != selected[selected.len() - 1]
            }
            None => {
                let mut ids: HashSet<_> = selected.iter().copied().collect();
                let fresh = ids.insert(id);
                seen = Some(ids);
                fresh
            }
        };
        if fresh {
            increasing &= selected.last().is_none_or(|&last| id > last);
            selected.push(id);
        }
    }
    Ok(selected)
}

/// The `bcftools consensus` options that shape every sequence; which
/// alleles a sequence takes is chosen per query ([`Alleles`]).
#[derive(Debug, Clone, Default)]
pub struct ConsensusOptions {
    /// `-I`: IUPAC codes, of REF and ALT without samples, or of the
    /// sample's genotype with a haplotype.
    pub iupac_codes: bool,
    /// `-M`: write this for a missing genotype instead of leaving the
    /// reference.
    pub missing: Option<u8>,
    /// `-a`: write this where no record covers the reference.
    pub absent: Option<u8>,
    /// `--mark-del`: write this for each deleted base instead of deleting.
    pub mark_del: Option<u8>,
    /// `--mark-ins`.
    pub mark_ins: Option<Mark>,
    /// `--mark-snv`.
    pub mark_snv: Option<Mark>,
    /// `-m` with `--mask-with`.
    pub mask: Option<Arc<Mask>>,
    /// `--regions-overlap`.
    pub regions_overlap: Overlap,
    /// `-c`: also build the chain from reference to consensus.
    pub chain: bool,
}

/// How inserted or substituted bases are marked.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Mark {
    /// `uc`
    Upper,
    /// `lc`
    Lower,
    /// This character over each base.
    Char(u8),
}

/// Which alleles a sequence takes (`-s`, `-H`).
#[derive(Debug, Clone, Copy)]
pub enum Alleles<'a> {
    /// The first ALT of every record, samples ignored (`-s -`); IUPAC codes
    /// of REF and ALT with [`ConsensusOptions::iupac_codes`].
    AllAlt,
    /// IUPAC codes of these samples' genotypes (`-s a,b` without `-H`).
    Samples(&'a [&'a str]),
    /// One haplotype of one sample: `-s x -H` with `1`, `2`, ..., `R`, `A`,
    /// `LR`, `LA`, `SR`, `SA`, `I` or `<N>pIu`.
    Haplotype { sample: &'a str, haplotype: &'a str },
}

/// What every context of one variant file shares: sample lookup and the
/// compiled `-i`/`-e` expression (`true`: include).
struct Shared {
    sample_ids: HashMap<String, usize>,
    filter: Option<(Filter, bool)>,
}

/// A loaded region: its reference, masked, and the records it takes.
#[derive(Default)]
struct Loaded {
    chrom: String,
    start: i64,
    reference: Vec<u8>,
    records: Records,
    /// Per record: left out by the mask.
    masked: Vec<bool>,
}

/// Builds consensus sequences from a BCF read by region through its CSI.
/// [`Consensus::load`] reads a region's reference and decodes the records
/// the index gives for it, each once; each [`Consensus::query`] then builds
/// one sequence of it (one sample's haplotype, say), held until the next
/// query and copied out by [`Consensus::write`]. Builders for other threads
/// come from [`Consensus::clone_context`] or [`Consensus::try_clone`].
pub struct Consensus<S = FileSource> {
    variants: BcfReader<S>,
    header: Arc<Header>,
    shared: Arc<Shared>,
    options: ConsensusOptions,
    record: Record,
    loaded: Option<Loaded>,
    seq: Vec<u8>,
    ready: bool,
    chain: String,
    chain_builder: Option<Chain>,
    scratch: Scratch,
    spare: Vec<u8>,
}

impl<S: ReadAt> Consensus<S> {
    /// Builds from `variants`, applying the records `filter` passes. An
    /// expression that does not compile is an error here.
    pub fn new(
        variants: BcfReader<S>,
        filter: Option<RecordFilter<'_>>,
        options: ConsensusOptions,
    ) -> Result<Consensus<S>> {
        let header = Arc::clone(variants.shared_header());
        let filter = match filter {
            None => None,
            Some(RecordFilter::Include(e)) => Some((Filter::new(&header, e)?, true)),
            Some(RecordFilter::Exclude(e)) => Some((Filter::new(&header, e)?, false)),
        };
        let sample_ids = header
            .samples
            .iter()
            .enumerate()
            .map(|(i, s)| (s.clone(), i))
            .collect();
        let shared = Arc::new(Shared { sample_ids, filter });
        Ok(Self::with_parts(variants, header, shared, options))
    }

    fn with_parts(
        variants: BcfReader<S>,
        header: Arc<Header>,
        shared: Arc<Shared>,
        options: ConsensusOptions,
    ) -> Consensus<S> {
        let chain_builder = options.chain.then(|| Chain::new(0, 0));
        Consensus {
            variants,
            header,
            shared,
            options,
            record: Record::default(),
            loaded: None,
            seq: Vec::new(),
            ready: false,
            chain: String::new(),
            chain_builder,
            scratch: Scratch::default(),
            spare: Vec::new(),
        }
    }

    pub fn variants(&self) -> &BcfReader<S> {
        &self.variants
    }

    /// The variant reader, e.g. to attach a block cache.
    pub fn variants_mut(&mut self) -> &mut BcfReader<S> {
        &mut self.variants
    }

    /// The samples, in header order.
    pub fn samples(&self) -> &[String] {
        &self.header.samples
    }

    pub fn sample_index(&self, name: &str) -> Option<usize> {
        self.shared.sample_ids.get(name).copied()
    }

    pub fn options(&self) -> &ConsensusOptions {
        &self.options
    }

    /// Shares the header, filter and options with an independent execution
    /// context. Loaded reference, records, output and scratch start empty.
    pub fn clone_context(&self) -> Self
    where
        S: Clone,
    {
        self.fork(self.variants.clone_context())
    }

    fn fork(&self, variants: BcfReader<S>) -> Self {
        Self::with_parts(
            variants,
            Arc::clone(&self.header),
            Arc::clone(&self.shared),
            self.options.clone(),
        )
    }

    /// Reads the reference of `region` (which must lie within its
    /// chromosome) and decodes the records it takes, for the queries that
    /// follow. A chromosome the BCF does not declare has no records, as
    /// bcftools only warns of it.
    pub fn load(&mut self, reference: &mut impl SequenceRead, region: Region<'_>) -> Result<()> {
        self.ready = false;
        let mut loaded = self.loaded.take().unwrap_or_default();
        region.check()?;
        region.check_within(reference.chrom_len(region.chrom)?)?;
        loaded.reference.clear();
        loaded.reference.resize(region.len() as usize, 0);
        reference.sequence(&Locus::from(region), &mut loaded.reference)?;
        loaded.chrom.clear();
        loaded.chrom.push_str(region.chrom);
        loaded.start = region.start as i64;
        loaded.masked.clear();
        if let Some(mask) = &self.options.mask {
            mask.apply(region.chrom, loaded.start, &mut loaded.reference);
        }
        self.read_records(region, &mut loaded.records)?;
        let skips = self.options.mask.as_ref().filter(|m| m.skips_records());
        if let Some(mask) = skips {
            mask.overlapping_records(
                region.chrom,
                loaded.records.sites.iter().map(|s| (s.pos, s.last())),
                &mut loaded.masked,
            );
        } else {
            loaded.masked.resize(loaded.records.sites.len(), false);
        }
        self.loaded = Some(loaded);
        Ok(())
    }

    /// Decodes the records of `region` that `--regions-overlap` takes and
    /// the filter passes, in file order: the records htslib's synced reader
    /// gives bcftools for the region. The filter starts fresh per region,
    /// as in a `bcftools consensus` run over that one region.
    fn read_records(&mut self, region: Region<'_>, records: &mut Records) -> Result<()> {
        records.clear();
        if region.is_empty() || self.header.contig_id(region.chrom).is_none() {
            return Ok(());
        }
        let Consensus {
            variants,
            header,
            shared,
            options,
            record,
            ..
        } = self;
        let header: &Header = header;
        let n_samples = header.samples.len();
        let (start, end) = (region.start as i64, region.end as i64 - 1);
        let overlap = options.regions_overlap;
        let mut filter = shared.filter.clone();
        variants.visit(region, record, |rec| {
            if !overlap.takes(rec, start, end) {
                return Ok(());
            }
            let mut error = None;
            if let Some((filter, include)) = &mut filter {
                match filter.test(header, rec) {
                    Ok(pass) if pass != *include => return Ok(()),
                    Ok(_) => {}
                    Err(e) => error = Some(e),
                }
            }
            records.push(rec, n_samples, error)
        })
    }

    fn mode(&self, alleles: Alleles<'_>) -> Result<Mode> {
        let sample = |name: &str| {
            self.sample_index(name)
                .ok_or_else(|| Error::invalid("consensus", format!("no sample {name:?}")))
        };
        Ok(match alleles {
            Alleles::AllAlt if self.options.iupac_codes => Mode::IupacRefAlt,
            Alleles::AllAlt => Mode::AllAlt,
            Alleles::Samples([]) => return Err(Error::invalid("consensus", "no samples")),
            Alleles::Samples(names) => Mode::IupacSamples(select_samples(names, sample)?),
            Alleles::Haplotype {
                sample: name,
                haplotype,
            } => Mode::Sample {
                sample: sample(name)?,
                haplotype: Haplotype::parse(haplotype).ok_or_else(|| {
                    Error::invalid("consensus", format!("bad haplotype {haplotype:?}"))
                })?,
                iupac: self.options.iupac_codes,
            },
        })
    }

    /// Builds the consensus of the loaded region with `alleles` and holds
    /// it; returns its length, which differs from the region's where indels
    /// apply.
    pub fn query(&mut self, alleles: Alleles<'_>) -> Result<usize> {
        self.ready = false;
        let mode = self.mode(alleles)?;
        let Consensus {
            options: o,
            loaded,
            seq,
            chain: chain_text,
            chain_builder,
            scratch,
            spare,
            ..
        } = self;
        let loaded = loaded
            .as_ref()
            .ok_or_else(|| Error::invalid("consensus", "query before load"))?;
        let records = &loaded.records;
        let marks = Marks {
            missing: o.missing,
            absent: o.absent,
            mark_del: o.mark_del,
            mark_ins: o.mark_ins,
            mark_snv: o.mark_snv,
        };
        let mut chain = chain_builder.as_mut();
        if let Some(chain) = &mut chain {
            chain.reset(loaded.start, loaded.reference.len() as i64);
        }
        let mut edit = Edit::new(&loaded.reference, loaded.start, std::mem::take(seq));
        for (site, &masked) in records.sites.iter().zip(&loaded.masked) {
            if let Some(e) = records.error(site) {
                *seq = edit.abandon();
                return Err(e.clone());
            }
            if let Some(a) = o.absent {
                edit.absent(site.pos, a);
            }
            if (site.n_allele == 1 && o.missing.is_none() && o.absent.is_none()) || masked {
                continue;
            }
            let applied =
                choose(records, site, &mode, o.missing.is_some(), scratch).and_then(|choice| {
                    edit.apply(
                        site,
                        choice,
                        records.allele(site, 0),
                        &mut scratch.allele,
                        spare,
                        &marks,
                        chain.as_deref_mut(),
                    )
                });
            if let Err(e) = applied {
                *seq = edit.abandon();
                return Err(e);
            }
        }
        if let Some(a) = o.absent {
            edit.absent(i64::MAX, a);
        }
        *seq = edit.finish();
        if let Some(chain) = chain {
            chain.render(&loaded.chrom, chain_text);
        }
        self.ready = true;
        Ok(self.seq.len())
    }

    /// Copies the last query's sequence into `out` (its length), 5' to 3'
    /// of `strand`: the reverse complement (IUPAC codes complemented, case
    /// kept) for `-`, the sequence itself otherwise. The chain stays on the
    /// forward strand.
    pub fn write(&self, strand: Strand, out: &mut [u8]) -> Result<()> {
        if !self.ready {
            return Err(Error::NotQueried);
        }
        check_output("consensus sequence", out.len(), self.seq.len())?;
        if strand.is_reverse() {
            for (o, &b) in out.iter_mut().zip(self.seq.iter().rev()) {
                *o = complement(b);
            }
        } else {
            out.copy_from_slice(&self.seq);
        }
        Ok(())
    }

    /// The last query's chain (UCSC chain format), with
    /// [`ConsensusOptions::chain`].
    pub fn chain(&self) -> Option<&str> {
        (self.ready && self.options.chain).then_some(self.chain.as_str())
    }
}

impl Consensus<FileSource> {
    /// Another execution context on a duplicated file handle, sharing the
    /// header, index, filter and options.
    pub fn try_clone(&self) -> Result<Self> {
        Ok(self.fork(self.variants.try_clone()?))
    }
}

#[cfg(test)]
mod tests {
    use super::select_samples;

    #[test]
    fn sample_union_resolves_all_names_and_keeps_first_occurrence_order() {
        for (unique, reverse) in [1, 4, 5, 1024]
            .into_iter()
            .flat_map(|n| [(n, false), (n, true)])
        {
            let names: Vec<_> = (0..unique)
                .map(|i| if reverse { unique - i - 1 } else { i })
                .map(|i| i.to_string())
                .collect();
            let requested: Vec<_> = names
                .iter()
                .cycle()
                .take(unique * 17)
                .map(String::as_str)
                .collect();
            let mut lookups = 0;
            let selected = select_samples(&requested, |s| {
                lookups += 1;
                Ok(s.parse().unwrap())
            })
            .unwrap();
            assert_eq!(lookups, requested.len());
            assert_eq!(
                &selected[..],
                &names
                    .iter()
                    .map(|s| s.parse::<usize>().unwrap())
                    .collect::<Vec<_>>()
            );
            assert_eq!(selected.spilled(), unique > 4);
        }
        let mut lookups = 0;
        assert!(
            select_samples(&["0", "0", "bad", "1"], |s| {
                lookups += 1;
                s.parse()
                    .map_err(|_| omicsweave_core::Error::UnknownChrom(s.to_owned()))
            })
            .is_err()
        );
        assert_eq!(lookups, 3);
    }
}
