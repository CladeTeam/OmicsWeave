//! [`Consensus`]: personal sequences of a region, as `bcftools consensus`
//! builds them, from a resident [`VcfStore`] and a reference read through
//! [`SequenceRead`].

use std::sync::Arc;

use omicsweave_core::{Error, Region, Result, SequenceRead, check_output};

use crate::chain::Chain;
use crate::edit::{Edit, Marks};
use crate::mask::Mask;
use crate::select::{Haplotype, Mode, Scratch, choose};
use crate::store::{Overlap, VcfStore};

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

/// A loaded region: its reference, masked, and the records it takes.
struct Loaded {
    chrom: String,
    start: i64,
    reference: Vec<u8>,
    contig: Option<usize>,
    sites: Vec<usize>,
    /// Per record: left out by the mask.
    masked: Vec<bool>,
}

/// Builds consensus sequences. [`Consensus::load`] reads a region's
/// reference and finds its records once; each [`Consensus::query`] then
/// builds one sequence of it (one sample's haplotype, say), held until the
/// next query and copied out by [`Consensus::write`]. Many builders, one
/// per thread, can share a store.
pub struct Consensus {
    store: Arc<VcfStore>,
    options: ConsensusOptions,
    loaded: Option<Loaded>,
    seq: Vec<u8>,
    ready: bool,
    chain: String,
    chain_builder: Option<Chain>,
    scratch: Scratch,
    spare: Vec<u8>,
}

impl Consensus {
    pub fn new(store: Arc<VcfStore>, options: ConsensusOptions) -> Consensus {
        let chain_builder = options.chain.then(|| Chain::new(0, 0));
        Consensus {
            store,
            options,
            loaded: None,
            seq: Vec::new(),
            ready: false,
            chain: String::new(),
            chain_builder,
            scratch: Scratch::default(),
            spare: Vec::new(),
        }
    }

    pub fn store(&self) -> &VcfStore {
        &self.store
    }

    pub fn options(&self) -> &ConsensusOptions {
        &self.options
    }

    /// Reads the reference of `region` (which must lie within its
    /// chromosome) and finds the records it takes, for the queries that
    /// follow.
    pub fn load(&mut self, reference: &mut impl SequenceRead, region: Region<'_>) -> Result<()> {
        self.ready = false;
        region.check()?;
        let mut loaded = self.loaded.take().unwrap_or(Loaded {
            chrom: String::new(),
            start: 0,
            reference: Vec::new(),
            contig: None,
            sites: Vec::new(),
            masked: Vec::new(),
        });
        loaded.reference.clear();
        loaded.reference.resize(region.len() as usize, 0);
        reference.sequence(region, &mut loaded.reference)?;
        loaded.chrom.clear();
        loaded.chrom.push_str(region.chrom);
        loaded.start = region.start as i64;
        loaded.contig = self.store.contig_id(region.chrom);
        loaded.sites.clear();
        loaded.masked.clear();
        if let Some(mask) = &self.options.mask {
            mask.apply(region.chrom, loaded.start, &mut loaded.reference);
        }
        if let (Some(contig), false) = (loaded.contig, region.is_empty()) {
            let end = region.end as i64 - 1;
            let overlap = self.options.regions_overlap;
            self.store
                .overlapping(contig, loaded.start, end, overlap, &mut loaded.sites);
            let skips = self.options.mask.as_ref().filter(|m| m.skips_records());
            if let Some(mask) = skips {
                mask.overlapping_records(
                    region.chrom,
                    loaded.sites.iter().map(|&i| {
                        let site = self.store.site(contig, i);
                        (site.pos, site.last())
                    }),
                    &mut loaded.masked,
                );
            } else {
                loaded.masked.resize(loaded.sites.len(), false);
            }
        }
        self.loaded = Some(loaded);
        Ok(())
    }

    fn mode(&self, alleles: Alleles<'_>) -> Result<Mode> {
        let sample = |name: &str| {
            self.store
                .sample_index(name)
                .ok_or_else(|| Error::invalid("consensus", format!("no sample {name:?}")))
        };
        Ok(match alleles {
            Alleles::AllAlt if self.options.iupac_codes => Mode::IupacRefAlt,
            Alleles::AllAlt => Mode::AllAlt,
            Alleles::Samples([]) => return Err(Error::invalid("consensus", "no samples")),
            Alleles::Samples(names) => {
                Mode::IupacSamples(names.iter().map(|n| sample(n)).collect::<Result<_>>()?)
            }
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
            store,
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
        if let Some(contig) = loaded.contig {
            for (&i, &masked) in loaded.sites.iter().zip(&loaded.masked) {
                let site = store.site(contig, i);
                if let Some(e) = store.error(site) {
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
                    choose(store, site, &mode, o.missing.is_some(), scratch).and_then(|choice| {
                        edit.apply(
                            site,
                            choice,
                            store.allele(site, 0),
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

    /// Copies the last query's sequence into `out` (its length).
    pub fn write(&self, out: &mut [u8]) -> Result<()> {
        if !self.ready {
            return Err(Error::NotQueried);
        }
        check_output("consensus sequence", out.len(), self.seq.len())?;
        out.copy_from_slice(&self.seq);
        Ok(())
    }

    /// The last query's chain (UCSC chain format), with
    /// [`ConsensusOptions::chain`].
    pub fn chain(&self) -> Option<&str> {
        (self.ready && self.options.chain).then_some(self.chain.as_str())
    }
}
