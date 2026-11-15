//! Allele counts at known SNVs (cellsnp-lite mode 1, vartrix): variant ×
//! row matrices of the units supporting each site's REF, ALT or another
//! base, a per-variant table, and gene-level allele-specific expression
//! from one sample's phased genotypes.
//!
//! Definitions:
//! - The **sites** are the biallelic SNVs of a BCF (with its CSI), every
//!   other record skipped; all of them are columns, in file order.
//! - A read passing the filter **observes** a site when an aligned base
//!   (`M` / `=` / `X`) lies on it with quality `min_quality` or more (all
//!   bases of a record without qualities): `ref` when it is the REF base,
//!   `alt` the ALT base, `other` any other of `A`, `C`, `G`, `T` (`N`
//!   observes nothing). Deletions and introns over a site observe nothing.
//! - **Units**: with `umi`, distinct UMIs (exact values of `umi_tag`; reads
//!   without one are not counted) per row, site and allele; otherwise
//!   distinct read names, so both mates of a pair count once. A unit whose
//!   reads disagree counts in each allele it shows.
//! - **Rows**: a tag (`CB`), the read group or the input file, optionally
//!   a whitelist.
//! - **Gene-level ASE** (`genes`): at sites where the phasing sample is
//!   heterozygous and phased (`0|1`, `1|0`), a unit observing haplotype
//!   1's allele is evidence for haplotype 1, the other allele for
//!   haplotype 2; a unit is counted in the gene its read's aligned blocks
//!   overlap the exons of (exactly one gene, strand by the library) as
//!   `hap1`, `hap2`, or `both` when it shows both.
//!
//! Products: `<name>/` (10x directory with layers `ref`, `alt`, `other`:
//! `ref.cpk`, `alt.cpk`, `other.cpk`), `<name>/variants.tsv.gz` (per site:
//! REF / ALT / other units over all rows and `alt / (ref + alt)`), and
//! `<name>/genes/` (layers `hap1`, `hap2`, `both`).

use std::path::PathBuf;
use std::sync::Arc;

use omicsweave_bam::flag;
use omicsweave_convert::stage::Sorter;
use omicsweave_core::{Error, Result};

use super::counts::{
    Entry, Fold, GlobalRows, LocalRows, Reduce, RowKey, Rows, VELOCITY_LAYER, collapse,
};
use super::coverage::Tally;
use super::velocity::{EXONIC, INTRONIC};
use super::{Accumulator, BatchInfo, Delivery, MemoryNeed, Needs, PlanInfo, RunCtx};
use crate::assign::{Query, SpaceId};
use crate::model::{BarcodeTable, Cells, CountMatrix, FeatureList, FeatureTable, VariantTable};
use crate::prelude::{Filter, ReadView, Strandedness, fnv, umi_key};
use crate::reference::SitePanel;
use crate::sink::json::Json;

/// Gene-level allele-specific expression of an [`Alleles`] product.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct AlleleGenes {
    /// The sample whose phased genotypes give the haplotypes (`None`: the
    /// first of the BCF).
    pub sample: Option<String>,
    /// How reads relate to the genes' strand.
    pub strand: Strandedness,
}

/// An allele-count product; see the module documentation.
#[derive(Debug, Clone)]
pub struct Alleles {
    /// The output directory (`alleles`).
    pub name: String,
    /// The SNVs (BCF with its CSI).
    pub sites: PathBuf,
    pub rows: RowKey,
    pub whitelist: Option<PathBuf>,
    /// Count UMIs; otherwise read names.
    pub umi: bool,
    pub umi_tag: [u8; 2],
    pub filter: Filter,
    pub min_quality: u8,
    /// Write `variants.tsv.gz`.
    pub table: bool,
    pub genes: Option<AlleleGenes>,
}

impl Alleles {
    /// Single-cell (cellsnp-lite): sites × `CB`, UMIs, primary mapped
    /// records that are not QC failures, MAPQ 20 or more.
    pub fn single_cell(sites: PathBuf) -> Alleles {
        Alleles {
            name: "alleles".into(),
            sites,
            rows: RowKey::Tag(*b"CB"),
            whitelist: None,
            umi: true,
            umi_tag: *b"UB",
            filter: Filter {
                require: 0,
                exclude: flag::UNMAPPED | flag::SECONDARY | flag::QC_FAIL,
                min_mapq: 20,
                unique: false,
            },
            min_quality: 20,
            table: true,
            genes: None,
        }
    }

    /// Bulk: sites × input file, read pairs (by name), duplicates excluded
    /// too.
    pub fn bulk(sites: PathBuf) -> Alleles {
        Alleles {
            rows: RowKey::File,
            umi: false,
            filter: Filter {
                exclude: flag::UNMAPPED | flag::SECONDARY | flag::QC_FAIL | flag::DUPLICATE,
                ..Alleles::single_cell(PathBuf::new()).filter
            },
            ..Alleles::single_cell(sites)
        }
    }
}

/// Why reads were or were not counted.
#[derive(Debug, Clone, Copy, Default)]
pub struct AlleleStatus {
    /// Reads observing at least one site, and the observations.
    pub reads: u64,
    pub observations: u64,
    pub no_row: u64,
    pub no_umi: u64,
    /// Observations of a base below `min_quality`.
    pub low_quality: u64,
    /// ASE: units with haplotype evidence over no gene, several genes.
    pub ase_no_gene: u64,
    pub ase_multi_gene: u64,
}

impl AlleleStatus {
    fn add(&mut self, o: &AlleleStatus) {
        self.reads += o.reads;
        self.observations += o.observations;
        self.no_row += o.no_row;
        self.no_umi += o.no_umi;
        self.low_quality += o.low_quality;
        self.ase_no_gene += o.ase_no_gene;
        self.ase_multi_gene += o.ase_multi_gene;
    }
}

/// An [`Alleles`] product resolved against a task.
pub struct AlleleStage {
    pub(crate) spec: Alleles,
    pub(crate) panel: Arc<SitePanel>,
    pub(crate) rows: Rows,
    /// Genes by their exons, for ASE.
    pub(crate) genes: Option<(SpaceId, Arc<FeatureTable>)>,
    pub(crate) table_bytes: u64,
}

pub struct AlleleBatch {
    variants: Vec<Entry>,
    genes: Vec<Entry>,
    local: LocalRows,
    tally: Tally,
    status: AlleleStatus,
}

pub struct AlleleTotal<'p> {
    variants: Sorter<'p, Entry>,
    genes: Sorter<'p, Entry>,
    rows: GlobalRows,
    tally: Tally,
    status: AlleleStatus,
}

impl AlleleStage {
    fn reduce(&self) -> Reduce {
        Reduce::Distinct
    }

    fn gene_fold(&self) -> Fold {
        Fold {
            spliced: Some(0),
            ..Fold::plain(self.reduce())
        }
    }
}

impl Accumulator for AlleleStage {
    type Batch = AlleleBatch;
    type Total<'p> = AlleleTotal<'p>;
    type Product = (CountMatrix, Option<CountMatrix>, Tally, AlleleStatus, u64);

    fn name(&self) -> &str {
        &self.spec.name
    }

    fn needs(&self) -> Needs {
        let mut n = Needs::CIGAR | Needs::TAGS | Needs::SEQUENCE | Needs::QUALITY;
        if self.genes.is_some() {
            n |= Needs::ASSIGN;
        }
        n
    }

    fn estimate(&self, info: &PlanInfo) -> MemoryNeed {
        let sums = if self.spec.table {
            self.panel.sites.len() as u64 * 24
        } else {
            0
        };
        MemoryNeed {
            fixed: self.panel.heap_bytes()
                + self.table_bytes
                + sums
                + self.genes.as_ref().map_or(0, |g| g.1.heap_bytes()),
            per_batch: info.batch_bytes / 8,
            min_buffer: 2 << 20,
            weight: 4,
        }
    }

    fn batch(&self, _info: BatchInfo) -> AlleleBatch {
        AlleleBatch {
            variants: Vec::new(),
            genes: Vec::new(),
            local: LocalRows::default(),
            tally: Tally::default(),
            status: AlleleStatus::default(),
        }
    }

    fn observe(&self, b: &mut AlleleBatch, r: &ReadView<'_>) -> Result<()> {
        b.tally.seen += 1;
        if let Some(why) = self.spec.filter.reject(r) {
            b.tally.reject(why);
            return Ok(());
        }
        let (first, sites) = self.panel.sites_of(r.input as usize, r.tid());
        let (lo, hi) = (r.pos(), r.end());
        let at = sites.partition_point(|s| s.pos < lo);
        if sites.get(at).is_none_or(|s| s.pos >= hi) {
            return Ok(());
        }
        let Some(row) = self.rows.of(r, &mut b.local) else {
            b.status.no_row += 1;
            return Ok(());
        };
        let sub = if self.spec.umi {
            match r.tag_bytes(self.spec.umi_tag) {
                Some(u) => umi_key(u),
                None => {
                    b.status.no_umi += 1;
                    return Ok(());
                }
            }
        } else {
            fnv(r.name())
        };
        let qual = r.rec.qual();
        let no_qual = qual.first() == Some(&0xff);
        let mut haps = 0u32;
        let mut observed = 0u64;
        for (k, s) in sites[at..].iter().enumerate() {
            if s.pos >= hi {
                break;
            }
            let Some(q) = r.query_at(s.pos) else { continue };
            if !no_qual && qual[q as usize] < self.spec.min_quality {
                b.status.low_quality += 1;
                continue;
            }
            let base = match r.base(q as usize) {
                b'=' => s.reference,
                x => x,
            };
            let allele = if base == s.reference {
                0
            } else if base == s.alt {
                1
            } else if matches!(base, b'A' | b'C' | b'G' | b'T') {
                2
            } else {
                continue;
            };
            observed += 1;
            b.variants.push(Entry {
                row,
                col: (first + at + k) as u32,
                layer: allele,
                sub,
                value: 1,
            });
            if let (Some(h), true) = (s.hap1, allele < 2) {
                haps |= if allele as u8 == h { EXONIC } else { INTRONIC };
            }
        }
        if observed > 0 {
            b.status.reads += 1;
            b.status.observations += observed;
        }
        if let (Some((space, _)), true) = (&self.genes, haps != 0) {
            let strand = match self.spec.genes.as_ref().map(|g| g.strand) {
                Some(Strandedness::Unstranded) | None => None,
                Some(lib) => Some(r.reverse_for(lib)),
            };
            let hits = r.assign(*space, Query::Blocks, strand);
            match hits.len() {
                0 => b.status.ase_no_gene += 1,
                1 => b.genes.push(Entry {
                    row,
                    col: hits[0],
                    layer: VELOCITY_LAYER,
                    sub,
                    value: haps,
                }),
                _ => b.status.ase_multi_gene += 1,
            }
        }
        Ok(())
    }

    fn seal(&self, b: &mut AlleleBatch) -> Result<()> {
        collapse(&mut b.variants, Fold::direct(self.reduce()));
        collapse(&mut b.genes, self.gene_fold());
        Ok(())
    }

    fn total<'p>(&'p self, rc: &RunCtx<'p>) -> Result<AlleleTotal<'p>> {
        let genes = if self.genes.is_some() {
            rc.buffer / 4
        } else {
            0
        };
        Ok(AlleleTotal {
            variants: Sorter::new(rc.pool, rc.buffer - genes, &rc.tmp),
            genes: Sorter::new(rc.pool, genes.max(1 << 20), &rc.tmp),
            rows: GlobalRows::default(),
            tally: Tally::default(),
            status: AlleleStatus::default(),
        })
    }

    fn merge<'p>(
        &'p self,
        t: &mut AlleleTotal<'p>,
        b: AlleleBatch,
        _frontier: u64,
        rc: &RunCtx<'p>,
    ) -> Result<()> {
        t.tally.add(&b.tally);
        t.status.add(&b.status);
        let AlleleBatch {
            mut variants,
            mut genes,
            local,
            ..
        } = b;
        if self.rows.local() {
            let map = t.rows.remap(local);
            for e in variants.iter_mut().chain(genes.iter_mut()) {
                e.row = map[e.row as usize];
            }
            if t.rows.bytes() > rc.buffer.max(64 << 20) {
                return Err(Error::invalid(
                    "budget",
                    format!(
                        "{}: {} row keys outgrow the memory budget; give a whitelist or raise \
                         the memory",
                        self.spec.name,
                        t.rows.names.len()
                    ),
                ));
            }
        }
        t.variants.extend(variants)?;
        t.genes.extend(genes)
    }

    fn finish<'p>(
        &'p self,
        t: AlleleTotal<'p>,
        _rc: &RunCtx<'p>,
    ) -> Result<(CountMatrix, Option<CountMatrix>, Tally, AlleleStatus, u64)> {
        let spilled = t.variants.spilled() + t.genes.spilled();
        let names = match self.rows.fixed_names() {
            Some(n) => n.to_vec(),
            None => t.rows.names,
        };
        let variants = CountMatrix {
            rows: BarcodeTable {
                names: names.clone(),
            },
            features: Arc::new(self.panel.features()) as Arc<dyn FeatureList>,
            layers: vec!["ref".into(), "alt".into(), "other".into()],
            named_layers: true,
            cells: Cells {
                merge: Some(t.variants.finish()?),
                fold: Fold::direct(self.reduce()),
                head: None,
            },
        };
        let genes = match &self.genes {
            Some((_, features)) => Some(CountMatrix {
                rows: BarcodeTable { names },
                features: Arc::clone(features) as Arc<dyn FeatureList>,
                layers: vec!["hap1".into(), "hap2".into(), "both".into()],
                named_layers: true,
                cells: Cells {
                    merge: Some(t.genes.finish()?),
                    fold: self.gene_fold(),
                    head: None,
                },
            }),
            None => None,
        };
        Ok((variants, genes, t.tally, t.status, spilled))
    }

    fn deliver(
        &self,
        (variants, genes, tally, status, spilled): (
            CountMatrix,
            Option<CountMatrix>,
            Tally,
            AlleleStatus,
            u64,
        ),
        d: &mut Delivery<'_>,
    ) -> Result<Json> {
        let name = &self.spec.name;
        let n_rows = variants.rows.len();
        let w = crate::sink::matrix::write_with(variants, name, false, self.spec.table, d)?;
        let mut files: Vec<Json> = w.files.iter().cloned().map(Json::from).collect();
        let nonzero = w.nonzero;
        let mut covered = 0u64;
        if self.spec.table {
            let n = self.panel.sites.len();
            let counts: Vec<[u64; 3]> = (0..n)
                .map(|i| {
                    [
                        w.column_sums[0][i],
                        w.column_sums[1][i],
                        w.column_sums[2][i],
                    ]
                })
                .collect();
            covered = counts.iter().filter(|c| c.iter().any(|&v| v > 0)).count() as u64;
            let table = VariantTable {
                contigs: Arc::clone(&self.panel.contigs),
                sites: Arc::clone(&self.panel.sites),
                counts,
            };
            let rel = format!("{name}/variants.tsv.gz");
            let path = d.path(&rel)?;
            crate::sink::table::write_variants(&table, &path, d)?;
            files.push(Json::from(rel));
        }
        let mut gene_fields = None;
        if let Some(g) = genes {
            let n_genes = g.features.len();
            let w = crate::sink::matrix::write(g, &format!("{name}/genes"), false, d)?;
            files.extend(w.files.into_iter().map(Json::from));
            gene_fields = Some(Json::obj([
                ("genes", Json::from(n_genes)),
                ("nonzero", Json::from(w.nonzero)),
                (
                    "sample",
                    Json::from(self.panel.sample.clone().unwrap_or_default()),
                ),
                ("phased_sites", Json::from(self.panel.phased)),
                ("unassigned_no_gene", Json::from(status.ase_no_gene)),
                ("unassigned_multi_gene", Json::from(status.ase_multi_gene)),
            ]));
        }
        let mut fields = vec![
            ("product", Json::from("alleles")),
            (
                "unit",
                Json::from(if self.spec.umi {
                    "molecule"
                } else {
                    "fragment"
                }),
            ),
            ("rows", Json::from(n_rows)),
            ("sites", Json::from(self.panel.sites.len())),
            (
                "sites_skipped_not_biallelic_snv",
                Json::from(self.panel.skipped),
            ),
            ("sites_covered", Json::from(covered)),
            ("nonzero", Json::from(nonzero)),
            ("reads_observing", Json::from(status.reads)),
            ("observations", Json::from(status.observations)),
            ("observations_low_quality", Json::from(status.low_quality)),
            ("unassigned_no_row_key", Json::from(status.no_row)),
            ("unassigned_no_umi", Json::from(status.no_umi)),
            ("spilled_bytes", Json::from(spilled)),
        ];
        if let Some(g) = gene_fields {
            fields.push(("genes", g));
        }
        fields.extend(tally.json());
        fields.push(("files", Json::Arr(files)));
        Ok(Json::obj(fields))
    }
}
