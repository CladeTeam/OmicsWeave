//! The reference genome and the known variants a task reads alongside its
//! BAMs: opened once at plan time and matched to each input's references,
//! then read by every batch through a window that slides along its sorted
//! records, so each base is fetched about once per batch. A panel of SNVs
//! to count alleles at ([`SitePanel`]) is read whole at plan time.

use std::cell::RefCell;
use std::path::{Path, PathBuf};
use std::sync::{Arc, Mutex};

use omicsweave_bgzf::{BinIndex, locate_index};
use omicsweave_core::{Error, FileSource, Locus, ReadAt, Region, Result, SequenceRead};
use omicsweave_fasta::{Fai, FastaReader, Gzi, TWOBIT_MAGIC, TwoBitReader, fai_path, gzi_path};
use omicsweave_gtf::{Aliases, match_names};
use omicsweave_vcf::{BcfReader, Record as BcfRecord, gt_allele, gt_is_missing, gt_is_phased};

use crate::model::{FeatureList, Site};

use crate::assign::{BatchSpaces, Query, SpaceId};
use crate::prelude::{Context, Hits, ReadView};
use crate::source::Input;

/// Bases fetched beyond what a read asks for, so that the reads after it
/// are served from the same fetch.
const AHEAD: u32 = 256 << 10;

/// A sequence reader of either kind.
enum Sequences {
    Fasta(FastaReader<FileSource>),
    TwoBit(TwoBitReader<FileSource>),
}

impl Sequences {
    fn try_clone(&self) -> Result<Sequences> {
        Ok(match self {
            Sequences::Fasta(r) => Sequences::Fasta(r.try_clone()?),
            Sequences::TwoBit(r) => Sequences::TwoBit(r.try_clone()?),
        })
    }

    fn fetch(&mut self, chrom: &str, start: u32, end: u32, out: &mut [u8]) -> Result<()> {
        let locus = Locus::from(Region::new(chrom, start, end));
        match self {
            Sequences::Fasta(r) => r.sequence(&locus, out),
            Sequences::TwoBit(r) => r.sequence(&locus, out),
        }
    }
}

/// A reference genome (bgzip or plain FASTA with its `.fai` / `.gzi`, or
/// 2bit) matched to the references of the inputs.
pub struct Reference {
    pub path: PathBuf,
    proto: Mutex<Sequences>,
    /// `chroms[input][tid]`: the reference's sequence for the BAM's
    /// reference, with its length.
    chroms: Vec<Vec<Option<(String, u32)>>>,
}

/// The index of `path` at `at`, which must exist.
fn index<T>(at: PathBuf, path: &Path, parse: fn(&[u8]) -> Result<T>) -> Result<T> {
    match std::fs::read(&at) {
        Ok(bytes) => parse(&bytes),
        Err(e) if e.kind() == std::io::ErrorKind::NotFound => Err(Error::Io {
            kind: e.kind(),
            detail: format!(
                "the reference {} has no index {}; write it with \
                 `ow convert {} ref.fa.gz` (bgzip FASTA with .fai and .gzi)",
                path.display(),
                at.display(),
                path.display()
            ),
        }),
        Err(e) => Err(e.into()),
    }
}

impl Reference {
    /// Opens `path` and matches its sequences to every input's references
    /// (by name, `aliases` and the built-in `chr1` ↔ `1` rules); a matched
    /// sequence of another length than the BAM's reference is an error.
    pub fn open(path: &Path, inputs: &[Input], aliases: Option<&Aliases>) -> Result<Reference> {
        let source = FileSource::open(path).map_err(|e| match e {
            Error::Io { kind, detail } => Error::Io {
                kind,
                detail: format!("reference {}: {detail}", path.display()),
            },
            other => other,
        })?;
        let mut head = [0u8; 4];
        let n = (source.len() as usize).min(4);
        source.read_exact_at(&mut head[..n], 0)?;
        let magic = u32::from_le_bytes(head);
        let (reader, names): (Sequences, Vec<(String, u32)>) =
            if n == 4 && (magic == TWOBIT_MAGIC || magic.swap_bytes() == TWOBIT_MAGIC) {
                let r = TwoBitReader::open(source, false)?;
                let names: Vec<String> = r.names().map(str::to_owned).collect();
                let names = names
                    .into_iter()
                    .map(|c| {
                        let len = r.chrom_len(&c)?;
                        Ok((c, len))
                    })
                    .collect::<Result<Vec<_>>>()?;
                (Sequences::TwoBit(r), names)
            } else {
                let fai = index(fai_path(path), path, Fai::parse)?;
                let gzi = if head[..2] == [0x1f, 0x8b] {
                    Some(index(gzi_path(path), path, Gzi::parse)?)
                } else {
                    None
                };
                let names = fai
                    .entries()
                    .iter()
                    .map(|e| (e.name.clone(), e.len))
                    .collect();
                (Sequences::Fasta(FastaReader::new(source, fai, gzi)), names)
            };
        let own: Vec<&str> = names.iter().map(|n| n.0.as_str()).collect();
        let mut chroms = Vec::with_capacity(inputs.len());
        for input in inputs {
            let refs: Vec<&str> = input.refs.iter().map(|r| r.0.as_str()).collect();
            let matched = match_names(&own, &refs, aliases);
            let mut per = Vec::with_capacity(refs.len());
            for (tid, m) in matched.into_iter().enumerate() {
                per.push(match m {
                    None => None,
                    Some(i) => {
                        let (name, len) = &names[i as usize];
                        let want = input.refs[tid].1;
                        if *len != want {
                            return Err(Error::invalid(
                                "reference",
                                format!(
                                    "{} of {} has {len} bases, but {} of {} has {want}: \
                                     the BAM was aligned to another genome",
                                    name,
                                    path.display(),
                                    input.refs[tid].0,
                                    input.path.display()
                                ),
                            ));
                        }
                        Some((name.clone(), *len))
                    }
                });
            }
            chroms.push(per);
        }
        Ok(Reference {
            path: path.to_owned(),
            proto: Mutex::new(reader),
            chroms,
        })
    }

    fn reader(&self) -> Result<Sequences> {
        self.proto
            .lock()
            .unwrap_or_else(|e| e.into_inner())
            .try_clone()
    }

    /// The reference's sequence for reference `tid` of input `input`.
    pub fn chrom(&self, input: usize, tid: u32) -> Option<&(String, u32)> {
        self.chroms.get(input)?.get(tid as usize)?.as_ref()
    }

    /// A reader of single positions, for lookups outside the pass.
    pub fn lookup(&self, input: usize) -> Result<Lookup<'_>> {
        Ok(Lookup {
            reference: self,
            input,
            reader: self.reader()?,
        })
    }
}

/// Reads a few bases at a time from a [`Reference`].
pub struct Lookup<'a> {
    reference: &'a Reference,
    input: usize,
    reader: Sequences,
}

impl Lookup<'_> {
    /// The uppercase bases of `[start, end)` on `tid` (`N` past its end),
    /// or `None` without a sequence for `tid`.
    pub fn bases(&mut self, tid: u32, start: u32, end: u32) -> Result<Option<Vec<u8>>> {
        let Some((name, len)) = self.reference.chrom(self.input, tid) else {
            return Ok(None);
        };
        let mut out = vec![b'N'; end.saturating_sub(start) as usize];
        let hi = end.min(*len);
        if start < hi {
            self.reader
                .fetch(name, start, hi, &mut out[..(hi - start) as usize])?;
        }
        out.make_ascii_uppercase();
        Ok(Some(out))
    }
}

/// The reference bases of one batch: a window over one reference at a
/// time that moves forward with the records.
pub(crate) struct RefWindow<'a> {
    reference: &'a Reference,
    input: usize,
    reader: Option<Sequences>,
    tid: u32,
    start: u32,
    bases: Vec<u8>,
    /// One bounded disjoint window keeps alternating splice arms from
    /// replacing each other's read-ahead on every overlapping read.
    previous: Option<(u32, u32, Vec<u8>)>,
    #[cfg(test)]
    fetched: u64,
}

impl<'a> RefWindow<'a> {
    pub fn new(reference: &'a Reference, input: usize) -> RefWindow<'a> {
        RefWindow {
            reference,
            input,
            reader: None,
            tid: u32::MAX,
            start: 0,
            bases: Vec::new(),
            previous: None,
            #[cfg(test)]
            fetched: 0,
        }
    }

    fn get(&mut self, tid: u32, start: u32, end: u32, with: &mut dyn FnMut(&[u8])) -> Result<bool> {
        let Some((name, len)) = self.reference.chrom(self.input, tid) else {
            return Ok(false);
        };
        let have = self.start + self.bases.len() as u32;
        if (tid != self.tid || start < self.start || end > have)
            && let Some((other_tid, other_start, bases)) = &mut self.previous
            && *other_tid == tid
            && start >= *other_start
            && start <= *other_start + bases.len() as u32
        {
            std::mem::swap(&mut self.tid, other_tid);
            std::mem::swap(&mut self.start, other_start);
            std::mem::swap(&mut self.bases, bases);
            if bases.capacity() > (2 * AHEAD) as usize {
                self.previous = None;
            }
        }
        let have = self.start + self.bases.len() as u32;
        if tid != self.tid || start < self.start || start > have {
            if !self.bases.is_empty() && self.bases.capacity() <= (2 * AHEAD) as usize {
                let mut fresh = self.previous.take().map(|(_, _, b)| b).unwrap_or_default();
                fresh.clear();
                self.previous = Some((
                    self.tid,
                    self.start,
                    std::mem::replace(&mut self.bases, fresh),
                ));
            } else {
                self.bases.clear();
            }
            self.tid = tid;
            self.start = start;
        } else if start - self.start > AHEAD {
            self.bases.drain(..(start - self.start) as usize);
            self.start = start;
        }
        let have = self.start + self.bases.len() as u32;
        if end > have {
            let to = end.max(have.saturating_add(AHEAD)).min(*len).max(have);
            if to > have {
                if self.reader.is_none() {
                    self.reader = Some(self.reference.reader()?);
                }
                let at = self.bases.len();
                self.bases.resize(at + (to - have) as usize, 0);
                if let Err(e) = self.reader.as_mut().expect("opened above").fetch(
                    name,
                    have,
                    to,
                    &mut self.bases[at..],
                ) {
                    // Only successfully fetched bytes may satisfy a
                    // later query or move into the reusable window.
                    self.bases.truncate(at);
                    return Err(e);
                }
                #[cfg(test)]
                {
                    self.fetched += (to - have) as u64;
                }
                self.bases[at..].make_ascii_uppercase();
            }
        }
        let lo = (start - self.start) as usize;
        let hi = ((end.min(*len)).max(start) - self.start) as usize;
        if end <= *len {
            with(&self.bases[lo..hi]);
        } else {
            let mut padded = self.bases[lo..hi].to_vec();
            padded.resize((end - start) as usize, b'N');
            with(&padded);
        }
        Ok(true)
    }
}

/// Known variant sites (a BGZF BCF with its CSI): every reference base a
/// record's REF covers is masked.
pub struct Variants {
    pub path: PathBuf,
    proto: Mutex<BcfReader<FileSource>>,
    /// `contigs[input][tid]`: the BCF's contig for the BAM's reference.
    contigs: Vec<Vec<Option<String>>>,
}

impl Variants {
    pub fn open(path: &Path, inputs: &[Input], aliases: Option<&Aliases>) -> Result<Variants> {
        let at = locate_index(path).ok_or_else(|| Error::Io {
            kind: std::io::ErrorKind::NotFound,
            detail: format!(
                "no CSI index of {} found (.csi beside it); write it with \
                 `ow convert in.vcf.gz out.bcf --index`",
                path.display()
            ),
        })?;
        let reader = BcfReader::open(FileSource::open(path)?, BinIndex::read_file(&at)?)?;
        let names: Vec<String> = reader
            .header()
            .contigs
            .iter()
            .map(|c| c.name.clone())
            .collect();
        let own: Vec<&str> = names.iter().map(String::as_str).collect();
        let contigs = inputs
            .iter()
            .map(|input| {
                let refs: Vec<&str> = input.refs.iter().map(|r| r.0.as_str()).collect();
                match_names(&own, &refs, aliases)
                    .into_iter()
                    .map(|m| m.map(|i| names[i as usize].clone()))
                    .collect()
            })
            .collect();
        Ok(Variants {
            path: path.to_owned(),
            proto: Mutex::new(reader),
            contigs,
        })
    }

    fn reader(&self) -> Result<BcfReader<FileSource>> {
        self.proto
            .lock()
            .unwrap_or_else(|e| e.into_inner())
            .try_clone()
    }

    fn contig(&self, input: usize, tid: u32) -> Option<&str> {
        self.contigs.get(input)?.get(tid as usize)?.as_deref()
    }
}

/// The masked positions of one batch, a window moving forward like
/// [`RefWindow`].
pub(crate) struct MaskWindow<'a> {
    variants: &'a Variants,
    input: usize,
    reader: Option<BcfReader<FileSource>>,
    rec: BcfRecord,
    tid: u32,
    /// Positions are known for `[start, covered)`.
    start: u32,
    covered: u32,
    positions: Vec<u32>,
}

impl<'a> MaskWindow<'a> {
    pub fn new(variants: &'a Variants, input: usize) -> MaskWindow<'a> {
        MaskWindow {
            variants,
            input,
            reader: None,
            rec: BcfRecord::default(),
            tid: u32::MAX,
            start: 0,
            covered: 0,
            positions: Vec::new(),
        }
    }

    fn get(&mut self, tid: u32, start: u32, end: u32, out: &mut Vec<u32>) -> Result<()> {
        let Some(contig) = self.variants.contig(self.input, tid) else {
            return Ok(());
        };
        if tid != self.tid || start < self.start || start > self.covered {
            self.tid = tid;
            self.start = start;
            self.covered = start;
            self.positions.clear();
        } else if start - self.start > AHEAD {
            let keep = self.positions.partition_point(|&p| p < start);
            self.positions.drain(..keep);
            self.start = start;
        }
        if end > self.covered {
            let to = end.max(self.covered.saturating_add(AHEAD));
            if self.reader.is_none() {
                self.reader = Some(self.variants.reader()?);
            }
            let from = self.covered;
            let positions = &mut self.positions;
            let at = positions.len();
            if let Err(e) = self.reader.as_mut().expect("opened above").visit(
                Region::new(contig, from, to),
                &mut self.rec,
                |r| {
                    let s = r.pos.max(0) as u32;
                    let e = (r.pos + r.rlen.max(1)).clamp(0, u32::MAX as i64) as u32;
                    positions.extend(s.max(from)..e.min(to));
                    Ok(())
                },
            ) {
                positions.truncate(at);
                return Err(e);
            }
            positions.sort_unstable();
            positions.dedup();
            self.covered = to;
        }
        let lo = self.positions.partition_point(|&p| p < start);
        let hi = self.positions.partition_point(|&p| p < end);
        out.extend_from_slice(&self.positions[lo..hi]);
        Ok(())
    }
}

/// What a batch's records consult: the sweep cursors of the feature spaces,
/// and, when the task has them and some product reads them, the reference
/// and variant windows.
pub(crate) struct BatchCtx<'a> {
    pub spaces: BatchSpaces<'a>,
    pub reference: Option<RefCell<RefWindow<'a>>>,
    pub variants: Option<RefCell<MaskWindow<'a>>>,
}

impl Context for BatchCtx<'_> {
    fn assign(&self, id: SpaceId, r: &ReadView<'_>, query: Query, strand: Option<bool>) -> Hits {
        self.spaces.assign(id, r, query, strand)
    }

    fn assign_at(&self, id: SpaceId, r: &ReadView<'_>, start: u32, end: u32) -> Hits {
        self.spaces.assign_at(id, r, start, end)
    }

    fn reference(
        &self,
        tid: u32,
        start: u32,
        end: u32,
        with: &mut dyn FnMut(&[u8]),
    ) -> Result<bool> {
        match &self.reference {
            Some(w) => w.borrow_mut().get(tid, start, end, with),
            None => Ok(false),
        }
    }

    fn variants(&self, tid: u32, start: u32, end: u32, out: &mut Vec<u32>) -> Result<()> {
        match &self.variants {
            Some(w) => w.borrow_mut().get(tid, start, end, out),
            None => Ok(()),
        }
    }
}

/// Biallelic SNVs of a BCF (with its CSI), every site held in memory in
/// the file's order (contig, then position), with the phasing of one
/// sample: the columns of an allele matrix.
pub struct SitePanel {
    pub path: PathBuf,
    pub contigs: Arc<Vec<String>>,
    pub sites: Arc<Vec<Site>>,
    /// The sites of each contig: `sites[ranges[c].0..ranges[c].1]`.
    ranges: Vec<(usize, usize)>,
    /// `map[input][tid]`: the panel's contig for the BAM's reference.
    map: Vec<Vec<Option<u32>>>,
    /// Records that are not biallelic SNVs (`A`, `C`, `G`, `T`).
    pub skipped: u64,
    /// The phasing sample and its phased heterozygous sites.
    pub sample: Option<String>,
    pub phased: u64,
}

impl SitePanel {
    /// Reads every site of `path`; with `phasing`, the genotypes of that
    /// sample (`Some("")`: the first) give each site's haplotype 1 allele.
    pub fn open(
        path: &Path,
        inputs: &[Input],
        aliases: Option<&Aliases>,
        phasing: Option<&str>,
    ) -> Result<SitePanel> {
        let at = locate_index(path).ok_or_else(|| Error::Io {
            kind: std::io::ErrorKind::NotFound,
            detail: format!(
                "no CSI index of {} found (.csi beside it); write it with \
                 `ow convert in.vcf.gz out.bcf --index`",
                path.display()
            ),
        })?;
        let mut reader = BcfReader::open(FileSource::open(path)?, BinIndex::read_file(&at)?)?;
        let header = reader.header().clone();
        let sample = match phasing {
            None => None,
            Some(name) => {
                let i = if name.is_empty() {
                    (!header.samples.is_empty()).then_some(0)
                } else {
                    header.samples.iter().position(|s| s == name)
                };
                match i {
                    Some(i) => Some(i),
                    None => {
                        return Err(Error::invalid(
                            "sites",
                            format!(
                                "{} has no sample {:?} to phase alleles by",
                                path.display(),
                                name
                            ),
                        ));
                    }
                }
            }
        };
        let contigs: Vec<String> = header.contigs.iter().map(|c| c.name.clone()).collect();
        let mut sites = Vec::new();
        let mut ranges = Vec::with_capacity(contigs.len());
        let mut skipped = 0u64;
        let mut phased = 0u64;
        let mut rec = BcfRecord::default();
        let base = |a: &[u8]| -> Option<u8> {
            match a {
                [b] if matches!(b.to_ascii_uppercase(), b'A' | b'C' | b'G' | b'T') => {
                    Some(b.to_ascii_uppercase())
                }
                _ => None,
            }
        };
        for (c, name) in contigs.iter().enumerate() {
            let start = sites.len();
            reader.visit(Region::new(name, 0, u32::MAX), &mut rec, |r| {
                let alleles = (r.n_alleles() == 2).then(|| (base(r.allele(0)), base(r.allele(1))));
                let (Some((Some(reference), Some(alt))), true) =
                    (alleles, r.pos >= 0 && r.pos < u32::MAX as i64)
                else {
                    skipped += 1;
                    return Ok(());
                };
                let hap1 = sample.and_then(|i| {
                    let g = r.genotypes()?;
                    let v = g.sample(i);
                    if v.len() != 2 || gt_is_missing(v[0]) || gt_is_missing(v[1]) {
                        return None;
                    }
                    let (a, b) = (gt_allele(v[0]), gt_allele(v[1]));
                    let het = a != b && (0..=1).contains(&a) && (0..=1).contains(&b);
                    (het && gt_is_phased(v[1])).then_some(a as u8)
                });
                phased += hap1.is_some() as u64;
                sites.push(Site {
                    contig: c as u32,
                    pos: r.pos as u32,
                    reference,
                    alt,
                    hap1,
                });
                Ok(())
            })?;
            ranges.push((start, sites.len()));
        }
        let own: Vec<&str> = contigs.iter().map(String::as_str).collect();
        let map = inputs
            .iter()
            .map(|input| {
                let refs: Vec<&str> = input.refs.iter().map(|r| r.0.as_str()).collect();
                match_names(&own, &refs, aliases)
            })
            .collect();
        Ok(SitePanel {
            path: path.to_owned(),
            sample: sample.map(|i| header.samples[i].clone()),
            contigs: Arc::new(contigs),
            sites: Arc::new(sites),
            ranges,
            map,
            skipped,
            phased,
        })
    }

    /// The sites on reference `tid` of input `input`, with the column of
    /// the first.
    pub fn sites_of(&self, input: usize, tid: u32) -> (usize, &[Site]) {
        match self
            .map
            .get(input)
            .and_then(|m| m.get(tid as usize))
            .copied()
            .flatten()
        {
            Some(c) => {
                let (s, e) = self.ranges[c as usize];
                (s, &self.sites[s..e])
            }
            None => (0, &[]),
        }
    }

    pub fn heap_bytes(&self) -> u64 {
        (self.sites.len() * std::mem::size_of::<Site>()
            + self.contigs.iter().map(|c| c.len() + 24).sum::<usize>()
            + self.ranges.len() * 16) as u64
    }

    /// The sites as matrix columns.
    pub fn features(&self) -> SiteFeatures {
        SiteFeatures {
            contigs: Arc::clone(&self.contigs),
            sites: Arc::clone(&self.sites),
        }
    }
}

/// Sites as the columns of a matrix: `chrom:pos:REF>ALT` (1-based), type
/// `Variant`, spanning their base.
pub struct SiteFeatures {
    contigs: Arc<Vec<String>>,
    sites: Arc<Vec<Site>>,
}

impl FeatureList for SiteFeatures {
    fn len(&self) -> usize {
        self.sites.len()
    }

    fn id(&self, i: usize, out: &mut Vec<u8>) {
        let s = &self.sites[i];
        out.extend_from_slice(self.contigs[s.contig as usize].as_bytes());
        out.push(b':');
        crate::sink::push_u64(out, s.pos as u64 + 1);
        out.extend_from_slice(&[b':', s.reference, b'>', s.alt]);
    }

    fn line(&self, i: usize, out: &mut Vec<u8>) {
        let s = &self.sites[i];
        self.id(i, out);
        out.push(b'\t');
        self.id(i, out);
        out.extend_from_slice(b"\tVariant\t");
        out.extend_from_slice(self.contigs[s.contig as usize].as_bytes());
        out.push(b'\t');
        crate::sink::push_u64(out, s.pos as u64);
        out.push(b'\t');
        crate::sink::push_u64(out, s.pos as u64 + 1);
    }
}

#[cfg(test)]
#[path = "../tests/common/mod.rs"]
mod test_common;

#[cfg(test)]
mod reference_window_tests {
    use super::*;

    #[test]
    fn alternating_splice_arms_reuse_bounded_real_fasta_windows() {
        let dir = tempfile::tempdir().unwrap();
        let path = dir.path().join("reference.fa");
        let len = 1_500_000u32;
        let mut bytes = b">chr1\n".to_vec();
        bytes.extend((0..len).map(|i| b"ACGT"[(i % 4) as usize]));
        bytes.push(b'\n');
        std::fs::write(&path, &bytes).unwrap();
        std::fs::write(
            fai_path(&path),
            format!("chr1\t{len}\t6\t{len}\t{}\n", len + 1),
        )
        .unwrap();
        let input = Input {
            path: dir.path().join("input.bam"),
            refs: vec![("chr1".into(), len)],
            read_groups: vec![],
            len: 0,
        };
        let reference = Reference::open(&path, &[input], None).unwrap();
        let mut window = RefWindow::new(&reference, 0);
        for offset in 0..1000 {
            for s in [100 + offset, 1_000_100 + offset] {
                assert!(
                    window
                        .get(0, s, s + 64, &mut |b| {
                            let expected: Vec<_> =
                                (s..s + 64).map(|p| b"ACGT"[(p % 4) as usize]).collect();
                            assert_eq!(b, expected);
                        })
                        .unwrap()
                );
            }
        }
        // Counts bytes actually passed to the real FastaReader fetch,
        // including read-ahead, rather than just callback payload sizes.
        assert_eq!(window.fetched, 2 * AHEAD as u64);
        let spare = window.previous.as_ref().unwrap();
        assert!(window.bases.capacity() + spare.2.capacity() <= (4 * AHEAD) as usize);

        let successful = window.bases.len();
        let s = window.start + AHEAD - 4;
        std::fs::write(&path, b">chr1\n").unwrap();
        assert!(window.get(0, s, s + 64, &mut |_| {}).is_err());
        assert_eq!(window.bases.len(), successful);
        let fetched = window.fetched;
        window
            .get(0, 1_000_100, 1_000_164, &mut |b| {
                assert_eq!(&b[..4], b"ACGT")
            })
            .unwrap();
        assert_eq!(window.fetched, fetched);
        std::fs::write(&path, &bytes).unwrap();
        window
            .get(0, s, s + 64, &mut |b| {
                let expected: Vec<_> = (s..s + 64).map(|p| b"ACGT"[(p % 4) as usize]).collect();
                assert_eq!(b, expected);
            })
            .unwrap();
        assert!(window.fetched > fetched);

        window.get(0, 0, len, &mut |_| {}).unwrap();
        assert!(window.bases.capacity() > (2 * AHEAD) as usize);
        window
            .get(0, len + 100, len + 110, &mut |b| {
                assert_eq!(b, b"NNNNNNNNNN")
            })
            .unwrap();
        window.get(0, 1_000_100, 1_000_164, &mut |_| {}).unwrap();
        assert!(window.previous.is_none());
        assert!(window.bases.capacity() <= (2 * AHEAD) as usize);
    }

    #[test]
    fn failed_variant_extension_keeps_sorted_prefix_and_retries() {
        let dir = tempfile::tempdir().unwrap();
        let refs = vec![("chr1".into(), 1_000_000)];
        let sites: Vec<_> = (0..5000).map(|i| (0, 100 + i * 100, b'A', b'C')).collect();
        let path = super::test_common::write_bcf(dir.path(), "sites", &refs, &sites);
        let input = Input {
            path: dir.path().join("input.bam"),
            refs,
            read_groups: vec![],
            len: 0,
        };
        let variants = Variants::open(&path, &[input], None).unwrap();
        let mut window = MaskWindow::new(&variants, 0);
        window.get(0, 0, 64, &mut Vec::new()).unwrap();
        let prefix = window.positions.clone();
        let covered = window.covered;
        assert!(!prefix.is_empty());
        let bytes = std::fs::read(&path).unwrap();
        let mut corrupt = bytes.clone();
        let (mut at, mut last) = (0, 0);
        while at < bytes.len() {
            let len = u16::from_le_bytes([bytes[at + 16], bytes[at + 17]]) as usize + 1;
            if bytes[at + len - 4..at + len] != [0, 0, 0, 0] {
                last = at;
            }
            at += len;
        }
        let len = u16::from_le_bytes([bytes[last + 16], bytes[last + 17]]) as usize + 1;
        corrupt[last + len - 8] ^= 1;
        std::fs::write(&path, corrupt).unwrap();
        let mut diagnostic = variants.reader().unwrap();
        let mut visited = 0;
        assert!(
            diagnostic
                .visit(
                    Region::new("chr1", covered, 600_000),
                    &mut BcfRecord::default(),
                    |_| {
                        visited += 1;
                        Ok(())
                    }
                )
                .is_err()
        );
        assert!(visited > 0);
        assert!(window.get(0, 0, 600_000, &mut Vec::new()).is_err());
        assert_eq!(window.positions, prefix);
        assert_eq!(window.covered, covered);
        std::fs::write(&path, &bytes).unwrap();
        let mut out = Vec::new();
        window.get(0, 0, 600_000, &mut out).unwrap();
        let expected: Vec<_> = sites.iter().map(|s| s.1).collect();
        assert_eq!(out, expected);
    }
}
