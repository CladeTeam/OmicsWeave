//! [`BcfReader`]: a BGZF-compressed BCF file with its CSI index, read by
//! region as htslib's `bcf_itr_queryi` reads it, and the genotype matrix
//! of a region ([`GenotypeMatrix`]).

use std::sync::Arc;

use omicsweave_bgzf::{BinIndex, ChunkWalk, IndexKind, SeekReader, SeekStats};
use omicsweave_core::cache::{BlockCache, BlockKey, SourceToken};
use omicsweave_core::{Error, FileSource, ReadAt, Region, Result};

use crate::bcf::{Fixed, parse_record_fields};
use crate::fields::{FormatField, Typed, Value, type_width};
use crate::header::{Header, VCF44, ValueType};
use crate::record::{GT_MISSING_VALUE, GT_VECTOR_END, Record, gt_allele, gt_is_phased};
use crate::table::{RecordFilter, RecordQuery, TableBuilder, VariantTable};

/// A BCF file and its CSI index. Opening reads the header; each query reads
/// only the chunks the index gives for its region.
pub struct BcfReader<S> {
    source: S,
    header: Arc<Header>,
    index: Arc<BinIndex>,
    stream: SeekReader,
}

impl<S: ReadAt> BcfReader<S> {
    /// Reads the header of `source`, a BGZF BCF file that `index` (a CSI)
    /// indexes.
    pub fn open(source: S, index: BinIndex) -> Result<BcfReader<S>> {
        if index.kind() != IndexKind::Csi {
            return Err(Error::invalid("bcf index", "a BCF is indexed by a CSI"));
        }
        let mut stream = SeekReader::new();
        stream.seek(&source, 0, 0)?;
        let magic = stream.take(&source, 5, "bcf header")?;
        let magic = stream.slice(magic);
        if magic[..4] != *b"BCF\x02" || !matches!(magic[4], 1 | 2) {
            return Err(Error::invalid(
                "bcf header",
                "not a BCF 2 file (an indexed BCF is BGZF-compressed)",
            ));
        }
        let len = stream.take(&source, 4, "bcf header")?;
        let len = u32::from_le_bytes(stream.slice(len).try_into().unwrap()) as usize;
        let text = stream.take(&source, len, "bcf header")?;
        let text = stream.slice(text);
        if text.last() != Some(&0) {
            return Err(Error::invalid("bcf header", "missing terminating NUL"));
        }
        let text = std::str::from_utf8(text)
            .map_err(|_| Error::invalid("bcf header", "header is not UTF-8"))?;
        let header = Header::parse(text)?;
        Ok(BcfReader {
            source,
            header: Arc::new(header),
            index: Arc::new(index),
            stream,
        })
    }

    pub fn header(&self) -> &Header {
        &self.header
    }

    /// The header, as shared with every context cloned from this reader.
    pub fn shared_header(&self) -> &Arc<Header> {
        &self.header
    }

    pub fn index(&self) -> &BinIndex {
        &self.index
    }

    /// Attaches externally managed BGZF blocks. Source contents and indexes
    /// must remain unchanged for this reader and its shared contexts.
    pub fn set_cache(&mut self, cache: Option<Arc<BlockCache>>) {
        let namespace = self.stream.cache_namespace();
        self.stream.set_cache(cache, namespace);
    }

    pub fn cache_namespace(&self) -> SourceToken {
        self.stream.cache_namespace()
    }

    /// Actual positional reads and BGZF inflations in this context.
    pub fn stream_stats(&self) -> SeekStats {
        self.stream.stats()
    }

    pub fn reset_stats(&mut self) {
        self.stream.reset_stats();
    }

    /// BGZF dependencies of the CSI chunks; no record decoding is performed.
    pub fn describe_cache(&mut self, region: Region<'_>, out: &mut Vec<BlockKey>) -> Result<()> {
        region.check()?;
        let rid = self
            .header
            .contig_id(region.chrom)
            .ok_or_else(|| Error::UnknownChrom(region.chrom.to_owned()))?;
        if region.is_empty() {
            out.clear();
            return Ok(());
        }
        let chunks = self
            .index
            .chunks(rid as usize, region.start as i64, region.end as i64);
        self.stream.describe_cache(&self.source, &chunks, out)
    }

    /// Populates only external block storage, without constructing records or
    /// a genotype matrix. The caller chooses when to prefetch and evict.
    pub fn prefetch_cache(&mut self, region: Region<'_>) -> Result<()> {
        region.check()?;
        let rid = self
            .header
            .contig_id(region.chrom)
            .ok_or_else(|| Error::UnknownChrom(region.chrom.to_owned()))?;
        if region.is_empty() {
            return self.stream.prefetch_cache(&self.source, &[]);
        }
        let chunks = self
            .index
            .chunks(rid as usize, region.start as i64, region.end as i64);
        self.stream.prefetch_cache(&self.source, &chunks)
    }

    /// Shares an immutable source session and metadata, retaining independent
    /// decompression and query scratch in the new execution context.
    pub fn clone_context(&self) -> Self
    where
        S: Clone,
    {
        Self {
            source: self.source.clone(),
            header: Arc::clone(&self.header),
            index: Arc::clone(&self.index),
            stream: self.stream.clone_context(),
        }
    }

    /// Decodes into `rec` every record whose span `[pos, pos + rlen)`
    /// overlaps `region` and calls `visit` with it, in file order: the
    /// records htslib's `bcf_itr_queryi` iterator returns. The span is the
    /// stored `rlen`, as htslib's iterator reads it; `rec.rlen` holds the
    /// span htslib decodes. Reading stops at the first record of another
    /// contig or starting at or after the region's end. `rec` is cleared
    /// on an error or when the query yields no records.
    pub fn visit(
        &mut self,
        region: Region<'_>,
        rec: &mut Record,
        visit: impl FnMut(&Record) -> Result<()>,
    ) -> Result<()> {
        rec.clear();
        let mut visit = visit;
        let result = self.visit_fields::<true>(region, rec, |rec, _| visit(rec));
        if result.is_err() {
            rec.clear();
        }
        result
    }

    fn visit_fields<const UNPACK_FIELDS: bool>(
        &mut self,
        region: Region<'_>,
        rec: &mut Record,
        mut visit: impl FnMut(&Record, Option<FormatField<'_>>) -> Result<()>,
    ) -> Result<()> {
        self.stream.reset_query();
        region.check()?;
        let rid = self
            .header
            .contig_id(region.chrom)
            .ok_or_else(|| Error::UnknownChrom(region.chrom.to_owned()))?;
        if region.is_empty() {
            return Ok(());
        }
        let (beg, end) = (region.start as i64, region.end as i64);
        let mut walk = ChunkWalk::new(self.index.chunks(rid as usize, beg, end));
        let BcfReader {
            source,
            header,
            stream,
            ..
        } = self;
        while walk.next(stream, source)? {
            if stream.fill(source, 1)? == 0 {
                return Ok(());
            }
            let lens = stream.take(source, 8, "bcf record")?;
            let lens = stream.slice(lens);
            let l_shared = u32::from_le_bytes(lens[..4].try_into().unwrap()) as usize;
            let l_indiv = u32::from_le_bytes(lens[4..].try_into().unwrap()) as usize;
            if l_shared < 24 {
                return Err(Error::invalid("bcf record", "shared block too short"));
            }
            let body_len = l_shared
                .checked_add(l_indiv)
                .ok_or_else(|| Error::invalid("bcf record", "record length overflows usize"))?;
            let body = stream.take(source, body_len, "bcf record")?;
            let body = stream.slice(body);
            let fixed = Fixed::parse(body[..24].try_into().unwrap());
            let pos = if fixed.pos == u32::MAX {
                -1
            } else {
                fixed.pos as i64
            };
            if fixed.rid != rid as i32 || pos >= end {
                return Ok(());
            }
            if pos + fixed.rlen as i64 > beg && end > pos {
                let (shared, indiv) = body[24..].split_at(l_shared - 24);
                let gt = parse_record_fields::<UNPACK_FIELDS>(header, &fixed, shared, indiv, rec)?;
                visit(rec, gt)?;
            }
        }
        Ok(())
    }

    fn check_samples(&self, samples: &[usize]) -> Result<()> {
        let n_sample = self.header.samples.len();
        if let Some(&bad) = samples.iter().find(|&&s| s >= n_sample) {
            return Err(Error::invalid(
                "samples",
                format!("sample index {bad} is not below the {n_sample} samples"),
            ));
        }
        Ok(())
    }

    /// The records `visit` returns for `region` that `filter` keeps, as
    /// columns: the genotypes of `samples` (header sample indexes, in the
    /// order given), QUAL, FILTER and the fields of `query`; one pass,
    /// each record decoded once.
    pub fn records(
        &mut self,
        region: Region<'_>,
        samples: &[usize],
        query: &RecordQuery,
        mut filter: Option<&mut RecordFilter>,
    ) -> Result<VariantTable> {
        self.check_samples(samples)?;
        let header = Arc::clone(&self.header);
        let mut genotypes = Builder::new(samples);
        let mut table = TableBuilder::new(&header, query, samples)?;
        let mut rec = Record::default();
        // Unpacked records hold GT among their FORMAT fields, its phasing
        // already brought up to VCF 4.4 as htslib reads older files.
        let gt_key = header.key_id("GT").filter(|&id| {
            header
                .key(id)
                .is_some_and(|key| key.format == Some(ValueType::String))
        });
        let result = self.visit_fields::<true>(region, &mut rec, |rec, _| {
            if let Some(f) = filter.as_deref_mut()
                && !f.keeps(&header, rec)?
            {
                return Ok(());
            }
            genotypes.push(rec, gt_key.and_then(|k| rec.format_field(k)), false)?;
            table.push(&header, rec)
        });
        rec.clear();
        result?;
        table.finish(genotypes.finish()?)
    }

    /// The genotypes of `samples` (header sample indexes, in the order
    /// given) at every record `visit` returns for `region`.
    pub fn genotypes(&mut self, region: Region<'_>, samples: &[usize]) -> Result<GenotypeMatrix> {
        self.check_samples(samples)?;
        let mut builder = Builder::new(samples);
        let mut rec = Record::default();
        let update_phase = self.header.version < VCF44;
        self.visit_fields::<false>(region, &mut rec, |rec, gt| {
            builder.push(rec, gt, update_phase)
        })?;
        builder.finish()
    }
}

impl BcfReader<FileSource> {
    /// Clones the open file handle and shares its immutable header/index/cache
    /// session; no file data is reread to create the context.
    pub fn try_clone(&self) -> Result<Self> {
        Ok(Self {
            source: self.source.try_clone()?,
            header: Arc::clone(&self.header),
            index: Arc::clone(&self.index),
            stream: self.stream.clone_context(),
        })
    }
}

/// Columns of a string per record, in Arrow's `LargeString` layout.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Strings {
    /// `len + 1` offsets into `data`, starting at 0.
    pub offsets: Vec<i64>,
    pub data: Vec<u8>,
}

impl Default for Strings {
    fn default() -> Self {
        Self::new()
    }
}

impl Strings {
    fn new() -> Strings {
        Strings {
            offsets: vec![0],
            data: Vec::new(),
        }
    }

    pub(crate) fn end(&mut self) {
        self.offsets.push(self.data.len() as i64);
    }

    pub fn len(&self) -> usize {
        self.offsets.len().saturating_sub(1)
    }

    pub fn is_empty(&self) -> bool {
        self.len() == 0
    }

    pub fn get(&self, i: usize) -> &[u8] {
        &self.data[self.offsets[i] as usize..self.offsets[i + 1] as usize]
    }
}

/// Allele index of a genotype matrix entry for a missing allele (`.`).
pub const GT_MISSING: i8 = -1;
/// Entry past the end of a sample's alleles: a sample of lower ploidy than
/// the matrix, or a record without genotypes.
pub const GT_ABSENT: i8 = -2;

/// The variants of a region and the genotypes of chosen samples at them.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct GenotypeMatrix {
    /// 0-based position of each record's first REF base.
    pub start: Vec<i64>,
    /// `start` plus the reference span htslib decodes.
    pub end: Vec<i64>,
    pub id: Strings,
    pub reference: Strings,
    /// The ALT alleles, comma-separated; empty when there are none.
    pub alt: Strings,
    /// Samples per record.
    pub samples: usize,
    /// Alleles per sample: the most any record of the region gives.
    pub ploidy: usize,
    /// `[record][sample][allele]` allele indexes (0 = REF), with
    /// [`GT_MISSING`] and [`GT_ABSENT`].
    pub gt: Vec<i8>,
    /// `[record][sample]`: every allele after the first is phased; a sample
    /// of one allele takes that allele's phase, as htslib decodes it.
    pub phased: Vec<bool>,
}

impl GenotypeMatrix {
    pub fn len(&self) -> usize {
        self.start.len()
    }

    pub fn is_empty(&self) -> bool {
        self.start.is_empty()
    }
}

/// Collects a matrix one record at a time. Each record's genotypes are
/// stored at its own ploidy and laid out at the largest one once the region
/// ends; when every record has the same ploidy that layout is the stored one.
struct Builder<'a> {
    samples: &'a [usize],
    m: GenotypeMatrix,
    widths: Vec<usize>,
}

impl<'a> Builder<'a> {
    fn new(samples: &'a [usize]) -> Builder<'a> {
        Builder {
            samples,
            m: GenotypeMatrix {
                start: Vec::new(),
                end: Vec::new(),
                id: Strings::new(),
                reference: Strings::new(),
                alt: Strings::new(),
                samples: samples.len(),
                ploidy: 0,
                gt: Vec::new(),
                phased: Vec::new(),
            },
            widths: Vec::new(),
        }
    }

    fn push(
        &mut self,
        rec: &Record,
        gt: Option<FormatField<'_>>,
        update_phase: bool,
    ) -> Result<()> {
        let m = &mut self.m;
        m.start.push(rec.pos);
        m.end.push(
            rec.pos
                .checked_add(rec.rlen)
                .ok_or_else(|| Error::invalid("genotype matrix", "record end overflows i64"))?,
        );
        m.id.data.extend_from_slice(rec.id());
        m.id.end();
        m.reference.data.extend_from_slice(rec.allele(0));
        m.reference.end();
        for i in 1..rec.n_alleles() {
            if i > 1 {
                m.alt.data.push(b',');
            }
            m.alt.data.extend_from_slice(rec.allele(i));
        }
        m.alt.end();
        let Some(gts) = gt else {
            self.widths.push(0);
            m.phased
                .extend(std::iter::repeat_n(false, self.samples.len()));
            return Ok(());
        };
        let width = gts.n;
        self.widths.push(width);
        m.ploidy = m.ploidy.max(width);
        for &s in self.samples {
            let values = gts.sample(s);
            // Pre-4.4 normalization changes the first stored byte, including
            // sentinels; apply it before interpreting that selected value.
            let first = if update_phase && width > 0 && values.bytes[0] & 1 == 0 {
                let w = type_width(values.ty);
                let phase = if width == 1 {
                    u8::from(values.bytes[0] != 0)
                } else {
                    values
                        .bytes
                        .chunks_exact(w)
                        .skip(1)
                        .fold(1, |all, v| all & v[0])
                        & 1
                };
                let mut bytes = [0; 4];
                bytes[..w].copy_from_slice(&values.bytes[..w]);
                bytes[0] |= phase;
                Some(
                    Typed {
                        ty: values.ty,
                        bytes: &bytes[..w],
                    }
                    .value(0),
                )
            } else {
                None
            };
            let mut alleles = 0;
            let mut first_phased = false;
            let mut rest_phased = true;
            let mut ended = false;
            for j in 0..width {
                let v = if ended {
                    GT_VECTOR_END
                } else {
                    match first.filter(|_| j == 0).unwrap_or_else(|| values.value(j)) {
                        Value::Missing => GT_MISSING_VALUE,
                        Value::VectorEnd => {
                            ended = true;
                            GT_VECTOR_END
                        }
                        Value::Int(v) => i32::try_from(v)
                            .map_err(|_| Error::invalid("genotype", "GT integer exceeds int32"))?,
                        _ => return Err(Error::invalid("genotype", "GT is not an integer vector")),
                    }
                };
                let entry = match v {
                    GT_VECTOR_END => GT_ABSENT,
                    GT_MISSING_VALUE => GT_MISSING,
                    v => {
                        let a = gt_allele(v);
                        i8::try_from(a).map_err(|_| {
                            Error::invalid(
                                "genotype",
                                format!("allele index {a} does not fit the int8 matrix"),
                            )
                        })?
                    }
                };
                if entry != GT_ABSENT {
                    alleles += 1;
                    if j == 0 {
                        first_phased = gt_is_phased(v);
                    } else if !gt_is_phased(v) {
                        rest_phased = false;
                    }
                }
                m.gt.push(entry);
            }
            m.phased.push(match alleles {
                0 => false,
                1 => first_phased,
                _ => rest_phased,
            });
        }
        Ok(())
    }

    fn finish(mut self) -> Result<GenotypeMatrix> {
        let m = &mut self.m;
        let ploidy = m.ploidy;
        if self.widths.iter().all(|&w| w == ploidy) {
            return Ok(self.m);
        }
        let count = m
            .len()
            .checked_mul(m.samples)
            .and_then(|n| n.checked_mul(ploidy))
            .ok_or_else(|| Error::invalid("genotype matrix", "matrix length overflows usize"))?;
        let mut gt = Vec::new();
        gt.try_reserve_exact(count)
            .map_err(|e| Error::invalid("genotype matrix", e.to_string()))?;
        gt.resize(count, GT_ABSENT);
        let mut from = 0;
        for (r, &w) in self.widths.iter().enumerate() {
            for s in 0..m.samples {
                let to = (r * m.samples + s) * ploidy;
                gt[to..to + w].copy_from_slice(&m.gt[from..from + w]);
                from += w;
            }
        }
        m.gt = gt;
        Ok(self.m)
    }
}
