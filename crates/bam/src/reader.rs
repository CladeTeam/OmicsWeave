//! [`BamReader`]: a BAM file with its index, read by region.

use std::ops::Range;
use std::sync::Arc;

use omicsweave_bgzf::{BinIndex, ChunkWalk, SeekReader};
use omicsweave_core::cache::{BlockCache, BlockKey, SourceToken};
use omicsweave_core::{
    AlignmentEntry, AlignmentFields, AlignmentRead, AlignmentSize, AlignmentsOut, Error,
    ModificationEntry, ReadAt, Region, Result, StringColumn,
};

use crate::header::BamHeader;
use crate::mods::{ModDecoder, ModSeq, decode_record};
use crate::record::{Record, unpack_seq};
use crate::stream::{Reading, next_record};

/// A BAM file and its BAI or CSI index. Opening reads the header; each
/// query reads only the chunks the index gives for its region.
pub struct BamReader<S> {
    source: S,
    header: Arc<BamHeader>,
    index: Arc<BinIndex>,
    stream: SeekReader,
    staged: Staged,
}

/// What a query keeps for [`AlignmentRead::write`]: the fixed fields, and
/// the record bytes behind the variable fields asked for.
#[derive(Default)]
struct Staged {
    queried: bool,
    fields: AlignmentFields,
    size: AlignmentSize,
    entries: Vec<AlignmentEntry>,
    bytes: Vec<u8>,
    records: Vec<Kept>,
    decoder: ModDecoder,
    mods: Vec<ModificationEntry>,
    mod_ends: Vec<i64>,
}

/// Where a kept record's fields lie in [`Staged::bytes`].
struct Kept {
    // Selected fields are contiguous in NAME, CIGAR, SEQ, QUAL, AUX order;
    // the previous record's aux end is this record's name start.
    name_end: usize,
    cigar_end: usize,
    seq_end: usize,
    l_seq: usize,
    qual_end: usize,
    aux_end: usize,
}

impl<S: ReadAt> BamReader<S> {
    /// Reads the header of `source`, a BAM file that `index` indexes.
    pub fn open(source: S, index: BinIndex) -> Result<BamReader<S>> {
        let mut stream = SeekReader::new();
        stream.seek(&source, 0, 0)?;
        let header = BamHeader::read(&mut Reading {
            stream: &mut stream,
            source: &source,
        })?;
        Ok(BamReader {
            source,
            header: Arc::new(header),
            index: Arc::new(index),
            stream,
            staged: Staged::default(),
        })
    }

    pub fn header(&self) -> &BamHeader {
        &self.header
    }

    pub fn index(&self) -> &BinIndex {
        &self.index
    }

    pub fn source(&self) -> &S {
        &self.source
    }

    /// Attaches caller-owned BGZF blocks. The source and index must remain
    /// unchanged for this source session and all its execution contexts.
    pub fn set_cache(&mut self, cache: Option<Arc<BlockCache>>) {
        let namespace = self.stream.cache_namespace();
        self.stream.set_cache(cache, namespace);
    }

    pub fn cache_namespace(&self) -> SourceToken {
        self.stream.cache_namespace()
    }

    fn chunks(&self, region: Region<'_>) -> Result<Vec<omicsweave_bgzf::Chunk>> {
        region.check()?;
        let tid = self
            .header
            .tid(region.chrom)
            .ok_or_else(|| Error::UnknownChrom(region.chrom.to_owned()))?;
        Ok(if region.is_empty() {
            Vec::new()
        } else {
            self.index
                .chunks(tid, region.start as i64, region.end as i64)
        })
    }

    /// Describes the complete physical blocks behind the indexed chunks,
    /// without inflating or constructing alignment output.
    pub fn describe_cache(&mut self, region: Region<'_>, out: &mut Vec<BlockKey>) -> Result<()> {
        let chunks = self.chunks(region)?;
        self.stream.describe_cache(&self.source, &chunks, out)
    }

    /// Loads indexed dependencies without constructing alignment output.
    pub fn prefetch_cache(&mut self, region: Region<'_>) -> Result<()> {
        let chunks = self.chunks(region)?;
        self.stream.prefetch_cache(&self.source, &chunks)
    }

    /// Physical input and decoding work, independent of output formatting.
    pub fn stats(&self) -> omicsweave_bgzf::SeekStats {
        self.stream.stats()
    }

    pub fn reset_stats(&mut self) {
        self.stream.reset_stats();
    }

    fn context_with_source(&self, source: S) -> Self {
        Self {
            source,
            header: Arc::clone(&self.header),
            index: Arc::clone(&self.index),
            stream: self.stream.clone_context(),
            staged: Staged::default(),
        }
    }

    /// Shares the source session, index and external cache; execution scratch
    /// and staged results start empty and are independent.
    pub fn clone_context(&self) -> Self
    where
        S: Clone,
    {
        self.context_with_source(self.source.clone())
    }

    /// Calls `visit` with every record whose reference span `[pos,
    /// bam_endpos)` overlaps `region`, in file order: the records htslib's
    /// `sam_itr_queryi` iterator returns. Reading stops at the first record
    /// of another reference or starting at or after the region's end.
    pub fn visit(
        &mut self,
        region: Region<'_>,
        mut visit: impl FnMut(&Record<'_>) -> Result<()>,
    ) -> Result<()> {
        region.check()?;
        let tid = self
            .header
            .tid(region.chrom)
            .ok_or_else(|| Error::UnknownChrom(region.chrom.to_owned()))?;
        if region.is_empty() {
            return Ok(());
        }
        let (beg, end) = (region.start as i64, region.end as i64);
        let mut walk = ChunkWalk::new(self.index.chunks(tid, beg, end));
        let n_refs = self.header.references().len();
        let BamReader { source, stream, .. } = self;
        stream.reset_query();
        while walk.next(stream, source)? {
            let Some(range) = next_record(stream, source)? else {
                return Ok(());
            };
            let rec = Record::parse(stream.slice(range), n_refs)?;
            if rec.tid as i64 != tid as i64 || rec.pos as i64 >= end {
                return Ok(());
            }
            if rec.end() > beg && end > rec.pos as i64 {
                visit(&rec)?;
            }
        }
        Ok(())
    }
}

impl BamReader<omicsweave_core::FileSource> {
    pub fn try_clone(&self) -> Result<Self> {
        Ok(self.context_with_source(self.source.try_clone()?))
    }
}

/// Copies one kept part of every record into a string column.
fn fill_strings(
    column: &mut StringColumn<'_>,
    bytes: &[u8],
    records: &[Kept],
    part: impl Fn(usize, &Kept) -> Range<usize>,
) {
    column.offsets[0] = 0;
    let mut at = 0;
    let mut start = 0;
    for (i, k) in records.iter().enumerate() {
        let p = &bytes[part(start, k)];
        column.data[at..at + p.len()].copy_from_slice(p);
        at += p.len();
        column.offsets[i + 1] = at as i64;
        start = k.aux_end;
    }
}

impl<S: ReadAt> AlignmentRead for BamReader<S> {
    fn chrom_len(&self, chrom: &str) -> Result<u32> {
        self.header.chrom_len(chrom)
    }

    fn query_filtered(
        &mut self,
        region: Region<'_>,
        fields: AlignmentFields,
        modification_codes: &[i32],
    ) -> Result<AlignmentSize> {
        let mut staged = std::mem::take(&mut self.staged);
        staged.queried = false;
        staged.fields = fields;
        staged.entries.clear();
        staged.bytes.clear();
        staged.records.clear();
        staged.mods.clear();
        staged.mod_ends.clear();
        let want_mods = fields.contains(AlignmentFields::MODIFICATIONS);
        let want_record_fields = fields.bits() & !AlignmentFields::MODIFICATIONS.bits() != 0;
        let mut size = AlignmentSize::default();
        let result = self.visit(region, |rec| {
            staged.entries.push(rec.entry());
            size.records += 1;
            if fields.is_empty() {
                return Ok(());
            }
            if want_mods {
                let (a, b) = rec.aux();
                decode_record(
                    &mut staged.decoder,
                    &[a, b],
                    ModSeq::Packed {
                        bytes: rec.packed_seq(),
                        len: rec.l_seq(),
                    },
                    rec.flag,
                    rec.pos as i64,
                    rec.cigar(),
                    rec.name(),
                    modification_codes,
                    &mut staged.mods,
                )?;
                staged.mod_ends.push(staged.mods.len() as i64);
            }
            if !want_record_fields {
                return Ok(());
            }
            let bytes = &mut staged.bytes;
            let mut keep = |part: &[u8]| {
                bytes.extend_from_slice(part);
                bytes.len()
            };
            let name = if fields.contains(AlignmentFields::NAME) {
                size.name_bytes += rec.name().len();
                keep(rec.name())
            } else {
                keep(&[])
            };
            let cigar = if fields.contains(AlignmentFields::CIGAR) {
                size.cigar_ops += rec.n_cigar();
                keep(rec.cigar_bytes())
            } else {
                keep(&[])
            };
            let seq = if fields.contains(AlignmentFields::SEQ) {
                keep(rec.packed_seq())
            } else {
                keep(&[])
            };
            let qual = if fields.contains(AlignmentFields::QUAL) {
                keep(rec.qual())
            } else {
                keep(&[])
            };
            if fields.contains(AlignmentFields::SEQ) || fields.contains(AlignmentFields::QUAL) {
                size.seq_bytes += rec.l_seq();
            }
            let aux = if fields.contains(AlignmentFields::AUX) {
                size.aux_bytes += rec.aux_len();
                let (a, b) = rec.aux();
                keep(a);
                keep(b)
            } else {
                keep(&[])
            };
            staged.records.push(Kept {
                name_end: name,
                cigar_end: cigar,
                seq_end: seq,
                l_seq: rec.l_seq(),
                qual_end: qual,
                aux_end: aux,
            });
            Ok(())
        });
        size.modifications = staged.mods.len();
        staged.size = size;
        staged.queried = result.is_ok();
        self.staged = staged;
        result.map(|()| size)
    }

    fn write(&self, mut out: AlignmentsOut<'_>) -> Result<()> {
        let s = &self.staged;
        if !s.queried {
            return Err(Error::NotQueried);
        }
        out.check(s.fields, s.size)?;
        out.entries.copy_from_slice(&s.entries);
        let bytes = &s.bytes;
        if let Some(c) = &mut out.names {
            fill_strings(c, bytes, &s.records, |start, k| start..k.name_end);
        }
        if let Some(c) = &mut out.quals {
            fill_strings(c, bytes, &s.records, |_, k| k.seq_end..k.qual_end);
        }
        if let Some(c) = &mut out.aux {
            fill_strings(c, bytes, &s.records, |_, k| k.qual_end..k.aux_end);
        }
        if let Some(c) = &mut out.seqs {
            c.offsets[0] = 0;
            let mut at = 0;
            for (i, k) in s.records.iter().enumerate() {
                unpack_seq(
                    &bytes[k.cigar_end..k.seq_end],
                    k.l_seq,
                    &mut c.data[at..at + k.l_seq],
                );
                at += k.l_seq;
                c.offsets[i + 1] = at as i64;
            }
        }
        if let Some(m) = &mut out.modifications {
            m.offsets[0] = 0;
            m.offsets[1..].copy_from_slice(&s.mod_ends);
            m.entries.copy_from_slice(&s.mods);
        }
        if let Some(c) = &mut out.cigars {
            c.offsets[0] = 0;
            let mut at = 0;
            for (i, k) in s.records.iter().enumerate() {
                let packed = &bytes[k.name_end..k.cigar_end];
                let count = packed.len() / 4;
                #[cfg(target_endian = "little")]
                if count > 4 {
                    // SAFETY: query only stages whole four-byte CIGAR words,
                    // and write checks the destination against their total count.
                    // Copying bytes requires no source alignment; every u32 bit
                    // pattern is valid and little-endian words need no conversion.
                    // Reader-owned bytes cannot alias caller-owned mutable output.
                    unsafe {
                        std::ptr::copy_nonoverlapping(
                            packed.as_ptr(),
                            c.ops.as_mut_ptr().add(at).cast::<u8>(),
                            packed.len(),
                        );
                    }
                } else {
                    for (dst, word) in c.ops[at..at + count].iter_mut().zip(packed.chunks_exact(4))
                    {
                        *dst = u32::from_le_bytes(word.try_into().unwrap());
                    }
                }
                #[cfg(target_endian = "big")]
                for (dst, word) in c.ops[at..at + count].iter_mut().zip(packed.chunks_exact(4)) {
                    *dst = u32::from_le_bytes(word.try_into().unwrap());
                }
                at += count;
                c.offsets[i + 1] = at as i64;
            }
        }
        Ok(())
    }
}
