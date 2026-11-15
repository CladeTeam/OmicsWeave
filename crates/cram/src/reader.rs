//! [`CramReader`]: a CRAM file with its `.crai`, read by region.

use std::sync::Arc;

use omicsweave_core::{
    AlignmentEntry, AlignmentFields, AlignmentRead, AlignmentSize, AlignmentsOut, Error, ReadAt,
    Region, Result, SequenceRead,
};

use crate::comp::{CompressionHeader, ds};
use crate::container::{Block, content};
use crate::header::{CramHeader, read_container_body, read_container_header};
use crate::index::CramIndex;
use crate::record::{BamFields, Naming};
use crate::reference::{ByName, NoReference};
use crate::slice::{Decode, Decoded, References, SliceInput};

/// How records are decoded.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct CramOptions {
    /// Check block CRC32s and each slice's reference MD5 against the
    /// reference given, as htslib does unless told to ignore them.
    pub verify_checksums: bool,
    /// Generate the MD and NM tags of records that store none, as htslib
    /// does by default.
    pub decode_md: bool,
    /// Prefix of the names given to records stored without one (htslib:
    /// the file's base name).
    pub name_prefix: Vec<u8>,
}

impl Default for CramOptions {
    fn default() -> CramOptions {
        CramOptions {
            verify_checksums: true,
            decode_md: true,
            name_prefix: Vec::new(),
        }
    }
}

/// The data series every query decodes: what the fixed fields of
/// [`AlignmentEntry`] and the CIGAR are computed from.
pub(crate) const FIXED: u32 = ds::BF
    | ds::RI
    | ds::AP
    | ds::MQ
    | ds::CIGAR
    | ds::CF
    | ds::NF
    | ds::NS
    | ds::NP
    | ds::TS
    | ds::MF;

/// The data series for `fields`, and whether bases are reconstructed.
pub(crate) fn series(
    fields: AlignmentFields,
    comp: &CompressionHeader,
    decode_md: bool,
) -> (u32, bool) {
    let mut wanted = FIXED;
    if comp.has_cg {
        wanted |= ds::TL | ds::AUX;
    }
    if fields.contains(AlignmentFields::NAME) {
        wanted |= ds::RN;
    }
    if fields.contains(AlignmentFields::SEQ) {
        wanted |= ds::SEQ;
    }
    if fields.contains(AlignmentFields::QUAL) {
        wanted |= ds::QUAL;
    }
    if fields.contains(AlignmentFields::AUX) {
        wanted |= ds::RG | ds::TL | ds::AUX;
        if decode_md {
            wanted |= ds::SEQ;
        }
    }
    let bases = fields.contains(AlignmentFields::SEQ)
        || (fields.contains(AlignmentFields::AUX) && decode_md);
    (wanted, bases)
}

/// What a query keeps for [`AlignmentRead::write`].
#[derive(Default)]
struct Staged {
    queried: bool,
    fields: AlignmentFields,
    size: AlignmentSize,
    entries: Vec<AlignmentEntry>,
    names: Vec<u8>,
    name_ends: Vec<i64>,
    seqs: Vec<u8>,
    quals: Vec<u8>,
    seq_ends: Vec<i64>,
    cigars: Vec<u32>,
    cigar_ends: Vec<i64>,
    aux: Vec<u8>,
    aux_ends: Vec<i64>,
}

impl Staged {
    fn clear(&mut self, fields: AlignmentFields) {
        self.queried = false;
        self.fields = fields;
        self.size = AlignmentSize::default();
        self.entries.clear();
        for v in [
            &mut self.names,
            &mut self.seqs,
            &mut self.quals,
            &mut self.aux,
        ] {
            v.clear();
        }
        self.cigars.clear();
        for v in [
            &mut self.name_ends,
            &mut self.seq_ends,
            &mut self.cigar_ends,
            &mut self.aux_ends,
        ] {
            v.clear();
        }
    }
}

/// Reusable decoding buffers.
#[derive(Default)]
struct Work {
    head: Vec<u8>,
    /// The last container read: its file offset, bytes after the header
    /// and compression header, and the slices decoded from it.
    body_at: Option<u64>,
    body: Vec<u8>,
    comp: Option<CompressionHeader>,
    slices: Vec<CachedSlice>,
    scratch: Vec<u8>,
    fields: BamFields,
}

/// A decoded slice of the last container, kept for following queries.
struct CachedSlice {
    index: usize,
    wanted: u32,
    bases: bool,
    /// The query reference a multi-reference slice was decoded for.
    range_ref: i32,
    decoded: Decoded,
}

/// Counters of the work queries did.
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub struct CramStats {
    pub containers: u64,
    pub slices_decoded: u64,
    pub records_decoded: u64,
    pub bytes_read: u64,
}

/// A CRAM file, its `.crai` index and, when its slices need one, the
/// reference they were compressed against. Opening reads the header; each
/// query reads the containers from the one the index gives for its start.
/// The last container read and the slices decoded from it are kept, so
/// queries close to each other decode each slice once.
pub struct CramReader<S, R = NoReference> {
    source: S,
    reference: Option<R>,
    header: Arc<CramHeader>,
    index: Arc<CramIndex>,
    options: CramOptions,
    staged: Staged,
    work: Work,
    stats: CramStats,
}

/// The range of a query, as htslib's `cram_range`: 1-based start and
/// inclusive end.
#[derive(Clone, Copy)]
struct Range {
    refid: i32,
    start: i64,
    end: i64,
}

impl<S: ReadAt, R: SequenceRead> CramReader<S, R> {
    /// Reads the header of `source`, a CRAM file that `index` indexes.
    /// `reference` supplies the bases of slices that refer to an external
    /// reference; without it only fields that need no bases can be read
    /// from them.
    pub fn open(
        source: S,
        index: CramIndex,
        reference: Option<R>,
        options: CramOptions,
    ) -> Result<CramReader<S, R>> {
        let header = CramHeader::read(&source)?;
        Ok(CramReader {
            source,
            reference,
            header: Arc::new(header),
            index: Arc::new(index),
            options,
            staged: Staged::default(),
            work: Work::default(),
            stats: CramStats::default(),
        })
    }

    pub fn header(&self) -> &CramHeader {
        &self.header
    }

    pub fn index(&self) -> &CramIndex {
        &self.index
    }

    pub fn source(&self) -> &S {
        &self.source
    }

    /// The reference given at open.
    pub fn reference(&self) -> Option<&R> {
        self.reference.as_ref()
    }

    pub fn options(&self) -> &CramOptions {
        &self.options
    }

    pub fn stats(&self) -> CramStats {
        self.stats
    }

    pub fn reset_stats(&mut self) {
        self.stats = CramStats::default();
    }

    /// A reader of the same file, index and options with another source
    /// handle and reference; staged results start empty.
    pub fn with_parts<S2: ReadAt, R2: SequenceRead>(
        &self,
        source: S2,
        reference: Option<R2>,
    ) -> CramReader<S2, R2> {
        CramReader {
            source,
            reference,
            header: Arc::clone(&self.header),
            index: Arc::clone(&self.index),
            options: self.options.clone(),
            staged: Staged::default(),
            work: Work::default(),
            stats: CramStats::default(),
        }
    }

    /// Calls `visit` with the decoded slice and record index of every
    /// record of `region`, in file order: the records htslib's
    /// `sam_itr_queryi` returns for a CRAM. `fields` chooses what is
    /// decoded besides the fixed fields.
    fn visit(
        &mut self,
        region: Region<'_>,
        fields: AlignmentFields,
        mut visit: impl FnMut(&Decoded, usize) -> Result<()>,
    ) -> Result<()> {
        region.check()?;
        let tid = self
            .header
            .tid(region.chrom)
            .ok_or_else(|| Error::UnknownChrom(region.chrom.to_owned()))?;
        if region.is_empty() {
            return Ok(());
        }
        let range = Range {
            refid: tid as i32,
            start: region.start as i64 + 1,
            end: region.end as i64,
        };
        let Some(mut offset) = self.index.query(range.refid, range.start) else {
            return Ok(());
        };
        let offset_u =
            |o: i64| u64::try_from(o).map_err(|_| Error::invalid("CRAM index", "negative offset"));
        let CramReader {
            source,
            reference,
            header,
            options,
            work,
            stats,
            ..
        } = self;
        let mut first = true;
        loop {
            let Some(c) = read_container_header(&*source, offset_u(offset)?, &mut work.head)?
            else {
                return Ok(());
            };
            let next = offset + c.header.header_len as i64 + c.header.length.max(0) as i64;
            if c.header.length == 0 {
                offset = next;
                continue;
            }
            let h = &c.header;
            if first {
                if h.ref_seq_id != -2
                    && (h.ref_seq_id < range.refid
                        || (range.refid >= 0
                            && h.ref_seq_id == range.refid
                            && h.ref_seq_end() < range.start))
                {
                    offset = next;
                    continue;
                }
                first = false;
                if h.ref_seq_id != -2 && h.ref_seq_id != range.refid {
                    return Ok(());
                }
            } else if h.ref_seq_id != -2 {
                if h.ref_seq_id != range.refid {
                    return Ok(());
                }
                if range.refid != -1 && h.ref_seq_start > range.end {
                    return Ok(());
                }
                if range.refid != -1 && h.ref_seq_end() < range.start {
                    offset = next;
                    continue;
                }
            }
            let Work {
                body_at,
                body,
                comp,
                slices,
                scratch,
                ..
            } = &mut *work;
            if *body_at != Some(c.offset) {
                *body_at = None;
                slices.clear();
                read_container_body(&*source, &c, body)?;
                stats.containers += 1;
                stats.bytes_read += (c.header.header_len + body.len()) as u64;
                let mut at = 0;
                let comp_block = Block::parse(body, &mut at)?;
                if comp_block.content_type != content::COMPRESSION_HEADER {
                    return Err(Error::invalid("CRAM container", "no compression header"));
                }
                *comp = Some(CompressionHeader::parse(
                    comp_block.uncompress(scratch, options.verify_checksums)?,
                )?);
                *body_at = Some(c.offset);
            }
            let comp = comp.as_ref().unwrap();
            if c.header.num_records == 0 {
                offset = next;
                continue;
            }
            let mut at = 0;
            Block::parse(body, &mut at)?;
            let (wanted, bases) = series(fields, comp, options.decode_md);
            for index in 0..c.header.landmarks.len() {
                let slice = SliceInput::read(body, &mut at, scratch)?;
                let sh = &slice.header;
                if sh.ref_seq_id != -2 {
                    if sh.ref_seq_id != range.refid {
                        return Ok(());
                    }
                    if range.refid != -1 && sh.ref_seq_start > range.end {
                        return Ok(());
                    }
                    if range.refid != -1 && sh.ref_seq_end() < range.start {
                        continue;
                    }
                }
                let range_ref = if sh.ref_seq_id == -2 { range.refid } else { -3 };
                let found = slices.iter().position(|k| {
                    k.index == index
                        && k.range_ref == range_ref
                        && k.wanted & wanted == wanted
                        && (k.bases || !bases)
                });
                let at_cache = match found {
                    Some(i) => i,
                    None => {
                        let decode = Decode {
                            comp,
                            ref_lens: &header.ref_lens,
                            unknown_rg: header.unknown_rg,
                            wanted,
                            bases,
                            decode_md: options.decode_md,
                            mates: true,
                            holds_refs: true,
                            range_ref: range.refid,
                            verify: options.verify_checksums,
                        };
                        let mut by_name;
                        let refs: Option<&mut dyn References> = match reference.as_mut() {
                            Some(r) => {
                                by_name = ByName {
                                    reader: r,
                                    refs: header.references(),
                                };
                                Some(&mut by_name)
                            }
                            None => None,
                        };
                        let mut decoded = Decoded::default();
                        decode.run(&slice, refs, &mut decoded)?;
                        stats.slices_decoded += 1;
                        stats.records_decoded += decoded.recs.len() as u64;
                        slices.retain(|k| k.index != index || k.range_ref != range_ref);
                        slices.push(CachedSlice {
                            index,
                            wanted,
                            bases,
                            range_ref,
                            decoded,
                        });
                        slices.len() - 1
                    }
                };
                let dec = &slices[at_cache].decoded;
                for (i, cr) in dec.recs.iter().enumerate() {
                    if cr.ref_id < range.refid && cr.ref_id != -1 {
                        continue;
                    }
                    if cr.ref_id != range.refid {
                        return Ok(());
                    }
                    if range.refid != -1 && cr.apos > range.end {
                        return Ok(());
                    }
                    if range.refid != -1 && cr.aend < range.start {
                        continue;
                    }
                    visit(dec, i)?;
                }
            }
            offset = next;
        }
    }
}

impl<S: ReadAt, R: SequenceRead> AlignmentRead for CramReader<S, R> {
    fn chrom_len(&self, chrom: &str) -> Result<u32> {
        self.header.chrom_len(chrom)
    }

    fn query(&mut self, region: Region<'_>, fields: AlignmentFields) -> Result<AlignmentSize> {
        let mut staged = std::mem::take(&mut self.staged);
        staged.clear(fields);
        let mut f = std::mem::take(&mut self.work.fields);
        let header = Arc::clone(&self.header);
        let prefix = self.options.name_prefix.clone();
        let naming = Naming {
            prefix: &prefix,
            read_groups: &header.read_groups,
        };
        let want_name = fields.contains(AlignmentFields::NAME);
        let want_aux = fields.contains(AlignmentFields::AUX);
        let want_seq = fields.contains(AlignmentFields::SEQ);
        let want_qual = fields.contains(AlignmentFields::QUAL);
        let want_cigar = fields.contains(AlignmentFields::CIGAR);
        let result = self.visit(region, fields, |dec, i| {
            dec.bam_fields(i, &naming, want_name, want_aux || dec.has_aux, &mut f)?;
            staged.entries.push(f.entry);
            if want_name {
                staged.names.extend_from_slice(&f.name);
                staged.name_ends.push(staged.names.len() as i64);
            }
            if want_cigar {
                staged.cigars.extend_from_slice(&f.cigar);
                staged.cigar_ends.push(staged.cigars.len() as i64);
            }
            if want_seq {
                dec.push_seq(i, f.l_seq, &mut staged.seqs);
            }
            if want_qual {
                dec.push_qual(i, f.l_seq, &mut staged.quals);
            }
            if want_seq || want_qual {
                let prev = staged.seq_ends.last().copied().unwrap_or(0);
                staged.seq_ends.push(prev + f.l_seq as i64);
            }
            if want_aux {
                staged.aux.extend_from_slice(&f.aux);
                staged.aux_ends.push(staged.aux.len() as i64);
            }
            Ok(())
        });
        self.work.fields = f;
        staged.size = AlignmentSize {
            records: staged.entries.len(),
            name_bytes: staged.names.len(),
            seq_bytes: staged.seq_ends.last().copied().unwrap_or(0) as usize,
            cigar_ops: staged.cigars.len(),
            aux_bytes: staged.aux.len(),
        };
        staged.queried = result.is_ok();
        let size = staged.size;
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
        let fill = |offsets: &mut [i64], ends: &[i64]| {
            offsets[0] = 0;
            offsets[1..].copy_from_slice(ends);
        };
        if let Some(c) = &mut out.names {
            fill(c.offsets, &s.name_ends);
            c.data.copy_from_slice(&s.names);
        }
        if let Some(c) = &mut out.seqs {
            fill(c.offsets, &s.seq_ends);
            c.data.copy_from_slice(&s.seqs);
        }
        if let Some(c) = &mut out.quals {
            fill(c.offsets, &s.seq_ends);
            c.data.copy_from_slice(&s.quals);
        }
        if let Some(c) = &mut out.cigars {
            fill(c.offsets, &s.cigar_ends);
            c.ops.copy_from_slice(&s.cigars);
        }
        if let Some(c) = &mut out.aux {
            fill(c.offsets, &s.aux_ends);
            c.data.copy_from_slice(&s.aux);
        }
        Ok(())
    }
}
