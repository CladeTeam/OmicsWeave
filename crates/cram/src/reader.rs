//! [`CramReader`]: a CRAM file with its `.crai`, read by region.

use std::cell::Cell;
use std::sync::Arc;

use omicsweave_core::cache::{BlockCache, BlockKey, SourceToken};

use omicsweave_bam::mods::{ModDecoder, ModSeq, decode_record};
use omicsweave_core::{
    AlignmentEntry, AlignmentFields, AlignmentRead, AlignmentSize, AlignmentsOut, Error,
    ModificationEntry, ReadAt, Region, Result, SequenceRead,
};

use crate::cache::{self, ContainerUse, Meta, PreparedContainer, PreparedSlice, SliceUse};
use crate::comp::{CompressionHeader, ds};
use crate::header::CramHeader;
use crate::index::CramIndex;
use crate::record::{BamFields, Naming};
use crate::reference::{ByName, NoReference};
use crate::slice::{Decode, Decoded, References};

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
    // Modifications are counted along the bases, from the MM / ML tags.
    let mods = fields.contains(AlignmentFields::MODIFICATIONS);
    if mods {
        wanted |= ds::SEQ | ds::TL | ds::AUX;
    }
    let bases = fields.contains(AlignmentFields::SEQ)
        || mods
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
    decoder: ModDecoder,
    mods: Vec<ModificationEntry>,
    mod_ends: Vec<i64>,
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
        self.mods.clear();
        for v in [
            &mut self.name_ends,
            &mut self.seq_ends,
            &mut self.cigar_ends,
            &mut self.aux_ends,
            &mut self.mod_ends,
        ] {
            v.clear();
        }
    }
}

/// Reusable decoding buffers.
#[derive(Default)]
struct Work {
    head: Vec<u8>,
    body: Vec<u8>,
    decoded: Decoded,
    scratch: Vec<u8>,
    fields: BamFields,
}

/// Counters of the work queries did.
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub struct CramStats {
    pub read_calls: u64,
    pub containers: u64,
    pub slices_decoded: u64,
    pub records_decoded: u64,
    pub records_considered: u64,
    pub bytes_read: u64,
}

/// A CRAM file, its `.crai` index and, when its slices need one, the
/// reference they were compressed against. Opening reads the header; each
/// query reads the containers from the one the index gives for its start.
/// Decoded slices are reused only when a caller explicitly attaches a
/// BlockCache. Query scratch and staged output are independent per context.
pub struct CramReader<S, R = NoReference> {
    source: S,
    reference: Option<R>,
    header: Arc<CramHeader>,
    index: Arc<CramIndex>,
    options: CramOptions,
    staged: Staged,
    work: Work,
    stats: CramStats,
    cache: Option<Arc<BlockCache>>,
    namespace: SourceToken,
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
            cache: None,
            namespace: SourceToken::new(),
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

    /// Attaches caller-managed immutable containers and decoded slices.
    /// The CRAM, reference, index and options must remain unchanged throughout
    /// this source session and all contexts cloned from it.
    pub fn set_cache(&mut self, cache: Option<Arc<BlockCache>>) {
        self.cache = cache;
    }

    pub fn cache_namespace(&self) -> SourceToken {
        self.namespace
    }

    /// Describes physical container and slice dependencies for these fields,
    /// without decoding records or constructing alignment output. Reading
    /// compression and slice headers may require metadata decompression.
    pub fn describe_cache(
        &mut self,
        region: Region<'_>,
        fields: AlignmentFields,
        out: &mut Vec<BlockKey>,
    ) -> Result<()> {
        self.visit(region, fields, Some(out), false, |_, _| Ok(()))
    }

    /// Loads the dependencies of the query without constructing its output.
    /// Without an attached cache the work is performed but not retained.
    pub fn prefetch_cache(&mut self, region: Region<'_>, fields: AlignmentFields) -> Result<()> {
        self.visit(region, fields, None, true, |_, _| Ok(()))
    }

    /// A new source/reference session sharing only immutable header and index.
    /// These parts must match that metadata. Changed reference contents are
    /// isolated even when the same external cache remains attached.
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
            cache: self.cache.clone(),
            namespace: SourceToken::new(),
        }
    }

    /// Clones handles of this same immutable source and semantically identical
    /// reference. The factories must duplicate these handles, not replace their
    /// contents; use `with_parts` for a new source/reference session.
    pub fn try_clone_with<C, F>(&self, source: C, reference: F) -> Result<Self>
    where
        C: FnOnce(&S) -> Result<S>,
        F: FnOnce(&R) -> Result<R>,
    {
        let mut cloned = self.with_parts(
            source(&self.source)?,
            self.reference.as_ref().map(reference).transpose()?,
        );
        cloned.namespace = self.namespace;
        Ok(cloned)
    }

    /// An independent execution context of a Clone source/reference session.
    pub fn clone_context(&self) -> Self
    where
        S: Clone,
        R: Clone,
    {
        self.try_clone_with(|s| Ok(s.clone()), |r| Ok(r.clone()))
            .expect("infallible handle clones")
    }

    /// Calls `visit` with the decoded slice and record index of every
    /// record of `region`, in file order: the records htslib's
    /// `sam_itr_queryi` returns for a CRAM. `fields` chooses what is
    /// decoded besides the fixed fields.
    fn visit(
        &mut self,
        region: Region<'_>,
        fields: AlignmentFields,
        mut keys: Option<&mut Vec<BlockKey>>,
        prepare: bool,
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
        let Some(offset) = self.index.query(range.refid, range.start) else {
            return Ok(());
        };
        let mut offset =
            u64::try_from(offset).map_err(|_| Error::invalid("CRAM index", "negative offset"))?;
        let calls = Cell::new(0u64);
        let bytes = Cell::new(0u64);
        let CramReader {
            source,
            reference,
            header,
            options,
            work,
            stats,
            cache: attached,
            namespace,
            ..
        } = self;
        let source = Counted {
            source,
            calls: &calls,
            bytes: &bytes,
        };
        let result = (|| {
            let mut first = true;
            loop {
                let key = cache::meta_key(*namespace, offset);
                if let Some(keys) = &mut keys {
                    keys.push(key);
                }
                let meta = cache::load(attached, key, || {
                    Meta::read(&source, offset, &mut work.head)
                })?;
                let Some(h) = &meta.header else {
                    return Ok(());
                };
                let next = offset
                    .checked_add(h.header_len as u64)
                    .and_then(|n| n.checked_add(h.length as u64))
                    .ok_or_else(|| Error::invalid("CRAM container", "offset overflow"))?;
                if h.length == 0 {
                    offset = next;
                    continue;
                }
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
                    if h.ref_seq_id != range.refid
                        || (range.refid != -1 && h.ref_seq_start > range.end)
                    {
                        return Ok(());
                    }
                    if range.refid != -1 && h.ref_seq_end() < range.start {
                        offset = next;
                        continue;
                    }
                }
                let key = cache::body_key(*namespace, offset);
                if let Some(keys) = &mut keys {
                    keys.push(key);
                }
                let container = cache::load(attached, key, || {
                    PreparedContainer::read(
                        &source,
                        offset,
                        h,
                        &mut work.scratch,
                        options.verify_checksums,
                        stats,
                        std::mem::take(&mut work.body),
                    )
                })?;
                let container = ContainerUse {
                    loaded: container,
                    spare: &mut work.body,
                };
                let (wanted, bases) = series(fields, &container.comp, options.decode_md);
                for (index, part) in container.parts.iter().enumerate() {
                    let sh = &part.header;
                    if sh.ref_seq_id != -2 {
                        if sh.ref_seq_id != range.refid
                            || (range.refid != -1 && sh.ref_seq_start > range.end)
                        {
                            return Ok(());
                        }
                        if range.refid != -1 && sh.ref_seq_end() < range.start {
                            continue;
                        }
                    }
                    let key =
                        cache::slice_key(*namespace, offset, index, range.refid, fields, options);
                    if let Some(keys) = &mut keys {
                        keys.push(key);
                    }
                    if !prepare {
                        continue;
                    }
                    let decoded = cache::load(attached, key, || {
                        let input = container.slice(part)?;
                        let decode = Decode {
                            comp: &container.comp,
                            ref_lens: &header.ref_lens,
                            unknown_rg: header.unknown_rg,
                            wanted,
                            bases,
                            // MM/ML access needs stored AUX fields, but only
                            // requested AUX output needs absent MD/NM tags
                            // synthesized. Stored placeholders still rebuild.
                            decode_md: options.decode_md && fields.contains(AlignmentFields::AUX),
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
                        let mut decoded = std::mem::take(&mut work.decoded);
                        decode.run(&input, refs, &mut decoded)?;
                        stats.slices_decoded += 1;
                        stats.records_decoded += decoded.recs.len() as u64;
                        Ok(PreparedSlice::new(decoded, attached))
                    })?;
                    let decoded = SliceUse {
                        loaded: decoded,
                        spare: &mut work.decoded,
                    };
                    let candidates = decoded.candidates(range.refid, range.start, range.end);
                    for i in candidates.clone() {
                        stats.records_considered += 1;
                        let cr = &decoded.decoded.recs[i];
                        if cr.ref_id < range.refid && cr.ref_id != -1 {
                            continue;
                        }
                        if cr.ref_id != range.refid || (range.refid != -1 && cr.apos > range.end) {
                            return Ok(());
                        }
                        if range.refid != -1 && cr.aend < range.start {
                            continue;
                        }
                        visit(&decoded.decoded, i)?;
                    }
                    if candidates.end < decoded.decoded.recs.len() {
                        return Ok(());
                    }
                }
                offset = next;
            }
        })();
        stats.read_calls += calls.get();
        stats.bytes_read += bytes.get();
        result
    }
}

impl<S: ReadAt, R: SequenceRead> AlignmentRead for CramReader<S, R> {
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
        let want_mods = fields.contains(AlignmentFields::MODIFICATIONS);
        let result = self.visit(region, fields, None, true, |dec, i| {
            // Only CG conversion and requested AUX output need an owned copy;
            // modification tags can borrow the immutable decoded slice.
            dec.bam_fields(
                i,
                &naming,
                want_name,
                want_aux || dec.recs[i].has_cg,
                &mut f,
            )?;
            staged.entries.push(f.entry);
            if want_mods {
                let e = &f.entry;
                decode_record(
                    &mut staged.decoder,
                    &[dec.aux(i)],
                    ModSeq::Ascii(dec.seq(i, f.l_seq)),
                    e.flag,
                    e.start,
                    f.cigar.iter().map(|&c| ((c & 0xf) as u8, c >> 4)),
                    // Names are decoded only when asked for.
                    if f.name.is_empty() {
                        b"a record"
                    } else {
                        &f.name
                    },
                    modification_codes,
                    &mut staged.mods,
                )?;
                staged.mod_ends.push(staged.mods.len() as i64);
            }
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
            modifications: staged.mods.len(),
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
        if let Some(m) = &mut out.modifications {
            fill(m.offsets, &s.mod_ends);
            m.entries.copy_from_slice(&s.mods);
        }
        Ok(())
    }
}

struct Counted<'a, S> {
    source: &'a S,
    calls: &'a Cell<u64>,
    bytes: &'a Cell<u64>,
}
impl<S: ReadAt> ReadAt for Counted<'_, S> {
    fn len(&self) -> u64 {
        self.source.len()
    }
    fn read_exact_at(&self, out: &mut [u8], offset: u64) -> Result<()> {
        self.calls.set(self.calls.get().saturating_add(1));
        self.bytes
            .set(self.bytes.get().saturating_add(out.len() as u64));
        self.source.read_exact_at(out, offset)
    }
}
