//! The per-record work shared by every accumulator, done once and only when
//! asked for: FLAG / MAPQ / NH filters, optional-field lookup, library
//! strandedness, aligned blocks, junctions, fragment spans, Tn5 cut sites
//! and the read's bases against the reference ([`ReadView`]), and mate
//! pairing ([`Pairs`], [`Mates`]).

use std::cell::{OnceCell, RefCell};
use std::cmp::Reverse;
use std::collections::{BinaryHeap, HashMap};

use omicsweave_bam::{Record, flag, op};
use omicsweave_core::Result;
use smallvec::SmallVec;

use crate::assign::{Query, SpaceId};

/// Aligned reference intervals of a read, ascending.
pub type Blocks = SmallVec<[(u32, u32); 4]>;

/// Columns (feature ids) a read or fragment is assigned to, ascending and
/// distinct.
pub type Hits = SmallVec<[u32; 4]>;

/// Aligned segments of a read: `(reference start, query start, length)` of
/// each `M` / `=` / `X` operation, in alignment order.
pub type Segments = SmallVec<[(u32, u32, u32); 4]>;

/// Reference bases fetched beyond each end of a read, so that the sequence
/// context of its first and last bases is known.
pub const REFERENCE_PAD: u32 = 2;

/// Reference sequence only around aligned segments, with nearby flanks merged.
struct ReferenceContext {
    ranges: SmallVec<[(u32, u32, usize); 4]>,
    bases: Vec<u8>,
}

/// Which records an accumulator counts.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Filter {
    /// FLAG bits a record must all have.
    pub require: u16,
    /// FLAG bits a record must have none of.
    pub exclude: u16,
    pub min_mapq: u8,
    /// Only records whose `NH` tag is 1 (a record without `NH` counts as
    /// unique).
    pub unique: bool,
}

impl Filter {
    /// Unmapped, secondary, QC-failed and duplicate records excluded
    /// (`samtools depth`).
    pub const DEFAULT_EXCLUDE: u16 =
        flag::UNMAPPED | flag::SECONDARY | flag::QC_FAIL | flag::DUPLICATE;

    /// [`Filter::DEFAULT_EXCLUDE`], any MAPQ, multi-mappers included.
    pub fn new() -> Filter {
        Filter {
            require: 0,
            exclude: Filter::DEFAULT_EXCLUDE,
            min_mapq: 0,
            unique: false,
        }
    }

    /// Why `r` is not counted, or `None` when it is.
    #[inline]
    pub fn reject(&self, r: &ReadView<'_>) -> Option<Rejected> {
        let f = r.flag();
        if f & self.exclude != 0 || f & self.require != self.require {
            return Some(Rejected::Flag);
        }
        if r.mapq() < self.min_mapq {
            return Some(Rejected::Mapq);
        }
        if self.unique && r.nh() > 1 {
            return Some(Rejected::MultiMapping);
        }
        None
    }
}

impl Default for Filter {
    fn default() -> Filter {
        Filter::new()
    }
}

/// Why a [`Filter`] turned a record away.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Rejected {
    Flag,
    Mapq,
    MultiMapping,
}

/// How reads relate to the transcribed strand.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum Strandedness {
    /// No relation (featureCounts `-s 0`).
    #[default]
    Unstranded,
    /// Read 1 (or a single-end read) lies on the transcript's strand
    /// (`-s 1`, ligation / SMART-seq3 style).
    Forward,
    /// Read 1 lies on the opposite strand (`-s 2`, dUTP / TruSeq Stranded,
    /// 10x 3' read 2 on the transcript strand).
    Reverse,
}

impl Strandedness {
    /// `0`, `1` or `2`, as featureCounts numbers them.
    pub fn from_code(code: u8) -> Option<Strandedness> {
        Some(match code {
            0 => Strandedness::Unstranded,
            1 => Strandedness::Forward,
            2 => Strandedness::Reverse,
            _ => return None,
        })
    }
}

/// An optional field's value.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum Aux<'a> {
    Char(u8),
    Int(i64),
    Float(f64),
    /// `Z` (without its NUL) or `H`.
    Str(&'a [u8]),
    /// `B`: element type and the raw little-endian elements.
    Array(u8, &'a [u8]),
}

impl<'a> Aux<'a> {
    /// The value as text: a string, or a one-byte character.
    pub fn bytes(&self) -> Option<&'a [u8]> {
        match *self {
            Aux::Str(s) => Some(s),
            _ => None,
        }
    }

    pub fn int(&self) -> Option<i64> {
        match *self {
            Aux::Int(v) => Some(v),
            _ => None,
        }
    }
}

fn aux_width(t: u8) -> Option<usize> {
    Some(match t {
        b'A' | b'c' | b'C' => 1,
        b's' | b'S' => 2,
        b'i' | b'I' | b'f' => 4,
        b'd' => 8,
        _ => return None,
    })
}

/// The value of `tag` in `aux` (optional fields as stored); a malformed
/// field ends the search.
pub fn find_aux(aux: &[u8], tag: [u8; 2]) -> Option<Aux<'_>> {
    let mut at = 0usize;
    while at + 3 <= aux.len() {
        let (key, t) = ([aux[at], aux[at + 1]], aux[at + 2]);
        let v = at + 3;
        let (value, next) = match t {
            b'Z' | b'H' => {
                let len = aux[v..].iter().position(|&c| c == 0)?;
                (Aux::Str(&aux[v..v + len]), v + len + 1)
            }
            b'B' => {
                let sub = *aux.get(v)?;
                let n = u32::from_le_bytes(aux.get(v + 1..v + 5)?.try_into().ok()?) as usize;
                let w = aux_width(sub)?;
                let end = v.checked_add(5)?.checked_add(n.checked_mul(w)?)?;
                (Aux::Array(sub, aux.get(v + 5..end)?), end)
            }
            _ => {
                let w = aux_width(t)?;
                let b = aux.get(v..v + w)?;
                let value = match t {
                    b'A' => Aux::Char(b[0]),
                    b'c' => Aux::Int(b[0] as i8 as i64),
                    b'C' => Aux::Int(b[0] as i64),
                    b's' => Aux::Int(i16::from_le_bytes([b[0], b[1]]) as i64),
                    b'S' => Aux::Int(u16::from_le_bytes([b[0], b[1]]) as i64),
                    b'i' => Aux::Int(i32::from_le_bytes(b.try_into().ok()?) as i64),
                    b'I' => Aux::Int(u32::from_le_bytes(b.try_into().ok()?) as i64),
                    b'f' => Aux::Float(f32::from_le_bytes(b.try_into().ok()?) as f64),
                    _ => Aux::Float(f64::from_le_bytes(b.try_into().ok()?)),
                };
                (value, v + w)
            }
        };
        if key == tag {
            return Some(value);
        }
        at = next;
    }
    None
}

/// Positions of a Tn5 insertion relative to a read's 5' end: +4 on the
/// forward strand, -5 on the reverse strand (10x Genomics / ArchR).
pub const TN5_FORWARD_SHIFT: u32 = 4;
/// See [`TN5_FORWARD_SHIFT`].
pub const TN5_REVERSE_SHIFT: u32 = 5;

/// A fragment of a read pair, 0-based half-open on the reference.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, PartialOrd, Ord)]
pub struct Fragment {
    pub start: u32,
    pub end: u32,
}

/// What a batch offers the records it holds: feature assignment through its
/// sweep cursors, and the reference bases and known variants of the region
/// its records cover.
pub(crate) trait Context {
    fn assign(&self, id: SpaceId, r: &ReadView<'_>, query: Query, strand: Option<bool>) -> Hits;
    /// The columns of space `id` that `[start, end)` overlaps, for an
    /// interval within the read's span (at or after its position).
    fn assign_at(&self, id: SpaceId, r: &ReadView<'_>, start: u32, end: u32) -> Hits;
    /// Calls `with` with the reference bases of `[start, end)` on `tid`
    /// (uppercase; `N` past the sequence's end); `false` when the task has
    /// no reference for `tid`.
    fn reference(
        &self,
        tid: u32,
        start: u32,
        end: u32,
        with: &mut dyn FnMut(&[u8]),
    ) -> Result<bool>;
    /// Appends the positions in `[start, end)` of `tid` that a known
    /// variant covers, ascending.
    fn variants(&self, tid: u32, start: u32, end: u32, out: &mut Vec<u32>) -> Result<()>;
}

/// A read base over a reference base: one position of an `M`, `=` or `X`
/// operation.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct AlignedBase {
    pub ref_pos: u32,
    /// Position in the read as stored (soft clips counted, hard clips not).
    pub qpos: u32,
    /// The read's base, uppercase IUPAC.
    pub read: u8,
    /// The reference base, uppercase IUPAC (`N` when unknown).
    pub reference: u8,
    /// Phred quality; 255 when the record has none.
    pub qual: u8,
    /// A known variant covers `ref_pos`.
    pub masked: bool,
}

/// `=ACMGRSVTWYHKDBN`: the 4-bit base codes of BAM.
const NT16: &[u8; 16] = b"=ACMGRSVTWYHKDBN";

/// One placed record as accumulators see it. Derived values (blocks, tags,
/// NH, assignments, aligned bases) are computed on first use and shared by
/// every accumulator that asks.
pub struct ReadView<'a> {
    /// The input BAM (0-based, as given).
    pub input: u32,
    /// The record's number in its BAM, from 0.
    pub ordinal: u64,
    pub rec: Record<'a>,
    ctx: &'a (dyn Context + 'a),
    blocks: OnceCell<Blocks>,
    blocks_del: OnceCell<Blocks>,
    nh: OnceCell<u32>,
    segments: OnceCell<Segments>,
    reference: OnceCell<Option<(u32, Vec<u8>)>>,
    reference_context: OnceCell<Option<ReferenceContext>>,
    mask: OnceCell<Vec<u32>>,
    aligned: OnceCell<Option<Vec<AlignedBase>>>,
    /// (space, query, strand) → hits, for this record.
    hits: RefCell<SmallVec<[(u64, Hits); 2]>>,
}

impl<'a> ReadView<'a> {
    pub(crate) fn new(
        input: u32,
        ordinal: u64,
        rec: Record<'a>,
        ctx: &'a (dyn Context + 'a),
    ) -> ReadView<'a> {
        ReadView {
            input,
            ordinal,
            rec,
            ctx,
            blocks: OnceCell::new(),
            blocks_del: OnceCell::new(),
            nh: OnceCell::new(),
            segments: OnceCell::new(),
            reference: OnceCell::new(),
            reference_context: OnceCell::new(),
            mask: OnceCell::new(),
            aligned: OnceCell::new(),
            hits: RefCell::new(SmallVec::new()),
        }
    }

    /// Reference id (the record is placed).
    #[inline]
    pub fn tid(&self) -> u32 {
        self.rec.tid as u32
    }

    /// 0-based leftmost position.
    #[inline]
    pub fn pos(&self) -> u32 {
        self.rec.pos.max(0) as u32
    }

    /// End of the alignment on the reference (`bam_endpos`).
    #[inline]
    pub fn end(&self) -> u32 {
        self.rec.end().clamp(0, u32::MAX as i64) as u32
    }

    #[inline]
    pub fn flag(&self) -> u16 {
        self.rec.flag
    }

    #[inline]
    pub fn mapq(&self) -> u8 {
        self.rec.mapq
    }

    #[inline]
    pub fn is_reverse(&self) -> bool {
        self.rec.flag & flag::REVERSE != 0
    }

    #[inline]
    pub fn is_paired(&self) -> bool {
        self.rec.flag & flag::PAIRED != 0
    }

    /// Read 2 of a pair.
    #[inline]
    pub fn is_read2(&self) -> bool {
        self.rec.flag & (flag::PAIRED | flag::READ2) == flag::PAIRED | flag::READ2
    }

    pub fn name(&self) -> &'a [u8] {
        self.rec.name()
    }

    /// The optional field `tag` (a `CG` tag holding the CIGAR is not one).
    pub fn tag(&self, tag: [u8; 2]) -> Option<Aux<'a>> {
        let (a, b) = self.rec.aux();
        find_aux(a, tag).or_else(|| find_aux(b, tag))
    }

    /// A string field (`Z`, `H`) or a character (`A`) as bytes.
    pub fn tag_bytes(&self, tag: [u8; 2]) -> Option<&'a [u8]> {
        let (a, b) = self.rec.aux();
        let raw = |aux: &'a [u8]| -> Option<&'a [u8]> {
            match find_aux(aux, tag)? {
                Aux::Str(s) => Some(s),
                Aux::Char(_) => {
                    // Point at the character inside the record's bytes.
                    let at = locate(aux, tag)?;
                    aux.get(at + 3..at + 4)
                }
                _ => None,
            }
        };
        raw(a).or_else(|| raw(b))
    }

    /// The `NH` tag: alignments reported for the read; 1 without the tag.
    pub fn nh(&self) -> u32 {
        *self.nh.get_or_init(|| {
            self.tag(*b"NH")
                .and_then(|a| a.int())
                .map_or(1, |v| v.clamp(0, u32::MAX as i64) as u32)
        })
    }

    /// The reference intervals of the read's `M` / `=` / `X` operations,
    /// adjacent ones merged; with `deletions`, `D` operations join the
    /// blocks they separate. `N` (and, without `deletions`, `D`) separates
    /// blocks; insertions and clips take no reference bases.
    pub fn blocks(&self, deletions: bool) -> &[(u32, u32)] {
        let cell = if deletions {
            &self.blocks_del
        } else {
            &self.blocks
        };
        cell.get_or_init(|| {
            let mut out = Blocks::new();
            let mut at = self.pos();
            for (o, len) in self.rec.cigar() {
                match o {
                    op::MATCH | op::EQUAL | op::DIFF => {
                        push_block(&mut out, at, at.saturating_add(len));
                        at = at.saturating_add(len);
                    }
                    op::DEL if deletions => {
                        push_block(&mut out, at, at.saturating_add(len));
                        at = at.saturating_add(len);
                    }
                    op::DEL | op::REF_SKIP => at = at.saturating_add(len),
                    _ => {}
                }
            }
            out
        })
    }

    /// The introns of a spliced alignment: the reference intervals of its
    /// `N` operations, ascending.
    pub fn junctions(&self) -> Blocks {
        let mut out = Blocks::new();
        let mut at = self.pos();
        for (o, len) in self.rec.cigar() {
            if o == op::REF_SKIP {
                out.push((at, at.saturating_add(len)));
            }
            if omicsweave_bam::consumes_reference(o) {
                at = at.saturating_add(len);
            }
        }
        out
    }

    /// The reference position of the read's 5' end.
    pub fn five_prime(&self) -> u32 {
        if self.is_reverse() {
            self.end().saturating_sub(1)
        } else {
            self.pos()
        }
    }

    /// The reference position of the read's 3' end.
    pub fn three_prime(&self) -> u32 {
        if self.is_reverse() {
            self.pos()
        } else {
            self.end().saturating_sub(1)
        }
    }

    /// Whether this record stands for its pair's fragment: the leftmost
    /// mate of a pair with both mates mapped to the same reference (at the
    /// same position, read 1).
    pub fn is_fragment_leader(&self) -> bool {
        let f = self.rec.flag;
        if f & flag::PAIRED == 0 || f & (flag::UNMAPPED | flag::MATE_UNMAPPED) != 0 {
            return false;
        }
        if self.rec.mate_tid != self.rec.tid || self.rec.tlen == 0 {
            return false;
        }
        self.rec.pos < self.rec.mate_pos
            || (self.rec.pos == self.rec.mate_pos && f & flag::READ1 != 0)
    }

    /// The fragment of a pair, from its leftmost mate: `[pos, pos + |TLEN|)`;
    /// `None` for any other record and for fragments longer than
    /// `max_len` (0: no limit).
    pub fn fragment(&self, max_len: u32) -> Option<Fragment> {
        if !self.is_fragment_leader() {
            return None;
        }
        let len = self.rec.tlen.unsigned_abs();
        if max_len > 0 && len > max_len {
            return None;
        }
        let start = self.pos();
        Some(Fragment {
            start,
            end: start.saturating_add(len),
        })
    }

    /// The Tn5-corrected fragment of a pair (10x `fragments.tsv`): the
    /// leftmost mate's 5' end + 4 to the rightmost mate's 5' end - 5, from
    /// the leftmost mate. Its Tn5 insertion sites are `start` and `end - 1`.
    /// `None` when [`ReadView::fragment`] is, or when the correction leaves
    /// nothing.
    pub fn atac_fragment(&self, max_len: u32) -> Option<Fragment> {
        let f = self.fragment(0)?;
        let start = f.start.saturating_add(TN5_FORWARD_SHIFT);
        let end = f.end.checked_sub(TN5_REVERSE_SHIFT)?;
        if end <= start || (max_len > 0 && end - start > max_len) {
            return None;
        }
        Some(Fragment { start, end })
    }

    /// Whether the read (as its fragment) lies on the reverse strand of the
    /// transcript: read 2 of a pair counts as read 1's strand, and a
    /// reverse-stranded library flips the result. Unstranded libraries are
    /// read as forward-stranded.
    pub fn reverse_for(&self, lib: Strandedness) -> bool {
        let fragment_reverse = self.is_reverse() ^ self.is_read2();
        fragment_reverse ^ (lib == Strandedness::Reverse)
    }

    /// The columns of feature space `space` that `query` of this read
    /// overlaps (by at least one base), keeping only features on `strand`
    /// when given (`Some(true)`: reverse); features without a strand match
    /// either. Computed once per record and query.
    pub fn assign(&self, space: SpaceId, query: Query, strand: Option<bool>) -> Hits {
        let key = (space.0 as u64) << 16
            | (query.code() as u64) << 2
            | strand.map_or(0, |s| 1 + s as u64);
        if let Some((_, h)) = self.hits.borrow().iter().find(|(k, _)| *k == key) {
            return h.clone();
        }
        let hits = self.ctx.assign(space, self, query, strand);
        self.hits.borrow_mut().push((key, hits.clone()));
        hits
    }

    /// The columns of `space` that `[start, end)` overlaps, for an interval
    /// within the read's span; not cached (positions of single bases).
    pub fn assign_at(&self, space: SpaceId, start: u32, end: u32) -> Hits {
        self.ctx.assign_at(space, self, start, end)
    }

    /// The read's `M` / `=` / `X` operations as `(reference start, query
    /// start, length)`, computed once.
    pub fn segments(&self) -> &[(u32, u32, u32)] {
        self.segments.get_or_init(|| {
            let mut out = Segments::new();
            let mut at = self.pos();
            let mut q = 0u32;
            for (o, len) in self.rec.cigar() {
                match o {
                    op::MATCH | op::EQUAL | op::DIFF => {
                        if len > 0 {
                            out.push((at, q, len));
                        }
                        at = at.saturating_add(len);
                        q = q.saturating_add(len);
                    }
                    op::INS | op::SOFT_CLIP => q = q.saturating_add(len),
                    op::DEL | op::REF_SKIP => at = at.saturating_add(len),
                    _ => {}
                }
            }
            out
        })
    }

    /// The query position aligned to reference position `p`, if a `M` /
    /// `=` / `X` operation covers it.
    pub fn query_at(&self, p: u32) -> Option<u32> {
        let s = self.segments();
        let i = s.partition_point(|x| x.0 <= p).checked_sub(1)?;
        let (r, q, len) = s[i];
        (p - r < len).then(|| q.checked_add(p - r)).flatten()
    }

    /// The reference position query base `q` is aligned to, if it lies in
    /// a `M` / `=` / `X` operation.
    pub fn reference_at(&self, q: u32) -> Option<u32> {
        let s = self.segments();
        let i = s.partition_point(|x| x.1 <= q).checked_sub(1)?;
        let (r, qs, len) = s[i];
        (q - qs < len).then(|| r.checked_add(q - qs)).flatten()
    }

    /// The task's reference bases under the read, [`REFERENCE_PAD`] bases
    /// beyond each end (`N` past the sequence), uppercase, fetched once:
    /// `(first position, bases)`; `None` without a reference for the read's
    /// sequence.
    pub fn reference_bases(&self) -> Result<Option<(u32, &[u8])>> {
        if let Some(v) = self.reference.get() {
            return Ok(v.as_ref().map(|(s, b)| (*s, b.as_slice())));
        }
        let lo = self.pos().saturating_sub(REFERENCE_PAD);
        let hi = self.end().max(self.pos() + 1).saturating_add(REFERENCE_PAD);
        let mut bases = Vec::new();
        let have = self
            .ctx
            .reference(self.tid(), lo, hi, &mut |b| bases.extend_from_slice(b))?;
        let v = self.reference.get_or_init(|| have.then_some((lo, bases)));
        Ok(v.as_ref().map(|(s, b)| (*s, b.as_slice())))
    }

    /// The reference base at `p` from [`ReadView::reference_bases`] (`N`
    /// outside them), or `None` without a reference.
    pub fn reference_base(&self, p: u32) -> Result<Option<u8>> {
        Ok(self.reference_bases()?.map(|(s, b)| {
            p.checked_sub(s)
                .and_then(|i| b.get(i as usize).copied())
                .unwrap_or(b'N')
        }))
    }

    fn reference_context(&self) -> Result<Option<&ReferenceContext>> {
        if let Some(v) = self.reference_context.get() {
            return Ok(v.as_ref());
        }
        let mut ranges: SmallVec<[(u32, u32, usize); 4]> = SmallVec::new();
        for &(start, _, len) in self.segments() {
            let lo = start.saturating_sub(REFERENCE_PAD);
            let hi = start.saturating_add(len).saturating_add(REFERENCE_PAD);
            if let Some(last) = ranges.last_mut()
                && lo <= last.1
            {
                last.1 = last.1.max(hi);
            } else {
                ranges.push((lo, hi, 0));
            }
        }
        let mut bases = Vec::new();
        let mut have = true;
        if ranges.is_empty() {
            have = self
                .ctx
                .reference(self.tid(), self.pos(), self.pos(), &mut |_| {})?;
        }
        for (lo, hi, offset) in &mut ranges {
            *offset = bases.len();
            if !self
                .ctx
                .reference(self.tid(), *lo, *hi, &mut |b| bases.extend_from_slice(b))?
            {
                have = false;
                break;
            }
        }
        let value = have.then_some(ReferenceContext { ranges, bases });
        Ok(self.reference_context.get_or_init(|| value).as_ref())
    }

    /// Whether the reference can supply aligned bases and their short context.
    /// This does not fetch introns or long deletions between aligned segments.
    pub(crate) fn has_reference_context(&self) -> Result<bool> {
        match self.reference.get() {
            Some(v) => Ok(v.is_some()),
            None => Ok(self.reference_context()?.is_some()),
        }
    }

    /// A base from an explicit full reference or aligned-segment context.
    /// Positions outside the available ranges read N.
    pub(crate) fn reference_context_base(&self, p: u32) -> Result<Option<u8>> {
        if let Some(v) = self.reference.get() {
            return Ok(v.as_ref().map(|(s, b)| {
                p.checked_sub(*s)
                    .and_then(|i| b.get(i as usize).copied())
                    .unwrap_or(b'N')
            }));
        }
        Ok(self.reference_context()?.map(|context| {
            let i = context.ranges.partition_point(|r| r.0 <= p);
            i.checked_sub(1)
                .and_then(|i| {
                    let &(lo, hi, offset) = &context.ranges[i];
                    (p < hi)
                        .then(|| offset.checked_add((p - lo) as usize))
                        .flatten()
                })
                .and_then(|i| context.bases.get(i).copied())
                .unwrap_or(b'N')
        }))
    }

    /// Positions in the read's span that a known variant covers, ascending,
    /// looked up once.
    pub fn masked(&self) -> Result<&[u32]> {
        if let Some(v) = self.mask.get() {
            return Ok(v);
        }
        let mut mask = Vec::new();
        self.ctx
            .variants(self.tid(), self.pos(), self.end(), &mut mask)?;
        Ok(self.mask.get_or_init(|| mask))
    }

    /// The read's base at query position `q`, uppercase IUPAC (`=` as
    /// stored).
    #[inline]
    pub fn base(&self, q: usize) -> u8 {
        let b = self.rec.packed_seq()[q / 2];
        NT16[((b >> (4 * (1 - q % 2))) & 0xf) as usize]
    }

    /// Hard-clipped bases before and after the stored read.
    pub fn hard_clips(&self) -> (u32, u32) {
        let mut ops = self.rec.cigar();
        let n = ops.len();
        let mut lead = 0;
        let mut tail = 0;
        if let Some((o, l)) = ops.next()
            && o == op::HARD_CLIP
        {
            lead = l;
        }
        if n > 1
            && let Some((o, l)) = self.rec.cigar().nth(n - 1)
            && o == op::HARD_CLIP
        {
            tail = l;
        }
        (lead, tail)
    }

    /// Every aligned base of the read with the reference base under it,
    /// computed once: the reference comes from the task's reference genome
    /// when it has the read's sequence, else from the `MD` tag. `None`
    /// when neither is there (or the `MD` tag does not fit the CIGAR), or
    /// the record has no bases.
    pub fn aligned(&self) -> Result<Option<&[AlignedBase]>> {
        if let Some(v) = self.aligned.get() {
            return Ok(v.as_deref());
        }
        let walked = self.walk()?;
        Ok(self.aligned.get_or_init(|| walked).as_deref())
    }

    fn walk(&self) -> Result<Option<Vec<AlignedBase>>> {
        let l_seq = self.rec.l_seq();
        if l_seq == 0 || self.flag() & flag::UNMAPPED != 0 {
            return Ok(None);
        }
        let lo = self.pos();
        let cached = self.reference.get().and_then(|v| v.as_ref());
        let have_reference = match self.reference.get() {
            Some(v) => v.is_some(),
            None => self.ctx.reference(self.tid(), lo, lo, &mut |_| {})?,
        };
        let mut md = if have_reference {
            None
        } else {
            match self.tag_bytes(*b"MD") {
                Some(t) => Some(Md::new(t)),
                None => return Ok(None),
            }
        };
        let mut reference = Vec::new();
        let mask = self.masked()?;
        let packed = self.rec.packed_seq();
        let qual = self.rec.qual();
        let no_qual = qual.first() == Some(&0xff);
        let mut out = Vec::with_capacity(l_seq);
        let mut at = lo;
        let mut q = 0usize;
        let mut m = 0usize;
        for (o, len) in self.rec.cigar() {
            match o {
                op::MATCH | op::EQUAL | op::DIFF => {
                    let Some(qend) = q.checked_add(len as usize) else {
                        return Ok(None);
                    };
                    if qend > l_seq || qend > qual.len() || qend.div_ceil(2) > packed.len() {
                        return Ok(None);
                    }
                    let Some(rend) = at.checked_add(len) else {
                        return Ok(None);
                    };
                    // Fetch aligned segments separately: an intron can
                    // be much larger than all the read's aligned bases.
                    // A caller's explicit full-span fetch is reused.
                    if cached.is_none() && have_reference {
                        reference.clear();
                        self.ctx.reference(self.tid(), at, rend, &mut |b| {
                            reference.extend_from_slice(b)
                        })?;
                    }
                    for (qi, pos) in (q..qend).zip(at..rend) {
                        let (byte, quality) = (packed[qi / 2], qual[qi]);
                        let read = NT16[((byte >> (4 * (1 - qi % 2))) & 0xf) as usize];
                        let refb = match &mut md {
                            None => match cached {
                                Some((s, b)) => pos
                                    .checked_sub(*s)
                                    .and_then(|i| b.get(i as usize).copied())
                                    .unwrap_or(b'N'),
                                None => reference.get((pos - at) as usize).copied().unwrap_or(b'N'),
                            },
                            Some(md) => match md.next_base(read) {
                                Some(b) => b,
                                None => return Ok(None),
                            },
                        };
                        let read = if read == b'=' { refb } else { read };
                        while m < mask.len() && mask[m] < pos {
                            m += 1;
                        }
                        out.push(AlignedBase {
                            ref_pos: pos,
                            qpos: qi as u32,
                            read,
                            reference: refb,
                            qual: if no_qual { 255 } else { quality },
                            masked: m < mask.len() && mask[m] == pos,
                        });
                    }
                    q = qend;
                    at = rend;
                }
                op::INS | op::SOFT_CLIP => q += len as usize,
                op::DEL => {
                    if let Some(md) = &mut md
                        && !md.deletion(len)
                    {
                        return Ok(None);
                    }
                    at = at.saturating_add(len);
                }
                op::REF_SKIP => at = at.saturating_add(len),
                _ => {}
            }
        }
        Ok(Some(out))
    }
}

/// A walk along an `MD` tag: runs of matching bases, mismatched reference
/// bases, and `^`-prefixed deleted bases, in alignment order.
struct Md<'a> {
    tag: &'a [u8],
    at: usize,
    /// Matching bases left in the current run.
    run: u32,
}

impl<'a> Md<'a> {
    fn new(tag: &'a [u8]) -> Md<'a> {
        Md { tag, at: 0, run: 0 }
    }

    fn number(&mut self) {
        while self.at < self.tag.len() && self.tag[self.at].is_ascii_digit() {
            self.run = self
                .run
                .saturating_mul(10)
                .saturating_add((self.tag[self.at] - b'0') as u32);
            self.at += 1;
        }
    }

    /// The reference base of the next aligned position, whose read base
    /// is `read`; `None` when the tag does not fit.
    fn next_base(&mut self, read: u8) -> Option<u8> {
        if self.run == 0 {
            self.number();
        }
        if self.run > 0 {
            self.run -= 1;
            return Some(if read == b'=' { b'N' } else { read });
        }
        let b = *self.tag.get(self.at)?;
        if !b.is_ascii_alphabetic() {
            return None;
        }
        self.at += 1;
        Some(b.to_ascii_uppercase())
    }

    /// Steps over a deletion of `len` bases.
    fn deletion(&mut self, len: u32) -> bool {
        if self.run == 0 {
            self.number();
        }
        if self.run > 0 || self.tag.get(self.at) != Some(&b'^') {
            return false;
        }
        self.at += 1;
        for _ in 0..len {
            match self.tag.get(self.at) {
                Some(b) if b.is_ascii_alphabetic() => self.at += 1,
                _ => return false,
            }
        }
        true
    }
}

/// The complement of an uppercase base (IUPAC codes other than `ACGT`
/// unchanged).
#[inline]
pub fn complement(b: u8) -> u8 {
    match b {
        b'A' => b'T',
        b'C' => b'G',
        b'G' => b'C',
        b'T' => b'A',
        other => other,
    }
}

/// Where `tag`'s field starts in `aux`.
fn locate(aux: &[u8], tag: [u8; 2]) -> Option<usize> {
    let mut at = 0usize;
    while at + 3 <= aux.len() {
        if aux[at..at + 2] == tag {
            return Some(at);
        }
        let t = aux[at + 2];
        let v = at + 3;
        at = match t {
            b'Z' | b'H' => v + aux[v..].iter().position(|&c| c == 0)? + 1,
            b'B' => {
                let w = aux_width(*aux.get(v)?)?;
                let n = u32::from_le_bytes(aux.get(v + 1..v + 5)?.try_into().ok()?) as usize;
                v + 5 + n * w
            }
            _ => v + aux_width(t)?,
        };
    }
    None
}

#[inline]
fn push_block(out: &mut Blocks, s: u32, e: u32) {
    if let Some(last) = out.last_mut()
        && last.1 == s
    {
        last.1 = e;
        return;
    }
    if e > s {
        out.push((s, e));
    }
}

/// Identifies a mate pair: the read name with both mates' positions, so
/// that secondary pairs of the same name stay apart. `pending` is the key a
/// record files itself under, `wanted` the one its mate files itself under.
#[derive(Debug, Clone, PartialEq, Eq, Hash, PartialOrd, Ord)]
pub struct MateKey {
    name: SmallVec<[u8; 64]>,
    /// (tid, pos) of the record filed, then of its mate.
    at: (u64, u64),
}

/// `(tid, pos)` packed to compare in coordinate order.
#[inline]
pub fn coord_key(tid: i32, pos: i32) -> u64 {
    ((tid as u32 as u64) << 32) | (pos as i64 + 1).clamp(0, u32::MAX as i64) as u64
}

impl MateKey {
    /// The key `r` files itself under, and the one its mate does.
    pub fn of(r: &ReadView<'_>) -> (MateKey, MateKey) {
        let me = coord_key(r.rec.tid, r.rec.pos);
        let mate = coord_key(r.rec.mate_tid, r.rec.mate_pos);
        let name: SmallVec<[u8; 64]> = SmallVec::from_slice(r.name());
        (
            MateKey {
                name: name.clone(),
                at: (me, mate),
            },
            MateKey {
                name,
                at: (mate, me),
            },
        )
    }

    /// Where the mate of the record filed under this key is placed.
    pub fn mate_at(&self) -> u64 {
        self.at.1
    }
}

/// Whether `r` has a mate that should appear in the BAM: paired, both
/// mates mapped and placed.
pub fn expects_mate(r: &ReadView<'_>) -> bool {
    let f = r.flag();
    f & flag::PAIRED != 0 && f & (flag::UNMAPPED | flag::MATE_UNMAPPED) == 0 && r.rec.mate_tid >= 0
}

/// Mate pairing within one batch: the first mate of a pair waits by key;
/// halves whose mate is not in the batch are handed to [`Mates`] at merge.
pub struct Pairs<T> {
    waiting: HashMap<MateKey, T>,
}

impl<T> Default for Pairs<T> {
    fn default() -> Pairs<T> {
        Pairs {
            waiting: HashMap::new(),
        }
    }
}

impl<T> Pairs<T> {
    /// Offers `r`'s half: returns the mate's half when it is waiting,
    /// otherwise files this one.
    pub fn offer(&mut self, r: &ReadView<'_>, half: T) -> Option<(T, T)> {
        let (me, wanted) = MateKey::of(r);
        match self.waiting.remove(&wanted) {
            Some(mate) => Some((mate, half)),
            None => {
                self.waiting.insert(me, half);
                None
            }
        }
    }

    /// The halves left without their mate, in a deterministic order.
    pub fn into_dangling(self) -> Vec<(MateKey, T)> {
        let mut v: Vec<(MateKey, T)> = self.waiting.into_iter().collect();
        v.sort_unstable_by(|a, b| (a.0.at, &a.0.name).cmp(&(b.0.at, &b.0.name)));
        v
    }

    pub fn len(&self) -> usize {
        self.waiting.len()
    }

    pub fn is_empty(&self) -> bool {
        self.waiting.is_empty()
    }
}

/// Mate pairing across batches, on the merging thread: halves wait by key
/// until their mate arrives, or until every record at or before the mate's
/// position has been merged without it.
pub struct Mates<T> {
    waiting: HashMap<MateKey, T>,
    /// (mate position, sequence) of each half filed, to expire them in
    /// order.
    expiry: BinaryHeap<Reverse<(u64, u64, MateKey)>>,
    filed: u64,
    /// Bytes per waiting half beyond its key, for the budget.
    pub half_bytes: usize,
}

impl<T> Mates<T> {
    pub fn new(half_bytes: usize) -> Mates<T> {
        Mates {
            waiting: HashMap::new(),
            expiry: BinaryHeap::new(),
            filed: 0,
            half_bytes,
        }
    }

    /// Takes a batch's dangling half: returns the pair when its mate was
    /// waiting, otherwise files it.
    pub fn offer(&mut self, key: MateKey, half: T) -> Option<(T, T)> {
        let wanted = MateKey {
            name: key.name.clone(),
            at: (key.at.1, key.at.0),
        };
        if let Some(mate) = self.waiting.remove(&wanted) {
            return Some((mate, half));
        }
        self.filed += 1;
        self.expiry
            .push(Reverse((key.mate_at(), self.filed, key.clone())));
        self.waiting.insert(key, half);
        None
    }

    /// The halves whose mate would sit before `frontier` (all records
    /// before it are merged): their mates never came.
    pub fn expire(&mut self, frontier: u64, mut alone: impl FnMut(T)) {
        while let Some(Reverse((at, _, _))) = self.expiry.peek() {
            if *at >= frontier {
                break;
            }
            let Reverse((_, _, key)) = self.expiry.pop().expect("peeked");
            if let Some(half) = self.waiting.remove(&key) {
                alone(half);
            }
        }
    }

    /// Every half still waiting, in filing order.
    pub fn drain(&mut self, mut alone: impl FnMut(T)) {
        let mut all: Vec<(u64, MateKey)> = self
            .expiry
            .drain()
            .map(|Reverse((_, seq, key))| (seq, key))
            .collect();
        all.sort_unstable_by_key(|(seq, _)| *seq);
        for (_, key) in all {
            if let Some(half) = self.waiting.remove(&key) {
                alone(half);
            }
        }
    }

    pub fn len(&self) -> usize {
        self.waiting.len()
    }

    pub fn is_empty(&self) -> bool {
        self.waiting.is_empty()
    }

    /// Approximate bytes held.
    pub fn bytes(&self) -> u64 {
        (self.waiting.len() as u64)
            * (std::mem::size_of::<MateKey>() as u64 * 2 + self.half_bytes as u64 + 32)
    }
}

/// Packs a UMI of `A`, `C`, `G`, `T` (up to 31 bases) into a distinct
/// integer; other UMIs get a 63-bit hash with the top bit set.
pub fn umi_key(umi: &[u8]) -> u64 {
    if umi.len() <= 31 {
        let mut v = 1u64;
        let mut ok = true;
        for &b in umi {
            let code = match b {
                b'A' => 0,
                b'C' => 1,
                b'G' => 2,
                b'T' => 3,
                _ => {
                    ok = false;
                    break;
                }
            };
            v = v << 2 | code;
        }
        if ok {
            return v;
        }
    }
    (1 << 63) | (fnv(umi) >> 1)
}

/// FNV-1a.
pub fn fnv(b: &[u8]) -> u64 {
    let mut h = 0xcbf2_9ce4_8422_2325u64;
    for &c in b {
        h ^= c as u64;
        h = h.wrapping_mul(0x100_0000_01b3);
    }
    h
}

#[cfg(test)]
mod tests {
    use super::*;

    struct NoReference;

    impl Context for NoReference {
        fn assign(&self, _: SpaceId, _: &ReadView<'_>, _: Query, _: Option<bool>) -> Hits {
            Hits::new()
        }
        fn assign_at(&self, _: SpaceId, _: &ReadView<'_>, _: u32, _: u32) -> Hits {
            Hits::new()
        }
        fn reference(&self, _: u32, _: u32, _: u32, _: &mut dyn FnMut(&[u8])) -> Result<bool> {
            Ok(false)
        }
        fn variants(&self, _: u32, _: u32, _: u32, _: &mut Vec<u32>) -> Result<()> {
            Ok(())
        }
    }

    fn raw_record(cigar: &[(u8, u32)], l_seq: usize) -> Vec<u8> {
        let mut bytes = vec![0; 32];
        bytes[4..8].copy_from_slice(&100i32.to_le_bytes());
        bytes[8] = 2;
        bytes[12..14].copy_from_slice(&(cigar.len() as u16).to_le_bytes());
        bytes[16..20].copy_from_slice(&(l_seq as i32).to_le_bytes());
        bytes[20..24].copy_from_slice(&(-1i32).to_le_bytes());
        bytes[24..28].copy_from_slice(&(-1i32).to_le_bytes());
        bytes.extend_from_slice(b"x\0");
        for &(op, len) in cigar {
            bytes.extend_from_slice(&((len << 4) | op as u32).to_le_bytes());
        }
        bytes.extend(std::iter::repeat_n(0x12, l_seq.div_ceil(2)));
        bytes.extend((0..l_seq).map(|q| 10 + q as u8));
        bytes.extend_from_slice(format!("MDZ{l_seq}\0").as_bytes());
        bytes
    }

    #[test]
    fn aligned_bases_keep_odd_offsets_qualities_and_hard_clips() {
        for lead in 0..4 {
            for len in 1..17 {
                let cigar = [
                    (op::HARD_CLIP, 3),
                    (op::SOFT_CLIP, lead),
                    (op::MATCH, len),
                    (op::HARD_CLIP, 7),
                ];
                let bytes = raw_record(&cigar, (lead + len) as usize);
                let rec = Record::parse(&bytes, 1).unwrap();
                let view = ReadView::new(0, 0, rec, &NoReference);
                assert_eq!(view.hard_clips(), (3, 7));
                let bases = view.aligned().unwrap().unwrap();
                assert_eq!(bases.len(), len as usize);
                for (i, b) in bases.iter().enumerate() {
                    let q = i + lead as usize;
                    assert_eq!(b.qpos, q as u32);
                    assert_eq!(b.ref_pos, 100 + i as u32);
                    assert_eq!(b.read, if q.is_multiple_of(2) { b'A' } else { b'C' });
                    assert_eq!(b.reference, b.read);
                    assert_eq!(b.qual, 10 + q as u8);
                }
            }
        }
        for (cigar, want) in [
            (vec![], (0, 0)),
            (vec![(op::HARD_CLIP, 5)], (5, 0)),
            (vec![(op::MATCH, 1)], (0, 0)),
        ] {
            let len = cigar
                .iter()
                .filter(|x| x.0 == op::MATCH)
                .map(|x| x.1)
                .sum::<u32>();
            let bytes = raw_record(&cigar, len as usize);
            let rec = Record::parse(&bytes, 1).unwrap();
            assert_eq!(ReadView::new(0, 0, rec, &NoReference).hard_clips(), want);
        }
    }

    #[test]
    fn spliced_walk_fetches_segments_and_matches_explicit_full_reference() {
        #[derive(Default)]
        struct Reference {
            ranges: RefCell<Vec<(u32, u32)>>,
        }
        impl Context for Reference {
            fn assign(&self, _: SpaceId, _: &ReadView<'_>, _: Query, _: Option<bool>) -> Hits {
                Hits::new()
            }
            fn assign_at(&self, _: SpaceId, _: &ReadView<'_>, _: u32, _: u32) -> Hits {
                Hits::new()
            }
            fn reference(
                &self,
                _: u32,
                s: u32,
                e: u32,
                with: &mut dyn FnMut(&[u8]),
            ) -> Result<bool> {
                self.ranges.borrow_mut().push((s, e));
                let bases: Vec<_> = (s..e).map(|p| b"ACGT"[(p % 4) as usize]).collect();
                with(&bases);
                Ok(true)
            }
            fn variants(&self, _: u32, _: u32, _: u32, out: &mut Vec<u32>) -> Result<()> {
                out.extend([103, 1_000_112]);
                Ok(())
            }
        }
        let cigar = [
            (op::HARD_CLIP, 3),
            (op::SOFT_CLIP, 2),
            (op::MATCH, 3),
            (op::INS, 1),
            (op::EQUAL, 4),
            (op::DEL, 2),
            (op::DIFF, 2),
            (op::REF_SKIP, 1_000_000),
            (op::MATCH, 5),
            (op::SOFT_CLIP, 1),
            (op::HARD_CLIP, 7),
        ];
        let bytes = raw_record(&cigar, 18);
        let ctx = Reference::default();
        for reverse in [false, true] {
            let mut record = Record::parse(&bytes, 1).unwrap();
            if reverse {
                record.flag |= flag::REVERSE;
            }
            ctx.ranges.borrow_mut().clear();
            let view = ReadView::new(0, 0, record.clone(), &ctx);
            let bases = view.aligned().unwrap().unwrap();
            assert_eq!(bases.len(), 14);
            assert_eq!(
                ctx.ranges.borrow().iter().map(|(s, e)| e - s).sum::<u32>(),
                14
            );
            for b in bases {
                assert_eq!(b.reference, b"ACGT"[(b.ref_pos % 4) as usize]);
                assert_eq!(b.qual, 10 + b.qpos as u8);
                assert_eq!(b.masked, [103, 1_000_112].contains(&b.ref_pos));
            }
            ctx.ranges.borrow_mut().clear();
            assert!(view.has_reference_context().unwrap());
            assert_eq!(
                ctx.ranges.borrow().iter().map(|(s, e)| e - s).sum::<u32>(),
                24
            );
            assert_eq!(view.reference_context().unwrap().unwrap().ranges.len(), 2);
            for b in bases {
                for p in b.ref_pos.saturating_sub(2)..=b.ref_pos.saturating_add(2) {
                    assert_eq!(
                        view.reference_context_base(p).unwrap(),
                        Some(b"ACGT"[(p % 4) as usize])
                    );
                }
            }
            assert_eq!(view.reference_context_base(500_000).unwrap(), Some(b'N'));
            let explicit = ReadView::new(0, 0, record, &ctx);
            let (start, reference) = explicit.reference_bases().unwrap().unwrap();
            assert_eq!(start, 98);
            assert_eq!(reference.len(), 1_000_020);
            let calls = ctx.ranges.borrow().len();
            assert_eq!(explicit.aligned().unwrap().unwrap(), bases);
            assert!(explicit.has_reference_context().unwrap());
            assert_eq!(
                explicit.reference_context_base(500_000).unwrap(),
                Some(b'A')
            );
            assert_eq!(ctx.ranges.borrow().len(), calls);
        }
    }

    #[test]
    fn aux_fields_are_found_by_tag() {
        let aux = b"NMi\x05\x00\x00\x00CBZAACG-1\x00UBZTTT\x00XAAzHCC\x00";
        assert_eq!(find_aux(aux, *b"NM"), Some(Aux::Int(5)));
        assert_eq!(find_aux(aux, *b"CB"), Some(Aux::Str(b"AACG-1")));
        assert_eq!(find_aux(aux, *b"UB"), Some(Aux::Str(b"TTT")));
        assert_eq!(find_aux(aux, *b"XA"), Some(Aux::Char(b'z')));
        assert_eq!(find_aux(aux, *b"QQ"), None);
        assert_eq!(locate(aux, *b"XA"), Some(24));
    }

    #[test]
    fn md_tags_give_reference_bases() {
        let mut md = Md::new(b"2A0C1^GT1");
        let read = b"ACGTAC";
        let mut got = Vec::new();
        for &r in &read[..5] {
            got.push(md.next_base(r).unwrap());
        }
        assert!(md.deletion(2));
        got.push(md.next_base(read[5]).unwrap());
        assert_eq!(&got, b"ACACAC");
        assert_eq!(md.next_base(b'A'), None);
        assert!(!Md::new(b"3").deletion(1));
    }

    #[test]
    fn umis_pack_without_collisions() {
        assert_ne!(umi_key(b"A"), umi_key(b"AA"));
        assert_ne!(umi_key(b"ACGT"), umi_key(b"ACGG"));
        assert!(umi_key(b"ACGN") >> 63 == 1);
        assert!(umi_key(b"ACGT") >> 63 == 0);
    }
}
