//! Decoding a slice into records, as htslib's `cram_decode_slice`,
//! `cram_decode_seq`, `cram_decode_aux` and `cram_decode_slice_xref` do,
//! restricted to the data series a caller needs.

use omicsweave_core::{Error, Result};

use crate::codec::Streams;
use crate::comp::{CompressionHeader, SliceHeader, ds, s};
use crate::container::{Block, content, method};
use crate::md5::md5;

/// BAM flag bits.
const PAIRED: i32 = 0x1;
const UNMAPPED: i32 = 0x4;
const MATE_UNMAPPED: i32 = 0x8;
const REVERSE: i32 = 0x10;
const MATE_REVERSE: i32 = 0x20;
const READ1: i32 = 0x40;

/// CRAM record flags (CF).
const PRESERVE_QUAL: i32 = 1;
const DETACHED: i32 = 2;
const MATE_DOWNSTREAM: i32 = 4;
const NO_SEQ: i32 = 8;
const EXPLICIT_TLEN: i32 = 16;

/// Mate flags (MF).
const M_REVERSE: i32 = 1;
const M_UNMAP: i32 = 2;

/// CIGAR operations.
const C_MATCH: u32 = 0;
const C_INS: u32 = 1;
const C_DEL: u32 = 2;
const C_SKIP: u32 = 3;
const C_SOFT: u32 = 4;
const C_HARD: u32 = 5;
const C_PAD: u32 = 6;

/// Longest read name BAM holds.
pub(crate) const MAX_NAME: i32 = 254;

/// Reference bases for a slice, as a caller provides them.
pub trait References {
    /// The length of reference `tid`.
    fn ref_len(&mut self, tid: usize) -> Result<i64>;

    /// Writes the bases `start..=end` (1-based) of reference `tid`,
    /// uppercase, into `out`; `end` is within the reference.
    fn bases(&mut self, tid: usize, start: i64, end: i64, out: &mut Vec<u8>) -> Result<()>;
}

/// A slice as stored: its header and its blocks.
pub(crate) struct SliceInput<'a> {
    pub header: SliceHeader,
    pub blocks: Vec<Block<'a>>,
}

impl<'a> SliceInput<'a> {
    /// Reads the slice at `buf[*at..]` (`cram_read_slice`), advancing `at`.
    pub(crate) fn read(buf: &'a [u8], at: &mut usize, scratch: &mut Vec<u8>) -> Result<Self> {
        let head = Block::parse(buf, at)?;
        if head.content_type != content::MAPPED_SLICE
            && head.content_type != content::UNMAPPED_SLICE
        {
            return Err(Error::invalid(
                "CRAM slice",
                "a slice does not start with its header",
            ));
        }
        let data = if head.method == method::RAW {
            head.data
        } else {
            head.uncompress(scratch, true)?
        };
        let header = SliceHeader::parse(head.content_type, data)?;
        let mut blocks = Vec::with_capacity(header.num_blocks as usize);
        for _ in 0..header.num_blocks {
            blocks.push(Block::parse(buf, at)?);
        }
        Ok(SliceInput { header, blocks })
    }
}

/// What a slice decode works from, besides the slice.
pub(crate) struct Decode<'c> {
    pub comp: &'c CompressionHeader,
    /// Reference lengths of the SAM header, by reference id.
    pub ref_lens: &'c [i64],
    /// The read group htslib treats as none (a last `@RG ID:UNKNOWN`).
    pub unknown_rg: i32,
    /// The data series asked for, before [`CompressionHeader::closure`].
    pub wanted: u32,
    /// Reconstruct bases from the reference.
    pub bases: bool,
    /// Generate MD and NM where the record has none.
    pub decode_md: bool,
    /// Resolve mates within the slice; when false, mate fields are zero.
    pub mates: bool,
    /// Decode as htslib does when it reads every field, which loads the
    /// reference of every slice that has one (this decides where an
    /// alignment running past its reference ends); false as when it
    /// builds an index.
    pub holds_refs: bool,
    /// The reference a query reads (-2 for all): multi-reference slices
    /// load only its sequence.
    pub range_ref: i32,
    /// Check block CRC32s and the slice's reference MD5.
    pub verify: bool,
}

/// One decoded record (`cram_record`). Ranges index the buffers of
/// [`Decoded`].
#[derive(Debug, Clone, Default)]
pub(crate) struct Rec {
    pub flags: i32,
    pub cram_flags: i32,
    pub ref_id: i32,
    pub len: i32,
    /// 1-based alignment start and end.
    pub apos: i64,
    pub aend: i64,
    pub rg: i32,
    pub name: u32,
    pub name_len: u32,
    pub mate_flags: i32,
    pub mate_ref_id: i32,
    pub mate_pos: i64,
    pub mate_line: i32,
    pub tlen: i64,
    pub explicit_tlen: Option<i64>,
    pub mqual: i32,
    pub cigar: u32,
    pub ncigar: u32,
    pub seq: u32,
    pub qual: u32,
    pub aux: u32,
    pub aux_size: u32,
    pub has_cg: bool,
}

/// A decoded slice.
#[derive(Debug, Default)]
pub(crate) struct Decoded {
    pub recs: Vec<Rec>,
    pub record_counter: i64,
    pub names: Vec<u8>,
    pub cigars: Vec<u32>,
    pub seqs: Vec<u8>,
    pub quals: Vec<u8>,
    pub aux: Vec<u8>,
    /// Whether bases and qualities were decoded.
    pub has_seq: bool,
    pub has_qual: bool,
}

/// `fd->L1`: A, C, G, T to 0-3, anything else to 4.
#[inline]
fn l1(c: u8) -> usize {
    match c {
        b'A' | b'a' => 0,
        b'C' | b'c' => 1,
        b'G' | b'g' => 2,
        b'T' | b't' => 3,
        _ => 4,
    }
}

fn bad(detail: impl Into<String>) -> Error {
    Error::invalid("CRAM record", detail)
}

fn beyond() -> Error {
    bad("CIGAR extends beyond the slice's reference extent")
}

/// `append_uint32`: decimal digits.
fn push_uint(out: &mut Vec<u8>, v: u32) {
    let mut buf = [0u8; 10];
    let mut i = buf.len();
    let mut v = v;
    loop {
        i -= 1;
        buf[i] = b'0' + (v % 10) as u8;
        v /= 10;
        if v == 0 {
            break;
        }
    }
    out.extend_from_slice(&buf[i..]);
}

/// The reference a record reads: bases from 1-based `start`, with the
/// slice's extent end.
struct RefView<'r> {
    bases: &'r [u8],
    start: i64,
    end: i64,
}

impl RefView<'_> {
    /// The base at 0-based `pos`.
    #[inline]
    fn at(&self, pos: i64) -> u8 {
        self.bases[(pos - self.start + 1) as usize]
    }

    fn range(&self, pos: i64, n: i64) -> &[u8] {
        let i = (pos - self.start + 1) as usize;
        &self.bases[i..i + n as usize]
    }
}

/// Per-record state of `cram_decode_seq`.
struct SeqState<'x, 'r> {
    reference: Option<&'x RefView<'r>>,
    /// Whether htslib would hold a reference for this record (it changes
    /// how an alignment past the reference end ends).
    holds_ref: bool,
    ref_len: i64,
}

impl Decode<'_> {
    /// Decodes every record of `slice`.
    pub(crate) fn run(
        &self,
        slice: &SliceInput<'_>,
        mut refs: Option<&mut dyn References>,
        out: &mut Decoded,
    ) -> Result<()> {
        let comp = self.comp;
        let h = &slice.header;
        out.recs.clear();
        out.names.clear();
        out.cigars.clear();
        out.seqs.clear();
        out.quals.clear();
        out.aux.clear();
        out.record_counter = h.record_counter;
        let present: Vec<(u8, i32)> = slice
            .blocks
            .iter()
            .map(|b| (b.content_type, b.content_id))
            .collect();
        let d = comp.closure(self.wanted, &present);
        out.has_seq = self.bases;
        out.has_qual = d & ds::QS != 0;
        let Some(core_block) = slice.blocks.first() else {
            return Err(Error::invalid("CRAM slice", "no blocks"));
        };
        if core_block.content_type != content::CORE {
            return Err(Error::invalid("CRAM slice", "the first block is not CORE"));
        }
        let mut core_buf = Vec::new();
        let core = core_block.uncompress(&mut core_buf, self.verify)?;
        let mut st = Streams::new(core, &slice.blocks, self.verify);

        // The slice's reference.
        let ref_id = h.ref_seq_id;
        let embed = h.ref_base_id >= 0;
        let mut ref_bytes: Vec<u8> = Vec::new();
        let mut ref_extent: Option<(i64, i64)> = None;
        let mut single_holds = false;
        if ref_id >= 0 {
            if embed {
                let b = st.whole(h.ref_base_id)?.ok_or_else(|| {
                    Error::invalid("CRAM slice", "embedded reference block missing")
                })?;
                if h.ref_seq_span > b.len() as i64 {
                    return Err(Error::invalid(
                        "CRAM slice",
                        "embedded reference is too small",
                    ));
                }
                ref_bytes.extend_from_slice(b);
                ref_extent = Some((h.ref_seq_start, h.ref_seq_end()));
                single_holds = true;
            } else if !comp.no_ref {
                single_holds = self.holds_refs;
                if self.bases {
                    let refs = refs.as_deref_mut().ok_or_else(|| {
                        Error::invalid(
                            "CRAM slice",
                            "the slice needs its reference sequence; give a reference",
                        )
                    })?;
                    let len = refs.ref_len(ref_id as usize)?;
                    let start = h.ref_seq_start.max(0);
                    let mut end = h.ref_seq_end();
                    if h.ref_seq_start >= 1 {
                        let fetch_end = end.min(len).max(h.ref_seq_start);
                        refs.bases(ref_id as usize, h.ref_seq_start, fetch_end, &mut ref_bytes)?;
                    }
                    if end > len {
                        end = len;
                    }
                    ref_extent = Some((start, end));
                }
            }
        }
        if self.bases && self.verify && ref_id >= 0 && h.md5 != [0; 16] && !comp.no_ref {
            let digest = match ref_extent {
                Some((rs, re)) => {
                    let mut start = h.ref_seq_start - rs;
                    if h.ref_seq_start < rs {
                        start = 0;
                    }
                    let have = re - rs + 1;
                    let mut len = if h.ref_seq_span <= have {
                        h.ref_seq_span
                    } else {
                        have
                    };
                    if start + len > have {
                        len = have - start;
                    }
                    let (start, len) = (start.max(0) as usize, len.max(0) as usize);
                    md5(ref_bytes.get(start..start + len).unwrap_or_default())
                }
                None => [0; 16],
            };
            if digest != h.md5 {
                return Err(Error::invalid(
                    "CRAM slice",
                    "MD5 checksum reference mismatch; check the reference given is correct",
                ));
            }
        }
        let single_view = ref_extent.map(|(start, end)| RefView {
            bases: &ref_bytes,
            start,
            end,
        });

        // Multi-reference slices: whole references, by id.
        let mut multi: Vec<(i32, Vec<u8>, i64)> = Vec::new();
        let mut last_ref_id = -9;
        let mut multi_at: Option<usize> = None;

        let n = h.num_records.max(0) as usize;
        out.recs.reserve(n);
        let mut last_apos = h.ref_seq_start;
        for rec in 0..n {
            let mut cr = Rec {
                mate_line: -1,
                mate_ref_id: -1,
                ..Rec::default()
            };
            let bf = if d & ds::BF != 0 {
                let v = comp
                    .codec(s::BF)
                    .ok_or_else(|| bad("no BF codec"))?
                    .int(&mut st)?;
                if !(0..0x1000).contains(&v) {
                    return Err(bad(format!("flag {v}")));
                }
                v
            } else {
                UNMAPPED
            };
            cr.flags = bf;
            let cf = if d & ds::CF != 0 {
                comp.codec(s::CF)
                    .ok_or_else(|| bad("no CF codec"))?
                    .int(&mut st)?
            } else {
                0
            };
            cr.cram_flags = cf;
            if ref_id == -2 {
                if d & ds::RI != 0 {
                    cr.ref_id = comp
                        .codec(s::RI)
                        .ok_or_else(|| bad("no RI codec"))?
                        .int(&mut st)?;
                    if cr.ref_id < -1 || cr.ref_id as i64 >= self.ref_lens.len() as i64 {
                        return Err(bad(format!("unknown reference id {}", cr.ref_id)));
                    }
                    if cr.ref_id >= 0 && cr.ref_id != last_ref_id {
                        multi_at = None;
                        let need =
                            !comp.no_ref && (self.range_ref == -2 || cr.ref_id == self.range_ref);
                        if need && self.bases {
                            let at = match multi.iter().position(|m| m.0 == cr.ref_id) {
                                Some(at) => at,
                                None => {
                                    let refs = refs.as_deref_mut().ok_or_else(|| {
                                        Error::invalid(
                                            "CRAM slice",
                                            "the slice needs its reference sequences; give a reference",
                                        )
                                    })?;
                                    let len = refs.ref_len(cr.ref_id as usize)?;
                                    let mut b = Vec::new();
                                    if len > 0 {
                                        refs.bases(cr.ref_id as usize, 1, len, &mut b)?;
                                    }
                                    multi.push((cr.ref_id, b, len));
                                    multi.len() - 1
                                }
                            };
                            multi_at = Some(at);
                        }
                        last_ref_id = cr.ref_id;
                    }
                } else {
                    cr.ref_id = -1;
                }
            } else {
                cr.ref_id = ref_id;
            }
            if cr.ref_id < -1 || cr.ref_id as i64 >= self.ref_lens.len() as i64 {
                return Err(bad(format!("unknown reference id {}", cr.ref_id)));
            }
            if d & ds::RL != 0 {
                cr.len = comp
                    .codec(s::RL)
                    .ok_or_else(|| bad("no RL codec"))?
                    .int(&mut st)?;
                if cr.len < 0 {
                    return Err(bad("negative read length"));
                }
            }
            if d & ds::AP != 0 {
                let v = comp
                    .codec(s::AP)
                    .ok_or_else(|| bad("no AP codec"))?
                    .int(&mut st)? as i64;
                cr.apos = if comp.ap_delta { v + last_apos } else { v };
                last_apos = cr.apos;
                if h.ref_seq_id >= 0 && cr.apos < h.ref_seq_start {
                    return Err(bad("alignment position before the slice"));
                }
            }
            cr.rg = -1;
            if d & ds::RG != 0 {
                cr.rg = comp
                    .codec(s::RG)
                    .ok_or_else(|| bad("no RG codec"))?
                    .int(&mut st)?;
                if cr.rg == self.unknown_rg {
                    cr.rg = -1;
                }
            }
            let name = |st: &mut Streams<'_>, cr: &mut Rec, out: &mut Decoded| -> Result<()> {
                cr.name = out.names.len() as u32;
                if d & ds::RN != 0 {
                    let n = comp.codec(s::RN).ok_or_else(|| bad("no RN codec"))?.block(
                        st,
                        1,
                        &mut out.names,
                    )?;
                    if n > MAX_NAME {
                        return Err(bad("read name too long"));
                    }
                    cr.name_len = n as u32;
                }
                Ok(())
            };
            if comp.read_names_included {
                name(&mut st, &mut cr, out)?;
            }
            cr.explicit_tlen = None;
            let mut tlen: Option<i64> = None;
            if d & ds::CF != 0 && cf & DETACHED != 0 {
                if d & ds::MF != 0 {
                    cr.mate_flags = comp
                        .codec(s::MF)
                        .ok_or_else(|| bad("no MF codec"))?
                        .int(&mut st)?;
                }
                if !comp.read_names_included {
                    name(&mut st, &mut cr, out)?;
                }
                if d & ds::NS != 0 {
                    cr.mate_ref_id = comp
                        .codec(s::NS)
                        .ok_or_else(|| bad("no NS codec"))?
                        .int(&mut st)?;
                    if cr.mate_ref_id < -1 || cr.mate_ref_id as i64 >= self.ref_lens.len() as i64 {
                        return Err(bad(format!("unknown mate reference id {}", cr.mate_ref_id)));
                    }
                }
                if d & ds::NP != 0 {
                    cr.mate_pos = comp
                        .codec(s::NP)
                        .ok_or_else(|| bad("no NP codec"))?
                        .int(&mut st)? as i64;
                }
                if d & ds::TS != 0 {
                    tlen = Some(
                        comp.codec(s::TS)
                            .ok_or_else(|| bad("no TS codec"))?
                            .int(&mut st)? as i64,
                    );
                }
            } else if d & ds::CF != 0 && cf & MATE_DOWNSTREAM != 0 {
                if d & ds::NF != 0 {
                    let nf = comp
                        .codec(s::NF)
                        .ok_or_else(|| bad("no NF codec"))?
                        .int(&mut st)?;
                    cr.mate_line = nf.wrapping_add(rec as i32 + 1);
                    cr.mate_ref_id = -1;
                    cr.mate_pos = 0;
                } else {
                    cr.mate_flags = 0;
                }
                if cf & EXPLICIT_TLEN != 0 {
                    if d & ds::TS != 0 {
                        cr.explicit_tlen = Some(
                            comp.codec(s::TS)
                                .ok_or_else(|| bad("no TS codec"))?
                                .int(&mut st)? as i64,
                        );
                    } else {
                        cr.mate_flags = 0;
                    }
                }
            } else if d & ds::CF != 0 && cf & EXPLICIT_TLEN != 0 {
                if d & ds::TS != 0 {
                    cr.explicit_tlen = Some(
                        comp.codec(s::TS)
                            .ok_or_else(|| bad("no TS codec"))?
                            .int(&mut st)? as i64,
                    );
                } else {
                    cr.mate_flags = 0;
                }
            } else {
                cr.mate_flags = 0;
            }
            // `None` stands for htslib's INT64_MIN "not yet known".
            cr.tlen = tlen.unwrap_or(i64::MIN);

            // Tags.
            let (mut has_md, mut has_nm) = (0i64, 0i64);
            self.decode_aux(&mut st, d, &mut cr, out, &mut has_md, &mut has_nm)?;

            cr.seq = out.seqs.len() as u32;
            cr.qual = out.quals.len() as u32;
            let len = cr.len as usize;
            if d & ds::RL != 0 {
                if self.bases {
                    out.seqs.resize(out.seqs.len() + len, b'=');
                }
                if out.has_qual {
                    out.quals.resize(out.quals.len() + len, 0);
                }
            }

            // The reference of this record.
            let (view, holds) = if ref_id == -2 {
                let need = !comp.no_ref
                    && self.holds_refs
                    && cr.ref_id >= 0
                    && (self.range_ref == -2 || cr.ref_id == self.range_ref);
                let view = multi_at.map(|at| RefView {
                    bases: &multi[at].1,
                    start: 1,
                    end: multi[at].2,
                });
                (view, need)
            } else {
                (
                    single_view.as_ref().map(|v| RefView {
                        bases: v.bases,
                        start: v.start,
                        end: v.end,
                    }),
                    single_holds,
                )
            };
            let sq = SeqState {
                reference: view.as_ref(),
                holds_ref: holds,
                ref_len: if cr.ref_id >= 0 {
                    self.ref_lens[cr.ref_id as usize]
                } else {
                    0
                },
            };
            if view.is_some() && self.bases {
                let at = cr.seq as usize;
                out.seqs[at..at + len].fill(0);
            }

            if bf & UNMAPPED == 0 {
                if d & ds::AP != 0 && cr.apos <= 0 {
                    return Err(bad(format!(
                        "alignment position {} without the unmapped flag",
                        cr.apos
                    )));
                }
                if d & (ds::SEQ | ds::MQ) != 0 {
                    self.decode_seq(&mut st, d, &mut cr, out, &sq, has_md, has_nm)?;
                } else {
                    cr.aend = cr.apos;
                }
            } else {
                cr.aend = cr.apos;
                if d & ds::BA != 0 && cr.len > 0 {
                    let codec = comp.codec(s::BA).ok_or_else(|| bad("no BA codec"))?;
                    let at = cr.seq as usize;
                    let seq = self.bases.then(|| &mut out.seqs[at..at + len]);
                    codec.bytes(&mut st, len, seq)?;
                }
                if d & ds::CF != 0 && cf & PRESERVE_QUAL != 0 {
                    if d & ds::QS != 0 {
                        let codec = comp.codec(s::QS).ok_or_else(|| bad("no QS codec"))?;
                        let at = cr.qual as usize;
                        let qual = out.has_qual.then(|| &mut out.quals[at..at + len]);
                        codec.bytes(&mut st, len, qual)?;
                    }
                } else if out.has_qual {
                    let at = cr.qual as usize;
                    out.quals[at..at + len].fill(255);
                }
            }
            if !comp.qs_seq_orient && d & ds::QS != 0 && cr.flags & REVERSE != 0 {
                let at = cr.qual as usize;
                out.quals[at..at + len].reverse();
            }
            out.recs.push(cr);
        }
        if self.mates {
            xref(&mut out.recs)?;
        } else {
            for cr in &mut out.recs {
                cr.tlen = 0;
                cr.mate_pos = 0;
                cr.mate_ref_id = -1;
            }
        }
        Ok(())
    }

    /// `cram_decode_aux`: the record's tags in tag-line order. `has_md` /
    /// `has_nm` become 0 when absent, positive when stored, and minus the
    /// value position of a `*` placeholder.
    fn decode_aux(
        &self,
        st: &mut Streams<'_>,
        d: u32,
        cr: &mut Rec,
        out: &mut Decoded,
        has_md: &mut i64,
        has_nm: &mut i64,
    ) -> Result<()> {
        let comp = self.comp;
        cr.has_cg = false;
        if d & (ds::TL | ds::AUX) == 0 {
            cr.aux = 0;
            cr.aux_size = 0;
            return Ok(());
        }
        let tl = comp
            .codec(crate::comp::s::TL)
            .ok_or_else(|| bad("no TL codec"))?
            .int(st)?;
        if tl < 0 || tl as usize >= comp.tag_lines.len() {
            return Err(bad(format!("tag line {tl}")));
        }
        cr.aux = out.aux.len() as u32;
        cr.aux_size = 0;
        let line = &comp.tag_lines[tl as usize];
        if d & ds::AUX == 0 {
            cr.has_cg = line.iter().any(|k| k[0] == b'C' && k[1] == b'G');
            return Ok(());
        }
        for key in line {
            let at = (out.aux.len() + 3) as i64;
            if key[0] == b'M' && key[1] == b'D' {
                *has_md = if key[2] == b'*' { -at } else { at };
            }
            if key[0] == b'N' && key[1] == b'M' {
                *has_nm = if key[2] == b'*' { -at } else { at };
            }
            if key[0] == b'C' && key[1] == b'G' {
                cr.has_cg = true;
            }
            let id = (key[0] as u32) << 16 | (key[1] as u32) << 8 | key[2] as u32;
            let codec = comp
                .tags
                .get(&id)
                .ok_or_else(|| bad(format!("no codec for tag {}", String::from_utf8_lossy(key))))?;
            out.aux.extend_from_slice(key);
            let n = codec.block(st, 1, &mut out.aux)?;
            cr.aux_size += n as u32 + 3;
            if key == b"cFC" && n == 1 {
                let flags = *out.aux.last().unwrap();
                out.aux.truncate(out.aux.len() - 4);
                cr.aux_size -= 4;
                if flags & 1 != 0 && *has_md == 0 {
                    *has_md = 1;
                }
                if flags & 2 != 0 && *has_nm == 0 {
                    *has_nm = 1;
                }
            }
        }
        Ok(())
    }

    /// `cram_decode_seq`: bases, qualities, CIGAR, alignment end, mapping
    /// quality, and the MD / NM htslib generates.
    #[allow(clippy::too_many_arguments)]
    fn decode_seq(
        &self,
        st: &mut Streams<'_>,
        d: u32,
        cr: &mut Rec,
        out: &mut Decoded,
        sq: &SeqState<'_, '_>,
        has_md: i64,
        mut has_nm: i64,
    ) -> Result<()> {
        let comp = self.comp;
        let cf = cr.cram_flags;
        let len = cr.len as i64;
        let reference = sq.reference;
        let has_ref = reference.is_some();
        let mut decode_md = has_ref
            && cr.ref_id >= 0
            && ((self.decode_md && has_md == 0) || has_md < 0)
            && d & ds::AUX != 0;
        let mut decode_nm = has_ref
            && cr.ref_id >= 0
            && ((self.decode_md && has_nm == 0) || has_nm < 0)
            && d & ds::AUX != 0;
        if cf & NO_SEQ != 0 {
            decode_md = false;
            decode_nm = false;
        }
        let seq_at = cr.seq as usize;
        let qual_at = cr.qual as usize;
        let seq_on = self.bases;
        let qual_on = out.has_qual;
        if d & ds::QS != 0 && cf & PRESERVE_QUAL == 0 && qual_on {
            out.quals[qual_at..qual_at + len as usize].fill(255);
        }
        let mut md = Vec::new();
        let mut md_dist: i32 = 0;
        let mut nm: u32 = 0;
        let add_md = |md: &mut Vec<u8>, c: u8, md_dist: &mut i32| {
            if decode_md {
                push_uint(md, *md_dist as u32);
                md.push(c);
                *md_dist = 0;
            }
        };

        let fn_ = if d & ds::FN != 0 {
            comp.codec(s::FN)
                .ok_or_else(|| bad("no FN codec"))?
                .int(st)?
        } else {
            0
        };
        let mut ref_pos = cr.apos - 1;
        cr.cigar = out.cigars.len() as u32;
        let cigar = &mut out.cigars;
        let mut cig_len: u32 = 0;
        let mut cig_op = C_MATCH;
        let mut seq_pos: i32 = 1;
        let mut prev_pos: i32 = 0;
        let ref_len = sq.ref_len;
        macro_rules! flush_if {
            ($op:expr) => {
                if cig_len != 0 && cig_op != $op {
                    cigar.push(cig_len << 4 | cig_op);
                    cig_len = 0;
                }
            };
        }

        if d & (ds::FC | ds::FP) != 0 {
            if fn_ > 0
                && ((d & ds::FC != 0 && comp.codec(s::FC).is_none())
                    || (d & ds::FP != 0 && comp.codec(s::FP).is_none()))
            {
                return Err(bad("no feature codecs"));
            }
            for _ in 0..fn_ {
                let mut op = 0u8;
                if d & ds::FC != 0 {
                    let mut b = [0u8];
                    comp.codec(s::FC).unwrap().bytes(st, 1, Some(&mut b))?;
                    op = b[0];
                }
                if d & ds::FP == 0 {
                    continue;
                }
                let pos = comp.codec(s::FP).unwrap().int(st)?.wrapping_add(prev_pos);
                if pos <= 0 {
                    return Err(bad(format!("feature position {pos} before the read")));
                }
                if cr.len != 0 && pos > cr.len {
                    let valid_end = if matches!(op, b'N' | b'P' | b'H' | b'D') {
                        cr.len + 1
                    } else {
                        cr.len
                    };
                    if pos > valid_end {
                        return Err(bad(format!("feature position {pos} after the read")));
                    }
                }
                if pos > seq_pos {
                    let frag = (pos - seq_pos) as i64;
                    if let (Some(r), true) = (reference, cr.ref_id >= 0) {
                        if ref_pos + frag > ref_len {
                            let rlen = ref_len - ref_pos;
                            if rlen > 0 {
                                if ref_pos + rlen > r.end {
                                    return Err(beyond());
                                }
                                if cr.len != 0 && seq_on {
                                    let at = seq_at + seq_pos as usize - 1;
                                    out.seqs[at..at + rlen as usize]
                                        .copy_from_slice(r.range(ref_pos, rlen));
                                    if frag - rlen > 0 {
                                        out.seqs[at + rlen as usize..at + frag as usize].fill(b'N');
                                    }
                                }
                            } else if cr.len != 0 && seq_on {
                                let at = seq_at + seq_pos as usize - 1;
                                out.seqs[at..seq_at + len as usize].fill(b'N');
                            }
                            if md_dist >= 0 {
                                md_dist += frag as i32;
                            }
                        } else {
                            if ref_pos + frag > r.end {
                                return Err(beyond());
                            }
                            let refp = r.range(ref_pos, frag);
                            if decode_md || decode_nm {
                                if memchr::memchr(b'N', refp).is_some() {
                                    for &base in refp {
                                        if base == b'N' {
                                            add_md(&mut md, b'N', &mut md_dist);
                                            nm += 1;
                                        } else {
                                            md_dist += 1;
                                        }
                                    }
                                } else {
                                    md_dist += frag as i32;
                                }
                            }
                            if cr.len != 0 && seq_on {
                                let at = seq_at + seq_pos as usize - 1;
                                out.seqs[at..at + frag as usize].copy_from_slice(refp);
                            }
                        }
                    }
                    flush_if!(C_MATCH);
                    cig_op = C_MATCH;
                    cig_len = cig_len.wrapping_add(frag as u32);
                    ref_pos += frag;
                    seq_pos = pos;
                }
                prev_pos = pos;
                if d & ds::FC == 0 {
                    break;
                }
                let room = if cr.len != 0 { cr.len - (pos - 1) } else { 1 };
                let seq_out = || {
                    (cr.len != 0 && seq_on)
                        .then(|| (seq_at + pos as usize - 1, seq_at + len as usize))
                };
                match op {
                    b'S' => {
                        if cig_len != 0 {
                            cigar.push(cig_len << 4 | cig_op);
                            cig_len = 0;
                        }
                        if d & ds::SC != 0 {
                            let n = match comp.codec(s::SC) {
                                Some(c) => {
                                    let w = seq_out();
                                    let o = w.map(|(a, b)| &mut out.seqs[a..b]);
                                    c.array(st, room, o)?
                                }
                                None => {
                                    if let Some((a, _)) = seq_out() {
                                        out.seqs[a] = b'N';
                                    }
                                    1
                                }
                            };
                            cigar.push((n as u32) << 4 | C_SOFT);
                            cig_op = C_SOFT;
                            seq_pos += n;
                        }
                    }
                    b'X' => {
                        flush_if!(C_MATCH);
                        if d & ds::BS != 0 {
                            let c = comp.codec(s::BS).ok_or_else(|| bad("no BS codec"))?;
                            let mut b = [0u8];
                            c.bytes(st, 1, Some(&mut b))?;
                            let code = b[0] as usize;
                            if code > 3 {
                                return Err(bad(format!("substitution code {code}")));
                            }
                            let no_ref = cr.ref_id < 0 || ref_pos >= ref_len || !has_ref;
                            if no_ref {
                                if pos - 1 < cr.len && seq_on {
                                    out.seqs[seq_at + pos as usize - 1] = comp.sub_matrix[4][code];
                                }
                                if decode_md || decode_nm {
                                    if md_dist >= 0 && decode_md {
                                        push_uint(&mut md, md_dist as u32);
                                    }
                                    md_dist = -1;
                                    nm = nm.wrapping_sub(1);
                                }
                            } else {
                                let r = reference.unwrap();
                                let call = if ref_pos < r.end { r.at(ref_pos) } else { b'N' };
                                if pos - 1 < cr.len && seq_on {
                                    out.seqs[seq_at + pos as usize - 1] =
                                        comp.sub_matrix[l1(call)][code];
                                }
                                add_md(&mut md, call, &mut md_dist);
                            }
                        }
                        cig_op = C_MATCH;
                        nm = nm.wrapping_add(1);
                        cig_len += 1;
                        seq_pos += 1;
                        ref_pos += 1;
                    }
                    b'D' => {
                        flush_if!(C_DEL);
                        if d & ds::DL != 0 {
                            let n = comp
                                .codec(s::DL)
                                .ok_or_else(|| bad("no DL codec"))?
                                .int(st)?;
                            if n < 0 {
                                return Err(beyond());
                            }
                            let n64 = n as i64;
                            if decode_md || decode_nm {
                                let r = reference.unwrap();
                                if ref_pos + n64 > r.end {
                                    return Err(beyond());
                                }
                                if md_dist >= 0 && decode_md {
                                    push_uint(&mut md, md_dist as u32);
                                }
                                if ref_pos + n64 <= ref_len {
                                    if decode_md {
                                        md.push(b'^');
                                        md.extend_from_slice(r.range(ref_pos, n64));
                                        md_dist = 0;
                                    }
                                    nm = nm.wrapping_add(n as u32);
                                } else {
                                    let dlen = if ref_len >= ref_pos {
                                        if decode_md {
                                            md.push(b'^');
                                            md.extend_from_slice(
                                                r.range(ref_pos, ref_len - ref_pos),
                                            );
                                            push_uint(&mut md, 0);
                                        }
                                        n64 - (ref_len - ref_pos)
                                    } else {
                                        n64
                                    };
                                    nm = nm.wrapping_add((n64 - dlen) as u32);
                                    md_dist = -1;
                                }
                            }
                            cig_op = C_DEL;
                            cig_len = cig_len.wrapping_add(n as u32);
                            ref_pos += n64;
                        }
                    }
                    b'I' => {
                        flush_if!(C_INS);
                        if d & ds::IN != 0 {
                            let c = comp.codec(s::IN).ok_or_else(|| bad("no IN codec"))?;
                            let w = seq_out();
                            let o = w.map(|(a, b)| &mut out.seqs[a..b]);
                            let n = c.array(st, room, o)?;
                            cig_op = C_INS;
                            cig_len = cig_len.wrapping_add(n as u32);
                            seq_pos += n;
                            nm = nm.wrapping_add(n as u32);
                        }
                    }
                    b'i' => {
                        flush_if!(C_INS);
                        if d & ds::BA != 0 {
                            let c = comp.codec(s::BA).ok_or_else(|| bad("no BA codec"))?;
                            let w = seq_out();
                            let o = w.map(|(a, _)| &mut out.seqs[a..a + 1]);
                            c.bytes(st, 1, o)?;
                        }
                        cig_op = C_INS;
                        cig_len += 1;
                        seq_pos += 1;
                        nm = nm.wrapping_add(1);
                    }
                    b'b' => {
                        let mut n = if cr.len != 0 { cr.len - (pos - 1) } else { 1 };
                        flush_if!(C_MATCH);
                        if d & ds::BB != 0 {
                            let c = comp.codec(s::BB).ok_or_else(|| bad("no BB codec"))?;
                            let w = seq_out();
                            let o = w.map(|(a, b)| &mut out.seqs[a..b]);
                            n = c.array(st, n, o)?;
                            if decode_md || decode_nm {
                                if md_dist >= 0 && decode_md {
                                    push_uint(&mut md, md_dist as u32);
                                }
                                let mut x = 0i64;
                                while x < n as i64 {
                                    if x != 0 && decode_md {
                                        push_uint(&mut md, 0);
                                    }
                                    if ref_pos + x >= ref_len || !has_ref {
                                        break;
                                    }
                                    if decode_md {
                                        let r = reference.unwrap();
                                        if ref_pos + x >= r.end {
                                            return Err(beyond());
                                        }
                                        md.push(r.at(ref_pos + x));
                                    }
                                    x += 1;
                                }
                                nm = nm.wrapping_add(x as u32);
                                md_dist = 0;
                            }
                        }
                        cig_op = C_MATCH;
                        cig_len = cig_len.wrapping_add(n as u32);
                        seq_pos += n;
                        ref_pos += n as i64;
                    }
                    b'q' => {
                        let mut n = if cr.len != 0 { cr.len - (pos - 1) } else { 1 };
                        flush_if!(C_MATCH);
                        if d & ds::QQ != 0 {
                            let c = comp.codec(s::QQ).ok_or_else(|| bad("no QQ codec"))?;
                            if d & ds::QS != 0
                                && cf & PRESERVE_QUAL == 0
                                && cr.len > 0
                                && qual_on
                                && out.quals[qual_at] == 255
                            {
                                out.quals[qual_at..qual_at + len as usize].fill(30);
                            }
                            let o = (cr.len != 0 && qual_on).then(|| {
                                &mut out.quals[qual_at + pos as usize - 1..qual_at + len as usize]
                            });
                            n = c.array(st, n, o)?;
                        }
                        let _ = n;
                        cig_op = C_MATCH;
                    }
                    b'B' => {
                        flush_if!(C_MATCH);
                        if d & ds::BA != 0 {
                            let c = comp.codec(s::BA).ok_or_else(|| bad("no BA codec"))?;
                            let w = seq_out();
                            let o = w.map(|(a, _)| &mut out.seqs[a..a + 1]);
                            c.bytes(st, 1, o)?;
                            if decode_md || decode_nm {
                                if md_dist >= 0 && decode_md {
                                    push_uint(&mut md, md_dist as u32);
                                }
                                if ref_pos >= ref_len || !has_ref {
                                    md_dist = -1;
                                } else {
                                    if decode_md {
                                        let r = reference.unwrap();
                                        if ref_pos >= r.end {
                                            return Err(beyond());
                                        }
                                        md.push(r.at(ref_pos));
                                    }
                                    nm = nm.wrapping_add(1);
                                    md_dist = 0;
                                }
                            }
                        }
                        if d & ds::QS != 0 {
                            let c = comp.codec(s::QS).ok_or_else(|| bad("no QS codec"))?;
                            if cf & PRESERVE_QUAL == 0
                                && cr.len > 0
                                && qual_on
                                && out.quals[qual_at] == 255
                            {
                                out.quals[qual_at..qual_at + len as usize].fill(30);
                            }
                            let o = (cr.len != 0 && qual_on).then(|| {
                                let a = qual_at + pos as usize - 1;
                                &mut out.quals[a..a + 1]
                            });
                            c.bytes(st, 1, o)?;
                        }
                        cig_op = C_MATCH;
                        cig_len += 1;
                        seq_pos += 1;
                        ref_pos += 1;
                    }
                    b'Q' => {
                        if d & ds::QS != 0 {
                            let c = comp.codec(s::QS).ok_or_else(|| bad("no QS codec"))?;
                            if cf & PRESERVE_QUAL == 0
                                && cr.len > 0
                                && qual_on
                                && out.quals[qual_at] == 255
                            {
                                out.quals[qual_at..qual_at + len as usize].fill(30);
                            }
                            let o = (cr.len != 0 && qual_on).then(|| {
                                let a = qual_at + pos as usize - 1;
                                &mut out.quals[a..a + 1]
                            });
                            c.bytes(st, 1, o)?;
                        }
                    }
                    b'H' => {
                        flush_if!(C_HARD);
                        if d & ds::HC != 0 {
                            let n = comp
                                .codec(s::HC)
                                .ok_or_else(|| bad("no HC codec"))?
                                .int(st)?;
                            if n < 0 {
                                return Err(beyond());
                            }
                            cig_op = C_HARD;
                            cig_len = cig_len.wrapping_add(n as u32);
                        }
                    }
                    b'P' => {
                        flush_if!(C_PAD);
                        if d & ds::PD != 0 {
                            let n = comp
                                .codec(s::PD)
                                .ok_or_else(|| bad("no PD codec"))?
                                .int(st)?;
                            if n < 0 {
                                return Err(beyond());
                            }
                            cig_op = C_PAD;
                            cig_len = cig_len.wrapping_add(n as u32);
                        }
                    }
                    b'N' => {
                        flush_if!(C_SKIP);
                        if d & ds::RS != 0 {
                            let n = comp
                                .codec(s::RS)
                                .ok_or_else(|| bad("no RS codec"))?
                                .int(st)?;
                            if n < 0 {
                                return Err(beyond());
                            }
                            cig_op = C_SKIP;
                            cig_len = cig_len.wrapping_add(n as u32);
                            ref_pos += n as i64;
                        }
                    }
                    other => {
                        return Err(bad(format!("unknown feature code {:?}", other as char)));
                    }
                }
            }
            if d & ds::FC != 0 && d & ds::FN != 0 && len >= seq_pos as i64 {
                let rest = len - seq_pos as i64 + 1;
                if sq.holds_ref && cr.ref_id >= 0 {
                    if ref_pos + rest > ref_len {
                        let rlen = ref_len - ref_pos;
                        if let Some(r) = reference
                            && seq_on
                        {
                            let at = seq_at + seq_pos as usize - 1;
                            if rlen > 0 {
                                if ref_pos + rlen > r.end {
                                    return Err(beyond());
                                }
                                if (seq_pos as i64 - 1 + rlen) < len {
                                    out.seqs[at..at + rlen as usize]
                                        .copy_from_slice(r.range(ref_pos, rlen));
                                }
                                if rest - rlen > 0 {
                                    out.seqs[at + rlen as usize..seq_at + len as usize].fill(b'N');
                                }
                            } else if rest > 0 {
                                out.seqs[at..seq_at + len as usize].fill(b'N');
                            }
                        }
                        if md_dist >= 0 {
                            md_dist += rest as i32;
                        }
                    } else {
                        if let Some(r) = reference
                            && rest > 0
                        {
                            if ref_pos + rest > r.end {
                                return Err(beyond());
                            }
                            let refp = r.range(ref_pos, rest);
                            if decode_md || decode_nm {
                                match memchr::memchr(b'N', refp) {
                                    None => md_dist += rest as i32,
                                    Some(first) => {
                                        md_dist += first as i32;
                                        for &base in &refp[first..] {
                                            if base == b'N' {
                                                add_md(&mut md, b'N', &mut md_dist);
                                                nm += 1;
                                            } else {
                                                md_dist += 1;
                                            }
                                        }
                                    }
                                }
                            }
                            if seq_on {
                                let at = seq_at + seq_pos as usize - 1;
                                out.seqs[at..at + rest as usize].copy_from_slice(refp);
                            }
                        }
                        ref_pos += rest;
                    }
                } else if cr.ref_id >= 0 {
                    ref_pos += rest;
                }
                flush_if!(C_MATCH);
                cig_op = C_MATCH;
                cig_len = cig_len.wrapping_add(rest as u32);
            }
        }
        if d & ds::FN != 0 && decode_md && md_dist >= 0 {
            push_uint(&mut md, md_dist as u32);
        }
        if cig_len != 0 {
            cigar.push(cig_len << 4 | cig_op);
        }
        cr.ncigar = cigar.len() as u32 - cr.cigar;
        cr.aend = if ref_pos > cr.apos { ref_pos } else { cr.apos };

        cr.mqual = if d & ds::MQ != 0 {
            comp.codec(s::MQ)
                .ok_or_else(|| bad("no MQ codec"))?
                .int(st)?
        } else {
            40
        };
        if d & ds::QS != 0 && cf & PRESERVE_QUAL != 0 {
            let c = comp.codec(s::QS).ok_or_else(|| bad("no QS codec"))?;
            let o = qual_on.then(|| &mut out.quals[qual_at..qual_at + len as usize]);
            c.bytes(st, len as usize, o)?;
        }
        if cf & NO_SEQ != 0 {
            cr.len = 0;
        }

        if decode_md {
            md.push(0);
            let sz = md.len();
            if has_md < 0 {
                let at = (-has_md) as usize;
                out.aux.splice(at..at, md.iter().copied());
                if -has_nm > -has_md {
                    has_nm -= sz as i64;
                }
            } else {
                out.aux.extend_from_slice(b"MDZ");
                out.aux.extend_from_slice(&md);
                cr.aux_size += 3;
            }
            cr.aux_size += sz as u32;
        }
        if decode_nm {
            if has_nm == 0 {
                out.aux.extend_from_slice(b"NM");
                if nm <= u8::MAX as u32 {
                    out.aux.push(b'C');
                    out.aux.push(nm as u8);
                    cr.aux_size += 4;
                } else if nm <= u16::MAX as u32 {
                    out.aux.push(b'S');
                    out.aux.extend_from_slice(&(nm as u16).to_le_bytes());
                    cr.aux_size += 5;
                } else {
                    out.aux.push(b'I');
                    out.aux.extend_from_slice(&nm.to_le_bytes());
                    cr.aux_size += 7;
                }
            } else {
                let at = (-has_nm) as usize;
                if at + 4 <= out.aux.len() {
                    out.aux[at..at + 4].copy_from_slice(&nm.to_le_bytes());
                }
            }
        }
        Ok(())
    }
}

/// `cram_decode_slice_xref`: mate fields and template lengths of mates
/// within the slice.
fn xref(recs: &mut [Rec]) -> Result<()> {
    let n = recs.len();
    for rec in 0..n {
        if recs[rec].mate_line >= 0 {
            let ml = recs[rec].mate_line as usize;
            if ml < n {
                if recs[rec].tlen == i64::MIN {
                    let (id1, mut id2) = (rec, rec);
                    let (mut aleft, mut aright) = (recs[rec].apos, recs[rec].aend);
                    let mut reference = recs[rec].ref_id;
                    let (mut left_cnt, mut right_cnt) = (0, 0);
                    loop {
                        let r2 = &recs[id2];
                        if aleft > r2.apos {
                            aleft = r2.apos;
                            left_cnt = 1;
                        } else if aleft == r2.apos {
                            left_cnt += 1;
                        }
                        if aright < r2.aend {
                            aright = r2.aend;
                            right_cnt = 1;
                        } else if aright == r2.aend {
                            right_cnt += 1;
                        }
                        if r2.mate_line == -1 {
                            recs[id2].mate_line = rec as i32;
                            break;
                        }
                        if (r2.mate_line as i64) <= id2 as i64 || r2.mate_line as usize >= n {
                            return Err(bad("mate line out of order"));
                        }
                        id2 = r2.mate_line as usize;
                        if recs[id2].ref_id != reference {
                            reference = -1;
                        }
                        if id2 == id1 {
                            break;
                        }
                    }
                    if reference != -1 {
                        let mut tlen = aright - aleft + 1;
                        let r = &recs[rec];
                        if r.apos == aleft && (r.aend < aright || left_cnt <= 1) {
                            recs[rec].tlen = tlen;
                            tlen = -tlen;
                        } else if r.apos == aleft
                            && r.aend == aright
                            && left_cnt > 1
                            && right_cnt > 1
                        {
                            if r.flags & READ1 != 0 {
                                recs[rec].tlen = tlen;
                                tlen = -tlen;
                            } else {
                                recs[rec].tlen = -tlen;
                            }
                        } else {
                            recs[rec].tlen = -tlen;
                        }
                        let mut id2 = recs[rec].mate_line as usize;
                        while id2 != rec {
                            recs[id2].tlen = tlen;
                            id2 = recs[id2].mate_line as usize;
                        }
                    } else {
                        recs[rec].tlen = 0;
                        let mut id2 = recs[rec].mate_line as usize;
                        while id2 != rec {
                            recs[id2].tlen = 0;
                            id2 = recs[id2].mate_line as usize;
                        }
                    }
                }
                let ml = recs[rec].mate_line as usize;
                let (mpos, mref, mflags) = (recs[ml].apos, recs[ml].ref_id, recs[ml].flags);
                let cr = &mut recs[rec];
                cr.mate_pos = mpos;
                cr.mate_ref_id = mref;
                cr.flags |= PAIRED;
                if mflags & UNMAPPED != 0 {
                    cr.flags |= MATE_UNMAPPED;
                    cr.tlen = 0;
                }
                if cr.flags & UNMAPPED != 0 {
                    cr.tlen = 0;
                }
                if mflags & REVERSE != 0 {
                    cr.flags |= MATE_REVERSE;
                }
            }
        } else {
            let cr = &mut recs[rec];
            if cr.mate_flags & M_REVERSE != 0 {
                cr.flags |= PAIRED | MATE_REVERSE;
            }
            if cr.mate_flags & M_UNMAP != 0 {
                cr.flags |= MATE_UNMAPPED;
            }
            if cr.flags & PAIRED == 0 {
                cr.mate_ref_id = -1;
            }
        }
        if recs[rec].tlen == i64::MIN {
            recs[rec].tlen = 0;
        }
    }
    for cr in recs.iter_mut() {
        if let Some(t) = cr.explicit_tlen {
            cr.tlen = t;
        }
    }
    Ok(())
}
