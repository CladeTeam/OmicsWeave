//! The compression header of a container and the header of a slice, as
//! htslib's `cram_decode_compression_header` and
//! `cram_decode_slice_header` read them; which data series a decode needs
//! (`cram_dependent_data_series`).

use std::collections::HashMap;

use omicsweave_core::{Error, Result};

use crate::codec::{Codec, ValueType};
use crate::container::content;
use crate::varint::{itf8, ltf8};

/// Data series, by the bit htslib gives each in a slice's set of series to
/// decode (`enum cram_fields`).
pub(crate) mod ds {
    pub const BF: u32 = 1 << 0;
    pub const AP: u32 = 1 << 1;
    pub const FP: u32 = 1 << 2;
    pub const RL: u32 = 1 << 3;
    pub const DL: u32 = 1 << 4;
    pub const NF: u32 = 1 << 5;
    pub const BA: u32 = 1 << 6;
    pub const QS: u32 = 1 << 7;
    pub const FC: u32 = 1 << 8;
    pub const FN: u32 = 1 << 9;
    pub const BS: u32 = 1 << 10;
    pub const IN: u32 = 1 << 11;
    pub const RG: u32 = 1 << 12;
    pub const MQ: u32 = 1 << 13;
    pub const TL: u32 = 1 << 14;
    pub const RN: u32 = 1 << 15;
    pub const NS: u32 = 1 << 16;
    pub const NP: u32 = 1 << 17;
    pub const TS: u32 = 1 << 18;
    pub const MF: u32 = 1 << 19;
    pub const CF: u32 = 1 << 20;
    pub const RI: u32 = 1 << 21;
    pub const RS: u32 = 1 << 22;
    pub const PD: u32 = 1 << 23;
    pub const HC: u32 = 1 << 24;
    pub const SC: u32 = 1 << 25;
    pub const BB: u32 = 1 << 26;
    pub const QQ: u32 = 1 << 28;
    pub const AUX: u32 = 1 << 30;
    pub const ALL: u32 = 0x7fff_ffff;

    pub const CIGAR: u32 = FN | FP | FC | DL | IN | SC | HC | PD | RS | RL | BF;
    pub const SEQ: u32 = CIGAR | BA | BS | RL | AP | BB;
    pub const QUAL: u32 = CIGAR | RL | AP | QS | QQ;
}

/// Index of a data series in [`CompressionHeader::codecs`].
pub(crate) type Series = usize;

/// The two-letter keys of the record encoding map, their data series bit
/// and value type. The index in this table is the series' [`Series`].
pub(crate) const SERIES: [([u8; 2], u32, ValueType); 29] = [
    (*b"BF", ds::BF, ValueType::Int),
    (*b"CF", ds::CF, ValueType::Int),
    (*b"RI", ds::RI, ValueType::Int),
    (*b"RL", ds::RL, ValueType::Int),
    (*b"AP", ds::AP, ValueType::Int),
    (*b"RG", ds::RG, ValueType::Int),
    (*b"MF", ds::MF, ValueType::Int),
    (*b"NS", ds::NS, ValueType::Int),
    (*b"NP", ds::NP, ValueType::Int),
    (*b"TS", ds::TS, ValueType::Int),
    (*b"NF", ds::NF, ValueType::Int),
    (*b"TC", 0, ValueType::Byte),
    (*b"TN", 0, ValueType::Int),
    (*b"FN", ds::FN, ValueType::Int),
    (*b"FC", ds::FC, ValueType::Byte),
    (*b"FP", ds::FP, ValueType::Int),
    (*b"BS", ds::BS, ValueType::Byte),
    (*b"IN", ds::IN, ValueType::ByteArray),
    (*b"SC", ds::SC, ValueType::ByteArray),
    (*b"DL", ds::DL, ValueType::Int),
    (*b"BA", ds::BA, ValueType::Byte),
    (*b"BB", ds::BB, ValueType::ByteArray),
    (*b"RS", ds::RS, ValueType::Int),
    (*b"PD", ds::PD, ValueType::Int),
    (*b"HC", ds::HC, ValueType::Int),
    (*b"MQ", ds::MQ, ValueType::Int),
    (*b"RN", ds::RN, ValueType::ByteArrayBlock),
    (*b"QS", ds::QS, ValueType::Byte),
    (*b"QQ", ds::QQ, ValueType::ByteArray),
];

/// Series indexes, in [`SERIES`] order.
pub(crate) mod s {
    use super::Series;
    pub const BF: Series = 0;
    pub const CF: Series = 1;
    pub const RI: Series = 2;
    pub const RL: Series = 3;
    pub const AP: Series = 4;
    pub const RG: Series = 5;
    pub const MF: Series = 6;
    pub const NS: Series = 7;
    pub const NP: Series = 8;
    pub const TS: Series = 9;
    pub const NF: Series = 10;
    pub const FN: Series = 13;
    pub const FC: Series = 14;
    pub const FP: Series = 15;
    pub const BS: Series = 16;
    pub const IN: Series = 17;
    pub const SC: Series = 18;
    pub const DL: Series = 19;
    pub const BA: Series = 20;
    pub const BB: Series = 21;
    pub const RS: Series = 22;
    pub const PD: Series = 23;
    pub const HC: Series = 24;
    pub const MQ: Series = 25;
    pub const RN: Series = 26;
    pub const QS: Series = 27;
    pub const QQ: Series = 28;
    pub const TL: Series = 29;
}

/// Slots in [`CompressionHeader::codecs`]: [`SERIES`] and TL.
const SLOTS: usize = SERIES.len() + 1;

/// A container's compression header.
#[derive(Debug, Clone)]
pub struct CompressionHeader {
    pub(crate) read_names_included: bool,
    pub(crate) ap_delta: bool,
    pub(crate) no_ref: bool,
    pub(crate) qs_seq_orient: bool,
    /// Base for (reference base A, C, G, T, N; substitution code 0-3).
    pub(crate) sub_matrix: [[u8; 4]; 5],
    /// Tag lines (TD): each a list of `tag + type` keys.
    pub(crate) tag_lines: Vec<Vec<[u8; 3]>>,
    /// Whether any tag line holds a `CG` tag.
    pub(crate) has_cg: bool,
    pub(crate) codecs: Vec<Option<Codec>>,
    /// Tag codecs by `tag[0] << 16 | tag[1] << 8 | type`.
    pub(crate) tags: HashMap<u32, Codec>,
}

fn invalid(detail: impl Into<String>) -> Error {
    Error::invalid("CRAM compression header", detail)
}

/// The tag dictionary (`cram_decode_TD`): NUL-separated lines of 3-byte
/// keys.
fn tag_lines(data: &[u8]) -> Vec<Vec<[u8; 3]>> {
    let mut data = data.to_vec();
    if data.last().is_some_and(|&c| c != 0) {
        data.push(0);
    }
    let mut lines = Vec::new();
    let mut i = 0;
    while i < data.len() {
        let end = i + data[i..]
            .iter()
            .position(|&c| c == 0)
            .unwrap_or(data.len() - i);
        lines.push(
            data[i..end]
                .chunks_exact(3)
                .map(|k| [k[0], k[1], k[2]])
                .collect(),
        );
        i = end + 1;
    }
    lines
}

impl CompressionHeader {
    /// Parses the uncompressed bytes of a compression header block.
    pub(crate) fn parse(data: &[u8]) -> Result<CompressionHeader> {
        let what = "CRAM compression header";
        let mut h = CompressionHeader {
            read_names_included: false,
            ap_delta: true,
            no_ref: false,
            qs_seq_orient: true,
            sub_matrix: [*b"CGTN", *b"AGTN", *b"ACTN", *b"ACGN", *b"ACGT"],
            tag_lines: Vec::new(),
            has_cg: false,
            codecs: vec![None; SLOTS],
            tags: HashMap::new(),
        };
        let mut cp = data;

        // Preservation map.
        let size = itf8(&mut cp, what)?;
        let start = cp.len();
        let count = itf8(&mut cp, what)?;
        for _ in 0..count {
            if cp.len() < 3 {
                return Err(Error::Truncated {
                    what,
                    need: 3,
                    have: cp.len() as u64,
                });
            }
            let key = [cp[0], cp[1]];
            cp = &cp[2..];
            match &key {
                b"MI" | b"UI" | b"PI" => cp = &cp[1..],
                b"RN" => {
                    h.read_names_included = cp[0] != 0;
                    cp = &cp[1..];
                }
                b"AP" => {
                    h.ap_delta = cp[0] != 0;
                    cp = &cp[1..];
                }
                b"RR" => {
                    h.no_ref = cp[0] == 0;
                    cp = &cp[1..];
                }
                b"QO" => {
                    h.qs_seq_orient = cp[0] != 0;
                    cp = &cp[1..];
                }
                b"SM" => {
                    if cp.len() < 5 {
                        return Err(invalid("substitution matrix cut short"));
                    }
                    const ROWS: [&[u8; 4]; 5] = [b"CGTN", b"AGTN", b"ACTN", b"ACGN", b"ACGT"];
                    for (r, row) in ROWS.iter().enumerate() {
                        let b = cp[r];
                        for (k, &base) in row.iter().enumerate() {
                            h.sub_matrix[r][((b >> (6 - 2 * k)) & 3) as usize] = base;
                        }
                    }
                    cp = &cp[5..];
                }
                b"TD" => {
                    let n = itf8(&mut cp, what)?;
                    if n < 0 || cp.len() < n as usize {
                        return Err(invalid("tag dictionary cut short"));
                    }
                    h.tag_lines = tag_lines(&cp[..n as usize]);
                    cp = &cp[n as usize..];
                }
                _ => cp = &cp[1..],
            }
        }
        if start - cp.len() != size as usize {
            return Err(invalid("preservation map size"));
        }
        h.has_cg = h
            .tag_lines
            .iter()
            .any(|l| l.iter().any(|k| k[0] == b'C' && k[1] == b'G'));

        // Data series encodings.
        let size = itf8(&mut cp, what)?;
        let start = cp.len();
        let count = itf8(&mut cp, what)?;
        for _ in 0..count {
            if cp.len() < 4 {
                return Err(invalid("data series encoding cut short"));
            }
            let key = [cp[0], cp[1]];
            cp = &cp[2..];
            let encoding = itf8(&mut cp, what)?;
            let len = itf8(&mut cp, what)?;
            if encoding == 0 {
                continue;
            }
            if len < 0 || cp.len() < len as usize {
                return Err(invalid("data series encoding cut short"));
            }
            let params = &cp[..len as usize];
            let slot = if &key == b"TL" {
                Some((s::TL, ValueType::Int))
            } else {
                SERIES
                    .iter()
                    .position(|(k, _, _)| *k == key)
                    .map(|i| (i, SERIES[i].2))
            };
            if let Some((slot, ty)) = slot {
                h.codecs[slot] = Some(Codec::parse(encoding, params, ty)?);
            }
            cp = &cp[len as usize..];
        }
        if start - cp.len() != size as usize {
            return Err(invalid("data series encoding map size"));
        }

        // Tag encodings.
        let size = itf8(&mut cp, what)?;
        let start = cp.len();
        let count = itf8(&mut cp, what)?;
        for _ in 0..count {
            if cp.len() < 6 {
                return Err(invalid("tag encoding cut short"));
            }
            let key = itf8(&mut cp, what)? as u32;
            let encoding = itf8(&mut cp, what)?;
            let len = itf8(&mut cp, what)?;
            if len < 0 || cp.len() < len as usize {
                return Err(invalid("tag encoding cut short"));
            }
            let codec = Codec::parse(encoding, &cp[..len as usize], ValueType::ByteArrayBlock)?;
            h.tags.insert(key & 0xff_ffff, codec);
            cp = &cp[len as usize..];
        }
        if start - cp.len() != size as usize {
            return Err(invalid("tag encoding map size"));
        }
        Ok(h)
    }

    pub(crate) fn codec(&self, series: Series) -> Option<&Codec> {
        self.codecs[series].as_ref()
    }

    /// The data series to decode for `wanted`, closed over what decoding
    /// them needs: series their values depend on, and every series (or tag)
    /// that reads a block one of them reads, since a block's values are
    /// consumed in order.
    pub(crate) fn closure(&self, wanted: u32, blocks: &[(u8, i32)]) -> u32 {
        let mut d = wanted;
        let bit_of = |i: usize| if i == s::TL { ds::TL } else { SERIES[i].1 };
        let ids = |c: &Codec| {
            let (a, b) = c.blocks();
            [a, b]
        };
        let present = |id: i32| {
            id == -1
                || blocks
                    .iter()
                    .any(|&(ty, cid)| ty == content::EXTERNAL && cid == id)
        };
        loop {
            let before = d;
            for (dep, needs) in [
                (
                    ds::RS | ds::PD | ds::HC | ds::QS | ds::IN | ds::SC | ds::BS,
                    ds::FC | ds::FP,
                ),
                (ds::DL | ds::BA | ds::BB | ds::QQ, ds::FC | ds::FP),
                (ds::SEQ | ds::CIGAR, ds::RL),
                (ds::FP, ds::FC),
                (ds::FC, ds::FN),
                (ds::AUX, ds::TL),
                (ds::MF, ds::CF),
                (ds::MQ, ds::BF),
                (ds::BS, ds::RI),
                (ds::MF | ds::NS | ds::NP | ds::TS | ds::NF, ds::CF),
                (ds::BA | ds::QS | ds::BB | ds::QQ, ds::BF | ds::CF | ds::RL),
                (ds::FN, ds::SC | ds::IN | ds::BB),
            ] {
                if d & dep != 0 {
                    d |= needs;
                }
            }
            if !self.read_names_included && d & ds::RN != 0 {
                d |= ds::CF | ds::NF;
            }
            let mut used: Vec<i32> = Vec::new();
            for (i, c) in self.codecs.iter().enumerate() {
                if let Some(c) = c
                    && d & bit_of(i) != 0
                {
                    used.extend(ids(c).into_iter().filter(|&id| id != -2 && present(id)));
                }
            }
            if d & ds::AUX != 0 {
                for c in self.tags.values() {
                    used.extend(ids(c).into_iter().filter(|&id| id != -2 && present(id)));
                }
            }
            for (i, c) in self.codecs.iter().enumerate() {
                if let Some(c) = c
                    && ids(c).iter().any(|id| *id != -2 && used.contains(id))
                {
                    d |= bit_of(i);
                }
            }
            if self
                .tags
                .values()
                .any(|c| ids(c).iter().any(|id| *id != -2 && used.contains(id)))
            {
                d |= ds::AUX;
            }
            if d == before {
                return d;
            }
        }
    }
}

/// A slice header.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct SliceHeader {
    /// Reference id; -1 unmapped, -2 several.
    pub ref_seq_id: i32,
    pub ref_seq_start: i64,
    pub ref_seq_span: i64,
    pub num_records: i32,
    pub record_counter: i64,
    pub num_blocks: i32,
    /// External block holding the slice's reference, if embedded.
    pub ref_base_id: i32,
    pub md5: [u8; 16],
}

impl SliceHeader {
    pub(crate) fn parse(content_type: u8, data: &[u8]) -> Result<SliceHeader> {
        let what = "CRAM slice header";
        let mut cp = data;
        let mapped = content_type == content::MAPPED_SLICE;
        if !mapped && content_type != content::UNMAPPED_SLICE {
            return Err(Error::invalid(what, "not a slice header block"));
        }
        let (mut ref_seq_id, mut ref_seq_start, mut ref_seq_span) = (-1, 0, 0);
        if mapped {
            ref_seq_id = itf8(&mut cp, what)?;
            ref_seq_start = itf8(&mut cp, what)? as i64;
            ref_seq_span = itf8(&mut cp, what)? as i64;
            if ref_seq_start < 0 || ref_seq_span < 0 {
                return Err(Error::invalid(what, "negative start or span"));
            }
        }
        let num_records = itf8(&mut cp, what)?;
        let record_counter = ltf8(&mut cp, what)?;
        let num_blocks = itf8(&mut cp, what)?;
        let num_ids = itf8(&mut cp, what)?;
        if !(1..10000).contains(&num_ids) || !(1..10000).contains(&num_blocks) {
            return Err(Error::invalid(what, "block counts out of range"));
        }
        for _ in 0..num_ids {
            itf8(&mut cp, what)?;
        }
        let ref_base_id = if mapped { itf8(&mut cp, what)? } else { -1 };
        if cp.len() < 16 {
            return Err(Error::invalid(what, "reference MD5 cut short"));
        }
        let md5 = cp[..16].try_into().unwrap();
        Ok(SliceHeader {
            ref_seq_id,
            ref_seq_start,
            ref_seq_span,
            num_records,
            record_counter,
            num_blocks,
            ref_base_id,
            md5,
        })
    }

    pub fn ref_seq_end(&self) -> i64 {
        self.ref_seq_start + self.ref_seq_span - 1
    }
}
