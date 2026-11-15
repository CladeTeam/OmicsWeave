//! The CRAM 3 data series codecs, decoded as htslib's `cram_codecs.c`
//! decodes them: EXTERNAL, HUFFMAN, BETA, SUBEXP, GAMMA, BYTE_ARRAY_LEN and
//! BYTE_ARRAY_STOP (GOLOMB and GOLOMB_RICE have no decoder there either).
//! Bit codecs read the slice's CORE block, most significant bit first;
//! the others read external blocks, which are decompressed when first
//! read.

use std::borrow::Cow;
use std::collections::HashMap;

use omicsweave_core::{Error, Result};

use crate::container::{Block, content, method};
use crate::varint::itf8;

/// How a data series' values are serialized (`enum cram_external_type`).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) enum ValueType {
    Int,
    Byte,
    ByteArray,
    /// Byte arrays appended to a growing buffer (read names, tags).
    ByteArrayBlock,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) struct HuffmanCode {
    symbol: i64,
    len: i32,
    code: i32,
    /// `code - index` of the first code of this length.
    p: i32,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub(crate) enum Codec {
    External {
        id: i32,
        ty: ValueType,
    },
    Huffman {
        codes: Box<[HuffmanCode]>,
        ty: ValueType,
    },
    Beta {
        offset: i32,
        nbits: i32,
        ty: ValueType,
    },
    Subexp {
        offset: i32,
        k: i32,
    },
    Gamma {
        offset: i32,
    },
    ByteArrayLen {
        len: Box<Codec>,
        val: Box<Codec>,
    },
    ByteArrayStop {
        stop: u8,
        id: i32,
        ty: ValueType,
    },
}

fn malformed(codec: &'static str) -> Error {
    Error::invalid("CRAM codec", format!("malformed {codec} parameters"))
}

/// Encoding ids.
mod encoding {
    pub const EXTERNAL: i32 = 1;
    pub const HUFFMAN: i32 = 3;
    pub const BYTE_ARRAY_LEN: i32 = 4;
    pub const BYTE_ARRAY_STOP: i32 = 5;
    pub const BETA: i32 = 6;
    pub const SUBEXP: i32 = 7;
    pub const GAMMA: i32 = 9;
}

impl Codec {
    /// The codec of `encoding` with parameters `data`, for values of type
    /// `ty` (`cram_decoder_init`).
    pub(crate) fn parse(encoding: i32, data: &[u8], ty: ValueType) -> Result<Codec> {
        let mut cp = data;
        let done = |cp: &[u8]| cp.is_empty();
        match encoding {
            encoding::EXTERNAL => {
                if data.is_empty() {
                    return Err(malformed("EXTERNAL"));
                }
                let id = itf8(&mut cp, "CRAM codec").map_err(|_| malformed("EXTERNAL"))?;
                if !done(cp) {
                    return Err(malformed("EXTERNAL"));
                }
                Ok(Codec::External { id, ty })
            }
            encoding::HUFFMAN => huffman(data, ty),
            encoding::BYTE_ARRAY_LEN => {
                let mut sub = |ty: ValueType| -> Result<Codec> {
                    let enc = itf8(&mut cp, "CRAM codec").unwrap_or(0);
                    let size = itf8(&mut cp, "CRAM codec").unwrap_or(0);
                    if size < 0 || cp.len() < size as usize {
                        return Err(malformed("BYTE_ARRAY_LEN"));
                    }
                    let (params, rest) = cp.split_at(size as usize);
                    let codec = Codec::parse(enc, params, ty)?;
                    cp = rest;
                    Ok(codec)
                };
                let len = Box::new(sub(ValueType::Int)?);
                let val = Box::new(sub(ty)?);
                if !done(cp) {
                    return Err(malformed("BYTE_ARRAY_LEN"));
                }
                Ok(Codec::ByteArrayLen { len, val })
            }
            encoding::BYTE_ARRAY_STOP => {
                if data.len() < 2 {
                    return Err(malformed("BYTE_ARRAY_STOP"));
                }
                if !matches!(ty, ValueType::ByteArray | ValueType::ByteArrayBlock) {
                    return Err(Error::invalid(
                        "CRAM codec",
                        "BYTE_ARRAY_STOP only codes byte arrays",
                    ));
                }
                let stop = cp[0];
                cp = &cp[1..];
                let id = itf8(&mut cp, "CRAM codec").map_err(|_| malformed("BYTE_ARRAY_STOP"))?;
                if !done(cp) {
                    return Err(malformed("BYTE_ARRAY_STOP"));
                }
                Ok(Codec::ByteArrayStop { stop, id, ty })
            }
            encoding::BETA => {
                if ty == ValueType::ByteArrayBlock {
                    return Err(Error::invalid(
                        "CRAM codec",
                        "BETA does not code byte arrays",
                    ));
                }
                let offset = itf8(&mut cp, "CRAM codec").unwrap_or(0);
                let mut nbits = -1;
                if !cp.is_empty() {
                    nbits = itf8(&mut cp, "CRAM codec").unwrap_or(-1);
                }
                if !done(cp) || !(0..=32).contains(&nbits) {
                    return Err(malformed("BETA"));
                }
                Ok(Codec::Beta { offset, nbits, ty })
            }
            encoding::SUBEXP => {
                if ty != ValueType::Int {
                    return Err(Error::invalid("CRAM codec", "SUBEXP only codes integers"));
                }
                let offset = itf8(&mut cp, "CRAM codec").unwrap_or(0);
                let k = itf8(&mut cp, "CRAM codec").unwrap_or(-1);
                if !done(cp) || !(0..=31).contains(&k) {
                    return Err(malformed("SUBEXP"));
                }
                Ok(Codec::Subexp { offset, k })
            }
            encoding::GAMMA => {
                if ty != ValueType::Int {
                    return Err(Error::invalid("CRAM codec", "GAMMA only codes integers"));
                }
                if data.is_empty() {
                    return Err(malformed("GAMMA"));
                }
                let offset = itf8(&mut cp, "CRAM codec").map_err(|_| malformed("GAMMA"))?;
                if !done(cp) {
                    return Err(malformed("GAMMA"));
                }
                Ok(Codec::Gamma { offset })
            }
            _ => Err(Error::invalid(
                "CRAM codec",
                format!("unimplemented encoding {encoding}"),
            )),
        }
    }

    /// The blocks the codec reads (`cram_codec_to_id`): an external
    /// content id, -1 for the CORE block, -2 for none; the second is the
    /// value codec of BYTE_ARRAY_LEN.
    pub(crate) fn blocks(&self) -> (i32, i32) {
        match self {
            Codec::Huffman { codes, .. } => (if codes.len() == 1 { -2 } else { -1 }, -2),
            Codec::Beta { .. } | Codec::Subexp { .. } | Codec::Gamma { .. } => (-1, -2),
            Codec::External { id, .. } => (*id, -2),
            Codec::ByteArrayLen { len, val } => (len.blocks().0, val.blocks().0),
            Codec::ByteArrayStop { id, .. } => (*id, -2),
        }
    }
}

/// `cram_huffman_decode_init`: canonical codes ordered by length, then
/// symbol.
fn huffman(data: &[u8], ty: ValueType) -> Result<Codec> {
    if ty == ValueType::ByteArrayBlock {
        return Err(Error::invalid(
            "CRAM codec",
            "HUFFMAN does not code byte arrays",
        ));
    }
    let bad = || malformed("HUFFMAN");
    let mut cp = data;
    let n = itf8(&mut cp, "CRAM codec").map_err(|_| bad())?;
    if !(0..1_000_000).contains(&n) {
        return Err(bad());
    }
    let mut codes = Vec::with_capacity(n as usize);
    for _ in 0..n {
        let symbol = match ty {
            ValueType::Int | ValueType::Byte => {
                itf8(&mut cp, "CRAM codec").map_err(|_| bad())? as i64
            }
            _ => return Err(bad()),
        };
        codes.push(HuffmanCode {
            symbol,
            len: 0,
            code: 0,
            p: 0,
        });
    }
    if itf8(&mut cp, "CRAM codec").map_err(|_| bad())? != n {
        return Err(bad());
    }
    if n == 0 {
        return Ok(Codec::Huffman {
            codes: Box::new([]),
            ty,
        });
    }
    let mut max_len = 0;
    for c in &mut codes {
        c.len = itf8(&mut cp, "CRAM codec").map_err(|_| bad())?;
        if c.len < 0 {
            return Err(bad());
        }
        max_len = max_len.max(c.len);
    }
    if !cp.is_empty() || max_len >= n || max_len > 31 {
        return Err(bad());
    }
    codes.sort_by(|a, b| a.len.cmp(&b.len).then(a.symbol.cmp(&b.symbol)));
    let (mut val, mut last_len, mut max_val) = (-1i64, 0i32, 0u64);
    for c in &mut codes {
        val += 1;
        if val as u64 > max_val {
            return Err(bad());
        }
        if c.len > last_len {
            val <<= c.len - last_len;
            last_len = c.len;
            max_val = (1u64 << c.len) - 1;
        }
        c.code = val as i32;
    }
    let (mut last_len, mut j) = (0, 0);
    for (i, c) in codes.iter_mut().enumerate() {
        if c.len > last_len {
            j = c.code - i as i32;
            last_len = c.len;
        }
        c.p = j;
    }
    Ok(Codec::Huffman {
        codes: codes.into_boxed_slice(),
        ty,
    })
}

/// The CORE block, read a bit at a time, most significant bit first.
#[derive(Debug, Default)]
pub(crate) struct Bits<'a> {
    data: &'a [u8],
    byte: usize,
    /// Next bit of the current byte, 7 (the MSB) to 0.
    bit: i32,
}

impl<'a> Bits<'a> {
    pub(crate) fn new(data: &'a [u8]) -> Bits<'a> {
        Bits {
            data,
            byte: 0,
            bit: 7,
        }
    }

    /// `cram_not_enough_bits`.
    #[inline]
    fn short(&self, nbits: u64) -> bool {
        let left = self.data.len() as i64 - self.byte as i64;
        (left <= 0 && nbits > 0)
            || (left <= (i32::MAX / 8 + 1) as i64 && left * 8 + self.bit as i64 - 7 < nbits as i64)
    }

    #[inline]
    fn bit(&mut self) -> u32 {
        let v = (self.data[self.byte] >> self.bit) as u32 & 1;
        self.bit -= 1;
        if self.bit < 0 {
            self.bit = 7;
            self.byte += 1;
        }
        v
    }

    fn bits(&mut self, n: i32) -> u64 {
        let mut v = 0u64;
        for _ in 0..n {
            v = v << 1 | self.bit() as u64;
        }
        v
    }

    /// `get_one_bits_MSB`: the number of 1 bits before a 0.
    fn ones(&mut self) -> Option<i32> {
        if self.byte >= self.data.len() {
            return None;
        }
        let mut n = 0;
        loop {
            let b = self.data[self.byte] >> self.bit;
            self.bit -= 1;
            if self.bit == -1 {
                self.bit = 7;
                self.byte += 1;
                if self.byte == self.data.len() && b & 1 != 0 {
                    return None;
                }
            }
            n += 1;
            if b & 1 == 0 {
                return Some(n - 1);
            }
        }
    }

    /// `get_zero_bits_MSB`: the number of 0 bits before a 1.
    fn zeros(&mut self) -> Option<i32> {
        if self.byte >= self.data.len() {
            return None;
        }
        let mut n = 0;
        loop {
            let b = self.data[self.byte] >> self.bit;
            self.bit -= 1;
            if self.bit == -1 {
                self.bit = 7;
                self.byte += 1;
                if self.byte == self.data.len() && b & 1 == 0 {
                    return None;
                }
            }
            n += 1;
            if b & 1 != 0 {
                return Some(n - 1);
            }
        }
    }
}

/// An external block and the read position of its decoders.
struct External<'a> {
    block: Block<'a>,
    data: Option<Cow<'a, [u8]>>,
    idx: usize,
}

/// The blocks of a slice as its decoders read them: the CORE bit stream
/// and every external block by content id.
pub(crate) struct Streams<'a> {
    pub(crate) core: Bits<'a>,
    ext: Vec<External<'a>>,
    /// Content ids below 256: the last external block of that id.
    small: [u16; 256],
    /// Larger ids by `256 + id % 251`, as htslib's table keys them.
    large: HashMap<u32, usize>,
    crc: bool,
}

fn exhausted() -> Error {
    Error::invalid(
        "CRAM slice",
        "a data series reads past the end of its block",
    )
}

fn missing(id: i32) -> Error {
    Error::invalid("CRAM slice", format!("no external block {id}"))
}

impl<'a> Streams<'a> {
    /// The streams of a slice's blocks: `core` uncompressed, the external
    /// ones decompressed on first use (checking their CRC32 when `crc`).
    pub(crate) fn new(core: &'a [u8], blocks: &[Block<'a>], crc: bool) -> Streams<'a> {
        let mut s = Streams {
            core: Bits::new(core),
            ext: Vec::new(),
            small: [0; 256],
            large: HashMap::new(),
            crc,
        };
        for b in blocks
            .iter()
            .filter(|b| b.content_type == content::EXTERNAL)
        {
            let at = s.ext.len();
            s.ext.push(External {
                block: *b,
                data: None,
                idx: 0,
            });
            let v = b.content_id as u32;
            if v < 256 {
                s.small[v as usize] = at as u16 + 1;
            } else {
                s.large.insert(256 + v % 251, at);
            }
        }
        s
    }

    fn find(&self, id: i32) -> Option<usize> {
        let v = id as u32;
        if v < 256 {
            return (self.small[v as usize] as usize).checked_sub(1);
        }
        if let Some(&at) = self.large.get(&(256 + v % 251))
            && self.ext[at].block.content_id == id
        {
            return Some(at);
        }
        self.ext.iter().position(|e| e.block.content_id == id)
    }

    /// The external block `id`, decompressed; `None` when the slice has
    /// none.
    fn external(&mut self, id: i32) -> Result<Option<(&[u8], &mut usize, bool)>> {
        let Some(at) = self.find(id) else {
            return Ok(None);
        };
        let crc = self.crc;
        let e = &mut self.ext[at];
        if e.data.is_none() {
            if crc {
                e.block.check_crc()?;
            }
            let b = &e.block;
            e.data = Some(if b.uncomp_size == 0 {
                Cow::Borrowed(&[][..])
            } else if b.method == method::RAW {
                Cow::Borrowed(b.data)
            } else {
                let mut out = Vec::new();
                crate::container::decompress(b.method, b.data, b.uncomp_size, &mut out)?;
                Cow::Owned(out)
            });
        }
        let tok3 = e.block.method == method::TOK3;
        Ok(Some((e.data.as_deref().unwrap(), &mut e.idx, tok3)))
    }

    /// The uncompressed bytes of external block `id`, whole (embedded
    /// references).
    pub(crate) fn whole(&mut self, id: i32) -> Result<Option<&[u8]>> {
        Ok(self.external(id)?.map(|(d, _, _)| d))
    }

    /// `cram_extract_block`: the next `n` bytes of block `id`.
    fn extract(&mut self, id: i32, n: usize) -> Result<&[u8]> {
        let (data, idx, _) = self.external(id)?.ok_or_else(|| missing(id))?;
        let at = *idx;
        *idx += n;
        if *idx > data.len() {
            return Err(exhausted());
        }
        Ok(&data[at..at + n])
    }
}

impl Codec {
    fn huffman_symbol(codes: &[HuffmanCode], bits: &mut Bits<'_>) -> Result<i64> {
        let (mut idx, mut val, mut len, mut last_len) = (0usize, 0i32, 0i32, 0i32);
        loop {
            let dlen = codes[idx].len - last_len;
            if bits.short(dlen.max(0) as u64) {
                return Err(exhausted());
            }
            len += dlen;
            last_len = len;
            for _ in 0..dlen {
                val = val << 1 | bits.bit() as i32;
            }
            let i = val as i64 - codes[idx].p as i64;
            if i < 0 || i >= codes.len() as i64 {
                return Err(Error::invalid("CRAM slice", "invalid HUFFMAN code"));
            }
            idx = i as usize;
            if codes[idx].code == val && codes[idx].len == len {
                return Ok(codes[idx].symbol);
            }
        }
    }

    /// One integer of an integer data series.
    pub(crate) fn int(&self, s: &mut Streams<'_>) -> Result<i32> {
        match self {
            Codec::External {
                id,
                ty: ValueType::Int,
            } => {
                let (data, idx, _) = s.external(*id)?.ok_or_else(|| missing(*id))?;
                let mut cp = &data[*idx..];
                let v = itf8(&mut cp, "CRAM external integer")?;
                *idx = data.len() - cp.len();
                Ok(v)
            }
            Codec::Huffman { codes, ty } if *ty != ValueType::Byte => {
                if codes.is_empty() {
                    return Err(Error::invalid("CRAM slice", "a data series has no codes"));
                }
                if codes[0].len == 0 {
                    return Ok(codes[0].symbol as i32);
                }
                Ok(Self::huffman_symbol(codes, &mut s.core)? as i32)
            }
            Codec::Beta { offset, nbits, ty } if *ty != ValueType::Byte => {
                if *nbits == 0 {
                    return Ok(offset.wrapping_neg());
                }
                if s.core.short(*nbits as u64) {
                    return Err(exhausted());
                }
                Ok((s.core.bits(*nbits) as i64 - *offset as i64) as i32)
            }
            Codec::Subexp { offset, k } => {
                let k = *k;
                let u = s.core.ones().ok_or_else(exhausted)?;
                if s.core
                    .short(if u > 0 { (u + k - 1) as u64 } else { k as u64 })
                {
                    return Err(exhausted());
                }
                let val = if u > 0 {
                    if u > 31 - k {
                        return Err(Error::invalid("CRAM slice", "SUBEXP value too large"));
                    }
                    let tail = u + k - 1;
                    (s.core.bits(tail) as i32).wrapping_add(1 << tail)
                } else {
                    s.core.bits(k) as i32
                };
                Ok(val.wrapping_sub(*offset))
            }
            Codec::Gamma { offset } => {
                let nz = s.core.zeros().ok_or_else(exhausted)?;
                if s.core.short(nz as u64) {
                    return Err(exhausted());
                }
                let mut val: i32 = 1;
                for _ in 0..nz {
                    val = val.wrapping_shl(1) | s.core.bit() as i32;
                }
                Ok(val.wrapping_sub(*offset))
            }
            _ => Err(Error::invalid(
                "CRAM slice",
                "an integer data series has a byte codec",
            )),
        }
    }

    /// `n` values of a byte data series (or of the values of a byte
    /// array), into `out` when given.
    pub(crate) fn bytes(
        &self,
        s: &mut Streams<'_>,
        n: usize,
        out: Option<&mut [u8]>,
    ) -> Result<()> {
        match self {
            Codec::External { id, .. } => {
                let got = s.extract(*id, n)?;
                if let Some(out) = out {
                    out[..n].copy_from_slice(got);
                }
                Ok(())
            }
            Codec::Huffman { codes, .. } => {
                if codes.is_empty() {
                    return Err(Error::invalid("CRAM slice", "a data series has no codes"));
                }
                if codes[0].len == 0 {
                    if let Some(out) = out {
                        out[..n].fill(codes[0].symbol as u8);
                    }
                    return Ok(());
                }
                let mut out = out;
                for i in 0..n {
                    let v = Self::huffman_symbol(codes, &mut s.core)? as u8;
                    if let Some(o) = out.as_deref_mut() {
                        o[i] = v;
                    }
                }
                Ok(())
            }
            Codec::Beta { offset, nbits, .. } => {
                if *nbits == 0 {
                    if let Some(out) = out {
                        out[..n].fill(offset.wrapping_neg() as u8);
                    }
                    return Ok(());
                }
                if s.core.short(*nbits as u64 * n as u64) {
                    return Err(exhausted());
                }
                let mut out = out;
                for i in 0..n {
                    let v = (s.core.bits(*nbits) as i64 - *offset as i64) as u8;
                    if let Some(o) = out.as_deref_mut() {
                        o[i] = v;
                    }
                }
                Ok(())
            }
            _ => Err(Error::invalid(
                "CRAM slice",
                "a byte data series has an integer or array codec",
            )),
        }
    }

    /// One byte array of a data series of type BYTE_ARRAY (IN, SC, BB,
    /// QQ), at most `room` bytes when `out` is given; returns its length.
    pub(crate) fn array(
        &self,
        s: &mut Streams<'_>,
        room: i32,
        out: Option<&mut [u8]>,
    ) -> Result<i32> {
        match self {
            Codec::ByteArrayLen { len, val } => {
                let n = len.int(s)?;
                let unbounded = matches!(
                    **val,
                    Codec::External {
                        ty: ValueType::ByteArrayBlock,
                        ..
                    }
                );
                if n < 0
                    || (n > room && !unbounded)
                    || out.as_ref().is_some_and(|o| n as usize > o.len())
                {
                    return Err(Error::invalid(
                        "CRAM slice",
                        "byte array overflows its read",
                    ));
                }
                val.bytes(s, n as usize, out)?;
                Ok(n)
            }
            Codec::ByteArrayStop { stop, id, .. } => {
                let stop = *stop;
                let (data, idx, _) = s.external(*id)?.ok_or_else(|| missing(*id))?;
                if *idx >= data.len() {
                    return Err(exhausted());
                }
                let rest = &data[*idx..];
                let found = match &out {
                    // At most `room` bytes, then the stop byte.
                    Some(_) => rest[..rest.len().min(room.max(0) as usize + 1)]
                        .iter()
                        .position(|&c| c == stop),
                    None => memchr::memchr(stop, rest),
                };
                let Some(n) = found else {
                    return Err(Error::invalid("CRAM slice", "byte array has no stop byte"));
                };
                if let Some(out) = out {
                    out[..n].copy_from_slice(&rest[..n]);
                }
                *idx += n + 1;
                Ok(n as i32)
            }
            _ => {
                let n = room.max(0) as usize;
                self.bytes(s, n, out)?;
                Ok(room)
            }
        }
    }

    /// One byte array of a data series of type BYTE_ARRAY_BLOCK (read
    /// names, tags), appended to `out`; `n` is the length a plain EXTERNAL
    /// codec reads. Returns its length.
    pub(crate) fn block(&self, s: &mut Streams<'_>, n: i32, out: &mut Vec<u8>) -> Result<i32> {
        match self {
            Codec::External { id, .. } => {
                let got = s.extract(*id, n.max(0) as usize)?;
                out.extend_from_slice(got);
                Ok(n)
            }
            Codec::ByteArrayLen { len, val } => {
                let len_ = len.int(s)?;
                let unbounded = matches!(
                    **val,
                    Codec::External {
                        ty: ValueType::ByteArrayBlock,
                        ..
                    }
                );
                if len_ < 0 || (len_ > n && !unbounded) {
                    return Err(Error::invalid("CRAM slice", "byte array overflows"));
                }
                val.block(s, len_, out)?;
                Ok(len_)
            }
            Codec::ByteArrayStop { stop, id, .. } => {
                let (data, idx, tok3) = s.external(*id)?.ok_or_else(|| missing(*id))?;
                if *idx >= data.len() {
                    return Err(exhausted());
                }
                // The name tokeniser always ends names with NUL.
                let stop = if tok3 { 0 } else { *stop };
                let rest = &data[*idx..];
                let n = memchr::memchr(stop, rest).unwrap_or(rest.len());
                out.extend_from_slice(&rest[..n]);
                *idx += n + 1;
                Ok(n as i32)
            }
            _ => Err(Error::invalid(
                "CRAM slice",
                "a byte-array data series has an unsupported codec",
            )),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn huffman_canonical_codes() {
        // Symbols 1, 2, 3 with lengths 1, 2, 2: codes 0, 10, 11.
        let params = [3, 1, 2, 3, 3, 1, 2, 2];
        let c = Codec::parse(encoding::HUFFMAN, &params, ValueType::Int).unwrap();
        let core = [0b0101_1000];
        let mut s = Streams::new(&core, &[], false);
        let got: Vec<i32> = (0..3).map(|_| c.int(&mut s).unwrap()).collect();
        assert_eq!(got, [1, 2, 3]);
    }

    #[test]
    fn single_symbol_huffman_reads_nothing() {
        let c = Codec::parse(encoding::HUFFMAN, &[1, 65, 1, 0], ValueType::Byte).unwrap();
        assert_eq!(c.blocks(), (-2, -2));
        let mut s = Streams::new(&[], &[], false);
        let mut out = [0u8; 3];
        c.bytes(&mut s, 3, Some(&mut out)).unwrap();
        assert_eq!(&out, b"AAA");
    }

    #[test]
    fn beta_and_gamma() {
        let beta = Codec::parse(encoding::BETA, &[0, 3], ValueType::Int).unwrap();
        let gamma = Codec::parse(encoding::GAMMA, &[0], ValueType::Int).unwrap();
        // beta 101 = 5, gamma 00101 = 5
        let core = [0b1010_0101];
        let mut s = Streams::new(&core, &[], false);
        assert_eq!(beta.int(&mut s).unwrap(), 5);
        assert_eq!(gamma.int(&mut s).unwrap(), 5);
    }
}
