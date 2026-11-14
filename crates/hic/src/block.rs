//! Contact blocks: zlib streams of records, decoded into (x bin, y bin,
//! value) triples, and encoded from them.

use flate2::{Compress, Compression, FlushCompress, Status};
use omicsweave_core::{Error, Result};

/// One cell of a map, in bins: `x` on the map's first chromosome, `y` on
/// its second.
///
/// numpy: `np.dtype([("x", "<u4"), ("y", "<u4"), ("value", "<f4")])`.
#[repr(C)]
#[derive(Debug, Clone, Copy, Default, PartialEq)]
pub struct Record {
    pub x: u32,
    pub y: u32,
    pub value: f32,
}

const _: () = assert!(std::mem::size_of::<Record>() == 12);

/// Reusable zlib state; keep one per thread.
pub struct Inflater {
    inner: libdeflater::Decompressor,
}
impl Default for Inflater {
    fn default() -> Self {
        Self::new()
    }
}
impl Inflater {
    pub fn new() -> Self {
        Self {
            inner: libdeflater::Decompressor::new(),
        }
    }

    /// Inflates a whole block, growing initialized caller storage as needed.
    /// Previous output length seeds the next block's capacity estimate.
    pub fn inflate(&mut self, raw: &[u8], out: &mut Vec<u8>) -> Result<()> {
        let mut room = out.len().max(raw.len().saturating_mul(4)).max(4096);
        loop {
            if out.len() < room {
                out.try_reserve(room - out.len())
                    .map_err(|_| Error::Inflate("output allocation failed"))?;
                out.resize(room, 0);
            }
            match self.inner.zlib_decompress(raw, out) {
                Ok(n) => {
                    out.truncate(n);
                    return Ok(());
                }
                Err(libdeflater::DecompressionError::BadData) => {
                    out.clear();
                    return Err(Error::Inflate("corrupt or truncated zlib stream"));
                }
                Err(libdeflater::DecompressionError::InsufficientSpace) => {
                    room = room
                        .checked_mul(2)
                        .ok_or(Error::Inflate("output size overflows usize"))?;
                }
            }
        }
    }
}

/// Reusable zlib compression state.
pub struct Deflater {
    z: Compress,
}

impl Deflater {
    /// `level` 0-9; the Juicer tools use 6 (zlib's default).
    pub fn new(level: u32) -> Deflater {
        Deflater {
            z: Compress::new(Compression::new(level.min(9)), true),
        }
    }

    /// Compresses `data` into `out`, replacing its contents.
    pub fn deflate(&mut self, data: &[u8], out: &mut Vec<u8>) -> Result<()> {
        self.z.reset();
        out.clear();
        out.reserve(data.len() / 2 + 64);
        loop {
            let read = self.z.total_in() as usize;
            let status = self
                .z
                .compress_vec(&data[read..], out, FlushCompress::Finish)
                .map_err(|_| Error::invalid("hic block", "compression failed"))?;
            match status {
                Status::StreamEnd => return Ok(()),
                _ => out.reserve(out.capacity().max(64)),
            }
        }
    }
}

/// Little-endian reads with bounds checks over a decompressed block.
struct Reader<'a> {
    buf: &'a [u8],
    at: usize,
}

impl<'a> Reader<'a> {
    fn records(&mut self, n: usize, width: usize) -> Result<&'a [u8]> {
        let end = n
            .checked_mul(width)
            .and_then(|len| self.at.checked_add(len))
            .ok_or_else(|| Error::invalid("hic block", "record array overflows usize"))?;
        let out = self.buf.get(self.at..end).ok_or(Error::Truncated {
            what: "hic block",
            need: end as u64,
            have: self.buf.len() as u64,
        })?;
        self.at = end;
        Ok(out)
    }

    fn take<const N: usize>(&mut self) -> Result<[u8; N]> {
        let end = self.at + N;
        let bytes = self.buf.get(self.at..end).ok_or(Error::Truncated {
            what: "hic block",
            need: end as u64,
            have: self.buf.len() as u64,
        })?;
        self.at = end;
        Ok(bytes.try_into().unwrap())
    }

    fn u8(&mut self) -> Result<u8> {
        Ok(self.take::<1>()?[0])
    }

    fn i16(&mut self) -> Result<i16> {
        Ok(i16::from_le_bytes(self.take()?))
    }

    fn i32(&mut self) -> Result<i32> {
        Ok(i32::from_le_bytes(self.take()?))
    }

    /// A value stored as `int` when `wide`, else as `short`.
    fn int(&mut self, wide: bool) -> Result<i32> {
        if wide {
            self.i32()
        } else {
            Ok(self.i16()? as i32)
        }
    }

    fn count(&mut self, wide: bool) -> Result<usize> {
        let n = self.int(wide)?;
        usize::try_from(n).map_err(|_| Error::invalid("hic block", format!("negative count {n}")))
    }
}

/// # Safety
/// ptr points to at least N readable bytes, with no alignment requirement.
#[inline]
unsafe fn field<const N: usize>(ptr: *const u8) -> [u8; N] {
    // SAFETY: The caller supplies N bytes. Byte arrays have alignment one
    // and all bit patterns are valid; no reference into input escapes.
    unsafe { ptr.cast::<[u8; N]>().read_unaligned() }
}

fn coordinate(offset: i32, delta: i32) -> Result<i32> {
    offset
        .checked_add(delta)
        .ok_or_else(|| Error::invalid("hic block", "bin overflows i32"))
}

#[inline(always)]
fn sparse<const X: bool, const Y: bool, const FLOAT: bool>(
    r: &mut Reader<'_>,
    x0: i32,
    y0: i32,
    f: &mut impl FnMut(i32, i32, f32),
) -> Result<()> {
    let xw = if X { 4 } else { 2 };
    let width = xw + if FLOAT { 4 } else { 2 };
    for _ in 0..r.count(Y)? {
        let y = coordinate(y0, r.int(Y)?)?;
        let count = r.count(X)?;
        let records = r.records(count, width)?;
        for i in 0..count {
            // SAFETY: records checked count*width bytes; every immutable chunk
            // contains the complete X/FLOAT-selected fields, even unaligned.
            let (x, value) = unsafe {
                let p = records.as_ptr().add(i * width);
                let x = if X {
                    i32::from_le_bytes(field(p))
                } else {
                    i16::from_le_bytes(field(p)) as i32
                };
                let value = if FLOAT {
                    f32::from_le_bytes(field(p.add(xw)))
                } else {
                    i16::from_le_bytes(field(p.add(xw))) as f32
                };
                (x, value)
            };
            f(coordinate(x0, x)?, y, value);
        }
    }
    Ok(())
}

#[inline(always)]
fn dense<const FLOAT: bool>(
    r: &mut Reader<'_>,
    x0: i32,
    y0: i32,
    f: &mut impl FnMut(i32, i32, f32),
) -> Result<()> {
    let points = r.count(true)?;
    let width = r.i16()?;
    if width <= 0 && points > 0 {
        return Err(Error::invalid("hic block", "nonpositive dense width"));
    }
    if points == 0 {
        return Ok(());
    }
    let width = width as usize;
    let vw = if FLOAT { 4 } else { 2 };
    let records = r.records(points, vw)?;
    coordinate(x0, (points.min(width) - 1) as i32)?;
    coordinate(y0, ((points - 1) / width) as i32)?;
    for (row, bytes) in records.chunks(width * vw).enumerate() {
        for col in 0..bytes.len() / vw {
            // SAFETY: Each chunk contains exactly one complete value in the
            // selected representation. All float/integer byte patterns are valid.
            let value = unsafe {
                let value = bytes.as_ptr().add(col * vw);
                if FLOAT {
                    f32::from_le_bytes(field(value))
                } else {
                    let v = i16::from_le_bytes(field(value));
                    if v == i16::MIN {
                        continue;
                    }
                    v as f32
                }
            };
            if !value.is_nan() {
                f(x0 + col as i32, y0 + row as i32, value);
            }
        }
    }
    Ok(())
}

/// Decodes a v8 or v9 block, skipping dense missing-value sentinels.
/// Malformed counts, truncation and coordinate overflow return errors;
/// callbacks may have received records before an error in a later row.
#[inline]
pub fn decode_block(version: i32, block: &[u8], mut f: impl FnMut(i32, i32, f32)) -> Result<()> {
    if !(8..=9).contains(&version) {
        return Err(Error::Unsupported(".hic versions other than 8 and 9"));
    }
    let mut r = Reader { buf: block, at: 0 };
    r.count(true)?;
    let (x0, y0) = (r.i32()?, r.i32()?);
    let float = r.u8()? != 0;
    let (x, y) = if version > 8 {
        (r.u8()? != 0, r.u8()? != 0)
    } else {
        (false, false)
    };
    match r.u8()? {
        1 => match (x, y, float) {
            (false, false, false) => sparse::<false, false, false>(&mut r, x0, y0, &mut f),
            (false, false, true) => sparse::<false, false, true>(&mut r, x0, y0, &mut f),
            (false, true, false) => sparse::<false, true, false>(&mut r, x0, y0, &mut f),
            (false, true, true) => sparse::<false, true, true>(&mut r, x0, y0, &mut f),
            (true, false, false) => sparse::<true, false, false>(&mut r, x0, y0, &mut f),
            (true, false, true) => sparse::<true, false, true>(&mut r, x0, y0, &mut f),
            (true, true, false) => sparse::<true, true, false>(&mut r, x0, y0, &mut f),
            (true, true, true) => sparse::<true, true, true>(&mut r, x0, y0, &mut f),
        },
        2 if float => dense::<true>(&mut r, x0, y0, &mut f),
        2 => dense::<false>(&mut r, x0, y0, &mut f),
        t => Err(Error::invalid(
            "hic block",
            format!("unknown representation {t}"),
        )),
    }
}

/// Widths of the narrow (`short`) block fields.
const SHORT_LIMIT: u32 = i16::MAX as u32;

fn short_value(v: f32) -> bool {
    v.fract() == 0.0 && v > i16::MIN as f32 && v <= i16::MAX as f32
}

/// Encodes one block of a v8 or v9 file, uncompressed, into `out` (replacing its
/// contents). `records` must be sorted by (y, x) with unique cells and
/// finite values. Values are stored as `short` when all are integers that
/// fit, positions as `short` where the block is narrow enough (v9 can widen
/// them; v8 cannot), and the smaller of the list-of-rows and dense
/// layouts is chosen.
pub fn encode_block(version: i32, records: &[Record], out: &mut Vec<u8>) -> Result<()> {
    if !(8..=9).contains(&version) {
        return Err(Error::Unsupported(".hic versions other than 8 and 9"));
    }
    out.clear();
    let (Some(first), Some(last)) = (records.first(), records.last()) else {
        return Err(Error::invalid("hic block", "no records"));
    };
    if records.iter().any(|r| !r.value.is_finite())
        || records
            .windows(2)
            .any(|w| (w[0].y, w[0].x) >= (w[1].y, w[1].x))
    {
        return Err(Error::invalid(
            "hic block",
            "records must be unique, sorted by (y,x), and finite",
        ));
    }
    let (mut x0, mut x_max) = (u32::MAX, 0);
    for r in records {
        x0 = x0.min(r.x);
        x_max = x_max.max(r.x);
    }
    let y0 = first.y;
    if x_max > i32::MAX as u32 || last.y > i32::MAX as u32 {
        return Err(Error::invalid("hic block", "bin exceeds i32"));
    }
    let (w, h) = (x_max - x0 + 1, last.y - y0 + 1);
    let shorts = records.iter().all(|r| short_value(r.value));
    let (wide_x, wide_y) = (w >= SHORT_LIMIT, h >= SHORT_LIMIT);
    if version < 9 && (wide_x || wide_y) {
        return Err(Error::invalid(
            "hic block",
            format!("{w} x {h} bins exceed a v{version} block"),
        ));
    }
    let n = i32::try_from(records.len())
        .map_err(|_| Error::invalid("hic block", "more than 2^31 records"))?;

    let value = if shorts { 2 } else { 4 };
    let (xw, yw) = (if wide_x { 4 } else { 2 }, if wide_y { 4 } else { 2 });
    let mut rows = 0u64;
    let mut previous = None;
    for r in records {
        if previous != Some(r.y) {
            rows += 1;
            previous = Some(r.y);
        }
    }
    let rows_size = yw + rows * (yw + xw) + records.len() as u64 * (xw + value);
    let points = (last.y - y0) as u64 * w as u64 + (last.x - x0) as u64 + 1;
    let dense_size = if w < SHORT_LIMIT && points <= i32::MAX as u64 {
        6 + points * value
    } else {
        u64::MAX
    };
    let dense = dense_size < rows_size;

    out.reserve(14 + rows_size.min(dense_size) as usize);
    out.extend_from_slice(&n.to_le_bytes());
    out.extend_from_slice(&(x0 as i32).to_le_bytes());
    out.extend_from_slice(&(y0 as i32).to_le_bytes());
    out.push(!shorts as u8);
    if version > 8 {
        out.push(wide_x as u8);
        out.push(wide_y as u8);
    }
    let put_value = |out: &mut Vec<u8>, v: f32| {
        if shorts {
            out.extend_from_slice(&(v as i16).to_le_bytes());
        } else {
            out.extend_from_slice(&v.to_le_bytes());
        }
    };
    let put_int = |out: &mut Vec<u8>, wide: bool, v: u32| {
        if wide {
            out.extend_from_slice(&(v as i32).to_le_bytes());
        } else {
            out.extend_from_slice(&(v as i16).to_le_bytes());
        }
    };
    if dense {
        out.push(2);
        out.extend_from_slice(&(points as i32).to_le_bytes());
        out.extend_from_slice(&(w as i16).to_le_bytes());
        let mut next = 0u64;
        for r in records {
            let at = (r.y - y0) as u64 * w as u64 + (r.x - x0) as u64;
            for _ in next..at {
                if shorts {
                    out.extend_from_slice(&i16::MIN.to_le_bytes());
                } else {
                    out.extend_from_slice(&f32::NAN.to_le_bytes());
                }
            }
            put_value(out, r.value);
            next = at + 1;
        }
    } else {
        out.push(1);
        put_int(out, wide_y, rows as u32);
        for row in records.chunk_by(|a, b| a.y == b.y) {
            put_int(out, wide_y, row[0].y - y0);
            put_int(out, wide_x, row.len() as u32);
            for r in row {
                put_int(out, wide_x, r.x - x0);
                put_value(out, r.value);
            }
        }
    }
    Ok(())
}
