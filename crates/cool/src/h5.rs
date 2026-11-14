//! Metadata access through libhdf5: attributes, small whole datasets, and
//! the storage facts (layout, filters, chunk addresses) the parallel pixel
//! path needs. Everything here runs under hdf5-metno's global lock, so it
//! is serial; it is only used when opening a file and by the library
//! fallback for pixel columns the fast path cannot decode.

use hdf5_metno as hdf5;
use hdf5_metno_sys as sys;

use hdf5::Location;
use hdf5::types::{FloatSize, IntSize, TypeDescriptor, VarLenAscii, VarLenUnicode};
use omicsweave_core::{Error, Result};

/// Wraps a libhdf5 failure with what was being read.
pub(crate) fn err(what: &'static str) -> impl FnOnce(hdf5::Error) -> Error {
    move |e| Error::invalid(what, e.to_string())
}

/// Opens a named member dataset, naming it in the error when missing.
pub(crate) fn dataset(group: &hdf5::Group, path: &'static str) -> Result<hdf5::Dataset> {
    if !group.link_exists(path) {
        return Err(Error::invalid(
            "cool layout",
            format!("missing dataset {path}"),
        ));
    }
    group.dataset(path).map_err(err(path))
}

/// A 1-D dataset's length.
pub(crate) fn len_1d(ds: &hdf5::Dataset, what: &'static str) -> Result<u64> {
    let shape = ds.shape();
    match shape.as_slice() {
        [n] => Ok(*n as u64),
        _ => Err(Error::invalid(
            what,
            format!("expected a 1-D dataset, shape is {shape:?}"),
        )),
    }
}

fn has_attr(loc: &Location, name: &str) -> Result<bool> {
    let names = loc.attr_names().map_err(err("hdf5 attributes"))?;
    Ok(names.iter().any(|n| n == name))
}

/// A string attribute (variable- or fixed-length, scalar or the first
/// element of an array), or `None` when absent.
pub(crate) fn attr_string(loc: &Location, name: &'static str) -> Result<Option<String>> {
    if !has_attr(loc, name)? {
        return Ok(None);
    }
    let attr = loc.attr(name).map_err(err(name))?;
    let dtype = attr.dtype().map_err(err(name))?;
    let desc = dtype.to_descriptor().map_err(err(name))?;
    let string = match desc {
        TypeDescriptor::VarLenUnicode => attr
            .read_raw::<VarLenUnicode>()
            .map_err(err(name))?
            .first()
            .map(|s| s.as_str().to_owned())
            .ok_or_else(|| Error::invalid(name, "empty string attribute"))?,
        TypeDescriptor::VarLenAscii => attr
            .read_raw::<VarLenAscii>()
            .map_err(err(name))?
            .first()
            .map(|s| s.as_str().to_owned())
            .ok_or_else(|| Error::invalid(name, "empty string attribute"))?,
        TypeDescriptor::FixedAscii(n) | TypeDescriptor::FixedUnicode(n) => {
            let space = attr.space().map_err(err(name))?;
            let _lock = sys::LOCK.lock();
            // SAFETY: `space` owns a valid dataspace id throughout the query.
            let count = unsafe { sys::h5s::H5Sget_simple_extent_npoints(space.id()) };
            let count = usize::try_from(count)
                .map_err(|_| Error::invalid(name, "invalid string attribute size"))?;
            if count == 0 {
                return Err(Error::invalid(name, "empty string attribute"));
            }
            let shape_count = space
                .shape()
                .iter()
                .try_fold(1usize, |n, &dim| n.checked_mul(dim))
                .ok_or_else(|| Error::invalid(name, "string attribute shape overflow"))?;
            if count != shape_count {
                return Err(Error::invalid(name, "inconsistent string attribute size"));
            }
            let mut buf = fixed_buffer(count, n, name)?;
            // SAFETY: checked multiplication and successful allocation in
            // fixed_buffer give `count * n` initialized bytes. The memory
            // datatype is the attribute's own fixed-size `n`-byte type and
            // the queried dataspace has exactly `count` elements.
            let status =
                unsafe { sys::h5a::H5Aread(attr.id(), dtype.id(), buf.as_mut_ptr().cast()) };
            if status < 0 {
                return Err(Error::invalid(name, "H5Aread failed"));
            }
            fixed_string(&buf[..n])?
        }
        other => {
            return Err(Error::invalid(
                name,
                format!("expected a string attribute, found {other:?}"),
            ));
        }
    };
    Ok(Some(string))
}

/// An integer attribute (scalar or the first element of an array), or
/// `None` when absent.
pub(crate) fn attr_int(loc: &Location, name: &'static str) -> Result<Option<i64>> {
    if !has_attr(loc, name)? {
        return Ok(None);
    }
    let attr = loc.attr(name).map_err(err(name))?;
    let desc = attr
        .dtype()
        .and_then(|t| t.to_descriptor())
        .map_err(err(name))?;
    if !matches!(
        desc,
        TypeDescriptor::Integer(_) | TypeDescriptor::Unsigned(_)
    ) {
        return Err(Error::invalid(
            name,
            format!("expected an integer attribute, found {desc:?}"),
        ));
    }
    let v = attr.read_raw::<i64>().map_err(err(name))?;
    v.first()
        .copied()
        .map(Some)
        .ok_or_else(|| Error::invalid(name, "empty integer attribute"))
}

fn fixed_buffer(count: usize, size: usize, what: &'static str) -> Result<Vec<u8>> {
    if size == 0 {
        return Err(Error::invalid(what, "zero-width fixed string"));
    }
    let bytes = count
        .checked_mul(size)
        .ok_or_else(|| Error::invalid(what, "fixed string buffer size overflow"))?;
    let mut buf = Vec::new();
    buf.try_reserve_exact(bytes)
        .map_err(|_| Error::invalid(what, "fixed string buffer allocation failed"))?;
    buf.resize(bytes, 0);
    Ok(buf)
}

fn fixed_string(s: &[u8]) -> Result<String> {
    let end = s.iter().position(|&b| b == 0).unwrap_or(s.len());
    let s = s[..end].trim_ascii_end();
    String::from_utf8(s.to_vec())
        .map_err(|_| Error::invalid("hdf5 string", "fixed-length string is not UTF-8"))
}

/// A whole 1-D string dataset (variable- or fixed-length strings).
pub(crate) fn read_strings(ds: &hdf5::Dataset, what: &'static str) -> Result<Vec<String>> {
    let n = len_1d(ds, what)? as usize;
    let dtype = ds.dtype().map_err(err(what))?;
    match dtype.to_descriptor().map_err(err(what))? {
        TypeDescriptor::VarLenUnicode => Ok(ds
            .read_raw::<VarLenUnicode>()
            .map_err(err(what))?
            .iter()
            .map(|s| s.as_str().to_owned())
            .collect()),
        TypeDescriptor::VarLenAscii => Ok(ds
            .read_raw::<VarLenAscii>()
            .map_err(err(what))?
            .iter()
            .map(|s| s.as_str().to_owned())
            .collect()),
        TypeDescriptor::FixedAscii(size) | TypeDescriptor::FixedUnicode(size) => {
            let mut buf = fixed_buffer(n, size, what)?;
            if n > 0 {
                let _lock = sys::LOCK.lock();
                // SAFETY: fixed_buffer checked `n * size` and allocated that
                // many initialized bytes. This 1-D dataset has `n` elements
                // in its own fixed-size `size`-byte memory datatype.
                let status = unsafe {
                    sys::h5d::H5Dread(
                        ds.id(),
                        dtype.id(),
                        sys::h5s::H5S_ALL,
                        sys::h5s::H5S_ALL,
                        sys::h5p::H5P_DEFAULT,
                        buf.as_mut_ptr().cast(),
                    )
                };
                if status < 0 {
                    return Err(Error::invalid(what, "H5Dread failed"));
                }
            }
            buf.chunks_exact(size).map(fixed_string).collect()
        }
        other => Err(Error::invalid(
            what,
            format!("expected strings, found {other:?}"),
        )),
    }
}

/// A whole 1-D integer (or integer-based enum) dataset, as stored values.
pub(crate) fn read_ints(ds: &hdf5::Dataset, what: &'static str) -> Result<Vec<i64>> {
    len_1d(ds, what)?;
    let desc = ds
        .dtype()
        .and_then(|t| t.to_descriptor())
        .map_err(err(what))?;
    if !matches!(
        desc,
        TypeDescriptor::Integer(_) | TypeDescriptor::Unsigned(_) | TypeDescriptor::Enum(_)
    ) {
        return Err(Error::invalid(
            what,
            format!("expected integers, found {desc:?}"),
        ));
    }
    ds.read_raw::<i64>().map_err(err(what))
}

/// The enum members of a dataset's type, `(name, value)`, or `None` when
/// the type is not an enum.
pub(crate) fn enum_members(
    ds: &hdf5::Dataset,
    what: &'static str,
) -> Result<Option<Vec<(String, i64)>>> {
    let desc = ds
        .dtype()
        .and_then(|t| t.to_descriptor())
        .map_err(err(what))?;
    let TypeDescriptor::Enum(e) = desc else {
        return Ok(None);
    };
    let bits = int_bytes(e.size) * 8;
    let members = e
        .members
        .iter()
        .map(|m| {
            // Members carry the raw bits of their value; sign-extend
            // narrow signed values.
            let v = if e.signed && bits < 64 {
                ((m.value << (64 - bits)) as i64) >> (64 - bits)
            } else {
                m.value as i64
            };
            (m.name.clone(), v)
        })
        .collect();
    Ok(Some(members))
}

fn int_bytes(s: IntSize) -> usize {
    match s {
        IntSize::U1 => 1,
        IntSize::U2 => 2,
        IntSize::U4 => 4,
        IntSize::U8 => 8,
    }
}

/// Float widths decoded outside libhdf5 (half floats, when the bindings
/// know them, are left to libhdf5).
#[allow(unreachable_patterns)]
fn float_bytes(s: FloatSize) -> Option<usize> {
    match s {
        FloatSize::U4 => Some(4),
        FloatSize::U8 => Some(8),
        _ => None,
    }
}

/// How a column's stored elements are laid out, for decoding them outside
/// libhdf5.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) struct Elem {
    pub kind: ElemKind,
    pub size: usize,
    pub big_endian: bool,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) enum ElemKind {
    Signed,
    Unsigned,
    Float,
}

/// A column's element type, or `None` when it is not a plain integer or
/// IEEE float of 1–8 bytes in a fixed byte order using all its bits (such
/// columns are read through libhdf5, which converts them).
pub(crate) fn elem(ds: &hdf5::Dataset, what: &'static str) -> Result<Option<Elem>> {
    let dtype = ds.dtype().map_err(err(what))?;
    let desc = dtype.to_descriptor().map_err(err(what))?;
    let (kind, size) = match desc {
        TypeDescriptor::Integer(s) => (ElemKind::Signed, int_bytes(s)),
        TypeDescriptor::Unsigned(s) => (ElemKind::Unsigned, int_bytes(s)),
        TypeDescriptor::Float(s) => match float_bytes(s) {
            Some(n) => (ElemKind::Float, n),
            None => return Ok(None),
        },
        TypeDescriptor::Enum(_) => return Ok(None),
        other => {
            return Err(Error::invalid(
                what,
                format!("expected numbers, found {other:?}"),
            ));
        }
    };
    let big_endian = match dtype.byte_order() {
        hdf5::datatype::ByteOrder::LittleEndian => false,
        hdf5::datatype::ByteOrder::BigEndian => true,
        _ => return Ok(None),
    };
    let _lock = sys::LOCK.lock();
    // SAFETY: plain queries on a valid datatype id.
    let (precision, offset) = unsafe {
        (
            sys::h5t::H5Tget_precision(dtype.id()),
            sys::h5t::H5Tget_offset(dtype.id()),
        )
    };
    if precision != size * 8 || offset != 0 {
        return Ok(None);
    }
    Ok(Some(Elem {
        kind,
        size,
        big_endian,
    }))
}

#[cfg(test)]
mod tests {
    use super::*;
    use hdf5::types::FixedAscii;

    #[test]
    fn fixed_buffer_rejects_overflow_before_ffi() {
        assert!(fixed_buffer(usize::MAX, 2, "test").is_err());
        assert!(fixed_buffer(isize::MAX as usize + 1, 1, "test").is_err());
        assert!(fixed_buffer(1, 0, "test").is_err());
        assert_eq!(fixed_buffer(3, 4, "test").unwrap(), vec![0; 12]);
    }

    #[test]
    fn fixed_attribute_scalar_array_and_empty() {
        let dir = tempfile::tempdir().unwrap();
        let file = hdf5::File::create(dir.path().join("attrs.h5")).unwrap();
        let first = FixedAscii::<8>::from_ascii("first").unwrap();
        let second = FixedAscii::<8>::from_ascii("second").unwrap();
        file.new_attr::<FixedAscii<8>>()
            .create("scalar")
            .unwrap()
            .write_scalar(&first)
            .unwrap();
        file.new_attr::<FixedAscii<8>>()
            .shape((1, 2))
            .create("array")
            .unwrap()
            .write_raw(&[first, second])
            .unwrap();
        file.new_attr::<FixedAscii<8>>()
            .shape(0)
            .create("empty")
            .unwrap();
        assert_eq!(
            attr_string(&file, "scalar").unwrap().as_deref(),
            Some("first")
        );
        assert_eq!(
            attr_string(&file, "array").unwrap().as_deref(),
            Some("first")
        );
        assert!(attr_string(&file, "empty").is_err());
    }
}
