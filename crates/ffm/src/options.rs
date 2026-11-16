use crate::error::{Result, invalid};

/// Allocation limits for file-controlled lengths, in bytes. These are reader
/// policy, not format limits; increase them explicitly for larger collections.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct ReaderOptions {
    pub max_metadata_bytes: usize,
    pub max_block_bytes: usize,
    pub max_stream_bytes: usize,
    /// Logical resident payload (arrays, descriptors and metadata). Temporary
    /// decode buffers, Vec spare capacity and the input source are additional.
    pub max_resident_bytes: usize,
}

impl Default for ReaderOptions {
    fn default() -> Self {
        Self {
            max_metadata_bytes: 256 << 20,
            max_block_bytes: 256 << 20,
            max_stream_bytes: 256 << 20,
            max_resident_bytes: 1 << 30,
        }
    }
}

pub(crate) fn reserve<T>(len: usize, limit: usize, what: &'static str) -> Result<Vec<T>> {
    let bytes = len
        .checked_mul(std::mem::size_of::<T>())
        .filter(|&n| n <= limit && n <= isize::MAX as usize)
        .ok_or_else(|| invalid(what, "decoded size exceeds allocation limit"))?;
    let mut out = Vec::new();
    out.try_reserve_exact(len)
        .map_err(|_| invalid(what, format!("cannot allocate {bytes} bytes")))?;
    Ok(out)
}

pub(crate) fn allocate<T: Default + Clone>(
    len: usize,
    limit: usize,
    what: &'static str,
) -> Result<Vec<T>> {
    let mut out = reserve(len, limit, what)?;
    out.resize(len, T::default());
    Ok(out)
}
