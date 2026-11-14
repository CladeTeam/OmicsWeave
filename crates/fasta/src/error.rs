pub use omicsweave_core::{Error, Result};

pub(crate) fn invalid(what: &'static str, detail: impl Into<String>) -> Error {
    Error::invalid(what, detail)
}

pub(crate) fn require(buf: &[u8], need: usize, what: &'static str) -> Result<()> {
    if buf.len() < need {
        return Err(Error::Truncated {
            what,
            need: need as u64,
            have: buf.len() as u64,
        });
    }
    Ok(())
}

#[inline]
pub(crate) fn u16_at(buf: &[u8], at: usize) -> u16 {
    u16::from_le_bytes([buf[at], buf[at + 1]])
}

#[inline]
pub(crate) fn u32_at(buf: &[u8], at: usize) -> u32 {
    u32::from_le_bytes(buf[at..at + 4].try_into().unwrap())
}

#[inline]
pub(crate) fn u64_at(buf: &[u8], at: usize) -> u64 {
    u64::from_le_bytes(buf[at..at + 8].try_into().unwrap())
}
