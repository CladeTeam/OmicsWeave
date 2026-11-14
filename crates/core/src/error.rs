/// Every error a format reader reports.
#[derive(Debug, Clone, PartialEq, Eq, thiserror::Error)]
pub enum Error {
    /// A buffer ends before the structure it holds.
    #[error("{what}: need {need} bytes, buffer has {have}")]
    Truncated {
        what: &'static str,
        need: u64,
        have: u64,
    },

    /// A magic number does not identify the expected structure.
    #[error("{what}: bad magic {found:#010x}")]
    BadMagic { what: &'static str, found: u32 },

    /// A byte-swapped (big-endian) file.
    #[error("big-endian files are not supported")]
    BigEndian,

    /// A valid input this implementation does not handle.
    #[error("unsupported: {0}")]
    Unsupported(&'static str),

    /// A file offset or buffer position outside the bytes it must address.
    #[error("{what}: offset {offset} is outside [{lo}, {hi})")]
    OutOfRange {
        what: &'static str,
        offset: u64,
        lo: u64,
        hi: u64,
    },

    /// A structurally invalid field or record.
    #[error("{what}: {detail}")]
    Invalid { what: &'static str, detail: String },

    /// A tree deeper than any well-formed file has (or a cycle).
    #[error("tree deeper than {0} levels")]
    TooDeep(usize),

    /// Text that cannot be parsed, with its 1-based line.
    #[error("line {line}: {detail}")]
    Format { line: u64, detail: String },

    /// Compressed data that does not decompress or check out.
    #[error("decompression: {0}")]
    Inflate(&'static str),

    /// The source has no chromosome of this name.
    #[error("unknown chromosome {0:?}")]
    UnknownChrom(String),

    /// A region with `start > end`.
    #[error("invalid region {start}..{end}")]
    InvalidRegion { start: u32, end: u32 },

    /// A region that ends past its sequence, where output must cover it.
    #[error("region {start}..{end} outside sequence of length {len}")]
    RegionOutOfRange { start: u32, end: u32, len: u32 },

    /// An output buffer of the wrong size.
    #[error("{what}: {have} elements, need {need}")]
    OutputSize {
        what: &'static str,
        have: usize,
        need: usize,
    },

    /// Output was requested before any query staged it.
    #[error("write before query")]
    NotQueried,

    /// A sidecar built for a different source (or an earlier version of it);
    /// rebuild it from the source.
    #[error("sidecar does not belong to this source")]
    StaleSidecar,

    /// An I/O failure of the byte source.
    #[error("I/O ({kind:?}): {detail}")]
    Io {
        kind: std::io::ErrorKind,
        detail: String,
    },
}

impl Error {
    pub fn invalid(what: &'static str, detail: impl Into<String>) -> Error {
        Error::Invalid {
            what,
            detail: detail.into(),
        }
    }
}

impl From<std::io::Error> for Error {
    fn from(e: std::io::Error) -> Error {
        Error::Io {
            kind: e.kind(),
            detail: e.to_string(),
        }
    }
}

pub type Result<T> = std::result::Result<T, Error>;

/// Fails unless an output buffer has exactly `need` elements.
pub fn check_output(what: &'static str, have: usize, need: usize) -> Result<()> {
    if have != need {
        return Err(Error::OutputSize { what, have, need });
    }
    Ok(())
}
