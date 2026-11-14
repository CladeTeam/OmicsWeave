pub use omicsweave_core::{Error, Result};

pub(crate) fn invalid(what: &'static str, detail: impl Into<String>) -> Error {
    Error::invalid(what, detail)
}
