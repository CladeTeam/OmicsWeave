//! Plain, gzip and BGZF inputs exercise the same conversion semantics.

use std::io::Write;
use std::path::Path;

/// The three storages an input can have.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Storage {
    Plain,
    Gzip,
    Bgzf,
}

pub const STORAGES: [Storage; 3] = [Storage::Plain, Storage::Gzip, Storage::Bgzf];

/// Writes `text` to `path` in `storage`.
pub fn write_input(path: &Path, text: &[u8], storage: Storage) {
    let bytes = match storage {
        Storage::Plain => text.to_vec(),
        Storage::Gzip => {
            let mut e = flate2::write::GzEncoder::new(Vec::new(), flate2::Compression::new(3));
            e.write_all(text).unwrap();
            e.finish().unwrap()
        }
        Storage::Bgzf => {
            let mut w = omicsweave_bgzf::Writer::new(
                Vec::new(),
                1,
                Box::new(omicsweave_core::exec::Inline::new()),
            )
            .unwrap();
            w.write(text).unwrap();
            w.finish().unwrap().out
        }
    };
    std::fs::write(path, bytes).unwrap();
}

/// The text of a BGZF file.
pub fn inflate(path: &Path) -> Vec<u8> {
    use std::io::Read;
    let raw = std::fs::read(path).unwrap();
    let mut out = Vec::new();
    omicsweave_bgzf::Reader::new(&raw[..], Box::new(omicsweave_core::exec::Inline::new()))
        .read_to_end(&mut out)
        .unwrap();
    out
}
