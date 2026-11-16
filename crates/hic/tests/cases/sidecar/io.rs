use std::io::{self, Cursor, Write};

use omicsweave_core::sidecar::SourceIdentity;
use omicsweave_core::{Error, ReadAt, Result};
use omicsweave_hic::{Chrom, HicHeader, HicSidecar, HicWriter};

struct ShortWriter {
    bytes: Vec<u8>,
    remaining: usize,
}
impl Write for ShortWriter {
    fn write(&mut self, bytes: &[u8]) -> io::Result<usize> {
        if self.remaining == 0 {
            return Err(io::Error::other("injected sidecar write failure"));
        }
        let len = bytes.len().min(self.remaining).min(3);
        self.bytes.extend_from_slice(&bytes[..len]);
        self.remaining -= len;
        Ok(len)
    }
    fn flush(&mut self) -> io::Result<()> {
        Ok(())
    }
}

struct BrokenSource;
impl ReadAt for BrokenSource {
    fn len(&self) -> u64 {
        100
    }
    fn read_exact_at(&self, _: &mut [u8], _: u64) -> Result<()> {
        Err(Error::from(io::Error::other(
            "injected identity read failure",
        )))
    }
}

#[test]
fn write_to_matches_encoding_handles_short_writes_and_propagates_io_failures() {
    let header = HicHeader::new(
        "g",
        vec![Chrom {
            name: "a".into(),
            len: 100,
        }],
        vec![1],
    );
    let (out, sidecar) = HicWriter::create(Cursor::new(Vec::new()), header, 1)
        .unwrap()
        .finish_with_sidecar()
        .unwrap();
    let source = out.into_inner();
    let expected = sidecar.to_bytes(SourceIdentity::of(source.as_slice()).unwrap());
    let mut short = ShortWriter {
        bytes: Vec::new(),
        remaining: usize::MAX,
    };
    sidecar.write_to(source.as_slice(), &mut short).unwrap();
    assert_eq!(short.bytes, expected);
    assert_eq!(
        HicSidecar::parse(&short.bytes, SourceIdentity::of(source.as_slice()).unwrap()).unwrap(),
        sidecar
    );
    let mut broken = ShortWriter {
        bytes: Vec::new(),
        remaining: 7,
    };
    assert!(sidecar.write_to(source.as_slice(), &mut broken).is_err());
    assert_eq!(broken.bytes, expected[..7]);
    let mut untouched = Vec::new();
    assert!(sidecar.write_to(&BrokenSource, &mut untouched).is_err());
    assert!(untouched.is_empty(), "identity failure must precede output");
    let mut retry = Vec::new();
    sidecar.write_to(source.as_slice(), &mut retry).unwrap();
    assert_eq!(retry, expected);
}
