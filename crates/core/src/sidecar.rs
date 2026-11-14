//! Sidecars: compact binary companions of a source file (`<file>.sidecar`)
//! holding facts a loader needs before reading any data block, such as the
//! decompressed size of every block, that the source format does not store.
//!
//! A sidecar is a derived product: it can be deleted and rebuilt from its
//! source at any time, and it is bound to the source it was built from by a
//! [`SourceIdentity`]. Each format defines its own body; this module frames
//! it.
//!
//! # Container layout
//!
//! Every field is little-endian. The body starts 8-byte aligned and, by
//! convention, holds fixed-width records whose 8-byte fields stay 8-byte
//! aligned, so a reader can index its tables in place.
//!
//! | Offset | Type | Field |
//! |---|---|---|
//! | 0 | `[u8; 8]` | [`SIDECAR_MAGIC`] |
//! | 8 | `u16` | container version ([`SIDECAR_VERSION`]) |
//! | 10 | `u16` | [`SidecarKind`] of the body |
//! | 12 | `u32` | body version, defined by the format |
//! | 16 | `u64` | source length in bytes |
//! | 24 | `u64` | source fingerprint ([`SourceIdentity::of`]) |
//! | 32 | `u64` | body length in bytes |
//! | 40 | `u64` | body checksum ([`checksum`]) |
//! | 48 | | body |

use std::fs::OpenOptions;
use std::io::Write;
use std::path::{Path, PathBuf};
use std::sync::atomic::{AtomicU64, Ordering};

use crate::error::{Error, Result};
use crate::source::ReadAt;

/// First bytes of every sidecar.
pub const SIDECAR_MAGIC: [u8; 8] = *b"OWSIDECR";
/// Version of the container framing.
pub const SIDECAR_VERSION: u16 = 1;
/// Bytes of the container header; the body follows.
pub const SIDECAR_HEADER_LEN: usize = 48;
/// Appended to a source path to name its sidecar.
pub const SIDECAR_SUFFIX: &str = ".sidecar";

/// Which format's body a sidecar holds.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
#[repr(u16)]
pub enum SidecarKind {
    /// bigWig and bigBed (`omicsweave-bbi`).
    Bbi = 1,
    /// Juicer `.hic` (`omicsweave-hic`).
    Hic = 2,
}

impl SidecarKind {
    fn from_u16(v: u16) -> Option<SidecarKind> {
        match v {
            1 => Some(SidecarKind::Bbi),
            2 => Some(SidecarKind::Hic),
            _ => None,
        }
    }
}

/// What a sidecar records about its source: the length and a fingerprint
/// of sampled content (the first and last 64 KiB and 14 evenly spaced
/// 4 KiB pages). It depends on content only, so a source copied or moved
/// together with its sidecar stays valid; a rewrite that keeps the length
/// and changes no sampled byte goes unnoticed.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub struct SourceIdentity {
    pub len: u64,
    pub fingerprint: u64,
}

const EDGE: u64 = 64 * 1024;
const PAGE: u64 = 4096;
const PAGES: u64 = 14;

impl SourceIdentity {
    /// Reads the sampled ranges of `source` and fingerprints them.
    pub fn of<S: ReadAt + ?Sized>(source: &S) -> Result<SourceIdentity> {
        let len = source.len();
        let mut ranges = [(0, 0); PAGES as usize + 2];
        ranges[0] = (0, len.min(EDGE));
        ranges[1] = (len.saturating_sub(EDGE), len);
        for k in 1..=PAGES {
            let at = (len / (PAGES + 1)).saturating_mul(k);
            ranges[k as usize + 1] = (at, (at + PAGE).min(len));
        }
        let mut h = Hasher::new(len);
        let mut buf = Vec::new();
        if len <= 2 * EDGE {
            // The two edge samples cover the whole source, so one read
            // serves every sample without fetching any unsampled bytes.
            buf.resize(len as usize, 0);
            source.read_exact_at(&mut buf, 0)?;
            for (lo, hi) in ranges {
                h.write_u64(lo);
                h.write(&buf[lo as usize..hi as usize]);
            }
        } else {
            for (lo, hi) in ranges {
                buf.resize((hi - lo) as usize, 0);
                source.read_exact_at(&mut buf, lo)?;
                h.write_u64(lo);
                h.write(&buf);
            }
        }
        Ok(SourceIdentity {
            len,
            fingerprint: h.finish(),
        })
    }
}

/// The sidecar path of `source`: its path with [`SIDECAR_SUFFIX`] appended
/// (`a.bw` -> `a.bw.sidecar`).
pub fn sidecar_path(source: impl AsRef<Path>) -> PathBuf {
    let mut name = source.as_ref().as_os_str().to_owned();
    name.push(SIDECAR_SUFFIX);
    PathBuf::from(name)
}

/// The 64-bit checksum the container stores for its body.
pub fn checksum(bytes: &[u8]) -> u64 {
    let mut h = Hasher::new(bytes.len() as u64);
    h.write(bytes);
    h.finish()
}

/// A word-at-a-time multiply-rotate hash with a splitmix64 finalizer: fast,
/// deterministic across platforms, not cryptographic.
struct Hasher(u64);

impl Hasher {
    fn new(seed: u64) -> Hasher {
        Hasher(0x243f_6a88_85a3_08d3 ^ seed)
    }

    fn write_u64(&mut self, w: u64) {
        self.0 = (self.0 ^ w)
            .wrapping_mul(0x9e37_79b9_7f4a_7c15)
            .rotate_left(31);
    }

    fn write(&mut self, bytes: &[u8]) {
        let mut words = bytes.chunks_exact(8);
        for w in &mut words {
            self.write_u64(u64::from_le_bytes(w.try_into().unwrap()));
        }
        let rest = words.remainder();
        let mut last = [0u8; 8];
        last[..rest.len()].copy_from_slice(rest);
        self.write_u64(u64::from_le_bytes(last) ^ (rest.len() as u64) << 56);
    }

    fn finish(&self) -> u64 {
        let mut z = self.0;
        z = (z ^ (z >> 30)).wrapping_mul(0xbf58_476d_1ce4_e5b9);
        z = (z ^ (z >> 27)).wrapping_mul(0x94d0_49bb_1331_11eb);
        z ^ (z >> 31)
    }
}

/// Frames `body` as a sidecar of `kind` bound to `source`.
pub fn seal(kind: SidecarKind, body_version: u32, source: SourceIdentity, body: &[u8]) -> Vec<u8> {
    let mut out = Vec::with_capacity(SIDECAR_HEADER_LEN + body.len());
    out.extend_from_slice(&SIDECAR_MAGIC);
    out.extend_from_slice(&SIDECAR_VERSION.to_le_bytes());
    out.extend_from_slice(&(kind as u16).to_le_bytes());
    out.extend_from_slice(&body_version.to_le_bytes());
    out.extend_from_slice(&source.len.to_le_bytes());
    out.extend_from_slice(&source.fingerprint.to_le_bytes());
    out.extend_from_slice(&(body.len() as u64).to_le_bytes());
    out.extend_from_slice(&checksum(body).to_le_bytes());
    out.extend_from_slice(body);
    out
}

/// Checks a sidecar's framing and checksum, and that it was built for
/// `source` ([`Error::StaleSidecar`] otherwise); returns the body version
/// and the body.
pub fn unseal(bytes: &[u8], kind: SidecarKind, source: SourceIdentity) -> Result<(u32, &[u8])> {
    let mut r = BodyReader::new(bytes, "sidecar header");
    let magic = r.bytes(8)?;
    if magic != SIDECAR_MAGIC {
        return Err(Error::BadMagic {
            what: "sidecar",
            found: u32::from_le_bytes(magic[..4].try_into().unwrap()),
        });
    }
    let version = r.u16()?;
    if version != SIDECAR_VERSION {
        return Err(Error::invalid(
            "sidecar",
            format!("container version {version}, expected {SIDECAR_VERSION}"),
        ));
    }
    let found = r.u16()?;
    if SidecarKind::from_u16(found) != Some(kind) {
        return Err(Error::invalid(
            "sidecar",
            format!("kind {found}, expected {}", kind as u16),
        ));
    }
    let body_version = r.u32()?;
    let identity = SourceIdentity {
        len: r.u64()?,
        fingerprint: r.u64()?,
    };
    let body_len = r.u64()?;
    let sum = r.u64()?;
    let body = &bytes[SIDECAR_HEADER_LEN..];
    if body.len() as u64 != body_len {
        return Err(Error::Truncated {
            what: "sidecar body",
            need: body_len,
            have: body.len() as u64,
        });
    }
    if checksum(body) != sum {
        return Err(Error::invalid("sidecar", "body checksum mismatch"));
    }
    if identity != source {
        return Err(Error::StaleSidecar);
    }
    Ok((body_version, body))
}

/// Writes `bytes` to `path` through a temporary file in the same directory,
/// so a reader never sees a partial sidecar.
pub fn write_file(path: &Path, bytes: &[u8]) -> Result<()> {
    static NEXT_TEMP: AtomicU64 = AtomicU64::new(0);

    // create_new reserves each writer's inode exclusively, including when
    // another process writes the same sidecar or leaves a temporary file.
    for _ in 0..128 {
        let tmp = path.with_file_name(format!(
            ".omicsweave-sidecar-{}-{}.tmp",
            std::process::id(),
            NEXT_TEMP.fetch_add(1, Ordering::Relaxed)
        ));
        let mut file = match OpenOptions::new().write(true).create_new(true).open(&tmp) {
            Ok(file) => file,
            Err(error) if error.kind() == std::io::ErrorKind::AlreadyExists => continue,
            Err(error) => return Err(error.into()),
        };
        let pending = PendingSidecar(tmp);
        file.write_all(bytes)?;
        drop(file);
        std::fs::rename(&pending.0, path)?;
        return Ok(());
    }
    Err(std::io::Error::new(
        std::io::ErrorKind::AlreadyExists,
        "could not reserve a sidecar temporary file",
    )
    .into())
}

struct PendingSidecar(PathBuf);

impl Drop for PendingSidecar {
    fn drop(&mut self) {
        // A successful rename has removed this name; on failure only this
        // writer's temporary file is removed, leaving the target untouched.
        let _ = std::fs::remove_file(&self.0);
    }
}

/// Little-endian appends for building sidecar bodies.
#[derive(Debug, Default)]
pub struct BodyWriter(pub Vec<u8>);

impl BodyWriter {
    pub fn u8(&mut self, v: u8) {
        self.0.push(v);
    }
    pub fn u16(&mut self, v: u16) {
        self.0.extend_from_slice(&v.to_le_bytes());
    }
    pub fn u32(&mut self, v: u32) {
        self.0.extend_from_slice(&v.to_le_bytes());
    }
    pub fn u64(&mut self, v: u64) {
        self.0.extend_from_slice(&v.to_le_bytes());
    }
}

/// Bounds-checked little-endian reads of a sidecar body.
#[derive(Debug)]
pub struct BodyReader<'a> {
    buf: &'a [u8],
    at: usize,
    what: &'static str,
}

impl<'a> BodyReader<'a> {
    pub fn new(buf: &'a [u8], what: &'static str) -> BodyReader<'a> {
        BodyReader { buf, at: 0, what }
    }

    /// Bytes not read yet.
    pub fn remaining(&self) -> usize {
        self.buf.len() - self.at
    }

    pub fn bytes(&mut self, n: usize) -> Result<&'a [u8]> {
        let end = self.at.checked_add(n).filter(|&e| e <= self.buf.len());
        let Some(end) = end else {
            return Err(Error::Truncated {
                what: self.what,
                need: (self.at as u64).saturating_add(n as u64),
                have: self.buf.len() as u64,
            });
        };
        let out = &self.buf[self.at..end];
        self.at = end;
        Ok(out)
    }

    pub fn u8(&mut self) -> Result<u8> {
        Ok(self.bytes(1)?[0])
    }
    pub fn u16(&mut self) -> Result<u16> {
        Ok(u16::from_le_bytes(self.bytes(2)?.try_into().unwrap()))
    }
    pub fn u32(&mut self) -> Result<u32> {
        Ok(u32::from_le_bytes(self.bytes(4)?.try_into().unwrap()))
    }
    pub fn u64(&mut self) -> Result<u64> {
        Ok(u64::from_le_bytes(self.bytes(8)?.try_into().unwrap()))
    }

    /// A record count whose records take at least `record_len` bytes each,
    /// checked against the bytes left before anything is allocated for it.
    pub fn count(&mut self, record_len: usize) -> Result<usize> {
        let n = self.u64()?;
        match usize::try_from(n) {
            Ok(n)
                if n.checked_mul(record_len)
                    .is_some_and(|b| b <= self.remaining()) =>
            {
                Ok(n)
            }
            _ => Err(Error::invalid(self.what, format!("{n} records cannot fit"))),
        }
    }

    /// Fails unless every byte was read.
    pub fn finish(&self) -> Result<()> {
        if self.at != self.buf.len() {
            return Err(Error::invalid(
                self.what,
                format!("{} trailing bytes", self.buf.len() - self.at),
            ));
        }
        Ok(())
    }
}
