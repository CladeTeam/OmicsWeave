//! Writing UCSC 2bit.
//!
//! A record's header lists its N and mask blocks before the packed bases,
//! and the directory of record offsets precedes all records, so a sequence
//! is encoded in pieces ([`TwoBitPiece`], independent of each other and
//! cheap to run in parallel), their runs joined in order
//! ([`TwoBitRecord`]), and the record written once complete
//! ([`TwoBitWriter`]); the directory is filled in at the end.

use std::io::{Read, Seek, SeekFrom, Write};

use crate::error::{Result, invalid};
use crate::twobit::{HEADER_LEN, TWOBIT_MAGIC};

/// 2bit code of each byte (T 0, C 1, A 2, G 3); every other byte is stored
/// as N, packed as T as the UCSC tools pack it.
const fn code_table() -> [u8; 256] {
    let mut t = [0u8; 256];
    t[b'C' as usize] = 1;
    t[b'c' as usize] = 1;
    t[b'A' as usize] = 2;
    t[b'a' as usize] = 2;
    t[b'G' as usize] = 3;
    t[b'g' as usize] = 3;
    t
}
static CODE: [u8; 256] = code_table();

fn is_acgt(c: u8) -> bool {
    matches!(c, b'A' | b'C' | b'G' | b'T' | b'a' | b'c' | b'g' | b't')
}

/// Runs of `[start, end)` positions where `hit` holds, offset by `base`.
fn runs(bases: &[u8], base: u64, hit: impl Fn(u8) -> bool, out: &mut Vec<(u64, u64)>) {
    let mut i = 0;
    while i < bases.len() {
        if !hit(bases[i]) {
            i += 1;
            continue;
        }
        let start = i;
        while i < bases.len() && hit(bases[i]) {
            i += 1;
        }
        out.push((base + start as u64, base + i as u64));
    }
}

/// Bases `[start, start + len)` of one sequence, packed with their N and
/// soft-mask (lowercase) runs. Bases other than A, C, G and T (either case)
/// become N; IUPAC codes cannot be stored in 2bit.
#[derive(Debug, Clone, Default, PartialEq, Eq)]
pub struct TwoBitPiece {
    pub start: u64,
    pub len: u64,
    pub packed: Vec<u8>,
    pub n_runs: Vec<(u64, u64)>,
    pub mask_runs: Vec<(u64, u64)>,
}

impl TwoBitPiece {
    /// Encodes `bases` starting at sequence position `start`, which must be
    /// a multiple of 4 so that pieces pack to whole bytes; only a
    /// sequence's last piece may have a length that is not.
    pub fn encode(start: u64, bases: &[u8]) -> Result<TwoBitPiece> {
        if !start.is_multiple_of(4) {
            return Err(invalid("2bit piece", "start is not a multiple of 4"));
        }
        start
            .checked_add(bases.len() as u64)
            .ok_or_else(|| invalid("2bit piece", "end overflows u64"))?;
        let mut packed = Vec::with_capacity(bases.len().div_ceil(4));
        let mut quads = bases.chunks_exact(4);
        for q in &mut quads {
            packed.push(
                CODE[q[0] as usize] << 6
                    | CODE[q[1] as usize] << 4
                    | CODE[q[2] as usize] << 2
                    | CODE[q[3] as usize],
            );
        }
        let rest = quads.remainder();
        if !rest.is_empty() {
            let mut b = 0u8;
            for (i, &c) in rest.iter().enumerate() {
                b |= CODE[c as usize] << (6 - 2 * i);
            }
            packed.push(b);
        }
        let mut n_runs = Vec::new();
        let mut mask_runs = Vec::new();
        runs(bases, start, |c| !is_acgt(c), &mut n_runs);
        runs(bases, start, |c| c.is_ascii_lowercase(), &mut mask_runs);
        Ok(TwoBitPiece {
            start,
            len: bases.len() as u64,
            packed,
            n_runs,
            mask_runs,
        })
    }
}

/// The header of one sequence's record, gathered from its pieces in order.
#[derive(Debug, Clone, Default, PartialEq, Eq)]
pub struct TwoBitRecord {
    len: u64,
    n_runs: Vec<(u64, u64)>,
    mask_runs: Vec<(u64, u64)>,
}

fn join(all: &mut Vec<(u64, u64)>, more: &[(u64, u64)]) {
    let mut more = more.iter();
    if let (Some(last), Some(&first)) = (all.last_mut(), more.as_slice().first())
        && last.1 == first.0
    {
        last.1 = first.1;
        more.next();
    }
    all.extend(more);
}

impl TwoBitRecord {
    pub fn new() -> TwoBitRecord {
        TwoBitRecord::default()
    }

    /// Adds the runs of the next piece, which must start where the previous
    /// one ended; runs that meet at the boundary become one block. The
    /// caller keeps `piece.packed`, in order, for [`TwoBitWriter::add`].
    /// Rejects inconsistent packed lengths or runs outside the piece,
    /// without changing this record.
    pub fn push(&mut self, piece: &TwoBitPiece) -> Result<()> {
        if piece.start != self.len || !self.len.is_multiple_of(4) {
            return Err(invalid(
                "2bit record",
                "pieces must follow each other and all but the last span a multiple of 4 bases",
            ));
        }
        let end = self
            .len
            .checked_add(piece.len)
            .ok_or_else(|| invalid("2bit record", "length overflows u64"))?;
        if piece.packed.len() as u64 != piece.len.div_ceil(4) {
            return Err(invalid(
                "2bit piece",
                "packed length differs from base length",
            ));
        }
        for list in [&piece.n_runs, &piece.mask_runs] {
            let mut previous_end = piece.start;
            for &(start, stop) in list {
                if start < previous_end || start >= stop || stop > end {
                    return Err(invalid(
                        "2bit piece",
                        "runs must be nonempty, ordered and within the piece",
                    ));
                }
                previous_end = stop;
            }
        }
        join(&mut self.n_runs, &piece.n_runs);
        join(&mut self.mask_runs, &piece.mask_runs);
        self.len = end;
        Ok(())
    }

    /// Bases so far.
    pub fn len(&self) -> u64 {
        self.len
    }

    pub fn is_empty(&self) -> bool {
        self.len == 0
    }

    /// Bytes of the packed bases.
    pub fn packed_len(&self) -> u64 {
        self.len.div_ceil(4)
    }

    /// Bytes of the record header.
    pub fn header_len(&self) -> u64 {
        16 + 8 * (self.n_runs.len() + self.mask_runs.len()) as u64
    }

    fn header(&self) -> Result<Vec<u8>> {
        let len = u32::try_from(self.len)
            .map_err(|_| invalid("2bit record", "sequence longer than 2^32 - 1 bases"))?;
        let mut h = Vec::with_capacity(self.header_len() as usize);
        h.extend_from_slice(&len.to_le_bytes());
        for list in [&self.n_runs, &self.mask_runs] {
            h.extend_from_slice(&(list.len() as u32).to_le_bytes());
            for &(s, _) in list.iter() {
                h.extend_from_slice(&(s as u32).to_le_bytes());
            }
            for &(s, e) in list.iter() {
                h.extend_from_slice(&((e - s) as u32).to_le_bytes());
            }
        }
        h.extend_from_slice(&0u32.to_le_bytes());
        Ok(h)
    }
}

/// Size of a 2bit file of `seqs` (name, bases) with no N or mask blocks; a
/// lower bound that decides whether 32-bit offsets can address the file.
pub fn twobit_base_size<'a>(seqs: impl IntoIterator<Item = (&'a str, u64)>) -> u64 {
    let mut size = HEADER_LEN as u64;
    for (name, len) in seqs {
        size += 1 + name.len() as u64 + 8 + 16 + len.div_ceil(4);
    }
    size
}

/// Writes a 2bit file whose sequence names are known up front (the
/// directory precedes the records). Records are added in the order of the
/// names; the directory is written by [`TwoBitWriter::finish`].
pub struct TwoBitWriter<W: Write + Seek> {
    out: W,
    version: u32,
    names: Vec<String>,
    offsets: Vec<u64>,
    pos: u64,
    broken: bool,
}

impl<W: Write + Seek> TwoBitWriter<W> {
    /// Writes the header and a directory placeholder to `out`, which must
    /// start empty at position 0. `version` 1 stores 64-bit offsets, needed
    /// for files of 4 GiB or more ([`twobit_base_size`]).
    pub fn create(mut out: W, names: Vec<String>, version: u32) -> Result<TwoBitWriter<W>> {
        if version > 1 {
            return Err(invalid("2bit writer", format!("unknown version {version}")));
        }
        let count = u32::try_from(names.len())
            .map_err(|_| invalid("2bit writer", "more than 2^32 - 1 sequences"))?;
        let mut seen = std::collections::HashSet::with_capacity(names.len());
        for n in &names {
            if n.is_empty() || n.len() > 255 {
                return Err(invalid(
                    "2bit writer",
                    format!("sequence name {n:?} must have 1 to 255 bytes"),
                ));
            }
            if !seen.insert(n.as_str()) {
                return Err(invalid("2bit writer", format!("duplicate name {n:?}")));
            }
        }
        let width = if version == 1 { 8 } else { 4 };
        let mut head = Vec::with_capacity(HEADER_LEN);
        head.extend_from_slice(&TWOBIT_MAGIC.to_le_bytes());
        head.extend_from_slice(&version.to_le_bytes());
        head.extend_from_slice(&count.to_le_bytes());
        head.extend_from_slice(&0u32.to_le_bytes());
        for n in &names {
            head.push(n.len() as u8);
            head.extend_from_slice(n.as_bytes());
            head.extend(std::iter::repeat_n(0u8, width));
        }
        out.write_all(&head)?;
        Ok(TwoBitWriter {
            out,
            version,
            offsets: Vec::with_capacity(names.len()),
            names,
            pos: head.len() as u64,
            broken: false,
        })
    }

    /// Name of the next record to add.
    pub fn next_name(&self) -> Option<&str> {
        self.names.get(self.offsets.len()).map(String::as_str)
    }

    /// Writes the next record: its header from `record` and its packed
    /// bases, `record.packed_len()` bytes read from `packed`.
    /// A failed write or short input leaves an incomplete record in the
    /// sink; subsequent writes and [`TwoBitWriter::finish`] then fail.
    pub fn add(&mut self, name: &str, record: &TwoBitRecord, packed: &mut impl Read) -> Result<()> {
        self.check_intact()?;
        match self.next_name() {
            Some(next) if next == name => {}
            Some(next) => {
                return Err(invalid(
                    "2bit writer",
                    format!("record {name:?} out of order; expected {next:?}"),
                ));
            }
            None => {
                return Err(invalid(
                    "2bit writer",
                    format!("unexpected record {name:?}"),
                ));
            }
        }
        if self.version == 0 && self.pos > u32::MAX as u64 {
            return Err(invalid(
                "2bit writer",
                "file exceeds 4 GiB; write version 1 (64-bit offsets)",
            ));
        }
        let header = record.header()?;
        let want = record.packed_len();
        let end = self
            .pos
            .checked_add(header.len() as u64)
            .and_then(|pos| pos.checked_add(want))
            .ok_or_else(|| invalid("2bit writer", "file position overflows u64"))?;
        let result = (|| {
            self.out.write_all(&header)?;
            let copied = std::io::copy(&mut packed.take(want), &mut self.out)?;
            if copied != want {
                return Err(invalid(
                    "2bit writer",
                    format!("{name}: {copied} packed bytes, need {want}"),
                ));
            }
            Ok(())
        })();
        if let Err(error) = result {
            self.broken = true;
            return Err(error);
        }
        self.offsets.push(self.pos);
        self.pos = end;
        Ok(())
    }

    fn check_intact(&self) -> Result<()> {
        if self.broken {
            return Err(invalid("2bit writer", "an earlier write failed"));
        }
        Ok(())
    }

    /// Fills in the directory once every named record has been added.
    pub fn finish(mut self) -> Result<W> {
        self.check_intact()?;
        if self.offsets.len() != self.names.len() {
            return Err(invalid(
                "2bit writer",
                format!(
                    "{} of {} records written",
                    self.offsets.len(),
                    self.names.len()
                ),
            ));
        }
        let wide = self.version == 1;
        let mut dir = Vec::new();
        for (n, &off) in self.names.iter().zip(&self.offsets) {
            dir.push(n.len() as u8);
            dir.extend_from_slice(n.as_bytes());
            if wide {
                dir.extend_from_slice(&off.to_le_bytes());
            } else {
                dir.extend_from_slice(&(off as u32).to_le_bytes());
            }
        }
        self.out.seek(SeekFrom::Start(HEADER_LEN as u64))?;
        self.out.write_all(&dir)?;
        self.out.seek(SeekFrom::Start(self.pos))?;
        self.out.flush()?;
        Ok(self.out)
    }
}

/// Encodes a whole sequence as one record.
pub fn encode_twobit(seq: &[u8]) -> Result<(TwoBitRecord, Vec<u8>)> {
    let piece = TwoBitPiece::encode(0, seq)?;
    let mut record = TwoBitRecord::new();
    record.push(&piece)?;
    Ok((record, piece.packed))
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::TwoBitReader;
    use omicsweave_core::{Locus, Region, SequenceRead};
    use std::io::Cursor;

    /// What the reader returns: N blocks read as `N` whatever their mask
    /// (as lib2bit reads them), other bases keep their case.
    fn expected(seq: &[u8]) -> Vec<u8> {
        seq.iter()
            .map(|&c| if is_acgt(c) { c } else { b'N' })
            .collect()
    }

    fn write(seqs: &[(&str, &[u8])], piece: usize, version: u32) -> Vec<u8> {
        let names = seqs.iter().map(|(n, _)| n.to_string()).collect();
        let mut w = TwoBitWriter::create(Cursor::new(Vec::new()), names, version).unwrap();
        for (name, seq) in seqs {
            let mut rec = TwoBitRecord::new();
            let mut packed = Vec::new();
            for (i, chunk) in seq.chunks(piece).enumerate() {
                let p = TwoBitPiece::encode((i * piece) as u64, chunk).unwrap();
                rec.push(&p).unwrap();
                packed.extend_from_slice(&p.packed);
            }
            w.add(name, &rec, &mut &packed[..]).unwrap();
        }
        w.finish().unwrap().into_inner()
    }

    #[test]
    fn pieces_round_trip_through_the_reader() {
        let a: Vec<u8> = (0..1003u32)
            .map(|i| b"ACGTacgtNNnnRyAC-"[(i * 7 % 17) as usize])
            .collect();
        let b = b"nnnnACGT".to_vec();
        let c = b"G".to_vec();
        for version in [0, 1] {
            for piece in [4, 8, 400, 2000] {
                let file = write(
                    &[("a", &a[..]), ("b", &b[..]), ("c", &c[..])],
                    piece,
                    version,
                );
                let mut reader = TwoBitReader::open(file, true).unwrap();
                for (name, seq) in [("a", &a), ("b", &b), ("c", &c)] {
                    let mut out = vec![0u8; seq.len()];
                    reader
                        .sequence(
                            &Locus::from(Region::new(name, 0, seq.len() as u32)),
                            &mut out,
                        )
                        .unwrap();
                    assert_eq!(out, expected(seq), "{name} piece {piece} v{version}");
                }
            }
        }
    }

    #[test]
    fn piece_boundaries_do_not_change_the_file() {
        let a: Vec<u8> = (0..999u32)
            .map(|i| b"nnNNacgtACGT"[(i / 3 % 12) as usize])
            .collect();
        let one = write(&[("a", &a[..])], 10_000, 0);
        for piece in [4, 12, 100] {
            assert_eq!(write(&[("a", &a[..])], piece, 0), one);
        }
    }
}
