//! Blocks: consecutive records stored as a structure stream (`u32`), a value
//! stream (`f64`) and an annotation stream (bytes).

use std::sync::OnceLock;

use crate::view::{FactorView, MotifView};

use crate::error::{Error, Result, invalid};
use crate::layout::Cursor;
use crate::model::{FactorKind, Motif, check_factor, check_header, check_support, table_len};
use crate::schema::Schema;
use crate::stream::{Codec, Levels, Stream, next_stream, push_bytes, push_numbers};

pub(crate) const FLAG_SITES: u32 = 1;
pub(crate) const FLAG_GAPS: u32 = 2;

/// Where a record starts in each stream of its block.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub(crate) struct RecordOffsets {
    pub structure: u32,
    pub values: u32,
    pub annotation: u32,
}

/// Records of the block being written.
#[derive(Default)]
pub(crate) struct BlockBuilder {
    pub structure: Vec<u32>,
    pub values: Vec<f64>,
    pub annotation: Vec<u8>,
    pub records: u32,
}

impl BlockBuilder {
    /// Appends a validated record and returns where it starts.
    pub fn push(
        &mut self,
        motif: &Motif,
        alphabet_len: usize,
        annotation: Option<&[&str]>,
    ) -> RecordOffsets {
        self.push_parts(
            motif.width,
            motif.support.codes().map_or(alphabet_len, <[u8]>::len),
            motif
                .support
                .codes()
                .into_iter()
                .flatten()
                .map(|&c| c as u32),
            motif.sites,
            motif
                .gaps
                .as_ref()
                .map(|g| (g.deletions.as_slice(), g.insertions.as_slice())),
            motif.factors.iter().map(|f| FactorView {
                kind: f.kind,
                positions: &f.positions,
                table: &f.table,
            }),
            annotation,
        )
    }

    pub fn push_view(
        &mut self,
        motif: MotifView<'_>,
        annotation: Option<&[&str]>,
    ) -> RecordOffsets {
        self.push_parts(
            motif.width(),
            motif.support().map_or(motif.alphabet_len(), <[u32]>::len),
            motif.support().into_iter().flatten().copied(),
            motif.sites(),
            motif.gaps().map(|g| (g.deletions(), g.insertions())),
            motif.factors(),
            annotation,
        )
    }

    #[allow(clippy::too_many_arguments)]
    fn push_parts<'a>(
        &mut self,
        width: u32,
        support_count: usize,
        support: impl Iterator<Item = u32>,
        sites: Option<f64>,
        gaps: Option<(&[f64], &[f64])>,
        factors: impl ExactSizeIterator<Item = FactorView<'a>>,
        annotation: Option<&[&str]>,
    ) -> RecordOffsets {
        let starts = RecordOffsets {
            structure: self.structure.len() as u32,
            values: self.values.len() as u32,
            annotation: self.annotation.len() as u32,
        };
        let s = &mut self.structure;
        s.push(width);
        s.push(support_count as u32);
        s.extend(support);
        let mut flags = 0;
        if sites.is_some() {
            flags |= FLAG_SITES;
        }
        if gaps.is_some() {
            flags |= FLAG_GAPS;
        }
        s.push(flags);
        s.push(factors.len() as u32);
        if let Some(sites) = sites {
            self.values.push(sites);
        }
        if let Some((deletions, insertions)) = gaps {
            self.values.extend_from_slice(deletions);
            self.values.extend_from_slice(insertions);
        }
        for f in factors {
            s.push(f.kind.id());
            s.push(f.positions.len() as u32);
            s.extend_from_slice(f.positions);
            self.values.extend_from_slice(f.table);
        }
        if let Some(values) = annotation {
            for (i, v) in values.iter().enumerate() {
                if i > 0 {
                    self.annotation.push(b'\t');
                }
                self.annotation.extend_from_slice(v.as_bytes());
            }
            self.annotation.push(0);
        }
        self.records += 1;
        starts
    }

    pub fn clear(&mut self) {
        self.structure.clear();
        self.values.clear();
        self.annotation.clear();
        self.records = 0;
    }

    /// Encodes without changing the builder, so compression failure is retryable.
    pub fn encode(&self, codecs: &BlockCodecs, levels: Levels) -> Result<Vec<u8>> {
        let mut out = Vec::new();
        out.extend_from_slice(&self.records.to_le_bytes());
        out.extend_from_slice(&0u32.to_le_bytes());
        push_numbers(&mut out, &self.structure, codecs.structure, levels)?;
        push_numbers(&mut out, &self.values, codecs.values, levels)?;
        push_bytes(&mut out, &self.annotation, codecs.annotation, levels)?;
        Ok(out)
    }
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) struct BlockCodecs {
    pub structure: Codec,
    pub values: Codec,
    pub annotation: Codec,
}

/// Checked stream starts and shape of one record.
#[derive(Clone, Copy)]
pub(crate) struct BlockRecord {
    pub starts: RecordOffsets,
    pub is_pfm: bool,
}

/// A decoded block. Structure and values are decoded and checked when the
/// block is read; the annotation stays compressed until first asked for.
pub(crate) struct Block {
    pub structure: Vec<u32>,
    pub values: Vec<f64>,
    annotation_codec: Codec,
    annotation_len: usize,
    annotation_stored: Vec<u8>,
    pub records: Vec<BlockRecord>,
    annotation: OnceLock<Box<str>>,
    max_stream_bytes: usize,
}

impl Block {
    /// Decodes a block read whole, checking every record against the index
    /// rows `starts` of its motifs.
    pub fn decode(
        bytes: &[u8],
        starts: &[RecordOffsets],
        alphabet_len: usize,
        max_stream_bytes: usize,
    ) -> Result<Block> {
        const WHAT: &str = "FFM block";
        let mut r = Cursor::new(bytes, WHAT);
        let records = r.u32()? as usize;
        r.zeros(4)?;
        if records != starts.len() {
            return Err(invalid(WHAT, "record count differs from the directory"));
        }
        let structure = next_stream(&mut r, WHAT)?.numbers::<u32>(WHAT, max_stream_bytes)?;
        let values = next_stream(&mut r, WHAT)?.numbers::<f64>(WHAT, max_stream_bytes)?;
        let Stream {
            codec: annotation_codec,
            decoded_len: annotation_len,
            data,
        } = next_stream(&mut r, WHAT)?;
        if annotation_codec == Codec::Pco {
            return Err(invalid(WHAT, "pcodec annotation stream"));
        }
        if annotation_len > max_stream_bytes {
            return Err(invalid(WHAT, "annotation size exceeds allocation limit"));
        }
        let annotation_stored = data.to_vec();
        r.finish()?;

        let mut s = 0usize;
        let mut v = 0usize;
        let mut records = crate::options::reserve(starts.len(), max_stream_bytes, WHAT)?;
        for (i, st) in starts.iter().enumerate() {
            if st.structure as usize != s || st.values as usize != v {
                return Err(invalid(
                    WHAT,
                    format!("record {i} does not start where the index says"),
                ));
            }
            let (next_s, next_v, is_pfm) = check_record(&structure, &values, s, v, alphabet_len)?;
            records.push(BlockRecord {
                starts: *st,
                is_pfm,
            });
            (s, v) = (next_s, next_v);
        }
        if s != structure.len() || v != values.len() {
            return Err(invalid(WHAT, "streams hold data past the last record"));
        }
        Ok(Block {
            structure,
            values,
            annotation_codec,
            annotation_len,
            annotation_stored,
            records,
            annotation: OnceLock::new(),
            max_stream_bytes,
        })
    }

    /// Heap bytes the block holds once its annotation is decoded.
    pub fn heap_bytes(&self) -> usize {
        [
            self.structure.capacity().saturating_mul(4),
            self.values.capacity().saturating_mul(8),
            self.annotation_stored.capacity(),
            self.records
                .capacity()
                .saturating_mul(std::mem::size_of::<BlockRecord>()),
            self.annotation_len,
        ]
        .into_iter()
        .fold(0, usize::saturating_add)
    }

    /// The decoded annotation stream, decompressed and checked on first use.
    pub fn annotation(&self, schema: &Schema) -> Result<&str> {
        if let Some(bytes) = self.annotation.get() {
            return Ok(bytes);
        }
        let stream = Stream {
            codec: self.annotation_codec,
            decoded_len: self.annotation_len,
            data: &self.annotation_stored,
        };
        let bytes = stream.bytes("FFM annotation", self.max_stream_bytes)?;
        check_annotation(&bytes, &self.records, schema)?;
        let text = String::from_utf8(bytes).expect("annotation checked as UTF-8");
        // Admission reserves annotation_len before lazy decode. Retaining
        // only a successful boxed slice keeps this immutable payload within
        // that charge, including concurrent first readers and error paths.
        let _ = self.annotation.set(text.into_boxed_str());
        Ok(self.annotation.get().expect("annotation published"))
    }

    /// The `\0`-terminated annotation of record `j`, without the terminator.
    pub fn annotation_of(&self, schema: &Schema, j: usize) -> Result<&str> {
        let text = self.annotation(schema)?;
        if schema.is_empty() {
            return Ok("");
        }
        let start = self.records[j].starts.annotation as usize;
        let end = match self.records.get(j + 1) {
            Some(e) => e.starts.annotation as usize,
            None => text.len(),
        };
        Ok(&text[start..end - 1])
    }
}

fn check_annotation(bytes: &[u8], records: &[BlockRecord], schema: &Schema) -> Result<()> {
    const WHAT: &str = "FFM annotation";
    if schema.is_empty() {
        if !bytes.is_empty() || records.iter().any(|r| r.starts.annotation != 0) {
            return Err(invalid(WHAT, "annotation without schema keys"));
        }
        return Ok(());
    }
    let text = std::str::from_utf8(bytes).map_err(|_| invalid(WHAT, "not UTF-8"))?;
    let mut pos = 0usize;
    for (j, record) in records.iter().enumerate() {
        if record.starts.annotation as usize != pos {
            return Err(invalid(
                WHAT,
                format!("record {j} does not start where the index says"),
            ));
        }
        let Some(len) = text[pos..].find('\0') else {
            return Err(invalid(WHAT, format!("record {j} is not NUL-terminated")));
        };
        let record = &text[pos..pos + len];
        let mut count = 0;
        for (i, value) in record.split('\t').enumerate() {
            if i >= schema.len() {
                return Err(invalid(WHAT, format!("record {j} has too many values")));
            }
            schema.check_value(i, value)?;
            count += 1;
        }
        if count != schema.len() {
            return Err(invalid(WHAT, format!("record {j} has too few values")));
        }
        pos += len + 1;
    }
    if pos != text.len() {
        return Err(invalid(WHAT, "data past the last record"));
    }
    Ok(())
}

/// Checks the record starting at `s` / `v` and returns where the next one
/// starts.
fn check_record(
    structure: &[u32],
    values: &[f64],
    s: usize,
    v: usize,
    alphabet_len: usize,
) -> Result<(usize, usize, bool)> {
    let mut r = StructureCursor {
        data: structure,
        pos: s,
    };
    let width = r.next()?;
    let support_count = r.next()? as usize;
    let support = if support_count == alphabet_len {
        None
    } else {
        let codes = r.take(support_count)?;
        if codes.iter().any(|&c| c >= alphabet_len as u32) {
            return Err(invalid("FFM motif", "support code outside the alphabet"));
        }
        Some(codes)
    };
    let flags = r.next()?;
    if flags & !(FLAG_SITES | FLAG_GAPS) != 0 {
        return Err(invalid("FFM motif", format!("unknown flags {flags:#x}")));
    }
    let mut vals = ValueCursor {
        data: values,
        pos: v,
    };
    let sites = if flags & FLAG_SITES != 0 {
        Some(vals.take(1)?[0])
    } else {
        None
    };
    check_header(width, None, sites, alphabet_len)?;
    let gaps = if flags & FLAG_GAPS != 0 {
        let w = width as usize;
        Some((vals.take(w)?, vals.take(w - 1)?))
    } else {
        None
    };
    if let Some(codes) = support {
        check_support(codes, alphabet_len)?;
    }
    let ar = support.map_or(alphabet_len, |codes| codes.len() + 1);
    let factor_count = r.next()?;
    if factor_count == 0 {
        return Err(invalid("FFM motif", "no factors"));
    }
    let mut shapes = Vec::new();
    let mut is_pfm = factor_count == width;
    for i in 0..factor_count as usize {
        let kind_id = r.next()?;
        let kind = FactorKind::from_id(kind_id)
            .ok_or_else(|| invalid("FFM factor", format!("unknown kind {kind_id}")))?;
        let order = r.next()? as usize;
        if order == 0 || order > crate::model::MAX_ORDER {
            return Err(invalid("FFM factor", format!("order {order}")));
        }
        let positions = r.take(order)?;
        let table = vals.take(table_len(ar, order)?)?;
        let deleted = gaps.is_some_and(|(d, _)| d.get(i).is_some_and(|&v| v > 0.0));
        check_factor(kind, positions, table, width, ar, deleted)?;
        is_pfm &= kind == FactorKind::Joint && positions == [i as u32];
        if gaps.is_some() {
            shapes.push((kind, positions));
        }
    }
    if let Some((deletions, insertions)) = gaps {
        crate::model::check_gaps(width, deletions, insertions, shapes.into_iter())?;
    }
    Ok((r.pos, vals.pos, is_pfm))
}

struct StructureCursor<'a> {
    data: &'a [u32],
    pos: usize,
}

impl<'a> StructureCursor<'a> {
    fn next(&mut self) -> Result<u32> {
        Ok(self.take(1)?[0])
    }

    fn take(&mut self, n: usize) -> Result<&'a [u32]> {
        let end = self.pos.checked_add(n).filter(|&e| e <= self.data.len());
        let Some(end) = end else {
            return Err(Error::Truncated {
                what: "FFM structure stream",
                need: (self.pos as u64).saturating_add(n as u64).saturating_mul(4),
                have: (self.data.len() as u64).saturating_mul(4),
            });
        };
        let s = &self.data[self.pos..end];
        self.pos = end;
        Ok(s)
    }
}

struct ValueCursor<'a> {
    data: &'a [f64],
    pos: usize,
}

impl<'a> ValueCursor<'a> {
    fn take(&mut self, n: usize) -> Result<&'a [f64]> {
        let end = self.pos.checked_add(n).filter(|&e| e <= self.data.len());
        let Some(end) = end else {
            return Err(Error::Truncated {
                what: "FFM value stream",
                need: (self.pos as u64).saturating_add(n as u64).saturating_mul(8),
                have: (self.data.len() as u64).saturating_mul(8),
            });
        };
        let s = &self.data[self.pos..end];
        self.pos = end;
        Ok(s)
    }
}
