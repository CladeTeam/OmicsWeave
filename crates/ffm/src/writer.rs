//! Writes an FFM file: motifs are encoded into compressed blocks as they are
//! pushed; `finish` writes the header, the metadata and then the blocks.

use std::io::Write;

use crate::alphabet::Alphabet;
use crate::block::{BlockBuilder, BlockCodecs, RecordOffsets};
use crate::error::{Result, invalid};
use crate::layout::{DirEntry, HEADER_LEN, Header, MINOR, Span, seal};
use crate::meta::{IdTable, index_columns};
use crate::model::Motif;
use crate::schema::Schema;
use crate::stream::{Codec, Levels, push_bytes, push_numbers};
use crate::view::MotifView;

/// Block sizing, codecs and compression effort.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct WriterOptions {
    /// A block is closed before its value stream would pass this many cells.
    pub block_values: usize,
    /// A block is closed once it holds this many records.
    pub block_records: usize,
    pub structure_codec: Codec,
    pub values_codec: Codec,
    pub annotation_codec: Codec,
    pub index_codec: Codec,
    pub ids_codec: Codec,
    pub levels: Levels,
}

impl Default for WriterOptions {
    fn default() -> Self {
        WriterOptions {
            block_values: 4096,
            block_records: 4096,
            structure_codec: Codec::Pco,
            values_codec: Codec::Pco,
            annotation_codec: Codec::Zstd,
            index_codec: Codec::Pco,
            ids_codec: Codec::Zstd,
            levels: Levels {
                pco: pco::DEFAULT_COMPRESSION_LEVEL,
                zstd: 19,
            },
        }
    }
}

impl WriterOptions {
    fn check(&self) -> Result<()> {
        let numeric = |c: Codec| matches!(c, Codec::Raw | Codec::Pco);
        let bytes = |c: Codec| matches!(c, Codec::Raw | Codec::Zstd);
        if !(numeric(self.structure_codec)
            && numeric(self.values_codec)
            && numeric(self.index_codec)
            && bytes(self.annotation_codec)
            && bytes(self.ids_codec))
        {
            return Err(invalid(
                "FFM writer",
                "numeric streams take raw or pcodec, byte streams raw or zstd",
            ));
        }
        if self.block_values == 0 || self.block_records == 0 {
            return Err(invalid("FFM writer", "block limits must be positive"));
        }
        if [self.structure_codec, self.values_codec, self.index_codec].contains(&Codec::Pco)
            && self.levels.pco > 12
        {
            return Err(invalid(
                "FFM writer",
                "pcodec compression level must be in 0..=12",
            ));
        }
        if [self.annotation_codec, self.ids_codec].contains(&Codec::Zstd)
            && !zstd::compression_level_range().contains(&self.levels.zstd)
        {
            return Err(invalid(
                "FFM writer",
                "zstd compression level outside supported range",
            ));
        }
        Ok(())
    }
}

/// Writes one FFM file. Compressed blocks are kept in memory until
/// [`FfmWriter::finish`], because the metadata ahead of them depends on every
/// motif.
pub struct FfmWriter<W: Write> {
    out: W,
    alphabet: Alphabet,
    schema: Schema,
    opts: WriterOptions,
    builder: BlockBuilder,
    blocks: Vec<u8>,
    directory: Vec<(u32, u32)>,
    ids: IdTable,
    motif_blocks: Vec<u32>,
    motif_starts: Vec<RecordOffsets>,
}

impl<W: Write> FfmWriter<W> {
    pub fn new(out: W, alphabet: Alphabet, schema: Schema, opts: WriterOptions) -> Result<Self> {
        opts.check()?;
        Ok(FfmWriter {
            out,
            alphabet,
            schema,
            opts,
            builder: BlockBuilder::default(),
            blocks: Vec::new(),
            directory: Vec::new(),
            ids: IdTable::new(),
            motif_blocks: Vec::new(),
            motif_starts: Vec::new(),
        })
    }

    pub fn alphabet(&self) -> &Alphabet {
        &self.alphabet
    }

    pub fn schema(&self) -> &Schema {
        &self.schema
    }

    /// Motifs pushed so far.
    pub fn len(&self) -> usize {
        self.ids.len()
    }

    pub fn is_empty(&self) -> bool {
        self.ids.len() == 0
    }

    /// Adds a motif with its ID and one value per schema key (empty strings
    /// for missing values).
    pub fn push(&mut self, id: &str, motif: &Motif, annotation: &[&str]) -> Result<()> {
        let n = self.alphabet.len();
        motif.validate(n)?;
        let size = model_size(
            motif.support.codes().map_or(0, <[u8]>::len),
            motif.sites,
            motif
                .gaps
                .as_ref()
                .map_or(0, |g| g.deletions.len() + g.insertions.len()),
            motif
                .factors
                .iter()
                .map(|f| (f.positions.len(), f.table.len())),
        )?;
        self.prepare(id, size, annotation)?;
        let annotation = (!self.schema.is_empty()).then_some(annotation);
        let starts = self.builder.push(motif, n, annotation);
        self.record(starts);
        Ok(())
    }

    /// Copies a checked, borrowed model directly into the block streams.
    /// Codes must use this writer's alphabet in the same order.
    pub fn push_view(&mut self, id: &str, motif: MotifView<'_>, annotation: &[&str]) -> Result<()> {
        if motif.alphabet_len() != self.alphabet.len() {
            return Err(invalid(
                "FFM writer",
                "view and writer alphabet sizes differ",
            ));
        }
        let size = model_size(
            motif.support().map_or(0, <[u32]>::len),
            motif.sites(),
            motif.gaps().map_or(0, |g| g.values().len()),
            motif.factors().map(|f| (f.positions.len(), f.table.len())),
        )?;
        self.prepare(id, size, annotation)?;
        let annotation = (!self.schema.is_empty()).then_some(annotation);
        let starts = self.builder.push_view(motif, annotation);
        self.record(starts);
        Ok(())
    }

    fn prepare(&mut self, id: &str, size: (usize, usize), annotation: &[&str]) -> Result<()> {
        if annotation.len() != self.schema.len() {
            return Err(invalid(
                "FFM annotation",
                format!(
                    "{} values for {} schema keys",
                    annotation.len(),
                    self.schema.len()
                ),
            ));
        }
        for (i, v) in annotation.iter().enumerate() {
            self.schema.check_value(i, v)?;
        }
        self.ids.check_push(id)?;
        let (structure, cells) = size;
        let ann_bytes = annotation
            .iter()
            .try_fold(annotation.len(), |len, value| len.checked_add(value.len()))
            .filter(|&n| n <= u32::MAX as usize)
            .ok_or_else(|| invalid("FFM writer", "annotation exceeds stream limit"))?;
        let b = &self.builder;
        if b.records > 0
            && (b.records as usize >= self.opts.block_records
                || b.records == u32::MAX
                || b.values.len().saturating_add(cells) > self.opts.block_values
                || b.values.len().saturating_add(cells) > u32::MAX as usize
                || b.structure.len().saturating_add(structure) > u32::MAX as usize
                || b.annotation.len().saturating_add(ann_bytes) > u32::MAX as usize)
        {
            self.flush()?;
        }
        self.ids.push(id)?;
        Ok(())
    }

    fn record(&mut self, starts: RecordOffsets) {
        self.motif_blocks.push(self.directory.len() as u32);
        self.motif_starts.push(starts);
    }

    fn flush(&mut self) -> Result<()> {
        if self.builder.records == 0 {
            return Ok(());
        }
        let records = self.builder.records;
        let codecs = BlockCodecs {
            structure: self.opts.structure_codec,
            values: self.opts.values_codec,
            annotation: self.opts.annotation_codec,
        };
        let bytes = self.builder.encode(&codecs, self.opts.levels)?;
        let len = u32::try_from(bytes.len())
            .map_err(|_| invalid("FFM writer", "block exceeds 2^32 - 1 bytes"))?;
        self.blocks.extend_from_slice(&bytes);
        self.directory.push((len, records));
        self.builder.clear();
        Ok(())
    }

    /// Writes the file and returns the sink.
    pub fn finish(mut self) -> Result<W> {
        self.flush()?;
        let levels = self.opts.levels;

        let mut alphabet = Vec::new();
        self.alphabet.write_to(&mut alphabet);
        seal(&mut alphabet);
        let mut schema = Vec::new();
        self.schema.write_to(&mut schema);
        seal(&mut schema);
        let mut ids = Vec::new();
        push_bytes(&mut ids, &self.ids.stored(), self.opts.ids_codec, levels)?;
        let mut index = Vec::new();
        let columns = index_columns(&self.motif_blocks, &self.motif_starts);
        push_numbers(&mut index, &columns, self.opts.index_codec, levels)?;
        let directory_len = self.directory.len() * crate::layout::DIR_ENTRY_LEN + 8;

        let mut at = HEADER_LEN as u64;
        let mut span = |len: usize| {
            let s = Span {
                offset: at,
                len: len as u64,
            };
            at += len as u64;
            s
        };
        let header = Header {
            minor: MINOR,
            motif_count: self.ids.len() as u64,
            block_count: self.directory.len() as u64,
            alphabet: span(alphabet.len()),
            schema: span(schema.len()),
            ids: span(ids.len()),
            directory: span(directory_len),
            index: span(index.len()),
        };
        let mut directory = Vec::with_capacity(directory_len);
        let mut offset = at;
        for &(length, record_count) in &self.directory {
            DirEntry {
                offset,
                length,
                record_count,
            }
            .write_to(&mut directory);
            offset += length as u64;
        }
        seal(&mut directory);

        for part in [
            &header.encode()[..],
            &alphabet,
            &schema,
            &ids,
            &directory,
            &index,
            &self.blocks,
        ] {
            self.out.write_all(part)?;
        }
        self.out.flush()?;
        Ok(self.out)
    }
}

pub(crate) fn model_size(
    support_len: usize,
    sites: Option<f64>,
    gap_values: usize,
    mut factors: impl Iterator<Item = (usize, usize)>,
) -> Result<(usize, usize)> {
    factors.try_fold(
        (
            4usize + support_len,
            usize::from(sites.is_some()) + gap_values,
        ),
        |(s, v), (order, cells)| {
            let s = s.checked_add(2 + order).filter(|&n| n <= u32::MAX as usize);
            let v = v.checked_add(cells).filter(|&n| n <= u32::MAX as usize);
            s.zip(v)
                .ok_or_else(|| invalid("FFM writer", "model exceeds stream limits"))
        },
    )
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::FfmReader;

    #[test]
    fn compression_failure_preserves_records_and_does_not_claim_id() {
        let opts = WriterOptions {
            block_records: 1,
            ..WriterOptions::default()
        };
        let mut writer =
            FfmWriter::new(Vec::new(), Alphabet::dna(), Schema::default(), opts).unwrap();
        let model = Motif::pfm(&[vec![1.0; 4]], None);
        writer.push("first", &model, &[]).unwrap();
        writer.opts.levels.pco = 13;
        assert!(writer.push("second", &model, &[]).is_err());
        assert_eq!(writer.len(), 1);
        assert_eq!(writer.builder.records, 1);
        writer.opts.levels.pco = opts.levels.pco;
        writer.push("second", &model, &[]).unwrap();
        let reader = FfmReader::open(writer.finish().unwrap()).unwrap();
        assert_eq!(reader.len(), 2);
        for i in 0..2 {
            assert_eq!(reader.motif(i).unwrap().view().to_motif(), model);
        }
    }

    #[test]
    fn sizes_are_checked_without_truncation() {
        assert!(model_size(0, None, 0, [(1, u32::MAX as usize), (1, 1)].into_iter()).is_err());
        assert!(model_size(0, None, 0, [(u32::MAX as usize, 1)].into_iter()).is_err());
        let mut opts = WriterOptions::default();
        opts.levels.pco = 13;
        assert!(FfmWriter::new(Vec::new(), Alphabet::dna(), Schema::default(), opts).is_err());
    }
}
