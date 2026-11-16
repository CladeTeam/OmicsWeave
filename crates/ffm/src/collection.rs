//! Resident models: shared ownership, contiguous arenas and borrowed views.

use std::io::Write;
use std::ops::Range;
use std::sync::Arc;

use omicsweave_core::ReadAt;

use crate::error::{Result, invalid};
use crate::meta::IdTable;
use crate::view::FactorDescriptor;
use crate::{
    Alphabet, FactorKind, FactorView, FfmReader, FfmWriter, Motif, MotifView, ReaderOptions,
    Schema, WriterOptions,
};

pub(crate) struct MotifDescriptor {
    pub width: u32,
    pub sites: Option<f64>,
    pub gaps: Option<Range<usize>>,
    pub support: Range<usize>,
    pub factors: FactorRange,
    pub values: Range<usize>,
    pub annotation: Range<usize>,
}

pub(crate) enum FactorRange {
    Pfm,
    Indexed(Range<usize>),
}

#[derive(Default)]
pub(crate) struct MotifArena {
    pub motifs: Vec<MotifDescriptor>,
    pub factors: Vec<FactorDescriptor>,
    factor_count: usize,
    pub positions: Vec<u32>,
    pub pfm_positions: Vec<u32>,
    pub values: Vec<f64>,
    pub gap_values: Vec<f64>,
    pub support: Vec<u32>,
    annotation_text: String,
    field_offsets: Vec<usize>,
}

struct MotifParts<'a, F> {
    width: u32,
    support: Option<&'a [u32]>,
    sites: Option<f64>,
    gaps: Option<(&'a [f64], &'a [f64])>,
    factors: F,
}

impl MotifArena {
    pub(crate) fn push_view(
        &mut self,
        model: MotifView<'_>,
        annotation: &str,
        keys: usize,
        limit: usize,
    ) -> Result<()> {
        self.append(
            MotifParts {
                width: model.width(),
                support: model.support(),
                sites: model.sites(),
                gaps: model.gaps().map(|g| (g.deletions(), g.insertions())),
                factors: model.factors(),
            },
            annotation.split('\t').take(keys),
            keys,
            limit,
        )
    }

    fn append<'a>(
        &mut self,
        model: MotifParts<'_, impl ExactSizeIterator<Item = FactorView<'a>> + Clone>,
        annotation: impl Iterator<Item = &'a str> + Clone,
        keys: usize,
        limit: usize,
    ) -> Result<()> {
        let MotifParts {
            width,
            support,
            sites,
            gaps,
            factors,
        } = model;
        let gap_count = gaps.map_or(0, |(d, i)| d.len() + i.len());
        let factor_count = factors.len();
        let mut position_count = 0usize;
        let mut value_count = 0usize;
        let mut is_pfm = factor_count == width as usize;
        for (p, f) in factors.clone().enumerate() {
            is_pfm &= f.kind == FactorKind::Joint && f.positions == [p as u32];
            position_count = position_count
                .checked_add(f.positions.len())
                .ok_or_else(|| invalid("FFM collection", "too many positions"))?;
            value_count = value_count
                .checked_add(f.table.len())
                .ok_or_else(|| invalid("FFM collection", "too many values"))?;
        }
        let annotation_bytes = annotation
            .clone()
            .try_fold(keys, |n, v| n.checked_add(v.len()))
            .ok_or_else(|| invalid("FFM collection", "annotation too long"))?;
        let field_count = if keys == 0 { 0 } else { keys + 1 };
        let support_len = support.map_or(0, <[u32]>::len);
        let descriptor_count = if is_pfm { 0 } else { factor_count };
        let new_pfm_positions = if is_pfm {
            (width as usize).saturating_sub(self.pfm_positions.len())
        } else {
            0
        };
        let position_count = if is_pfm { 0 } else { position_count };
        let sizes = [
            (self.motifs.len(), 1, std::mem::size_of::<MotifDescriptor>()),
            (
                self.factors.len(),
                descriptor_count,
                std::mem::size_of::<FactorDescriptor>(),
            ),
            (self.positions.len(), position_count, 4),
            (self.pfm_positions.len(), new_pfm_positions, 4),
            (self.values.len(), value_count, 8),
            (self.gap_values.len(), gap_count, 8),
            (self.support.len(), support_len, 4),
            (self.annotation_text.len(), annotation_bytes, 1),
            (
                self.field_offsets.len(),
                field_count,
                std::mem::size_of::<usize>(),
            ),
        ];
        let size = sizes
            .into_iter()
            .try_fold(0usize, |total, (n, add, width)| {
                n.checked_add(add)
                    .and_then(|n| n.checked_mul(width))
                    .and_then(|n| total.checked_add(n))
            });
        if size.is_none_or(|n| n > limit) {
            return Err(invalid(
                "FFM collection",
                "resident data exceeds allocation limit",
            ));
        }
        let allocation = |_| invalid("FFM collection", "cannot reserve resident data");
        self.motifs.try_reserve(1).map_err(allocation)?;
        self.factors
            .try_reserve(descriptor_count)
            .map_err(allocation)?;
        self.positions
            .try_reserve(position_count)
            .map_err(allocation)?;
        self.pfm_positions
            .try_reserve(new_pfm_positions)
            .map_err(allocation)?;
        self.values.try_reserve(value_count).map_err(allocation)?;
        self.gap_values.try_reserve(gap_count).map_err(allocation)?;
        self.support.try_reserve(support_len).map_err(allocation)?;
        self.annotation_text
            .try_reserve(annotation_bytes)
            .map_err(allocation)?;
        self.field_offsets
            .try_reserve(field_count)
            .map_err(allocation)?;

        let gaps_range = gaps.map(|(d, i)| {
            let start = self.gap_values.len();
            self.gap_values.extend_from_slice(d);
            self.gap_values.extend_from_slice(i);
            start..self.gap_values.len()
        });
        let support_start = self.support.len();
        self.support.extend_from_slice(support.unwrap_or_default());
        let factor_start = self.factors.len();
        let values_start = self.values.len();
        if is_pfm {
            self.pfm_positions
                .extend(self.pfm_positions.len() as u32..width);
        }
        for f in factors {
            if !is_pfm {
                self.factors.push(FactorDescriptor {
                    table_offset: self.values.len() - values_start,
                    position_offset: self.positions.len(),
                    table_len: f.table.len() as u32,
                    order: f.positions.len() as u8,
                    kind: f.kind,
                });
                self.positions.extend_from_slice(f.positions);
            }
            self.values.extend_from_slice(f.table);
        }
        let annotation_start = self.annotation_text.len();
        for value in annotation {
            self.field_offsets.push(self.annotation_text.len());
            self.annotation_text.push_str(value);
            self.annotation_text.push('\t');
        }
        if keys > 0 {
            self.annotation_text.pop();
            self.annotation_text.push('\0');
            self.field_offsets.push(self.annotation_text.len());
        }
        self.factor_count += factor_count;
        self.motifs.push(MotifDescriptor {
            width,
            sites,
            gaps: gaps_range,
            support: support_start..self.support.len(),
            factors: if is_pfm {
                FactorRange::Pfm
            } else {
                FactorRange::Indexed(factor_start..self.factors.len())
            },
            values: values_start..self.values.len(),
            annotation: annotation_start..self.annotation_text.len(),
        });
        Ok(())
    }

    fn heap_bytes(&self) -> usize {
        self.motifs.capacity() * std::mem::size_of::<MotifDescriptor>()
            + self.factors.capacity() * std::mem::size_of::<FactorDescriptor>()
            + (self.positions.capacity() + self.pfm_positions.capacity() + self.support.capacity())
                * 4
            + (self.values.capacity() + self.gap_values.capacity()) * 8
            + self.annotation_text.capacity()
            + self.field_offsets.capacity() * std::mem::size_of::<usize>()
    }

    fn shrink_to_fit(&mut self) {
        self.motifs.shrink_to_fit();
        self.factors.shrink_to_fit();
        self.positions.shrink_to_fit();
        self.pfm_positions.shrink_to_fit();
        self.values.shrink_to_fit();
        self.gap_values.shrink_to_fit();
        self.support.shrink_to_fit();
        self.annotation_text.shrink_to_fit();
        self.field_offsets.shrink_to_fit();
    }
}

struct MotifSetInner {
    alphabet: Alphabet,
    schema: Schema,
    ids: IdTable,
    arena: MotifArena,
}

/// An immutable, fully decoded collection, independent of the byte source.
/// Cloning shares the arenas; lookups borrow without locks or allocations.
#[derive(Clone)]
pub struct MotifSet {
    inner: Arc<MotifSetInner>,
}

impl MotifSet {
    /// Reads, decodes and validates all models and annotations.
    pub fn load(source: impl ReadAt) -> Result<Self> {
        FfmReader::open(source)?.into_set()
    }

    pub fn load_with_options(source: impl ReadAt, options: ReaderOptions) -> Result<Self> {
        FfmReader::open_with_options(source, options)?.into_set()
    }

    pub(crate) fn from_parts(
        alphabet: Alphabet,
        schema: Schema,
        ids: IdTable,
        mut arena: MotifArena,
    ) -> Self {
        arena.shrink_to_fit();
        Self {
            inner: Arc::new(MotifSetInner {
                alphabet,
                schema,
                ids,
                arena,
            }),
        }
    }

    pub fn len(&self) -> usize {
        self.inner.arena.motifs.len()
    }
    pub fn is_empty(&self) -> bool {
        self.len() == 0
    }
    pub fn alphabet(&self) -> &Alphabet {
        &self.inner.alphabet
    }
    pub fn schema(&self) -> &Schema {
        &self.inner.schema
    }
    pub fn id(&self, i: usize) -> Option<&str> {
        self.inner.ids.get(i)
    }
    pub fn find(&self, id: &str) -> Option<usize> {
        self.inner.ids.find(id)
    }
    pub fn factor_count(&self) -> usize {
        self.inner.arena.factor_count
    }
    /// Table cells; `sites` is held in the motif descriptors separately.
    pub fn value_count(&self) -> usize {
        self.inner.arena.values.len()
    }
    /// Every factor table in motif order, with no `sites` cells or padding.
    pub fn table_values(&self) -> &[f64] {
        &self.inner.arena.values
    }

    pub fn get(&self, i: usize) -> Option<MotifEntry<'_>> {
        self.inner.arena.motifs.get(i)?;
        Some(MotifEntry {
            set: self,
            index: i,
        })
    }

    pub fn get_by_id(&self, id: &str) -> Option<MotifEntry<'_>> {
        self.get(self.find(id)?)
    }

    pub fn motif(&self, i: usize) -> Option<MotifView<'_>> {
        let record = self.inner.arena.motifs.get(i)?;
        Some(MotifView::from_arena(
            &self.inner.arena,
            record,
            self.alphabet().len(),
        ))
    }

    pub fn annotation(&self, i: usize) -> Option<AnnotationView<'_>> {
        let record = self.inner.arena.motifs.get(i)?;
        let keys = self.schema().len();
        let field_offsets = if keys == 0 {
            &[][..]
        } else {
            &self.inner.arena.field_offsets[i * (keys + 1)..(i + 1) * (keys + 1)]
        };
        Some(AnnotationView {
            text: &self.inner.arena.annotation_text,
            row: record.annotation.clone(),
            field_offsets,
            schema: self.schema(),
        })
    }

    pub fn iter(&self) -> impl ExactSizeIterator<Item = MotifEntry<'_>> + DoubleEndedIterator {
        (0..self.len()).map(|index| MotifEntry { set: self, index })
    }

    /// Bytes of owned heap storage, including capacities. Shared data is counted
    /// once per call, not once per clone; reference counts and allocator
    /// bookkeeping are excluded.
    pub fn heap_bytes(&self) -> usize {
        std::mem::size_of::<MotifSetInner>()
            + self.inner.arena.heap_bytes()
            + self.inner.ids.heap_bytes()
            + self.alphabet().heap_bytes()
            + self.schema().heap_bytes()
    }

    /// Writes borrowed views directly, without materializing per-factor Vecs.
    pub fn write_to<W: Write>(&self, out: W, options: WriterOptions) -> Result<W> {
        let mut writer =
            FfmWriter::new(out, self.alphabet().clone(), self.schema().clone(), options)?;
        let mut annotation = Vec::with_capacity(self.schema().len());
        for record in self.iter() {
            annotation.clear();
            annotation.extend(record.annotation().values());
            writer.push_view(record.id(), record.model(), &annotation)?;
        }
        writer.finish()
    }
}

#[derive(Clone, Copy)]
pub struct MotifEntry<'a> {
    set: &'a MotifSet,
    index: usize,
}

impl<'a> MotifEntry<'a> {
    pub fn index(&self) -> usize {
        self.index
    }
    pub fn id(&self) -> &'a str {
        self.set.id(self.index).expect("record index")
    }
    pub fn model(&self) -> MotifView<'a> {
        self.set.motif(self.index).expect("record index")
    }
    pub fn annotation(&self) -> AnnotationView<'a> {
        self.set.annotation(self.index).expect("record index")
    }
}

/// Indexed annotation fields. An empty string denotes missing data.
#[derive(Clone)]
pub struct AnnotationView<'a> {
    text: &'a str,
    row: Range<usize>,
    field_offsets: &'a [usize],
    schema: &'a Schema,
}

impl<'a> AnnotationView<'a> {
    pub fn as_str(&self) -> &'a str {
        if self.schema.is_empty() {
            ""
        } else {
            &self.text[self.row.start..self.row.end - 1]
        }
    }
    pub fn get(&self, i: usize) -> Option<&'a str> {
        if i >= self.schema.len() {
            return None;
        }
        Some(&self.text[self.field_offsets[i]..self.field_offsets[i + 1] - 1])
    }
    pub fn get_by_name(&self, name: &str) -> Option<&'a str> {
        self.get(self.schema.index(name)?)
    }
    pub fn values(&self) -> impl ExactSizeIterator<Item = &'a str> + DoubleEndedIterator + '_ {
        (0..self.schema.len()).map(|i| self.get(i).expect("key index"))
    }
    /// Parses according to the declared schema; missing and out-of-range values
    /// return None. Raw strings, including empty ones, remain available via get.
    pub fn value(&self, i: usize) -> Option<crate::AnnotationValue<'a>> {
        self.schema.value(i, self.get(i)?)
    }
}

/// Builds a collection directly from models, without an encode/decode round trip.
pub struct MotifSetBuilder {
    alphabet: Alphabet,
    schema: Schema,
    ids: IdTable,
    arena: MotifArena,
}

impl MotifSetBuilder {
    pub fn new(alphabet: Alphabet, schema: Schema) -> Self {
        Self {
            alphabet,
            schema,
            ids: IdTable::new(),
            arena: MotifArena::default(),
        }
    }
    pub fn len(&self) -> usize {
        self.ids.len()
    }
    pub fn is_empty(&self) -> bool {
        self.ids.len() == 0
    }
    pub fn alphabet(&self) -> &Alphabet {
        &self.alphabet
    }
    pub fn schema(&self) -> &Schema {
        &self.schema
    }

    /// Adds a whole-alphabet PFM from a flat position-major matrix. Width is
    /// inferred from values.len() / alphabet.len(); there are no per-row Vecs.
    pub fn push_pfm(
        &mut self,
        id: &str,
        values: &[f64],
        sites: Option<f64>,
        annotation: &[&str],
    ) -> Result<()> {
        let arity = self.alphabet.len();
        if !values.len().is_multiple_of(arity) {
            return Err(invalid(
                "FFM PFM",
                "matrix length is not a multiple of alphabet size",
            ));
        }
        let width = u32::try_from(values.len() / arity)
            .map_err(|_| invalid("FFM PFM", "width exceeds u32"))?;
        crate::model::check_header(width, None, sites, arity)?;
        crate::writer::model_size(0, sites, 0, (0..width).map(|_| (1, arity)))?;
        self.check(id, annotation)?;
        let positions: Vec<u32> = (0..width).collect();
        for (position, row) in positions.chunks_exact(1).zip(values.chunks_exact(arity)) {
            crate::model::check_factor(FactorKind::Joint, position, row, width, arity, false)?;
        }
        self.arena.append(
            MotifParts {
                width,
                support: None,
                sites,
                gaps: None,
                factors: positions
                    .chunks_exact(1)
                    .zip(values.chunks_exact(arity))
                    .map(|(positions, table)| FactorView {
                        kind: FactorKind::Joint,
                        positions,
                        table,
                    }),
            },
            annotation.iter().copied(),
            self.schema.len(),
            usize::MAX,
        )?;
        self.ids.push(id)?;
        Ok(())
    }

    pub fn push(&mut self, id: &str, model: &Motif, annotation: &[&str]) -> Result<()> {
        model.validate(self.alphabet.len())?;
        crate::writer::model_size(
            model.support.codes().map_or(0, <[u8]>::len),
            model.sites,
            model
                .gaps
                .as_ref()
                .map_or(0, |g| g.deletions.len() + g.insertions.len()),
            model
                .factors
                .iter()
                .map(|f| (f.positions.len(), f.table.len())),
        )?;
        self.check(id, annotation)?;
        let mut support = [0u32; crate::MAX_SYMBOLS];
        let codes = model.support.codes().map(|c| {
            for (o, &v) in support.iter_mut().zip(c) {
                *o = v as u32;
            }
            &support[..c.len()]
        });
        self.arena.append(
            MotifParts {
                width: model.width,
                support: codes,
                sites: model.sites,
                gaps: model
                    .gaps
                    .as_ref()
                    .map(|g| (g.deletions.as_slice(), g.insertions.as_slice())),
                factors: model.factors.iter().map(|f| FactorView {
                    kind: f.kind,
                    positions: &f.positions,
                    table: &f.table,
                }),
            },
            annotation.iter().copied(),
            self.schema.len(),
            usize::MAX,
        )?;
        self.ids.push(id)?;
        Ok(())
    }

    pub fn push_view(&mut self, id: &str, model: MotifView<'_>, annotation: &[&str]) -> Result<()> {
        if model.alphabet_len() != self.alphabet.len() {
            return Err(invalid(
                "FFM collection",
                "view and builder alphabet sizes differ",
            ));
        }
        self.check(id, annotation)?;
        self.arena.append(
            MotifParts {
                width: model.width(),
                support: model.support(),
                sites: model.sites(),
                gaps: model.gaps().map(|g| (g.deletions(), g.insertions())),
                factors: model.factors(),
            },
            annotation.iter().copied(),
            self.schema.len(),
            usize::MAX,
        )?;
        self.ids.push(id)?;
        Ok(())
    }

    fn check(&self, id: &str, annotation: &[&str]) -> Result<()> {
        self.ids.check_push(id)?;
        if annotation.len() != self.schema.len() {
            return Err(invalid("FFM annotation", "value count differs from schema"));
        }
        for (i, v) in annotation.iter().enumerate() {
            self.schema.check_value(i, v)?;
        }
        Ok(())
    }

    pub fn finish(self) -> MotifSet {
        MotifSet::from_parts(self.alphabet, self.schema, self.ids, self.arena)
    }
}
