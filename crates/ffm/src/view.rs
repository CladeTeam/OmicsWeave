//! Borrowed model views and local symbol encoding.

use crate::block::{Block, FLAG_GAPS, FLAG_SITES};
use crate::error::{Result, invalid};
use crate::model::{FactorKind, Motif, Support, table_len};

/// A checked model borrowed from a decoded block or a resident collection.
#[derive(Clone, Copy)]
pub struct MotifView<'a> {
    width: u32,
    alphabet_len: usize,
    support: Option<&'a [u32]>,
    sites: Option<f64>,
    /// Deletion then insertion counts of a gapped profile.
    gaps: Option<&'a [f64]>,
    factor_count: u32,
    is_pfm: bool,
    layout: FactorLayout<'a>,
    values: &'a [f64],
}

/// The gap counts of a gapped profile (see [`crate::Gaps`]).
#[derive(Clone, Copy, Debug)]
pub struct GapsView<'a> {
    values: &'a [f64],
    width: usize,
}

impl<'a> GapsView<'a> {
    /// Per column: sequences in which the column is deleted.
    pub fn deletions(&self) -> &'a [f64] {
        &self.values[..self.width]
    }

    /// Per junction after column `i`: residues inserted there.
    pub fn insertions(&self) -> &'a [f64] {
        &self.values[self.width..]
    }

    /// Deletions then insertions, as stored.
    pub fn values(&self) -> &'a [f64] {
        self.values
    }
}

#[derive(Clone, Copy)]
enum FactorLayout<'a> {
    Stream(&'a [u32]),
    PfmStream(&'a [u32]),
    Pfm(&'a [u32]),
    Indexed {
        factors: &'a [FactorDescriptor],
        positions: &'a [u32],
    },
}

/// Offsets are relative to the motif's values and the shared positions arena.
#[derive(Clone, Copy)]
pub(crate) struct FactorDescriptor {
    pub table_offset: usize,
    pub position_offset: usize,
    pub table_len: u32,
    pub order: u8,
    pub kind: FactorKind,
}

impl<'a> MotifView<'a> {
    /// Reads one record; its ranges and shape were checked when decoded.
    pub(crate) fn from_block(
        block: &'a Block,
        record: usize,
        alphabet_len: usize,
    ) -> MotifView<'a> {
        let meta = &block.records[record];
        let s = meta.starts.structure as usize;
        let v = meta.starts.values as usize;
        let (s_end, v_end) = block
            .records
            .get(record + 1)
            .map_or((block.structure.len(), block.values.len()), |next| {
                (next.starts.structure as usize, next.starts.values as usize)
            });
        let st = &block.structure[s..s_end];
        let width = st[0];
        let support_count = st[1] as usize;
        let (support, rest) = if support_count == alphabet_len {
            (None, &st[2..])
        } else {
            (Some(&st[2..2 + support_count]), &st[2 + support_count..])
        };
        let has_sites = rest[0] & FLAG_SITES != 0;
        let has_gaps = rest[0] & FLAG_GAPS != 0;
        let factor_count = rest[1];
        let values = &block.values[v..v_end];
        let after_sites = &values[usize::from(has_sites)..];
        let gap_len = if has_gaps { 2 * width as usize - 1 } else { 0 };
        MotifView {
            width,
            alphabet_len,
            support,
            sites: has_sites.then(|| values[0]),
            gaps: has_gaps.then(|| &after_sites[..gap_len]),
            factor_count,
            is_pfm: meta.is_pfm,
            layout: if meta.is_pfm {
                FactorLayout::PfmStream(&rest[2..])
            } else {
                FactorLayout::Stream(&rest[2..])
            },
            values: &after_sites[gap_len..],
        }
    }

    pub(crate) fn from_arena(
        arena: &'a crate::collection::MotifArena,
        record: &crate::collection::MotifDescriptor,
        alphabet_len: usize,
    ) -> Self {
        let (factor_count, layout) = match &record.factors {
            crate::collection::FactorRange::Pfm => (
                record.width,
                FactorLayout::Pfm(&arena.pfm_positions[..record.width as usize]),
            ),
            crate::collection::FactorRange::Indexed(factors) => (
                factors.len() as u32,
                FactorLayout::Indexed {
                    factors: &arena.factors[factors.clone()],
                    positions: &arena.positions,
                },
            ),
        };
        Self {
            width: record.width,
            alphabet_len,
            support: (!record.support.is_empty()).then(|| &arena.support[record.support.clone()]),
            sites: record.sites,
            gaps: record.gaps.clone().map(|r| &arena.gap_values[r]),
            factor_count,
            is_pfm: matches!(record.factors, crate::collection::FactorRange::Pfm),
            layout,
            values: &arena.values[record.values.clone()],
        }
    }

    pub fn width(&self) -> u32 {
        self.width
    }

    pub fn alphabet_len(&self) -> usize {
        self.alphabet_len
    }

    /// Global codes modelled individually; `None` when the motif models the
    /// whole alphabet.
    pub fn support(&self) -> Option<&'a [u32]> {
        self.support
    }

    /// Local codes per position.
    pub fn arity(&self) -> usize {
        match self.support {
            None => self.alphabet_len,
            Some(c) => c.len() + 1,
        }
    }

    pub fn sites(&self) -> Option<f64> {
        self.sites
    }

    /// The gap counts when the motif is a gapped profile.
    pub fn gaps(&self) -> Option<GapsView<'a>> {
        self.gaps.map(|values| GapsView {
            values,
            width: self.width as usize,
        })
    }

    pub fn factor_count(&self) -> usize {
        self.factor_count as usize
    }

    pub fn factors(&self) -> Factors<'a> {
        Factors {
            remaining: self.factor_count,
            arity: self.arity(),
            layout: self.layout,
            values: self.values,
        }
    }

    /// O(1) for a resident collection or a PFM; O(i) for other block models.
    pub fn factor(&self, i: usize) -> Option<FactorView<'a>> {
        self.factors().nth(i)
    }

    /// All factor tables back to back, excluding `sites`.
    pub fn table_values(&self) -> &'a [f64] {
        self.values
    }

    /// A position-major matrix when there is one single-position joint
    /// factor per position, in position order. Columns are local codes.
    pub fn pfm(&self) -> Option<PfmView<'a>> {
        self.is_pfm.then(|| PfmView {
            values: self.values,
            arity: self.arity(),
        })
    }

    /// Byte-to-local-code maps. Reverse codes complement the global symbol
    /// before applying this motif's support; reverse the positions separately.
    pub fn encoding(&self, alphabet: &crate::Alphabet) -> Result<MotifEncoding> {
        if alphabet.len() != self.alphabet_len {
            return Err(invalid("FFM encoding", "alphabet size differs from model"));
        }
        let mut local = [crate::INVALID; 256];
        if let Some(support) = self.support {
            local[..alphabet.len()].fill(support.len() as u8);
            for (i, &code) in support.iter().enumerate() {
                local[code as usize] = i as u8;
            }
        } else {
            for (i, code) in local[..alphabet.len()].iter_mut().enumerate() {
                *code = i as u8;
            }
        }
        let mut forward = [crate::INVALID; 256];
        let mut reverse = alphabet.complement().map(|_| [crate::INVALID; 256]);
        for (byte, &global) in alphabet.lut().iter().enumerate() {
            if global == crate::INVALID {
                continue;
            }
            forward[byte] = local[global as usize];
            if let (Some(comp), Some(out)) = (alphabet.complement(), reverse.as_mut()) {
                out[byte] = local[comp[global as usize] as usize];
            }
        }
        Ok(MotifEncoding { forward, reverse })
    }

    pub fn to_motif(&self) -> Motif {
        Motif {
            width: self.width,
            support: match self.support {
                None => Support::All,
                Some(c) => Support::Codes(c.iter().map(|&x| x as u8).collect()),
            },
            sites: self.sites,
            factors: self
                .factors()
                .map(|f| crate::model::Factor {
                    kind: f.kind,
                    positions: f.positions.to_vec(),
                    table: f.table.to_vec(),
                })
                .collect(),
            gaps: self.gaps().map(|g| crate::model::Gaps {
                deletions: g.deletions().to_vec(),
                insertions: g.insertions().to_vec(),
            }),
        }
    }
}

/// One factor of a [`MotifView`].
#[derive(Clone, Copy, Debug)]
pub struct FactorView<'a> {
    pub kind: FactorKind,
    pub positions: &'a [u32],
    pub table: &'a [f64],
}

impl FactorView<'_> {
    /// Writes normalized frequencies using scaled sums, so finite counts near
    /// f64::MAX can be normalized without overflow. Potentials are not normalized.
    pub fn normalize(&self, arity: usize, out: &mut [f64]) -> Result<()> {
        omicsweave_core::check_output("normalized factor", out.len(), self.table.len())?;
        if table_len(arity, self.positions.len())? != self.table.len()
            || self.table.iter().any(|v| !v.is_finite() || *v < 0.0)
        {
            return Err(invalid("FFM factor", "invalid table for normalization"));
        }
        let group = match self.kind {
            FactorKind::Joint => self.table.len(),
            FactorKind::Conditional => arity,
            FactorKind::Potential => {
                return Err(invalid(
                    "FFM factor",
                    "potential has no normalization groups",
                ));
            }
        };
        if self
            .table
            .chunks(group)
            .any(|g| !g.iter().any(|&v| v > 0.0))
        {
            return Err(invalid("FFM factor", "normalization group has zero mass"));
        }
        for (input, output) in self.table.chunks(group).zip(out.chunks_mut(group)) {
            let max = input.iter().copied().fold(0.0f64, f64::max);
            let sum: f64 = input.iter().map(|v| v / max).sum();
            for (&v, o) in input.iter().zip(output) {
                *o = (v / max) / sum;
            }
        }
        Ok(())
    }
}

#[derive(Clone, Copy, Debug)]
pub struct PfmView<'a> {
    values: &'a [f64],
    arity: usize,
}

impl<'a> PfmView<'a> {
    pub fn width(&self) -> usize {
        self.values.len() / self.arity
    }
    pub fn arity(&self) -> usize {
        self.arity
    }
    pub fn values(&self) -> &'a [f64] {
        self.values
    }
    pub fn row(&self, position: usize) -> Option<&'a [f64]> {
        if position >= self.width() {
            return None;
        }
        // SAFETY: position < width and private constructors validate a whole
        // number of arity-sized rows, so both products fit values.len().
        Some(unsafe {
            self.values
                .get_unchecked(position * self.arity..(position + 1) * self.arity)
        })
    }
    pub fn rows(&self) -> std::slice::ChunksExact<'a, f64> {
        self.values.chunks_exact(self.arity)
    }
}

#[derive(Clone, Debug, PartialEq, Eq)]
pub struct MotifEncoding {
    forward: [u8; 256],
    reverse: Option<[u8; 256]>,
}

impl MotifEncoding {
    pub fn forward(&self) -> &[u8; 256] {
        &self.forward
    }
    pub fn reverse(&self) -> Option<&[u8; 256]> {
        self.reverse.as_ref()
    }
}

#[derive(Clone)]
pub struct Factors<'a> {
    remaining: u32,
    arity: usize,
    layout: FactorLayout<'a>,
    values: &'a [f64],
}

impl<'a> Iterator for Factors<'a> {
    type Item = FactorView<'a>;

    fn next(&mut self) -> Option<FactorView<'a>> {
        if self.remaining == 0 {
            return None;
        }
        self.remaining -= 1;
        if let FactorLayout::Pfm(positions) = self.layout {
            // SAFETY: private constructors validate one position and arity
            // cells per remaining PFM factor; each step consumes one row.
            let (position, rest) = unsafe { positions.split_at_unchecked(1) };
            self.layout = FactorLayout::Pfm(rest);
            // SAFETY: each remaining validated PFM row has arity cells.
            let (table, rest) = unsafe { self.values.split_at_unchecked(self.arity) };
            self.values = rest;
            return Some(FactorView {
                kind: FactorKind::Joint,
                positions: position,
                table,
            });
        }
        if let FactorLayout::PfmStream(structure) = self.layout {
            // SAFETY: block validation established each PFM factor's three
            // words (Joint, order=1, position), and one arity-sized table.
            let (position, rest, table, values) = unsafe {
                (
                    structure.get_unchecked(2..3),
                    structure.get_unchecked(3..),
                    self.values.get_unchecked(..self.arity),
                    self.values.get_unchecked(self.arity..),
                )
            };
            self.layout = FactorLayout::PfmStream(rest);
            self.values = values;
            return Some(FactorView {
                kind: FactorKind::Joint,
                positions: position,
                table,
            });
        }
        if let FactorLayout::Indexed { factors, positions } = self.layout {
            // SAFETY: remaining > 0 and private arena construction stores
            // exactly remaining validated factor descriptors in this slice.
            let (f, rest) = unsafe { (factors.get_unchecked(0), factors.get_unchecked(1..)) };
            self.layout = FactorLayout::Indexed {
                factors: rest,
                positions,
            };
            return Some(FactorView {
                kind: f.kind,
                // SAFETY: arena descriptors are made from checked models;
                // positions use global arena offsets, tables motif-local ones.
                positions: unsafe {
                    positions.get_unchecked(f.position_offset..f.position_offset + f.order as usize)
                },
                table: unsafe {
                    self.values
                        .get_unchecked(f.table_offset..f.table_offset + f.table_len as usize)
                },
            });
        }
        let FactorLayout::Stream(structure) = self.layout else {
            unreachable!()
        };
        // SAFETY: block validation consumes every factor's header, positions
        // and table before this private iterator can be constructed.
        let (kind_id, order) = unsafe {
            (
                *structure.get_unchecked(0),
                *structure.get_unchecked(1) as usize,
            )
        };
        let kind = FactorKind::from_id(kind_id).expect("checked kind");
        // SAFETY: the checked order fits this factor's structure slice.
        let (positions, rest) = unsafe {
            (
                structure.get_unchecked(2..2 + order),
                structure.get_unchecked(2 + order..),
            )
        };
        self.layout = FactorLayout::Stream(rest);
        let cells = self.arity.pow(order as u32);
        // SAFETY: table_len validated this power without overflow, and the
        // checked record contains this factor's complete table.
        let (table, rest) = unsafe { self.values.split_at_unchecked(cells) };
        self.values = rest;
        Some(FactorView {
            kind,
            positions,
            table,
        })
    }

    fn nth(&mut self, n: usize) -> Option<Self::Item> {
        if n >= self.remaining as usize {
            self.remaining = 0;
            return None;
        }
        match self.layout {
            FactorLayout::Pfm(positions) => {
                // SAFETY: n < remaining; a PFM owns one position and arity
                // cells per remaining factor. Products fit validated lengths.
                unsafe {
                    self.layout = FactorLayout::Pfm(positions.get_unchecked(n..));
                    self.values = self.values.get_unchecked(n * self.arity..);
                }
            }
            FactorLayout::PfmStream(structure) => {
                // SAFETY: each checked PFM stream factor is three words;
                // n < remaining, so both skips fit the record's slices.
                unsafe {
                    self.layout = FactorLayout::PfmStream(structure.get_unchecked(n * 3..));
                    self.values = self.values.get_unchecked(n * self.arity..);
                }
            }
            FactorLayout::Indexed { factors, positions } => {
                self.layout = FactorLayout::Indexed {
                    // SAFETY: n < remaining and descriptor count is exact.
                    factors: unsafe { factors.get_unchecked(n..) },
                    positions,
                };
            }
            FactorLayout::Stream(mut structure) => {
                for _ in 0..n {
                    // SAFETY: n < remaining and block validation proved each
                    // skipped header, order, table length and position range.
                    unsafe {
                        let order = *structure.get_unchecked(1) as usize;
                        self.values = self.values.get_unchecked(self.arity.pow(order as u32)..);
                        structure = structure.get_unchecked(2 + order..);
                    }
                }
                self.layout = FactorLayout::Stream(structure);
            }
        }
        self.remaining -= n as u32;
        self.next()
    }

    fn size_hint(&self) -> (usize, Option<usize>) {
        (self.remaining as usize, Some(self.remaining as usize))
    }
}

impl ExactSizeIterator for Factors<'_> {}
