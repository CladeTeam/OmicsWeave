//! [`HicBuilder`]: a `.hic` file from contacts, binned at every resolution
//! as the Juicer tools' `pre` bins them. No expected or normalization
//! vectors are computed.

use std::collections::HashMap;
use std::io::{Seek, Write};

use omicsweave_core::exec::Exec;
use omicsweave_core::{ContactUnit, Error, Result};

use crate::block::Record;
use crate::sidecar::HicSidecar;
use crate::writer::{BlockLayout, HicHeader, HicWriter};

/// Records of one map, merged (sorted by (x, y), cells unique) whenever the
/// unmerged tail grows as large as the merged part.
#[derive(Default)]
struct Accumulator {
    records: Vec<Record>,
    merged: usize,
}

impl Accumulator {
    const MIN_TAIL: usize = 1 << 16;

    fn push(&mut self, record: Record) {
        self.records.push(record);
        if self.records.len() - self.merged >= self.merged.max(Self::MIN_TAIL) {
            self.merge();
        }
    }

    fn merge(&mut self) {
        self.records.sort_by_key(|r| (r.x, r.y));
        self.records.dedup_by(|b, a| {
            let same = a.x == b.x && a.y == b.y;
            if same {
                a.value += b.value;
            }
            same
        });
        self.merged = self.records.len();
    }
}

/// How [`HicBuilder::finish`] compresses and cuts maps into blocks.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct BuildOptions {
    /// zlib level of blocks.
    pub level: u32,
    pub blocks: BlockLayout,
}

impl Default for BuildOptions {
    /// Level 6, as the Juicer tools default to, and the default
    /// [`BlockLayout`].
    fn default() -> BuildOptions {
        BuildOptions {
            level: 6,
            blocks: BlockLayout::default(),
        }
    }
}

/// Collects contacts in memory, then writes the file.
pub struct HicBuilder {
    header: HicHeader,
    // Every contact populates every declared resolution, so a pair owns
    // one contiguous list in header order (BP, then FRAG).
    maps: HashMap<(u32, u32), Vec<Accumulator>>,
}

impl HicBuilder {
    pub fn new(header: HicHeader) -> Result<HicBuilder> {
        header.check()?;
        Ok(HicBuilder {
            header,
            maps: HashMap::new(),
        })
    }

    pub fn header(&self) -> &HicHeader {
        &self.header
    }

    /// Adds `count` contacts between base `pos1` of `chr1` and base `pos2`
    /// of `chr2` (chromosome indexes) at every resolution. Positions past a
    /// chromosome's end count at its end. Fragment bins come from the
    /// header's sites as the Juicer tools assign them: fragment `i` runs
    /// from site `i - 1` (inclusive) to site `i`.
    pub fn add(&mut self, chr1: u32, pos1: u32, chr2: u32, pos2: u32, count: f32) -> Result<()> {
        let n = self.header.chroms.len();
        if chr1 as usize >= n || chr2 as usize >= n {
            return Err(Error::invalid("contact", "chromosome index out of range"));
        }
        if !count.is_finite() {
            return Err(Error::invalid("contact", "non-finite count"));
        }
        let (c1, p1, c2, p2) = if chr1 <= chr2 {
            (chr1, pos1, chr2, pos2)
        } else {
            (chr2, pos2, chr1, pos1)
        };
        let clamp =
            |c: u32, p: u32| p.min(self.header.chroms[c as usize].len.min(u32::MAX as u64) as u32);
        let (p1, p2) = (clamp(c1, p1), clamp(c2, p2));
        let resolutions = self.header.bp_resolutions.len() + self.header.frag_resolutions.len();
        if resolutions == 0 {
            return Ok(());
        }
        let maps = self.maps.entry((c1, c2)).or_insert_with(|| {
            std::iter::repeat_with(Accumulator::default)
                .take(resolutions)
                .collect()
        });
        let mut maps = maps.iter_mut();
        let intra = c1 == c2;
        let mut put = |a: u32, b: u32, res: u32| {
            let (x, y) = (a / res, b / res);
            let (x, y) = if intra && x > y { (y, x) } else { (x, y) };
            // The list has exactly one accumulator per header resolution.
            maps.next().unwrap().push(Record { x, y, value: count });
        };
        for &res in &self.header.bp_resolutions {
            put(p1, p2, res);
        }
        if !self.header.frag_resolutions.is_empty() {
            let fragment =
                |c: u32, p: u32| self.header.sites[c as usize].partition_point(|&s| s <= p) as u32;
            let (f1, f2) = (fragment(c1, p1), fragment(c2, p2));
            for &res in &self.header.frag_resolutions {
                put(f1, f2, res);
            }
        }
        Ok(())
    }

    /// Writes the file: every map, resolution by resolution.
    pub fn finish<W: Write + Seek>(self, out: W, options: &BuildOptions) -> Result<W> {
        Ok(self.finish_with_sidecar(out, options)?.0)
    }

    /// [`HicBuilder::finish`], also returning the sidecar of the file
    /// written ([`HicWriter::finish_with_sidecar`]).
    pub fn finish_with_sidecar<W: Write + Seek>(
        self,
        out: W,
        options: &BuildOptions,
    ) -> Result<(W, HicSidecar)> {
        self.finish_exec(out, options, Exec::Inline)
    }

    /// [`HicBuilder::finish_with_sidecar`], encoding and compressing blocks
    /// on `exec` ([`HicWriter::create_exec`]); the file and sidecar are the
    /// same for every executor.
    pub fn finish_exec<W: Write + Seek>(
        self,
        out: W,
        options: &BuildOptions,
        exec: Exec<'_>,
    ) -> Result<(W, HicSidecar)> {
        let HicBuilder { header, maps } = self;
        let resolutions: Vec<_> = [ContactUnit::Bp, ContactUnit::Frag]
            .into_iter()
            .flat_map(|unit| header.resolutions(unit).iter().map(move |&res| (unit, res)))
            .collect();
        let mut writer = HicWriter::create_exec(out, header, options.level, options.blocks, exec)?;
        let mut pairs: Vec<_> = maps.into_iter().collect();
        // All pairs have the same resolution list, so sort their shared
        // order once and consume each accumulator in resolution order.
        pairs.sort_unstable_by_key(|&(pair, _)| pair);
        for (index, (unit, res)) in resolutions.into_iter().enumerate() {
            for &mut ((c1, c2), ref mut maps) in &mut pairs {
                let mut acc = std::mem::take(&mut maps[index]);
                acc.merge();
                writer.write_matrix(c1, c2, unit, res, &acc.records)?;
            }
        }
        writer.finish_with_sidecar()
    }
}
