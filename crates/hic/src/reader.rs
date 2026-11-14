//! [`HicReader`]: standard [`ContactRead`] access to a `.hic` file over a
//! caller-provided byte source.

use std::collections::HashMap;
use std::ops::Range;

use omicsweave_core::{
    ContactEntry, ContactQuery, ContactRead, ContactUnit, Error, ReadAt, Result, check_output,
};

use crate::block::{Inflater, decode_block};
use crate::cursor::Cursor;
use crate::matrix::{BlockEntry, ZoomData};
use crate::meta::{Depth, Footer, Header};

type ZoomKey = (u32, u32, ContactUnit, u32);

/// A `.hic` file (v8 or v9). Opening reads the header and the footer's
/// master index; each chromosome pair's block index and each fragment site
/// list is read the first time a query needs it and kept. A query reads and
/// decodes each block it overlaps once. Values are the stored contact
/// counts; stored normalization and expected vectors are not read.
pub struct HicReader<S> {
    source: S,
    header: Header,
    footer: Footer,
    depth: Depth,
    chrom_ids: HashMap<String, u32>,
    zooms: HashMap<ZoomKey, Option<ZoomData>>,
    sites: HashMap<u32, Vec<u32>>,
    inflater: Option<Inflater>,
    raw: Vec<u8>,
    block: Vec<u8>,
    blocks: Vec<BlockEntry>,
    entries: Vec<ContactEntry>,
    entries_ready: bool,
}

/// A query resolved against the file: which stored matrix to read, which
/// bins of its x (chr1) and y (chr2) axes, and how cells map onto the tile.
struct Plan {
    chr1: u32,
    chr2: u32,
    x: Range<u32>,
    y: Range<u32>,
    /// The tile's rows run along y (rows lie on the higher-index chromosome).
    transposed: bool,
}

fn missing(what: &str, detail: String) -> Error {
    Error::invalid("hic query", format!("no {what} {detail}"))
}

impl<S: ReadAt> HicReader<S> {
    pub fn open(source: S) -> Result<HicReader<S>> {
        let mut cursor = Cursor::new(&source, 0);
        let header = Header::parse(&mut cursor)?;
        let footer = Footer::parse(&mut cursor, &header, source.len())?;
        let mut chrom_ids = HashMap::new();
        for (i, chrom) in header.chroms.iter().enumerate() {
            if chrom_ids.insert(chrom.name.clone(), i as u32).is_some() {
                return Err(Error::invalid("hic header", "duplicate chromosome name"));
            }
        }
        Ok(HicReader {
            source,
            depth: header.depth(),
            header,
            footer,
            chrom_ids,
            zooms: HashMap::new(),
            sites: HashMap::new(),
            inflater: Some(Inflater::new()),
            raw: Vec::new(),
            block: Vec::new(),
            blocks: Vec::new(),
            entries: Vec::new(),
            entries_ready: false,
        })
    }

    pub fn header(&self) -> &Header {
        &self.header
    }

    pub fn footer(&self) -> &Footer {
        &self.footer
    }

    /// Index of a chromosome by name.
    pub fn chrom_index(&self, name: &str) -> Result<u32> {
        self.chrom_ids
            .get(name)
            .copied()
            .ok_or_else(|| Error::UnknownChrom(name.to_owned()))
    }

    /// A chromosome's fragment count (its number of restriction sites).
    pub fn fragment_count(&self, chrom: &str) -> Result<u32> {
        let c = self.chrom_index(chrom)?;
        self.header
            .sites
            .get(c as usize)
            .map(|s| s.count)
            .ok_or_else(|| missing("fragment sites", format!("for {chrom}")))
    }

    /// A chromosome's restriction site positions, read once.
    pub fn sites(&mut self, chrom: u32) -> Result<&[u32]> {
        if !self.sites.contains_key(&chrom) {
            let list = *self
                .header
                .sites
                .get(chrom as usize)
                .ok_or_else(|| missing("fragment sites", format!("for chromosome {chrom}")))?;
            let buf = self.read_range(list.position, list.count as u64 * 4, "fragment sites")?;
            let sites = buf
                .chunks_exact(4)
                .map(|b| {
                    u32::try_from(i32::from_le_bytes(b.try_into().unwrap()))
                        .map_err(|_| Error::invalid("fragment sites", "negative position"))
                })
                .collect::<Result<Vec<_>>>()?;
            if sites.windows(2).any(|w| w[0] > w[1]) {
                return Err(Error::invalid("fragment sites", "unsorted positions"));
            }
            self.sites.insert(chrom, sites);
        }
        Ok(&self.sites[&chrom])
    }

    /// The block index of a chromosome pair (`chr1 <= chr2`) at `bin_size`
    /// in `unit`, read once; `None` when the file has no such map.
    pub fn zoom(
        &mut self,
        chr1: u32,
        chr2: u32,
        unit: ContactUnit,
        bin_size: u32,
    ) -> Result<Option<&ZoomData>> {
        let key = (chr1, chr2, unit, bin_size);
        if !self.zooms.contains_key(&key) {
            let zoom = match self.footer.matrices.get(&(chr1, chr2)) {
                Some(&position) => ZoomData::parse(
                    &mut Cursor::new(&self.source, position),
                    position,
                    unit,
                    bin_size,
                )?,
                None => None,
            };
            if zoom
                .as_ref()
                .is_some_and(|z| (z.chr1, z.chr2) != (chr1, chr2))
            {
                return Err(Error::invalid(
                    "hic matrix",
                    "chromosome pair disagrees with master index",
                ));
            }
            self.zooms.insert(key, zoom);
        }
        Ok(self.zooms[&key].as_ref())
    }

    fn read_range(&self, position: u64, nbytes: u64, what: &'static str) -> Result<Vec<u8>> {
        if position
            .checked_add(nbytes)
            .is_none_or(|end| end > self.source.len())
        {
            return Err(Error::invalid(what, "range outside source"));
        }
        let n = usize::try_from(nbytes).map_err(|_| Error::invalid(what, "range exceeds usize"))?;
        let mut buf = vec![0u8; n];
        self.source.read_exact_at(&mut buf, position)?;
        Ok(buf)
    }

    /// Resolves the query and reads its block index.
    fn prepare(&mut self, query: &ContactQuery<'_>) -> Result<Plan> {
        query.check()?;
        let (unit, res) = (query.unit, query.resolution);
        let (r, c) = (
            self.chrom_index(query.rows.chrom)?,
            self.chrom_index(query.cols.chrom)?,
        );
        let plan = if r <= c {
            Plan {
                chr1: r,
                chr2: c,
                x: query.row_bins(),
                y: query.col_bins(),
                transposed: false,
            }
        } else {
            Plan {
                chr1: c,
                chr2: r,
                x: query.col_bins(),
                y: query.row_bins(),
                transposed: true,
            }
        };
        let stored = self.zoom(plan.chr1, plan.chr2, unit, res)?.is_some();
        if !stored && !self.header.resolutions(unit).contains(&res) {
            return Err(missing("resolution", format!("{res} {}", unit.as_str())));
        }
        Ok(plan)
    }

    /// Calls `f(row_bin, col_bin, count)` for every stored cell of the tile
    /// with a finite, nonzero count, decoding each overlapping block once.
    fn visit(
        &mut self,
        query: &ContactQuery<'_>,
        plan: &Plan,
        mut f: impl FnMut(u32, u32, f32),
    ) -> Result<()> {
        if plan.x.is_empty() || plan.y.is_empty() {
            return Ok(());
        }
        let (unit, res) = (query.unit, query.resolution);
        let version = self.header.version;
        let mut blocks = std::mem::take(&mut self.blocks);
        if let Some(zoom) = &self.zooms[&(plan.chr1, plan.chr2, unit, res)] {
            zoom.blocks_for_into(
                version,
                self.depth,
                plan.x.clone(),
                plan.y.clone(),
                &mut blocks,
            );
        } else {
            blocks.clear();
        }
        let mut inflater = self.inflater.take().unwrap_or_default();
        let mut raw = std::mem::take(&mut self.raw);
        let mut block = std::mem::take(&mut self.block);
        let result = (|| {
            let intra = plan.chr1 == plan.chr2;
            let mut emit = |x: i32, y: i32, v: f32| {
                let (Ok(x), Ok(y)) = (u32::try_from(x), u32::try_from(y)) else {
                    return;
                };
                if !plan.x.contains(&x) || !plan.y.contains(&y) {
                    return;
                }
                if !v.is_finite() || v == 0.0 {
                    return;
                }
                if plan.transposed {
                    f(y, x, v)
                } else {
                    f(x, y, v)
                }
            };
            let mut first = 0;
            while first < blocks.len() {
                let start = blocks[first].position;
                let end_of = |entry: &BlockEntry| {
                    entry
                        .position
                        .checked_add(entry.size as u64)
                        .filter(|&end| end <= self.source.len())
                        .ok_or_else(|| Error::invalid("hic block", "file range outside source"))
                };
                let mut end = end_of(&blocks[first])?;
                let mut last = first + 1;
                while last < blocks.len() && blocks[last].position == end {
                    let next = end_of(&blocks[last])?;
                    if next - start > 1 << 20 {
                        break;
                    }
                    end = next;
                    last += 1;
                }
                let len = usize::try_from(end - start)
                    .map_err(|_| Error::invalid("hic block", "read exceeds usize"))?;
                raw.resize(len, 0);
                self.source.read_exact_at(&mut raw, start)?;
                for entry in &blocks[first..last] {
                    let lo = (entry.position - start) as usize;
                    inflater.inflate(&raw[lo..lo + entry.size as usize], &mut block)?;
                    decode_block(version, &block, |x, y, v| {
                        emit(x, y, v);
                        if intra && x != y {
                            emit(y, x, v);
                        }
                    })?;
                }
                first = last;
            }
            Ok(())
        })();
        self.blocks = blocks;
        self.inflater = Some(inflater);
        self.raw = raw;
        self.block = block;
        result
    }
}

impl<S: ReadAt> ContactRead for HicReader<S> {
    fn chrom_len(&self, chrom: &str) -> Result<u32> {
        let len = self.header.chroms[self.chrom_index(chrom)? as usize].len;
        u32::try_from(len)
            .map_err(|_| Error::invalid("chromosome", format!("{chrom} longer than 2^32 - 1")))
    }

    fn resolutions(&self, unit: ContactUnit) -> &[u32] {
        self.header.resolutions(unit)
    }

    fn matrix(&mut self, query: ContactQuery<'_>, out: &mut [f32]) -> Result<()> {
        query.check()?;
        let (rows, cols) = query.shape();
        let cells = rows
            .checked_mul(cols)
            .ok_or_else(|| Error::invalid("contact matrix", "shape overflows usize"))?;
        check_output("contact matrix", out.len(), cells)?;
        let plan = self.prepare(&query)?;
        if cells == 0 {
            return Ok(());
        }
        let (r0, c0) = (query.row_bins().start, query.col_bins().start);
        out.fill(0.0);
        self.visit(&query, &plan, |r, c, v| {
            let at = (r - r0) as usize * cols + (c - c0) as usize;
            // SAFETY: visit clips emitted bins to the query. The checked
            // product rows*cols equals out.len(), so at is in range.
            unsafe {
                *out.get_unchecked_mut(at) = v;
            }
        })
    }

    fn query_entries(&mut self, query: ContactQuery<'_>) -> Result<usize> {
        self.entries_ready = false;
        let plan = self.prepare(&query)?;
        let mut staged = std::mem::take(&mut self.entries);
        staged.clear();
        let res = query.resolution;
        let result = self.visit(&query, &plan, |r, c, value| {
            staged.push(ContactEntry {
                row: r.saturating_mul(res),
                col: c.saturating_mul(res),
                value,
            })
        });
        self.entries = staged;
        result?;
        let staged = &mut self.entries;
        staged.sort_unstable_by_key(|e| (e.row, e.col));
        let n = staged.len();
        self.entries_ready = true;
        Ok(n)
    }

    fn write_entries(&self, out: &mut [ContactEntry]) -> Result<()> {
        if !self.entries_ready {
            return Err(Error::NotQueried);
        }
        let staged = &self.entries;
        check_output("contact entries", out.len(), staged.len())?;
        out.copy_from_slice(staged);
        Ok(())
    }
}
